# -*- coding: utf-8 -*-
# The template for batch run of DFTTK
import argparse
import datetime
from pymatgen.ext.matproj import MPRester, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.inputs import Potcar
#from dfttk.wflows import get_wf_gibbs, get_wf_EV_bjb, get_wf_gibbs_robust
from dfttk.wflows import get_wf_gibbs
from dfttk.utils import recursive_glob, sort_x_by_y
from dfttk.structure_builders.parse_anrl_prototype import multi_replace
from monty.serialization import loadfn, dumpfn
import dfttk.pythelec as pythelec
from dfttk.pythelec import thelecMDB
import warnings
import copy
import os
import sys
import json
import subprocess
import shutil
import numpy as np
from fireworks.fw_config import config_to_dict
from monty.serialization import loadfn
from atomate.vasp.database import VaspCalcDb
from dfttk.analysis.ywutils import formula2composition, reduced_formula, get_used_pot, get_Magnetic_State

def findjobdir(jobpath, metatag):
    try:
        if not os.path.isdir(jobpath): return None
    except:
        return None

    for dir in os.listdir(jobpath):
        dirpath = jobpath+'/'+dir
        if os.path.isdir(dirpath):
            for file in os.listdir(dirpath):
                if file.startswith("METADATA"):
                    with open(os.path.join(dirpath, file),'r') as fp:
                        lines = fp.readlines()
                        for line in lines:
                            try:
                                line.index(metatag)
                                return dir
                            except:
                                pass
    return None

class thfindMDB ():
    """
    find the metadata tag that has finished.

    Parameters
        STR_FOLDER = args.STRUCTURE_FOLDER
            folder/file containing structures
        MATCH_PATTERN = args.MATCH_PATTERN
            Match patterns for structure file, e.g. *POSCAR
        RECURSIVE = args.RECURSIVE
            recursive or not
        WORKFLOW = args.WORKFLOW
            workflow, current only get_wf_gibbs
    """
    def __init__(self, args, vasp_db):
        self.check = args.check
        self.remove = args.remove

        self.tags = []
        self._Yphon = []
        self.within = []
        self.containall = []
        self.containany = []
        self.excludeall = []
        self.excludeany = []
        self.nV = args.nV
        self.metatag = args.metatag
        self.get = args.get
        self.supercellN = args.supercellN
        self.t0 = args.t0
        self.t1 = args.t1
        self.td = args.td
        self.jobpath = args.jobpath
        self.findbandgap = args.findbandgap
        self.db_repair = args.db_repair
        self.db_renew = args.db_renew
        if args.within is not None: self.within, tmp = formula2composition(args.within)
        if args.containall is not None: self.containall, tmp = formula2composition(args.containall)
        if args.containany is not None: self.containany, tmp = formula2composition(args.containany)
        if args.excludeall is not None: self.excludeall, tmp = formula2composition(args.excludeall)
        if args.excludeany is not None: self.excludeany, tmp = formula2composition(args.excludeany)

        self.plotonly = args.plotonly
        if args.qhamode is not None:
            self.qhamode = args.qhamode
        else:
            self.qhamode = 'phonon'
        if args.qhamode == 'debye' : self.qhamode = 'qha'

        search_condition = [{'output.structure.lattice.volume': {'$exists': True}}]
        if len(self.containall)!=0:
            search_condition.append({"elements":{"$all":self.containall}})
        if len(self.containany)!=0:
            search_condition.append({"elements":{"$in":self.containany}})
        if len(self.excludeany)!=0:
            search_condition.append({"elements":{"$nin":self.excludeany}})

        scondition = []
        if len(search_condition) > 1:
            try:
                metadata_list =  vasp_db.collection.find({'$and':search_condition},{'metadata':1})
                metadata_list = list(set([i['metadata']['tag'] for i in metadata_list]))
                metadata_list = [{'tag':i} for i in metadata_list]
                scondition.append({"metadata":{"$in":metadata_list}})
            except:
                pass

        if not self.plotonly:
            try:
                self.vasp_db = vasp_db
                if self.qhamode=='phonon':
                    scondition.append({"S_vib": { "$exists": True } })
                    self.items = list((self.vasp_db).db['phonon'].find({'$and':scondition},\
                        {'metadata':1, 'unitcell':1, 'volume':1, 'supercell_matrix':1}))
                else:
                    scondition.append({"debye": { "$exists": True } })
                    self.items = list((self.vasp_db).db['qha'].find({'$and':scondition},\
                        {'metadata':1, 'structure':1}))
            except:
                self.vasp_db = None
                warnings.warn("\n*********WARNING: CANNOT get MongoDB service, so I will proceed using local data")

    def skipby(self, phase, metatag):
        if self.metatag!=None:
            if self.metatag!=metatag: return True
        try:
            els,tmp = formula2composition(phase.split('_')[0])
        except:
            return True
        if len (self.within) != 0:
            for e in els:
                if e not in self.within: return True
        if len (self.excludeall) != 0:
            for e in self.excludeall:
                rx = e not in els
                if rx: break
            if not rx: return True
        if len (self.excludeany) != 0:
            for e in self.excludeany:
                return True
        if len (self.containall) != 0:
            for e in self.containall:
               if e not in els: return True
            return False
        if len (self.containany) != 0:
            for e in self.containany:
                if e in els: return False
            return True
        return False


    def find_plotfiles(self):
        if self.jobpath!=None:
            jobpath = os.listdir(self.jobpath)
        else:
            self.jobpath='./'
            jobpath = os.listdir('./')

        print("\nfound complete local calculations ...\n")
        for _dir in jobpath:
            dir = self.jobpath+_dir
            thermofile = dir+"/fvib_ele"
            if not os.path.exists(thermofile): continue
            volumes = None
            energies = None
            ss = [s for s in dir.split('/') if s!=""]
            formula = ss[-1].split('_')[0]
            #print ("eeeeeeee", dir, ss, formula)
            metatag = None
            readme = dir+"/readme"
            if os.path.exists(readme):
                with open(readme, "r") as fp:
                    info = json.load(fp)
                    try:
                        metatag = info['METADATA']['tag']
                        volumes = info['E-V']['volumes']
                        energies = info['E-V']['energies']
                    except:
                        #print (readme, "FILED")
                        pass
            if self.skipby(formula, metatag): continue
            #print (self.metatag,metatag)
            sys.stdout.write('{}, dir: {}, formula: {}\n'.format(metatag, dir, formula))
            self.tags.append([metatag, thermofile, volumes, energies, dir, formula])


    def run_console(self):
        if self.plotonly: self.find_plotfiles()
        elif self.vasp_db==None: self.find_plotfiles()
        elif self.check: self.check_find()
        elif self.db_repair or self.db_renew: self.qha_renew()
        elif self.qhamode=='phonon': self.phonon_find()
        elif self.qhamode=='qha': self.debye_find()
        return self.tags


    def phonon_find(self):
        hit = []
        count = []
        phases = []
        volumes = []
        ITEMS = []
        self.supercellsize = []
        for i in self.items:
            mm = i['metadata']
            if mm in hit:
                if i['volume'] not in volumes[hit.index(mm)]:
                    volumes[hit.index(mm)].append(i['volume'])
                    count[hit.index(mm)] += 1
            else:
                ITEMS.append(i)
                hit.append(mm)
                count.append(1)
                volumes.append([i['volume']])
                structure = Structure.from_dict(i['unitcell'])
                natoms = len(structure.sites)
                supercell_matrix = i['supercell_matrix']
                self.supercellsize.append(natoms*int(np.linalg.det(np.array(supercell_matrix))+.5))
                formula_pretty = structure.composition.reduced_formula
                try:
                    formula2composition(formula_pretty)
                except:
                    formula_pretty = reduced_formula(structure.composition.alphabetical_formula)

                sa = SpacegroupAnalyzer(structure)
                phasename = formula_pretty+'_'\
                    + sa.get_space_group_symbol().replace('/','.')+'_'+str(sa.get_space_group_number())
                if phasename in phases:
                    for jj in range (10000):
                        nphasename = phasename + "#" + str(jj)
                        if nphasename in phases: continue
                        phasename = nphasename
                        break
                phases.append(phasename)

        print("\nfound complete calculations in the collection:", self.qhamode, "\n")
        total = 0
        total_qha_phonon = 0
        all_static_calculations = list((self.vasp_db).db['tasks'].\
            find({'$and':[{'metadata': { "$exists": True }}, {'adopted': True} ]},\
            {'metadata':1, 'output':1, 'input':1, 'orig_inputs':1}))
        all_qha_calculations = list((self.vasp_db).db['qha'].\
            find({'$and':[{'metadata': { "$exists": True }},{'has_phonon':True}]}, {'metadata':1, 'temperatures':1}))
        all_qha_phonon_calculations = list((self.vasp_db).db['qha_phonon'].\
            find({'$and':[{'metadata': { "$exists": True }},{'has_phonon':True}]}, {'metadata':1, 'temperatures':1}))
        for i,m in enumerate(hit):
            if self.skipby(phases[i], m['tag']): continue
            total += 1
            static_calculations = [f for f in all_static_calculations if f['metadata']['tag']==m['tag']]
            qha_calculations = [f for f in all_qha_calculations if f['metadata']['tag']==m['tag']]
            qha_phonon_calculations = [f for f in all_qha_phonon_calculations if f['metadata']['tag']==m['tag']]
            """
            static_calculations = [f for f in all_static_calculations if f['metadata']==m]
            qha_calculations = [f for f in all_qha_calculations if f['metadata']==m]
            qha_phonon_calculations = [f for f in all_qha_phonon_calculations if f['metadata']==m]
            """
            qha_phonon_success =True
            if len(qha_calculations) > 0:
                total_qha_phonon += 1
            elif len(qha_phonon_calculations) > 0:
                total_qha_phonon += 1
            else:
                qha_phonon_success = False

            gapfound = False
            potsoc = None
            _volumes = []
            emin = 1.e36
            for ii, calc in enumerate(static_calculations):
                vol = calc['output']['structure']['lattice']['volume']
                ene = calc['output']['energy']
                if ene < emin: _calc = calc
                if vol not in _volumes: _volumes.append(vol)
                bandgap = calc['output']['bandgap']
                if not gapfound: gapfound = float(bandgap) > 0.0
            potsoc = get_used_pot(_calc)
            Mag_State = get_Magnetic_State(_calc)
            pname = phases[i].split('#')
            if len(pname)>1: phases[i] = pname[0]+potsoc+Mag_State+'#'+pname[1]
            else: phases[i] = pname[0]+potsoc+Mag_State
            nS = len(_volumes)
            if self.findbandgap:
                if gapfound: sys.stdout.write('{}, phonon: {:>2}, static: {:>2}, supercellsize: {:>3}, {}\n'.format(m, count[i], nS, self.supercellsize[i], phases[i]))
            else:
                if count[i] < self.nV: continue
                if self.db_repair:
                    if qha_phonon_success and not self.db_renew: continue
                if self.supercellsize[i] < self.supercellN: continue
                jobpath = findjobdir(self.jobpath, m['tag'])
                if self.remove:
                    sys.stdout.write('dfttk db_remove --force -m all -tag {} phonon: {:>2}, static: {:>2}, SN: {:>3}, qha_phonon: {:<1.1s}, {}\n'\
                        .format(m['tag'], count[i], nS, self.supercellsize[i], str(qha_phonon_success), phases[i]))
                elif jobpath==None:
                    sys.stdout.write('{}, phonon: {:>2}, static: {:>2}, SN: {:>3}, qha_phonon: {:<1.1s}, {}\n'\
                        .format(m, count[i], nS, self.supercellsize[i], str(qha_phonon_success), phases[i]))
                else:
                    sys.stdout.write('{}, phonon: {:>2}, static: {:>2}, SN: {:>3}, qha_phonon: {:<1.1s}, {},{}\n'\
                        .format(m, count[i], nS, self.supercellsize[i], str(qha_phonon_success), phases[i],jobpath))
                #if count[i]>=5: self.tags.append({'tag':m['tag'],'phasename':phases[i]})
                self.tags.append({'tag':m['tag'],'phasename':phases[i]})
                #print(sorted(volumes[i]))
        sys.stdout.write ('\n{}/{} qha_phonon successful under the given searching conditions.\n'\
            .format(total_qha_phonon, total))


    def qha_renew(self):
        hit = []
        phases = []
        static_collection = (self.vasp_db).collection.find({'$and':[{'metadata': { "$exists": True }}, \
            {'adopted': True} ]})
        for i in static_collection:
            mm = i['metadata']
            if mm in hit: continue
            if len(mm)>1: continue
            else:
                hit.append(mm)
                structure = Structure.from_dict(i['output']['structure'])
                formula_pretty = structure.composition.reduced_formula
                try:
                    formula2composition(formula_pretty)
                except:
                    formula_pretty = reduced_formula(structure.composition.alphabetical_formula)

                sa = SpacegroupAnalyzer(structure)
                phasename = formula_pretty+'_'\
                    + sa.get_space_group_symbol().replace('/','.')+'_'+str(sa.get_space_group_number())
                if phasename in phases:
                    for jj in range (10000):
                        nphasename = phasename + "#" + str(jj)
                        if nphasename in phases: continue
                        phasename = nphasename
                        break
                phases.append(phasename)

        print("\nfound complete calculations in the task collection:\n")
        total = 0
        total_qha = 0
        total_qha_phonon = 0
        all_static_calculations = list((self.vasp_db).db['tasks'].\
            find({'$and':[{'metadata': { "$exists": True }}, {'adopted': True} ]},\
            {'metadata':1, 'output':1, 'input':1, 'orig_inputs':1}))
        all_qha_calculations = list((self.vasp_db).db['qha'].\
            find({'$and':[{'metadata': { "$exists": True }},{'has_phonon':True}]}, {'metadata':1, 'temperatures':1}))
        all_qha_phonon_calculations = list((self.vasp_db).db['qha_phonon'].\
            find({'$and':[{'metadata': { "$exists": True }},{'has_phonon':True}]}, {'metadata':1, 'temperatures':1}))
        for i,m in enumerate(hit):
            if self.skipby(phases[i], m['tag']): continue
            total += 1
            static_calculations = [f for f in all_static_calculations if f['metadata']['tag']==m['tag']]
            qha_calculations = [f for f in all_qha_calculations if f['metadata']['tag']==m['tag']]
            qha_phonon_calculations = [f for f in all_qha_phonon_calculations if f['metadata']['tag']==m['tag']]
            qha_phonon_success = len(qha_phonon_calculations) > 0
            if qha_phonon_success: total_qha_phonon += 1
            if len(qha_calculations) > 0 or qha_phonon_success: total_qha += 1

            potsoc = None
            volumes = []
            energies = []
            for ii, calc in enumerate(static_calculations):
                vol = calc['output']['structure']['lattice']['volume']
                if vol in volumes:
                    if len (calc['metadata']) > 1: continue
                    else:
                        ix = volumes.index(vol)
                        volumes.pop(ix)
                        energies.pop(ix)
                volumes.append(calc['output']['structure']['lattice']['volume'])
                energies.append(calc['output']['energy'])
                if potsoc is None:
                    potsoc = get_used_pot(calc)
                    pname = phases[i].split('#')
                    if len(pname)>1: phases[i] = pname[0]+potsoc+'#'+pname[1]
                    else: phases[i] = pname[0]+potsoc
            nS = len(volumes)
            if nS < 6: continue
            if qha_phonon_success and not self.db_renew: continue
            energies = sort_x_by_y(energies, volumes)
            volumes = sorted(volumes)
            volumes = np.array(volumes)
            energies = np.array(energies)
            val, idx = min((val, idx) for (idx, val) in enumerate(energies))
            if idx <2 or idx>nS-2: continue

            jobpath = findjobdir(self.jobpath, m['tag'])
            if jobpath==None:
                sys.stdout.write('{}, static: {:>2}, qha_phonon: {:<1.1s}, {}\n'\
                    .format(m, nS, str(qha_phonon_success), phases[i]))
            else:
                sys.stdout.write('{}, static: {:>2}, qha_phonon: {:<1.1s}, {},{}\n'\
                    .format(m, nS, str(qha_phonon_success), phases[i],jobpath))
            self.tags.append({'tag':m['tag'],'phasename':phases[i]})
        sys.stdout.write ('\n({},{})/{} (qha, qha_phonon) entries returned under the given searching conditions.\n'\
            .format(total_qha, total_qha_phonon, total))


    def debye_find(self):
        hit = []
        phases = []
        count = []
        for i in self.items:
            mm = i['metadata']
            if mm in hit:
                count[hit.index(mm)] += 1
            else:
                hit.append(mm)
                count.append(1)
                structure = Structure.from_dict(i['structure'])
                formula_pretty = structure.composition.reduced_formula
                try:
                    formula2composition(formula_pretty)
                except:
                    formula_pretty = reduced_formula(structure.composition.alphabetical_formula)
                sa = SpacegroupAnalyzer(structure)
                phasename = formula_pretty+'_'\
                    + sa.get_space_group_symbol().replace('/','.')+'_'+str(sa.get_space_group_number())
                if phasename in phases:
                    for jj in range (10000):
                        nphasename = phasename + "#" + str(jj)
                        if nphasename in phases: continue
                        phasename = nphasename
                        break
                phases.append(phasename)

        print("\nfound complete calculations in the collection:", self.qhamode, "\n")
        all_static_calculations = list((self.vasp_db).db['tasks'].\
            find({'$and':[{'metadata': { "$exists": True }}, {'adopted': True} ]},\
            {'metadata':1, 'output':1, 'input':1, 'orig_inputs':1}))
        for i,m in enumerate(hit):
            if self.skipby(phases[i], m['tag']): continue
            static_calculations = [f for f in all_static_calculations if f['metadata']['tag']==m['tag']]
            for ii, calc in enumerate(static_calculations):
                potsoc = get_used_pot(calc)
                if self.qhamode == 'qha' : potsoc += "_debye"
                pname = phases[i].split('#')
                if len(pname)>1: phases[i] = pname[0]+potsoc+'#'+pname[1]
                else: phases[i] = pname[0]+potsoc
                break
            print (m, ":", phases[i])
            self.tags.append({'tag':m['tag'],'phasename':phases[i]})


    def check_find(self):
        """
        hit = []
        relaxations_collection = (self.vasp_db).db['relaxations'].find({})
        for i in relaxations_collection:
            try:
                mm = i['metadata']['tag']
                hit.append(mm)
            except:
                continue

        lastupdated = [None] * len(hit)
        for i,mm in enumerate(hit):
            static_calc = (self.vasp_db).collection.\
                find({'$and':[ {'metadata.tag': mm} ]})
            for calc in static_calc:
                lnew = calc['last_updated']
                if lastupdated[i]!=None:
                    lold = lastupdated[i]
                    if lnew > lold: lastupdated[i] = lnew
                else:
                    lastupdated[i] = lnew

        """


        print("\nfinding tags for complete calculations in the static collection\n")
        hit = []
        lastupdated = []
        static_collection = (self.vasp_db).collection.find({})
        for i in static_collection:
            try:
                mm = i['metadata']['tag']
            except:
                continue
            if mm in hit:
                idx = hit.index(mm)
                lold = lastupdated[idx]
                lnew = i['last_updated']
                if lnew > lold: lastupdated[idx] = lnew
            else:
                lastupdated.append(i['last_updated'])
                hit.append(mm)


        print("\nfinding complete calculations in the phonon collection\n")
        phases =  [""] * len(hit)
        supercellsize =  [0] * len(hit)
        phonon_count = [0] * len(hit)
        phonon_calc = list((self.vasp_db).db['phonon'].find({"S_vib": { "$exists": True } },\
            {'metadata':1, 'unitcell':1, 'supercell_matrix':1}))

        for i,mm in enumerate(hit):
            for calc in phonon_calc:
                if calc['metadata']['tag']!=mm: continue
                phonon_count[i] += 1
                if phonon_count[i]==1:
                    structure = Structure.from_dict(calc['unitcell'])
                    natoms = len(structure.sites)
                    supercell_matrix = calc['supercell_matrix']
                    supercellsize[i]= (natoms*int(np.linalg.det(np.array(supercell_matrix))+.5))
                    formula_pretty = structure.composition.reduced_formula
                    try:
                        formula2composition(formula_pretty)
                    except:
                        formula_pretty = reduced_formula(structure.composition.alphabetical_formula)
                    sa = SpacegroupAnalyzer(structure)
                    phasename = formula_pretty+'_'\
                        + sa.get_space_group_symbol().replace('/','.')+'_'+str(sa.get_space_group_number())
                    if phasename in phases:
                        for jj in range (10000):
                            nphasename = phasename + "#" + str(jj)
                            if nphasename in phases: continue
                            phasename = nphasename
                            break
                    phases[i] = phasename


        print("\nfinding complete calculations in the static collection\n")
        static_count = [0] * len(hit)
        for i,mm in enumerate(hit):
            if self.skipby(phases[i], mm): continue
            static_calc = (self.vasp_db).collection.\
                find({'$and':[ {'metadata.tag': mm} ]})
            for calc in static_calc:
                static_count[i] += 1

        print("\nfinding complete calculations in the qha collection\n")

        qha_count = [0] * len(hit)
        for i,mm in enumerate(hit):
            qha_calc = (self.vasp_db).db['qha'].\
                find({'$and':[ {'metadata.tag': mm} ]})
            for calc in qha_calc:
                qha_count[i] += 1

        print("\nfinding complete calculations in the qha_phonon collection\n")

        qha_phonon_count = [0] * len(hit)
        for i,mm in enumerate(hit):
            try:
                qha_phonon_calculations = self.vasp_db.db['qha_phonon'].find({'metadata.tag': mm})
                T = qha_phonon_calculations[0]['phonon']['temperatures']
                qha_phonon_count[i] += 1
            except:
                try:
                    qha_phonon_calculations = self.vasp_db.db['qha'].find({'metadata.tag': mm})
                    T = qha_phonon_calculations[0]['phonon']['temperatures']
                    qha_phonon_count[i] += 1
                except:
                    pass

        print("\nfinding complete calculations in the relaxations collection\n")

        relaxations_count = [0] * len(hit)
        for i,mm in enumerate(hit):
            relaxations_calc = (self.vasp_db).db['relaxations'].\
                find({'$and':[ {'metadata.tag': mm} ]})
            for calc in relaxations_calc:
                relaxations_count[i] += 1

        nTBD = 0
        for i,mm in enumerate(hit):
                #dd = datetime.datetime.strptime(lastupdated[i], '%Y-%m-%d %H:%M:%S.%f').date()
                dd = lastupdated[i].date()
                now = datetime.datetime.now().date()
                #if supercellsize[i]>=16 and phonon_count[i]>=5: continue
                if supercellsize[i]>=self.supercellN and phonon_count[i]>=self.nV: continue
                if dd >now-datetime.timedelta(days=7): continue
                nTBD += 1
                sys.stdout.write('[{:>04}] relax: {:>2}, static: {:>2}, qha: {:>2}, qha_phonon: {:>2}, phonon: {:>2}, SN: {:>3}, phases: {}, {}\n'.format(i, relaxations_count[i], static_count[i], qha_count[i], qha_phonon_count[i], phonon_count[i], supercellsize[i], phases[i], dd))
                #sys.stdout.write('{}, static: {:>2}, qha: {:>2}, qha_phonon: {:>2}, phonon: {:>2}, SN: {:>3}, phases: {}, date: {}\n'.format(mm['tag'], static_count[i], qha_count[i], qha_phonon_count[i], phonon_count[i], supercellsize[i], phases[i], lastupdated[i]))
                self.tags.append({'tag':mm,'phasename':phases[i]})

        print("\n", nTBD,"/", len(hit), "recs to be removed\n")

        for t in self.tags:
            print(t)

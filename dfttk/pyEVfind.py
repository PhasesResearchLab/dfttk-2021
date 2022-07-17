# -*- coding: utf-8 -*-
# The template for batch run of DFTTK
import argparse
from pymatgen.ext.matproj import MPRester, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.inputs import Potcar
#from dfttk.wflows import get_wf_gibbs, get_wf_EV_bjb, get_wf_gibbs_robust
from dfttk.wflows import get_wf_gibbs
from dfttk.utils import recursive_glob
from dfttk.structure_builders.parse_anrl_prototype import multi_replace
from monty.serialization import loadfn, dumpfn
from dfttk.pythelec import thelecMDB, vol_within
import warnings
import copy
import os
import sys
import subprocess
import shutil
import numpy as np
from fireworks.fw_config import config_to_dict
from atomate.vasp.database import VaspCalcDb
from dfttk.analysis.ywutils import formula2composition, reduced_formula
from dfttk.analysis.ywplot import myjsonout, thermoplot
from dfttk.utils import sort_x_by_y
from dfttk.analysis.ywutils import get_rec_from_metatag, get_used_pot

class EVfindMDB ():
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
        self.vasp_db = vasp_db
        self.qhamode = args.qhamode
        #self.items = (self.vasp_db).collection.find({'adopted': True})
        #print(self.hit_condition)
        self.within = []
        self.containall = []
        self.containany = []
        self.excludeall = []
        self.excludeany = []
        self.nV = args.nV
        self.print = args.print
        self.plot = args.plot
        self.natoms = args.natoms
        self.findbandgap = args.findbandgap
        self.metatag = args.metatag
        if args.within is not None: self.within, tmp = formula2composition(args.within)
        if args.containall is not None: self.containall, tmp = formula2composition(args.containall)
        if args.containany is not None: self.containany, tmp = formula2composition(args.containany)
        if args.excludeall is not None: self.excludeall, tmp = formula2composition(args.excludeall)
        if args.excludeany is not None: self.excludeany, tmp = formula2composition(args.excludeany)

        search_condition = [{'output.structure.lattice.volume': {'$exists': True}}, {'metadata': {'$exists': True}}]
        if len(self.containall)!=0:
            search_condition.append({"elements":{"$all":self.containall}})
        if len(self.containany)!=0:
            search_condition.append({"elements":{"$in":self.containany}})
        if len(self.excludeany)!=0:
            search_condition.append({"elements":{"$nin":self.excludeany}})
        #print (search_condition)
        self.items = (self.vasp_db).collection.find({'$and':search_condition},\
            {'metadata':1, 'output':1, 'input':1, 'orig_inputs':1, 'elements':1})

        scondition = []
        if len(search_condition) > 1:
            metadata_list =  (self.vasp_db).collection.find({'$and':search_condition},{'metadata':1})
            metadata_list = list(set([i['metadata']['tag'] for i in metadata_list]))
            metadata_list = [{'tag':i} for i in metadata_list]
            scondition.append({"metadata":{"$in":metadata_list}})
        
        if self.qhamode=='phonon':
            scondition.extend([{'adopted': True}, {"S_vib": { "$exists": True } }])
            hit_condition = list((self.vasp_db).db['phonon'].find({'$and':scondition},{'metadata':1}))
        elif self.qhamode=='debye':
            scondition.append({"debye": { "$exists": True } })
            hit_condition = list((self.vasp_db).db['qha'].find({'$and':scondition},{'metadata':1}))
        else:
            hit_condition = None
            
        if hit_condition is not None:
            metadata = [i['metadata']['tag'] for i in hit_condition]
            self.hit_condition = list(set(metadata))
            self.hit_count = {i:metadata.count(i) for i in metadata}
        else: self.hit_condition = None

    def skipby(self, phase, metatag, els=None):
        if self.metatag!=None:
            if self.metatag!=metatag: return True
        if els is None:
            els,tmp = formula2composition(phase.split('_')[0])
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

    def run_console(self):
        self.EV_find()


    def EV_find(self):
        evdirhome = 'E-V'
        if not os.path.exists(evdirhome): os.mkdir(evdirhome)

        hit = []
        count = []
        phases = []
        volumes = []
        ITEMS = []
        potname = []
        fp_ev = open(os.path.join(evdirhome,"E-V.dat"),"w")
        for i in self.items:
            mm = i['metadata']['tag']
            els = i['elements']
            if self.skipby("", mm, els=els): continue
            if self.hit_condition is not None:
                if mm not in self.hit_condition: continue
            if mm in hit:
                volume = i['output']['structure']['lattice']['volume']
                if volume not in volumes[hit.index(mm)]:
                    volumes[hit.index(mm)].append(volume)
                    count[hit.index(mm)] += 1
            else:
                ITEMS.append(i)
                hit.append(mm)
                count.append(1)
                volumes.append([i['output']['structure']['lattice']['volume']])

                potsoc = get_used_pot(i)

                structure = Structure.from_dict(i['output']['structure'])
                formula_pretty = structure.composition.reduced_formula
                try:
                    formula2composition(formula_pretty)
                except:
                    formula_pretty = reduced_formula(structure.composition.alphabetical_formula)
                sa = SpacegroupAnalyzer(structure)
                phasename = formula_pretty+'_'\
                    + sa.get_space_group_symbol().replace('/','.')+'_'+str(sa.get_space_group_number())
                potname.append(potsoc)

                if phasename in phases:
                    for jj in range (10000):
                        nphasename = phasename + "#" + str(jj)
                        if nphasename in phases: continue
                        phasename = nphasename
                        break
                phases.append(phasename)

        blank_lines = False
        for i,mm in enumerate(hit):
            if self.hit_condition is not None:
                if mm not in self.hit_condition: continue
                if self.qhamode == 'phonon':
                    if self.hit_count[mm] < self.nV: continue
            if count[i]<self.nV: continue
            if self.skipby(phases[i], mm): continue
            EV, POSCAR, INCAR = get_rec_from_metatag(self.vasp_db, mm)
            metadata = {'tag':mm}
            pname = phases[i].split('#')
            if len(pname)>1: phases[i] = pname[0]+potname[i]+EV['MagState']+'#'+pname[1]
            else: phases[i] = pname[0]+potname[i]+EV['MagState']
            if EV['natoms'] < self.natoms: continue
            
            sys.stdout.write('{}, static: {:>2}, natoms: {:>3}, {}\n'.format(metadata, count[i], EV['natoms'], phases[i]))

            folder = os.path.join(evdirhome,phases[i])
            if not os.path.exists(folder): os.mkdir(folder)
            with open (os.path.join(folder,'POSCAR'), 'w') as fp:
                fp.write(POSCAR)
            readme = {}
            readme['E-V'] = EV
            readme['INCAR'] = INCAR
            readme['POSCAR'] = POSCAR
            natoms = readme['E-V']['natoms']
            with open (os.path.join(folder,'readme'), 'w') as fp:
                myjsonout(readme, fp, indent="", comma="")
            i_energies = np.array(EV['energies'])/natoms
            i_volumes = np.array(EV['volumes'])/natoms
            val, idx = min((val, idx) for (idx, val) in enumerate(i_energies))
            if blank_lines: print("\n", file=fp_ev)
            blank_lines=True
            print("#phase:", phases[i], file=fp_ev)
            print("#metadata:", mm, file=fp_ev)
            print("#natoms:", natoms, file=fp_ev)
            elist = ['Fe','Cu', 'Se', 'Al', 'Ni', 'Co', 'Pt', 'Ta', 'O']
            el = []
            nel = []
            if len(EV['magmoms'])>0:
                fp_ev.write("#magmoms:")
                magmoms = EV['magmoms'][idx]
                m0 = magmoms[0]
                n0 = 1
                for j in range(1,len(magmoms)):
                    if magmoms[j] == m0:
                        n0 += 1
                    else:
                        if n0==1: fp_ev.write('{},'.format(m0))
                        else: fp_ev.write('{}*{},'.format(n0,m0))
                        idx = len(el)%len(elist)
                        el.append(elist[idx])
                        nel.append(n0)
                        n0 = 1
                        m0 = magmoms[j]
                if n0==1: fp_ev.write('{}\n'.format(m0))
                else: fp_ev.write('{}*{}\n'.format(n0,m0))
                idx = len(el)%len(elist)
                el.append(elist[idx])
                nel.append(n0)

                lines = [l for l in POSCAR.split('\n') if l!=""]
                with open (os.path.join(folder,phases[i]+'.VASP'), 'w') as fp:
                    for j in range(0,5):
                        print(lines[j], file=fp)
                    for j in range(len(el)):
                        fp.write(' {}'.format(el[j]))
                    fp.write('\n')
                    for j in range(len(nel)):
                        fp.write(' {}'.format(nel[j]))
                    fp.write('\n')
                    print(lines[7], file=fp)
                    for j in range(8,len(lines)):
                        print(lines[j], float(list(magmoms[j-8].values())[0]), file=fp)

            for j in range(len(i_volumes)):
                print(i_volumes[j], i_energies[j],  file=fp_ev)
            thermoplot(folder,"0 K total energies (eV/atom)", i_volumes, i_energies)

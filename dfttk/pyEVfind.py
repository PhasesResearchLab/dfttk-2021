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
        self.items = (self.vasp_db).collection.find({'output.structure.lattice.volume': {'$exists': True} })
        
        if self.qhamode=='phonon':
            hit_condition = list((self.vasp_db).db['phonon'].find({'$and':[{'adopted': True}, {"S_vib": { "$exists": True } }]},\
                {'metadata':1, 'volume':1}))
        elif self.qhamode=='debye':
            hit_condition = list((self.vasp_db).db['qha'].find({"debye": { "$exists": True } },\
                        {'metadata':1}))
        else:
            hit_condition = None
        if hit_condition is not None:
            self.hit_condition = []
            self.hit_count = []
            volumes = []
            for i in hit_condition:
                mm = i['metadata']['tag']
                if mm in self.hit_condition:
                    if self.qhamode == "phonon":
                        if i['volume'] not in volumes[self.hit_condition.index(mm)]:
                            volumes[self.hit_condition.index(mm)].append(i['volume'])
                            self.hit_count[self.hit_condition.index(mm)] += 1
                    else: self.hit_count[self.hit_condition.index(mm)] += 1
                else:
                    if self.qhamode == "phonon": volumes.append([i['volume']])
                    self.hit_condition.append(mm)
                    self.hit_count.append(1)
        else: self.hit_condition = None
        #print(self.hit_condition)
        self.within = []
        self.containall = []
        self.containany = []
        self.excludeall = []
        self.excludeany = []
        self.nV = args.nV
        self.print = args.print
        self.plot = args.plot
        self.findbandgap = args.findbandgap
        self.metatag = args.metatag
        if args.within is not None: self.within, tmp = formula2composition(args.within)
        if args.containall is not None: self.containall, tmp = formula2composition(args.containall)
        if args.containany is not None: self.containany, tmp = formula2composition(args.containany)
        if args.excludeall is not None: self.excludeall, tmp = formula2composition(args.excludeall)

    def skipby(self, phase, metatag):
        if self.metatag!=None:
            if self.metatag!=metatag: return True
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
        hit = []
        count = []
        phases = []
        volumes = []
        ITEMS = []
        potname = []
        for i in self.items:
            try:
                mm = i['metadata']['tag']
            except:
                continue
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

        for i,mm in enumerate(hit):
            if self.hit_condition is not None:
                if mm not in self.hit_condition: continue
                if self.qhamode == 'phonon':
                    if self.hit_count[self.hit_condition.index(mm)] < self.nV: continue
            if count[i]<self.nV: continue
            if self.skipby(phases[i], mm): continue
            if self.qhamode == 'phonon':
                if self.hit_count[self.hit_condition.index(mm)] < self.nV: continue
            metadata = {'tag':mm}
            pname = phases[i].split('#')
            if len(pname)>1: phases[i] = pname[0]+potname[i]+'#'+pname[1]
            else: phases[i] = pname[0]+potname[i]
            sys.stdout.write('{}, static: {:>2}, {}\n'.format(metadata, count[i], phases[i]))
            EV, POSCAR, INCAR = get_rec_from_metatag(self.vasp_db, mm)

            evdir = 'E-V'
            if not os.path.exists(evdir): os.mkdir(evdir)
            folder = os.path.join(evdir,phases[i])
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
            thermoplot(folder,"0 K total energies (eV/atom)", \
                np.array(EV['volumes'])/natoms, np.array(EV['energies'])/natoms)

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
            
            folder = os.path.join(evdirhome,phases[i])
            if not os.path.exists(folder): os.mkdir(folder)
            self.num_Born = 0
            self.get_data(folder, EV, mm)
            sys.stdout.write('{}, static: {:>2}, natoms: {:>3}, Born: {}, {}\n'.format(metadata, \
                count[i], EV['natoms'], self.num_Born, phases[i]))
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

    def get_data(self, phasename, EV, tag):
        phdir = os.path.join(phasename,'Yphon')
        if not os.path.exists(phdir):
            os.mkdir(phdir)

        self.num_Born = self.get_dielecfij(phdir, tag)
        for i in (self.vasp_db).db['phonon'].find({'metadata.tag': tag}):
            try:
                self.force_constant_factor = i['force_constant_factor']
            except:
                if self.static_vasp_version[0:1] >= '6':
                    self.force_constant_factor = 0.004091649655126895

            if i['volume'] not in EV['volumes']: continue
            voldir = self.get_superfij(i, phdir, EV['volumes'], EV['energies'])


    def get_superfij(self,i, phdir, volumes, energies):
        try:
            structure = Structure.from_dict(i['unitcell'])
        except:
            print("\nit seems phonon for", i['volume'],"is not finished yet  and so it is discared\n")
            return None
        vol = 'V{:010.6f}'.format(float(i['volume']))
        voldir = os.path.join(phdir,vol)
        if not os.path.exists(voldir):
            os.mkdir(voldir)
        elif os.path.exists(os.path.join(voldir,'superfij.out')): return voldir

        poscar = structure.to(fmt="poscar")
        unitcell_l = str(poscar).split('\n')
        natom = len(structure.sites)

        supercell_matrix = i['supercell_matrix']
        supercell_structure = copy.deepcopy(structure)
        supercell_structure.make_supercell(supercell_matrix)

        sa = SpacegroupAnalyzer(supercell_structure)
        #reduced_structure = supercell_structure.get_reduced_structure(reduction_algo='niggli')
        #print ('niggli reduced structure', reduced_structure)
        #poscar = reduced_structure.to(fmt="poscar")
        primitive_unitcell_structure = sa.find_primitive()
        poscar = primitive_unitcell_structure.to(fmt="poscar")
        punitcell_l = str(poscar).split('\n')

        natoms = len(supercell_structure.sites)
        ##print(supercell_structure.sites)
        poscar = supercell_structure.to(fmt="poscar")
        supercell_l = str(poscar).split('\n')
        structure.to(filename=os.path.join(voldir,'POSCAR'))

        with open (os.path.join(voldir,'OSZICAR'),'w') as out:
            out.write('   1 F= xx E0= {}\n'.format(energies[(list(volumes)).index(i['volume'])]))
        with open (os.path.join(voldir,'superfij.out'),'w') as out:
            for line in range (2,5):
                out.write('{}\n'.format(unitcell_l[line]))
            for line in range (2,5):
                out.write('{}\n'.format(supercell_l[line]))
            out.write('{} {}\n'.format(natoms, natoms//natom))
            for line in range (7,natoms+8):
                out.write('{}\n'.format(supercell_l[line]))
            force_constant_matrix = np.array(i['force_constants'])
            hessian_matrix = np.empty((natoms*3, natoms*3), dtype=float)

            for ii in range(natoms):
               for jj in range(natoms):
                    for x in range(3):
                        for y in range(3):
                            hessian_matrix[ii*3+x, jj*3+y] = -force_constant_matrix[ii,jj,x,y]

            hessian_matrix *= self.force_constant_factor
            if self.force_constant_factor!=1.0: 
                print ("\n force constant matrix has been rescaled by :", self.force_constant_factor)

            for xx in range(natoms*3):
                for yy in range(natoms*3-1):
                    out.write('{} '.format(hessian_matrix[xx,yy]))
                out.write('{}\n'.format(hessian_matrix[xx,natoms*3-1]))
        return voldir

    def get_dielecfij(self, phdir, tag):

        volumes_b = (self.vasp_db).db['borncharge'].find({'metadata.tag': tag}, {'_id':0, 'volume':1})
        volumes_b = [i['volume'] for i in volumes_b]
        num_Born = len(volumes_b)      
        if num_Born>0:
            for i in (self.vasp_db).db['borncharge'].find({'metadata.tag': tag}):
                vol = 'V{:010.6f}'.format(float(i['volume']))
                voldir = phdir+'/'+vol
                if not os.path.exists(voldir):
                   os.mkdir(voldir)

                structure = Structure.from_dict(i['structure'])
                poscar = structure.to(fmt="poscar")
                poscar = str(poscar).split('\n')
                natom = len(structure.sites)
                with open (voldir+'/dielecfij.out','w') as out:
                    for line in range (2,5):
                        out.write('{}\n'.format(poscar[line]))
                    for line in range (8,natom+8):
                        out.write('{}\n'.format(poscar[line]))
                    dielectric_tensor = np.array(i['dielectric_tensor'])
                    for x in range(3):
                        for y in range(3):
                            out.write('{} '.format(dielectric_tensor[x,y]))
                        out.write('\n')
                    born_charge = np.array(i['born_charge'])
                    for ii in range(natom):
                        out.write(' ion  {}\n'.format(ii+1))
                        for x in range(3):
                            out.write('  {} '.format(x+1))
                            for y in range(3):
                                out.write('{} '.format(born_charge[ii,x,y]))
                            out.write('\n')
        return num_Born

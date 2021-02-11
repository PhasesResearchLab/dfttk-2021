import sys
import numpy as np
from scipy.optimize import linprog
from numpy.linalg import solve
from fractions import Fraction
from dfttk.utils import sort_x_by_y
import os
import json
from pymatgen import MPRester, Structure
from atomate.vasp.database import VaspCalcDb

MM_of_Elements = {'H': 1.00794, 'He': 4.002602, 'Li': 6.941, 'Be': 9.012182, 'B': 10.811, 'C': 12.0107, 'N': 14.0067,
              'O': 15.9994, 'F': 18.9984032, 'Ne': 20.1797, 'Na': 22.98976928, 'Mg': 24.305, 'Al': 26.9815386,
              'Si': 28.0855, 'P': 30.973762, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.0983, 'Ca': 40.078,
              'Sc': 44.955912, 'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938045,
              'Fe': 55.845, 'Co': 58.933195, 'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.409, 'Ga': 69.723, 'Ge': 72.64,
              'As': 74.9216, 'Se': 78.96, 'Br': 79.904, 'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.90585,
              'Zr': 91.224, 'Nb': 92.90638, 'Mo': 95.94, 'Tc': 98.9063, 'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.42,
              'Ag': 107.8682, 'Cd': 112.411, 'In': 114.818, 'Sn': 118.71, 'Sb': 121.760, 'Te': 127.6,
              'I': 126.90447, 'Xe': 131.293, 'Cs': 132.9054519, 'Ba': 137.327, 'La': 138.90547, 'Ce': 140.116,
              'Pr': 140.90465, 'Nd': 144.242, 'Pm': 146.9151, 'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25,
              'Tb': 158.92535, 'Dy': 162.5, 'Ho': 164.93032, 'Er': 167.259, 'Tm': 168.93421, 'Yb': 173.04,
              'Lu': 174.967, 'Hf': 178.49, 'Ta': 180.9479, 'W': 183.84, 'Re': 186.207, 'Os': 190.23, 'Ir': 192.217,
              'Pt': 195.084, 'Au': 196.966569, 'Hg': 200.59, 'Tl': 204.3833, 'Pb': 207.2, 'Bi': 208.9804,
              'Po': 208.9824, 'At': 209.9871, 'Rn': 222.0176, 'Fr': 223.0197, 'Ra': 226.0254, 'Ac': 227.0278,
              'Th': 232.03806, 'Pa': 231.03588, 'U': 238.02891, 'Np': 237.0482, 'Pu': 244.0642, 'Am': 243.0614,
              'Cm': 247.0703, 'Bk': 247.0703, 'Cf': 251.0796, 'Es': 252.0829, 'Fm': 257.0951, 'Md': 258.0951,
              'No': 259.1009, 'Lr': 262, 'Rf': 267, 'Db': 268, 'Sg': 271, 'Bh': 270, 'Hs': 269, 'Mt': 278,
              'Ds': 281, 'Rg': 281, 'Cn': 285, 'Nh': 284, 'Fl': 289, 'Mc': 289, 'Lv': 292, 'Ts': 294, 'Og': 294,
              'ZERO': 0}

periodictable = MM_of_Elements.keys() #""" list of all elements from the periodic table"""


"""check if value is a float number"""
def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False


"""convert a chemical formula into element list and composition
formula - chemical formula
return:
element list and composition list
"""
def formula2composition(formula, normalize=True):
  formula = formula.replace(" ",'').replace("-",'').replace(",",'')
  newc = ""
  """Follow the convention, elemental symbol must start from capital letter"""
  for c in formula:
    if c in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
      newc = newc + '|'
    newc = newc + c
  els = newc.split('|')
  els = [k for k in els if k != '']

  """now get the composition for each element"""
  ele = []
  com = []
  for el in els:
    newel = ""
    newcc = ""
    for c in el:
      if c.isalpha():
        newel = newel + c
      else:
        newcc = newcc + c

    if (newel not in periodictable):
      raise ValueError('"'+newel+'" is not an int number! your formula is wrong!')
    ele.append(newel)

    if (len(newcc)!=0):
      if (isfloat(newcc)):
        com.append(float(newcc))
      else:
        raise ValueError('"'+newcc+'" is not an int number! your formula is wrong!')
    else:
      com.append(1.0)
  com = np.array(list(map(float,com)))
  #print("eeeeee", formula, ele, com)
  if normalize: 
      if sum(com)==0.0: raise ValueError("divided by zero")
      com = com/sum(com)

  #sorted the sequence and merge the duplicate
  elist = sorted(set(ele))
  clist = np.zeros(len(elist), dtype=float)
  for j,el in enumerate(ele):
    ix = elist.index(el)
    clist[ix] += com[j]

  return elist,clist


"""make the reduced formula
els - a list of elements
natype - a list of number of elements
return:
"""
def reduced_formula(formula):
  _els, nat = formula2composition(formula.replace(' ',''), False)
  #_els=els.split(' ')
  #_els = [k for k in _els if k != '']
  _nat = np.array(list(map(int,nat)))
  els = sorted(set(_els))
  nat = np.zeros(len(els),dtype=int)
  for i,el in enumerate(_els):
    ix = els.index(el)
    nat[ix] += _nat[i]

  Nd = min(nat)
  for i in range(Nd,0,-1):
    out = True
    for j in range(len(nat)):
      if ((nat[j]//i)*i!=nat[j]):
        out = False
        break
    if out:
      break
  form = ""
  for j,el in enumerate(els):
    ix = nat[j]//i
    form = form+el
    if ix!=1:
      form = form+str(ix)
  #print("eeeeeeee", _els, _nat, formula.replace(' ',''), _nat,form)
  return form

"""exclude almost repeated volume data
vol - single volume to checked
volumes - list of found volumes
thr - threshold
return:
True: vol is in volumes; otherwise False
"""
def vol_within(vol, volumes, thr=0.001):
    for i,v in enumerate(volumes):
        if (abs(vol-v) < thr*vol): return True
    return False


"""get experimental data from a list of dict
expt - a list of dict containing experimental data
formula - chemical formula
return:
matched dict containing experimental data
"""
def get_expt(expt, formula):
    #print ('eeeeeeee', expt)
    if isinstance(expt, str):
        with open(expt, encoding='utf-8') as element_json_file:
            _expt = json.load(element_json_file)
    else: _expt = expt

    f0 = reduced_formula(formula)
    #print ("eeeeeee f0=", f0, formula)
    data = []
    for ee in _expt:
        #print ("eeeeeee", ee, f0)
        if reduced_formula(ee['Compound'])!=f0: continue
        #print ("eeeeeee ee['Compound']", ee['Compound'])
        data.append(ee)
    return data


"""get melting temperature form jANAF table
formula - chemical formula
expt - a dict containing heat capacity data
return:
if found, the highest temperature from the heat capacity data
otherwise, the up temperature limit defined in the argparse args 
"""
def get_melting_temperature_from_JANAF(expt=None, formula=None):

    if expt!=None and formula!=None:
        _expt =get_expt(expt, formula)
        for d in _expt:
            if "Malcolm W. Chase, Jr. NIST-JANAF Thermochemical Tables" in d['Author']:
                if 'heat capacity'==d['property']:
                    return max(d['data'][0::2])
    return None


"""get melting temperature form SSUB database
formula - chemical formula
expt - a dict containing heat capacity data
return:
if found the dict key 'melting', return its value
otherwise, None
"""
def get_melting_temperature(expt=None, formula=None):

    if expt!=None and formula!=None:
        _expt =get_expt(expt, formula)
        for d in _expt:
            #print ('eeeeeeeeee 0', d)
            if "Malcolm W. Chase, Jr. NIST-JANAF Thermochemical Tables" in d['Author']:
                if 'heat capacity'==d['property']:
                    return max(d['data'][0::2])
        for d in _expt:
            #print ('eeeeeeeeee 1', d)
            """
            if "Andersson(CALPHAD), Andersson J.O.," in d['Author']:
                if 'heat capacity'==d['property']:
            """
            try:
                return d['melting']
            except:
                pass
    return None


"""return E-V data information from MonggoDB database of static calculations
vasp_db - MonggoDB database connection
m - metadata tag value
return: E-V, strain, stresses, bandgap etc
"""
def get_rec_from_metatag(vasp_db,m):
    static_calculations = vasp_db.collection.\
        find({'$and':[ {'metadata.tag': m}, {'adopted': True} ]})
    gapfound = False
    energies = []
    volumes = []
    stresses = []
    lattices = []
    bandgaps = []
    pressures = []
    magmoms = []
    emin = 1.e36
    for calc in static_calculations:
        vol = calc['output']['structure']['lattice']['volume']
        if vol_within(vol, volumes): continue
        natoms = len(calc['output']['structure']['sites'])
        try:
            sites = calc['output']['structure']['sites']
            magmoms.append([{s['label']:s['properties']['magmom']} for s in sites])
        except:
            pass
        lat = calc['output']['structure']['lattice']['matrix']
        sts = calc['output']['stress']
        ene = calc['output']['energy']
        if ene < emin:
            structure = Structure.from_dict(calc['input']['structure'])
            POSCAR = structure.to(fmt="poscar")
            INCAR = calc['input']['incar']
        gap = calc['output']['bandgap']
        volumes.append(vol)
        energies.append(ene)
        stresses.append(sts)
        lattices.append(lat)
        bandgaps.append(gap)
        if sts!=None: pressures.append((sts[0][0]+sts[1][1]+sts[2][2])/3.)
        else: pressures.append(None)
        if not gapfound: gapfound = float(gap) > 0.0
    energies = sort_x_by_y(energies, volumes)
    pressures = sort_x_by_y(pressures, volumes)
    stresses = sort_x_by_y(stresses, volumes)
    lattices = sort_x_by_y(lattices, volumes)
    bandgaps = sort_x_by_y(bandgaps, volumes)
    try:
        magmoms = sort_x_by_y(magmoms, volumes)
    except:
        pass
    volumes = sort_x_by_y(volumes, volumes)
    EV = {}
    EV['metatag'] = m
    EV['natoms'] = natoms
    EV['volumes'] = volumes
    EV['stresses'] = stresses
    EV['energies'] = energies
    EV['pressures'] = pressures
    EV['bandgaps'] = bandgaps
    EV['lattices'] = lattices
    try:
        for volume in magmoms:
            for magmom in volume.values():
                if magmom!=0.:
                    EV['magmoms'] = magmoms
                    break
    except:
        pass
    return EV,POSCAR,INCAR

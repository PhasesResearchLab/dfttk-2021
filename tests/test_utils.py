#!python
import pytest

import dfttk.utils as dfttkutils
from dfttk.input_sets import RelaxSet
from pymatgen import Structure, SETTINGS

POSCAR_STR_check_symbol = """FCC_Fe_WithDiffMagMom
1.0
3.5391385555 0.0000000000 0.0000000000
0.0000000000 3.5391385555 0.0000000000
0.0000000000 0.0000000000 3.5391385555
Fe Fe
2 2
Direct
0.000000000 0.000000000 0.000000000
0.000000000 0.500000000 0.500000000
0.500000000 0.500000000 0.000000000
0.500000000 0.000000000 0.500000000"""

POSCAR_scalling_t1 = """Ti1 Pb1 O3
1.0
2.0 0.000000 0.000000
0.000000 2.0 0.000000
0.000000 0.000000 8.0
Ti Pb O
1 1 3
direct
0.500000 0.500000 0.520923 Ti
0.000000 0.000000 0.969212 Pb
0.500000 0.500000 0.142449 O
0.000000 0.500000 0.626658 O
0.500000 0.000000 0.626658 O"""

try:
    API_KEY = SETTINGS["PMG_MAPI_KEY"]
    PMG_VASP_PSP_DIR = SETTINGS["PMG_VASP_PSP_DIR"]
except Exception as e:
    print("Please provide the API_KEY.")
    API_KEY = None
    PMG_VASP_PSP_DIR = None

@pytest.mark.skipif(API_KEY is None, reason="MAPI_KEY required")
def test_mp_structures_from_ids():
    mp_ids = ["mp-66", "mp-22862"]  #66 for Diamond, 22862 for NaCl
    structs = dfttkutils.mp_structures_from_ids(mp_ids, API_KEY=API_KEY)
    assert(structs[0].composition.reduced_formula == "C")
    assert(structs[1].composition.reduced_formula == "NaCl")

@pytest.mark.skipif(API_KEY is None, reason="MAPI_KEY required")
def test_mp_structures_from_system():
    system = "Fe-Cr"
    structs = dfttkutils.mp_structures_from_system(system, API_KEY=API_KEY)
    formula = []
    for s in structs:
        formula.append(s.composition.reduced_formula)
    assert(formula == ['CrFe3', 'Cr2Fe', 'CrFe3', 'Cr3Fe', 'CrFe4', 'CrFe', 'Cr3Fe', 'Cr3Fe'])

@pytest.mark.skipif(API_KEY is None, reason="MAPI_KEY required")
def test_mp_sorted_structures_from_system():
    system = "Fe-Cr"
    sorted_structs = dfttkutils.mp_sorted_structures_from_system(system, API_KEY=API_KEY)
    formula = []
    for s in sorted_structs:
        formula.append(s.composition.reduced_formula)
    assert(formula == ['CrFe4', 'CrFe3', 'CrFe3', 'CrFe', 'Cr3Fe'])

def test_check_symbol():
    struc = Structure.from_str(POSCAR_STR_check_symbol, fmt="POSCAR")
    magmoms = [4.0, 4.0, -4.0, -4.0]
    struc.add_site_property("magmom", magmoms)
    InputSet = RelaxSet(struc)
    symbols, natom = dfttkutils.check_symbol(InputSet)
    assert(symbols == ["Fe", "Fe"])
    assert(natom == ["2", "2"])

def test_update_pos_by_symbols():
    struc = Structure.from_str(POSCAR_STR_check_symbol, fmt="POSCAR")
    magmoms = [4.0, 4.0, -4.0, -4.0]
    struc.add_site_property("magmom", magmoms)
    InputSet = RelaxSet(struc)
    poscar_str = dfttkutils.update_pos_by_symbols(InputSet, write_file=False)
    syms = poscar_str.split("\n")[5]
    natom = poscar_str.split("\n")[6]
    assert(syms == "Fe Fe")
    assert(natom == "2 2")

@pytest.mark.skipif(PMG_VASP_PSP_DIR is None, reason="PMG_VASP_PSP_DIR required")
def test_update_pot_by_symbols():
    struc = Structure.from_str(POSCAR_STR_check_symbol, fmt="POSCAR")
    magmoms = [4.0, 4.0, -4.0, -4.0]
    struc.add_site_property("magmom", magmoms)
    InputSet = RelaxSet(struc)
    potcar = dfttkutils.update_pot_by_symbols(InputSet, write_file=False)
    syms = potcar.symbols
    assert(syms == ["Fe", "Fe"])

def test_supercell_scaling_by_atom_lat_vol():
    min_atoms = 60
    max_atoms = 90
    lower_search_limit = -2
    upper_search_limit = 2
    target_shape = 'sc'
    #test for cubic
    stru1 = Structure.from_str(POSCAR_STR_check_symbol, fmt="POSCAR")
    #optimal_sc_shape = dfttkutils.supercell_scaling_by_target_atoms(stru1, min_atoms=min_atoms, max_atoms=max_atoms,
    #                                  target_shape=target_shape, lower_search_limit=lower_search_limit,
    #                                  upper_search_limit=upper_search_limit, verbose=False)
    stru2 = Structure.from_str(POSCAR_scalling_t1, fmt='POSCAR')
    optimal_sc_shape = dfttkutils.supercell_scaling_by_atom_lat_vol(stru2, min_obj=min_atoms, max_obj=max_atoms,
                                      target_shape=target_shape, lower_search_limit=lower_search_limit,
                                      upper_search_limit=upper_search_limit, verbose=False)
    print(optimal_sc_shape)

#test_supercell_scaling_by_target_atoms()

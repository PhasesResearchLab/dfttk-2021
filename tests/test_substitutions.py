#!python
#

import pytest
#import sys
#import os
#dfttkhome = os.path.abspath(os.path.join('..'))
#sys.path.append(dfttkhome)
import dfttk.structure_builders.substitutions as substitutions
import dfttk.utils
from pymatgen import Structure


POSCAR_STR = """FeNi3
1.0
3.5391385555 0.0000000000 0.0000000000
0.0000000000 3.5391385555 0.0000000000
0.0000000000 0.0000000000 3.5391385555
Fe Ni
1 3
Direct
0.000000000 0.000000000 0.000000000
0.000000000 0.500000000 0.500000000
0.500000000 0.500000000 0.000000000
0.500000000 0.000000000 0.500000000"""
struct = Structure.from_str(POSCAR_STR, fmt='POSCAR')


def test_canonicalize_config():
    configuration = [['Ni', 'Fe'], ['Fe']]
    occupancies = [[0.75, 0.25], [1.0]]
    (new_configuration, new_occupancies) = substitutions.canonicalize_config(configuration, occupancies)
    assert new_configuration == [['Fe', 'Ni'], ['Fe']]
    assert new_occupancies == [[0.25, 0.75], [1.0]]

def test_get_ele_list_from_struct():
    ele_list = substitutions.get_ele_list_from_struct(struct)
    assert ele_list == ['Fe', 'Ni', 'Ni', 'Ni']

def test_get_density_from_pt():
    den_dict = {'Nb': 8.57, 'Ti': 4.507}
    #test for the list of elements
    ele_list = ['Nb', 'Ti']
    density_dict = substitutions.get_density_from_pt(ele_list)
    for ele in ele_list:
        assert den_dict[ele] == density_dict[ele]
    #test for the dict of tlements
    ele_dict = {'Nb': 3, 'Ti': 1}
    density_dict = substitutions.get_density_from_pt(ele_dict)
    for ele in ele_dict:
        assert den_dict[ele] == density_dict[ele]

def test_scale_struct():
    struct_new = substitutions.scale_struct(struct)
    assert struct_new.lattice.a == 3.5443397446212437

def test_gen_replacement_dict():
    old_config = [['Fe', 'Cr'], ['Ni']]
    new_config = [['Ti', 'V'], ['Zr']]
    replacement_dict = substitutions.gen_replacement_dict(old_config, new_config)
    assert replacement_dict['Fe'] == 'Ti'
    assert replacement_dict['Cr'] == 'V'
    assert replacement_dict['Ni'] == 'Zr'

def test_substitute_configuration():
    template_config = [['Fe', 'Ni']]
    config = [['Al', 'Ti']]
    struct_new = substitutions.substitute_configuration(struct, template_config, config, check_sorting=True)
    assert struct_new.lattice.a == 4.11833824135907
    new_cofig = set()
    for ele in struct_new.species:
        new_cofig.add(str(ele))
    #Be cautions: this assert is right in current state, it may not true in partial substitute
    assert dfttk.utils.recursive_flatten(config) == sorted(list(new_cofig))

def test_substitute_configuration_with_metadata():
    template_config = [['Fe', 'Ni']]
    config = [['Al', 'Ti']]
    occupation = [[0.5, 0.5]]
    phase_name = "FCC_A1"
    site_ratios = [1.0]
    struct_new, metadata = substitutions.substitute_configuration_with_metadata(struct, template_config, config, occupation, phase_name, site_ratios)
    assert struct_new.lattice.a == 4.11833824135907
    assert metadata["phase_name"] == phase_name

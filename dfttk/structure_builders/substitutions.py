"""
Tools for substituting structures and generating metadata
"""

from copy import deepcopy
from dfttk.utils import sort_x_by_y

def canonicalize_config(configuration, occupancies):
    """
    Return canonicalized (sorted) configurations and occupancies.

    Parameters
    ----------
    configuration : list of lists
        DFTTK-style configuration
    occupancies :
        DFFTK-style occupancies

    Returns
    -------
    tuple
        Tuple of canonical (configuration, occupancies)

    Notes
    -----
    The sublattice ordering is preserved, but the species within a sublattice are sorted.

    Examples
    --------
    >>> canonicalize_config([['Fe', 'Ni'], ['Fe']], [[0.25, 0.75], [1.0]])  # already canonical
    ([['Fe', 'Ni'], ['Fe']], [[0.25, 0.75], [1.0]])
    >>> canonicalize_config([['Cu'], ['Mg']], [[1.0], [1.0]])  # already canonical
    ([['Cu'], ['Mg']], [[1.0], [1.0]])
    >>> canonicalize_config([['Cu', 'Mg']], [[0.9, 0.1]])  # already canonical
    ([['Cu', 'Mg']], [[0.9, 0.1]])
    >>> canonicalize_config([['Cu', 'Mg']], [[0.1, 0.9]])  # already canonical
    ([['Cu', 'Mg']], [[0.1, 0.9]])
    >>> canonicalize_config([['Ni', 'Fe'], ['Fe']], [[0.75, 0.25], [1.0]])
    ([['Fe', 'Ni'], ['Fe']], [[0.25, 0.75], [1.0]])
    >>> canonicalize_config([['Ni', 'Fe'], ['Fe', 'Cr', 'Ni']], [[0.75, 0.25], [0.1, 0.2, 0.7]])
    ([['Fe', 'Ni'], ['Cr', 'Fe', 'Ni']], [[0.25, 0.75], [0.2, 0.1, 0.7]])

    """
    new_occupancies = [sort_x_by_y(occ, config) for occ, config in zip(occupancies, configuration)]
    new_configuration = [sorted(config) for config in configuration]
    return (new_configuration, new_occupancies)

def get_density_from_pt(ele_list):
    """
    Get density(g/cm^3) from periodictable package

    Parameters
    ----------
        ele_list : list/dict
            The list of elements, e.g. ['Nb', 'Ti']/{'Nb': 3, 'Ti': 1}
    Returns
    -------
        density_dict : dict
            Dictionary of {element: density}, e.g. {'Nb': 8.57, 'Ti': 4.507}. 
    Examples
    --------
    >>> get_density_from_pt(['Nb', 'Ti'])
    {'Nb': 8.57, 'Ti': 4.507}
    """
    from pymatgen.core.periodic_table import Element
    #import periodictable as pt
    density_dict = {}
    for ele in ele_list:
        density_dict[ele] = float(Element(ele).density_of_solid)/1000.
    return density_dict

def get_ele_list_from_struct(struct):
    """
    Get elements list from pymatgen structure objective

    Parameters
    ----------
        struct : pymatgen objective
            The structure
    Returns
    -------
        ele_list : [str]
            The list of elements
    """
    ele_list = []
    for ele in struct.species:
        ele_list.append(str(ele))
    return ele_list

def scale_struct(struct):
    """Scale the structure according to the weighted average density of each element.

    Parameters
    ----------
    struct : pymatgen.Structure
    #density_dict : dict
    #    Dictionary of {element: density}, e.g. {'Fe': 9, 'Ti': 4}. The units do
    #    not matter as long as the densities in the dict are internally consistent.
    # This parameters is canceled by using periodictable module in pymatgen.core

    Returns
    -------
    pymatgen.Structure
        Modifies the structure in place, but also returns for convenience.

    """
    species_amnt_dict = struct.composition.get_el_amt_dict()  # dict of {'V': 10.0, 'Ni': 30.0}
    density_dict = get_density_from_pt(species_amnt_dict)
    # densities is dict of densities, {'V': 6.313, 'Ni': 9.03}
    expected_density = float(sum([density_dict[species]*amnt for species, amnt in species_amnt_dict.items()]))/struct.composition.num_atoms
    current_density = struct.density
    current_volume = struct.volume
    expected_volume = current_volume/expected_density*current_density
    struct.scale_lattice(float(expected_volume))
    return struct

def gen_replacement_dict(old_config, new_config):
    """Create a pymatgen replacement dict based on old and new sublattice configurations.

    Shapes of the config lists must match.

    Parameters
    ----------
    old_config : list
        DFTTK style configuration that will be replaced
    new_config : list
        DFTTK style configuration that the new structure will have

    Returns
    -------
    dict
        Dict of {element_to_replace: new_element}
    """
    replacement_dict = {}
    for new_subl, old_subl in zip(new_config, old_config):
        for new_atom, old_atom in zip(new_subl, old_subl):
            replacement_dict[old_atom] = new_atom
    return replacement_dict


def substitute_configuration(template_structure, template_config, config, check_sorting=True):
    """
    Replace the species in the template structure by switching the template_config elements for the config elements.

    Scales the structure according to the weighted element densities. Shapes of
    the config lists must match.

    Parameters
    ----------
    template_structure : pymatgen.Structure
        Structure with species that will be replaced
    template_config : list
        DFTTK style configuration that will be replaced
    config : list
        DFTTK style configuration that the new structure will have
    #density_dict : dict
    #    Dictionary of {element: density}, e.g. {'Fe': 9, 'Ti': 4}. The units do
    #    not matter as long as the densities in the dict are internally consistent.
    check_sorting : bool
        If True, will check that all sublattices are correctly sorted in the target config

    Returns
    -------
    pymatgen.Structure
        A new Structure object (the original is not modified so it can be reused in loops).
    """
    if check_sorting:
        for subl in config:
            if subl != list(sorted(subl)):
                raise ValueError("Configuration {} is not in sorted order. "
                                 "See information on DFTTK configurations in the docs.".format(config))
    struct = deepcopy(template_structure)
    struct.replace_species(gen_replacement_dict(template_config, config))
    scale_struct(struct)
    return struct



def substitute_configuration_with_metadata(template_structure, template_config, config, occupation, phase_name, site_ratios):
    """
    Replace the species in the template structure by switching the template_config elements for the config elements.

    Scales the structure according to the weighted element densities. Shapes of
    the config lists must match.

    Wrapper around subsitute_configuration that returns a tuple of the structure and metadata.

    Parameters
    ----------
    template_structure : pymatgen.Structure
        Structure with species that will be replaced
    template_config : list
        DFTTK style configuration that will be replaced
    config : list
        DFTTK style configuration that the new structure will have
    #density_dict : dict
    #    Dictionary of {element: density}, e.g. {'Fe': 9, 'Ti': 4}. The units do
    #    not matter as long as the densities in the dict are internally consistent.
    occupation : list
        DFTTK style occupancy fractions. Must match the shape to config and template_config.
    phase_name : str
        Name of the phase
    site_ratios : list
        Sublattice site ratios, 1-d list.

    Returns
    -------
    (pymatgen.Structure, dict)
        Tuple of a new Structure object (the original is not modified so it can be reused in loops) and a dict of metadata
    """
    struct = substitute_configuration(template_structure, template_config, config, check_sorting=False)
    config, occupation = canonicalize_config(config, occupation)
    metadata = {'phase_name': phase_name, 'sublattice': {'configuration': config, 'occupancies': occupation, 'site_ratios': site_ratios}}
    return struct, metadata

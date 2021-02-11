from dfttk.utils import recursive_flatten

def to_element_case(el):
    """Convert an uppercase elemnt to element case, e.g. FE to Fe, V to V."""
    return el[0].upper() + el[1:].lower()


def dfttk_config_to_espei(config):
    """Convert a DFTTK configuration, e.g. [['Fe', 'Ni'], ['Fe']] to ESPEI's configuration [['FE', 'NI'], 'FE']"""
    espei_config = []
    for subl in config:
        if len(subl) == 1:
            # this sublattice is an endmember
            espei_config.append(subl[0].upper())
        else:
            # this sublattice has an interaction
            espei_config.append([comp.upper() for comp in subl])
    return espei_config


def dfttk_occupancies_to_espei(occupancies):
    """Convert DFTTK occupancies, e.g. [[0.5 0.5], [1.0]] to ESPEI's configuration [[0.5, 0.5], 1.0]"""
    espei_occupancies = []
    for subl in occupancies:
        if len(subl) == 1:
            # this sublattice is an endmember
            espei_occupancies.append(1.0)
        else:
            # this sublattice has an interaction
            espei_occupancies.append([occ for occ in subl])
    return espei_occupancies


def espei_config_to_dfttk(config):
    """Convert a ESPEI configuration, e.g. [['FE', 'NI'], 'FE'] to DFTTK's configuration [['Fe', 'Ni'], ['Fe']]"""
    dfttk_config = []
    for subl in config:
        if isinstance(subl, str):
            # this sublattice is an endmember, e.g. 'FE'
            dfttk_config.append([to_element_case(subl)])
        else:
            # this sublattice is an interacting sublattice, e.g. ['FE', 'NI']
            dfttk_config.append([to_element_case(comp) for comp in subl])
    return dfttk_config


def espei_occupancies_to_dfttk(occupancies):
    """Convert a ESPEI configuration, e.g. [[0.5, 0.5], 1.0] to DFTTK's configuration [[0.5, 0.5], [1.0]]"""

    dfttk_occupancies = []
    for subl in occupancies:
        if isinstance(subl, float):
            # this sublattice is an endmember, e.g. 1.0 occupancy
            dfttk_occupancies.append([1.0])
        else:
            # this sublattice is an interacting sublattice, e.g. [0.5, 0.5]
            dfttk_occupancies.append([f for f in subl])
    return dfttk_occupancies


def make_dataset(phase_name, prop, subl_site_ratios, configurations,
                        conditions, values, occupancies=None, tag=None):
    """

    Parameters
    ----------
    phase_name : str
        Name of the phase, e.g. 'FCC_A1'
    prop : str
        Name of the property, e.g. 'HM_FORM', 'SM_MIX'
    tag : str
        UUID tag as a string, e.g. '389ng-n3nal3ih-andjladf'
    subl_site_ratios : list
        List of sublattice site ratios corresponding to the sublattice model, e.g. [1.0, 1.0, 3.0].
    configurations : list
        List of ESPEI-style sublattice configurations.
    conditions : dict
        Dictionary of condition lists. Conditions should be scalars or 1-d arrays (not broadcasted).
    values : numpy.array
        Multidimensional (broadcasted) array of values. Should conform to
        ESPEI datasets shape of (P, T, configurations).
    occupancies : optional, list
        List of ESPEI-style sublattice occupancies. Only required for
        configurations that have interactions (not endmembers).

    Returns
    -------
    dict
        ESPEI dataset dict

    Notes
    -----
    See description of ESPEI's single phase datasets for what the outputs should
    be: http://espei.org/en/latest/input_data.html#single-phase-data
    """
    components = sorted(set(recursive_flatten(configurations)))

    dataset_dict = {}
    dataset_dict['solver'] = {}
    dataset_dict['output'] = prop
    dataset_dict['components'] = components
    dataset_dict['phases'] = [phase_name]
    dataset_dict['conditions'] = conditions
    dataset_dict['values'] = values.tolist()
    dataset_dict['conditions'] = conditions
    dataset_dict['solver']['sublattice_site_ratios'] = subl_site_ratios
    dataset_dict['solver']['sublattice_configurations'] = configurations
    dataset_dict['solver']['mode'] = 'manual'
    dataset_dict['reference'] = "DFTTK calculation"
    if occupancies is not None:
        dataset_dict['solver']['sublattice_occupancies'] = occupancies
    if tag is not None:
        dataset_dict['comment'] = "DFTTK tag: >>{}<<".format(tag)
    return dataset_dict


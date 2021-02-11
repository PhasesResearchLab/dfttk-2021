"""
The sqs module handles converting abstract SQS Structure objects to concrete structures.

The SQS are regular pymatgen Structures with the species named according to sublattice and species type.
These species in pymatgen Structures are named to `Xab`, which corresponds to atom `B` in sublattice `a`.
"""

from __future__ import division

import copy
import itertools

import pymatgen as pmg
from pymatgen import Structure

from .prl_structure import PRLStructure


class AbstractSQS(Structure):
    """A pymatgen Structure with special features for SQS.
    """

    def __init__(self, *args, **kwargs):
        """Create a SQS object

        Parameters
        ----------
        args :
            args to pass to Structure
        sublattice_model : [[str]]
            Abstract sublattice model in the ESPEI style, e.g. `[['a', 'b'], ['a']]`.
        sublattice_names : [[str]]
            Names of the sublattices, or the second character in the species names, e.g. `['a', 'c']`.
        kwargs :
            kwargs to pass to Structure
        """
        self.sublattice_model = kwargs.pop('sublattice_model', None)
        self._sublattice_names = kwargs.pop('sublattice_names', None)
        super(AbstractSQS, self).__init__(*args, **kwargs)

    @property
    def normalized_sublattice_site_ratios(self):
        """Return normalized sublattice site ratio. E.g. [[0.25, 0.25], [0.1666, 0.1666, 0.1666]]
        """
        subl_model = self.sublattice_model
        subl_names = self._sublattice_names
        comp_dict = self.composition.as_dict()
        site_ratios = [[comp_dict['X'+name+e+'0+']/self.num_sites for e in subl] for subl, name in zip(subl_model, subl_names)]
        return site_ratios

    @property
    def sublattice_site_ratios(self):
        """Return normalized sublattice site ratio. E.g. [[0.25, 0.25], [0.1666, 0.1666, 0.1666]]
        """
        subl_model = self.sublattice_model
        subl_names = self._sublattice_names
        comp_dict = {k: int(v) for k, v in self.composition.reduced_composition.as_dict().items()}
        site_ratios = [[comp_dict['X'+name+e+'0+'] for e in subl] for subl, name in zip(subl_model, subl_names)]
        return site_ratios

    def get_concrete_sqs(self, subl_model, scale_volume=True):
        """Modify self to be a concrete SQS based on the sublattice model.

        Parameters
        ----------
        subl_model : [[str]]
            List of strings of species names. Must exactly match the shape of self.sublattice_model.
            **Note that order does matter!** [["Al", "Fe"]] and [["Fe", "Al"]] will produce different results!
        scale_volume : bool
            If True, scales the volume of the cell so the ions have at least their minimum atomic radii between them.
        """
        def _subl_error():
            raise ValueError('Concrete sublattice model {} does not match size of abstract sublattice model {}'.format(subl_model, self.sublattice_model))
        if len(subl_model) != len(self.sublattice_model):
            _subl_error()

        # build the replacement dictionary and the site ratios
        # we have to look up the sublattice names to build the replacement species names
        replacement_dict = {}
        site_occupancies = [] # list of [{'FE': 0.3333, 'NI': 0.6666}, {'FE': 1}] for [['FE', 'NI'], ['FE]]
        for abstract_subl, concrete_subl, subl_name, subl_ratios in zip(self.sublattice_model, subl_model, self._sublattice_names, self.sublattice_site_ratios):
            if len(abstract_subl) != len(concrete_subl):
                _subl_error()
            sublattice_ratio_sum = sum(subl_ratios)
            sublattice_occupancy_dict = {}
            for abstract_specie, concrete_specie, site_ratio in zip(abstract_subl, concrete_subl, subl_ratios):
                specie = 'X' + subl_name + abstract_specie
                replacement_dict[specie] = concrete_specie
                sublattice_occupancy_dict[concrete_specie] = sublattice_occupancy_dict.get(concrete_specie, 0) + site_ratio/sublattice_ratio_sum
            site_occupancies.append(sublattice_occupancy_dict)

        # create a copy of myself to make the transformations and make them
        self_copy = copy.deepcopy(self)
        self_copy.replace_species(replacement_dict)

        if scale_volume:
            fractional_comp = dict(self_copy.composition.fractional_composition)
            estimated_density = 0
            for component in self_copy.composition.elements :
                temp = pmg.Element(component).data['Density of solid']
                density = float(temp.split(' ')[0])
                estimated_density += (fractional_comp[component] * density)/1000
            self_copy.scale_lattice(float((self_copy.volume/estimated_density)*self_copy.density))

        # finally we will construct the SQS object and set the values for the canonicalized
        # sublattice configuration, site ratios, and site occupancies

        # first, canonicalize the sublattice model, e.g. [['FE', 'FE'], ['NI']] => [['FE'], ['NI']]
        sublattice_configuration = [sorted(set(subl)) for subl in subl_model]
        # construct the sublattice occupancies for the model
        sublattice_occupancies = [[occupancies[specie] for specie in subl] for occupancies, subl in zip(site_occupancies, sublattice_configuration)]
        # sum up the individual sublattice site ratios to the total sublattice ratios.
        # e.g [[0.25, 0.25], [0.1666, 0.1666, 0.1666]] => [0.5, 0.5]
        site_ratios = [sum(ratios) for ratios in self.sublattice_site_ratios]

        # create the SQS and add all of these properties to our SQS
        concrete_sqs = PRLStructure.from_sites(self_copy.sites)
        concrete_sqs.sublattice_configuration = sublattice_configuration
        concrete_sqs.sublattice_occupancies = sublattice_occupancies
        concrete_sqs.sublattice_site_ratios = site_ratios
        return concrete_sqs


    def get_endmember_space_group_info(self, symprec=1e-2, angle_tolerance=5.0):
        """
        Return endmember space group info..

        Args:
            symprec (float): Same definition as in SpacegroupAnalyzer.
                Defaults to 1e-2.
            angle_tolerance (float): Same definition as in SpacegroupAnalyzer.
                Defaults to 5 degrees.

        Returns:
            spacegroup_symbol, international_number
        """
        endmember_subl = [['X' + subl_name for _ in subl] for subl, subl_name in
                          zip(self.sublattice_model, self._sublattice_names)]
        # we need to replace the abstract names with real names of species.
        endmember_speices = {specie for subl in endmember_subl for specie in subl}
        real_species_dict = {abstract_specie: real_specie for abstract_specie, real_specie in
                             zip(endmember_speices, pmg.core.periodic_table._pt_data.keys())}
        # replace them
        endmember_subl = [[real_species_dict[specie] for specie in subl] for subl in endmember_subl]
        # get the structure and spacegroup info
        endmember_struct = self.get_concrete_sqs(endmember_subl, scale_volume=False)
        endmember_space_group_info = endmember_struct.get_space_group_info(symprec=symprec, angle_tolerance=angle_tolerance)
        return endmember_space_group_info

    def as_dict(self, verbosity=1, fmt=None, **kwargs):
        d = super(AbstractSQS, self).as_dict(verbosity=verbosity, fmt=fmt, **kwargs)
        d['sublattice_model'] = self.sublattice_model
        d['sublattice_names'] = self._sublattice_names
        d['sublattice_site_ratios'] = self.sublattice_site_ratios
        endmember_symmetry = self.get_endmember_space_group_info()
        d['symmetry'] = {'symbol': endmember_symmetry[0], 'number': endmember_symmetry[1]}
        return d

    @classmethod
    def from_dict(cls, d, fmt=None):
        sqs = super(AbstractSQS, cls).from_dict(d, fmt=fmt)
        sqs.sublattice_model = d.get('sublattice_model')
        sqs._sublattice_names = d.get('sublattice_names')
        return sqs


def enumerate_sqs(structure, subl_model, scale_volume=True, skip_on_failure=False):
    """
    Return a list of all of the concrete Structure objects from an abstract Structure and concrete sublattice model.
    Parameters
    ----------
    structure : AbstractSQS
        SQS object. Must be abstract.
    subl_model : [[str]]
        List of strings of species names, in the style of ESPEI `input.json`. This sublattice model
        can be of higher dimension than the SQS, e.g. a [["Al", "Fe", "Ni"]] for a fcc 75/25 binary SQS
        will generate the following Structures:
        Al0.75Fe0.25, Al0.75Ni0.25      Fe0.75Al0.25, Fe0.75Ni0.25      Ni0.75Al0.25, Ni0.75Fe0.25
        *Note that the ordering of species the sublattice model does not matter!*
    scale_volume : bool
        If True, scales the volume of the cell so the ions have at least their minimum atomic radii between them.
    skip_on_failure : bool
        If True, will skip if the sublattice model is lower order and return [] instead of raising

    Returns
    -------
    [PRLStructure]
        List of all concrete PRLStructure objects that can be created from the sublattice model.
    """
    if len(subl_model) != len(structure.sublattice_model):
        raise ValueError('Passed sublattice model ({}) does not agree with the passed structure ({})'.format(subl_model, structure.sublattice_model))
    possible_subls = []
    for subl, abstract_subl in zip(subl_model, structure.sublattice_model):
        subls = itertools.product(subl, repeat=len(abstract_subl))
        possible_subls.append(subls)
    unique_subl_models = itertools.product(*possible_subls)

    # create a list of unique concrete structures with the generated sublattice models
    unique_sqs = []
    unique_configurations_occupancies = []
    for model in unique_subl_models:
        proposed_sqs = structure.get_concrete_sqs(model, scale_volume)
        proposed_config_occupancy = (proposed_sqs.sublattice_configuration, proposed_sqs.sublattice_occupancies)
        if proposed_config_occupancy not in unique_configurations_occupancies:
            unique_configurations_occupancies.append(proposed_config_occupancy)
            unique_sqs.append(proposed_sqs)
    return unique_sqs

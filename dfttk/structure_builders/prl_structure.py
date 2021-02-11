from pymatgen import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

class PRLStructure(Structure):
    """A pymatgen Structure object, with some customizations for ESPEI.
    """

    def __init__(self, *args, **kwargs):
        """Create a Structure object, with some customizations for ESPEI

        Parameters
        ----------
        args :
            args to pass to Structure
        sublattice_configuration : [[str]]
            Sublattice configuration  e.g. `[['Fe', 'Ni'], ['Fe']]`.
        sublattice_occupancies : [[float]]
            Fraction of the sublattice each element in the configuration  has e.g. `[[0.3333, 0.6666], [1]]`.
        sublattice_site_ratios : [float]
            Ratios of sublattice multiplicity  e.g. `[3, 1]`.
        kwargs :
            kwargs to pass to Structure
        """
        self.sublattice_configuration = kwargs.pop('sublattice_configuration', None)
        self.sublattice_occupancies = kwargs.pop('sublattice_occupancies', None)
        self.sublattice_site_ratios = kwargs.pop('sublattice_site_ratios', None)
        super(PRLStructure, self).__init__(*args, **kwargs)

    def __eq__(self, other):
        """
        self and other are equivalent if the sublattice models are equal

        Parameters
        ----------
        other : PRLStructure
        """
        if not isinstance(other, PRLStructure):
            return False
        subl_config = self.sublattice_configuration == other.sublattice_configuration
        subl_site_ratios = self.sublattice_site_ratios == other.sublattice_site_ratios
        subl_occupancies = self.sublattice_occupancies == other.sublattice_occupancies
        return subl_config and subl_site_ratios and subl_occupancies

    @property
    def espei_sublattice_configuration(self):
        """
        Return ESPEI-formatted sublattice model [['a', 'b'], 'a'] for the concrete case
        """
        # short function to convert [['A', 'B'], ['A']] to [['A', 'B'], 'A'] as in ESPEI format
        canonicalize_sublattice = lambda sl: sl[0] if len(sl) == 1 else sl
        return [canonicalize_sublattice(sl) for sl in self.sublattice_configuration]

    @property
    def espei_sublattice_occupancies(self):
        """
        Return ESPEI-formatted sublattice occupancies [[0.3333, 0.6666], 1] for the concrete case
        """
        # short function to convert [[0.3333, 0.6666], [1]] to [[0.3333, 0.6666], 1] as in ESPEI format
        canonicalize_sublattice = lambda sl: sl[0] if len(sl) == 1 else sl
        return [canonicalize_sublattice(sl) for sl in self.sublattice_occupancies]

    def as_dict(self, verbosity=1, fmt=None, **kwargs):
        d = super(PRLStructure, self).as_dict(verbosity=verbosity, fmt=fmt, **kwargs)
        d['sublattice_configuration'] = self.sublattice_configuration
        d['sublattice_occupancies'] = self.sublattice_occupancies
        d['sublattice_site_ratios'] = self.sublattice_site_ratios
        return d

    @classmethod
    def from_dict(cls, d, fmt=None):
        struct = super(PRLStructure, cls).from_dict(d, fmt=fmt)
        struct.sublattice_configuration = d.get('sublattice_configuration')
        struct.sublattice_occupancies = d.get('sublattice_occupancies')
        struct.sublattice_site_ratios = d.get('sublattice_site_ratios')
        return struct


    @classmethod
    def from_structure(cls, structure, equivalent_wyckoff_sites=None):
        """

        Parameters
        ----------
        structure : pymatgen.Structure
        equivalent_wyckoff_sites : list of lists
            List of Wyckoff sites that are treated as the same sublattice, e.g. [['b', 'f']] will
            give combine Wyckoff site 'b' and Wyckoff site 'f' into one sublattice. Putting the same
            Wyckoff site in multiple equivalent groups will produce undefined results.

        Returns
        -------
        PRLStructure
        """
        struct = PRLStructure.from_dict(structure.as_dict())
        # normalize the input structure to a pure element to get Wyckoff sites
        structure = Structure.from_dict(structure.as_dict())
        structure.replace_species({sp.name: "H" for sp in structure.species})
        sga = SpacegroupAnalyzer(structure)
        wyckoff_sites = sga.get_symmetry_dataset()['wyckoffs']
        true_sublattices = sorted(set(wyckoff_sites))
        if equivalent_wyckoff_sites is not None:
            # transform the true sublattices by combining equivalent sites
            combined_sublattices = [''.join(sorted(sites)) for sites in equivalent_wyckoff_sites]
            def match_subl(candidate):
                for subl in combined_sublattices:
                    # if the candidate site is in the combined sublattice, return the combined sublattice
                    if candidate in subl:
                        return subl
                # no match found
                return candidate

            new_subl_model = sorted(set([match_subl(subl) for subl in true_sublattices]))
        else:
            new_subl_model = true_sublattices

        #ratios = [sum([1 if site in subl else 0 for site in wyckoff_sites]) for subl in new_subl_model]
        config = []
        occ = []
        ratios = []
        for subl in new_subl_model:
            species_frequency_dict = {}
            for site, wyckoff_site in zip(struct.sites, wyckoff_sites):
                if wyckoff_site in subl:
                    species = site.specie.name
                    species_frequency_dict[species] = species_frequency_dict.get(species, 0) + 1
            total_subl_occupation = sum(species_frequency_dict.values())
            subl_species = sorted(set(species_frequency_dict.keys()))
            subl_occpancy = [species_frequency_dict[sp]/total_subl_occupation for sp in subl_species]
            config.append(subl_species)
            occ.append(subl_occpancy)
            ratios.append(total_subl_occupation)
        #config = [sorted(set([site.specie.name for site, wyckoff in  if wyckoff in subl])) for subl in new_subl_model]

        struct.sublattice_configuration = config
        struct.sublattice_occupancies = occ
        struct.sublattice_site_ratios = ratios
        return struct

    @staticmethod
    def reindex_sublattice(new_indices, subl_model, subl_occupancies, subl_site_ratios):
        """
        Re-index the passed sublattice model, occupancies and site ratios according to the new index.

        Parameters
        ----------
        new_indices : [int]
            List of indicies corresponding to sublattices. There should be no duplicates. Specifically,
            sorted(new_indices) == list(range(len(subl_model))
        subl_model : [[str]]
            Sublattice configuration  e.g. `[['Fe', 'Ni'], ['Fe']]`.
        subl_occupancies : [[float]]
            Fraction of the sublattice each element in the configuration  has e.g. `[[0.3333, 0.6666], [1]]`.
        subl_site_ratios : [float]
            Ratios of sublattice multiplicity  e.g. `[3, 1]`.

        Returns
        -------
        tuple
            Tuple of (sublattice model, occupancies, site ratios) that have been re-indexed

        Examples
        --------

        >>> PRLStructure.reindex_sublattice([1, 0], [['Al', 'Ni'], ['Al']], [[0.333, 0.666], [1]], [3, 1])
        ([['Al'], ['Al', 'Ni']], [[1], [0.333, 0.666]], [1, 3])
        """
        if sorted(new_indices) != list(range(len(subl_model))):
            raise ValueError('Passed re-indexing indicies ({}) do not match the sublattice model indices ({}).'.format(new_indices, list(range(len(subl_model)))))
        new_subl_model = [subl_model[i] for i in new_indices]
        new_subl_occupancies = [subl_occupancies[i] for i in new_indices]
        new_subl_site_ratios = [subl_site_ratios[i] for i in new_indices]
        return (new_subl_model, new_subl_occupancies, new_subl_site_ratios)


    def reindex(self, new_indices):
        """
        Re-index the instance sublattice model, occupancies and site ratios according to the new index.

        Parameters
        ----------
        new_indices : [int]
            List of indicies corresponding to sublattices. There should be no duplicates. Specifically,
            sorted(new_indices) == list(range(len(subl_model))

        """
        self.sublattice_configuration, self.sublattice_occupancies, self.sublattice_site_ratios = \
            PRLStructure.reindex_sublattice(new_indices, self.sublattice_configuration, self.sublattice_occupancies, self.sublattice_site_ratios)

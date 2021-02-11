import numpy as np
from pymatgen import Structure
from dfttk.utils import eV_per_atom_to_J_per_mol, mget

def get_thermal_props(qha_result, phonon=True):
    """
    Return a dictionary of thermal properties in J/mol-atom from a QHA analysis.

    Parameters
    ----------
    qha_result : Dictionary of a QHA summary dict
    phonon : bool
        Whether to get the phonon data (True) or Debye data (False). Defaults to True.
        If the has_phonon key is False, it will fall back to Debye automatically.

    Returns
    -------
    dict
        Dictionary of thermal properties. Dictionary keys are GM, HM, SM, CPM, and T.
    """
    struct = Structure.from_dict(qha_result['structure'])
    if phonon and qha_result['has_phonon']:
        G = np.array(mget(qha_result, 'phonon.gibbs_free_energy')) * eV_per_atom_to_J_per_mol / struct.composition.num_atoms
        T = np.array(mget(qha_result, 'phonon.temperatures'))
    else:
        G = np.array(mget(qha_result, 'debye.gibbs_free_energy')) * eV_per_atom_to_J_per_mol / struct.composition.num_atoms
        T = np.array(mget(qha_result, 'debye.temperatures'))
    dT = T[1] - T[0]
    Cp = -T * np.gradient(np.gradient(G, dT), dT)
    S = - np.gradient(G, dT)
    H = G + T * S

    thermal_properties = {'GM': G, 'HM': H, 'SM': S, 'CPM': Cp, 'T': T}
    return thermal_properties


def get_formation_energy(tprops, struct, refstate_dict, prop, idx=None, thin=None):
    """

    Parameters
    ----------
    tprops : dict
        Dictionary of thermal properties containing a key of `prop` and a NumPy array as values
    struct : pymatgen.Structure
        pymatgen Structure object
    refstate_dict : dict
        Dictionary of thermal properties for each element chosen to be the reference.
        The dict should have keys of the element names, e.g. {'Cr', 'Ni'} and the
        values should be thermal properties dictionaries equivalent to the `tprops`
        argument above.
    prop : str
        String of the thermal property name, e.g. 'HM', 'SM', 'CPM'
    idx : optional, int
        Integer of the specific index to pull out from the thermal properties.
        E.g. the index corresponding to T=298K that can be used to get S298 and H298.
        Incompatible with thin.
    thin : int
        Integer stride length to thin the formation energies. Useful for CPM if
        there is a dense grid and only a few points are desired. Incompatible with idx.

    Returns
    -------
    np.array
        Array of formation properties. Array can be changed in size by idx and thin.
    """
    composition = struct.composition.get_el_amt_dict()
    ref_energy = np.zeros(tprops['T'].size)
    total_amount = np.sum([amnt for amnt in composition.values()])
    for comp, qty in composition.items():
        # surface of reference energy weighted by composition
        ref_energy += refstate_dict[comp][prop]*qty/total_amount
    delta_prop = tprops[prop] - ref_energy
    if (idx is not None) and (thin is not None):
        raise ValueError("Both `idx` and `thin` cannot be specified")
    if idx is not None:
        delta_prop = delta_prop[idx]
    if thin is not None:
        delta_prop = delta_prop[::thin]
    return delta_prop

"""
Phonon analysis using phonopy
"""

from phonopy import Phonopy
from phonopy.interface.vasp import Vasprun as PhonopyVasprun
from pymatgen.io.phonopy import get_phonopy_structure
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
from dfttk.utils import J_per_mol_to_eV_per_atom
from scipy.integrate import trapz
import numpy as np


def get_f_vib_phonopy(structure, supercell_matrix, vasprun_path,
                     qpoint_mesh=(50, 50, 50), t_min=5, t_step=5, t_max=2000.0,):
    """
    Return F_vib(T) for the unitcell in eV/atom

    Parameters
    ----------
    structure : pymatgen.Structure
        Unitcell (not supercell) of interest.
    supercell_matrix : numpy.ndarray
        3x3 matrix of the supercell deformation, e.g. [[3, 0, 0], [0, 3, 0], [0, 0, 3]].
    vasprun_path : str
        String pointing to a vasprun.xml file from a force constants run
    qpoint_mesh : list
        Mesh of q-points to calculate thermal properties on.
    t_min : float
        Minimum temperature
    t_step : float
        Temperature step size
    t_max : float
        Maximum temperature (inclusive)

    Returns
    -------
    tuple
        Tuple of (temperature, F_vib, S_vib, Cv_vib, force_constants)

    """
    # get the force constants from a vasprun.xml file
    vasprun = PhonopyVasprun(vasprun_path)
    force_constants, elements = vasprun.read_force_constants()

    ph_unitcell = get_phonopy_structure(structure)
    ph = Phonopy(ph_unitcell, supercell_matrix)
    # set the force constants we found
    ph.set_force_constants(force_constants)
    # calculate the thermal properties
    ph.run_mesh(qpoint_mesh)
    ph.run_thermal_properties(t_min=t_min, t_max=t_max, t_step=t_step)
    # the thermal properties are for the unit cell
    tp_dict = ph.get_thermal_properties_dict()
    temperatures = tp_dict['temperatures']
    # convert the units into our expected eV/atom-form (and per K)
    f_vib = tp_dict['free_energy'] * J_per_mol_to_eV_per_atom*1000
    s_vib = tp_dict['entropy'] * J_per_mol_to_eV_per_atom
    cv_vib = tp_dict['heat_capacity'] * J_per_mol_to_eV_per_atom
    return temperatures, f_vib, s_vib, cv_vib, ph.force_constants

def get_phonon_band(structure, supercell_matrix, force_constants, band_paths=None, npoints=51, labels=None,
                    save_data=False, save_fig=False):
    '''
    Return the phonon bandstructure

    Parameters
    ----------
    structure : pymatgen.Structure
        Unitcell (not supercell) of interest.
    supercell_matrix : numpy.ndarray
        3x3 matrix of the supercell deformation, e.g. [[3, 0, 0], [0, 3, 0], [0, 0, 3]].
    force_constants: list
        force constants
    band_paths :  list, multi dimention
        Sets of end points of paths, e.g. [[[0, 0, 0], [0.5, 0.5, 0], [0.5, 0.5, 0.5]], [[0.5, 0.25, 0.75], [0, 0, 0]]]
        If it equals None, it will determine the path automatically by phononpy
    npoints: int
        Number of q-points in each path including end points.
    labels: list of str
        The label of high symmetry points, if None, it will determine it automatically by phononpy
    save_data/save_fig: bool
        Determine if save the data/figure or not
    '''
    volume = structure.volume
    formula = structure.composition.reduced_formula
    filename = "{}-phonon-Vol{:.2f}".format(formula, volume)

    unitcell = get_phonopy_structure(structure)
    ph_band_obj = Phonopy(unitcell, supercell_matrix)
    ph_band_obj.set_force_constants(force_constants)

    if band_paths:
        qpoints, connections = get_band_qpoints_and_path_connections(band_paths, npoints=npoints)
        ph_band_obj.run_band_structure(qpoints, path_connections=connections, labels=labels)
    else:
        ph_band_obj.auto_band_structure()
    if save_fig:
        fig_band = ph_band_obj.plot_band_structure()
        fig_band.savefig(fname='{}-band.png'.format(filename))
        fig_band.close()
    if save_data:
        ph_band_obj.write_yaml_band_structure(filename='{}-band.yaml'.format(filename))
    return ph_band_obj

def get_phonon_dos(structure, supercell_matrix, force_constants, qpoint_mesh=(50, 50, 50), phonon_pdos=False,
                   save_data=False, save_fig=False):
    '''
    Return the phonon dos

    Parameters
    ----------
    structure : pymatgen.Structure
        Unitcell (not supercell) of interest.
    supercell_matrix : numpy.ndarray
        3x3 matrix of the supercell deformation, e.g. [[3, 0, 0], [0, 3, 0], [0, 0, 3]].
    force_constants: list
        force constants
    qpoint_mesh : list
        Mesh of q-points to calculate thermal properties on.
    phonon_pdos: bool
        Determine if calculate phonon pdos or not
    save_data/save_fig: bool
        Determine if save the data/figure or not
    '''
    volume = structure.volume
    formula = structure.composition.reduced_formula
    filename = "{}-phonon-Vol{:.2f}".format(formula, volume)

    unitcell = get_phonopy_structure(structure)
    ph_dos_obj = Phonopy(unitcell, supercell_matrix)
    ph_dos_obj.set_force_constants(force_constants)

    ph_dos_obj.run_mesh(qpoint_mesh)
    ph_dos_obj.run_total_dos()
    if save_fig:
        fig_dos = ph_dos_obj.plot_total_dos()
        fig_dos.savefig(fname='{}-dos.png'.format(filename))
        fig_dos.close()
    if save_data:
        ph_dos_obj.write_total_dos(filename='{}-dos.dat'.format(filename))

    if phonon_pdos:
        ph_dos_obj.run_mesh(qpoint_mesh, with_eigenvectors=True, is_mesh_symmetry=False)
        ph_dos_obj.run_projected_dos()
        if save_fig:
            ph_dos_obj.plot_projected_dos().savefig(fname='{}-pdos.png'.format(filename))
        if save_data:
            ph_dos_obj.write_projected_dos(filename='{}-pdos.dat'.format(filename))
    return ph_dos_obj

def get_phonon_band_dos(structure, supercell_matrix, force_constants, qpoint_mesh=(50, 50, 50), band_paths=None, 
                        npoints=51, labels=None, phonon_dos=True, phonon_band=True, phonon_pdos=False, 
                        save_data=False, save_fig=False):
    '''
    Return the phonon dos and band

    Parameters
    ----------
    structure : pymatgen.Structure
        Unitcell (not supercell) of interest.
    supercell_matrix : numpy.ndarray
        3x3 matrix of the supercell deformation, e.g. [[3, 0, 0], [0, 3, 0], [0, 0, 3]].
    force_constants: list
        force constants
    qpoint_mesh : list
        Mesh of q-points to calculate thermal properties on.
    band_paths :  list, multi dimention
        Sets of end points of paths, e.g. [[[0, 0, 0], [0.5, 0.5, 0], [0.5, 0.5, 0.5]], [[0.5, 0.25, 0.75], [0, 0, 0]]]
        If it equals None, it will determine the path automatically by phononpy
    npoints: int
        Number of q-points in each path including end points.
    labels: list of str
        The label of high symmetry points, if None, it will determine it automatically by phononpy
    phonon_dos/phonon_band/phonon_pdos: bool
        Determine if calculate dos/band/pdos or not
    save_data/save_fig: bool
        Determine if save the data/figure or not

    Returns
    -------
    '''
    ph_band_obj = None
    ph_dos_obj = None
    #for band
    if phonon_band:
        ph_band_obj = get_phonon_band(structure, supercell_matrix, force_constants, band_paths=band_paths,
                                      npoints=npoints, labels=labels, save_data=save_data, save_fig=save_fig)
    #for dos
    if phonon_dos:
        ph_dos_obj = get_phonon_dos(structure, supercell_matrix, force_constants, qpoint_mesh=qpoint_mesh,
                                    phonon_pdos=phonon_pdos, save_data=save_data, save_fig=save_fig)
    return (ph_band_obj, ph_dos_obj)

def phonon_stable(structure, supercell_matrix, force_constants, qpoint_mesh=(50, 50, 50), stable_tor=0.01):
    '''
    Judge the stability of structure from phonon

    Parameters
    ----------
    structure : pymatgen.Structure
        Unitcell (not supercell) of interest.
    supercell_matrix : numpy.ndarray
        3x3 matrix of the supercell deformation, e.g. [[3, 0, 0], [0, 3, 0], [0, 0, 3]].
    force_constants: list
        force constants
    qpoint_mesh : list
        Mesh of q-points to calculate thermal properties on.
    stable_tor: float
        The tolerance for the percentage of negative frequency. 
        If the percentage of negative frequency is lagre than the tor, then the structure is unstable

    Return
    ------
        structure_stability: bool
            True for stable, False for unstable
    '''
    structure_stability = True
    ph_dos_obj = get_phonon_dos(structure, supercell_matrix, force_constants, qpoint_mesh=qpoint_mesh)
    phonon_freq = ph_dos_obj._total_dos._frequency_points
    phonon_dos = ph_dos_obj._total_dos._dos
    freq_min = np.amin(phonon_freq)
    if freq_min < 0:
        integrate_full = trapz(phonon_dos, x=phonon_freq)
        ind_freq_le0 = np.where(phonon_freq < 0)
        integrate_le0 = trapz(phonon_dos[ind_freq_le0], x=phonon_freq[ind_freq_le0])
        if integrate_le0/integrate_full > stable_tor:
            structure_stability = False
    return structure_stability

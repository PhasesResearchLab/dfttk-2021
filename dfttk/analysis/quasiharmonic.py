"""
This module implements the Quasiharmonic approximation that can
be used to compute thermal properties.

It is based on pymatgen's QHA and further modified/refactored to abstract away the sources of
contributions to the Gibbs energy so that it may apply to the Debye models, phonon properties, etc.

"""
# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals, division, print_function

from collections import defaultdict

import numpy as np
from scipy.optimize import minimize
from pymatgen.analysis.eos import EOS, EOSError

from dfttk.analysis.thermal_electronic import calculate_thermal_electronic_contribution
#from dfttk.analysis.debye import DebyeModel
from dfttk.analysis.debye_ext import DebyeModel

__author__ = "Kiran Mathew, Brandon Bocklund"
__credits__ = "Cormac Toher"


class Quasiharmonic(object):
    """
    Class to perform quasiharmonic calculations.

    In principle, helps to abstract away where different energy contributions come from.

    Parameters
    ----------
    energies : list
        List of DFT energies in eV
    volumes : list
        List of volumes in Ang^3
    structure : pymatgen.Structure
        One of the structures on the E-V curve (can be any volume).
    dos_objects : list
        List of pymatgen Dos objects corresponding to the volumes. If passed, will enable the
        electronic contribution.
    f_vib : numpy.ndarray
        Array of F_vib(V,T) of shape (len(volumes), len(temperatures)). If absent, will use the Debye model.
    t_min : float
        Minimum temperature
    t_step : float
        Temperature step size
    t_max : float
        Maximum temperature (inclusive)
    eos : str
        Equation of state used for fitting the energies and the volumes.
        Options supported by pymatgen: "quadratic", "murnaghan", "birch", "birch_murnaghan",
        "pourier_tarantola", "vinet", "deltafactor", "numerical_eos". Default is "vinet".
    pressure : float
        Pressure to apply to the E-V curve/Gibbs energies in GPa. Defaults to 0.
    poisson : float
        Poisson ratio, defaults to 0.25. Only used in QHA
    bp2gru : float
        Fitting parameter for dBdP in the Gruneisen parameter. 2/3 is the high temperature
        value and 1 is the low temperature value. Defaults to 1.
    vib_kwargs : dict
        Additional keyword arguments to pass to the vibrational calculator
    """
    def __init__(self, energies, volumes, structure, dos_objects=None, F_vib=None, S_vib=None, C_vib=None,
                 t_min=5, t_step=5,
                 t_max=2000.0, eos="vinet", pressure=0.0, poisson=0.25,
                 bp2gru=1., vib_kwargs=None):
        self.energies = np.array(energies)
        self.volumes = np.array(volumes)
        self.natoms = len(structure)
        self.temperatures = np.arange(t_min, t_max+t_step, t_step)
        self.eos_name = eos
        self.pressure = pressure
        self.gpa_to_ev_ang = 1./160.21766208  # 1 GPa in ev/Ang^3
        self.eos = EOS(eos)

        # get the vibrational properties as a function of V and T
        #print ("xxxxxxxxxxxx volumes", self.volumes.shape, self.volumes)
        if F_vib is None:  # use the Debye model
            vib_kwargs = vib_kwargs or {}
            debye_model = DebyeModel(energies, volumes, structure, t_min=t_min, t_step=t_step,
                                     t_max=t_max, eos=eos, poisson=poisson, bp2gru=bp2gru, **vib_kwargs)
            self.F_vib = debye_model.F_vib  # vibrational free energy as a function of volume and temperature
            self.S_vib = debye_model.S_vib  # vibrational entropy as a function of volume and temperature
            self.C_vib = debye_model.C_vib  # vibrational heat capacity as a function of volume and temperature
            self.D_vib = debye_model.D_vib  # Debye temperature
        else:
            self.F_vib = F_vib
            self.S_vib = S_vib
            self.C_vib = C_vib
        #print ("xxxxxxxxxxxx volumes", F_vib.shape)


        # get the electronic properties as a function of V and T
        if dos_objects:
            # we set natom to 1 always because we want the property per formula unit here.
            thermal_electronic_props = [calculate_thermal_electronic_contribution(dos, t0=t_min, t1=t_max, td=t_step, natom=1) for dos in dos_objects]
            self.F_el = [p['free_energy'] for p in thermal_electronic_props]
        else:
            self.F_el = np.zeros((self.volumes.size, self.temperatures.size))

        # Set up the array of Gibbs energies
        # G = E_0(V) + F_vib(V,T) + F_el(V,T) + PV
        self.G = self.energies[:, np.newaxis] + self.F_vib + self.F_el + self.pressure * self.volumes[:, np.newaxis] * self.gpa_to_ev_ang

        # set up the final variables of the optimized Gibbs energies
        self.gibbs_free_energy = []  # optimized values, eV
        self.optimum_volumes = []  # in Ang^3
        self.optimize_gibbs_free_energy()

    def optimize_gibbs_free_energy(self):
        """
        Evaluate the gibbs free energy as a function of V, T and P i.e
        G(V, T, P), minimize G(V, T, P) wrt V for each T and store the
        optimum values.

        Note: The data points for which the equation of state fitting fails
            are skipped.
        """
        for temp_idx in range(self.temperatures.size):
            G_opt, V_opt = self.optimizer(temp_idx)
            self.gibbs_free_energy.append(float(G_opt))
            self.optimum_volumes.append(float(V_opt))

    def optimizer(self, temp_idx):
        """
        Evaluate G(V, T, P) at the given temperature(and pressure) and
        minimize it wrt V.

        1. Compute the  vibrational helmholtz free energy, A_vib.
        2. Compute the gibbs free energy as a function of volume, temperature
            and pressure, G(V,T,P).
        3. Preform an equation of state fit to get the functional form of
            gibbs free energy:G(V, T, P).
        4. Finally G(V, P, T) is minimized with respect to V.

        Args:
            temp_idx : int
            Index of the temperature of interest from self.temperatures

        Returns:
            float, float: G_opt(V_opt, T, P) in eV and V_opt in Ang^3.
        """
        G_V = self.G[:, temp_idx]

        # fit equation of state, G(V, T, P)
        try:
            eos_fit = self.eos.fit(self.volumes, G_V)
        except EOSError:
            return np.nan, np.nan
        # minimize the fit eos wrt volume
        # Note: the ref energy and the ref volume(E0 and V0) not necessarily
        # the same as minimum energy and min volume.
        volume_guess = eos_fit.volumes[np.argmin(eos_fit.energies)]
        min_wrt_vol = minimize(eos_fit.func, volume_guess)
        # G_opt=G(V_opt, T, P), V_opt
        return min_wrt_vol.fun, min_wrt_vol.x[0]


    def get_summary_dict(self):
        """
        Returns a dict with a summary of the computed properties.
        """
        d = defaultdict(list)
        d["pressure"] = self.pressure
        d["natoms"] = int(self.natoms)
        d["gibbs_free_energy"] = self.gibbs_free_energy
        d["temperatures"] = self.temperatures
        d["optimum_volumes"] = self.optimum_volumes
        d["volumes"] = self.volumes.tolist()
        d["helmholtz_energies"] = self.F_vib.tolist()

        try:
            d["entropies"] = self.S_vib.tolist()
        except:
            pass

        try:
            d["heat_capacities"] = self.C_vib.tolist()
        except:
            pass

        try:
            d["debye_temperatures"] = self.D_vib.tolist()
        except:
            pass

        return d

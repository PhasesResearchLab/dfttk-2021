"""
Calculate vibrational Helmholtz energies with the Debye model

Quasiharmonic Debye-Gruneisen model based on pymatgen's QHA and modified
to clean up and use the proper Gruneisen parameter and to factor out the Debye part from the QHA
"""
# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals, division, print_function

import numpy as np

from scipy.constants import physical_constants
from scipy.integrate import quadrature
from scipy.stats import gmean

from pymatgen.analysis.eos import EOS


__author__ = "Kiran Mathew, Brandon Bocklund"
__credits__ = "Cormac Toher"


class DebyeModel(object):
    """
    Calculate the vibrational free energy for volumes/temperatures using the Debye model

    Note that the properties are per unit formula!

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
        Poisson ratio, defaults to 0.363615, corresponding to the cubic scaling factor of 0.617 by Moruzzi
    gruneisen : bool
        Whether to use the Debye-Gruneisen model. Defaults to True.
    bp2gru : float
        Fitting parameter for dBdP in the Gruneisen parameter. 2/3 is the high temperature
        value and 1 is the low temperature value. Defaults to 1.
    mass_average_mode : str
        Either 'arithmetic' or 'geometric'. Default is 'arithmetic'

    """
    def __init__(self, energies, volumes, structure, t_min=5, t_step=5,
                 t_max=2000.0, eos="vinet", poisson=0.363615,
                 gruneisen=True, bp2gru=1., mass_average_mode='arithmetic'):
        self.energies = energies
        self.volumes = volumes
        self.structure = structure
        self.temperatures = np.arange(t_min, t_max+t_step, t_step)
        self.eos_name = eos
        self.poisson = poisson
        self.bp2gru = bp2gru
        self.gruneisen = gruneisen
        self.natoms = self.structure.composition.num_atoms
        self.kb = physical_constants["Boltzmann constant in eV/K"][0]
        # calculate the average masses
        masses = np.array([e.atomic_mass for e in self.structure.species]) * physical_constants["atomic mass constant"][0]
        if mass_average_mode == 'arithmetic':
            self.avg_mass = np.mean(masses)
        elif mass_average_mode == 'geometric':
            self.avg_mass = gmean(masses)
        else:
            raise ValueError("DebyeModel mass_average_mode must be either 'arithmetic' or 'geometric'")
        # fit E and V and get the bulk modulus(used to compute the Debye temperature)
        self.eos = EOS(eos)
        self.ev_eos_fit = self.eos.fit(volumes, energies)
        self.bulk_modulus = self.ev_eos_fit.b0_GPa  # in GPa

        self.calculate_F_el()

    def calculate_F_el(self):
        """
        Calculate the Helmholtz vibrational free energy

        """
        self.F_vib = np.zeros((len(self.volumes), self.temperatures.size ))
        for v_idx, vol in enumerate(self.volumes):
            for t_idx, temp in enumerate(self.temperatures):
                self.F_vib[v_idx, t_idx] = self.vibrational_free_energy(temp, vol)


    def vibrational_free_energy(self, temperature, volume):
        """
        Vibrational Helmholtz free energy, A_vib(V, T).
        Eq(4) in doi.org/10.1016/j.comphy.2003.12.001

        Args:
            temperature (float): temperature in K
            volume (float)

        Returns:
            float: vibrational free energy in eV
        """
        y = self.debye_temperature(volume) / temperature
        return self.kb * self.natoms * temperature * (9./8. * y + 3 * np.log(1 - np.exp(-y)) - self.debye_integral(y))


    def debye_temperature(self, volume):
        """
        Calculates the debye temperature.

        The anharmonic contribution is toggled by setting the anharmonic_contribution
        to True or False in the QuasiharmonicDebyeApprox constructor.

        Args:
            volume (float): in Ang^3

        Returns:
            float: debye temperature in K

        Notes
        -----
        The original code from pymatgen cites Toher [1], however the code here
        does not match the equation in that paper (the derivation seems
        incorrect). Chen and Sundman [2] have a clearer derivation in agreement
        with this work, except that their final equation is in terms of
        interatomic radii rather than volume.

        [1] Toher, Phys. Rev. B 90, 174107 (2014) doi:10.1016/j.comphy.2003.12.001
        [2] Chen and Sundman, Acta Materialia 49, 947--961 (2001) doi:10.1016/S1359-6454(01)00002-7

        """
        term1 = (2./3. * (1. + self.poisson) / (1. - 2. * self.poisson))**1.5
        term2 = (1./3. * (1. + self.poisson) / (1. - self.poisson))**1.5
        f = (3. / (2. * term1 + term2))**(1. / 3.)
        debye = 2.9772e-11 * (volume / self.natoms) ** (-1. / 6.) * f * np.sqrt(self.bulk_modulus/self.avg_mass)
        if self.gruneisen:
            # bp2gru should be the correction to the Gruneisen constant.
            # High temperature limit: 2/3
            # Low temperature limit: 1
            # take 0 K E-V curve properties
            dBdP = self.ev_eos_fit.b1  # bulk modulus/pressure derivative
            gamma = (1+dBdP)/2 - self.bp2gru  # 0K equilibrium Gruneisen parameter
            return debye * (self.ev_eos_fit.v0 / volume) ** (gamma)
        else:
            return debye


    @staticmethod
    def debye_integral(y):
        """
        Debye integral. Eq(5) in  doi.org/10.1016/j.comphy.2003.12.001

        Args:
            y (float): debye temperature/T, upper limit

        Returns:
            float: unitless
        """
        # floating point limit is reached around y=155, so values beyond that
        # are set to the limiting value(T-->0, y --> \infty) of
        # 6.4939394 (from wolfram alpha).
        factor = 3. / y ** 3
        if y < 155:
            integral = quadrature(lambda x: x ** 3 / (np.exp(x) - 1.), 0, y)
            return list(integral)[0] * factor
        else:
            return 6.493939 * factor

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
from scipy.interpolate import splrep, splev

from scipy.constants import physical_constants
import scipy.constants as scipy_constants

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
    def __init__(self, energies, volumes, structure, T=None, t_min=5, t_step=5,
                 t_max=2000.0, gruneisen_T0 = 0.0, debye_T0 = -1.0,
                 eos="vinet", poisson=0.363615,
                 gruneisen=True, bp2gru=2./3., mass_average_mode='arithmetic'):
        self.energies = energies
        self.volumes = volumes
        self.structure = structure
        if T is None: self.temperatures = np.arange(t_min, t_max+t_step, t_step)
        else: self.temperatures = T
        self.eos_name = eos
        self.poisson = poisson
        self.bp2gru = bp2gru
        self.gruneisen = gruneisen
        self.natoms = self.structure.composition.num_atoms
        self.kb = physical_constants["Boltzmann constant in eV/K"][0]
        self.gruneisen_T0 = gruneisen_T0
        if self.gruneisen_T0 is not None:
            if self.gruneisen_T0 < -99.0: self.gruneisen = False
        self.debye_T0 = debye_T0
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
        self.S_vib = np.zeros((len(self.volumes), self.temperatures.size ))
        self.C_vib = np.zeros((len(self.volumes), self.temperatures.size ))
        self.D_vib = np.zeros((len(self.volumes)))
        for v_idx, vol in enumerate(self.volumes):
            self.D_vib[v_idx] = self.debye_temperature(vol)
            for t_idx, temp in enumerate(self.temperatures):
                self.F_vib[v_idx, t_idx] = self.vibrational_free_energy(temp, vol)
                self.S_vib[v_idx, t_idx] = self.vibrational_entropy(temp, vol)
                self.C_vib[v_idx, t_idx] = self.vibrational_heat_capacity(temp, vol)


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
        if temperature <= 0.0: return self.kb * self.natoms * self.debye_temperature(volume) * 9./8.
        y = self.debye_temperature(volume) / temperature
        return self.kb * self.natoms * temperature * (9./8. * y + 3 * np.log(1 - np.exp(-y)) - self.debye_integral(y))


    def vibrational_entropy(self, temperature, volume):
        """
        Vibrational entropy, S_vib(V, T).

        Args:
            temperature (float): temperature in K
            volume (float)

        Returns:
            float: vibrational entropy in eV/K
        """
        if temperature <= 0.0: return 0.0
        y = self.debye_temperature(volume) / temperature
        return self.kb * self.natoms * ( -3 * np.log(1 - np.exp(-y)) + 4*self.debye_integral(y))


    def vibrational_heat_capacity(self, temperature, volume):
        """
        Vibrational heat capacity, C_vib(V, T).

        Args:
            temperature (float): temperature in K
            volume (float)

        Returns:
            float: vibrational heat capacity in eV/K
        """
        if temperature <= 0.0: return 0.0
        y = self.debye_temperature(volume) / temperature
        factor = 3. / y ** 3
        if y < 155:
            integral = quadrature(lambda x: x ** 4 *np.exp(x)/ (np.exp(x) - 1.)**2, 0, y)
            return 3*self.kb * self.natoms * list(integral)[0] * factor
        else:
            return self.kb * self.natoms * 4./5.*scipy_constants.pi**4 * factor


    def derivative(self, volume, der=1):
        vmin, vmax = min(self.volumes), max(self.volumes)
        vmin, vmax = (vmin - 0.01 * abs(vmin), vmax + 0.01 * abs(vmax))
        vfit = np.linspace(vmin, vmax, 100)

        spl = splrep(vfit,self.ev_eos_fit.func(vfit),k=5, s=3) # no smoothing, 3rd order spline
        ddy = splev(volume,spl,der=der) # use those knots to get second derivative
        return ddy

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

        hbar = scipy_constants.hbar #1.054571817e-34
        kB = scipy_constants.Boltzmann #1.38064852e-23
        pi = scipy_constants.pi
        A = (6*pi*pi)**(1/3)*hbar/kB #2.97721279650828e-11

        term1 = (2./3. * (1. + self.poisson) / (1. - 2. * self.poisson))**1.5
        term2 = (1./3. * (1. + self.poisson) / (1. - self.poisson))**1.5
        s = (3. / (2. * term1 + term2))**(1. / 3.) #0.6170015756491913 by default 

        if self.gruneisen:
            if self.gruneisen_T0 is not None:
                gamma = self.gruneisen_T0
            else:
                # take 0 K E-V curve properties
                # bp2gru should be the correction to the Gruneisen constant.
                # High temperature limit: 2/3
                # Low temperature limit: 1
                dBdP = self.ev_eos_fit.b1  # bulk modulus/pressure derivative
                gamma = (1+dBdP)/2 - self.bp2gru  # 0K equilibrium Gruneisen parameter

            if self.debye_T0 > 0.0:
                debye = self.debye_T0
            else:
                debye = s*A * (self.ev_eos_fit.v0*1.e-30/self.natoms) ** (1. / 6.) * np.sqrt(self.bulk_modulus*1e9/self.avg_mass)
            gamma *= volume/self.ev_eos_fit.v0
            debye = debye*(self.ev_eos_fit.v0 / volume) ** (gamma)
            return debye 
        else:
            bulk_modulus=160.2176621*volume*self.derivative(volume, der=2)
            pressure=-160.2176621*self.derivative(volume, der=1)
            #V^(1/3) [((1+λ)2)/3*∂E/∂V + V (∂^2 E)/(∂V^2 )]
            #V^(1/3) [((1+3*(x-1))2)/3*∂E/∂V + V (∂^2 E)/(∂V^2 )]
            shift_pressure = -(1+3*(self.bp2gru-1))*2/3*pressure
            if self.debye_T0 > 0.0: return self.debye_T0*(volume/self.ev_eos_fit.v0)**(1./6.)* \
                np.sqrt((shift_pressure+bulk_modulus)/self.bulk_modulus)
            return s*A * (volume*1.e-30/self.natoms) ** (1. / 6.) * np.sqrt((shift_pressure+bulk_modulus)*1e9/self.avg_mass)
                #t0 = s*A * (self.ev_eos_fit.v0*1.e-30/self.natoms) ** (1. / 6.) * np.sqrt(self.bulk_modulus*1e9/self.avg_mass)
                #debye *= self.debye_T0/t0


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

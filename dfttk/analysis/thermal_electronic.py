#!/usr/bin/env python

from __future__ import division
import math
import numpy as np
from scipy.constants import physical_constants
from scipy.optimize import brentq
from scipy.integrate import cumtrapz, trapz

k_B = physical_constants['Boltzmann constant in eV/K'][0]

# TODO: check and fix all XXX and TODOs

def getdos(dos, xdn, xup, dope, dos_grid_size, gaussian_grid_size): # Line 186
    """

    Parameters
    ----------
    dos : pymatgen.electronic_structure.dos.Dos
        DOS object from pymatgen
    xdn : float
        Minimum energy for integration
    xup : float
        Maximum energy for integration
    dope : float
        Number of electrons to dope (negative means to remove electrons, positive means add electrons)
    dos_grid_size : int
        Number of DOS points have in the energy/density grid
    gaussian_grid_size : int
        Number of Gaussian points to use in the grid mesh around the Fermi energy

    Returns
    -------
    tuple
        Tuple of a (float, float, array, array) of the number of electrons,
        Fermi level shift due to doping, and the arrays of energies and densities on a grid.
    """
    e_fermi = dos.efermi
    eup = np.max(dos.energies) - e_fermi
    edn = np.min(dos.energies) - e_fermi
    n_dos = len(dos.energies) # number of points in DOS
    vde = (eup - edn)/(n_dos-1) # change in energy per step
    # linearize: sometimes rounding errors in DOSCAR
    dos_energies = np.linspace(edn, eup, n_dos)
    grid_dos = np.zeros(dos_grid_size)

    vaspEdos = np.array(dos.get_densities())

    xdn = max(xdn,edn)
    xup = min(xup,eup)
    # TODO: from original code seems like an arbitrary number
    # set the minimum energy to the highest energy where the DOS is 0, as long as the energy is less than -15eV
    # easy area to optimize out by removal, but need a test case
    for i in range(dos_energies.size):
        if dos_energies[i] <-15.0 and vaspEdos[i]==0.0:
            xdn = dos_energies[i]
        if dos_energies[i] > -15.0:
            break
    grid_energies = np.linspace(xdn, xup, dos_grid_size, dtype=float)
    xde = (xup - xdn)/(dos_grid_size - 1)

    # TODO: What are iBoF and eBoF
    iBoF = -1
    for i,eband in enumerate(dos_energies):
        if eband > 0.0 and eband <=vde and vaspEdos[i]==0.0:
            iBoF = 1
            continue
        if iBoF == 1:
            if vaspEdos[i]!=0.0:
                iBoF = i-1
                break

    eBoF = -1.0
    if iBoF>0:
        eBoF = dos_energies[iBoF]
        espr = vaspEdos[iBoF+2]-vaspEdos[iBoF+1]
        if espr>0.0:
            espr = vaspEdos[iBoF+1]/espr*vde
            if (espr < vde):
                eBoF = dos_energies[iBoF+1] - espr

    if gaussian_grid_size != 0.0:
        grid_energies[0] = xdn
        # line 238
        sigma = (xup-xdn) / gaussian_grid_size
        fac = gaussian_grid_size / (math.sqrt(2.0 * math.pi))
        for i in range(1, dos_grid_size):
            if iBoF<0:
                f1 = fac*math.exp(-0.5*(grid_energies[i-1]/sigma)**2)
                yde = 2.0*xde/(1.0+f1)
            else:
                f1 = fac*math.exp(-0.5*(grid_energies[i-1]/sigma)**2)
                f2 = fac*math.exp(-0.5*((grid_energies[i-1]-eBoF)/sigma)**2)
                yde = 3.0*xde/(1.0+f1+f2)
            grid_energies[i] = grid_energies[i-1] + yde


    for i in range(0, dos_grid_size):
        tx = grid_energies[i]
        kx = int((tx-edn)/vde) # Converted to int, remember the type!
        kx = max([kx,0])
        kx = min([n_dos-2, kx])
        # handling near the top of valence band
        if vaspEdos[kx+1]==0.0 and dos_energies[kx+1]>0.0 and dos_energies[kx+1]<vde:
            if tx >= 0.0:
                grid_dos[i] = 0.0
            else:
                grid_dos[i] = vaspEdos[kx] * tx / dos_energies[kx]
        elif eBoF > 0.0 and vaspEdos[kx]==0.0 and dos_energies[kx+1]-eBoF<vde and dos_energies[kx+1]-eBoF>0.0:
            if tx <= eBoF:
                grid_dos[i] = 0.0
            else:
                grid_dos[i] = vaspEdos[kx + 1] * (tx - eBoF) / (dos_energies[kx + 1] - eBoF)
        else:
            grid_dos[i] = vaspEdos[kx] + (vaspEdos[kx + 1] - vaspEdos[kx]) / vde * (grid_energies[i] - dos_energies[kx])

    # find undoped number of electrons by integrating on the dos until the
    # energies change sign, then interpolate the sign change step
    ados = cumtrapz(grid_dos, grid_energies, initial=0.0)
    for i in range(0, dos_grid_size-1):
        if grid_energies[i]*grid_energies[i+1]<=0.0:
            n_electrons = ados[i] - grid_energies[i]/(grid_energies[i+1]-grid_energies[i])*(ados[i+1]-ados[i])
            break

    dF = 0.0
    if dope != 0.0:
        n_electrons += dope
        for i in range(0, dos_grid_size-1):
            if (ados[i] - n_electrons)*(ados[i+1]-n_electrons) < 0.0:
                if i == (dos_grid_size-1) or ados[i] == ados[i+1]:
                    # we are doping too much
                    raise ValueError('Too much doping')
                dF = (n_electrons-ados[i])/(ados[i+1] - ados[i])*(grid_energies[i+1] - grid_energies[i])+grid_energies[i]
                # dF is the shift in the Fermi energy due to doping
        grid_energies = grid_energies - dF # This is done in a loop (line 289), but I think we can do without

    if gaussian_grid_size != 0.0 and abs(dope)>0.0001:
        grid_energies[0] = xdn - dF
        sigma = (xup-xdn) / gaussian_grid_size
        fac = gaussian_grid_size / (math.sqrt(2.0 * math.pi))
        for i in range(1, dos_grid_size):
            if iBoF<0:
                f1 = fac*math.exp(-0.5*(grid_energies[i-1]/sigma)**2)
                yde = 2.0*xde/(1.0+f1)
            else:
                f1 = fac*math.exp(-0.5*(grid_energies[i-1]/sigma)**2)
                if dF < eBoF:
                    f2 = fac*math.exp(-0.5*((grid_energies[i-1]-eBoF+dF)/sigma)**2)
                else:
                    f2 = fac*math.exp(-0.5*((grid_energies[i-1]+dF)/sigma)**2)
                yde = 3.0*xde/(1.0+f1+f2)
            grid_energies[i] = grid_energies[i-1]+yde

        for i in range(0, dos_grid_size):
            tx = grid_energies[i] + dF
            kx = int((tx-edn)/vde) # Converted to int, remember the type!
            kx = max([kx,0])
            kx = min([n_dos-2, kx])
            if vaspEdos[kx+1]==0.0 and dos_energies[kx+1]>0.0 and dos_energies[kx+1]<vde:
                # handling near the Top of valence band
                if tx >= 0.0:
                    grid_dos[i] = 0.0
                else:
                    grid_dos[i] = vaspEdos[kx] * tx / dos_energies[kx]
            elif eBoF > 0.0 and vaspEdos[kx]==0.0 and dos_energies[kx+1]-eBoF<vde and dos_energies[kx+1]-eBoF>0.0:
                # handling near the bottom of conduction band
                if tx <= eBoF:
                    grid_dos[i] = 0.0
                else:
                    grid_dos[i] = vaspEdos[kx + 1] * (tx - eBoF) / (dos_energies[kx + 1] - eBoF)
            else:
                grid_dos[i] = vaspEdos[kx] + (vaspEdos[kx + 1] - vaspEdos[kx]) / vde * (tx - dos_energies[kx])

    ados = cumtrapz(grid_dos, grid_energies, initial=0.0)
    for i in range(0, dos_grid_size-1):
        if grid_energies[i]*grid_energies[i+1]<=0.0:
            n_electrons = ados[i] - grid_energies[i]/(grid_energies[i+1]-grid_energies[i])*(ados[i+1]-ados[i])
            break

    return n_electrons, dF, grid_energies, grid_dos


def gfind(mu_el, e, dos, n_electrons, beta):
    """
    Calculate the energy given chemical potential

    """
    tc = beta*(e-mu_el)
    tc = tc[np.where(tc<200)]
    k = len(tc)
    fn = dos[0:k]/(np.exp(tc[0:k])+1.0)
    return trapz(fn, e[0:k]) - n_electrons


def calculate_internal_energy(mu_el, energy, density, beta): # line 471
    tc = beta * (energy - mu_el)
    tc[tc>200] = 200
    fn = density * energy / (np.exp(tc) + 1.0)
#    fn[tc>200] = 0
    u = trapz(fn, energy)
    return u


def calculate_entropy(mu_el, energy, density, beta):
    tc = beta * (energy - mu_el)
    tc[tc>200] = 200
    tf = 1.0/(np.exp(tc)+1.0)  # To avoid RuntimeWarning: invalid value encountered in multiply
    tf1 = 1.0 - tf + 1.e-60
    fn = density * (tf * np.log(tf) + tf1 * np.log(tf1))
    fn[tc>200] = 0
    fn[tc<-200] = 0
    s = trapz(fn, energy)
    return -s*k_B


def calculate_thermal_electronic_contribution(dos, t0=0, t1=2000, td=5, xdn=-100, xup=100, ndosmx=10001, dope=0.0, natom=1, gaussian=1000):
    """
    Calculate thermal electronic contribution from pymatgen Dos objects

    Parameters
    ----------
    dos : pymatgen.electronic_structure.dos.Dos
        DOS object
    t0 : float
        Start temperature
    t1 : float
        Final temperature
    td : float
        Temperature step size
    xdn : float
        Minimum energy of the DOS to consider
    xup : float
        Maximum energy of the DOS to consider
    ndosmx : int
        Size of grid to interpolate the DOS on
    dope : float
        Doping level
    natom : int
        Number of atoms in the cell
    gaussian : int
        Number of grid points in the Gaussian mesh near the Fermi energy

    Returns
    -------

    """
    n_electrons, fermi_shift, e, dos = getdos(dos, xdn, xup, dope, ndosmx, gaussian)

    # for all temperatures
    nT = int(np.round((t1-t0)/td+1.0))
    T = np.arange(t0,t1+td,td)
    chemical_potential = np.zeros(nT)
    gmu0 = 0.0
    beta = 1.0/(T*k_B)
    beta[0] = 1.0e30
    for i, t in enumerate(T):
        chemical_potential[i] = brentq(gfind, gmu0-5.0, gmu0+5.0, args=(e, dos, n_electrons, beta[i]), maxiter=10000)
    U_el = calculate_internal_energy(chemical_potential[:, np.newaxis], e[np.newaxis, :], dos[np.newaxis, :], beta[:, np.newaxis])
    S_el = calculate_entropy(chemical_potential[:, np.newaxis], e[np.newaxis, :], dos[np.newaxis, :], beta[:, np.newaxis])
    C_el = np.gradient(U_el, td, edge_order=2)
    C_el[0] = 0

    # construct a dictionary of results
    results = {
        'temperature': T,
        'internal_energy': U_el/natom,
        'free_energy': (U_el-T*S_el-U_el[0])/natom,
        'entropy': S_el/natom,
        'heat_capacity': C_el/natom,
        'chemical_potential': chemical_potential,
        'n_electrons': n_electrons,
    }

    return results

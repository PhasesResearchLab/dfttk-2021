# -------- energy in eV, temperature in K
# assume every variable starting with a-h  and o-z are real numbers
# common block named comcon
from __future__ import division
import sys
import gzip
import os
from os import walk
import subprocess
import math
import copy
import json
import pickle
import numpy as np
from scipy.constants import physical_constants
import scipy.constants as scipy_constants 
from scipy.optimize import brentq, curve_fit
from scipy.integrate import cumtrapz, trapz, simps
from scipy.interpolate import interp1d, splev, splrep, BSpline
from scipy.integrate import quadrature
from scipy.interpolate import UnivariateSpline
from atomate.vasp.database import VaspCalcDb
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import dfttk.pyphon as ywpyphon
from dfttk.utils import sort_x_by_y
from dfttk.analysis.ywplot import myjsonout
from dfttk.analysis.ywutils import get_rec_from_metatag, get_used_pot
from dfttk.analysis.ywutils import formula2composition, reduced_formula, MM_of_Elements
import warnings


k_B = physical_constants['Boltzmann constant in eV/K'][0]


def substr(str1, str2, pos):
  try:
    if str1.index(str2)==pos:
        #print("idx=",str1.index(str2))
        return True
    else:
        return False
  except ValueError:
    return False

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

def isint(value):
  try:
    int(value)
    return True
  except ValueError:
    return False


# this is a FORTRAN function (e.g. 1 return value)
def pregetdos(f): # Line 186
    """
    to make the code can also handle WIEN2k dos in the unit of eV

    Parameters
    ----------
    f : file descriptor for the DOS file

    Returns
    -------
    xdn : lower energy to integrate over?
    xup : higher energy to integrate over?
    vde : band energy intercal
    e (array): band energy mesh the Fermi energy has been shifted to zero
    DOS (array) : e dos
    """
    # read the file
    lines = f.readlines() # read in all lines then determine is it is WIEN2k DOS (in the unit eV) file or VASP DOS file
    # now the first line should be the one with the data, lets remove the one into its own special line
    tmp = lines[0]
    if substr(tmp,"#  BAND", 0):
        tmp = lines[1]
        tmp1 = lines[2]
        if substr(tmp, "#EF=",0) and substr(tmp1, "# ENERGY",0):
            tmp1 = tmp[31:43].replace("NENRG=","")
            if isint(tmp1):
                n_dos = int(tmp1)
                tmp = lines[2]
                lines = lines[3:n_dos+3]
                wienEdos = np.zeros(n_dos)
                ve = np.zeros(n_dos)
                for i, l in enumerate(lines):
                    split_l = l.split(' ')
                    split_l = [k for k in split_l if k != '']
                    ve[i], wienEdos[i] = (float(split_l[0]), float(split_l[1]))
                edn = ve[0]
                eup = ve[n_dos-1]
                ve = np.linspace(edn, eup, n_dos)
                vde = (eup - edn)/(n_dos-1) # This appears to be the change of v per electron, so what is v? Voltage in eV?
                return edn, eup, vde, ve, wienEdos

    tmp = lines[5]
    data_line = tmp[0:32].split(' ') #n_dos >10000, no space left before it in VASP
    data_line.extend(tmp[32:].split(' '))
    # filter out empty spaces
    data_line = [k for k in data_line if k != '']
    #print (data_line)
    eup, edn, n_dos, eFermi = (float(data_line[0]),
                           float(data_line[1]),
                           int(data_line[2]),
                           float(data_line[3])) # we're leaving the last number behind
    lines = lines[6:n_dos+6]
    # line 197 goes to line 209

    eup = eup - eFermi
    edn = edn - eFermi
    vde = (eup - edn)/(n_dos-1) # This appears to be the change of v per electron, so what is v? Voltage in eV?

    # vectors
    ve = np.linspace(edn, eup, n_dos)
    vaspEdos = np.zeros(n_dos)

    for i, l in enumerate(lines):
        # why do we need to do this?
        split_l = l.split(' ')
        # filter again
        split_l = [k for k in split_l if k != '']
        if len(split_l)>=5: #spin polarized
            t, vaspEdos[i], y, vdos, x = (float(split_l[0]), float(split_l[1]), float(split_l[2]), float(split_l[3]), float(split_l[4]))
            vaspEdos[i] += y
        else:
            t, vaspEdos[i], vdos = (float(split_l[0]), float(split_l[1]), float(split_l[2]))
    _eFermi = CBMtoVBM(ve, vaspEdos)

    return edn-_eFermi, eup-_eFermi, vde, ve-_eFermi, vaspEdos


def CBMtoVBM(ve, vaspEdos):
    # move eFermi to VBM if it in CBM
    vde = ve[1] - ve[0]
    for i, dos in enumerate(vaspEdos):
        if ve[i] >= -vde: break
        if dos!=0.0:
            _eFermi = ve[i]
    if _eFermi < -3*vde:
        print ("Fermi energy shifted from CBM", 0.0, "to VBM", _eFermi)
        return _eFermi+vde
    else: return 0.0


def getdos(xdn, xup, dope, NEDOS, gaussian, edn, eup, vde, ve, tdos): # Line 186
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

    for i,energy in enumerate(ve):
      if energy <-15.0 and tdos[i]==0.0:
        xdn = energy

    n_dos = len(tdos)
    idx = closest(ve,0.0)
    for i in range(idx,n_dos):
        if tdos[i]!=0.0:
            iBoF = i-1
            break

    eBoF = -1.0
    if iBoF>=idx:
      eBoF = ve[iBoF]
      espr = tdos[iBoF+2]-tdos[iBoF+1]
      if espr>0.0:
        espr = tdos[iBoF+1]/espr*vde
        if (espr < vde):
          eBoF = ve[iBoF+1] - espr

    #print("eBoF=", eBoF)
    xdn = max(xdn,edn)
    xup = min(xup,eup)

    e = np.linspace(xdn,xup,NEDOS,dtype=float)
    if gaussian != 0.0:
      e = remesh(xdn, xup, gaussian, 0.0, eBoF, NEDOS)

    dos = refdos(eBoF, 0.0, vde, edn, e, ve, tdos)
    ados = cumtrapz(dos, e, initial=0.0)
    idx = closest(e,0.0)
    for idx1 in range(idx-1, 0, -1):
        if ados[idx1] != ados[idx] : break
    NELECTRONS = ados[idx] - e[idx]/(e[idx1]-e[idx])*(ados[idx1]-ados[idx])

    dF = 0.0
    if dope != 0.0:
        NELECTRONS = NELECTRONS+dope
        idx = closest(ados,NELECTRONS)
        for idx1 in range(idx-1, 0, -1):
            if ados[idx1] != ados[idx] : break
        #if idx == (NEDOS-1) or ados[idx] == ados[idx+1]:
        #print ("NELECTRONS=", NELECTRONS, "idx=", idx, ados[idx], "idx1=", idx1, ados[idx1], "NEDOS=", NEDOS)
        if idx1 <= 0 or idx >= (NEDOS-1) or ados[idx] == ados[idx1]:
            print ("NELECTRONS=", NELECTRONS, "idx=", idx, ados[idx], "idx1=", idx1, ados[idx1], "NEDOS=", NEDOS)
            # we are dopidxng too much
            raise ValueError('Too much doping')
        dF = (NELECTRONS-ados[idx])/(ados[idx1] - ados[idx])*(e[idx1] - e[idx])+e[idx]
                # dF is the shift in the Fermi energy due to doping
        e = e - dF # This is done in a loop (line 289), but I think we can do without

    if gaussian != 0.0 and abs(dope)>0.0001: # why did I do this ***********************
    #if gaussian != 0.0:
      e = remesh(xdn, xup, gaussian, dF, eBoF, NEDOS)

    dos = refdos(eBoF, dF, vde, edn, e, ve, tdos)
    edos = e*dos
    ados = cumtrapz(dos, e, initial=0.0)
    energy = cumtrapz(edos, e, initial=0.0)
    idx = closest(e,0.0)
    NELECTRONS = ados[idx] - e[idx]/(e[idx+1]-e[idx])*(ados[idx+1]-ados[idx])
    E0 = energy[idx] - e[idx]/(e[idx+1]-e[idx])*(energy[idx+1]-energy[idx])

    return NELECTRONS, E0, dF, e, dos, eBoF

def remesh(xdn, xup, gaussian, dF, eBoF, NEDOS):
    """
    refine the dos mesh by using denser mesh around the 0 K Fermi energy in order to decrease the numerical uncertainty
    Parameters
    ----------
    eBoF : Conduction band minimum
    dF : Fermi energy change due to doping
    ve : original e mesh
    NEDOS : original e mesh
    gaussian : parameter used to refine the e mesh near the Fermi energy

    Return
    ------
    e : refined e mesh
    """

    e = np.zeros(NEDOS)
    e[0] = xdn - dF
    xde = 2.0*(xup - xdn)/(NEDOS-1)
    if eBoF>0.0:
        xde = 3.0*(xup - xdn)/(NEDOS-1)
    sigma = -0.5*(gaussian/(xup-xdn))**2
    fac = gaussian/(math.sqrt(2.0*math.pi))
    for i in range(1,NEDOS):
        f1 = 1.0 + fac*math.exp(sigma*(e[i-1])**2)
        if eBoF>0.0:
          if dF < eBoF:
             f1 += fac*math.exp(sigma*((e[i-1]-eBoF+dF))**2)
          else:
             f1 += fac*math.exp(sigma*((e[i-1]+dF))**2)
        e[i] = e[i-1]+xde/f1
    return e

def refdos(eBoF, dF, vde, edn, e, ve, tdos):
    """
    refine the dos mesh by using denser mesh around the 0 K Fermi energy in order to decrease the numerical uncertainty
    Parameter
    ---------
    eBoF : Conduction band minimum
    dF : Fermi energy change due to doping
    e : refined e mesh
    ve : original e mesh
    tdos : original e dos
    Return
    ------
    dos : refined e dos
    """

    dos = np.zeros(len(e))
    n_dos = len(tdos)
    for i in range(0, len(e)):
        tx = e[i] + dF
        kx = int((tx-edn)/vde) # Converted to int, remember the type!
        kx = max([kx,0]) # XXX: is this translated correctly? What is the 1 in fortran?
        kx = min([n_dos-2, kx]) # TODO: the ndos-1 was here before. could be a source of error
        if tdos[kx+1]==0.0 and ve[kx+1]>0.0 and ve[kx+1]<vde:
          # handling near the Top of valence band
          if tx >= 0.0:
            dos[i] = 0.0
          else:
            dos[i] = tdos[kx]*tx/ve[kx]
            #dos[i] = tdos[kx]*(tx/ve[kx])**2
        elif eBoF > 0.0 and tdos[kx]==0.0 and ve[kx+1]-eBoF<vde and ve[kx+1]-eBoF>0.0:
          # handling near the bottom of conduction band
            if tx <= eBoF:
              dos[i] = 0.0
            else:
              dos[i] = tdos[kx+1]*(tx-eBoF)/(ve[kx+1]-eBoF)
        else:
          dos[i] = tdos[kx] + (tdos[kx+1] - tdos[kx])/vde*(tx - ve[kx])
    return dos

def closest(e,val):
    """
    find the index of the band energy which is the close to the energy val
    Parameters
    ----------
    e : float
    array of band energy for the e dos
    val : given value of band energy

    Return
    ------
    index of e that closest to the energy val
    """
    idx = np.abs(e-val).argmin()
    if e[idx] < val:
        idx = idx + 1
    return idx

def gfind(mu_el, pe, pdos, NELECTRONS, Beta, IntegrationFunc=trapz):
    """
    Calculate the number of electron difference from 0K given chemical potential. the purpose is the find the
    chemical potential to make zero of number of electron difference from 0K

    Parameters
    ----------
    mu_el : chemical potential, :math:`\mu`, in the Fermi distribution
    pe : eigenenergies
    pdos : density of states (:math:`n(\varepsilon) \varepsilon`
    NELECTRONS : Total number of electrons in the system at 0K
    Beta : :math:`\frac{1}{T*k_{B}}`

    Returns
    -------
    The number of electron difference from 0K given chemical potential
    """

    tc = Beta*(pe-mu_el)
    tc = tc[np.where(tc<200)]
    k = len(tc)
    fn = pdos[0:k]/(np.exp(tc[0:k])+1.0)
    return IntegrationFunc(fn, pe[0:k])- NELECTRONS


# line 363
def caclf(pe, pdos, NELECTRONS, Beta, mu_ref=0.0, dF=0.0, IntegrationFunc=trapz): #line 363
    """
    Calculate thermal free energy from electronic density of states (e DOS)

    Parameters
    ----------
    pe : band energy array
    pdos : e DOS
    NELECTRONS : total number of electrons
    Beta : 1/(kB*T)

    Returns
    -------
    electron chememical potential, internal energy, entropy, carrier amount, coefficient to cal Seebeck
    """

    #print ("dF=", dF)
    if 1==1:
        deltaE = 2
        for i in range(8):
            try:
                mu_el = brentq(gfind, mu_ref-deltaE, mu_ref+deltaE, args=(pe, pdos, NELECTRONS, Beta, IntegrationFunc), maxiter=10000)
                break
            except:
                deltaE *= 2
    else:
        t0 = mu_ref
        d0 = gfind(t0, pe, pdos, NELECTRONS, Beta, IntegrationFunc)
        if d0 > 0.0: td = -0.1
        elif d0 <0.0: td = 0.1
        else: return t0
        for i in range(999):
            t1 = t0 + td
            d1 = gfind(t1, pe, pdos, NELECTRONS, Beta, IntegrationFunc)
            if d1*d0 < 0.0: break
            elif d1*d0 == 0.0: break
            t0 = t1
            d0 = d1
            td = td + td
        for i in range(999):
            t2 = (t0 + t1)*0.5
            d2 = gfind(t2, pe, pdos, NELECTRONS, Beta, IntegrationFunc)
            if d2*d0 < 0.0:
                t1 = t2
                d1 = d2
            else:
                t0 = t2
                d0 = d2
            if abs(t1-t0) <1.e-8:
                mu_el = 0.5*(t0+t1)
                break
    tc = Beta*(pe-mu_el)
    tc = tc[np.where(tc<200)]
    k1 = len(tc)
    tf = 1.0/(np.exp(tc)+1.0)
    fn = pdos[0:k1]*pe[0:k1]*tf
    u = IntegrationFunc(fn, pe[0:k1])

    k0 = closest(tc,-200)
    tf0 = tf[k0:]
    pdos = pdos[k0:k1]
    pe = pe[k0:k1]
    tf1 = 1.0 - tf0 + 1.e-60 # 1.e-60 is used to avoid log exception
    fn = pdos*(tf0*np.log(tf0)+tf1*np.log(tf1))
    s = IntegrationFunc(fn, pe)

    tf = tf0*(1.0-tf0)
    fn = pdos*tf
    fn2 = pdos*tf*(pe-mu_el)
    Q_el = IntegrationFunc(fn, pe)
    Q_p = IntegrationFunc(fn[pe<=dF], pe[pe<=dF])
    Q_e = IntegrationFunc(fn[pe>dF], pe[pe>dF])
    Y_el = IntegrationFunc(fn2, pe)

    fn = pdos*(pe-mu_el)*tf
    if Q_el!=0.0:
        e_ = IntegrationFunc(fn, pe)/Q_el
        fn = pdos[0:k1]*(pe[0:k1]-mu_el-e_)**2*tf
        cv = IntegrationFunc(fn, pe[0:k1])
    else:
        cv = 0.0
    fn = pdos[0:k1]*(pe[0:k1]-mu_el)**2*tf
    c_mu = IntegrationFunc(fn, pe[0:k1])

#   hole/electron concentration by effective carrier
    tf = tf0*(1.0-tf0)
    fn = pdos*tf
    f2 = interp1d(pe, fn, kind='linear')
    fmu = f2(mu_el)
    x = np.hstack([pe[pe<mu_el],mu_el])
    y = np.hstack([fn[pe<mu_el],fmu])
    W_p = IntegrationFunc(y,x)
    x = np.hstack([mu_el, pe[pe>mu_el]])
    y = np.hstack([fmu, fn[pe>mu_el]])
    W_e = IntegrationFunc(y,x)
    #W_e = IntegrationFunc(fn[pe>mu_el], pe[pe>mu_el])
    #W_e = IntegrationFunc(fn[pe>dF], pe[pe>dF])

#   hole/electron concentration by alternative difination
    fn = pdos*(1.0-tf0)
    f2 = interp1d(pe, fn, kind='linear')
    #fmu = f2(mu_el)
    #x = np.hstack([pe[pe<mu_el],mu_el])
    #y = np.hstack([fn[pe<mu_el],fmu])
    try:
        fmu = f2(dF)
    except:
        fmu = 0.
    x = np.hstack([pe[pe<dF],dF])
    y = np.hstack([fn[pe<dF],fmu])
    Y_p = IntegrationFunc(y,x)

    fn = pdos*tf0
    f2 = interp1d(pe, fn, kind='linear')
    #fmu = f2(mu_el)
    #x = np.hstack([mu_el, pe[pe>mu_el]])
    #y = np.hstack([fmu, fn[pe>mu_el]])
    #print ("mu_el", mu_el, dF)
    try:
        fmu = f2(dF)
    except:
        fmu = 0.
    x = np.hstack([dF, pe[pe>dF]])
    y = np.hstack([fmu, fn[pe>dF]])
    Y_e = IntegrationFunc(y,x)

    return mu_el, u, -s*k_B, cv*k_B*Beta*Beta, Q_el, Y_el, Q_p, Q_e, c_mu*k_B*Beta*Beta, W_p, W_e, Y_p, Y_e


def T_remesh(t0, t1, td, _nT=-1):
    T = []
    if td > 0:
        for t in np.arange(t0,t1+td, td):
            T.append(t)
        return np.array(T)

    if _nT <= 0: nT = 51
    else: nT = _nT
    a = 100./nT
    dT_new = abs(td)/(1+(nT-1)*0.5*a)
    for i in range (nT):
      T.append(t0+i*dT_new*(1+i*a))
    T = np.array(T)
    p = (t1-t0)/(max(T)-t0)
    for i in range (nT):
      T[i] = round((T[i]-t0)*p+t0,2)
    return T


def runthelec(t0, t1, td, xdn, xup, dope, ndosmx, gaussian, natom,
    _T=[], dos=sys.stdin, fout=sys.stdout, vol=None, IntegrationFunc=trapz):
    """
    Calculate thermal free energy from electronic density of states (e DOS)

    Parameters
    ----------
    t0 : float
        Low temperature limit
    t1 : float
        High temperature limit
    td : float
        Temperature increment
    xdn : float
        Minimum energy for integration
    xup : float
        Maximum energy for integration
    dope : float
        Number of electrons to dope (negative means to remove electrons, positive means add electrons)
    ndosmx : int
        Refined number of DOS points for the energy/density grid
    gaussian_grid_size : int
        Gaussian parameter to refining the grid mesh around the Fermi energy
    natom : int
        Default 1. Number of atoms in the unit cell if one wants to renomalize
        the calculated properties in the unit of per atom
    dos : file description for the DOSCAR or pymatgen dos object
        Filename for VASP DOSCAR
    outf : file description
        Output file description for the calculated properties

    Return
    ------
    Tuple of 14 float array containing
    thermal electron free energy, entropy, specific heat, M_el, seebeck_coefficients,
    effective number of charge carrier, Q_p, Q_e, constant chemical potential specific heat, temperature.
    Other quantities are for researching purpose
    """

    if hasattr(dos, 'read'):
        edn, eup, vde, dos_energies, vaspEdos = pregetdos(dos) # Line 186
    else:
        e_fermi = dos.efermi
        eup = np.max(dos.energies) - e_fermi
        edn = np.min(dos.energies) - e_fermi
        n_dos = len(dos.energies) # number of points in DOS
        vde = (eup - edn)/(n_dos-1) # change in energy per step
        dos_energies = np.linspace(edn, eup, n_dos) # linearize: sometimes rounding errors in DOSCAR
        vaspEdos = np.array(dos.get_densities())
        _eFermi = CBMtoVBM(dos_energies, vaspEdos)
        eup -= _eFermi
        edn -= _eFermi
        dos_energies -= _eFermi
    NELECTRONS, E0, dF, e, dos, Eg = getdos(xdn, xup, dope, ndosmx, gaussian, edn, eup, vde, dos_energies, vaspEdos)

    if Eg < 0.0: Eg = 0.0
    if vol == None:
        fout.write('#Bandgap= {} eV. '.format(Eg))
    else:
        fout.write('#Bandgap= {} eV at volume= {} Angstrom^3/cell. '.format(Eg,vol))

    fout.write('Fermi energy was shifted {} due to doping of {} resulting Ne={} \n'.format(dF, dope, NELECTRONS))

    # for all temperatures
    if len(_T)!=0:
      T=copy.deepcopy(_T)
    elif td>0.0:
      T = np.arange(t0,t1+td,td) # temperature
    else:
      if self.debug:
        T = T_remesh(t0,t1,td,_nT=65)
      else:
        T = T_remesh(t0,t1,td,_nT=self.nT)
    nT = len(T)
    U_el = np.zeros(nT)
    S_el = np.zeros(nT)
    C_el = np.zeros(nT) # electronic specific heat
    C_mu = np.zeros(nT) # electronic specific heat at constant chemical potential
    M_el = np.zeros(nT) # electronic chemical potential, i.e., absolute thermal electric force
    Q_el = np.zeros(nT) # total number of thermal Carrier
    Y_el = np.zeros(nT)
    Q_p = np.zeros(nT)
    Q_e = np.zeros(nT)
    W_p = np.zeros(nT)
    W_e = np.zeros(nT)
    Y_p = np.zeros(nT)
    Y_e = np.zeros(nT)
    seebeck_coefficients = np.zeros(nT)
    U_el[0] = E0

    for i in range(0,nT):
        if T[i]==0.0: continue
        Beta = 1.0e0/(T[i]*k_B)
        M_el[i], U_el[i], S_el[i], C_el[i], Q_el[i],Y_el[i], Q_p[i],Q_e[i],  C_mu[i], W_p[i], W_e[i], Y_p[i], Y_e[i] = caclf(e, dos, NELECTRONS, Beta, M_el[i-1], dF=-dF, IntegrationFunc=IntegrationFunc)
        if Q_el[i]>0.0:
            seebeck_coefficients[i] = -1.0e6*Y_el[i]/Q_el[i]/T[i]

    F_el_atom = (U_el - T * S_el - E0) / natom  # electronic free energy per atom
    S_el_atom = S_el / natom  # entropy per atom
    #dU_dT = np.gradient(U_el, td) # gradient on U_el with step size of td
    #C_el_atom = dU_dT/natom # electronic heat capacity per atom
    C_el_atom = C_el/natom # electronic heat capacity per atom

    return F_el_atom, S_el_atom, C_el_atom, M_el, seebeck_coefficients, Q_el, Q_p, Q_e, C_mu, T, W_p, W_e, Y_p, Y_e


def thelecAPI(t0, t1, td, xdn, xup, dope, ndosmx, gaussian, natom, outf, doscar):
    """
    API to calculate the thermal electronic properties from DOSCAR

    Parameters
    ----------
    t0 : float
        Low temperature limit
    t1 : float
        High temperature limit
    td : float
        Temperature increment
    xdn : float
        Minimum energy for integration
    xup : float
        Maximum energy for integration
    dope : float
        Number of electrons to dope (negative means to remove electrons, positive means add electrons)
    ndosmx : int
        Refined number of DOS points for the energy/density grid
    gaussian_grid_size : int
        Gaussian parameter to refining the grid mesh around the Fermi energy
    natom : int
        Default 1. Number of atoms in the unit cell if one wants to renomalize
        the calculated properties in the unit of per atom
    doscar : str
        Filename for VASP DOSCAR
    poscar : str
        Filename for VASP POSCAR
    vdos   : str
        Filename for Yphon phonon DOS

    Output to (printed to outf)
    ---------------------------
    The properties in the order of temperature, thermal electron free energy, entropy,
    specific heat, seebeck_coefficients, Lorenz number,
    effective number of charge carrier, Q_p, Q_e, constant chemical potential specific heat
    """

    with open(doscar, 'r') as fp:
      with open(outf, 'w') as fvib:
        F_el_atom, S_el_atom, C_el_atom, M_el, seebeck_coefficients, Q_el, Q_p, Q_e, C_mu, T, W_p, W_e, Y_p, Y_e = runthelec(t0, t1, td, xdn, xup, dope, ndosmx, gaussian, natom, dos=fp, fout=fvib)
        fvib.write('#T, F_el_atom, S_el_atom, C_el_atom, M_el, seebeck_coefficients, Lorenz_number[WΩK−2], Q_el, Q_p, Q_e, C_mu, W_p, W_e, Y_p, Y_e\n')
        for i in range(T.size):
            L = 2.443004551768e-08 #1.380649e-23/1.60217662e-19x3.14159265359xx2/3
            if Q_el[i] != 0.0: L = C_el_atom[i]/Q_el[i]*k_B
            if Q_el[i] > 1.e-16: L = C_el_atom[i]/Q_el[i]*k_B
            fvib.write('{} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.format(T[i], F_el_atom[i], S_el_atom[i], C_el_atom[i], M_el[i], seebeck_coefficients[i], L, Q_el[i], Q_p[i], Q_e[i], C_mu[i], W_p[i], W_e[i], Y_p[i], Y_e[i]))

def BMvol(V,a):
  T = V**(-1./3)
  fval = a[0]+a[1]*T
  if len(a) > 2:
    fval += a[2]*T*T
  if len(a) > 3:
    fval += a[3]*T*T*T
  if len(a) > 4:
    fval += a[4]*T*T*T*T
  return(fval)

def BMvolP(V,a):
  T = V**(-1./3)
  fval = a[1]
  if len(a) > 2:
    fval += 2*a[2]*T
  if len(a) > 3:
    fval += 3*a[3]*T*T
  if len(a) > 4:
    fval += 4*a[4]*T*T*T
  return(fval)

def BMvolB(V,a):
  T = V**(-1./3)
  fval = 2*a[2]
  if len(a) > 3:
    fval += 2*3*a[3]*T
  if len(a) > 4:
    fval += 3*4*a[4]*T*T
  return(fval)

def BMvol4(T, a, b, c, d):
  return (BMvol(T, [a,b,c,d]))

def BMvol5(T, a, b, c, d, e):
  return (BMvol(T, [a,b,c,d,e]))

def alt_curve_fit(BMfunc, x, y):
  #it is found the python curve_fit can result in numerical instability
  if 1==0: return curve_fit(BMfunc, x, y)

  #change back to linear fitting to avoid numerical instability
  xx = np.array(x)**(-1/3)
  if BMfunc.__name__=="BMvol4":
    return np.polyfit(xx, y, 3)[::-1], 0
  elif BMfunc.__name__=="BMvol5":
    return np.polyfit(xx, y, 4)[::-1], 0

def BMfitB(V, x, y, BMfunc):
  f, pcov = alt_curve_fit(BMfunc, x, y)
  p = BMvolP(V, f)
  b = BMvolB(V, f)
  P = p*V**(-4./3)*(-1./3.)
  B = b*V**(-4./3)*(-1./3.)*V**(-4./3)*(-1./3.) + p*V**(-7./3)*(-1./3.)*(-4./3.)
  return B*V, P

def BMfitP(V, x, y, BMfunc):
  f, pcov = alt_curve_fit(BMfunc, x, y)
  p = BMvolP(V, f)
  P = p*V**(-4./3)*(-1./3.)
  return P

def BMfitF(V, x, y, BMfunc):
  f, pcov = alt_curve_fit(BMfunc, x, y)
  return BMvol(V,f)

def BMsmooth(_V, _E0, _Flat, _Fel, _Slat, _Sel, BMfunc, elmode):
    E0 = BMfitF(_V, _V, _E0, BMfunc)
    if elmode==1:
        Flat = UnivariateSpline(_V, _Flat)(_V)
        Slat = UnivariateSpline(_V, _Slat)(_V)
        Fel = UnivariateSpline(_V, _Fel)(_V)
        Sel = UnivariateSpline(_V, _Sel)(_V)
    elif elmode==2:
        Flat = _Flat
        Slat = _Slat
        Fel = UnivariateSpline(_V, _Fel)(_V)
        Sel = UnivariateSpline(_V, _Sel)(_V)
    elif elmode==3:
        Flat = _Flat
        Slat = _Slat
        p1 = np.poly1d(np.polyfit(_V, _Fel, 1))
        Fel = p1(_V)
        p1 = np.poly1d(np.polyfit(_V, _Sel, 1))
        Sel = p1(_V)
    elif elmode==4:
        p1 = np.poly1d(np.polyfit(_V, _Flat, 1))
        Flat = p1(_V)
        p1 = np.poly1d(np.polyfit(_V, _Slat, 1))
        Slat = p1(_V)
        p1 = np.poly1d(np.polyfit(_V, _Fel, 1))
        Fel = p1(_V)
        p1 = np.poly1d(np.polyfit(_V, _Sel, 1))
        Sel = p1(_V)
    else:
        Flat = _Flat
        Slat = _Slat
        Fel = _Fel
        Sel = _Sel
    return E0, Flat, Fel, Slat, Sel

def CenDif(v, vol, F, N=7,kind='cubic'):
    vn = min(vol)
    vx = max(vol)
    dV = 0.001*(vx-vn)
    if kind=='cubic':
        if 1==1:
            f2 = interp1d(vol, F, kind='cubic')
        else:
            f0 = splrep(vol, F)
            return splev(v, f0, der=1)
    else:
        f2 = UnivariateSpline(vol, F)

    try:
        nx = N//2
        result = 0.0
        for i in range(0, nx):
            result += (f2(v+dV*(i+1)) - f2(v-dV*(i+1)))/(2.0*(i+1)*dV)
        return result/nx
    except:
        return None


def CenDifB(vol, F, N=7,kind='cubic'):
    vn = min(vol)
    vx = max(vol)
    xx = np.linspace(vn,vx,1000)
    if kind=='cubic':
        f2 = interp1d(vol, F, kind='cubic')
    else:
        f2 = UnivariateSpline(vol, F)
    yy = f2(xx)
    val, idx = min((val, idx) for (idx, val) in enumerate(yy))
    if idx <N//2 or idx>=len(xx)-N//2-2:
        return -1.0, -1.0, 0.

    try:
        v = brentq(CenDif, xx[idx-1], xx[idx+1], args=(vol, F, N, kind), maxiter=10000)
        ff = interp1d(vol, F)(v)
    except:
        return -1.0, -1.0, 0.

    if 1==0:
        with open ("debug", "a") as fp:
            fp.write('#v={}\n'.format(v))
            for i,x in enumerate(vol):
                fp.write('{} {}\n'.format(x, F[i]))
            fp.write('\n\n')

            fp.write('#v={}\n'.format(v))
            for i,x in enumerate(xx):
                fp.write('{} {}\n'.format(x, yy[i]))
            fp.write('\n\n')

    n = N//2
    dV = xx[1]-xx[0]
    for nx in range(n, 0, -1):
        if v-dV*(nx+1) > vn:
            if v+dV*(nx+1) < vx: break
    if nx < 1: return -1.0, -1.0, 0.

    try:
        result = 0.0
        for i in range(0, nx):
            result += (CenDif(v+dV*(i+1), vol, F, kind=kind) - CenDif(v-dV*(i+1), vol, F, kind=kind))/(2.0*(i+1)*dV)
        return result*v/(nx), v, ff
    except:
        return -1.0, -1.0, 0.

def BMDifB(vol, F, BMfunc, N=7, _T=0):
    vn = min(vol)
    vx = max(vol)
    xx = np.linspace(vn,vx,1000)
    yy = BMfitF(xx, vol, F, BMfunc)
    val, idx = min((val, idx) for (idx, val) in enumerate(yy))
    if idx <N//2 or idx>=len(xx)-N//2-2:
        return -1.0, -1.0, 0.,0.
    v = brentq(BMfitP, xx[idx-1], xx[idx+1], args=(vol, F, BMfunc), maxiter=10000)
    ff = BMfitF(v, vol, F, BMfunc)
    bb, pp = BMfitB(v, vol, F, BMfunc)
    return bb, v, ff, pp


def debye_heat_capacity(temperature, debye_T, natoms):
    """
    debye Vibrational heat capacity, C_vib(V, T).
    Eq(4) in doi.org/10.1016/j.comphy.2003.12.001

    Args:
        temperature (float): temperature in K
        volume (float)

    Returns:
        float: vibrational heat capacity in eV
    """
    y = debye_T / temperature
    factor = 3. / y ** 3
    if y < 155:
        integral = quadrature(lambda x: x ** 4 *np.exp(x)/ (np.exp(x) - 1.)**2, 0, y)
        return 3*k_B * natoms * list(integral)[0] * factor
    else:
        return k_B * natoms * 4./5.*math.pi**4 * factor


def debye_phonon(x, temperature, natoms, clat):
    return debye_heat_capacity(temperature, x, natoms) - clat


def get_debye_T_from_phonon_Cv(temperature, clat, dlat, natoms, _td=50):
    if temperature <=0: return dlat
    t0 = dlat
    d0 = debye_phonon(t0, temperature, natoms, clat)
    if d0 > 0.0: td = _td
    elif d0 <0.0: td = -_td
    else: return t0
    for i in range(999):
        if t0 < 0.1 : return 0
        t1 = t0 + td
        t1 = max(t1,0.01)
        d1 = debye_phonon(t1, temperature, natoms, clat)
        if d1*d0 < 0.0: break
        t0 = t1
        d0 = d1
        td = td + td
    return brentq(debye_phonon, t0, t1, args=(temperature, natoms, clat), maxiter=10000)


def vol_within(vol, volumes, thr=0.001):
    for i,v in enumerate(volumes):
        if (abs(vol-v) < thr*vol): return True
    return False


def vol_closest(vol, volumes, thr=1.e-6):
    for i,v in enumerate(volumes):
        if (abs(vol-v) < thr*vol): return i
    return -1


def get_static_calculations(vasp_db, tag):

    # get the energies, volumes and DOS objects by searching for the tag
    if vasp_db.collection.count_documents({'$and':[ {'metadata.tag': tag}, {'adopted': True}, \
        {'output.structure.lattice.volume': {'$exists': True}}]}) <= 5:
        static_calculations = vasp_db.collection.find({'$and':[ {'metadata.tag': tag}, \
            {'output.structure.lattice.volume': {'$exists': True} }]})
    else:
        static_calculations = vasp_db.collection.find({'$and':[ {'metadata.tag': tag}, {'adopted': True} ]})
    energies = []
    volumes = []
    dos_objs = []  # pymatgen.electronic_structure.dos.Dos objects
    structure = None  # single Structure for QHA calculation
    _energies = []
    _volumes = []
    _dos_objs = []  # pymatgen.electronic_structure.dos.Dos objects
    emin = 1.e36
    for calc in static_calculations:
        ee = calc['output']['energy']
        if ee < emin : _calc = calc
        if len(calc['metadata'])==1:
            _vol = calc['output']['structure']['lattice']['volume']
            if np.any(abs(np.array(volumes)-_vol)<_vol*1.e-5): continue
            energies.append(ee)
            volumes.append(calc['output']['structure']['lattice']['volume'])
            dos_objs.append(vasp_db.get_dos(calc['task_id']))
        else:
            _energies.append(ee)
            _volumes.append(calc['output']['structure']['lattice']['volume'])
            _dos_objs.append(vasp_db.get_dos(calc['task_id']))

    tvolumes = np.array(sorted(volumes))
    dvolumes = tvolumes[1:-1] - tvolumes[0:-2]
    dvolumes = sorted(dvolumes)
    if abs(dvolumes[-1]-dvolumes[-2]) > 0.01*dvolumes[-1]:
        #adding useful contraint calculations if not calculated statically
        if len(_volumes)!=0:
            for i, _vol in enumerate(_volumes):
                if np.any(abs(np.array(volumes)-_vol)<_vol*1.e-5): continue
                volumes.append(_vol)
                energies.append(_energies[i])
                dos_objs.append(_dos_objs[i])

    # sort everything in volume order
    # note that we are doing volume last because it is the thing we are sorting by!
    energies = sort_x_by_y(energies, volumes)
    dos_objs = sort_x_by_y(dos_objs, volumes)
    volumes = sorted(volumes)
    volumes = np.array(volumes)
    energies = np.array(energies)
    return volumes, energies, dos_objs, _calc


def finished_calc():
    tags = []
    _, dirs, _ = next(walk("."))
    for calc in  dirs:
        readme = os.path.join(calc, 'readme')
        if not os.path.exists(readme) : continue
        with open(readme) as f:
            readme = json.load(f)
            tag = readme['METADATA']['tag']
            tags.append(tag)
    return tags


class thelecMDB():
    """
    API to calculate the thermal electronic properties from the saved dos and volume dependence in MongDB database

    Parameters
    ----------
    xdn : float
        Minimum energy for integration
    xup : float
        Maximum energy for integration
    dope : float
        Number of electrons to dope (negative means to remove electrons, positive means add electrons)
    ndosmx : int
        Refined number of DOS points for the energy/density grid
    gaussian_grid_size : int
        Gaussian parameter to refining the grid mesh around the Fermi energy
    natom : int
        Default 1. Number of atoms in the unit cell if one wants to renomalize
        the calculated properties in the unit of per atom
    everyT : int
        Default 1. number of temperature points skipped from QHA analysis for debug speeding purpose
    outf : str
        Output filename for the calculated properties
    vasp_db : obejct
        object containing the information to access the MongoDB database
    metatag : str
        metadata tag to access the to be calculated the compound
    qhamode : str
        Mode for the quasiharmonic calculations, according to it the T dependence volume to be extracted
    eqamode : int
        Mode to get LTC and the equilibrium volume
    doscar : str
        Filename for VASP DOSCAR
    poscar : str
        Filename for VASP POSCAR
    vdos   : str
        Filename for Yphon phonon DOS

    Output to (printed to outf)
    ---------------------------
    The properties in the order of temperature, thermal electron free energy, entropy,
    specific heat, seebeck_coefficients, Lorenz number,
    effective number of charge carrier, Q_p, Q_e, constant chemical potential specific heat, ...,
    and the equilibrium volume extracted from MongoDB in the last column
    """

    def __init__(self, t0, t1, td, xdn, xup, dope, ndosmx, gaussian, natfactor, outf, vasp_db=None,
                noel=False, everyT=1, metatag=None, qhamode=None, eqmode=0, elmode=1, smooth=False,
                phasename=None, pyphon=False, debug=False, renew=False, fitF=False, args=None):
        from atomate.vasp.database import VaspCalcDb
        from pymatgen.core import Structure
        self.vasp_db = vasp_db
        self.t0 = t0
        self.t1 = t1
        self.td = td
        self.xdn = xdn
        self.xup = xup
        self.dope = dope
        self.ndosmx = ndosmx
        self.gaussian = gaussian
        self.natfactor = natfactor
        self.outf = outf
        self.noel = noel
        self.everyT = everyT
        self.tag = metatag
        self.qhamode = qhamode
        self.eqmode = eqmode
        self.BMfunc = BMvol4
        if eqmode==5: self.BMfunc = BMvol5
        self.elmode = elmode
        self.smooth = smooth
        self.phasename=phasename
        self.pyphon=pyphon
        self.debug=debug
        self.renew=renew
        self.refresh=args.refresh
        self.fitF=fitF
        self.k_ph_mode = args.k_ph_mode
        self.force_constant_factor = 1.0

        if self.debug:
            if self.dope==0.0: self.dope=-1.e-5

        self.local=""
        if args!=None:
            self.nT = args.nT
            self.doscar=args.doscar
            self.poscar=args.poscar
            self.vdos=args.vdos
            self.local=args.local

        if self.vasp_db==None: self.pyphon=True
        #print ("iiiii=",len(self._Yphon))


    def get_superfij(self,i, phdir):
        if vol_within(float(i['volume']), self.Vlat):
            print("\nit seems a repeated phonon calculation for", i['volume'],"so it is discared\n")
            return None
        try:
            structure = Structure.from_dict(i['unitcell'])
        except:
            print("\nit seems phonon for", i['volume'],"is not finished yet  and so it is discared\n")
            return None
        vol = 'V{:010.6f}'.format(float(i['volume']))
        voldir = os.path.join(phdir,vol)
        if not os.path.exists(voldir):
            os.mkdir(voldir)
        elif os.path.exists(os.path.join(voldir,'superfij.out')): return voldir

        poscar = structure.to(fmt="poscar")
        unitcell_l = str(poscar).split('\n')
        natom = len(structure.sites)

        supercell_matrix = i['supercell_matrix']
        supercell_structure = copy.deepcopy(structure)
        supercell_structure.make_supercell(supercell_matrix)

        sa = SpacegroupAnalyzer(supercell_structure)
        #reduced_structure = supercell_structure.get_reduced_structure(reduction_algo='niggli')
        #print ('niggli reduced structure', reduced_structure)
        #poscar = reduced_structure.to(fmt="poscar")
        primitive_unitcell_structure = sa.find_primitive()
        poscar = primitive_unitcell_structure.to(fmt="poscar")
        punitcell_l = str(poscar).split('\n')

        natoms = len(supercell_structure.sites)
        ##print(supercell_structure.sites)
        poscar = supercell_structure.to(fmt="poscar")
        supercell_l = str(poscar).split('\n')
        structure.to(filename=os.path.join(voldir,'POSCAR'))

        with open (os.path.join(voldir,'metadata.json'),'w') as out:
            mm = i['metadata']
            mm['volume'] = i['volume']
            mm['energy'] = self.energies[(list(self.volumes)).index(i['volume'])]
            out.write('{}\n'.format(mm))


        if len(self.xmlvol)!=0:
            ii = vol_closest(mm['volume'],self.xmlvol)
            with open (os.path.join(voldir,'vasprun.xml.gz'),'wb') as out:
                out.write(self.xmlgz[ii])

            with open (os.path.join(voldir,'DOSCAR.gz'),'wb') as out:
                out.write(self.dosgz[ii])


        with open (os.path.join(voldir,'OSZICAR'),'w') as out:
            out.write('   1 F= xx E0= {}\n'.format(self.energies[(list(self.volumes)).index(i['volume'])]))
        with open (os.path.join(voldir,'superfij.out'),'w') as out:
            for line in range (2,5):
                out.write('{}\n'.format(unitcell_l[line]))
            for line in range (2,5):
                out.write('{}\n'.format(supercell_l[line]))
            out.write('{} {}\n'.format(natoms, natoms//natom))
            for line in range (7,natoms+8):
                out.write('{}\n'.format(supercell_l[line]))
            force_constant_matrix = np.array(i['force_constants'])
            hessian_matrix = np.empty((natoms*3, natoms*3), dtype=float)
            """
            PI2 = 2*3.141592653589793
            THz = 1e12*PI2
            AMU = 1.66053906660e-27
            eV = 1.602176634e-19
            M = 1.e10
            THz_to_eV=THz**2 *AMU/eV/M/M
            THz_to_eV = 0.004091649655126895
            """

            for ii in range(natoms):
               for jj in range(natoms):
                    for x in range(3):
                        for y in range(3):
                            hessian_matrix[ii*3+x, jj*3+y] = -force_constant_matrix[ii,jj,x,y]

            hessian_matrix *= self.force_constant_factor
            if self.force_constant_factor!=1.0: 
                print ("\n force constant matrix has been rescaled by :", self.force_constant_factor)

            for xx in range(natoms*3):
                for yy in range(natoms*3-1):
                    out.write('{} '.format(hessian_matrix[xx,yy]))
                out.write('{}\n'.format(hessian_matrix[xx,natoms*3-1]))
        return voldir


    def get_dielecfij(self, phdir):
        #volumes = (self.vasp_db).db['phonon'].find({'metadata.tag': self.tag}, {'_id':0, 'volume':1})
        #volumes_p = [i['volume'] for i in volumes]
        volumes = (self.vasp_db).db['phonon'].find({'metadata.tag': self.tag})
        volumes_p = []
        for i in volumes:
            try:
                structure = Structure.from_dict(i['unitcell'])
                volumes_p.append(i['volume'])
            except:
                pass
        try:
            volumes_b = (self.vasp_db).db['borncharge'].find({'metadata.tag': self.tag}, {'_id':0, 'volume':1})
            volumes_b = [i['volume'] for i in volumes_b]
            has_Born =  all(elem in volumes_b for elem in volumes_p)
        except:
            has_Born = False
        if has_Born:
            nV = 0
            for v in volumes_b:
                voldir = phdir+'/'+ 'V{:010.6f}'.format(float(v))
                if not os.path.exists(voldir+'/dielecfij.out'): break
                nV += 1
            if nV==len(volumes_b) and not self.refresh: return has_Born

            for i in (self.vasp_db).db['borncharge'].find({'metadata.tag': self.tag}):
                vol = 'V{:010.6f}'.format(float(i['volume']))
                voldir = phdir+'/'+vol
                if not os.path.exists(voldir):
                   os.mkdir(voldir)
                elif os.path.exists(voldir+'/dielecfij.out') and not self.refresh: continue

                structure = Structure.from_dict(i['structure'])
                poscar = structure.to(fmt="poscar")
                poscar = str(poscar).split('\n')
                natom = len(structure.sites)
                with open (voldir+'/dielecfij.out','w') as out:
                    for line in range (2,5):
                        out.write('{}\n'.format(poscar[line]))
                    for line in range (8,natom+8):
                        out.write('{}\n'.format(poscar[line]))
                    dielectric_tensor = np.array(i['dielectric_tensor'])
                    for x in range(3):
                        for y in range(3):
                            out.write('{} '.format(dielectric_tensor[x,y]))
                        out.write('\n')
                    born_charge = np.array(i['born_charge'])
                    for ii in range(natom):
                        out.write(' ion  {}\n'.format(ii+1))
                        for x in range(3):
                            out.write('  {} '.format(x+1))
                            for y in range(3):
                                out.write('{} '.format(born_charge[ii,x,y]))
                            out.write('\n')
        return has_Born


    def get_Cij_old(self, phdir):
        try:
            volumes_c = (self.vasp_db).db['elasticity'].find({'metadata.tag': self.tag}, \
                {'_id':0, 'elastic_tensor':1, 'initial_structure':1})
        except:
            volumes_c = []

        self.Cij = []
        self.VCij = []
        for i in volumes_c:
            vol  = float(i['initial_structure']['lattice']['volume'])
            if vol in self.VCij: continue
            if vol not in self.volumes: continue
            voldir = phdir+'/'+'V{:010.6f}'.format(vol)
            if not os.path.exists(voldir): os.mkdir(voldir)

            Cij = i['elastic_tensor']['ieee_format']
            self.Cij.append(Cij)
            self.VCij.append(vol)
            with open (voldir+'/Cij.out','w') as out:
                for x in range(6):
                    for y in range(6):
                        out.write('{:10.2f} '.format(Cij[x][y]))
                    out.write('\n')

        has_Cij = len(self.Cij)>0
        if has_Cij:
            self.Cij = np.array(sort_x_by_y(self.Cij, self.VCij))
            self.VCij = sort_x_by_y(self.VCij, self.VCij)
        return has_Cij


    def Cij_by_pinv(self,ii):
        all_stresses = np.array(ii['fitting_data']['cauchy_stresses'])
        all_strains = np.array(ii['fitting_data']['strains'])
        stresses = []
        strains = []
        eq_stress = None
        for jj, strain in enumerate(all_strains):
            if (abs(strain) < 1e-10).all():
                eq_stress = np.array(copy.deepcopy(all_stresses[jj]))
            else:
                stresses.append(all_stresses[jj])
                strains.append(all_strains[jj])
        stresses = np.array(stresses)
        strains = np.array(strains)
        if eq_stress is not None:
            stresses -= eq_stress
        stresses_v = np.zeros((len(strains),6))
        strains_v = np.zeros((len(strains),6))
        for i in range(len(strains)):
            #print (stresses.shape)
            stresses_v[i,0] = stresses[i, 0, 0]
            stresses_v[i,1] = stresses[i, 1, 1]
            stresses_v[i,2] = stresses[i, 2, 2]
            stresses_v[i,3] = (stresses[i, 1, 2] + stresses[i, 2, 1])*0.5
            stresses_v[i,4] = (stresses[i, 2, 0] + stresses[i, 0, 2])*0.5
            stresses_v[i,5] = (stresses[i, 0, 1] + stresses[i, 1, 0])*0.5
            strains_v[i,0] = strains[i, 0, 0]
            strains_v[i,1] = strains[i, 1, 1]
            strains_v[i,2] = strains[i, 2, 2]
            strains_v[i,3] = (strains[i, 1, 2] + strains[i, 2, 1])
            strains_v[i,4] = (strains[i, 2, 0] + strains[i, 0, 2])
            strains_v[i,5] = (strains[i, 0, 1] + strains[i, 1, 0])
        svd_strains_v = np.linalg.pinv(strains_v.T)
        Cij = np.matmul(stresses_v.T,svd_strains_v)
        for i in range(6):
            for j in range(i):
                Cij[i,j] = 0.5*(Cij[i,j]+Cij[j,i])
                Cij[j,i] = Cij[i,j]
        """
        print(svd_strains_v)
        print(stresses_v)
        print(strains_v)
        print(Cij)
        print(ii['elastic_tensor']['ieee_format'])
        sys.exit()
        """
        return Cij


    def vol_in_volumes(self, vol):
        for v in self.volumes:
            ncell = int(vol/v+1.e-12)
            if abs(vol-ncell*v)<1.e-12: return True, v
        return False, 0


    def get_Cij(self, phdir, pinv=True):
        try:
            volumes_c = (self.vasp_db).db['elasticity'].find({'metadata.tag': self.tag}, \
                {'_id':0, 'elastic_tensor':1, 'initial_structure':1, 'fitting_data':1})
        except:
            volumes_c = []

        self.Cij = []
        self.VCij = []
        self.Poisson_Ratio = []

        if not os.path.exists(phdir): os.mkdir(phdir)
        for i in volumes_c:
            vol  = float(i['initial_structure']['lattice']['volume'])
            if vol in self.VCij: continue
            #if vol not in self.volumes: continue
            volin, vol = self.vol_in_volumes(vol)
            if not volin: continue
            voldir = phdir+'/'+'V{:010.6f}'.format(vol)
            if not os.path.exists(voldir): os.mkdir(voldir)

            ngroup = self.space_group_number
            #if ngroup>=195 and ngroup<=230:
            if pinv:
                Cij = self.Cij_by_pinv(i)
            else:
                Cij = i['elastic_tensor']['ieee_format']
            """
            if ngroup>=195 and ngroup<=230: # for cubic system
                C11 = (Cij[0,0] + Cij[1,1] + Cij[2,2])/3.0
                C12 = (Cij[0,1] + Cij[1,2] + Cij[2,0])/3.0
                C44 = (Cij[3,3] + Cij[4,4] + Cij[5,5])/3.0
                for i in range(0,3):
                    Cij[i,i] = C11
                    Cij[i+3,i+3] = C44
                    for j in range(0,3):
                        if j==i: continue
                        Cij[i,j] = C12
            """

            self.Cij.append(Cij)
            self.VCij.append(vol)
            _,_,_,Poisson_Ratio=self.Cij_to_Moduli(np.stack(Cij))
            self.Poisson_Ratio.append(Poisson_Ratio)
            with open (voldir+'/Cij.out','w') as out:
                for x in range(6):
                    for y in range(6):
                        out.write('{:10.2f} '.format(Cij[x][y]))
                    out.write('\n')

        has_Cij = len(self.Cij)>0
        if has_Cij:
            self.Cij = np.array(sort_x_by_y(self.Cij, self.VCij))
            self.Poisson_Ratio = np.array(sort_x_by_y(self.Poisson_Ratio, self.VCij))
            self.VCij = sort_x_by_y(self.VCij, self.VCij)
        return has_Cij


    #For volume dependent thermoelectric calculations based on the output from BoltrzTraP2 code
    #btp2_dir is parent path containing thermoelectric results at several volumes
    def get_btp2(self, btp2_dir):
        volumes_uniform_tau = []
        volumes_uniform_lambda = []
        uniform_tau = []
        uniform_lambda = []
        yphondir = os.path.join(btp2_dir,"Yphon")
        if not os.path.exists(yphondir): yphondir = btp2_dir
        _, static_calculations, _ = next(walk(yphondir))
        for calc in static_calculations:
            poscar = os.path.join(yphondir, calc, 'POSCAR')
            if not os.path.exists(poscar) : continue
            tmp = os.path.join(yphondir, calc, 'interpolation.dope.trace.uniform_tau')
            if os.path.exists(tmp) :
                print ("Handling data in ", tmp)
                structure = Structure.from_file(poscar)
                vol = structure.volume
                if vol not in volumes_uniform_tau:
                    volumes_uniform_tau.append(vol)
                    uniform_tau.append(np.genfromtxt(tmp))
            tmp = os.path.join(yphondir, calc, 'interpolation.dope.trace.uniform_lambda')
            if os.path.exists(tmp) :
                print ("Handling data in ", tmp)
                structure = Structure.from_file(poscar)
                vol = structure.volume
                if vol not in volumes_uniform_lambda:
                    volumes_uniform_lambda.append(vol)
                    uniform_lambda.append(np.genfromtxt(tmp))

        has_btp2 = False
        from dfttk.utils import sort_x_by_y
        if len(volumes_uniform_tau)>0:
            has_btp2 = True
            self.uniform_tau = sort_x_by_y(uniform_tau, volumes_uniform_tau)
            self.volumes_uniform_tau = sort_x_by_y(volumes_uniform_tau,volumes_uniform_tau)
            print ("found volumes_uniform_tau volumes from static calculations:", self.volumes_uniform_tau)
        if len(volumes_uniform_lambda)>0:
            has_btp2 = True
            self.uniform_lambda = sort_x_by_y(uniform_lambda, volumes_uniform_lambda)
            self.volumes_uniform_lambda = sort_x_by_y(volumes_uniform_lambda,volumes_uniform_lambda)
            print ("found volumes_uniform_lambda volumes from static calculations:", self.volumes_uniform_lambda)
        return has_btp2


    def toYphon(self, _T=None, for_plot=False):
        if self.vasp_db==None or self.local!="":
            self.toYphon_without_DB(_T=_T, for_plot=for_plot)
            return

        self.Vlat = []
        self.Flat = []
        self.Clat = []
        self.Slat = []
        self.quality = []
        if _T is not None:
            self.T_vib = copy.deepcopy(_T)
        elif self.debug:
            self.T_vib = T_remesh(self.t0, self.t1, self.td, _nT=65)
        else:
            self.T_vib = T_remesh(self.t0, self.t1, self.td, _nT=self.nT)

        print ("extract the superfij.out used by Yphon ...")
        phdir = os.path.join(self.phasename,'Yphon')
        if not os.path.exists(phdir):
            os.mkdir(phdir)

        has_Born = self.get_dielecfij(phdir)
        for i in (self.vasp_db).db['phonon'].find({'metadata.tag': self.tag}):
            try:
                self.force_constant_factor = i['force_constant_factor']
            except:
                if self.static_vasp_version[0:1] >= '6':
                    self.force_constant_factor = 0.004091649655126895

            if i['volume'] not in self.volumes: continue
            voldir = self.get_superfij(i, phdir)
            if voldir is None: continue

            self.Vlat.append(float(i['volume']))
            cwd = os.getcwd()
            os.chdir( voldir )

            if os.path.exists('vdos.out'):
                updatevdos = False
                vdostime = os.path.getmtime('vdos.out')
                if os.path.exists('superfij.out'):
                    updatevdos = os.path.getmtime('superfij.out')>vdostime
                if os.path.exists('dielecfij.out'):
                    updatevdos = os.path.getmtime('dielecfij.out')>vdostime
            else: updatevdos = True

            #if not os.path.exists('vdos.out') or self.refresh:
            if updatevdos:
                _nqwave = ""
                if self.debug:
                    _nqwave = "-nqwave "+ str(1.e4)
                #md = "Yphon -tranI 2 -DebCut 0.5 " +_nqwave+ " <superfij.out"
                cmd = "Yphon -DebCut 0.5 " +_nqwave+ " <superfij.out"
                if has_Born: cmd += " -Born dielecfij.out"
                print(cmd, " at ", voldir)
                output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    universal_newlines=True)

            if len(self.Flat)==0:
                print ("Calling yphon to get f_vib, s_vib, cv_vib at ", phdir)
            with open("vdos.out", "r") as fp:
                f_vib, U_ph, s_vib, cv_vib, C_ph_n, Sound_ph, Sound_nn, N_ph, NN_ph, debyeT, quality, natoms \
                    = ywpyphon.vibrational_contributions(self.T_vib, dos_input=fp, energyunit='eV')
                self.quality.append(quality)

            self.Flat.append(f_vib)
            self.Slat.append(s_vib)
            self.Clat.append(cv_vib)
            self.quality.append(quality)
            os.chdir( cwd )

        if len(self.Vlat)<=0:
            print("\nFATAL ERROR! cannot find required data from phonon collection for metadata tag:", self.tag,"\n")
            raise ValueError()
            #sys.exit()
        self.Slat = np.array(sort_x_by_y(self.Slat, self.Vlat))
        self.Clat = np.array(sort_x_by_y(self.Clat, self.Vlat))
        self.Flat = np.array(sort_x_by_y(self.Flat, self.Vlat))
        self.Vlat = np.array(sort_x_by_y(self.Vlat, self.Vlat))
        self.Dlat = np.full((len(self.Vlat)), 400.)
        if for_plot: return
        self.volT = np.zeros(len(self.T_vib))
        self.GibT = np.zeros(len(self.T_vib))


    def toYphon_without_DB(self, _T=None, for_plot=False):
        self.Vlat = []
        self.Flat = []
        self.Clat = []
        self.Slat = []
        self.quality = []
        if _T is not None:
            self.T_vib = copy.deepcopy(_T)
        elif self.debug:
            self.T_vib = T_remesh(self.t0, self.t1, self.td, _nT=65)
        else:
            self.T_vib = T_remesh(self.t0, self.t1, self.td, _nT=self.nT)

        print ("extract the superfij.out used by Yphon ...")
        for i, vol in enumerate(self.volumes):
            if self.local!="":
                dir = self.dirs[i]
            else:
                dir = os.path.join(self.phasename,'Yphon','V{:010.6f}'.format(vol))
            cwd = os.getcwd()
            if not os.path.exists(dir): continue

            os.chdir( dir)
            if not os.path.exists('vdos.out'):
                _nqwave = ""
                if self.debug:
                    _nqwave = "-nqwave "+ str(1.e4)
                cmd = "Yphon -tranI 2 -DebCut 0.5 " +_nqwave+ " <superfij.out"
                print(cmd, " at ", dir)
                output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    universal_newlines=True)

            with open("vdos.out", "r") as fp:
                f_vib, U_ph, s_vib, cv_vib, C_ph_n, Sound_ph, Sound_nn, N_ph, NN_ph, debyeT, quality, natoms \
                    = ywpyphon.vibrational_contributions(self.T_vib, dos_input=fp, energyunit='eV')
                self.Vlat.append(float(self.volumes[i]))
                self.quality.append(quality)

            self.Flat.append(f_vib)
            self.Slat.append(s_vib)
            self.Clat.append(cv_vib)
            self.quality.append(quality)
            os.chdir( cwd )

        if len(self.Vlat)<=0:
            print("\nFATAL ERROR! cannot find required data from phonon collection for metadata tag:", self.tag,"\n")
            raise ValueError()
        self.Slat = np.array(sort_x_by_y(self.Slat, self.Vlat))
        self.Clat = np.array(sort_x_by_y(self.Clat, self.Vlat))
        self.Flat = np.array(sort_x_by_y(self.Flat, self.Vlat))
        self.Vlat = np.array(sort_x_by_y(self.Vlat, self.Vlat))
        self.Dlat = np.full((len(self.Vlat)), 400.)
        if for_plot: return
        self.volT = np.zeros(len(self.T_vib))
        self.GibT = np.zeros(len(self.T_vib))


    def check_vol(self):
        print ("\nChecking compatibility between qha/Yphon data and static calculation:\n")
        print (self.Vlat)
        for i,v in enumerate (self.Vlat):
            fvol = False
            for j,vol in enumerate(self.volumes):
                if abs(vol-v)<1.e-8:
                    print (v, self.energies[j])
                    fvol = True
            if not fvol:
                print("\nWarning! Not found v=",v,"\n")
        if len(self.Vlat)!=len(self.volumes):
            warnings.warn("\nWarning! The static/qha calculations are not inconsistent! Let me see if I can resolve it\n")

        _volumes = list(self.volumes)
        _v = []
        _e = []
        _d = []
        for i, vol in enumerate(list(self.volumes)):
            if vol not in self.Vlat:
                print ("data in static calculation with volume=", vol, "is discarded")
                continue
            _v.append(vol)
            _e.append(self.energies[i])
            _d.append(self.dos_objs[i])
        if len(_v)<len(self.volumes) :
             self.volumes = np.array(_v)
             self.energies = np.array(_e)
             self.dos_objs = _d

        good = True
        if len(self.Vlat)==len(self.volumes) and len(self.Vlat)>=5:
            val, idx = min((val, idx) for (idx, val) in enumerate(self.energies))
            if idx!=0 and idx!=len(self.energies)-1: print("\nOK, I found some data that I can try\n")
            else: good = False

        if not good:
            print("static volume:", self.volumes)
            print("phonon volume:", self.Vlat)
            print("\nFATAL ERROR! It appears that the calculations are not all done!\n")
        return good
            #raise ValueError("\nFATAL ERROR! It appears that the calculations are not all done!\n")


    # get the energies, volumes and DOS objects by searching for the tag
    def find_static_calculations(self):
        self.volumes, self.energies, self.dos_objs, _calc \
            = get_static_calculations(self.vasp_db, self.tag)
        structure = Structure.from_dict(_calc['output']['structure'])
        self.structure = structure
        print(structure)
        print ("\n")
        reduced_structure = structure.get_reduced_structure(reduction_algo='niggli')
        print ('niggli reduced structure', reduced_structure)
        print ("\n")
        self.formula_pretty = structure.composition.reduced_formula
        try:
            formula2composition(formula_pretty)
        except:
            self.formula_pretty = reduced_formula(reduced_structure.composition.alphabetical_formula)

        self.natoms = len(structure.sites)
        sa = SpacegroupAnalyzer(structure)

        potsoc = get_used_pot(_calc)
        self.space_group_number = sa.get_space_group_number()
        self.phase = sa.get_space_group_symbol().replace('/','.')+'_'+str(sa.get_space_group_number())+potsoc

        key_comments ={}
        tmp = _calc['input']['pseudo_potential']
        tmp['functional'] = potsoc
        key_comments['pseudo_potential'] = tmp
        key_comments['ENCUT'] = _calc['input']['incar']['ENCUT']
        key_comments['NEDOS'] = _calc['input']['incar']['NEDOS']
        try:
            key_comments['LSORBIT'] = _calc['input']['incar']['LSORBIT']
        except:
            key_comments['LSORBIT'] = False

        pot = _calc['orig_inputs']['kpoints']
        kpoints = {}
        kpoints['generation_style'] = pot['generation_style']
        kpoints['kpoints'] = pot['kpoints'][0]
        key_comments['kpoints'] = kpoints
        key_comments['bandgap'] = _calc['output']['bandgap']
        self.key_comments = key_comments

        print ("found volumes from static calculations:", self.volumes)

        self.xmlvol = []
        self.xmlgz = []
        self.dosgz = []
        xml = list(self.vasp_db.db['xmlgz'].find({'$and':[ {'metadata.tag': self.tag}, { 'vasprun_xml_gz': { '$exists': True } }]}))
        for x in xml:
            self.xmlgz.append(pickle.loads(x['vasprun_xml_gz']))
            self.dosgz.append(pickle.loads(x['DOSCAR_gz']))
            self.xmlvol.append(x['volume'])
            print ("found ", 'vasprun.xml.gz', 'DOSCAR.gz', "at", x['volume'])

        if self.phasename is None: self.phasename = self.formula_pretty+'_'+self.phase
        if not os.path.exists(self.phasename):
            os.mkdir(self.phasename)

        vasp_version = list(self.vasp_db.db['tasks'].find({'$and':[ {'metadata.tag': self.tag}, { 'calcs_reversed': { '$exists': True } }]}))
        self.static_vasp_version = None
        for i in vasp_version:
            v = i['calcs_reversed'][0]['vasp_version']
            if self.static_vasp_version is None: self.static_vasp_version = v
            elif v[0:1]!=self.static_vasp_version[0:1]:
                print("\n***********FETAL messing up calculation! please remove:", self.tag, "\n")
        if self.static_vasp_version is not None:
            print("\nvasp version for the static calculation is:", self.static_vasp_version, " for ", self.tag, "\n")


    # get the energies, volumes and DOS objects by searching for the tag
    def find_static_calculations_local(self):
        yphondir = os.path.join(self.local,"Yphon")
        if not os.path.exists(yphondir): yphondir = self.local

        _, static_calculations, _ = next(walk(yphondir))

        energies = []
        dirs = []
        volumes = []
        lattices = []
        _matrixs = []
        dos_objs = []  # pymatgen.electronic_structure.dos.Dos objects
        _structure = None  # single Structure for QHA calculation
        self.structure = None
        for calc in static_calculations:
            poscar = os.path.join(yphondir, calc, 'POSCAR')
            if not os.path.exists(poscar) : continue
            oszicar = os.path.join(yphondir, calc, 'OSZICAR')
            if not os.path.exists(oszicar) : continue
            print ("Handling data in ", os.path.join(yphondir,calc))
            structure = Structure.from_file(poscar)
            if self.structure == None:
                self.structure = structure
            vol = structure.volume
            sss = (structure.lattice.matrix).tolist()
            lattices.append(sss)

            if vol in volumes:
                print ("WARNING: skipped volume =", vol)
                continue
            volumes.append(vol)
            dirs.append(os.path.join(yphondir,calc))
            with open(oszicar,"r") as fp:
                lines = fp.readlines()
                for line in lines:
                    dat = [s for s in line.split() if s!=""]
                    if len(dat) < 5: continue
                    if dat[1]!="F=" or dat[3]!="E0=": continue
                    energies.append(float(dat[4]))
                    break

            # get a Structure. We only need one for the masses and number of atoms in the unit cell.
            if _structure is None:
                _structure = structure
                print(structure)
                print ("\n")
                reduced_structure = structure.get_reduced_structure(reduction_algo='niggli')
                print ('niggli reduced structure', reduced_structure)
                print ("\n")
                self.formula_pretty = structure.composition.reduced_formula
                try:
                    formula2composition(formula_pretty)
                except:
                    self.formula_pretty = reduced_formula(reduced_structure.composition.alphabetical_formula)
                self.natoms = len(structure.sites)
                sa = SpacegroupAnalyzer(structure)
                self.phase = sa.get_space_group_symbol().replace('/','.')+'_'+str(sa.get_space_group_number())
                if self.phasename==None:
                    self.phasename = self.formula_pretty+'_'+self.phase
                self.space_group_number = sa.get_space_group_number()

                key_comments ={}
                key_comments['comments'] = 'local calculations'
                self.key_comments = key_comments

            if os.path.exists(os.path.join(yphondir, calc, 'DOSCAR.gz')) : 
                dos_objs.append(os.path.join(yphondir, calc, 'DOSCAR.gz'))
            elif os.path.exists(os.path.join(yphondir, calc, 'DOSCAR')) : 
                dos_objs.append(os.path.join(yphondir, calc, 'DOSCAR'))
            else:
                continue

        # sort everything in volume order
        # note that we are doing volume last because it is the thing we are sorting by!

        from dfttk.utils import sort_x_by_y
        self.energies = sort_x_by_y(energies, volumes)
        self.dos_objs = sort_x_by_y(dos_objs, volumes)
        self.dirs = sort_x_by_y(dirs,volumes)
        self.volumes = sort_x_by_y(volumes,volumes)
        self.key_comments['E-V'] = {'lattices':sort_x_by_y(lattices, volumes),
            'volumes':self.volumes, 'energies':self.energies,
            'natoms':self.natoms}
        print ("found volumes from static calculations:", self.volumes)

        if self.phasename is None: self.phasename = self.formula_pretty+'_'+self.phase
        if not os.path.exists(self.phasename):
            os.mkdir(self.phasename)



    # get the energies, volumes and DOS objects by searching for the tag
    def find_static_calculations_without_DB(self):
        with open(self.phasename+'/readme') as f:
            readme = json.load(f)
            self.key_comments = readme
            self.natoms = readme['E-V']['natoms']
            self.volumes = np.array(readme['E-V']['volumes'])
            self.energies = np.array(readme['E-V']['energies'])
        self.dos_objs = []  # pymatgen.electronic_structure.dos.Dos objects
        self.structure = None
        for vol in self.volumes:
            dir = os.path.join(self.phasename,'Yphon','V{:010.6f}'.format(vol))
            if self.structure == None:
                poscar = os.path.join(dir, 'POSCAR')
                if not os.path.exists(poscar) : continue
                self.structure = Structure.from_file(poscar)
            if os.path.exists(os.path.join(dir,'DOSCAR.gz')):
                self.dos_objs.append(os.path.join(dir,'DOSCAR.gz'))
            elif os.path.exists(os.path.join(dir,'DOSCAR')):
                self.dos_objs.append(os.path.join(dir,'DOSCAR'))


    def get_lattice_conductivity(self, beta, blat, vol, clat, theta_D=None):   
        reduced_structure = self.structure.get_reduced_structure(reduction_algo='niggli')
        formula_pretty = self.structure.composition.reduced_formula
        try:
            formula2composition(formula_pretty)
        except:
            formula_pretty = reduced_formula(reduced_structure.composition.alphabetical_formula)
        gamma = beta*blat*vol/clat
        el, com = formula2composition(self.formula_pretty, normalize=True)
        M_average = 0.0
        for i, e in enumerate(el):
            M_average += com[i]*MM_of_Elements[e]
        natoms_prim = len(reduced_structure.sites)
        
        if theta_D is None:
            nT = len(self.T)
            for j in range(nT):
                clat = interp1d(self.volumes, self.Clat[:,j])(self.volT[j])
                d1 = get_debye_T_from_phonon_Cv(self.T[j], clat, 400., self.natoms)
                if d1 > self.T[j] :
                    d0 = d1 
                    t0 = self.T[j]
                else: 
                    t1 = self.T[j]
                    dt = t1-t0
                    d00 = abs(d0-t0)
                    d11 = abs(d1-t1)
                    f1 = d00/(d00+d11)
                    theta_D = (1-f1)*d0+f1*d1
                    break

        gamma = beta*blat*vol/clat
        fac0 = 0.849*3*4.**(1./3)
        fac1 = 20*np.pi**2*(1-0.514/gamma+0.228/gamma/gamma)
        k_B = physical_constants["Boltzmann constant"][0]
        hbar = scipy_constants.hbar
        p_mass = physical_constants["atomic mass constant"][0]
        theta_a = natoms_prim**(-1./3)*theta_D
        fac2 = (k_B*theta_a/hbar/gamma)**2*k_B*M_average*p_mass*theta_a* \
            (natoms_prim/self.natoms)**1./3*1.e-10/hbar
        return fac0/fac1*fac2, gamma, theta_D


    def get_lattice_conductivity_1(self,vol):   
        reduced_structure = self.structure.get_reduced_structure(reduction_algo='niggli')
        formula_pretty = self.structure.composition.reduced_formula
        try:
            formula2composition(formula_pretty)
        except:
            formula_pretty = reduced_formula(reduced_structure.composition.alphabetical_formula)
        el, com = formula2composition(self.formula_pretty, normalize=True)
        M_average = 0.0
        for i, e in enumerate(el):
            M_average += com[i]*MM_of_Elements[e]
        natoms_prim = len(reduced_structure.sites)
        
        debyeT = np.zeros(len(self.volumes))
        nT = len(self.T)
        for j in range(nT):
            tmp1 = np.array([get_debye_T_from_phonon_Cv(self.T[j], self.Clat[i,j], 400., self.natoms) \
                for i in range(len(self.volumes))])
            f1 = interp1d(self.volumes, tmp1)
            try:
                d1 = f1(self.volT[j])
            except:
                print("********WARNING: volume=", self.volT[j], "at T=", self.T[i], "is out of range of", self.volumes)
            print("d1=",d1)
            if d1 > self.T[j] :
                d0 = d1 
                t0 = self.T[j]
                tmp0 = copy.deepcopy(tmp1)
            else: 
                t1 = self.T[j]
                dt = t1-t0
                d00 = abs(d0-t0)
                d11 = abs(d1-t1)
                f1 = d00/(d00+d11)
                theta_D = (1-f1)*d0+f1*d1
                d00 = abs(theta_D-t0)
                d11 = abs(theta_D-t1)
                f1 = d00/(d00+d11)
                debyeT = (1-f1)*tmp0+f1*tmp1
                vol = (1-f1)*self.volT[j-1]+f1*self.volT[j]
                break

        logvol = np.log(self.volumes)
        logdeb = np.log(debyeT)
        gamma = -CenDif(np.log(vol), logvol, logdeb, kind='linear')
        fac0 = 0.849*3*4.**(1./3)
        fac1 = 20*np.pi**2*(1-0.514/gamma+0.228/gamma/gamma)
        k_B = physical_constants["Boltzmann constant"][0]
        hbar = scipy_constants.hbar
        p_mass = physical_constants["atomic mass constant"][0]
        theta_a = natoms_prim**(-1./3)*theta_D
        fac2 = (k_B*theta_a/hbar/gamma)**2*k_B*M_average*p_mass*theta_a* \
            np.power(natoms_prim/self.natoms, 1/3)*1.e-10/hbar
        return fac0/fac1*fac2, gamma


    def get_lattice_conductivity_0(self, _i, T, vol):   
        reduced_structure = self.structure.get_reduced_structure(reduction_algo='niggli')
        formula_pretty = self.structure.composition.reduced_formula
        try:
            formula2composition(formula_pretty)
        except:
            formula_pretty = reduced_formula(reduced_structure.composition.alphabetical_formula)
        el, com = formula2composition(self.formula_pretty, normalize=True)
        M_average = 0.0
        for i, e in enumerate(el):
            M_average += com[i]*MM_of_Elements[e]
        natoms_prim = len(reduced_structure.sites)
        debyeT = np.zeros(len(self.Clat[:,_i]))
        for i, c in enumerate(self.Clat[:,_i]):
            debyeT[i] = get_debye_T_from_phonon_Cv(T, c, 400., self.natoms)
        f1 = interp1d(self.volumes, debyeT)
        theta_D = f1(vol)
        logvol = np.log(self.volumes)
        logdeb = np.log(debyeT)
        gamma = -CenDif(np.log(vol), logvol, logdeb, kind='linear')
        fac0 = 0.849*3*4.**(1./3)
        fac1 = 20*np.pi**2*(1-0.514/gamma+0.228/gamma/gamma)
        k_B = physical_constants["Boltzmann constant"][0]
        hbar = scipy_constants.hbar
        p_mass = physical_constants["atomic mass constant"][0]
        theta_a = natoms_prim**(-1./3)*theta_D
        fac2 = (k_B*theta_a/hbar/gamma)**2*k_B*M_average*p_mass*theta_a* \
            np.power(natoms_prim/self.natoms, 1/3)*1.e-10/hbar
        return fac0/fac1*fac2, gamma


    def get_thermo_lattice_conductivity(self, beta, blat, vol, clat, theta_D):
        reduced_structure = self.structure.get_reduced_structure(reduction_algo='niggli')
        formula_pretty = self.structure.composition.reduced_formula
        try:
            formula2composition(formula_pretty)
        except:
            formula_pretty = reduced_formula(reduced_structure.composition.alphabetical_formula)
        el, com = formula2composition(self.formula_pretty, normalize=True)
        M_average = 0.0
        for i, e in enumerate(el):
            M_average += com[i]*MM_of_Elements[e]
        natoms_prim = len(reduced_structure.sites)
        gamma = beta*blat*vol/clat
        fac0 = 0.849*3*4.**(1./3)
        fac1 = 20*np.pi**2*(1-0.514/gamma+0.228/gamma/gamma)
        k_B = physical_constants["Boltzmann constant"][0]
        hbar = scipy_constants.hbar
        p_mass = physical_constants["atomic mass constant"][0]
        theta_a = natoms_prim**(-1./3)*theta_D
        fac2 = (k_B*theta_a/hbar/gamma)**2*k_B*M_average*p_mass*theta_a* \
            np.power(natoms_prim/self.natoms, 1/3)*1.e-10/hbar
        return fac0/fac1*fac2, gamma

    def calc_thermal_electron(self):
        t0 = min(self.T)
        t1 = max(self.T)
        td = (t1-t0)/(len(self.T)-1)
        if self.td < 0: td = self.td
        #print("xxxxx=", t0,t1,td)
        #theall = np.empty([len(prp_T), int((t1-t0)/td+1.5), len(self.volumes)])
        self.theall = np.empty([14, len(self.T), len(self.volumes)])
        for i,dos in enumerate(self.dos_objs):
            #print ("processing dos object at volume: ", self.volumes[i], " with nT =", len(self.T))
            if isinstance(dos, str):
                if dos.endswith(".gz"):
                    _dos = gzip.open (dos,'rt')
                else:
                    _dos = open (dos,'r')
                prp_vol = runthelec(t0, t1, td, self.xdn, self.xup, self.dope, self.ndosmx,
                    self.gaussian, self.natfactor, dos=_dos, _T=self.T, fout=sys.stdout, vol=self.volumes[i])
            else:
                prp_vol = runthelec(t0, t1, td, self.xdn, self.xup, self.dope, self.ndosmx,
                    self.gaussian, self.natfactor, dos=dos, _T=self.T, fout=sys.stdout, vol=self.volumes[i])
            """
            if 1==1:
                iFunc = trapz
            else:
                iFunc = simps
            F_el_atom, S_el_atom, C_el_atom, M_el, seebeck_coefficients, Q_el, Q_p, Q_e, C_mu, T, W_p, W_e, Y_p, Y_e = runthelec(t0, t1, td, xdn, xup, dope, ndosmx, gaussian, natfactor, dos=dos, fout=fvib, vol=volumes[i], IntegrationFunc=iFunc)
            Tuple of 14 float array containing
            thermal electron free energy, entropy, specific heat, M_el, seebeck_coefficients,
            effective number of charge carrier, Q_p, Q_e, constant chemical potential specific heat, temperature.
            Other quantities are for researching purpose
            """
            self.theall[:,:,i] = np.array( prp_vol ) # convert Tuple into array for the convenience of quasistatic interpolation

            """
            """
            if self.vasp_db==None: continue
            phdir = os.path.join(self.phasename,'Yphon')
            if not os.path.exists(phdir): os.mkdir(Yphon)
            vol = 'V{:010.6f}'.format(self.volumes[i])
            voldir = os.path.join(phdir, vol)
            if not os.path.exists(voldir):
                os.mkdir(voldir)
            doscar = os.path.join(voldir,'DOSCAR.gz')
            if not os.path.exists(doscar):
                with gzip.open (doscar,'wt') as out:
                    for j in range(5):
                        out.write('   {}\n'.format(j))
                    eup = np.max(dos.energies)
                    edn = np.min(dos.energies)
                    out.write('{:>16.8}{:>16.8}{:5}{:>16.8}{:>16.8}\n'.format(eup,edn,len(dos.energies),dos.efermi,1.0))
                    vaspEdos = np.array(dos.get_densities())
                    ados = cumtrapz(vaspEdos, dos.energies, initial=0.0)
                    for j,d in enumerate(dos.energies):
                        out.write('{:>11.3f} {:>11.4e} {:>11.4e}\n'.format(d,vaspEdos[j],ados[j]))


    def find_vibrational(self):
        if self.pyphon:
            self.toYphon()
            return
        self.from_phonon_collection = False
        if self.qhamode=="debye":
            self.qha_items = self.vasp_db.db['qha'].find({'metadata.tag': self.tag})
        elif self.qhamode=="phonon":
            self.qha_items = self.vasp_db.db['qha_phonon'].find({'metadata.tag': self.tag})
        else:
            try:
                self.qhamode='phonon'
                self.qha_items = self.vasp_db.db['qha_phonon'].find({'metadata.tag': self.tag})
            except:
                self.qhamode='debye'
                self.qha_items = self.vasp_db.db['qha'].find({'metadata.tag': self.tag})
        # check compatibility with vasp6
        if self.qhamode=='phonon':
            for i in (self.vasp_db).db['phonon'].find({'metadata.tag': self.tag}):
                try:
                    self.force_constant_factor = i['force_constant_factor']
                except:
                    if self.static_vasp_version[0:1] >= '6':
                        print("\n**************FETAL ERROR! force constant not compatible for :", self.tag, "by default phonopy with vasp6\n")
                        

        try:
            sys.stdout.write("\nTrying to get quasiharmonic mode : {}... \n".format(self.qhamode))
            self.T_vib = self.qha_items[0][self.qhamode]['temperatures'][::self.everyT]
            #print("xxxx=",self.T_vib)
        except:
            try:
                self.qha_items = self.vasp_db.db['qha_phonon'].find({'metadata': self.tag})
                self.T_vib = self.qha_items[0][self.qhamode]['temperatures'][::self.everyT]
            except:
                try:
                    self.qhamode = 'phonon'
                    self.qha_items = self.vasp_db.db[self.qhamode].find({'metadata.tag': self.tag})
                    self.T_vib = self.qha_items[0]['temperatures'][::self.everyT]
                    self.from_phonon_collection = True
                except:
                    warnings.warn ("\nWARNING! I cannot find required data from qha_phonon, am asking help from Yphon!\n")
                    self.pyphon = True
                    self.toYphon()


    def get_qha(self):
        if self.from_phonon_collection:
            _Vlat = []
            _Slat = []
            _Clat = []
            _Flat = []

            for i in self.qha_items:
                _Vlat.append(i['volume'])
                _Slat.append(i['S_vib'][::self.everyT])
                _Clat.append(i['CV_vib'][::self.everyT])
                _Flat.append(i['F_vib'][::self.everyT])
            self.volT = np.zeros(len(self.T_vib))
            self.GibT = np.zeros(len(self.T_vib))
            _Dlat = np.full((len(_Vlat)), 400.)
        elif self.qhamode=='debye':
            _Vlat = self.qha_items[0][self.qhamode]['volumes']
            _Slat = self.qha_items[0][self.qhamode]['entropies']
            _Clat = self.qha_items[0][self.qhamode]['heat_capacities']
            _Flat = self.qha_items[0][self.qhamode]['helmholtz_energies']
            self.volT = self.qha_items[0][self.qhamode]['optimum_volumes'][::self.everyT]
            self.GibT = self.qha_items[0][self.qhamode]['gibbs_free_energy'][::self.everyT]
            _Dlat = self.qha_items[0]['debye']['debye_temperatures']
        else:
            _Vlat = self.qha_items[0][self.qhamode]['volumes']
            _Slat = self.qha_items[0][self.qhamode]['entropies']
            _Clat = self.qha_items[0][self.qhamode]['heat_capacities']
            _Flat = self.qha_items[0][self.qhamode]['helmholtz_energies']
            self.volT = self.qha_items[0][self.qhamode]['optimum_volumes'][::self.everyT]
            self.GibT = self.qha_items[0][self.qhamode]['gibbs_free_energy'][::self.everyT]
            _Dlat = self.qha_items[0]['debye']['debye_temperatures']
        Vlat = []
        Slat = []
        Clat = []
        Flat = []
        Dlat = []
        for i, vol in enumerate(_Vlat):
            if vol_within(vol, Vlat): continue
            #if vol in Vlat: continue
            if vol not in self.volumes: continue
            Vlat.append(vol)
            Slat.append(_Slat[i])
            Clat.append(_Clat[i])
            Flat.append(_Flat[i])
            Dlat.append(_Dlat[i])

        Slat = sort_x_by_y(Slat, Vlat)
        Clat = sort_x_by_y(Clat, Vlat)
        Flat = sort_x_by_y(Flat, Vlat)
        Dlat = sort_x_by_y(Dlat, Vlat)
        Vlat = sort_x_by_y(Vlat, Vlat)
        self.Slat = np.array(Slat)[:,::self.everyT]
        self.Clat = np.array(Clat)[:,::self.everyT]
        self.Flat = np.array(Flat)[:,::self.everyT]
        self.Dlat = np.array(Dlat)
        self.Vlat = np.array(Vlat)
        if self.td < 0:
            s = []
            f = []
            c = []
            for i in range(len(self.volumes)):
                s.append(interp1d(self.T_vib, self.Slat[i,:])(self.T))
                f.append(interp1d(self.T_vib, self.Flat[i,:])(self.T))
                c.append(interp1d(self.T_vib, self.Clat[i,:])(self.T))
            self.Slat = np.array(s)
            self.Clat = np.array(c)
            self.Flat = np.array(f)
            self.volT = interp1d(self.T_vib, self.volT)(self.T)
            self.GibT = interp1d(self.T_vib, self.GibT)(self.T)
        self.Slat[np.isnan(self.Slat)] = 0


    def calc_TE_V_fitF(self):

        #self.energies_orig = copy.deepcopy(self.energies)
        self.energies = BMfitF(self.volumes, self.volumes, self.energies_orig, self.BMfunc)

        self.blat = []
        for i in range(len(self.T)):
            FF = self.Flat[:,i] + self.theall[0,i,:]
            #p2 = np.poly1d(np.polyfit(self.volumes, FF, 2))
            #FF = p2(self.volumes)
            p1 = np.poly1d(np.polyfit(self.volumes, FF, 1))
            FF = p1(self.volumes)
            blat, self.volT[i], newF = CenDifB(self.volumes, self.energies+FF, N=7,kind='cubic')
            self.blat.append(blat)
            if newF!=None: self.GibT[i] = newF
            if self.volT[i] < 0: break

        nT = len(self.blat)
        self.beta = np.zeros((nT), dtype=float)
        if self.blat[nT-1] < 0: nT -= 1
        if nT < 2:
            self.blat[0] = -1
            return

        f2 = splrep(self.T[0:nT], self.volT[0:nT])
        self.beta[0:nT] = splev(self.T[0:nT], f2, der=1)/self.volT[0:nT]
        if self.T[0] == 0: self.beta[0] = 0


    def calc_TE_V_general(self,i,kind='cubic'):
        #if self.elmode>=1:
        if self.smooth:
            E0, Flat, Fel, Slat, Sel = BMsmooth(self.volumes, self.energies, self.Flat[:,i],
                self.theall[0,i,:], self.Slat[:,i], self.theall[1,i,:], self.BMfunc, self.elmode)
        else:
            E0, Flat, Fel, Slat, Sel = self.energies, self.Flat[:,i], \
                self.theall[0,i,:], self.Slat[:,i], self.theall[1,i,:]

        #try:
        if True:
            FF = Flat+Fel

            if self.eqmode==4 or self.eqmode==5:
                #print ("iiii", self.T[i], self.volT[i], E0+Flat+Fel)
                blat, self.volT[i], self.GibT[i], P = BMDifB(self.volumes, E0+FF, self.BMfunc, N=7, _T=self.T[i])
                if blat < 0: return -1.0, 0.0
                if self.T[i]!=0.0: beta = (BMfitP(self.volT[i], self.volumes, E0+FF+self.T[i]*(Slat+Sel),
                    self.BMfunc) + P)/self.T[i]/blat
                else: beta = 0.0
            else:
                blat, self.volT[i], newF = CenDifB(self.volumes, E0+FF, N=7,kind=kind)
                if blat < 0: return -1.0, 0.0
                if newF!=None: self.GibT[i] = newF
                if self.T[i]!=0.0:
                    try:
                        beta = CenDif(self.volT[i], self.volumes, Slat+Sel, N=7,kind=kind)/blat
                    except:
                        return -1.0, 0.0
                else: beta = 0.0
            return blat, beta
        #except:
        #return -1.0, 0.0


    def dot(self,i,b):
        dx0 = self.T[i] - self.T[i-1]
        dx1 = self.T[i+1] - self.T[i]
        dy0 = b[i] - b[i-1]
        dy1 = b[i+1] - b[i]
        s0 = math.sqrt(dx0*dx0+dy0*dy0)
        s1 = math.sqrt(dx1*dx1+dy1*dy1)
        return (dx0*dx1+dy0*dy1)/s0/s1

    def calc_thermodynamics(self):
        Faraday_constant = physical_constants["Faraday constant"][0]
        electron_volt = physical_constants["electron volt"][0]
        angstrom = 1e-30
        toJmol = Faraday_constant/self.natoms
        toGPa = electron_volt/angstrom*1.e-9

        thermofile = self.phasename+'/'+self.outf
        with open(thermofile, 'w') as fvib:
            fvib.write('#Found quasiharmonic mode : {}\n'.format(self.qhamode))
            if self.hasSCF:
                fvib.write('#T(K), volume, F(eV), S(J/K), H(J/K), a(-6/K), Cp(J/mol), Cv, Cpion, Bt(GPa), T_ph-D(K), T-D(K), F_el_atom, S_el_atom, C_el_atom, M_el, Seebeck_coefficients(10**-6 V/K), Lorenz_number(WOK^{-2}), Q_el, Q_p, Q_e, C_mu, W_p, W_e, Y_p, Y_e\n')
            else:
                fvib.write('#T, F_el_atom, S_el_atom, C_el_atom, M_el, seebeck_coefficients, Lorenz_number[WΩK−2], Q_el, Q_p, Q_e, C_mu, W_p, W_e, Y_p, Y_e Vol Gibbs_energy\n')

            #round one
            if self.hasSCF:
                if self.fitF:
                    self.calc_TE_V_fitF()
                    nT = len(self.beta)
                    #_beta = copy.deepcopy(self.beta)
                else:
                    self.blat = np.zeros((len(self.T)), dtype=float)
                    self.beta = np.zeros((len(self.T)), dtype=float)
                    nT = len(self.T)
                    for i in range(len(self.T)):
                        if self.elmode>=1:
                            self.blat[i], self.beta[i] = self.calc_TE_V_general(i, kind='UnivariateSpline')
                        else:
                            self.blat[i], self.beta[i] = self.calc_TE_V_general(i, kind='cubic')
                        if self.blat[i] < 0:
                            nT = i
                            print ("\nblat<0! Perhaps it has reached the upvolume limit at T =", self.T[i], "\n")
                            break
                """
                    _beta = copy.deepcopy(self.beta)

                    if self.eqmode==4 or self.eqmode==5:
                        #f2 = splrep(self.T[0:nT], self.volT[0:nT], s=3, k=5)
                        f2 = splrep(self.T[0:nT], self.volT[0:nT])
                        if 1==0:
                            self.beta[0:nT] = BSpline(self.T[0:nT], f2, der=1)/self.volT[0:nT]
                        else:
                            self.beta[0:nT] = splev(self.T[0:nT], f2, der=1)/self.volT[0:nT]
                        if self.T[0] == 0: self.beta[0] = 0

                #optimized LTC
                _beused = copy.deepcopy(_beta)
                _b2 = copy.deepcopy(self.beta)
                _bsplev = copy.deepcopy(self.beta)
                for ii in range (32):
                    for i in range(1, len(self.T)-1):
                        if self.blat[i] < 0 or self.blat[i+1] < 0: break
                        t1 = self.dot(i, _beused)
                        t2 = self.dot(i, _b2)
                        #t1 = (_beused[i]-_beused[i-1])*(_beused[i+1]-_beused[i])
                        #t2 = (_b2[i]-_b2[i-1])*(_b2[i+1]-_b2[i])
                        if abs(t1) >=abs(t2):
                            tused = _beused[i]
                        else:
                            tused = _b2[i]
                        _beused[i]=tused
                        _b2[i]=tused
                self.beta = _beused

                LTCzigzag = 0
                for i in range(1, len(self.T)-1):
                    if self.blat[i] < 0 or self.blat[i+1] < 0: break
                    if self.beta[i] < 1.e-6: continue #check irregularity of along LTC curve
                    if (self.beta[i]-self.beta[i-1])*(self.beta[i+1]-self.beta[i]) < 0.0: LTCzigzag += 1
                self.key_comments['LTC quality'] = LTCzigzag
                """


            if nT <3: return np.array(self.volumes)/self.natoms, np.array(self.energies_orig)/self.natoms, thermofile
            self.TupLimit = self.T[-1]
            self.Cp = []
            self.Cv = []
            k_ph_fac = 0
            theta_D = None
            for i in range(len(self.T)):
                if self.hasSCF:
                    blat, beta = self.blat[i], self.beta[i]
                    if blat < 0:
                        self.TupLimit = self.T[i-1]
                        print ("\nblat<0! Perhaps it has reached the upvolume limit at T =", self.T[i], "\n")
                        break
                    try:
                        slat = interp1d(self.volumes, self.Slat[:,i])(self.volT[i])
                        clat = interp1d(self.volumes, self.Clat[:,i])(self.volT[i])
                        flat = interp1d(self.volumes, self.Flat[:,i])(self.volT[i])
                        dlat = interp1d(self.volumes, self.Dlat)(self.volT[i])
                        cplat = clat+beta*beta*blat*self.volT[i]*self.T[i]
                    except:
                        self.TupLimit = self.T[i-1]
                        print ("\nPerhaps it has reached the upvolume limit at T =", self.T[i], "\n")
                        break
                #print("xxxxxxxxxxxxxxx",self.T[i], blat, self.volT[i])
                prp_T = np.zeros((self.theall.shape[0]))
                for j in range(len(prp_T)):
                    prp_T[j] = interp1d(self.volumes, self.theall[j,i,:], kind='cubic')(self.volT[i])

                # 0 K limit for the Lorenz number
                L = 2.443004551768e-08 #(1.380649e-23/1.60217662e-19*3.14159265359)**2/3
                if prp_T[5] > 1.e-16: L = prp_T[2]/prp_T[5]*k_B #avoiding divided by zero

                if self.hasSCF:
                    if self.T[i]==0: gamma=0
                    else: gamma = beta*blat*self.volT[i]/clat
                    debyeT = get_debye_T_from_phonon_Cv(self.T[i], clat, dlat, self.natoms)
                    try:
                        if self.T[i] == 0:
                            k_ph_fac = 0
                            T_div = 1
                        elif self.k_ph_mode==1:
                        #elif True:
                            k_ph_fac, gamma = self.get_thermo_lattice_conductivity(beta, blat, self.volT[i], clat, debyeT)
                            T_div = self.T[i]
                        elif self.k_ph_mode==2:
                            k_ph_fac, gamma, theta_D = self.get_lattice_conductivity(beta, blat, self.volT[i], clat, theta_D=theta_D)
                            T_div = self.T[i]
                        elif k_ph_fac==0:
                            if self.k_ph_mode==0:
                                k_ph_fac, gamma = self.get_lattice_conductivity_1(self.volT[0])
                            else:
                                k_ph_fac, gamma = self.get_lattice_conductivity_0(i, self.T[i], self.volT[i])
                            T_div = self.T[i]
                        else:
                            T_div = self.T[i]
                    except:
                        k_ph_fac = 0
                        T_div = 1
                        pass
                    k_ph = k_ph_fac/T_div*self.volT[i]**(1./3)

                    self.Cp.append(cplat+prp_T[2])
                    self.Cv.append(clat+prp_T[2])
                    fvib.write('{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.
                    format(self.T[i], self.volT[i]/self.natoms, self.GibT[i]/self.natoms, (slat+prp_T[1])*toJmol,
                    (self.GibT[i]+self.T[i]*(slat+prp_T[1]))*toJmol,
                    beta/3., (cplat+prp_T[2])*toJmol, (clat+prp_T[2])*toJmol,
                    cplat*toJmol, blat*toGPa, debyeT, dlat,
                    prp_T[0]/self.natoms, prp_T[1]*toJmol, prp_T[2]*toJmol,
                    prp_T[3], prp_T[4], L, prp_T[5]/self.natoms, prp_T[6]/self.natoms,
                    prp_T[7]/self.natoms, prp_T[8]*toJmol, prp_T[10]/self.natoms,
                    prp_T[11]/self.natoms, prp_T[12]/self.natoms, prp_T[13]/self.natoms, k_ph, gamma))
                    #prp_T[7]/self.natoms, prp_T[8]*toJmol, _bsplev[i]/3.0, prp_T[10]/self.natoms,
                else:
                    #(T[i], F_el_atom[i], S_el_atom[i], C_el_atom[i], M_el[i], seebeck_coefficients[i], L,
                    #Q_el[i], Q_p[i], Q_e[i], C_mu[i], W_p[i], W_e[i], Y_p[i], Y_e[i])
                    fvib.write('{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.
                    format(self.T[i], prp_T[0], prp_T[1], prp_T[2], prp_T[3], prp_T[4], L,
                    prp_T[5], prp_T[6], prp_T[7], prp_T[8], prp_T[10], prp_T[11], prp_T[12], prp_T[13],
                    self.volT[i], self.GibT[i]))
        self.datasm(thermofile)
        return np.array(self.volumes)/self.natoms, np.array(self.energies_orig)/self.natoms, thermofile

    def datasm(self, fname):
        data = np.loadtxt(fname, comments="#", dtype=float)
        nT = data.shape[0]
        nF = data.shape[1]
        nSmooth = 11
        box = np.ones(nSmooth)/nSmooth
        if nT <nSmooth : return
        from scipy.signal import savgol_filter
        for i in range(1,nF):
            #data[:,i]=np.convolve(data[:,i], box, mode='same')
            data[:,i]=savgol_filter(data[:,i], nSmooth, 3)
        with open(fname+'_sm', 'w') as fout:
            with open(fname, 'r') as fin:
                lines = fin.readlines()
                for line in lines:
                    if line.startswith('#'): print(line.strip(),file=fout)

            for i in range(0,nT):
                for j in range(0,nF):
                    if j==nF-1: fout.write('{}\n'.format(data[i,j]))
                    else: fout.write('{} '.format(data[i,j])) 
    


    def add_comput_inf(self):
        if self.vasp_db!=None:
            self.key_comments['METADATA'] = {'tag':self.tag}
            self.key_comments['E-V'], self.key_comments['POSCAR'], self.key_comments['INCAR'] = \
                get_rec_from_metatag(self.vasp_db, self.tag)
            with open (self.phasename+'/POSCAR', 'w') as fp:
                fp.write(self.key_comments['POSCAR'])
        nT = min(len(self.volT),len(self.blat))
        for i in range(nT):
            if self.blat[i] < 0:
                nT = i
                break
        if nT < 3:
            self.key_comments['ERROR'] = "Fatal ERROR! Calculation corrupted due to certain reason! Perhaps very bad E-V curve!"
            return False

        vn = min(self.volT[0:nT])
        vx = max(self.volT[0:nT])
        for ix,vol in enumerate(self.volumes):
           if vol>vn: break
        try:
            n = 0
            q = 0.0
            for i in range(ix-1, len(self.volumes)):
                n = n+1
                q += self.quality[i]
                if self.volumes[i] > vx: break
            q /= n
            self.key_comments['phonon quality'] = '{:8.6}'.format(q)
        except:
            self.key_comments['phonon quality'] = '{:8.6}'.format(-1.0)
        return True


    def find_fitting_quality(self, xx, orig_points, fitted, readme):
        nT = len(self.volT)
        for i in range(nT):
            if self.blat[i] < 0:
                nT = i
                break
        vn = min(self.volT[0:nT])
        vx = max(self.volT[0:nT])
        for ix,vol in enumerate(self.volumes):
           if vol>vn: break
        #try:
        if True:
            n = 0
            q = 0.0
            qmax = 0.0
            qT = 0.0
            qmaxT = 0.0
            for i in range(ix-1, len(self.volumes)):
                n = n+1
                for k,t in enumerate(self.T):
                    for j in range(len(xx)-1):
                        if (xx[j]-self.volumes[i])*(xx[j+1]-self.volumes[i])<0:
                            fn0 = fitted[0][j]+(self.volumes[i]-xx[j])/(xx[j+1]-xx[j])*(fitted[0][j+1]-fitted[0][j])
                            fn = fitted[k][j]+(self.volumes[i]-xx[j])/(xx[j+1]-xx[j])*(fitted[k][j+1]-fitted[k][j])
                            q += abs(fn-orig_points[k][i])
                            qmax = max(qmax,abs(fn-orig_points[k][i]))
                            if k!=0:
                                #tmp = abs((fn-orig_points[k][i]-fn0+orig_points[0][i])/(fn-fn0))
                                tmp = abs((fn-orig_points[k][i]-fn0+orig_points[0][i])/(t))
                                qT += tmp
                                qmaxT = max(qmaxT,tmp)
                            break
                if self.volumes[i] > vx: break
            q /= n*len(self.T)*self.natoms
            qmax /= self.natoms
            qT /= n*(len(self.T)-1)*self.natoms
            qmaxT /= self.natoms
            readme['Helmholtz energy quality'] = '+-{:.1e} eV'.format(q)
            readme['Helmholtz energy max error'] = '{:.1e} eV'.format(qmax)
            readme['Entropy quality'] = '{:.1e} J/K'.format(qT*96484)
            readme['Entropy max error'] = '{:.1e} J/K'.format(qmaxT*96484)
        #except:
        #    readme['Helmholtz energy quality'] = '+-{} eV'.format(9999)


    def calc_eij(self):
        R = []
        #print (self.key_comments)
        lat = self.key_comments['E-V']['lattices']
        R_volumes = []
        for vol in self.volumes:
            if vol in self.key_comments['E-V']['volumes']:
                idx = self.key_comments['E-V']['volumes'].index(vol)
                R_volumes.append(vol)
                R.append(lat[idx])
        T = self.T[self.T <=self.TupLimit]
        nT = len(T)
        R = np.array(R)
        R_T = np.zeros((nT, 3, 3), dtype=float)
        for i in range(3):
            for j in range(3):
                f2 = splrep(R_volumes, R[:,i,j])
                R_T[:,i,j] = splev(self.volT[0:nT], f2)
        R_T_dT = np.zeros((nT, 3, 3), dtype=float)
        for i in range(3):
            for j in range(3):
                f2 = splrep(self.T[0:nT], R_T[:,i,j])
                R_T_dT[:,i,j] = splev(self.T[0:nT], f2, der=1)
        inv_R_T = np.linalg.inv(R_T)
        #eij_T = np.matmul(inv_R_T,R_T_dT)
        eij_T = np.matmul(R_T_dT, inv_R_T)

        with open (self.phasename+'/fvib_eij', 'w') as fp:
            fp.write('#T volume alpha alpha_xx alpha_yy alpha_zz alpha_yz alpha_zy alpha_xz alpha_zx alpha_xy alpha_yx\n')
            for i,m in enumerate(eij_T):
                fp.write('{} {:10.6f} {:12.4e} {:12.4e} {:12.4e} {:12.4e} {:12.4e} {:12.4e} '\
                    '{:12.4e} {:12.4e} {:12.4e} {:12.4e}\n'.format(T[i], self.volT[i]/self.natoms, self.beta[i]/3., \
                    m[0,0], m[1,1], m[2,2], m[1,2], m[2,1], m[0,2], m[2,0], m[0,1], m[1,0]))
        self.eij_T = eij_T


    def Cij_to_Moduli(self, C):
        A = (C[0,0] + C[1,1] + C[2,2])/3.
        B = (C[0,1] + C[0,2] + C[1,2])/3.
        C = (C[3,3] + C[4,4] + C[5,5])/3.
        Bv = (A + 2.*B)/3.
        Gv = (A - B + 3.*C)/5.
        Ev = 9.*Bv*Gv/(Gv+3.*Bv)
        Poisson_ration = 0.5*Ev/Gv - 1.0
        return Ev, Gv, Bv, Poisson_ration


    def calc_uniform(self, volumes, uniform, outf):
        if len(volumes) == 0: return
        T = uniform[0][:,1]
        nV = len(uniform)
        nT, nF = uniform[0].shape
        with open (outf, 'w') as fp:
          headerfmt = "#{:>11s} {:>9s}" + " ".join(18 * ["{:>25s}"])
          header = [
                    "mu-Ef[eV]", "T[K]", "N[e/uc]", "DOS(ef)[1/(Ha*uc)]", "S[V/K]",
                    "sigma/tau0[1/(ohm*m*s)]", "RH[m**3/C]", "kappae/tau0[W/(m*K*s)]",
                    "cv[J/(mole-atom*K)]", "chi[m**3/mol]",
                    "cv_x[J/(mole-atom*K)]", "S_x(V/K)", "N_x(e/cm^3)", "L(W*ohm/K**2)",
                    "L0_h", "L0_e", "M_h", "M_e", "N_h", "N_e"
                ]
          if outf.endswith("uniform_lambda"):
              header[5] = "sigma/lambda0[1/(ohm*m**2)]"
              header[7] = "kappae/lambda0[W/(m**2*K)]"
          elif outf.endswith("custom_tau"):
              header[5] = "sigma[1/(ohm*m)]"
              header[7] = "kappae[W/(m*K)]"
          print(headerfmt.format(*header).strip(), file=fp)
          rowfmt = "{:>14.8g} {:>9g}" + " ".join(18 * ["{:>25g}"])
          for i in range(nT):
            if T[i] <min(self.T): continue
            if T[i] >max(self.T): continue
            f2 = splrep(self.T, self.volT)
            vol = splev(T[i], f2)
            uniform_V = []
            for ii,u in enumerate(volumes):
                uniform_V.append(uniform[ii][i,:])
            uniform_V = np.array(uniform_V, dtype=float)
            values = []
            for j in range(nF):
              f2 = splrep(volumes, uniform_V[:,j])
              val = splev(vol, f2)
              values.append(val)
            print(rowfmt.format(*values), file=fp)
            #fp.write ("\n")


    def calc_Cij(self):
        if len(self.VCij) == 0: return
        T = self.T[self.T <=self.TupLimit]
        nT = len(T)
        if min(self.volT[0:nT]) < min(self.VCij) or max(self.volT[0:nT]) > max(self.VCij): return
        self.Cij_T = np.zeros((nT, 6, 6), dtype=float)
        electron_volt = physical_constants["electron volt"][0]
        angstrom = 1e-30
        toGPa = electron_volt/angstrom*1.e-9
        for i in range(6):
            for j in range(6):
                f2 = splrep(self.VCij, self.Cij[:,i,j])
                self.Cij_T[:,i,j] = splev(self.volT[0:nT], f2)
        """
        if True:
                    print ("db_file",self.VCij)
                    print (self.volT[0])
                    print (self.Cij_T[0,:,:])
                    sys.exit()
        """

        with open (self.phasename+'/fvib_Cij_T', 'w') as fp:
            ngroup = self.space_group_number
            self.Young_Modulus_Cij = []
            self.Shear_Modulus_Cij = []
            self.Bulk_Modulus_Cij = []
            self.Poisson_Ratio_Cij = []
            for i,m in enumerate(self.Cij_T):
                E,G,B,Poisson_Ratio = self.Cij_to_Moduli(m)
                correction_factor = self.blat[i]*toGPa/B
                #correction_factor = 1.0
                for j in range(3):
                    for k in range(3):
                        self.Cij_T[i,j,k] *=correction_factor
                E,G,B,Poisson_Ratio = self.Cij_to_Moduli(self.Cij_T[i,:,:])
                self.Young_Modulus_Cij.append(E)
                self.Shear_Modulus_Cij.append(G)
                #self.Bulk_Modulus_Cij.append(B)
                self.Bulk_Modulus_Cij.append(correction_factor)
                self.Poisson_Ratio_Cij.append(Poisson_Ratio)

            if ngroup>=1 and ngroup<=2: #for Triclinic system
                fp.write('# T(K) V(Ang^3) B(GPa) C11 C12 C13 C14 C15 C16 C22 C23 C24 C25 C26 C33 C34'\
                    ' C35 C36 C44 C45 C46 C55 C56 C66 E(GPa) G(GPa) B_correction_factor\n')
                fm = "{} {:10.6f}"+" {:10.4f}"*(21+4)
                for i,m in enumerate(self.Cij_T):
                    print(fm.format\
                    (T[i], self.volT[i]/self.natoms, self.blat[i]*toGPa, \
                    m[0,0], m[0,1], m[0,2], m[0,3], m[0,4], m[0,5], m[1,1], m[1,2], m[1,3], m[1,4], m[1,5], m[2,2], m[2,3], m[2,4], m[2,5], m[3,3], m[3,4], m[3,5], m[4,4], m[4,5], m[5,5], \
                    self.Young_Modulus_Cij[i], self.Shear_Modulus_Cij[i], self.Bulk_Modulus_Cij[i]),file=fp)
            elif ngroup>=3 and ngroup<=15: # for Monoclinic system
                fp.write('# T(K) V(Ang^3) B(GPa) C11 C12 C13 C16 C22 C23 C26 C33 C36 C44 C45 C55 C66'\
                    ' E(GPa) G(GPa) B_correction_factor\n')
                fm = "{} {:10.6f}"+" {:10.4f}"*(13+4)
                for i,m in enumerate(self.Cij_T):
                    print(fm.format\
                    (T[i], self.volT[i]/self.natoms, self.blat[i]*toGPa, \
                    m[0,0], m[0,1], m[0,2], m[0,5], m[1,1], m[1,2], m[1,5], m[2,2], m[2,5], m[3,3], m[3,4], m[4,4], m[5,5], \
                    self.Young_Modulus_Cij[i], self.Shear_Modulus_Cij[i], self.Bulk_Modulus_Cij[i]),file=fp)
            elif ngroup>=16 and ngroup<=74: # for Orthorhombic system
                fp.write('# T(K) V(Ang^3) B(GPa) C11 C12 C13 C22 C23 C33 C44 C55 C66 E(GPa) G(GPa) B_correction_factor\n')
                fm = "{} {:10.6f}"+" {:10.4f}"*(9+4)
                for i,m in enumerate(self.Cij_T):
                    print(fm.format\
                    (T[i], self.volT[i]/self.natoms, self.blat[i]*toGPa, \
                    m[0,0], m[0,1], m[0,2], m[1,1], m[1,2], m[2,2], m[3,3], m[4,4], m[5,5],\
                    self.Young_Modulus_Cij[i], self.Shear_Modulus_Cij[i], self.Bulk_Modulus_Cij[i]),file=fp)
            elif ngroup>=75 and ngroup<=142: # for Tetragonal system
                fp.write('# T(K) V(Ang^3) B(GPa) C11 C12 C13 C33 C44 C66 E(GPa) G(GPa) B_correction_factor\n')
                fm = "{} {:10.6f}"+" {:10.4f}"*(6+4)
                for i,m in enumerate(self.Cij_T):
                    print(fm.format\
                    (T[i], self.volT[i]/self.natoms, self.blat[i]*toGPa, \
                    m[0,0], m[0,1], m[0,2], m[2,2], m[3,3], m[5,5],\
                    self.Young_Modulus_Cij[i], self.Shear_Modulus_Cij[i], self.Bulk_Modulus_Cij[i]),file=fp)
            elif ngroup>=143 and ngroup<=194: # for Trigonal system
                fp.write('# T(K) V(Ang^3) B(GPa) C11 C12 C13 C33 C44 E(GPa) G(GPa) B_correction_factor\n')
                fm = "{} {:10.6f}"+" {:10.4f}"*(5+4)
                for i,m in enumerate(self.Cij_T):
                    print(fm.format\
                    (T[i], self.volT[i]/self.natoms, self.blat[i]*toGPa, \
                    m[0,0], m[0,1], m[0,2], m[2,2], m[3,3],\
                    self.Young_Modulus_Cij[i], self.Shear_Modulus_Cij[i], self.Bulk_Modulus_Cij[i]),file=fp)
            #elif ngroup>=168 and ngroup<=194: # for Hexagonal system
            elif ngroup>=195 and ngroup<=230: # for cubic system
                fp.write('# T(K) V(Ang^3) B(GPa) C11 C12 C44 E(GPa) G(GPa) B_correction_factor\n')
                fm = "{} {:10.6f}"+" {:10.4f}"*(3+4)
                for i,m in enumerate(self.Cij_T):
                    print(fm.format\
                    (T[i], self.volT[i]/self.natoms, self.blat[i]*toGPa, \
                    m[0,0], m[0,1], m[3,3], \
                    self.Young_Modulus_Cij[i], self.Shear_Modulus_Cij[i], self.Bulk_Modulus_Cij[i]),file=fp)

        with open (self.phasename+'/fvib_Mii_T', 'w') as fp:
            fp.write('# T(K) V(Ang^3) B(GPa) M1(GPa) M2(GPa) M3(GPa) M4(GPa) M5(GPa) M6(GPa)'
                ' E(GPa) G(GPa)\n')
            fm = "{} {:10.6f}"+" {:10.4f}"*(6+3)
            for i,m in enumerate(self.Cij_T):
                M = (m+m.T)/2.0
                w, v = np.linalg.eig(M)
                #w = np.real(w)
                w = sorted(w, reverse=True)
                print(fm.format\
                (T[i], self.volT[i]/self.natoms, self.blat[i]*toGPa, \
                w[0], w[1], w[2], w[3], w[4], w[5], \
                self.Young_Modulus_Cij[i], self.Shear_Modulus_Cij[i]),file=fp)


    def calc_Cij_S(self):
        if len(self.VCij) == 0: return
        T = self.T[self.T <=self.TupLimit]
        nT = len(T)
        #if min(self.volT) < min(self.VCij) or max(self.volT) > max(self.VCij): return
        if min(self.volT[0:nT]) < min(self.VCij) or max(self.volT[0:nT]) > max(self.VCij): return
        self.Cij_S = np.zeros((nT, 6, 6), dtype=float)
        electron_volt = physical_constants["electron volt"][0]
        angstrom = 1e-30
        toGPa = electron_volt/angstrom*1.e-9

        with open (self.phasename+'/fvib_Cij_S', 'w') as fp:
            ngroup = self.space_group_number
            self.Young_Modulus_Cij_S = []
            self.Shear_Modulus_Cij_S = []
            self.Bulk_Modulus_Cij_S = []
            self.Poisson_Ratio_Cij_S = []
            ev = np.zeros((6), dtype=float)
            for i,m in enumerate(self.Cij_T):
                eij = self.eij_T[i]
                ev[0] = eij[0,0]
                ev[1] = eij[1,1]
                ev[2] = eij[2,2]
                ev[3] = eij[1,2] + eij[1,2]
                ev[4] = eij[0,2] + eij[2,0]
                ev[5] = eij[0,1] + eij[1,0]
                ec = np.matmul(ev,m)/toGPa
                """
                if T[i] == 0.0:
                    print (m)
                    print (ec)
                    sys.exit()
                """
                for j in range(6):
                    for k in range(6):
                        if self.Cv[i] > 1.e-8:
                            self.Cij_S[i, j, k] = self.Cij_T[i, j, k] + T[i]*self.volT[i]/self.Cv[i]*ec[j]*ec[k]*toGPa
                        else:
                            self.Cij_S[i, j, k] = self.Cij_T[i, j, k]

                E,G,B,Poisson_Ratio = self.Cij_to_Moduli(self.Cij_S[i, :, :])
                self.Young_Modulus_Cij_S.append(E)
                self.Shear_Modulus_Cij_S.append(G)
                self.Bulk_Modulus_Cij_S.append(B)
                self.Poisson_Ratio_Cij_S.append(Poisson_Ratio)

            if ngroup>=1 and ngroup<=2: #for Triclinic system
                fp.write('# T(K) V(Ang^3) B(GPa) C11 C12 C13 C14 C15 C16 C22 C23 C24 C25 C26 C33 C34'\
                    ' C35 C36 C44 C45 C46 C55 C56 C66 E(GPa) G(GPa) B_correction_factor\n')
                fm = "{} {:10.6f}"+" {:10.4f}"*(21+4)
                for i,m in enumerate(self.Cij_S):
                    if self.Cv[i] > 1.e-8:
                        Bs = self.blat[i]*self.Cp[i]/self.Cv[i]
                    else:
                        Bs = self.blat[i]
                    print(fm.format\
                    (T[i], self.volT[i]/self.natoms, Bs*toGPa, \
                    m[0,0], m[0,1], m[0,2], m[0,3], m[0,4], m[0,5], m[1,1], m[1,2], m[1,3], m[1,4], m[1,5], m[2,2], m[2,3], m[2,4], m[2,5], m[3,3], m[3,4], m[3,5], m[4,4], m[4,5], m[5,5], \
                    self.Young_Modulus_Cij_S[i], self.Shear_Modulus_Cij_S[i], self.Bulk_Modulus_Cij[i]),file=fp)
            elif ngroup>=3 and ngroup<=15: # for Monoclinic system
                fp.write('# T(K) V(Ang^3) B(GPa) C11 C12 C13 C16 C22 C23 C26 C33 C36 C44 C45 C55 C66'\
                    ' E(GPa) G(GPa) B_correction_factor\n')
                fm = "{} {:10.6f}"+" {:10.4f}"*(13+4)
                for i,m in enumerate(self.Cij_S):
                    if self.Cv[i] > 1.e-8:
                        Bs = self.blat[i]*self.Cp[i]/self.Cv[i]
                    else:
                        Bs = self.blat[i]
                    print(fm.format\
                    (T[i], self.volT[i]/self.natoms, Bs*toGPa, \
                    m[0,0], m[0,1], m[0,2], m[0,5], m[1,1], m[1,2], m[1,5], m[2,2], m[2,5], m[3,3], m[3,4], m[4,4], m[5,5], \
                    self.Young_Modulus_Cij_S[i], self.Shear_Modulus_Cij_S[i], self.Bulk_Modulus_Cij[i]),file=fp)
            elif ngroup>=16 and ngroup<=74: # for Orthorhombic system
                fp.write('# T(K) V(Ang^3) B(GPa) C11 C12 C13 C22 C23 C33 C44 C55 C66 E(GPa) G(GPa) B_correction_factor\n')
                fm = "{} {:10.6f}"+" {:10.4f}"*(9+4)
                for i,m in enumerate(self.Cij_S):
                    if self.Cv[i] > 1.e-8:
                        Bs = self.blat[i]*self.Cp[i]/self.Cv[i]
                    else:
                        Bs = self.blat[i]
                    print(fm.format\
                    (T[i], self.volT[i]/self.natoms, Bs*toGPa, \
                    m[0,0], m[0,1], m[0,2], m[1,1], m[1,2], m[2,2], m[3,3], m[4,4], m[5,5],\
                    self.Young_Modulus_Cij_S[i], self.Shear_Modulus_Cij_S[i], self.Bulk_Modulus_Cij[i]),file=fp)
            elif ngroup>=75 and ngroup<=142: # for Tetragonal system
                fp.write('# T(K) V(Ang^3) B(GPa) C11 C12 C13 C33 C44 C66 E(GPa) G(GPa) B_correction_factor\n')
                fm = "{} {:10.6f}"+" {:10.4f}"*(6+4)
                for i,m in enumerate(self.Cij_S):
                    if self.Cv[i] > 1.e-8:
                        Bs = self.blat[i]*self.Cp[i]/self.Cv[i]
                    else:
                        Bs = self.blat[i]
                    print(fm.format\
                    (T[i], self.volT[i]/self.natoms, Bs*toGPa, \
                    m[0,0], m[0,1], m[0,2], m[2,2], m[3,3], m[5,5],\
                    self.Young_Modulus_Cij_S[i], self.Shear_Modulus_Cij_S[i], self.Bulk_Modulus_Cij[i]),file=fp)
            elif ngroup>=143 and ngroup<=194: # for Trigonal system
                fp.write('# T(K) V(Ang^3) B(GPa) C11 C12 C13 C33 C44 E(GPa) G(GPa) B_correction_factor\n')
                fm = "{} {:10.6f}"+" {:10.4f}"*(5+4)
                for i,m in enumerate(self.Cij_S):
                    if self.Cv[i] > 1.e-8:
                        Bs = self.blat[i]*self.Cp[i]/self.Cv[i]
                    else:
                        Bs = self.blat[i]
                    print(fm.format\
                    (T[i], self.volT[i]/self.natoms, Bs*toGPa, \
                    m[0,0], m[0,1], m[0,2], m[2,2], m[3,3],\
                    self.Young_Modulus_Cij_S[i], self.Shear_Modulus_Cij_S[i], self.Bulk_Modulus_Cij[i]),file=fp)
            #elif ngroup>=168 and ngroup<=194: # for Hexagonal system
            elif ngroup>=195 and ngroup<=230: # for cubic system
                fp.write('# T(K) V(Ang^3) B(GPa) C11 C12 C44 E(GPa) G(GPa) B_correction_factor\n')
                fm = "{} {:10.6f}"+" {:10.4f}"*(3+4)
                for i,m in enumerate(self.Cij_S):
                    if self.Cv[i] > 1.e-8:
                        Bs = self.blat[i]*self.Cp[i]/self.Cv[i]
                    else:
                        Bs = self.blat[i]
                    print(fm.format\
                    (T[i], self.volT[i]/self.natoms, Bs*toGPa, \
                    m[0,0], m[0,1], m[3,3], \
                    self.Young_Modulus_Cij_S[i], self.Shear_Modulus_Cij_S[i], self.Bulk_Modulus_Cij[i]),file=fp)


    def run_console(self):
        finished_tags = finished_calc()
        if not self.renew:
            if self.tag is not None:
                if self.tag in finished_tags:
                    print ("\nWARNING: previous calculation existed. supply -renew in the options to recalculate.\n")
                    return None, None, None, None

        if self.local!="": self.find_static_calculations_local()
        elif self.vasp_db!=None: self.find_static_calculations()
        else: self.find_static_calculations_without_DB()
        """
        if not self.fitF:
            if os.path.exists(self.phasename+'/fitF'):
                print ("\nWARNING: phase skipped due to file 'fitF' seen in ", self.phasename, "\n")
                return None, None, None, None
        else:
            if not os.path.exists(self.phasename+'/fitF'):
                print ("\nWARNING: phase skipped due to file 'fitF' not seen in ", self.phasename, "\n")
                return None, None, None, None
        """

        if not self.renew:
            pdis298 = self.phasename+'/figures/vdis298.15.png'
            if os.path.exists(pdis298):
                print ("\nWARNING: previous calculation existed. supply -renew in the options to recalculate.\n")
                return None, None, None, None

        self.find_vibrational()

        T = np.array(self.T_vib)
        if self.pyphon:
            self.T = self.T_vib
        elif self.td < 0:
            self.T = T_remesh(min(self.T_vib), min(self.t1,max(self.T_vib)), self.td, _nT=self.nT)
            #print ("xxxxx 2", len(self.T))
        else:
            self.T = T[T<=self.t1]
            #print ("xxxxx 3", len(self.T))

        if not self.pyphon: self.get_qha()
        self.hasSCF = True
        if not self.check_vol(): return None, None, None, None
        if self.local!="":
            self.has_Cij = self.get_Cij(self.phasename, pinv=False)
            self.has_btp2 = self.get_btp2(self.local)
        else:
            self.has_Cij = self.get_Cij(os.path.join(self.phasename,'Yphon'), pinv=False)
            self.has_btp2 = self.get_btp2(self.phasename)
        self.energies_orig = copy.deepcopy(self.energies)

        if self.noel : self.theall = np.zeros([14, len(self.T), len(self.volumes)])
        else : self.calc_thermal_electron()
        a,b,c = self.calc_thermodynamics()
        if self.add_comput_inf():
            self.calc_eij()
            if self.has_Cij:
                self.calc_Cij()
                self.calc_Cij_S()
            if self.has_btp2:
                self.calc_uniform(self.volumes_uniform_tau, self.uniform_tau,
                    os.path.join(self.local,'interpolation.dope.trace.uniform_tau'))
                self.calc_uniform(self.volumes_uniform_lambda, self.uniform_lambda,
                    os.path.join(self.local,'interpolation.dope.trace.uniform_lambda'))
            return a,b,c,self.key_comments
        else: return a,b,c,self.key_comments


    def calc_TE_V_general_for_plot(self,i,xx, kind='cubic'):
        if self.smooth:
            E0, Flat, Fel, Slat, Sel = BMsmooth(self.volumes, self.energies, self.Flat[:,i],
                self.theall[0,i,:], self.Slat[:,i], self.theall[1,i,:], self.BMfunc, self.elmode)
        else:
            E0, Flat, Fel, Slat, Sel = self.energies, self.Flat[:,i], \
                self.theall[0,i,:], self.Slat[:,i], self.theall[1,i,:]

        if True:
            FF = E0+Flat+Fel

            if self.eqmode==4 or self.eqmode==5:
                f, pcov = alt_curve_fit(self.BMfunc, self.volumes, FF)
                return BMvol(xx,f)
            else:
                f2 = interp1d(self.volumes, FF, kind=kind)
                return f2(xx)


    def calc_free_energy_for_plot(self, readme):
        vn = min(self.volumes)
        vx = max(self.volumes)
        xx = np.linspace(vn,vx,1000)
        fitted = []
        orig_points = []
        for i in range(len(self.T)):
            if self.fitF:
                FF = self.Flat[:,i] + self.theall[0,i,:]
                p1 = np.poly1d(np.polyfit(self.volumes, FF, 1))
                FF = p1(self.volumes)
                f2 = interp1d(self.volumes, self.energies+FF, kind='cubic')
                yy = f2(xx)
            elif self.elmode>=1:
                yy = self.calc_TE_V_general_for_plot(i, xx, kind='UnivariateSpline')
            else:
                yy = self.calc_TE_V_general_for_plot(i, xx, kind='cubic')
            fitted.append(yy)
            orig_points.append((self.energies_orig+self.Flat[:,i] + self.theall[0,i,:]))

        self.find_fitting_quality(xx, orig_points, fitted, readme)
        return np.array(self.volumes)/self.natoms, np.array(xx)/self.natoms, \
            np.array(orig_points)/self.natoms, np.array(fitted)/self.natoms


    def get_free_energy_for_plot(self, readme):
        if not self.hasSCF or not self.pyphon: return None
        self.T = np.linspace(0, self.TupLimit, 17, endpoint=True)
        self.toYphon(_T=self.T,for_plot=True)
        self.check_vol()
        self.energies_orig = copy.deepcopy(self.energies)

        if self.noel : self.theall = np.zeros([14, len(self.T), len(self.volumes)])
        else : self.calc_thermal_electron()
        v,x,o,f = self.calc_free_energy_for_plot(readme)
        return v, x, self.T, o, f


    def get_formula(self):
        if self.local!="":
            return self.formula_pretty
        elif self.vasp_db==None:
            ss = [s for s in self.phasename.split('/') if s!=""]
            return ss[-1].split('_')[0]
        else:
            return self.formula_pretty


    def calc_single(self):
        Faraday_constant = physical_constants["Faraday constant"][0]
        electron_volt = physical_constants["electron volt"][0]
        angstrom = 1e-30
        volume = 0
        if self.poscar!=None:
            structure = Structure.from_file(self.poscar)
            volume = structure.volume
            self.natoms = len(structure.sites)
            self.formula_pretty = structure.composition.reduced_formula
            try:
                formula2composition(formula_pretty)
            except:
                self.formula_pretty = reduced_formula(structure.composition.alphabetical_formula)
            sa = SpacegroupAnalyzer(structure)
            self.phase = sa.get_space_group_symbol().replace('/','.')+'_'+str(sa.get_space_group_number())
            if self.phasename==None:
                self.phasename = self.formula_pretty+'_'+self.phase
        else:
            volume = 0
            self.natoms = 1

        if self.phasename is None:
            self.phasename = 'unknown'

        self.theall = np.zeros([14, len(self.T)])
        if self.doscar!=None:
            with open(self.doscar, "r") as fp:
                prp_vol = runthelec(self.t0, self.t1, self.td, self.xdn, self.xup, self.dope, self.ndosmx,
                    self.gaussian, self.natfactor, dos=fp, _T=self.T, fout=sys.stdout, vol=volume)
                self.theall[:,:] = np.array( prp_vol ) # convert Tuple into array for the convenience of quasistatic interpolation

        if self.vdos is not None:
            with open(self.vdos, "r") as fp:
                Flat, U_ph, Slat, Clat, C_ph_n, Sound_ph, Sound_nn, N_ph, NN_ph, debyeT, quality, vdos_natoms \
                    = ywpyphon.vibrational_contributions(self.T, dos_input=fp, energyunit='eV')
        else:
            vdos_natoms = 1
            Flat = np.zeros((len(self.T)), type=float)
            Slat = np.zeros((len(self.T)), type=float)
            Clat = np.zeros((len(self.T)), type=float)

        try:
            self.key_comments
        except:
            self.key_comments = {}

        self.key_comments['phonon quality'] = '{:8.6}'.format(quality)
        print("\nnatoms in vdos.out=", vdos_natoms, "natoms in POSCAR=", self.natoms, "volume=", volume, "\n")

        toJmol = Faraday_constant/vdos_natoms
        elfac = self.natoms/vdos_natoms
        toJmolel = toJmol/elfac
        thermofile = self.phasename+'/'+self.outf
        if not os.path.exists(self.phasename):
            os.mkdir(self.phasename)

        with open(thermofile, 'w') as fvib:
            fvib.write('#T(K), volume, F(eV), S(J/K), H(J/K), a(-6/K), Cp(J/mol), Cv, Cpion, Bt(GPa), T_ph-D(K), T-D(K), F_el_atom, S_el_atom, C_el_atom, M_el, Seebeck_coefficients(μV/K), Lorenz_number(WΩK^{−2}), Q_el, Q_p, Q_e, C_mu, W_p, W_e, Y_p, Y_e\n')

            for i in range(len(self.T)):
                prp_T = self.theall[:,i]
                # 0 K limit for the Lorenz number
                L = 2.443004551768e-08 #(1.380649e-23/1.60217662e-19*3.14159265359)**2/3
                if prp_T[5] > 1.e-16: L = prp_T[2]/prp_T[5]*k_B #avoiding divided by zero
                f = Flat[i]+prp_T[0]/elfac
                s = Slat[i]+prp_T[1]/elfac
                c = Clat[i]+prp_T[2]/elfac
                debyeT = get_debye_T_from_phonon_Cv(self.T[i], Clat[i], 400., vdos_natoms)
                fvib.write('{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.
                format(self.T[i], volume/self.natoms, f/vdos_natoms, s*toJmol, (f+self.T[i]*s)*toJmol,
                0.0, c*toJmol, c*toJmol, Clat[i]*toJmol, 0, debyeT, 0,
                prp_T[0]/self.natoms, prp_T[1]*toJmolel, prp_T[2]*toJmolel,
                prp_T[3], prp_T[4], L, prp_T[5]/self.natoms, prp_T[6]/self.natoms,
                prp_T[7]/self.natoms, prp_T[8]*toJmolel, prp_T[10]/self.natoms,
                prp_T[11]/self.natoms, prp_T[12]/self.natoms, prp_T[13]/self.natoms))
        return thermofile


    def run_single(self):
        if self.td < 0:
            self.T = T_remesh(self.t0, self.t1, self.td, _nT=self.nT)
        else:
            nT = int((self.t1-self.t0)/self.td+1.5)
            self.T = np.linspace(self.t0, self.t1, nT, endpoint=True)
        thermofile = self.calc_single()
        return thermofile,self.key_comments,self.natoms


if __name__ == '__main__':
    # initialize temperatures
    t0 = 0  # low temperature
    t1 = 1300  # high temperature
    td = 10  #
    # both double precision
    xdn = -100  # what is this XXX
    xup = 100  # what is this XXX
    # back to reals
    ndosmx = 10001  # 10001 #aka n_dos
    dope = 0.0
    natom = 1
    gaussian = 1000. #densed mesh near Fermi energy
    outf = "fvib_ele"

    # handling the command line option
    # TODO: use proper argparse module for this
    count = 1
    while (count < len(sys.argv)):
      if (sys.argv[count] == "-T1"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        t1 = float(sys.argv[count])
      elif (sys.argv[count] == "-dT"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        td = float(sys.argv[count])
      elif (sys.argv[count] == "-dope"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        dope = float(sys.argv[count])
        if abs(dope)<5.e-9:
          ndosmx = 100001
          gaussian = 10000.
      elif (sys.argv[count] == "-gauss"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        gaussian = float(sys.argv[count])
      elif (sys.argv[count] == "-ne"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        ndosmx = int(sys.argv[count])
      elif (sys.argv[count] == "-natom"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        natom = int(sys.argv[count])
      elif (sys.argv[count] == "-outf"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        outf = str(sys.argv[count])
      count = count + 1

    with open(outf, 'w') as fvib:
        F_el_atom, S_el_atom, C_el_atom, M_el, seebeck_coefficients, Q_el, Q_p, Q_e, C_mu, T, W_p, W_e, Y_p, Y_e = runthelec(t0, t1, td, xdn, xup, dope, ndosmx, gaussian, natom)
        fvib.write('#T, F_el_atom, S_el_atom, C_el_atom, M_el, seebeck_coefficients, Lorenz_number[WΩK−2], Q_el, Q_p, Q_e, C_mu, W_p, W_e, Y_p, Y_e\n')
        for i in range(T.size):
            L = 2.443004551768e-08 #1.380649e-23/1.60217662e-19x3.14159265359xx2/3
            if Q_el[i] > 1.e-16: L = C_el_atom[i]/Q_el[i]*k_B
            fvib.write('{} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.format(T[i], F_el_atom[i], S_el_atom[i], C_el_atom[i], M_el[i], seebeck_coefficients[i], L, Q_el[i], Q_p[i], Q_e[i], C_mu[i], W_p[i], W_e[i], Y_p[i], Y_e[i]))

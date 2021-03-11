# -------- energy in eV, temperature in K
from __future__ import division
import sys
import math
import numpy as np
import scipy.constants
from scipy.constants import eV, physical_constants
from scipy.optimize import brentq
from scipy.integrate import cumtrapz, trapz, simps
from scipy.interpolate import interp1d
from dfttk.analysis.ywplot import myjsonout


kB_eV = physical_constants['Boltzmann constant in eV/K'][0]
kB = scipy.constants.k
h = scipy.constants.h
pi = scipy.constants.pi

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


# read phonon density of state from file f
def getdos(f): # Line 186
    """

    Parameters
    ----------
    f : file descripter for phonon dos from Yphon

    Returns
    -------
    freq : phonon frequency (array)
    pdos : phonon dos (array)
    quality : weigh of positive frequency
    """
    # read the file
    lines = f.readlines() # read in all lines then determine is it is WIEN2k DOS (in the unit eV) file or VASP DOS file
    freq = []
    pdos = []
    for line in lines:
        if line.startswith("#"): continue
        if line.strip()=="": continue
        data = [k for k in line.split(' ') if k != '']
        freq.append(float(data[0]))
        pdos.append(float(data[1]))

    freq = np.array(list(map(float,freq)))
    pdos = np.array(list(map(float,pdos)))
    quality = trapz(pdos,freq)
    pdos = pdos[freq>=0]
    freq = freq[freq>=0]
    NF = len(freq)
    Nlow = NF//100
    xfreq = freq[0:Nlow]
    yfreq = pdos[0:Nlow]
    cfreq = yfreq.sum()/(xfreq*xfreq).sum()
    A = freq[Nlow]
    B = A/Nlow
    a1 = B/10000
    d = 0.01*a1
    fnew = []
    pnew = []
    ff = 0
    i = 0
    while ff < A:
        fnew.append(ff)
        pnew.append(cfreq*ff*ff)
        i = i + 1
        ff += a1+i*d
    #print ("af=", A, B, d, a1, i)
    fnew.extend(freq[Nlow:])
    pnew.extend(pdos[Nlow:])
    good = trapz(pnew,fnew)
    natom = int(quality+0.5)//3
    quality = good/quality
    pnew = 3*natom/good*np.array(pnew)

    return np.array(list(map(float,fnew))), np.array(list(map(float,pnew))), quality, natom


def caclf(_freq, _pdos, T, dmu=0.0, energyunit='J'):
    """
    Calculate thermal free energy from phonon density of states (p DOS)

    Parameters

    _freq : phonon frequency
    -pdos : phobob DOS
    dmu : to be used external phonon chemical potential

    Returns

    f : vibrational free energy
    u : internal energy
    s : vibrational entropy
    cv : regular vibrational cv
    n : density of active phonons
    cv_n : vibrational cv of constant niumber of phonons
    sound_ph : phonon seebeck corfficient (freq/k)
    """

    hmu = h*_freq[np.where(_freq>0.0)]
    freq = _freq[np.where(_freq>0.0)]
    pdos = _pdos[np.where(_freq>0.0)]
    fn = pdos
    Nmode = trapz(fn, freq)
    fn = pdos*hmu*0.5
    u0 = trapz(fn, freq)
    #print ("u0=",u0, len(hmu))
    if T > 0: Beta = 1/(kB*T)
    else:
        if energyunit=='eV':
            #print ("pppppppppp", u0/eV)
            return u0/eV,u0/eV,0,0,0,0,0,0,0,0
        else:
            return u0,u0,0,0,0,0,0,0,0,0


    tc = Beta*(hmu-dmu)
    tc = tc[np.where(tc<500)]
    k1 = len(tc)
    tf = 1.0/(np.exp(tc)-1.0)
    hmu = hmu[0:k1]
    freq = freq[0:k1]
    pdos = pdos[0:k1]

    fn = pdos*tf*(hmu-dmu)
    tmp = trapz(fn, freq)
    active_freq = tmp/h
    lowT = active_freq/_freq.max() < 1.e-7
    if lowT:
        xfreq = freq[freq<1e-2*_freq.max()]
        yfreq = pdos[freq<1e-2*_freq.max()]
        cfreq = yfreq.sum()/(xfreq*xfreq).sum()
        #print ("af=", active_freq, active_freq/_freq.max(), cfreq, T)

    fn = pdos*tf
    nn = trapz(fn, freq)
    if lowT:
        nn = cfreq*(kB/h*T)**3*2.4041138064
    debye = nn/((kB/h*T)**3*2.4041138064)
    debye = (Nmode*3/debye)**(1/3)*h/kB
    #print ("debye=", debye, Nmode)

    fn = pdos*hmu*tf
    u_nn = trapz(fn, freq)
    u = u0+u_nn
    #print ("u=",u)

    #tf1 = tf + 1.e-60 # 1.e-60 is used to avoid log exception
    #fn = pdos*((1+tf)*np.log(1+tf)-tf1*np.log(tf1))
    fn = pdos*((1+tf)*np.log(1+tf)-tf*np.log(tf))
    s = trapz(fn, freq)*kB

    tf = tf*(1.0+tf)
    fn = pdos*tf
    n = trapz(fn, freq)
    if lowT:
        n = cfreq*(kB/h*T)**3*pi**2/3


    fn = pdos*tf*(hmu-dmu)
    tmp = trapz(fn, freq)
    u_n = tmp/n
    if lowT:
        u_n = kB*T*9*2.4041138064/pi**2
    sound_ph = u_n/h
    #print ("u_n=", u_n)

    fn = pdos*(hmu-dmu)**2*tf
    cv = trapz(fn, freq)/kB/T/T
    if lowT:
        cv = cfreq*kB*(kB/h*T)**3*4*pi**4/15
        s = cv/3.
        u = u0+cv*T/4.
        u_nn = cv*T/4.

    fn = pdos*(hmu-dmu-u_n)**2*tf
    cv_n = trapz(fn, freq)/kB/T/T
    if lowT:
        cv_n = cv - n*u_n*u_n/kB/T/T
        #print ("u_n=", cv, cv_n, u_n, cfreq, n*u_n*u_n/kB/T/T, T)


    #print ("cv_n=", cv_n, cv)

    #print (T, n, nn)
    if energyunit=='eV':
        return (u-T*s)/eV, u/eV, s/eV, cv/eV, cv_n/eV, sound_ph, u_nn/nn/h, n, nn, debye
    else:
        return u-T*s, u, s, cv, cv_n, sound_ph, u_nn/nn/h, n, nn, debye


def vibrational_contributions(T, dos_input=sys.stdin, _dmu=0.0, energyunit='J'):
    freq, pdos, quality, natom = getdos(dos_input)
    #print ("eeeeeeee", natom)
    nT = T.size
    F_ph = np.zeros(nT)
    U_ph = np.zeros(nT)
    S_ph = np.zeros(nT)
    C_ph_mu = np.zeros(nT) # phonon specific heat at constant mu
    NN_ph = np.zeros(nT) # total number of phonon
    N_ph = np.zeros(nT) # total number of thermal Carrier
    C_ph_n = np.zeros(nT) # phonon specific heat at constant N
    sound_ph = np.zeros(nT) # phonon seebeck corfficient (freq/k)
    sound_nn = np.zeros(nT) # averaged phonon frequency
    debyeT= np.zeros(nT) # averaged phonon frequency

    for i in range(nT):
        F_ph[i], U_ph[i], S_ph[i], C_ph_mu[i], C_ph_n[i], sound_ph[i], sound_nn[i], N_ph[i], NN_ph[i], debyeT[i] = caclf(freq, pdos, T[i], dmu=_dmu,energyunit=energyunit)
    #print ("eeeeee",C_ph_mu*96484)

    return F_ph, U_ph, S_ph, C_ph_mu, C_ph_n, sound_ph, sound_nn, N_ph, NN_ph, debyeT, quality, natom

if __name__ == '__main__':
    # initialize temperatures
    t0 = 0
    td = 10  #
    t1 = 1600  # high temperature
    #t0 = td  # low temperature
    natom = 1
    dmu = 0
    unit = 1

    # handling the command line option
    # TODO: use proper argparse module for this
    count = 1
    while (count < len(sys.argv)):
      if (sys.argv[count] == "-T0"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        t0 = float(sys.argv[count])
      elif (sys.argv[count] == "-T1"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        t1 = float(sys.argv[count])
      elif (sys.argv[count] == "-dT"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        td = float(sys.argv[count])
      elif (sys.argv[count] == "-dmu"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        dmu = float(sys.argv[count])
      elif (sys.argv[count] == "-natom"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        natom = int(sys.argv[count])
      elif (sys.argv[count] == "-moleatom"):
        unit = physical_constants['Avogadro constant'][0]
      count = count + 1

    unit = unit/natom
    # for all temperatures
    T = np.arange(t0,t1+td,td) # temperature
    F_ph, U_ph, S_ph, C_ph_mu, C_ph_n, Sound_ph, Sound_nn, N_ph, NN_ph, debyeT, quality, natom\
        = vibrational_contributions(T, dos_input=sys.stdin, _dmu=dmu, energyunit='J')

    sys.stderr.write ("\nThe phonon quality= {:08.6f}\n\n".format(quality))

    for i in range(T.size):
        tmp0 = 0.0
        tmp1 = 0.0
        if N_ph[i]!=0.:
            tmp0 = C_ph_mu[i]/N_ph[i]
            tmp1 = C_ph_n[i]/N_ph[i]

        sys.stdout.write('{:10.7g} {:10.7g} {:10.7g} {:10.7g} {:10.7g} {:10.7g} {:10.7g} \
        {:10.7g} {:10.7g} {:10.7g} {:10.7g} {:10.7g} {:10.7g}\n'.format(\
        T[i], F_ph[i]*unit, U_ph[i]*unit, S_ph[i]*unit, C_ph_mu[i]*unit, C_ph_n[i], \
        tmp0, tmp1, Sound_ph[i], Sound_nn[i], \
        N_ph[i], NN_ph[i], debyeT[i]))

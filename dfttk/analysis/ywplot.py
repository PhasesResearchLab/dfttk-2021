#!/usr/bin/python -x

import sys
from datetime import datetime

import os, fnmatch
import copy
import time
import datetime
import numpy as np
from scipy.optimize import linprog
from scipy.interpolate import interp1d
from scipy.interpolate import make_interp_spline, BSpline
from scipy import interpolate
from numpy.linalg import solve
from fractions import Fraction
from difflib import SequenceMatcher

#from PIL import Image
from scipy import misc
from scipy import ndimage as ndi
import math
import glob
from scipy.optimize import curve_fit
from scipy.constants import physical_constants
from scipy.optimize import brentq
from scipy.integrate import cumtrapz, trapz, simps
from pymatgen import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from dfttk.analysis.ywutils import get_expt, formula2composition, get_melting_temperature

import re
import json
import subprocess
from shutil import copyfile

from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw
from difflib import SequenceMatcher
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

#from elements import elements

MM_of_Elements = {'H': 1.00794, 'He': 4.002602, 'Li': 6.941, 'Be': 9.012182, 'B': 10.811, 'C': 12.0107, 'N': 14.0067,
              'O': 15.9994, 'F': 18.9984032, 'Ne': 20.1797, 'Na': 22.98976928, 'Mg': 24.305, 'Al': 26.9815386,
              'Si': 28.0855, 'P': 30.973762, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.0983, 'Ca': 40.078,
              'Sc': 44.955912, 'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938045,
              'Fe': 55.845, 'Co': 58.933195, 'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.409, 'Ga': 69.723, 'Ge': 72.64,
              'As': 74.9216, 'Se': 78.96, 'Br': 79.904, 'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.90585,
              'Zr': 91.224, 'Nb': 92.90638, 'Mo': 95.94, 'Tc': 98.9063, 'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.42,
              'Ag': 107.8682, 'Cd': 112.411, 'In': 114.818, 'Sn': 118.71, 'Sb': 121.760, 'Te': 127.6,
              'I': 126.90447, 'Xe': 131.293, 'Cs': 132.9054519, 'Ba': 137.327, 'La': 138.90547, 'Ce': 140.116,
              'Pr': 140.90465, 'Nd': 144.242, 'Pm': 146.9151, 'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25,
              'Tb': 158.92535, 'Dy': 162.5, 'Ho': 164.93032, 'Er': 167.259, 'Tm': 168.93421, 'Yb': 173.04,
              'Lu': 174.967, 'Hf': 178.49, 'Ta': 180.9479, 'W': 183.84, 'Re': 186.207, 'Os': 190.23, 'Ir': 192.217,
              'Pt': 195.084, 'Au': 196.966569, 'Hg': 200.59, 'Tl': 204.3833, 'Pb': 207.2, 'Bi': 208.9804,
              'Po': 208.9824, 'At': 209.9871, 'Rn': 222.0176, 'Fr': 223.0197, 'Ra': 226.0254, 'Ac': 227.0278,
              'Th': 232.03806, 'Pa': 231.03588, 'U': 238.02891, 'Np': 237.0482, 'Pu': 244.0642, 'Am': 243.0614,
              'Cm': 247.0703, 'Bk': 247.0703, 'Cf': 251.0796, 'Es': 252.0829, 'Fm': 257.0951, 'Md': 258.0951,
              'No': 259.1009, 'Lr': 262, 'Rf': 267, 'Db': 268, 'Sg': 271, 'Bh': 270, 'Hs': 269, 'Mt': 278,
              'Ds': 281, 'Rg': 281, 'Cn': 285, 'Nh': 284, 'Fl': 289, 'Mc': 289, 'Lv': 292, 'Ts': 294, 'Og': 294,
              'ZERO': 0}

periodictable = MM_of_Elements.keys() #""" list of all elements from the periodic table"""


from math import atan2,degrees
#Label line with line2D label data
#get from https://github.com/cphyc/matplotlib-label-lines
def labelLine(line,x,label=None,align=True,**kwargs):

    ax = line.axes
    xdata = line.get_xdata()
    ydata = line.get_ydata()

    if (x < xdata[0]) or (x > xdata[-1]):
        print('x label location is outside data range!')
        return

    #Find corresponding y co-ordinate and angle of the line
    ip = 1
    for i in range(len(xdata)):
        if x < xdata[i]:
            ip = i
            break

    y = ydata[ip-1] + (ydata[ip]-ydata[ip-1])*(x-xdata[ip-1])/(xdata[ip]-xdata[ip-1])

    if not label:
        label = line.get_label()

    if align:
        #Compute the slope
        dx = xdata[ip] - xdata[ip-1]
        dy = ydata[ip] - ydata[ip-1]
        ang = degrees(atan2(dy,dx))

        #Transform to screen co-ordinates
        pt = np.array([x,y]).reshape((1,2))
        trans_angle = ax.transData.transform_angles(np.array((ang,)),pt)[0]

    else:
        trans_angle = 0

    #Set a bunch of keyword arguments
    if 'color' not in kwargs:
        kwargs['color'] = line.get_color()

    if ('horizontalalignment' not in kwargs) and ('ha' not in kwargs):
        kwargs['ha'] = 'center'

    if ('verticalalignment' not in kwargs) and ('va' not in kwargs):
        kwargs['va'] = 'center'

    if 'backgroundcolor' not in kwargs:
        kwargs['backgroundcolor'] = ax.get_facecolor()

    if 'clip_on' not in kwargs:
        kwargs['clip_on'] = True

    if 'zorder' not in kwargs:
        kwargs['zorder'] = 2.5

    ax.text(x,y,label,rotation=trans_angle,**kwargs)


#get from https://github.com/cphyc/matplotlib-label-lines
def labelLines(lines,align=True,xvals=None,**kwargs):

    ax = lines[0].axes
    labLines = []
    labels = []

    #Take only the lines which have labels other than the default ones
    for line in lines:
        label = line.get_label()
        if "_line" not in label:
            labLines.append(line)
            labels.append(label)

    if xvals is None:
        xmin,xmax = ax.get_xlim()
        xvals = np.linspace(xmin,xmax,len(labLines)+2)[1:-1]

    for line,x,label in zip(labLines,xvals,labels):
        labelLine(line,x,label,align,**kwargs)


"""SGTE fitting using 
T - temperature
a - fitting parameters
"""
def SGTE(T,a):
  fval = a[0]+a[1]*T
  if len(a) > 2:
    fval += a[2]*T*np.log(T)
  if len(a) > 3:
    fval += a[3]*T*T
  if len(a) > 4:
    fval += a[4]*T*T*T
  if len(a) > 5:
    fval += a[5]/T
  return(fval)


"""SGTE fitting with two parameters"""
def SGTE2(T, a, b):
  return (SGTE(T, [a,b]))


"""SGTE fitting with three parameters"""
def SGTE3(T, a, b, c):
  return (SGTE(T, [a,b,c]))


"""SGTE fitting with four parameters"""
def SGTE4(T, a, b, c, d):
  return (SGTE(T, [a,b,c,d]))


"""SGTE fitting with five parameters"""
def SGTE5(T, a, b, c, d, e):
  return (SGTE(T, [a,b,c,d,e]))


"""SGTE fitting with six parameters"""
def SGTE6(T, a, b, c, d, e, f):
  return (SGTE(T, [a,b,c,d,e,f]))


"""SGTE fitting for heat capacity one parameter"""
def SGTEC1(T,a):
  return C_SGTE(T,[a])


"""SGTE fitting for heat capacity two parameters"""
def SGTEC2(T,a,b):
  return C_SGTE(T,[a,b])


"""SGTE fitting for heat capacity three parameters"""
def SGTEC3(T,a,b,c):
  return C_SGTE(T,[a,b,c])


"""SGTE fitting for heat capacity four parameters"""
def SGTEC4(T,a,b,c,d):
  return C_SGTE(T,[a,b,c,d])


"""SGTE fitting for heat capacity"""
def C_SGTE(T,a):
  fval = 0
  if len(a) > 0:
    fval += a[0]
  if len(a) > 1:
    fval += a[1]*T
  if len(a) > 2:
    fval += a[2]*T*T
  if len(a) > 3:
    fval += a[3]/T/T
  return(fval)


def SGTES(T,f):
  s = 0.0
  if len(f)>0:
    s += f[0]
  if len(f)>1:
    s += f[1]*np.log(T)
  if len(f)>2:
    s += f[2]*T
  if len(f)>3:
    s += f[3]*T*T
  if len(f)>4:
    s += f[4]/T/T
  return s

def SGTEH(T,f):
  h = 0.0
  if len(f)>0:
    h += f[0]
  if len(f)>1:
    h += f[1]*T
  if len(f)>2:
    h += f[2]*T*T
  if len(f)>3:
    h += f[3]*T*T*T
  if len(f)>4:
    h += f[4]/T
  return h

def SGTEC(T,f):
  s = 0.0
  if len(f)>0:
    s += f[0]
  if len(f)>1:
    s += f[1]*T
  if len(f)>2:
    s += f[2]*T*T
  if len(f)>3:
    s += f[3]/T/T
  return s


def CSGTEfit(f, x, y):
  popt,pcov = curve_fit(f, x, y)
  z = C_SGTE(x,popt)
  ferror=math.sqrt(((z-y)**2).sum()/len(z))
  return(popt,ferror)


def fitStoichiometricCp(x,y, thr=0.001):
  f,ferror = CSGTEfit(SGTEC2, x, y)
  if ferror > thr:
    f,ferror = CSGTEfit(SGTEC3, x, y)
  if ferror > thr:
    f,ferror = CSGTEfit(SGTEC4, x, y)
  return f,ferror


def H_SGTE(T,c):
  h = 0.
  if len(c)>0:
    h += c[0]*T
  if len(c)>1:
    h += c[1]/2*T*T
  if len(c)>2:
    h += c[2]/3*T*T*T
  if len(c)>3:
    h += -c[3]/T
  return h

def fitStoichiometricH(x,y,c):
  zz = H_SGTE(x,c)
  h = (y - zz).sum()/len(y)
  ferror=math.sqrt(((h+zz-y)**2).sum()/len(zz))
  h = [h]
  if len(c)>0:
    h.append(c[0])
  if len(c)>1:
    h.append(c[1]/2)
  if len(c)>2:
    h.append(c[2]/3)
  if len(c)>3:
    h.append(-c[3])
  return h,ferror


def S_SGTE(T,c):
  s = 0.
  if len(c)>0:
    s += c[0]+c[0]*np.log(T)
  if len(c)>1:
    s += c[1]*T
  if len(c)>2:
    s += c[2]/2*T*T
  if len(c)>3:
    s += -c[3]/2/T/T
  return s

def fitStoichiometricS(x,y,c):
  zz = S_SGTE(x,c)
  b = (y - zz).sum()/len(y)
  ferror=math.sqrt(((b+zz-y)**2).sum()/len(zz))
  s = []
  if len(c)>0:
    s.append(b+c[0])
    s.append(c[0])
  if len(c)>1:
    s.append(c[1])
  if len(c)>2:
    s.append(c[2]/2)
  if len(c)>3:
    s.append(-c[3]/2)
  return s,ferror

def fitStoichiometric(x,y, thr=1.0):
  f,ferror = SGTEfit(SGTE2, x, y)
  if ferror > thr:
    f,ferror = SGTEfit(SGTE3, x, y)
  if ferror > thr:
    f,ferror = SGTEfit(SGTE4, x, y)
  if ferror > thr:
    f,ferror = SGTEfit(SGTE5, x, y)
  if ferror > thr:
    f,ferror = SGTEfit(SGTE6, x, y)
  return f,ferror

def SGTEfit(f, x, y):
  popt,pcov = curve_fit(f, x, y)
  z = SGTE(x,popt)
  ferror=math.sqrt(((z-y)**2).sum()/len(z))
  return(popt,ferror)

def outexpressionG(f0):
    out = ""
    for i,f in enumerate(f0):
      if i==0:
        out += ' {:+g}'.format(f)
      elif i==1:
        out += ' {:+g}*T'.format(f)
      elif i==2:
        out += ' {:+g}*T*log(T)'.format(f)
      elif i==3:
        out += ' {:+g}*T*T'.format(f)
      elif i==4:
        out += ' {:+g}*T*T*T'.format(f)
      elif i==5:
        out += ' {:+g}/T'.format(f)
    return out 

def outexpressionS(f0):
    out = ""
    for i,f in enumerate(f0):
      if i==0:
        out += ' {:+g}'.format(f)
      elif i==1:
        out += ' {:+g}*log(T)'.format(f)
      elif i==2:
        out += ' {:+g}*T'.format(f)
      elif i==3:
        out += ' {:+g}*T*T'.format(f)
      elif i==4:
        out += ' {:+g}/T/T'.format(f)
    return out

def outexpressionH(f0):
    out = ""
    for i,f in enumerate(f0):
      if i==0:
        out += ' {:+g}'.format(f)
      elif i==1:
        out += ' {:+g}*T'.format(f)
      elif i==2:
        out += ' {:+g}*T*T'.format(f)
      elif i==3:
        out += ' {:+g}*T*T*T'.format(f)
      elif i==4:
        out += ' {:+g}/T'.format(f)
    return out

def outexpressionCp(f0):
    out = ""
    for i,f in enumerate(f0):
      if i==0:
        out += ' {:+g}'.format(f)
      elif i==1:
        out += ' {:+g}*T'.format(f)
      elif i==2:
        out += ' {:+g}*T*T'.format(f)
      elif i==3:
        out += ' {:+g}/T/T'.format(f)
    return out


def proStoichiometricG():
    #try:
      x = zthermo.get("temperature (K)")
      y = zthermo.get("Gibbs energy (eV/atom)")
      H298 = threcord.get("H298.15 (J/mol-atom)")
      x = np.array(list(map(float, x)))
      y = np.array(list(map(float, y)))*eVtoJ - H298
      i0 = 0
      for i,T in enumerate(x):
        if T < T0:
          i0 = i
      ifit0 = i0-15
      ifit0 = max(ifit0,0)
      
      f,ferror = fitStoichiometric(x[ifit0:],y[ifit0:])
      gout = 'G(T) =' + outexpressionG(f)
      #print(gout)
      s = []
      h = []
      c = []
      if len(f) >0:
        h.append(f[0])
      if len(f) >1:
        s.append(-f[1])
      if len(f) >2:
        s = []
        s.append(-f[1]-f[2])
        s.append(-f[2])
        h.append(-f[2])
        c.append(-f[2])
      if len(f) >3:
        s.append(-2.0*f[3])
        h.append(-f[3])
        c.append(-2.0*f[3])
      if len(f) >4:
        s.append(-3.0*f[4])
        h.append(-2.0*f[4])
        c.append(-6.0*f[4])
      if len(f) >5:
        s.append(f[5])
        h.append(2.0*f[5])
        c.append(-2.0*f[5])
      sout = 'S(T) =' + outexpressionS(s)
      hout = 'H(T) =' + outexpressionH(h)
      cout = 'Cp(T) =' + outexpressionCp(c)
      """
      print (sout)
      print (hout)
      print (cout)
      """
      uncertanty = {}
      SGTErec.update({"G-H298.15 (J/mol-atom)":gout})
      SGTErec.update({"H-H298.15 (J/mol-atom)":hout})
      SGTErec.update({"S (J/mol-atom/K)":sout})
      SGTErec.update({"Cp (J/mol-atom/K)":cout})
      SGTErec.update({"fitting uncertainty":round(ferror,1)})
      return(f,h,s,c,x[i0:])


def proStoichiometricCp():
    #try:
      uncertanty = {}
      x = zthermo.get("temperature (K)")
      y = zthermo.get("Cp (J/mol-atom/K)")
      H298 = threcord.get("H298.15 (J/mol-atom)")
      x = np.array(list(map(float, x)))
      y = np.array(list(map(float, y)))
      i0 = 0
      for i,T in enumerate(x):
        if T < T0:
          i0 = i
      ifit0 = i0
      ifit0 = max(ifit0,0)
      #print("xxxxxxxx=",x[ifit0:],T0) 
      c,cerror = fitStoichiometricCp(x[ifit0:],y[ifit0:])

      y = zthermo.get("enthalpy (J/mol-atom)")
      y = np.array(list(map(float, y))) - H298
      h,herror = fitStoichiometricH(x[ifit0:],y[ifit0:],c)

      y = zthermo.get("entropy (J/mol-atom/K)")
      y = np.array(list(map(float, y)))
      s,serror = fitStoichiometricS(x[ifit0:],y[ifit0:],c)

      f = [h[0]]
      if len(s) >0:
        f.append(-s[0]+c[0])
        f.append(-c[0])
      if len(c) >1:
        f.append(-c[1]/2)
      if len(c) >2:
        f.append(-c[2]/6)
      if len(c) >3:
        f.append(-c[3]/2)
      gout = 'G(T) =' + outexpressionG(f)
      #print (gout)

      SGTErec.update({"T":[x[ifit0], x[-1]]})
      SGTErec.update({"Cp (J/mol-atom/K)":[outexpressionCp(c),{"error":round(cerror,2)}]})
      SGTErec.update({"H-H298.15 (J/mol-atom)":[outexpressionH(h),{"error":round(herror,2)}]})
      SGTErec.update({"S (J/mol-atom/K)":[outexpressionS(s),{"error":round(serror,2)}]})
      SGTErec.update({"G-H298.15 (J/mol-atom)":[outexpressionG(f),{"error":round(herror,2)}]})
      return(f,h,s,c,x[i0:])


class thermoplot:
    def __init__(self, folder,thermodynamicproperty,x,y,reflin=None, yzero=None,fitted=None,xT=None,xlabel="T (K)", lp=False,
        ylabel=None, ytext=None, xlim=None, elonly=None, expt=None, CoT=False, label=None, single=False, plottitle=None):

        plt.rc('font', size=24)
        self.fig,self.ax=plt.subplots()
        self.fig.set_size_inches(12,9)
        self.ax.yaxis.set_ticks_position('both')
        self.cwd = os.getcwd()
        os.chdir( folder )

        self.folder = folder
        self.thermodynamicproperty = thermodynamicproperty
        self.x = np.array(x)
        self.y = np.array(y)
        self.reflin = reflin
        self.yzero = yzero
        self.fitted = fitted
        self.xT = xT
        self.ytext = ytext
        self.elonly = elonly
        self.expt = expt
        self.CoT = CoT
        self.single = single
        self.plottitle = plottitle

        self._xlabel = xlabel
        self.lp = lp
        self._ylabel = thermodynamicproperty
        if ylabel!=None: self._ylabel = ylabel
        self._label = self.thermodynamicproperty
        if label!=None: self._label = label
        self.fname = self.thermodynamicproperty.split('(')[0].strip().replace(' ','_')+".png"

        self.ax.set_xlim([0,np.array(list(map(float,x))).max()])
        self.plot_xlim = np.array(list(map(float,x))).max()
        self.xlim = xlim
        if xlim!=None:
            try: self.ax.set_xlim([0,xlim])
            except: self.ax.set_xlim(xlim)

        if self.thermodynamicproperty=="0 K total energies (eV/atom)": self.plot_EV()
        elif self.thermodynamicproperty=="Helmholtz energy (eV/atom)": self.plot_Helmholtz_energy_v0()
        elif self.thermodynamicproperty=="Helmholtz energy analysis (eV/atom)": self.plot_Helmholtz_energy_v1()
        elif self.thermodynamicproperty.lower()=="Effective charge carrier concentration ($e/cm^{3}$)".lower():
            self.plot_Effective_charge_carrier_concentration()
        elif self.thermodynamicproperty.lower()=="Electron DOS (States/Atom/eV)".lower(): self.plot_Electron_DOS()
        elif self.thermodynamicproperty.lower()=="Bulk modulus (GPa)".lower(): self.plot_Bulk_modulus()
        elif self.thermodynamicproperty.lower()=="LTC analysis (1/K)".lower(): self.plot_LTC_analysis()
        elif self.thermodynamicproperty=="Gamma point phonons": self.plot_Gamma_point_phonons()
        elif self.thermodynamicproperty.lower()!="heat capacities (J/mol-atom/K)".lower(): self.plot_default()
        else: self.plot_Heat_Capacity()

        if self.plottitle!=None: plt.title(self.plottitle)
        plt.xlabel(self._xlabel)
        plt.ylabel(self._ylabel)
        self.ax.legend(loc=0, prop={'size': 24})
        self.fig.savefig(self.fname,bbox_inches='tight')
        plt.close(self.fig)

        os.chdir( self.cwd )
        figures.update({self.thermodynamicproperty:folder.split('/')[-1]+'/'+self.fname})
    

    def plot_EV(self):
        self._xlabel = "atomic volume ($\AA^3$)"
        self._ylabel = "0 K total energies (eV/atom)"
        self.ax.set_xlim([min(self.x)*0.95,max(self.x)*1.05])
        if self.lp:
            self.ax.plot(self.x, self.y, marker='o', markersize=12, 
                color='r', linestyle='-', label=self._label)
        else:
            self.ax.plot(self.x, self.y, fillstyle='none', marker='o', markersize=12, 
                color='k', linestyle='None', label=self._label)
        xnew = np.linspace(min(self.x)*0.95,max(self.x)*1.05, 300)  
        from dfttk.pythelec import BMvol4, BMvol, alt_curve_fit
        f2, pcov = alt_curve_fit(BMvol4, self.x, self.y)
        ynew = BMvol(xnew, f2)
        self.ax.plot(xnew,ynew,'-',linewidth=1,color='b', label="BMvol4")


    def plot_Helmholtz_energy_v0(self):
        self.fig.set_size_inches(12,11)
        self._xlabel = "Atomic volume ($\AA^3$)"
        plt.ylabel("Helmholtz energy (eV/atom)")
        self.ax.plot(self.x, self.y, marker='o', markersize=4, color='k', linestyle=':')
        self.ax.set_xlim([min(self.x)*0.95,max(self.x)*1.05])
        fd = 0.05*(max(self.y)-min(self.y))
        self.ax.set_ylim([min(self.y)-fd,max(self.y)+2*fd])
        v,xx,t,o,f = self.reflin
        for i,T in enumerate(t):
            #self.ax.plot(v, o[i], fillstyle='none', marker='+', mew=2, markersize=12, color='b', linestyle='None')
            self.ax.plot(v, o[i], fillstyle='none', marker='+', markersize=12, color='b', linestyle='None')
            self.ax.plot(xx, f[i], color='b', linestyle='-')
        for i, l1 in enumerate(plt.gca().get_lines()):
            if i==0:
                x0 = 0.5*(min(self.x)+max(self.x))
                labelLine(l1,x0,label=r'$V_{eq}$',align = True)
            else:
                ii = (i-2)//4
                if (ii*4+2!=i): continue
                x0 = 0.90*0.95*min(self.x)+0.10*1.05*max(self.x)
                labelLine(l1,x0,label=r'${} K$'.format(int(t[ii*2])),align = True)
    

    def plot_Helmholtz_energy_v1(self):
        self.fig.set_size_inches(12,11)
        self._xlabel = "Atomic volume ($\AA^3$)"
        plt.ylabel("Helmholtz energy (eV/atom)")
        self.ax.plot(self.x, self.y, marker='o', markersize=4, color='k', linestyle=':')
        self.ax.set_xlim([min(self.x)*0.95,max(self.x)*1.05])
        fd = 0.05*(max(self.y)-min(self.y))
        self.ax.set_ylim([min(self.y)-fd,max(self.y)+2*fd])
        v,xx,t,o,f = self.reflin
        for i,T in enumerate(t):
            self.ax.plot(v, o[i], fillstyle='none', marker='+', mew=2, markersize=12, color='r', linestyle='None')
            #spl = make_interp_spline(v, o[i], k=3)  # type: BSpline
            spl = interp1d(v, o[i])
            power_smooth = spl(xx)
            self.ax.plot(xx, power_smooth, color='r', linestyle=':')
            self.ax.plot(xx, f[i], color='b', linestyle='-')
        for i, l1 in enumerate(plt.gca().get_lines()):
            if i==0:
                x0 = 0.5*(min(self.x)+max(self.x))
                labelLine(l1,x0,label=r'$V_{eq}$',align = True)
            else:
                ii = (i-3)//6
                if (ii*6+3!=i): continue
                x0 = 0.90*0.95*min(self.x)+0.10*1.05*max(self.x)
                labelLine(l1,x0,label=r'${} K$'.format(int(t[ii*2])),align = True)
    

    def plot_Effective_charge_carrier_concentration(self):
        self.ax.set_yscale('symlog')
        yy = self.y[self.x>0]
        xx = self.x[self.x>0]
        self.ax.plot(xx,yy,'-',linewidth=2,color='b', label=self._label)
        if self.xlim!=None:
            self.fname = self.thermodynamicproperty.split('(')[0].strip().replace(' ','_')+'_'+str(int(self.xlim))+".png"


    def plot_Electron_DOS(self):
        self.ax.plot(self.x,self.y,'-',linewidth=2,color='b', label=self._label)
        self.fname = self.thermodynamicproperty.split('(')[0].strip().replace(' ','_')+'_'+str(str(-self.xlim[0]))+"eV.png"


    def plot_Bulk_modulus(self):
        self.ax.plot(self.x,self.reflin,'--',linewidth=2,color='k', label=self._label+",$B_s$")
        self.ax.plot(self.x,self.y,'-',linewidth=2,color='b', label=self._label+",$B_T$")
        plot_expt(self.expt, 'bulk modulus', self.ax, xlim=self.xlim)


    def plot_LTC_analysis(self):
        self.ax.ticklabel_format(axis='y',style='sci',scilimits=(-2,4))
        self.ax.plot(self.x,self.y,'-',linewidth=2,color='b', label="dfttk")
        self.ax.plot(self.x,self.reflin,'--',linewidth=2,color='k', label="splev")


    def plot_Gamma_point_phonons(self):
        if self.reflin is not None:
            self.ax.plot(self.x,self.reflin,':',linewidth=1,color='k')
        self.ax.plot(self.x,self.y,'-',linewidth=2,color='b', label=self._label)
        self.ax.set_xlim([min(self.x)*1.05,max(self.x)*1.05])
        xx0 = np.array(self.ytext[0])
        yy0 = np.array(self.ytext[1])
        ss0 = self.ytext[2]
        for i in range (len(ss0)):
            self.ax.text(xx0[i], yy0[i], ss0[i], color='r', rotation=90, 
                horizontalalignment='left', verticalalignment='bottom')


    def plot_default(self):
        if self.thermodynamicproperty.split('(')[0].strip()=="Debye temperature":
            self.y=self.y[self.x>0]
            self.x=self.x[self.x>0]
        if self.yzero != None:
            y0 = np.nanmin(np.array(list(map(float,self.y))))
            y1 = np.nanmax(np.array(list(map(float,self.y))))
            self.ax.set_ylim([min(0.0, y0),y1*1.05])
            self.ax.ticklabel_format(axis='y',style='sci',scilimits=(-2,4))
        if self.reflin is not None:
            self.ax.plot(self.x,self.reflin,':',linewidth=1,color='k')
        self.ax.plot(self.x,self.y,'-',linewidth=2,color='b', label=self._label)
        if self.fitted!=None:
            self.ax.plot(self.xT[::5],self.fitted[::5],'--',fillstyle='none', marker='o', markersize=12, 
                linewidth=2,color='k', label="fitted")
        if self.xlim!=None: 
            self.ax.set_xlim([0.0,self.xlim])
            self.ax.set_ylim([0.98*min(self.y),1.02*max(self.y)])
            self.fname = self.thermodynamicproperty.split('(')[0].strip().replace(' ','_')+'_'+str(self.xlim)+".png"
        if self.thermodynamicproperty=="LTC (1/K)":
            plot_expt(self.expt, 'linear thermal expansion', self.ax, xlim=self.plot_xlim)
        elif self.thermodynamicproperty=="Entropy (J/mol-atom/K)":
            plot_expt(self.expt, 'entropy', self.ax, xlim=self.plot_xlim)
        elif self.thermodynamicproperty=="Enthalpy-H298 (J/mol-atom)":
            plot_expt(self.expt, 'enthalpy', self.ax, xlim=self.plot_xlim)
        elif self.thermodynamicproperty=="Lorenz number ($WΩK^{−2}$)":
            self.ax.set_ylim([min(2.e-8, np.array(self.y).min()),max(3.e-8,np.array(self.y).max())])


    def plot_Heat_Capacity(self):
        if self.fitted!=None:
            y = np.array(self.y)
            y0,y1,y2 = y[:,0], y[:,1], y[:,2]
            self.ax.set_ylim([0.0,np.array(list(map(float,y0))).max()*1.05])
            self.ax.plot(self.x,y0,'-',linewidth=2,color='b', label=self._label+",$C_p$")
            self.ax.plot(self.xT[::5],self.fitted[::5],'--',fillstyle='none', marker='o', markersize=12, 
                linewidth=2,color='k', label="fitted")
            self.ax.plot(self.x,y1,'--',linewidth=2,color='black', label="$C_v$")
            self.ax.plot(self.x,y2,':',linewidth=2,color='g', label="$C_{v,ion}$")
            y2 = np.array(list(map(float,y1))) - np.array(list(map(float,y2)))
            self.ax.plot(self.x,y2,'-.',linewidth=2,color='r', label="$C_{el}$")
            self.fname = self.thermodynamicproperty.split('(')[0].strip().replace(' ','_')+"_fitted.png"
        else:
            x = np.array(self.x)
            y = np.array(self.y)
            y0,y1 = y[:,0], y[:,1]
            y2 = y0 - y1
            if self.CoT:
                y0 = y0[x>0]
                y = y2[x>0]
                x = x[x>0]
                if self.elonly!=None:
                    self._xlabel = "$T (K)$"
                    self._ylabel = "$C_el/T$ (J/mol-atom/K/K)"
                    y = y[x<=self.elonly*1.2]
                    x = x[x<=self.elonly*1.2]
                    self.ax.set_xlim([0.0,self.elonly])
                    self.ax.plot(x,y/x,'-.',linewidth=2,color='r', label=self._label+",$C_{el}/T$")
                    ymax = plot_expt(self.expt, 'electronic heat capacity', self.ax, CoT=self.CoT, xlim=self.elonly)
                    self.ax.set_ylim([0.0,max(ymax,(y/x).max())*1.1])
                    self.fname = self.thermodynamicproperty.split('(')[0].strip().replace(' ','_')+\
                        '_'+str(self.elonly)+'_el_oT'".png"
                elif self.xlim!=None:
                    self._xlabel = "$T^2 (K^2)$"
                    self._ylabel = "$C/T$ (J/mol-atom/K/K)"
                    y = y0/x
                    x = x*x
                    y = y[x<self.xlim*1.2]
                    x = x[x<self.xlim*1.2]
                    self.ax.set_xlim([0.0,self.xlim])
                    if self.single:
                        self.ax.plot(x,y,'-',linewidth=2,color='b', label=self._label+",$C_{v,lat+el}/T$")
                    else:
                        self.ax.plot(x,y,'-',linewidth=2,color='b', label=self._label+",$C_{p,lat+el}/T$")
                    plot_expt(self.expt, 'heat capacity', self.ax, CoT=self.CoT, xlim=self.xlim)
                    self.fname = self.thermodynamicproperty.split('(')[0].strip().\
                        replace(' ','_')+'_'+str(self.xlim)+"_oT2.png"
            elif self.xlim!=None: 
                self.ax.set_xlim([0.0,self.xlim])
                if self.xlim==300: self.ax.set_ylim([0.0,25])
                y0 = y0[x<=self.xlim*1.1]
                y1 = y1[x<=self.xlim*1.1]
                y2 = y2[x<=self.xlim*1.1]
                x = x[x<=self.xlim*1.1]
                if self.single:
                    self.ax.plot(x,y0,'-',linewidth=2,color='b', label=self._label+",$C_{v,lat+el}$")
                    self.ax.plot(x,y1,'--',linewidth=2,color='black', label="$C_{v,lat}$")
                else:
                    self.ax.plot(x,y0,'-',linewidth=2,color='b', label=self._label+",$C_{p,lat+el}$")
                    self.ax.plot(x,y1,'--',linewidth=2,color='black', label="$C_{p,lat}$")
                plot_expt(self.expt, 'heat capacity', self.ax, xlim=self.xlim)
                """
                if y2.max() > 1.e-2:
                    self.ax.plot(x,y2,'-.',linewidth=2,color='r', label=self._label+",$C_{el}$")
                    plot_expt(self.expt, 'electronic heat capacity', self.ax, xlim=self.xlim)
                """
                self.fname = self.thermodynamicproperty.split('(')[0].strip().replace(' ','_')+'_'+str(self.xlim)+".png"
            elif self.elonly!=None:
                self.ax.set_xlim([0.0,self.elonly])
                y2 = y2[x<=self.xlim*1.1]
                x = x[x<=self.xlim*1.1]
                self.ax.plot(x,y2,'-.',linewidth=2,color='r', label=self._label+",$C_{el}$")
                plot_expt(self.expt, 'electronic heat capacity', self.ax, xlim=self.elonly)
                self.fname = self.thermodynamicproperty.split('(')[0].strip().replace(' ','_')+'_'+str(self.elonly)+'_el.png'
            else:
                if self.single:
                    self.ax.plot(x,y0,'-',linewidth=2,color='b', label=self._label+",$C_{v,lat+el}$")
                    self.ax.plot(x,y1,'--',linewidth=2,color='black', label="$C_{v,lat}$")
                else:
                    self.ax.plot(x,y0,'-',linewidth=2,color='b', label=self._label+",$C_{p,lat+el}$")
                    self.ax.plot(x,y1,'--',linewidth=2,color='black', label="$C_{p,lat}$")
                plot_expt(self.expt, 'heat capacity', self.ax, xlim=self.plot_xlim)
                """
                if y2.max() > 1.e-2:
                    self.ax.plot(x,y2,'-.',linewidth=2,color='r', label=self._label+",$C_{el}$")
                    plot_expt(self.expt, 'electronic heat capacity', self.ax, xlim=self.plot_xlim)
                """
                self.fname = self.thermodynamicproperty.split('(')[0].strip().replace(' ','_')+".png"

            plt.gca().set_ylim(bottom=0)


def plot_expt (expt, prp, ax, CoT=False, xlim=None):
    #global mindex
    mindex = 0
    ymax = 0.0
    if expt!=None:
        for rec in expt:
            if prp!=rec['property']: continue
            try:
                xval = np.array(rec['T'])
                yval = np.array(rec['val'])
            except:
                try:
                    lines = np.array(rec['data'])
                    xval = lines[0::2]
                    yval = lines[1::2]
                except:
                    continue
            Author = rec['Author']
            Unit = rec['Unit']
            natom = rec['natom']
            yval /= natom

            if Unit=='mJ/K' : yval /= 1000.
            elif Unit=='cal/K' : yval *= 4.184

            if CoT:
                if prp!='electronic heat capacity':
                    xx = xval[xval>0]*xval[xval>0]
                    yy = yval[xval>0]/xval[xval>0]
                    yy = yy[xx<xlim]
                    xx = xx[xx<xlim]
                else:
                    xx = xval[xval>0]
                    yy = yval[xval>0]/xval[xval>0]
                    yy = yy[xx<xlim]
                    xx = xx[xx<xlim]
            elif xlim!=None:
                yy = yval[xval<xlim]
                xx = xval[xval<xlim]
            else:
                yy = yval
                xx = xval

            if len(xx)>0:
                if Author.startswith('Andersson(CALPHAD)'):
                    ax.plot(xx,yy, marker='o', color='b',markersize=10, mew=2,
                        linestyle='None', fillstyle='none', label=Author.split(',')[0])
                elif Author.startswith('Andersson(CALPHAD)') or \
                    Author.startswith('Chase(JANAF)'):
                    ax.plot(xx,yy, marker='s', color='k',markersize=8, mew=2,
                        linestyle='None', fillstyle='none', label=Author.split(',')[0])
                    """
                    try:
                        ax.plot(xx,yy, marker=markers[mindex%len(markers)], markersize=8, mew=2,
                            linestyle='None', fillstyle='none', label=Author.split(',')[0])
                    except:
                        ax.plot(xx,yy, marker=markers[mindex%len(markers)], markersize=8,
                            linestyle='None', label=Author.split(',')[0])
                    """
                else:
                    ax.plot(xx,yy, marker=markers[mindex%len(markers)], markersize=8,
                        linestyle='None', label=Author.split(',')[0])
                ymax = max(yy.max(), ymax)
            mindex += 1
    return ymax

def myjsonout(data,fp,indent="",comma=""):
	#print (data)
	mj = ''
	if (isinstance(data,dict)):
		fp.write('{}\n'.format('{'))
			#sys.stdout.write('\n{}{}\n'.format(indent, '{'))
		nkey = 0
		for key in sorted(set(data.keys())):
			nkey += 1
			if nkey!=len(data):
				comma1 = ","
			else:
				comma1 = ""
			val = data[key]
			jval = json.dumps(val)
			jkey = json.dumps(key)
			#print (val)
			if (isinstance(val,dict)):
				fp.write('{}{}: '.format(indent+"    ",jkey))
				myjsonout(val,fp,indent+"    ",comma1)
			elif (isinstance(val,tuple)):
				#print (val)
				out = list(val)
				#print(out)
				fp.write('{}{}: {}{}\n'.format(indent + "    ", jkey, out, comma1))
			elif (isinstance(val,str)):
				if (indent == ""):
					fp.write('{}{}: {}{}\n'.format(indent + "    ", jkey, jval, comma1))
				else:
					fp.write('{}{}: {}{}\n'.format(indent + "    ", jkey, jval, comma1))
			else:
				if (indent==""):
					fp.write('{}{}: {}{}\n'.format(indent + "    ", jkey, jval, comma1))
				else:
					fp.write('{}{}: {}{}\n'.format(indent + "    ", jkey, jval, comma1))

				#print(val)
				"""
				if (nkey!=len(data)):
					sys.stdout.write('{}{}: {},\n'.format(indent+"    ", key, val))
				else:
					sys.stdout.write('{}{}: {}\n'.format(indent+"    ", key, val))
				"""
		if comma==',':
			fp.write('{}{}{}\n\n'.format(indent,'}', comma))
		else:
			fp.write('{}{}{}\n'.format(indent, '}', comma))


def Myjsonout(data,out):
    if (isinstance(data,dict)):
        myjsonout(data, out, indent="", comma="")
    elif (isinstance(data,list)):
        out.write("[\n")
        for j,rec in enumerate(data):
            if j!=len(data)-1: myjsonout(rec, out, indent="", comma=",")
            else: myjsonout(rec, out, indent="", comma="")
        out.write("]\n")


def similar(pp,pall):
  known = ["L12", "delta", "D022", "Gamma"]
  ii = -1
  for o in known:
    if pp.find(o)>-1:
      pname = o
      ii = 0
      break
  if ii == -1:
    return "unknown"

  s = 0.0
  for i,p in enumerate(pall):
    snew = SequenceMatcher ( None, pname, p ).ratio()
    if snew > s:
      ii = i
      s = snew
  print (pp, "= ", pall[ii], " by ", s)
  if s > 0.5:
    return pall[ii]
  else:
    return "unknown"


def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False


def formula2elist(formula):
  formula = formula.replace(" ",'').replace("-",'').replace(",",'')
  newc = ""
  """Follow the convention, elemental symbol must start from capital letter"""
  for c in formula:
    if c in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
      newc = newc + '|'
    newc = newc + c
  els = newc.split('|')
  els = [k for k in els if k != '']

  """now get the composition for each element"""
  ele = []
  com = []
  for el in els:
    newel = ""
    newcc = ""
    for c in el:
      if c.isalpha():
        newel = newel + c
      else:
        newcc = newcc + c

    if (newel not in periodictable):
      raise ValueError('"'+newel+'" is not an element! your formula is wrong!')
    ele.append(newel)

    if (len(newcc)!=0):
      if (isfloat(newcc)):
        com.append(float(newcc))
      else:
        raise ValueError('"'+newcc+'" is not an element! your formula is wrong!')
    else:
      com.append(1.0)
  com = np.array(list(map(float,com)))
  com = com/sum(com)
  #sorted the sequence and merge the duplicate
  elist = sorted(set(ele))
  clist = np.zeros(len(elist), dtype=float)
  for j,el in enumerate(ele):
    ix = elist.index(el)
    clist[ix] += com[j]

  return elist

def prety_formulaO(longphasename):
  puc = longphasename.split('|')[-1]
  _els,_nat=formula2composition(puc)

def prety_formula(_els,_nat):
  els = sorted(set(_els))
  nat = np.zeros(len(els),dtype=int)
  for i,el in enumerate(_els):
    ix = els.index(el)
    nat[ix] += _nat[i]

  Nd = min(nat)
  for i in range(Nd,0,-1):
    out = True
    for j in range(len(nat)):
      if ((nat[j]//i)*i!=nat[j]):
        out = False
        break
    if out:
      break
  form = ""
  for j,el in enumerate(els):
    ix = nat[j]//i
    form = form+el
    if ix!=1:
      form = form+str(ix)
  return form


def Genergy(thermofile,dir0):
  tmelt = 9999.
  ele = threcord.get("Elements")
  if ele!=None:
    if len(ele)==1:
      tmelt = ELEMENTS[ele[0]].tmelt

  folder = dir0+'/'+"figures"
  if not os.path.exists(folder):
    os.mkdir(folder)

  tmp = [s for s in thermofile.split('/') if s!=""]
  tmp[-1] = 'vdos_Cij'
  vdos_e_Cij = '/'.join(tmp)
  #print("Cij",vdos_e_Cij)
  if os.path.exists(vdos_e_Cij) :
    vdos_e_Cij = np.loadtxt(vdos_e_Cij, comments="#", dtype=np.float)
    ij = 0
    for i in range(1,7):
      for j in range(i,7):
        ij = ij + 2
        if abs(vdos_e_Cij[:,ij]).max() > 1.0:
          thermoplot(folder,"C"+"_"+str(i)+"_"+str(j),list(vdos_e_Cij[:,0]),list(vdos_e_Cij[:,ij]),yzero=0.0)
    thermoplot(folder,"B_v",list(vdos_e_Cij[:,0]),list(vdos_e_Cij[:,43]),yzero=0.0)
    thermoplot(folder,"G_v",list(vdos_e_Cij[:,0]),list(vdos_e_Cij[:,44]),yzero=0.0)
    thermoplot(folder,"E_v",list(vdos_e_Cij[:,0]),list(vdos_e_Cij[:,45]),yzero=0.0)
    sys.exit()


  thermo = np.loadtxt(thermofile, comments="#", dtype=np.float)
  thermo[np.isnan(thermo)] = 0.0
  for i,cp in enumerate(thermo[:,6]):
    if cp > CpMax: break
    elif thermo[i,0] > tmelt: break

  thermo = thermo[0:i,:]
    
  Vstack=interpolate.splrep(thermo[:,0], thermo[:,1])
  V298 = float(interpolate.splev(T0, Vstack))
  Hstack=interpolate.splrep(thermo[:,0], thermo[:,4])
  H298 = float(interpolate.splev(T0, Hstack))
  threcord.update({"H298.15 (J/mol-atom)":round(H298,4)})
  Sstack=interpolate.splrep(thermo[:,0], thermo[:,3])
  S298 = float(interpolate.splev(T0, Sstack))
  threcord.update({"S298.15 (J/mol-atom/K)":round(S298,6)})

  zthermo.update({"temperature (K)":list(thermo[:,0])})
  zthermo.update({"atomic volume ($\AA^3$)":list(thermo[:,1])})
  thermoplot(folder,"atomic volume ($\AA^3$)",list(thermo[:,0]),list(thermo[:,1]))
  zthermo.update({"Gibbs energy (eV/atom)":list(thermo[:,2])})
  zthermo.update({"enthalpy (J/mol-atom)":list(thermo[:,4])})
  zthermo.update({"entropy (J/mol-atom/K)":list(thermo[:,3])})
  zthermo.update({"Cp (J/mol-atom/K)":list(thermo[:,6])})

  if fitCp:
    g,h,s,c,x=proStoichiometricCp()
  else:
    g,h,s,c,x=proStoichiometricG()

  threcord.update({"SGTE fitting":SGTErec})
  thermoplot(folder,"Gibbs energy-H298 (J/mol-atom)",list(thermo[:,0]),list(thermo[:,2]*eVtoJ-H298),fitted=list(SGTE(x,g)), xT=list(x))
  thermoplot(folder,"enthalpy-H298 (J/mol-atom)",list(thermo[:,0]),list(thermo[:,4]-H298), fitted=list(SGTEH(x,h)), xT=list(x))
  #thermoplot(folder,"enthalpy-H298 (J/mol-atom)",list(thermo[:,0]),list(thermo[:,4]-H298), fitted=list(SGTE(x,g)+x*SGTES(x,s)), xT=list(x))
  thermoplot(folder,"entropy (J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,3]),yzero=0.0, fitted=list(SGTES(x,s)), xT=list(x))

  zthermo.update({"LTC (1/K)":list(thermo[:,5])})
  thermoplot(folder,"LTC (1/K)",list(thermo[:,0]),list(thermo[:,5]),yzero=0.0)
  zthermo.update({"Cv (J/mol-atom/K)":list(thermo[:,14])})
  zthermo.update({"Cv,ion (J/mol-atom/K)":list(thermo[:,7])})
  Cele = [round(c,6) for c in thermo[:,14]-thermo[:,7]]
  zthermo.update({"Cele (J/mol-atom/K)":Cele})
  ncols = [6,14,7]
  thermoplot(folder,"heat capacities ((J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]),fitted=list(SGTEC(x,c)), xT=list(x))
  ncols = [6,8]
  thermoplot(folder,"heat capacities ((J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), xT=list(x), expt=expt)
  thermoplot(folder,"heat capacities ((J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), xT=list(x), xlim=300,expt=expt)
  thermoplot(folder,"heat capacities ((J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), xT=list(x), xlim=70,expt=expt)
  thermoplot(folder,"heat capacities ((J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), xT=list(x), xlim=100,expt=expt, CoT=True)
  thermoplot(folder,"heat capacities ((J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), xT=list(x), xlim=1000,expt=expt, CoT=True)
  thermoplot(folder,"heat capacities ((J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), xT=list(x), xlim=10000,expt=expt, CoT=True)
  thermoplot(folder,"heat capacities ((J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), xT=list(x), elonly=300, expt=expt)
  thermoplot(folder,"heat capacities ((J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), xT=list(x), elonly=300, expt=expt, CoT=True)
  thermoplot(folder,"heat capacities ((J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), xT=list(x), elonly=70, expt=expt)
  thermoplot(folder,"heat capacities ((J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), xT=list(x), elonly=70, expt=expt, CoT=True)
  zthermo.update({"Debye temperature (K)":list(thermo[:,13])})
  thermoplot(folder,"Debye temperature (K)",list(thermo[:,0]),list(thermo[:,13]),yzero=0.0)
  thermoplot(folder,"Debye temperature (K)",list(thermo[:,0]),list(thermo[:,13]),yzero=0.0, xlim=70)
  zthermo.update({"bulk modulus (GPa)":list(thermo[:,15])})
  thermoplot(folder,"bulk modulus (GPa)",list(thermo[:,0]),list(thermo[:,15]),yzero=0.0)

  threcord.update({"zthermodynamic properies":zthermo})
  threcord.update({"Atomic volume at 298.15 K ($\AA^3$)":round(V298,6)})

  with open(vdos_e.replace('thermo/vdos_e','tplate/POSCAR'), 'r') as f:
    vvv = f.readlines()
  natom = sum([int(vv) for vv in vvv[6].split(' ') if vv!=""])
  structure.update({"number of atoms in POSCAR":natom})

  #natom = threcord.get("number of atoms in the primitive unit cell")

  with open(vdos_e.replace('thermo/vdos_e','thermo/data.in'), 'r') as f:
    vvv = f.readlines()
    Vfiles = []
    Pfiles = []
    volumes = []
    energies = []
    for vv in vvv[1:]:
      v = vv.split(' ')
      Pfiles.append("phonon/"+v[2].split('/')[-1].replace('\n', '').replace('"', ''))
      Vfiles.append(v[2].split('/')[-1].replace('\n', '').replace('"', ''))
      volumes.append(round(float(v[0])/natom,6))
      energies.append(round(float(v[1])/natom,6))
  structure.update({"Static vasp settings":Vfiles})
  structure.update({"phonon vasp settings and force constants":Pfiles})
  threcord.update({"volumes":volumes})
  threcord.update({"energies":energies})

  with open(dir0+'/E-V.dat','w') as f:
    for i,v in enumerate(volumes):
      f.write('{} {}\n'.format(v,energies[i]))
  #cmd = "YWfit -BMvol <"+dir0+'/E-V.dat | grep "f_expr(x) = "'
  ffun = "-Morse"
  cmd = "YWfit "+ffun+" <"+dir0+'/E-V.dat | grep "f_expr(x) = "'
  output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                      universal_newlines=True)
  cwd = os.getcwd()
  os.chdir( dir0+"/figures")
  fitF = output.stdout
  with open('E-V.plt','w') as f:
    f.write('set terminal postscript landscape enhanced color "Times_Roman" 20\n')
    f.write('set encoding iso_8859_1\n')
    f.write('set pointsize 1.2\n')
    f.write('set size 0.95,0.95\n')
    f.write('set output "E-V.eps"\n')
    f.write('{}\n'.format(fitF))
    f.write('set key right bottom\n')
    f.write('set xlabel "atomic volume ($\AA^3$)\n')
    f.write('set ylabel "static energy (eV/atom)\n')
    f.write('plot "../E-V.dat" title "calculated" w p pt 7, \\\n')
    f.write('     f_expr(x) title "'+ffun+'" w l lt -1\n')
  #cmd = "gnuplot E-V.plt; convert -fuzz 100% -transparent white -rotate 90 -density 120x120 E-V.eps E-V.png"
  cmd = "gnuplot E-V.plt; convert -background white -alpha remove -rotate 90 -density 120x120 E-V.eps E-V.png"
  output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                      universal_newlines=True)

  figures.update({"static E-V curve": "figures/E-V.png"})
  threcord.update({"figures":figures})
  os.chdir(cwd)
  return Vfiles,Pfiles,g

def BMfitP(x,z):
  p = 0.0
  N = len(z)
  for n in range(N):
    p = p + float(N-n-1)*z[n]*x**(N-n+0.5)*(-2.0/3.0)
  return (-p)

def BMfit(v,p,g, T):
  global tPmax
  v = np.array(list(map(float,v)))
  p = np.array(list(map(float,p)))
  g = np.array(list(map(float,g)))
  f = g - p*v
  x = v**(-2.0/3.0)
  z = np.polyfit(x,f,BMvol)
  gf = np.poly1d(z)
  if (Debug==1):
    for i, vv in enumerate(v):
      print(vv, p[i], BMfitP(x[i],z), f[i], gf(x[i]))
   
  xx = x[0]
  xd = (x[len(x)-1] - x[0])*0.02
  pp = []
  gg = []
  vv = []
  for i in range(999):
    ppxx = BMfitP(xx,z)
    if (ppxx > tPmax*1.1):
      break
    pp.append(ppxx)
    vv.append(xx**(-1.5))
    gg.append(gf(xx)+ppxx*xx**(-1.5))
    xx = xx + xd
  try:
    s = interpolate.splrep(pp, gg)
    sv = interpolate.splrep(pp, vv)
  except ValueError:
    print("*******fetal ERROR: BMvol of order: ", BMvol, "  fetal fitting error at T= ", T)
    print(pp)
    print(gg)
    sys.exit()
  gx =interpolate.splev(txx, s)
  vx =interpolate.splev(txx, sv)
  if (Debug==1):
    for i, pp in enumerate(txx):
      print(pp*eVtoGPa, gx[i])
    for i, pp in enumerate(p):
      print(pp*eVtoGPa, g[i])
    sys.exit()
  class result:
    G = gx
    V = vx
  return(result)

def mkDict(line):
  rec = {}
  skiprec = False
  ss = str(line)[0:].replace("'","").split()
  ss = [k for k in ss if k != '']
  for nc,el in enumerate(ss):
    if isfloat(el):
      break
  if len(within)!=0:
    for el in ss[0:nc]:
      if el not in within:
        return True, None, None, None

  _sideal = 0
  _PN = ""
  i = nc*2+2
  while i < len(ss):
    #print ("ncx=", nc, ss[i])
    if ss[i] == "PQ":
      try:
        _PQ = float(ss[i+1])
        #threcord.update({"amount of imaginary phonon mode":float('{:.6f}'.format(_PQ))})
        Uncertainty.update({"amount of imaginary phonon mode":round(_PQ,6)})
        skiprec = _PQ >= PQ
        i += 1
        if skiprec:
          if not paper: print (ss[nc*2+1],"skipped, PQ=", ss[i])
          break
      except:
        skiprec = True
        print ("********Wrong record", ss)
        break
    elif ss[i] == "EQ":
      try:
        _EQ =  float(ss[i+1])
        #threcord.update({"0 K energy uncertainty (eV/atom)":float('{:.6f}'.format(_EQ))})
        Uncertainty.update({"0 K energy uncertainty (eV/atom)":round(_EQ,6)})
        skiprec = _EQ >= EQ
        i += 1
        if skiprec:
          print (ss[nc*2+1],"skipped, EQ=", ss[i])
          break
      except:
        skiprec = True
        print ("********Wrong record", ss)
        break
    elif ss[i] == "PN":
        _PN = ss[i+1].strip("/")
        threcord.update({"Phase name":_PN})
        i += 1
        mpid = ""
        try:
          mp = _PN.index("mp-")
          mpid = _PN[mp:].split('_')[0]
        except:
          pass
        structure.update({"mpid":mpid})
    elif ss[i] == "E0":
        #threcord.update({"static energy (eV/atom)":float('{:.6f}'.format(float(ss[i+1].strip("/"))))})
        threcord.update({"Static energy (eV/atom)":round(float(ss[i+1]),6)})
        i += 1
    elif ss[i] == "TT":
        Tup = float(ss[i+1])
        threcord.update({"Tmax":Tup})
        i += 1
        if Tup < Tupmax:
          print (ss[nc*2+1],"skipped, Tmax=", ss[i])
          skiprec = True
          break
    elif isfloat(ss[i]):
      _sideal = float(ss[i])
    i += 1

  if skiprec:
    return True, None, None, None

  if nc!=0:
    threcord.update({"Uncertainty":Uncertainty})
    space = ss[nc*2].strip("/").split("|")
    structure.update({"space group":int(space[0])})
    structure.update({"point group symmetry":space[1]})
    structure.update({"space group symmetry":space[2]})
    structure.update({"primitive unit cell formula":space[3]})
    elist, clist = formula2composition(space[3])
    pnatom = sum(clist)
    structure.update({"number of atoms in the primitive unit cell":int(pnatom)})
  
    tComponents = ss[0:nc]
    tnComponents = np.array(list(map(int,ss[nc:nc+nc])))
    natom = sum(tnComponents)
    tnComponents = tnComponents/natom
  
    Components = sorted(set(tComponents))
    nComponents = np.zeros(len(Components))
    for i0,el in enumerate(tComponents):
      ix = Components.index(el)
      nComponents[ix] = nComponents[ix] +  tnComponents[i0]
  
    compositions = []
    for i in range(len(Components)):
      compositions.append(int(0.1+natom*nComponents[i]))
    threcord.update({"Elements":Components})
    threcord.update({"Occupancies":list(compositions)})
    
  
    i = nc*2+2
    while i < len(ss):
      if ss[i] == "disordered":
        if i+1>=len(ss):
          _sideal = -sum(nComponents*np.log(nComponents))
        elif isfloat(ss[i+1]):
          i += 1
          if float(ss[i])<0.0:
            _sideal = -sum(nComponents*np.log(nComponents))
          else:
            _sideal = float(ss[i])
        else:
          _sideal = -sum(nComponents*np.log(nComponents))
      i += 1
    threcord.update({"Ideal mixing entropy (kB/atom)":_sideal})

  keys = threcord.keys()
  if nc==0: nc=-1
  vdos_e = str(ss[nc+nc+1]).replace('//','/')
  threcord.update({"Calculation date":str(datetime.datetime.fromtimestamp(os.path.getmtime(vdos_e)))})
  #threcord.update({"Calculation date":str(date.fromtimestamp(os.path.getatime(vdos_e)))})
  if _PN=="":
    try:
      _PN = [s for s in vdos_e.split('/') if s!=""][-3]
    except:
      _PN = "unknown"

  dir0 = _PN
  idx = 1
  while True:
    if not os.path.exists(dir0): break
    recordfile = dir0+"/record.json"
    newdir = False
    try:
      if os.path.exists(recordfile):
        with open(recordfile) as jsonfile:
          orec = json.load (jsonfile)
        okeys = orec.keys()
        for k in keys:
          v = threcord.get(k)
          for ok in okeys:
            if ok != k: continue
            newdir = orec.get(ok) != v
            if newdir: break 
          if k == "Static energy":
            if k in okeys:
              if abs(float(v-okeys.get(k))) < THR0: newdir = False
    except:
      pass

    if not newdir: break
    idx += 1
    dir0 = _PN+"#"+str(idx)

  oldPN = dir0
  newPN = PhaseName.get(_PN)
  if newPN != None: dir0 = newPN
  #print (_PN,dir0,PhaseName)

  threcord.update({"Phase name":dir0})

  """
  if paper:
    pname = threcord.get("primitive unit cell formula")
    n0 = 1
    for px in papers:
      if pname == px.split('#')[0]:
        n0 += 1
    if n0!=1:
      pname = pname +'#'+str(n0)
    papers.append(pname)
    dir0 = pname
  global start
  print ("thermo files extracting cost", time.time()-start)
  start = time.time()
  """
  return False, vdos_e, dir0, oldPN
  

def VASPResults(dir0,vdos_e,Vfiles, Pfiles, phdft="phonon"):
  hdir = vdos_e.replace('thermo/vdos_e','')
  
  natom = structure.get("number of atoms in POSCAR")
  pdir = vdos_e.replace('thermo/vdos_e',phdft)+'/'
  phdir = dir0+'/phonon'
  if not os.path.exists(phdir):
    os.mkdir(phdir)
  for ff in Vfiles:
    vdir = dir0+"/"+ff
    if not os.path.exists(vdir):
      os.mkdir(vdir)
    pvdir = phdir+"/"+ff
    if not os.path.exists(pvdir):
      os.mkdir(pvdir)
       
    vdos = hdir+"/"+ff+"/vdos.out"
    copyfile(vdos,vdir+'/vdos.out')
    print(vdos)

    poscar = hdir+"/"+ff+"/CONTCAR"
    if not os.path.exists(poscar):
      poscar = hdir+"/"+ff+"/Static.CON"
      if not os.path.exists(poscar):
        poscar = hdir+"/"+ff+"/POSCAR"
    copyfile(poscar,vdir+'/POSCAR')

    outcar = hdir+ff+"/OUTCAR"
    if os.path.exists(outcar):
      output = subprocess.run("grep POTCAR "+outcar+" | sort -u", shell=True, stdout=subprocess.PIPE, 
                        universal_newlines=True)
      POTCAR = str(output.stdout)
    else:
      outcar = pdir+"/"+ff+"/OUTCAR.gz"
      if not os.path.exists(outcar):
        outcar = pdir+"/"+ff+"/Static.OUT.gz"
      output = subprocess.run("zgrep POTCAR "+outcar+" | sort -u", shell=True, stdout=subprocess.PIPE, 
                        universal_newlines=True)
      POTCAR = str(output.stdout)

    #print (outcar)
    with open(vdir+'/POTCAR', "w") as text_file:
      text_file.write(POTCAR)

    doscar = hdir+"/"+ff+"/DOSCAR"
    if not os.path.exists(doscar):
      doscar = hdir+"/"+ff+"/Static.DOS.gz"
    ddoscar = doscar.split('/')[-1].replace("Static.DOS", "DOSCAR")
    copyfile(doscar, vdir+'/'+ddoscar)

    oszicar = hdir+"/"+ff+"/OSZICAR"
    if not os.path.exists(doscar):
      oszicar = hdir+"/"+ff+"/Static.OSZ"
    output = subprocess.run("grep E0= "+oszicar+" | tail -1 | awk '{print $5}'", shell=True, stdout=subprocess.PIPE, 
                        universal_newlines=True)
    E0 = float(output.stdout)
    E0 /= natom
    with open(vdir+'/energy', "w") as text_file:
      text_file.write(str(E0))

    incars = fnmatch.filter(os.listdir(hdir+"/"+ff), 'INCAR_*')
    if len(incars)==0:
      INCAR = hdir+"/tplate/INCAR.Static"
    else:
      INCAR = hdir+"/"+ff+'/'+incars[0]
      d0 = os.path.getmtime(INCAR)
      for ii in range(1,len(incars)):
        d1 = os.path.getmtime(hdir+"/"+ff+'/'+incars[ii])
        if d1 > d0:
          INCAR = hdir+"/"+ff+'/'+incars[ii]
          d1 = d0
    copyfile(INCAR, vdir+'/INCAR')
    copyfile(hdir+"/tplate/KPOINTS", vdir+'/KPOINTS')

    sposcar = pdir+ff+"/POSCAR"
    copyfile(sposcar, pvdir+'/POSCAR')
    soutcar = pdir+ff+"/OUTCAR"
    if os.path.exists(soutcar):
      output = subprocess.run("grep POTCAR "+soutcar+" | sort -u", shell=True, stdout=subprocess.PIPE, 
                        universal_newlines=True)
      POTCAR = str(output.stdout)
    else:
      soutcar = pdir+"/"+ff+"/OUTCAR.gz"
      output = subprocess.run("zgrep POTCAR "+soutcar+" | sort -u", shell=True, stdout=subprocess.PIPE, 
                        universal_newlines=True)
      POTCAR = str(output.stdout)

    with open(pvdir+'/POTCAR', "w") as text_file:
      text_file.write(POTCAR)

    incars = fnmatch.filter(os.listdir(pdir+"/"+ff), 'INCAR_*')
    if len(incars)==0:
      INCAR = pdir+"/tplate/INCAR"
    else:
      INCAR = pdir+"/"+ff+'/'+incars[0]
      d0 = os.path.getmtime(INCAR)
      for ii in range(1,len(incars)):
        d1 = os.path.getmtime(pdir+"/"+ff+'/'+incars[ii])
        if d1 > d0:
          INCAR = pdir+"/"+ff+'/'+incars[ii]
          d1 = d0
    copyfile(INCAR, pvdir+'/INCAR')
    copyfile(pdir+"/tplate/KPOINTS", pvdir+'/KPOINTS')

    sposcar = pdir+ff+"/POSCAR"
    sxml = pdir+ff+"/vasprun.xml"
    if not os.path.exists(sxml):
      sxml = pdir+ff+"/vasprun.xml.gz"

    cwd = os.getcwd()
    os.chdir( pvdir )
    cmd = 'vasp_fij -outc '+soutcar+" -xml "+sxml+" -conc "+sposcar + " >& /dev/null"
    os.system(cmd)
    os.chdir( cwd )

  global start
  print ( round(time.time()-start,3), "Secs. costed in VASP files extracting")
  start = time.time()

def extractGph():
  phononmode = {}
  with open("symmetry.out", "r") as f:
    lines = f.readlines()
  i = 0
  while i < len(lines):
    ss = [s for s in lines[i].strip().split(' ') if s!='']
    if len(ss) >= 5:
      if ss[2] == "Modes" and ss[4] in ["silent_mode", "raman_active", "ir_active"]:
        mode = []
        for ii in range(int(ss[0])):
          i += 1
          mm = [s for s in lines[i].strip().replace('(',' ').replace(')',' ').split(' ') if s!='']
          mode.append(float(mm[3]))
        phononmode.update({ss[1]+" ( "+ss[4]+" )": sorted(mode)})
    i += 1
  threcord.update({"gamma point phonons (cm-1) ":phononmode})
    

class BornMix:
    def __init__(self, dir0, V0, V1, ff1, phdir298):
        F0 = dir0+'/'+V0+'/dielecfij.out'
        if not os.path.exists(F0): return
        F1 = dir0+'/'+V1+'/dielecfij.out'
        if not os.path.exists(F1): return
        with open (F0, 'r') as fp: data0 = fp.readlines()
        with open (F1, 'r') as fp: data1 = fp.readlines()
        out = phdir298+'/dielecfij.out'
        with open (out, 'w') as fp:
            for i, line in enumerate(data0):
                self.mix(line, data1[i], ff1, fp)

    def simplemix(self, ss0, ss1, ff1, fp):
        for i,s0 in enumerate(ss0):
            if s0 == ss1[i]: fp.write(' {}'.format(s0))
            else: fp.write(' {}'.format( float(s0)+(1.-float(ff1))*(float(ss1[i])-float(s0)) ))
        fp.write('\n')

    def sitemix(self, ss0, ss1, ff1, fp):
        for i,s0 in enumerate(ss0):
            if s0 == ss1[i]: fp.write(' {}'.format(s0))
            else: 
                change = float(ss1[i])-float(s0)
                if change>=0.5: change -= 1.
                if change<=-0.5: change += 1.
                fp.write(' {}'.format( float(s0)+(1.-float(ff1))*change ))
        fp.write('\n')

    def mix(self, line0, line1, ff1, fp):
        ss0 = [f.strip() for f in line0.split(' ') if f!='']
        ss1 = [f.strip() for f in line1.split(' ') if f!='']
        if len(ss0) >= 4: 
            if ss0[3] in periodictable and ss0[3] in periodictable: self.sitemix(ss0, ss1, ff1, fp)
            else: self.simplemix(ss0, ss1, ff1, fp)
        else:
            self.simplemix(ss0, ss1, ff1, fp)


def Phonon298(dir0, pvdos=False):
  V298 = threcord.get("Atomic volume at 298.15 K ($\AA^3$)")
  phdir298 = dir0 + '/phonon298.15K'
  if not os.path.exists(phdir298):
    os.mkdir(phdir298)
  volumes = threcord.get("volumes")
  i1 = 0
  for ii,vv in enumerate(volumes):
    if float(vv) < V298:
      i1 += 1
  i1 -= 1
  i1 = max(i1, 0)
  i1 = min(i1, len(volumes)-2)
  dV = float(volumes[i1+1]) - float(volumes[i1])
  ff1 = (float(volumes[i1+1]) - V298)/dV
  cmd = "Ymix -mlat -f "+str(ff1)+ " " + dir0+'/'+Pfiles[i1]+"/superfij.out " + " " + dir0+'/'+Pfiles[i1+1]+"/superfij.out >"+phdir298+"/superfij.out"
  output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                      universal_newlines=True)
  mix = BornMix(dir0, Pfiles[i1], Pfiles[i1+1], ff1, phdir298)

  cwd = os.getcwd()
  os.chdir( phdir298 )

  cmd = "Yphon -tranI 2 -eps -nqwave "+ str(nqwave)+ " <superfij.out"
  if os.path.exists('dielecfij.out') : cmd = cmd + ' -Born dielecfij.out'
  print(cmd)
  output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    universal_newlines=True)
  cmd = "gnuplot vdos.plt; convert background white -alpha remove -rotate 90 -density 120x120 vdos.eps vdos.png"
  output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    universal_newlines=True)
  figures = threcord.get("figures")
  figures.update({"phonon DOS at 298.15 K": "phonon298.15K/vdos.png"})

  if pvdos:
    cmd = "Yphon -tranI 2 -eps -pvdos -nqwave "+ str(nqwave/4)+ " <superfij.out"
    if os.path.exists('dielecfij.out') : cmd = cmd + ' -Born dielecfij.out'
    print(cmd)
    output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                      universal_newlines=True)
    cmd = "gnuplot pvdos.plt; convert background white -alpha remove -rotate 90 -density 120x120 pvdos.eps pvdos.png"
    output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                      universal_newlines=True)
    figures = threcord.get("figures")
    figures.update({"generalized phonon DOS at 298.15 K": "phonon298.15K/pvdos.png"})
  os.chdir( cwd )

  #ngroup = threcord.get("structure").get("space group")
  ngroup = structure.get("space group")
  dfile = ""
  if ngroup>=1 and ngroup<=2:
    dfile = home+"/bin/hbin/pycode/data/dfile.tri"
    structure.update({"crystal system": "Triclinic"})
  elif ngroup>=3 and ngroup<=15:
    dfile = home+"/bin/hbin/pycode/data/dfile.mon"
    structure.update({"crystal system": "Monoclinic"})
  elif ngroup>=16 and ngroup<=74:
    dfile = home+"/bin/hbin/pycode/data/dfile.oth"
    structure.update({"crystal system": "Orthorhombic"})
  elif ngroup>=75 and ngroup<=142:
    dfile = home+"/bin/hbin/pycode/data/dfile.tet"
    structure.update({"crystal system": "Tetragonal"})
  elif ngroup>=143 and ngroup<=167:
    dfile = home+"/bin/hbin/pycode/data/dfile.rho"
    structure.update({"crystal system": "Trigonal"})
  elif ngroup>=168 and ngroup<=194:
    dfile = home+"/bin/hbin/pycode/data/dfile.hcp"
    structure.update({"crystal system": "Hexagonal"})
  elif ngroup>=195 and ngroup<=220:
    dfile = home+"/bin/hbin/pycode/data/dfile.scc"
    structure.update({"crystal system": "Cubic"})
  elif ngroup>=221 and ngroup<=224:
    dfile = home+"/bin/hbin/pycode/data/dfile.bcc"
    structure.update({"crystal system": "Cubic({bcc})"})
  elif ngroup>=225 and ngroup<=230:
    dfile = home+"/bin/hbin/pycode/data/dfile.fcc"
    structure.update({"crystal system": "Cubic({fcc})"})

  if dfile != "":
    dfile0 = dfile.split('/')[-1]
    copyfile(dfile,phdir298+'/'+dfile0)
    cwd = os.getcwd()
    os.chdir( phdir298 )
    cmd = 'timeout 6 pos2s Symmetry.pos -THR 3.e-4 >&symmetry.out'
    output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    universal_newlines=True)

    Gph = os.path.exists("symmetry.mode")
    if Gph:
      cmd = "Yphon -Gfile symmetry.mode -tranI 2 -eps -pdis "+dfile0+ " <superfij.out >symmetry.out"
    else:
      cmd = "Yphon -tranI 2 -eps -pdis "+dfile0+ " <superfij.out >symmetry.out"
    if os.path.exists('dielecfij.out') : cmd = cmd + ' -Born dielecfij.out'

    output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                      universal_newlines=True)
    if Gph:
      extractGph()

    cmd = "gnuplot vdis.plt; convert -background white -alpha remove -rotate 90 -density 120x120 vdis.eps vdis.png"
    output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                      universal_newlines=True)
    figures.update({"phonon dispersion at 298.15 K": "phonon298.15K/vdis.png"})
    os.chdir( cwd )

  threcord.update({"figures":figures})
  
  cwd = os.getcwd()
  os.chdir( phdir298 )
  cfile = ["findsym.log","run.log","exactQ.out","run","symmetry.out"]
  for f in cfile:
    if os.path.exists(f):
      os.remove(f)
  os.chdir( cwd )
    
  global start
  print (round(time.time()-start,3), "Secs. costed in calculations of phonon properties at 298.15K")
  start = time.time()

#outs = ["space group", "point group symmetry", "space group symmetry"]
outs = ["space group","space group symmetry"]
def addpapers(g,formula,pname):
  g[0] += float(threcord.get("H298.15 (J/mol-atom)"))
  if threcord.get("Ideal mixing entropy (kB/atom)")!=None:
    g[1] -= R*float(threcord.get("Ideal mixing entropy (kB/atom)"))

  sys.stdout.write("{},{}".format(pname,formula))
  for ss in outs:
    sys.stdout.write(",{}".format(structure.get(ss)))
  sys.stdout.write(",{:.6g}".format(g[0]))
  sys.stdout.write(",{:.6g}".format(g[1]))
  sys.stdout.write(",{:.5g}".format(g[2]))
  sys.stdout.write(",{:.3e}".format(g[3]))
  sys.stdout.write(",{:.3e}".format(g[4]))
  sys.stdout.write(",{:.3e}".format(g[5]))

  try:
      pq = 100.*threcord.get("Uncertainty").get("amount of imaginary phonon mode")
      eq = eVtoJ*threcord.get("Uncertainty").get("0 K energy uncertainty (eV/atom)")
  except:
      pq = 0.0
      eq = 0.0
  warning = SGTErec.get("G-H298.15 (J/mol-atom)")[1].get("error")
  #print(pq,eq,warning)
  #sys.stdout.write(",{:.1f},{:.0f},{:.0f},".format(pq,eq,warning))
  sys.stdout.write(",{:.1f},".format(pq))

  if warning > 99:
    sys.stdout.write("   **********WARNING, fittit error is too large! {}".format(warning)) 
  if abs(g[2]) > 32:
    sys.stdout.write("   **********WARNING, Cp at room condition is abnormal! {}".format(abs(g[2]))) 
  sys.stdout.write("\n")


markers=['o', 'v', 'd', '^', '<', '>', 's', '*', 'x', '+', '1', '2']

home = "/global/u2/y/yiwang62/"
k_B = 8.6173303e-5
R = 8.3144598

eVtoGPa = 160.21766208 
eVtoJ = 96486.9
THRE0 = 1.e-5
nqwave = 2.e6
nqwave = 1.e6
nqwave = 4.e6

T0 = 298.15
update = True

input_within = False #""" key to cotrol the within input"""
formula_within = "" #"""chemical formula"""
within = []

PQ = 0.075
EQ = 0.015
PQ = 0.01
EQ = 0.01
CpMax = 50.
Tupmax = 2000.0
start = time.time()
gamma_phonons = {}
threcord = {}
figures = {}
zthermo = {}
SGTErec = {}
structure = {}
expt = None
Uncertainty = {}
nphases = 0
debug = False
fitCp = False
paper = True
phdft = "phonon"
expt = None
xlim = None
pvdos = False

phases = []
papers = []
PhaseName = {}

def getdoslim(e, dos, xlim):
    xx, yy = [], []
    for i,energy in enumerate(e):
        if energy >xlim[0] and energy <xlim[1]:
            xx.append(energy)
            yy.append(dos[i])
    return xx, yy

def plotAPI(readme, thermofile, volumes=None, energies=None, expt=None, xlim=None, _fitCp=True,
    formula=None, debug=False, vtof=None, poscar=None, vdos=None, doscar=None, natoms=1, plotlabel=None):
  if plotlabel!=None:
      if plotlabel.lower().startswith("find_or_"):
          if "pseudo_potential" in readme.keys():
              plotlabel = readme["pseudo_potential"]['functional']
          else:
              plotlabel = plotlabel.replace("find_or_", "")
  else:
      plotlabel = 'DFT'
  if expt!=None: expt =get_expt(expt, formula)


  global fitCp
  fitCp = _fitCp
  phasedir = [substr for substr in thermofile.split('/') if substr!=""]
  plottitle = phasedir[-2]
  phasedir = ('/').join(phasedir[0:-1])
  if phasedir=="": phasedir="."
  folder = phasedir+"/figures/"
  print("All figures will be outputed into: ", folder, "  with T uplimt:", xlim, "\n\nEnjoy!\n")
  if not os.path.exists(folder):
    os.mkdir(folder)
  if volumes is not None: thermoplot(folder,"0 K total energies (eV/atom)",volumes, energies, plottitle=plottitle)

  thermo = np.loadtxt(thermofile, comments="#", dtype=np.float)
  thermo[np.isnan(thermo)] = 0.0
  _single = len(set(thermo[:,1])) == 1
  if len (thermo) < 1:
      print("\nCorrupted thermofile for", thermofile, "Please check it!")
      return False

  if vtof is not None:
    thermoplot(folder,"Helmholtz energy (eV/atom)",list(thermo[:,1]),list(thermo[:,2]), reflin=vtof,plottitle=plottitle)
    thermoplot(folder,"Helmholtz energy analysis (eV/atom)",list(thermo[:,1]),list(thermo[:,2]), reflin=vtof,plottitle=plottitle)
    return
  
           
  if xlim is not None:
      for i,x in enumerate(thermo[:,0]):
          if x>=xlim and x!=thermo[-1,0]:
              thermo = thermo[0:i+1,:]
              xlim = None
              break

  if expt!=None: 
      meltingT = get_melting_temperature(expt, formula)
      if meltingT!=None:
          for i,x in enumerate(thermo[:,0]):
              if x>=meltingT and x!=thermo[-1,0]:
                  thermo = thermo[0:i+1,:]
                  break

  for i,cp in enumerate(thermo[:,6]):
    if cp > CpMax: 
      thermo = thermo[0:i,:]
      break
    
  """
  f2=interpolate.splrep(thermo[:,0], thermo[:,1])
  V298 = float(interpolate.splev(T0, Vstack))
  Hstack=interpolate.splrep(thermo[:,0], thermo[:,4])
  H298 = float(interpolate.splev(T0, Hstack))
  Sstack=interpolate.splrep(thermo[:,0], thermo[:,3])
  S298 = float(interpolate.splev(T0, Sstack))
  """
  if T0 <= thermo[-1,0] : 
    f2=interp1d(thermo[:,0], thermo[:,1])
    V298 = f2(T0)
    f2=interp1d(thermo[:,0], thermo[:,4])
    H298 = f2(T0)
    f2=interp1d(thermo[:,0], thermo[:,3])
    S298 = f2(T0)
  else:
    print ("\nWarning! T0=", T0, "is higher than the T up limit:", thermo[-1,0], \
    " no SGTE fitting will be performed\n")

  if volumes is not None: 
    if T0 <= thermo[-1,0] :
      T = thermo[:,0]
      V = thermo[:,1]
      A = thermo[:,5]
      B = thermo[:,9]
      C = thermo[:,7]
      G = 3*A[T>T0]*B[T>T0]*physical_constants['Avogadro constant'][0]*1e-21*V[T>T0]/C[T>T0]
      v = V[T>T0]
      g = V298*G/v
      g = sum(g)/len(g)
      readme['Gruneisen parameter']= round(g,3)
      Gmax = 1.2*max(G)
      g = 3*A[T>0]*B[T>0]*physical_constants['Avogadro constant'][0]*1e-21*V[T>0]/C[T>0]
      t = T[T>0]
      Gmin = min(g)
      if Gmax>0: Gmin = max(Gmin,-Gmax)
      if Gmin>0: Gmin = 0
      ix = 0
      for i in range (len(g)-2,0,-1):
          if g[i]>Gmax or g[i] <Gmin:
              ix = i+1
              break
      g = g[ix:]
      t = t[ix:]
      thermoplot(folder,"Gruneisen coefficient",list(t),list(g), yzero=Gmin, expt=expt, xlim=xlim, label=plotlabel, single=vdos!=None,plottitle=plottitle)
      Plot298(folder, V298, volumes, debug=debug, plottitle=plottitle)
    else:
      print ("\nWarning! T0=", T0, "is higher than the T up limit:", thermo[-1,0], \
      " phonon perperties will be reported at 0 K\n")
      f2=interp1d(thermo[:,0], thermo[:,1])
      V0 = f2(0)
      Plot298(folder, V0, volumes, debug=debug)
  elif vdos!=None:
      #check if superfij.out file exist there. if yes, do phonon properties
      PlotVol(folder, vdos)

  if T0 <= thermo[-1,0] : 
    threcord.update({"H298.15 (J/mol-atom)":H298})
    threcord.update({"S298.15 (J/mol-atom/K)":S298})

  zthermo.update({"temperature (K)":list(thermo[:,0])})
  zthermo.update({"atomic volume ($\AA^3$)":list(thermo[:,1])})
  zthermo.update({"Gibbs energy (eV/atom)":list(thermo[:,2])})
  zthermo.update({"enthalpy (J/mol-atom)":list(thermo[:,4])})
  zthermo.update({"entropy (J/mol-atom/K)":list(thermo[:,3])})
  zthermo.update({"Cp (J/mol-atom/K)":list(thermo[:,6])})
  if T0 < thermo[-1,0] : 
    if fitCp:
      proStoichiometricCp()
    else:
      proStoichiometricG()
    with open(folder + '/../record.json', 'w') as fp:
      myjsonout(SGTErec, fp, indent="", comma="")
    myjsonout(SGTErec, sys.stdout, indent="", comma="")

  if volumes is not None:
    thermoplot(folder,"Atomic volume ($\AA^3$)",list(thermo[:,0]),list(thermo[:,1]), xlim=xlim, label=plotlabel,plottitle=plottitle)
  if T0 <= thermo[-1,0] : 
    thermoplot(folder,"Gibbs energy-H298 (J/mol-atom)",list(thermo[:,0]),list(thermo[:,2]*eVtoJ-H298), xlim=xlim,plottitle=plottitle)
    thermoplot(folder,"Enthalpy-H298 (J/mol-atom)",list(thermo[:,0]),list(thermo[:,4]-H298), 
      expt=expt, xlim=xlim,plottitle=plottitle)
  thermoplot(folder,"Entropy (J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,3]),yzero=0.0, expt=expt, 
      xlim=xlim,plottitle=plottitle)

  if volumes is not None:
    thermoplot(folder,"LTC (1/K)",list(thermo[:,0]),list(thermo[:,5]),yzero=0.0, expt=expt, xlim=xlim, label=plotlabel,plottitle=plottitle)
    #thermoplot(folder,"LTC analysis (1/K)",list(thermo[:,0]),list(thermo[:,5]),reflin=list(thermo[:,22]), yzero=0.0, xlim=xlim, label=plotlabel,plottitle=plottitle)
  ncols = [6,8]
  thermoplot(folder,"Heat capacities (J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), expt=expt, xlim=xlim, label=plotlabel, single=_single,plottitle=plottitle)
  thermoplot(folder,"Heat capacities (J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), xlim=300,expt=expt, label=plotlabel, single=_single,plottitle=plottitle)
  thermoplot(folder,"Heat capacities (J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), xlim=100,expt=expt, CoT=True, label=plotlabel, single=_single,plottitle=plottitle)
  thermoplot(folder,"Heat capacities (J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), xlim=1000,expt=expt, CoT=True, label=plotlabel, single=_single,plottitle=plottitle)
  tmp = 0.0
  for i,v in enumerate(thermo[:,0]):
    if v >300: break
    tmp = max(tmp, thermo[i,6]-thermo[i,8])
  if tmp>1.e-2:
    thermoplot(folder,"Heat capacities (J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), elonly=300, expt=expt, CoT=True, label=plotlabel,plottitle=plottitle)
  thermoplot(folder,"Debye temperature (K)",list(thermo[:,0]),list(thermo[:,10]),yzero=0.0, xlim=xlim, label=plotlabel,plottitle=plottitle)
  thermoplot(folder,"Debye temperature (K)",list(thermo[:,0]),list(thermo[:,10]),yzero=0.0, xlim=70, label=plotlabel,plottitle=plottitle)
  if volumes is not None:
    bs = np.ones((len(thermo[:,9])), dtype=float)
    bs[1:] = thermo[1:,6]/thermo[1:,7]*thermo[1:,9]
    thermoplot(folder,"Bulk modulus (GPa)",list(thermo[:,0]),list(thermo[:,9]), reflin=list(bs) , expt=expt, yzero=0.0,xlim=xlim, label=plotlabel,plottitle=plottitle)
  thermoplot(folder,"Seebeck coefficients (μV/K)",list(thermo[:,0]),list(thermo[:,16]),xlim=xlim, label=plotlabel,plottitle=plottitle)
  thermoplot(folder,"Lorenz number ($WΩK^{−2}$)",list(thermo[:,0]),list(thermo[:,17]),xlim=xlim, label=plotlabel,plottitle=plottitle)
  thermoplot(folder,"Absolute thermal electric force (V)",list(thermo[:,0]),list(thermo[:,15]), xlim=xlim, label=plotlabel,plottitle=plottitle)
  thermoplot(folder,"Effective charge carrier concentration ($e/cm^{3}$)",list(thermo[:,0]),
      list(thermo[:,18]/thermo[:,1]*1e24), label=plotlabel,plottitle=plottitle)
  thermoplot(folder,"Effective charge carrier concentration ($e/cm^{3}$)",list(thermo[:,0]),
      list(thermo[:,18]/thermo[:,1]*1e24), xlim=100, label=plotlabel,plottitle=plottitle)
  if len(gamma_phonons)!=0: readme['gamma phonons (cm^{-1})']= gamma_phonons
  if doscar!=None:
      from dfttk.pythelec import pregetdos, getdos
      with open (doscar, "r") as fp:
          edn, eup, vde, dos_energies, vaspEdos = pregetdos(fp) # Line 186
          NELECTRONS, E0, dF, e, dos, Eg =\
              getdos(-15, 15, 0.0, 10001, 1000., edn, eup, vde, dos_energies, vaspEdos)
          if Eg <0.0: Eg=0.
          xlim = [-0.1, Eg+0.1]
          xx, yy = getdoslim(dos_energies, vaspEdos, xlim)
          thermoplot(folder,"Electron DOS (States/Atom/eV)",list(xx),list(np.array(yy)/natoms), xlim=xlim,
              xlabel="Band energy (eV)", label=plotlabel,plottitle=plottitle)
          xlim = [-0.2, Eg+0.2]
          xx, yy = getdoslim(dos_energies, vaspEdos, xlim)
          thermoplot(folder,"Electron DOS (States/Atom/eV)",list(xx),list(np.array(yy)/natoms), xlim=xlim,
              xlabel="Band energy (eV)", label=plotlabel,plottitle=plottitle)
          xlim = [-0.5, Eg+0.5]
          xx, yy = getdoslim(dos_energies, vaspEdos, xlim)
          thermoplot(folder,"Electron DOS (States/Atom/eV)",list(xx),list(np.array(yy)/natoms), xlim=xlim,
              xlabel="Band energy (eV)", label=plotlabel,plottitle=plottitle)
          xlim = [-1.0, Eg+1.0]
          xx, yy = getdoslim(dos_energies, vaspEdos, xlim)
          thermoplot(folder,"Electron DOS (States/Atom/eV)",list(xx),list(np.array(yy)/natoms), xlim=xlim,
              xlabel="Band energy (eV)", label=plotlabel,plottitle=plottitle)
          xlim = [-2.0, Eg+2.0]
          xx, yy = getdoslim(dos_energies, vaspEdos, xlim)
          thermoplot(folder,"Electron DOS (States/Atom/eV)",list(xx),list(np.array(yy)/natoms), xlim=xlim,
              xlabel="Band energy (eV)", label=plotlabel,plottitle=plottitle)
          xlim = [-5.0, Eg+5.0]
          xx, yy = getdoslim(dos_energies, vaspEdos, xlim)
          thermoplot(folder,"Electron DOS (States/Atom/eV)",list(xx),list(np.array(yy)/natoms), xlim=xlim,
              xlabel="Band energy (eV)", label=plotlabel,plottitle=plottitle)
          xlim = [-10., Eg+10.]
          xx, yy = getdoslim(dos_energies, vaspEdos, xlim)
          thermoplot(folder,"Electron DOS (States/Atom/eV)",list(xx),list(np.array(yy)/natoms), xlim=xlim,
              xlabel="Band energy (eV)", label=plotlabel,plottitle=plottitle)
  return True


def plotCMD(thermofile, volumes=None, energies=None, expt=None, xlim=None, _fitCp=True, 
    poscar=None, vdos=None, doscar=None, natoms=1, plotlabel=None):
  global fitCp
  fitCp = _fitCp
  #print(expt)
  phasedir = [substr for substr in thermofile.split('/') if substr!=""]
  phasedir = ('/').join(phasedir[0:-1])
  if phasedir=="": phasedir="."
  folder = phasedir+"/figures/"
  print("All figures have been outputed into: ", folder, "  with T uplimt:", xlim, "\n\nEnjoy!\n")
  if not os.path.exists(folder):
    os.mkdir(folder)
  if volumes is not None: thermoplot(folder,"0 K total energies (eV/atom)",volumes, energies)

  thermo = np.loadtxt(thermofile, comments="#", dtype=np.float)
  thermo[np.isnan(thermo)] = 0.0
  if len (thermo) < 1:
      print("\nCorrupted thermofile for", thermofile, "Please check it!")
      return False

  if xlim is not None:
      for i,x in enumerate(thermo[:,0]):
          if x>=xlim:
              thermo = thermo[0:i,:]
              xlim = None
              break
           
  for i,cp in enumerate(thermo[:,6]):
    if cp > CpMax: 
      thermo = thermo[0:i,:]
      break
    
  """
  f2=interpolate.splrep(thermo[:,0], thermo[:,1])
  V298 = float(interpolate.splev(T0, Vstack))
  Hstack=interpolate.splrep(thermo[:,0], thermo[:,4])
  H298 = float(interpolate.splev(T0, Hstack))
  Sstack=interpolate.splrep(thermo[:,0], thermo[:,3])
  S298 = float(interpolate.splev(T0, Sstack))
  """
  f2=interp1d(thermo[:,0], thermo[:,1])
  V298 = f2(T0)
  f2=interp1d(thermo[:,0], thermo[:,4])
  H298 = f2(T0)
  f2=interp1d(thermo[:,0], thermo[:,3])
  S298 = f2(T0)

  #print(H298,V298,S298)

  if volumes is not None: Plot298(folder, V298, volumes)

  threcord.update({"H298.15 (J/mol-atom)":H298})
  threcord.update({"S298.15 (J/mol-atom/K)":S298})

  zthermo.update({"temperature (K)":list(thermo[:,0])})
  zthermo.update({"atomic volume ($\AA^3$)":list(thermo[:,1])})
  zthermo.update({"Gibbs energy (eV/atom)":list(thermo[:,2])})
  zthermo.update({"enthalpy (J/mol-atom)":list(thermo[:,4])})
  zthermo.update({"entropy (J/mol-atom/K)":list(thermo[:,3])})
  zthermo.update({"Cp (J/mol-atom/K)":list(thermo[:,6])})
  if fitCp:
    proStoichiometricCp()
  else:
    proStoichiometricG()
  with open(folder + '/../record.json', 'w') as fp:
    myjsonout(SGTErec, fp, indent="", comma="")
  myjsonout(SGTErec, sys.stdout, indent="", comma="")

  thermoplot(folder,"Atomic volume ($\AA^3$)",list(thermo[:,0]),list(thermo[:,1]), xlim=xlim, label=plotlabel)
  thermoplot(folder,"Gibbs energy-H298 (J/mol-atom)",list(thermo[:,0]),list(thermo[:,2]*eVtoJ-H298), xlim=xlim, label=plotlabel)
  #print(thermo[:,4]-H298)
  thermoplot(folder,"Enthalpy-H298 (J/mol-atom)",list(thermo[:,0]),list(thermo[:,4]-H298), 
    expt=expt, xlim=xlim, label=plotlabel)
  thermoplot(folder,"Entropy (J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,3]),yzero=0.0, xlim=xlim, label=plotlabel)

  thermoplot(folder,"LTC (1/K)",list(thermo[:,0]),list(1.e06*thermo[:,5]),yzero=0.0, xlim=xlim, label=plotlabel)
  ncols = [6,8]
  #print('eeeeeeee', plotlabel, expt)
  thermoplot(folder,"Heat capacities (J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), expt=expt, xlim=xlim, label=plotlabel)
  thermoplot(folder,"Heat capacities (J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), xlim=300,expt=expt, label=plotlabel)
  thermoplot(folder,"Heat capacities (J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), xlim=100,expt=expt, CoT=True, label=plotlabel)
  tmp = 0.0
  for i,v in enumerate(thermo[:,0]):
    if v >300: break
    tmp = max(tmp, thermo[i,6]-thermo[i,8])
  if tmp>1.e-2:
    thermoplot(folder,"Heat capacities (J/mol-atom/K)",list(thermo[:,0]),list(thermo[:,ncols]), elonly=300, expt=expt, CoT=True, label=plotlabel)
  thermoplot(folder,"Debye temperature (K)",list(thermo[:,0]),list(thermo[:,13]),yzero=0.0, xlim=xlim, label=plotlabel)
  thermoplot(folder,"Debye temperature (K)",list(thermo[:,0]),list(thermo[:,13]),yzero=0.0, xlim=70, label=plotlabel)
  #thermoplot(folder,"Bulk modulus (GPa)",list(thermo[:,0]),list(thermo[:,15]),yzero=0.0,xlim=xlim, label=plotlabel)
  bs = copy.deepcopy(thermo[:,9])
  for i,Cv in enumerate(thermo[1:,7]):
     if Cv>0.0: bs[i] = thermo[i,6]/Cv*thermo[i,9]
  thermoplot(folder,"Bulk modulus (GPa)",list(thermo[:,0]),list(thermo[:,9]), reflin=list(bs) , expt=expt, yzero=0.0,xlim=xlim, label=plotlabel)
  T = copy.deepcopy(thermo[:,0])
  t22 = copy.deepcopy(thermo[:,22])
  for i,tval in enumerate(t22):
      if T[i] <=0.0 : T[i]=1.e-8
      if t22[i] <=0.0 : t22[i]=1.e-8
  Lfactor = physical_constants['Boltzmann constant'][0]/physical_constants['atomic unit of charge'][0]**2/physical_constants['Avogadro constant'][0]
  thermoplot(folder,"Seebeck coefficients (μV/K)",list(thermo[:,0]),list(thermo[:,21]/t22/T),xlim=xlim, label=plotlabel)
  thermoplot(folder,"Lorenz number ($WΩK^{−2}$)",list(thermo[:,0]),list((thermo[:,6]-thermo[:,8])/t22*Lfactor),xlim=xlim, label=plotlabel)
  thermoplot(folder,"Absolute thermal electric force (V)",list(thermo[:,0]),list(thermo[:,19]), xlim=xlim, label=plotlabel)
  thermoplot(folder,"Effective charge carrier concentration ($e/cm^{3}$)",list(thermo[:,0]),
      list(thermo[:,22]/thermo[:,1]*1e24), label=plotlabel)
  thermoplot(folder,"Effective charge carrier concentration ($e/cm^{3}$)",list(thermo[:,0]),
      list(thermo[:,22]/thermo[:,1]*1e24), xlim=100, label=plotlabel)


def addvdos(x,y,f,w,h):
  for i,v in enumerate(x):
    dx = v - f
    y[i] += math.exp(-(dx/w)**2)*h


def plotRaman(folder, fp, vdos, plottitle=None):
  lines=fp.readlines()
  for i,line in enumerate(lines):
    if line.startswith("Setting workspace & pre-optimizing : Section time "):
      lines = lines[i+2:]
      break
  for i,line in enumerate(lines):
    ff = [f for f in line.strip().split(" ") if f!=""]
    if len(ff) < 3: continue
    if ff[2]=="Modes":
      lines = lines[i:]
      break
  
  I = []
  F_lo = []
  M = []
  F = []
  A = []
  G = {}
  global gamma_phonons
  gamma_phonons = {}
  for i,line in enumerate(lines):
    #if line.startswith("Handling symmetry : Section time "): break
    #if line.startswith("Handling symmetry : Section time "): break
    ff = [f for f in line.strip().split(" ") if f!=""]
    if len(ff) < 3: continue
    if ff[2]=="Modes":
      active = ff[4]
      continue

    if line.startswith(" No irrep        THz"): continue
    try:
      int(ff[0])
    except:
      break

    M.append(ff[1])
    I.append(ff[0])
    F.append(float(ff[2])/0.0299792458)
    A.append(active)
    kk = '{} {:05}'.format(ff[1],int(ff[0]))
    #THz = round(float(ff[2]),3)
    #THz = str(round(float(ff[2]),3))+" THz",
    #cm = str(round(float(ff[2])/0.0299792458,1))+" cm-1"
    cm = round(float(ff[2])/0.0299792458,1)
    try:
      F_lo.append(float(ff[3])/0.0299792458)
      if cm!=0 and active=="ir_active":
        cm = str(round(float(ff[2])/0.0299792458,1))+"(TO)+"+str(round(float(ff[3])/0.0299792458,1))+"(LO)"
    except:
      pass
    if cm!=0: gamma_phonons[kk] = [cm, active]
  #print(I,M,F,A)  
  x = vdos[:,0]*1.e-12/0.0299792458
  y = vdos[:,1]*1.e+12*0.0299792458
  yy = np.zeros((len(y)), dtype=float)
  w = max(x)*0.001
  h = 0.1*max(y)
  x0 = []
  y0 = []
  s0 = []
  for i,f in enumerate(F):
    if M[i].lower().startswith("e"): hh = h*2
    elif M[i].lower().startswith("t"): hh = h*3
    else: hh = h
    if float(f)<1.e-3: continue
    if len(F_lo)!=0:
      if A[i]=="ir_active":
        if M[i].lower().startswith("e"): hh = h
        elif M[i].lower().startswith("t"): hh = h*2
        else: hh = h/2
        addvdos(x,yy,float(f),w,hh)
        x0.append(float(f)-18*w)
        y0.append(hh)
        s0.append(M[i]+"(TO)")
        if hh==h/2: h = hh
        addvdos(x,yy,float(F_lo[i]),w,h)
        x0.append(float(F_lo[i])-18*w)
        y0.append(h)
        s0.append(M[i]+"(LO)")
      else:
        addvdos(x,yy,float(f),w,hh)
        x0.append(float(f)-18*w)
        y0.append(hh)
        s0.append(M[i])
    else:
      addvdos(x,yy,float(f),w,hh)
      x0.append(float(f)-18*w)
      y0.append(hh)
      s0.append(M[i])
  ix = sorted(range(len(x0)), key=lambda k: x0[k])
  _M = []
  _x0 = []
  _y0 = []
  _s0 = []
  for i in range(len(ix)):
    _M.append(s0[ix[i]])
    _x0.append(x0[ix[i]])
    _y0.append(y0[ix[i]])
    ss = s0[ix[i]]
    if len(ss)==1:
      ss = '$'+ss+'$'
    elif len(ss)>0:
      if ss[1].isdigit() or ss[1].isalpha():
        aa = ss[1:len(ss)].split('(')
        if len(aa)>1:
          ss = '$'+ss[0]+'_{'+aa[0]+'}^{('+aa[1]+'}$'
        else:
          ss = '$'+ss[0]+'_{'+aa[0]+'}$'
      else: ss = '$'+ss+'$'
    _s0.append(ss)
  M = _M
  x0 = _x0
  y0 = _y0
  s0 = _s0
  
  #print(x0)
  #print(M)
  nx0 = []
  ny0 = []
  ns0 = []
  nn0 = []
  #adjust y0 for overlapping text
  for i,v in enumerate(x0):
    ff = False
    for j in range(len(nx0)):
      if abs(x0[i] - nx0[j]) < 30*w:
        ns0[j] = ns0[j]+'+'+s0[i]
        nx0[j]= (nx0[j]*nn0[j]+x0[i])/(nn0[j]+1)
        ny0[j]= max(ny0[j], y0[i])
        nn0[j] += 1
        ff = True
        break
    if ff: continue
    nx0.append(x0[i])
    ny0.append(y0[i])
    ns0.append(s0[i])
    nn0.append(1)

  if len(s0)>0:
    for i,s in enumerate(ns0):
      ss = [f for f in s.split('+') if f!=""]
      if len(ss)>7:
        ns0[i] = '+'.join(ss[0:7])+'+...'
    thermoplot("./","Gamma point phonons",list(x),list(yy), 
      reflin=list(y), xlabel="Phonon frequency($cm^{-1}$)", ytext=[nx0,ny0,ns0], ylabel="Phonon DOS ($states.cm$)",plottitle=plottitle)
      #reflin=list(y), xlabel="Phonon frequency(THz)", ytext=[nx0,ny0,ns0], ylabel="Phonon DOS ($THz^{-1}$)")
    fn = "Gamma_point_phonons.png"
    os.rename(fn, folder+'/'+fn)

def Plot298(folder, V298, volumes, debug=False, plottitle=None):
  import dfttk.scripts.config_dfttk as dfttkconfig
  PATH_TO_STORE_CONFIG = dfttkconfig.default_path()
  plotdatabase = dfttkconfig.get_abspath(PATH_TO_STORE_CONFIG)+'/analysis/database/'
  #print (plotdatabase, folder)
  ydir = folder+'/../Yphon/'
  
  structure = None
  for root, dirs, files in os.walk(ydir):
    for dir in dirs:
      poscar = ydir+dir+'/POSCAR'
      if os.path.exists(poscar):
        structure = Structure.from_file(poscar)
        break
    if structure is not None: break
    #structure = Structure.from_file(ydir+dirs[len(dirs)//2]+'/POSCAR')
    #break

  try:
    natom = len(structure.sites)
    sa = SpacegroupAnalyzer(structure)
    ngroup = sa.get_space_group_number()
  except:
    return
      
  #print(natom,ngroup)

  i1 = 0
  for ii,vv in enumerate(volumes):
    if float(vv) < V298:
      i1 += 1
  i1 -= 1
  i1 = max(i1, 0)
  i1 = min(i1, len(volumes)-2)
  dV = float(volumes[i1+1]) - float(volumes[i1])
  ff1 = (float(volumes[i1+1]) - V298)/dV

  file1 = ydir+'/V{:010.6f}/superfij.out'.format(float(natom*volumes[i1]))
  if not os.path.exists(file1):
    print ("\nWARNING! I cannot find file :", file1, " so that I will not do phonon298.15 for you!\n")
    return
  file2 = ydir+'/V{:010.6f}/superfij.out'.format(float(natom*volumes[i1+1]))
  if not os.path.exists(file2):
    print ("\nWARNING! I cannot find file :", file2, " so that I will not do phonon298.15 for you!\n")
    return

  phdir298 = ydir+'/Phonon298.15'
  if not os.path.exists(phdir298):
      os.mkdir(phdir298)
  cmd = "Ymix -mlat -f "+str(ff1)+ " "+file1+ " "+file2 +" >"+phdir298+"/superfij.out"
  print(cmd)
  output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                      universal_newlines=True)

  mix = BornMix(ydir, 'V{:010.6f}'.format(float(natom*volumes[i1])), \
                      'V{:010.6f}'.format(float(natom*volumes[i1+1])), ff1, phdir298)

  cwd = os.getcwd()
  os.chdir( phdir298 )

  if debug: cmd = "Yphon -tranI 2 -eps <superfij.out"
  else: cmd = "Yphon -tranI 2 -eps -nqwave "+ str(nqwave)+ " <superfij.out"
  if os.path.exists('dielecfij.out') : cmd = cmd + ' -Born dielecfij.out'
  #cmd = "Yphon -tranI 2 -eps " + " <superfij.out"
  if not (debug and os.path.exists('vdos.out')):
      print(cmd)
      output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    universal_newlines=True)
      cmd = "gnuplot vdos.plt; convert -flatten -rotate 90 -density 120x120 vdos.eps vdos.png"
      #print(cmd)
      output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    universal_newlines=True)
  
      #copyfile("vdos.png", folder+'/vdos298.15.png')
      os.rename("vdos.eps", cwd+'/'+folder+'/vdos298.15.eps')
      os.rename("vdos.png", cwd+'/'+folder+'/vdos298.15.png')

  if not os.path.exists('symmetry.mode'):
      cmd = "pos2s Symmetry.pos -THR 0.001"
      print(cmd)
      output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    universal_newlines=True)
  """
  #temp for debug
  """
  #temp for debug

  if os.path.exists("vdos.out") :
    cmd = "Yphon -tranI 2 -eps -nqwave 100 -Gfile symmetry.mode <superfij.out >Raman.mode"
    if os.path.exists('dielecfij.out') : cmd = cmd + ' -Born dielecfij.out'
    os.rename("vdos.out", 'vdos.sav')
    output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    universal_newlines=True)
    os.rename("vdos.sav", 'vdos.out')
    vdos = np.loadtxt("vdos.out", comments="#", dtype=np.float)
    if os.path.exists("Raman.mode") :
      with open ("Raman.mode", "r") as fp:
        plotRaman(cwd+'/'+folder, fp, vdos, plottitle=plottitle)      

  dfile = ""
  if ngroup>=1 and ngroup<=2:
    dfile = plotdatabase+"/dfile.tri"
  elif ngroup>=3 and ngroup<=15:
    dfile = plotdatabase+"/dfile.mon"
  elif ngroup>=16 and ngroup<=74:
    dfile = plotdatabase+"/dfile.oth"
  elif ngroup>=75 and ngroup<=142:
    dfile = plotdatabase+"/dfile.tet"
  elif ngroup>=143 and ngroup<=167:
    dfile = plotdatabase+"/dfile.rho"
  elif ngroup>=168 and ngroup<=194:
    dfile = plotdatabase+"/dfile.hcp"
  elif ngroup>=195 and ngroup<=220:
    dfile = plotdatabase+"/dfile.scc"
  elif ngroup>=221 and ngroup<=224:
    dfile = plotdatabase+"/dfile.bcc"
  elif ngroup>=225 and ngroup<=230:
    dfile = plotdatabase+"/dfile.fcc"
    
  if dfile != "":
    dfile0 = dfile.split('/')[-1]
    copyfile(dfile,dfile0)
    cmd = "Yphon -tranI 2 -eps -pdis "+dfile0+ " <superfij.out"
    if os.path.exists('dielecfij.out') : cmd = cmd + ' -Born dielecfij.out -bvec'
    print(cmd)
    output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                      universal_newlines=True)
    cmd = "gnuplot vdis.plt; convert -flatten -rotate 90 -density 120x120 vdis.eps vdis.png"
    output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                      universal_newlines=True)
    os.rename("vdis.eps", cwd+'/'+folder+'/vdis298.15.eps')
    os.rename("vdis.png", cwd+'/'+folder+'/vdis298.15.png')
  os.chdir( cwd )


def PlotVol(folder, vdos):
  import dfttk.scripts.config_dfttk as dfttkconfig
  PATH_TO_STORE_CONFIG = dfttkconfig.default_path()
  plotdatabase = dfttkconfig.get_abspath(PATH_TO_STORE_CONFIG)+'/analysis/database/'
  #print (plotdatabase, folder)
  vdosdir = [substr for substr in vdos.split('/') if substr!=""]
  if vdos.startswith('/'):
    vdosdir = '/'+('/').join(vdosdir[0:-1])
  else:
    vdosdir = ('/').join(vdosdir[0:-1])
  if not os.path.exists(vdosdir+'/superfij.out') : return

  for ff in ['POSCAR', 'CONTCAR', 'Symmetry.pos']:
    if os.path.exists(vdosdir+'/'+ff) :
      structure = Structure.from_file(vdosdir+'/'+ff)
      break
  if structure==None: return

  try:
    natom = len(structure.sites)
    sa = SpacegroupAnalyzer(structure)
    ngroup = sa.get_space_group_number()
  except:
    return
      
  #print(natom,ngroup)

  cwd = os.getcwd()
  os.chdir( vdosdir )
  cmd = "Yphon -tranI 2 -eps -nqwave "+ str(nqwave)+ " <superfij.out"
  if os.path.exists('dielecfij.out') : cmd = cmd + ' -Born dielecfij.out'
  print(cmd)
  output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    universal_newlines=True)
  cmd = "gnuplot vdos.plt; convert -flatten -rotate 90 -density 120x120 vdos.eps vdos.png"
  #print(cmd)
  output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    universal_newlines=True)
  os.rename("vdos.eps", cwd+'/'+folder+'/vdos.eps')
  os.rename("vdos.png", cwd+'/'+folder+'/vdos.png')

  if not os.path.exists('symmetry.mode'):
      cmd = "pos2s Symmetry.pos -THR 0.001"
      print(cmd)
      output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    universal_newlines=True)
  """
  #temp for debug
  """
  #temp for debug

  if os.path.exists("vdos.out") :
    os.rename("vdos.out", 'vdos.sav')
    cmd = "Yphon -tranI 2 -eps -nqwave 100 -Gfile symmetry.mode <superfij.out >Raman.mode"
    if os.path.exists('dielecfij.out') : cmd = cmd + ' -Born dielecfij.out'
    output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    universal_newlines=True)
    os.rename("vdos.sav", 'vdos.out')
    vdos = np.loadtxt("vdos.out", comments="#", dtype=np.float)
    if os.path.exists("Raman.mode") :
      with open ("Raman.mode", "r") as fp:
        plotRaman(cwd+'/'+folder, fp, vdos)      

  dfile = ""
  if ngroup>=1 and ngroup<=2:
    dfile = plotdatabase+"/dfile.tri"
  elif ngroup>=3 and ngroup<=15:
    dfile = plotdatabase+"/dfile.mon"
  elif ngroup>=16 and ngroup<=74:
    dfile = plotdatabase+"/dfile.oth"
  elif ngroup>=75 and ngroup<=142:
    dfile = plotdatabase+"/dfile.tet"
  elif ngroup>=143 and ngroup<=167:
    dfile = plotdatabase+"/dfile.rho"
  elif ngroup>=168 and ngroup<=194:
    dfile = plotdatabase+"/dfile.hcp"
  elif ngroup>=195 and ngroup<=220:
    dfile = plotdatabase+"/dfile.scc"
  elif ngroup>=221 and ngroup<=224:
    dfile = plotdatabase+"/dfile.bcc"
  elif ngroup>=225 and ngroup<=230:
    dfile = plotdatabase+"/dfile.fcc"
    
  if dfile != "":
    dfile0 = dfile.split('/')[-1]
    copyfile(dfile,dfile0)
    cmd = "Yphon -tranI 2 -eps -pdis "+dfile0+ " <superfij.out"
    if os.path.exists('dielecfij.out') : cmd = cmd + ' -Born dielecfij.out'
    output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                      universal_newlines=True)
    cmd = "gnuplot vdis.plt; convert -flatten -rotate 90 -density 120x120 vdis.eps vdis.png"
    output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                      universal_newlines=True)
    os.rename("vdis.eps", cwd+'/'+folder+'/vdis.eps')
    os.rename("vdis.png", cwd+'/'+folder+'/vdis.png')
  os.chdir( cwd )


"""make the reduced formula
els - a list of elements
natype - a list of number of elements
return:
"""
def reduced_formula(formula):
  _els, _nat = formula2composition(formula)
  els = sorted(set(_els))
  nat = np.zeros(len(els),dtype=int)
  for i,el in enumerate(_els):
    ix = els.index(el)
    nat[ix] += _nat[i]

  Nd = min(nat)
  for i in range(Nd,0,-1):
    out = True
    for j in range(len(nat)):
      if ((nat[j]//i)*i!=nat[j]):
        out = False
        break
    if out:
      break
  form = ""
  for j,el in enumerate(els):
    ix = nat[j]//i
    form = form+el
    if ix!=1:
      form = form+str(ix)
  return form

if __name__ == '__main__':
    formula = None
    expt = None
    plotlabel = None
    count = 1
    while (count < len(sys.argv)):
      if (input_within):
        if sys.argv[count].startswith('-'):
          input_within = False
        else:
          formula_within = formula_within+sys.argv[count]
          count = count + 1
          if (count > len(sys.argv)):
            break
          continue
      if (sys.argv[count] == "-pvdos"):
        pvdos = True
      if (sys.argv[count] == "-within"):
        input_within = True
      elif (sys.argv[count] == "-formula"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        formula = sys.argv[count]
      elif (sys.argv[count] == "-expt"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        expt = sys.argv[count]
      elif (sys.argv[count] == "-phdft"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        phdft = sys.argv[count]
      elif (sys.argv[count] == "-xlim"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        xlim = float(sys.argv[count])
      elif (sys.argv[count] == "-T0"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        T0 = float(sys.argv[count])
      elif (sys.argv[count] == "-phasename"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        phasename = str(sys.argv[count])
      elif (sys.argv[count] == "-cpmax"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        CpMax = float(sys.argv[count])
      elif (sys.argv[count] == "-THRE0"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        THRE0 = float(sys.argv[count])
      elif (sys.argv[count] == "-PQ"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        PQ = float(sys.argv[count])
      elif (sys.argv[count] == "-EQ"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        EQ = float(sys.argv[count])
      elif (sys.argv[count] == "-plot"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        plotlabel = sys.argv[count]
      elif (sys.argv[count] == "-Tupmax"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        Tupmax = float(sys.argv[count])
      elif (sys.argv[count] == "-nqwave"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        nqwave = float(sys.argv[count])
      elif (sys.argv[count] == "-fitG"):
        fitCp = False
      elif (sys.argv[count] == "-fitCp"):
        fitCp = True
      elif (sys.argv[count] == "-debug"):
        debug = True
      elif (os.path.exists(sys.argv[count])):
        justplot=sys.argv[count]
      else:
        print ("*******Unknown option", sys.argv[count])
      count = count + 1
    
    if formula_within!="":
      within = formula2elist(formula_within)
      print ("data to be extracted within ",within)
    
    #print (phasename)
    #if True:
    try:
      with open (phasename,'r') as f:
        lines = f.readlines()
      for ll in lines:
        line = ll.strip('\n').replace(',', ' ').replace(':', ' ')
        ss = [s.strip() for s in line.split(' ') if s.strip()!='']
        if len(ss)>1: PhaseName.update({ss[0]:ss[1]})
    except:
      pass
    """
    """
    
    if expt!=None: expt=get_expt(expt, formula)

    if justplot==None: lines = sys.stdin.readlines()
    else: 
        plotCMD(justplot, volumes=None, energies=None, expt=expt, xlim=xlim, _fitCp = fitCp, plotlabel=plotlabel)
        sys.exit()
    
    sys.stdout.write("G(T)=a+b*T+c*T*Ln(T)+d*T*T+e*T*T*T+f/T (J/mol-atom)\n")
    sys.stdout.write("Phase,comp")
    for ss in outs:
      sys.stdout.write(",{}".format(ss))
    #sys.stdout.write(",a,b,c,d,e,f,PQ(%),EQ(J),GQ(J),\n".format(ss))
    sys.stdout.write(",a,b,c,d,e,f,PQ(%),\n".format(ss))
    
    for line in lines:
      if line.strip()=="": continue
      skip,vdos_e,dir0,oldPN = mkDict(line)
      if skip: continue
      #phases.append([threcord.get("phase name"),threcord.get("mpid")])
      phases.append(threcord.get("Phase name"))
      nphases += 1
      #print (line.strip())
    
      if not os.path.exists(dir0):
        os.mkdir(dir0)
      else:
        if not update: continue
    
      Vfiles,Pfiles,g = Genergy(vdos_e,dir0)
      addpapers(g,dir0,oldPN)
    
      if not debug:
        VASPResults(dir0,vdos_e,Vfiles, Pfiles, phdft=phdft)
        try:
          Phonon298(dir0, pvdos=pvdos)
        except:
          pass
      threcord.update({"structure":structure})
      #threcord.delete({"number of atoms in the primitive unit cell")
    
      with open(dir0 + '/record.json', 'w') as fp:
        myjsonout(threcord, fp, indent="", comma="")
    
    print ("\n", phases, "\n")
    print ("\n", nphases, "phases extracted\n")

#!python
#
import pytest

import warnings
import datetime
import subprocess
import os
import json
import numpy as np
import copy
import six
import shlex
from phonopy.interface.vasp import Vasprun as PhonopyVasprun
from pymatgen.core import Structure
from pymatgen.io.vasp.inputs import Incar
from pymatgen.io.vasp.outputs import Vasprun, Outcar
from custodian.custodian import Custodian
from custodian.vasp.handlers import VaspErrorHandler, AliasingErrorHandler, \
    MeshSymmetryErrorHandler, UnconvergedErrorHandler, PotimErrorHandler, \
    FrozenJobErrorHandler, NonConvergingErrorHandler, PositiveEnergyErrorHandler, \
    StdErrHandler, DriftErrorHandler
from custodian.vasp.jobs import VaspJob
from pymatgen.analysis.eos import Vinet, EOS
from fireworks import explicit_serialize, FiretaskBase, FWAction
from atomate.utils.utils import load_class, env_chk
from atomate.vasp.database import VaspCalcDb
from dfttk.analysis.phonon import get_f_vib_phonopy, get_phonon_band_dos, phonon_stable
from dfttk.analysis.relaxing import get_non_isotropic_strain, get_bond_distance_change
from dfttk.analysis.quasiharmonic import Quasiharmonic
from dfttk.utils import sort_x_by_y, update_pos_by_symbols, update_pot_by_symbols, check_symmetry
from dfttk.custodian_jobs import ATATWalltimeHandler, ATATInfDetJob
from atomate import __version__ as atomate_ver
from dfttk import __version__ as dfttk_ver
from pymatgen.core import __version__ as pymatgen_ver

from dfttk.EVcheck_QHA import *
from pymatgen.core import Structure
from pymatgen.analysis.eos import EOS
from fireworks import Firework
from atomate.vasp.config import VASP_CMD, DB_FILE
import os
from dfttk.ftasks import QHAAnalysis
from dfttk.analysis.debye import DebyeModel

POSCAR_STR = """fcc Al
4.10
0 .5 .5
.5 0 .5
.5 .5 0
Al
1
direct
0.000000 0.000000 0.000000 Al"""
structure = Structure.from_str(POSCAR_STR, fmt='POSCAR')
energies = [-3.62297814, -3.69619945, -3.73411839, -3.74537143, -3.73636044, -3.71218262, -3.676796, -3.63402384]
volumes = [14.004168009146943, 14.827947688018714, 15.651725983766324, 16.47548996484071, 17.29925450808183, 18.12305110099911, 18.946820774727207, 19.754109316522513]
D_vib = [534.74656038, 476.9830257, 428.09760996, 386.36052023, 350.44243067, 319.30871742, 292.14776627, 268.75843483]

@pytest.mark.DebyeModel
def test_DebyeModel():
    debye_model = DebyeModel(energies, volumes, structure, t_min=5, t_step=5,
        t_max=1000, eos="vinet", poisson=0.363615, bp2gru=2./3.)
    print (debye_model.D_vib)
    assert np.isclose(debye_model.D_vib, D_vib).all()

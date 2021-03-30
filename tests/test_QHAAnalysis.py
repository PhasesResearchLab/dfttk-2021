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

head,tail = os.path.split(__file__)
db_file = os.path.join(head,"db.json")

@pytest.mark.QHAAnalysis
def test_QHAAnalysis():
    tags = ['ec77b415-8e36-440a-997c-1c3d512099ce',
            '7e74e496-90b1-41c1-8113-38eec730b9f2']
    for tag in tags:
        vasp_db = VaspCalcDb.from_db_file(db_file, admin=False)
        phonon_calculations = list(vasp_db.db['phonon'].find({'$and':[ {'metadata.tag': tag}, {'adopted': True} ]}))
        T = phonon_calculations[0]['temperatures']
        t_min = min(T)
        everyT = 10
        t_step = (T[1]-T[0])*everyT
        t_max = max(T)
        nT = int((t_max-t_min)/t_step)
        t_max = t_min + nT*t_step
        print ("testing tag=", tag)
        proc = QHAAnalysis(phonon=True, t_min=t_min, t_max=t_max,
            t_step=t_step, everyT=everyT, db_file=db_file, test=True,
            tag=tag)
        proc.run_task("")


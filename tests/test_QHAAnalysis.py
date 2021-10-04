#!python
#
import pytest

import warnings
import datetime
import subprocess
import os
import sys
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
from dfttk.pythelec import get_static_calculations
from pymatgen.core import Structure
from pymatgen.analysis.eos import EOS
from fireworks import Firework
from atomate.vasp.config import VASP_CMD, DB_FILE
import os
from dfttk.ftasks import QHAAnalysis
from dfttk.analysis.ywutils import get_rec_from_metatag, get_used_pot
from monty.serialization import loadfn, dumpfn

FW_CONFIG_FILE = os.getenv('FW_CONFIG_FILE')
FW_CONFIG_FILE = loadfn(FW_CONFIG_FILE)

CONFIG_FILE_DIR = FW_CONFIG_FILE.get('CONFIG_FILE_DIR',None)
if CONFIG_FILE_DIR is not None:
    db_file = os.path.join(CONFIG_FILE_DIR, 'db.json')
else:
    head,tail = os.path.split(__file__)
    db_file = os.path.join(head,"db.json")

@pytest.mark.QHAAnalysis
def test_QHAAnalysis():
    tags = ['5e8b5a18-b2b9-4dcd-81a7-0bd75ee361ec']
    print(db_file)
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

@pytest.mark.check_points
def test_check_points():
    proc = EVcheck_QHA()
    volumes = [504.4854780252552, 534.1610477526048, 563.8366733894931, 593.512427132093, 623.1879003769909, 652.8636686613459, 682.5393342833485, 711.6214250632727, 740.703761424463, 769.7858059616457, 798.8678864247061]
    energies = [-342.98343642, -344.38954914, -344.89492019, -344.80547287, -344.28500957, -343.47527914, -342.47536783, -341.31384393, -340.16711161, -339.67940075, -338.41924052]
    proc.check_points("", "", 0.005, 14, 0.3, volumes, energies, True)

@pytest.mark.check_points_1
def test_check_points_1():
    tag = '19c9e217-4159-4bfe-9c3a-940fb40e023e'
    tag = 'd054780c-f051-4450-a611-d374d41d1884'
    tag = 'ed85a69b-4054-41d4-a724-7373934cdcc6'
    proc = EVcheck_QHA()
    volumes, energies, _ = proc.get_orig_EV(db_file, tag)
    proc.check_points("", "", 0.005, 14, 0.3, volumes, energies, True)
    #assert False

@pytest.mark.check_points_2
def test_check_points_2():
    tag = '19c9e217-4159-4bfe-9c3a-940fb40e023e'
    tag = 'd054780c-f051-4450-a611-d374d41d1884'
    tag = 'ed85a69b-4054-41d4-a724-7373934cdcc6'
    tag = '0267947f-b5a0-4993-8b5e-58f9dfb4e362'
    vasp_db = VaspCalcDb.from_db_file(db_file, admin=False)
    EV, POSCAR, INCAR = get_rec_from_metatag(vasp_db, tag, test=True)
    structure = Structure.from_str(POSCAR, fmt='POSCAR')
    proc = EVcheck_QHA(db_file=db_file, metadata={'tag':tag}, structure=structure, deformations = np.linspace(0.94,1.06,7), test=True)
    proc.run_task({})
    #assert False

@pytest.mark.EVcheck_QHA_2
def test_EVcheck_QHA_2():
    tag = '19c9e217-4159-4bfe-9c3a-940fb40e023e'
    wf = Firework(EVcheck_QHA(db_file=db_file,vasp_cmd=">>vasp_cmd<<",tag="test",metadata={'tag':tag}, test=True))
    print(wf.as_dict())
    #assert False

@pytest.mark.get_static_calculations
def test_get_static_calculations():
    tag = '0267947f-b5a0-4993-8b5e-58f9dfb4e362'
    vasp_db = VaspCalcDb.from_db_file(db_file=db_file, admin = True)
    volumes, energies, dos_objs, _ = get_static_calculations(vasp_db,tag)
    #print(volumes, energies)
    assert False

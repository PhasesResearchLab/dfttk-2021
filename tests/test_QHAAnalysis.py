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

from dfttk.wflows import _get_deformations
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


if os.getenv("GITHUB_ACTIONS") == "true":
    print("Skipping the test when running inside GitHub actions and Yi Wang's MongoDB is not accessible")
else:
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
        tags = ['700dd9c9-afb0-4a09-a39d-ad77622b3f08']
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
        tag = '8e7b216d-c6b6-4da4-905e-e7afd44195aa'
        vasp_db = VaspCalcDb.from_db_file(db_file, admin=False)
        EV, POSCAR, INCAR = get_rec_from_metatag(vasp_db, tag, test=True)
        structure = Structure.from_str(POSCAR, fmt='POSCAR')
        proc = EVcheck_QHA(db_file=db_file, metadata={'tag':tag}, structure=structure, deformations = np.linspace(0.94,1.06,7), test=True)
        proc.run_task({})
        #assert False

    POSCAR_STR = """STa1
    4.10
    0 .5 .5
    .5 0 .5
    .5 .5 0
    Al
    1
    direct
    0.000000 0.000000 0.000000 Al"""


    @pytest.mark.check_points_3
    def test_check_points_3():
        settings = {'deformation_fraction': [-0.05, 0.05], 'num_deformations': 3, 'deformation_scheme': 'volume', 'run_isif': 4, 'phonon': True, 'phonon_supercell_matrix': [[-2, 2, 2], [2, -2, 2], [2, 2, -2]], 'override_default_vasp_params': {'user_incar_settings': {'EDIFF_PER_ATOM': 1e-07, 'Relax_settings': {'PREC': 'HIGH', 'grid_density': 1000}, 'Static_settings': {'PREC': 'HIGH', 'grid_density': 1000}}, 'user_kpoints_settings': {'grid_density': 1000}}, 'metadata': {'tag': '67783198-6d5a-40c2-8d49-4f03e50ac130'}}

        ################ PARAMETERS FOR WF #############################
        #str, the absolute path of db.json file, e.g. /storage/home/mjl6505/atomate/config/db.json
        #  If None, it will use the configuration in fireworks
        #db_file = settings.get('db_file', DB_FILE)
        db_file = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]
        #str, the vasp command, if None then find in the FWorker configuration
        vasp_cmd = settings.get('vasp_cmd', VASP_CMD)
        #dict, metadata to be included, this parameter is useful for filter the data, e.g. metadata={"phase": "BCC_A2", "tag": "AFM"}
        metadata = settings.get('metadata',{})
        tag = metadata.get('tag', '{}'.format(str(uuid4())))
        metadata.update({'tag': tag})

        #list/tuple(min, max) or float(-max, max), the maximum amplitude of deformation, e.g. (-0.15, 0.15) means (0.95, 1.1) in volume
        deformation_fraction = settings.get('deformation_fraction', (-0.15, 0.20))
        #int, the number of initial deformations, e.g. 7
        num_deformations = settings.get('num_deformations', 8)
        if num_deformations==1:
            deformation_fraction[1] = deformation_fraction[0]
        #bool, run phonon(True) or not(False)
        phonon = settings.get('phonon', False)
        #list(3x3), the supercell matrix for phonon, e.g. [[2.0, 0, 0], [0, 2.0, 0], [0, 0, 2.0]]
        phonon_supercell_matrix = settings.get('phonon_supercell_matrix', None)
        phonon_supercell_matrix_min = settings.get('phonon_supercell_matrix_min', None)
        phonon_supercell_matrix_max = settings.get('phonon_supercell_matrix_max', None)
        optimize_sc = settings.get('optimize_sc', False)
        #The tolerance for phonon stable
        stable_tor = settings.get('stable_tor', 0.01)
        #float, the mimimum of temperature in QHA process, e.g. 5
        t_min = settings.get('t_min', 5)
        #float, the maximum of temperature in QHA process, e.g. 2000
        t_max = settings.get('t_max', 2000)
        #float, the step of temperature in QHA process, e.g. 5
        t_step = settings.get('t_step', 5)
        #float, acceptable value for average RMS, recommend >= 0.005
        eos_tolerance = settings.get('eos_tolerance', 0.01)

        #Global settings for all vasp job, e.g.
        #override_default_vasp_params = {'user_incar_settings': {}, 'user_kpoints_settings': {}, 'user_potcar_functional': str}
        #If some value in 'user_incar_settings' is set to None, it will use vasp's default value
        override_default_vasp_params = settings.get('override_default_vasp_params', {})

        #dict, dict of class ModifyIncar with keywords in Workflow name. e.g.
        modify_incar_params = settings.get('modify_incar_params', {})

        #check if fworker_name is assigned
        powerups = settings.get('powerups', {})
        if len(powerups)>0:
            if 'user_incar_settings' not in override_default_vasp_params:
                override_default_vasp_params.update({'user_incar_settings':{}})
            override_default_vasp_params['user_incar_settings'].update({'powerups':powerups})
            modify_incar_params.update({'powerups':powerups})
        
        #dict, dict of class ModifyKpoints with keywords in Workflow name, similar with modify_incar_params
        modify_kpoints_params = settings.get('modify_kpoints_params', {})
        #bool, print(True) or not(False) some informations, used for debug
        verbose = settings.get('verbose', False)
        #Save the volume data or not ("chgcar", "aeccar0", "aeccar2", "elfcar", "locpot")
        store_volumetric_data = settings.get('store_volumetric_data', False)

        if phonon:
            if isinstance(phonon_supercell_matrix, str):
                if phonon_supercell_matrix=='Yphon':
                    phonon_supercell_matrix = supercell_scaling_by_Yphon(structure, 
                                                supercellsize=phonon_supercell_matrix_max)
                else:
                    phonon_supercell_matrix = supercell_scaling_by_atom_lat_vol(structure, min_obj=phonon_supercell_matrix_min,
                                                max_obj=phonon_supercell_matrix_max, scale_object=phonon_supercell_matrix,
                                                target_shape='sc', lower_search_limit=-2, upper_search_limit=2,
                                                verbose=verbose, sc_tolerance=1e-5, optimize_sc=optimize_sc)

        _deformations = _get_deformations(deformation_fraction, num_deformations)
        if num_deformations > 1: vol_spacing = _deformations[1]-_deformations[0]
        else: vol_spacing=0.05

        structure = Structure.from_str(POSCAR_STR, fmt='POSCAR')
        t_kwargs = {'t_min': t_min, 't_max': t_max, 't_step': t_step}
        common_kwargs = {'vasp_cmd': vasp_cmd, 'db_file': db_file, "metadata": metadata, "tag": tag,
                        'override_default_vasp_params': override_default_vasp_params}
        vasp_kwargs = {'modify_incar_params': modify_incar_params, 'modify_kpoints_params': modify_kpoints_params}
        eos_kwargs = {'deformations': _deformations, 'vol_spacing': vol_spacing, 'eos_tolerance': eos_tolerance, 'threshold': 14}
        a_kwargs = {"structure":structure, "settings":settings, "eos_kwargs":eos_kwargs,
            "static": True, "phonon":phonon, "phonon_supercell_matrix":phonon_supercell_matrix}

        proc = Crosscom_EVcheck_QHA(verbose=verbose, stable_tor=stable_tor,
                run_num = 0,
                store_volumetric_data=store_volumetric_data, a_kwargs=a_kwargs,
                **eos_kwargs, **vasp_kwargs, **t_kwargs, **common_kwargs, test=True)
        proc.run_task({})
        assert False

    @pytest.mark.EVcheck_QHA_2
    def test_EVcheck_QHA_2():
        tag = '19c9e217-4159-4bfe-9c3a-940fb40e023e'
        wf = Firework(EVcheck_QHA(db_file=db_file,vasp_cmd=">>vasp_cmd<<",tag="test",metadata={'tag':tag}, test=True))
        print(wf.as_dict())
        #assert False

    @pytest.mark.get_static_calculations
    def test_get_static_calculations():
        tag = '8e7b216d-c6b6-4da4-905e-e7afd44195aa'
        vasp_db = VaspCalcDb.from_db_file(db_file=db_file, admin = True)
        volumes, energies, dos_objs, _ = get_static_calculations(vasp_db,tag)
        #print(volumes, energies)
        assert False

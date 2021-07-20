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
from dfttk.pythelec import get_static_calculations

head,tail = os.path.split(__file__)
db_file = os.path.join(head,"db.json")


@explicit_serialize
class QHAAnalysis_renew(FiretaskBase):
    """
    Do the quasiharmonic calculation from either phonon or Debye.

    Required params
    ---------------
    tag : str
        Tag to search the database for static calculations (energies, volumes, eDOS) from this job.
    db_file : str
        Points to the database JSON file. If None (the default) is passed, the path will be looked up in the FWorker.
    phonon : bool
        True if f_vib comes from phonon calculations (in the spec). If False, it is calculated by the Debye model.
    t_min : float
        Minimum temperature
    t_step : float
        Temperature step size
    t_max : float
        Maximum temperature (inclusive)

    Optional params
    ---------------
    poisson : float
        Poisson ratio, defaults to 0.25. Only used in Debye
    bp2gru : float
        Debye model fitting parameter for dBdP in the Gruneisen parameter. 2/3 is the high temperature
        value and 1 is the low temperature value. Defaults to 1.

    Notes
    -----
    Heterogeneity in the sources of E_0/F_el and F_vib is solved by sorting them according to increasing volume.
    """

    required_params = ["phonon", "db_file", "t_min", "t_max", "t_step", "tag"]

    optional_params = ["poisson", "bp2gru", "metadata", "test_failure", "admin"]

    def run_task(self):
        tag = self["tag"]
        admin = self.get('admin', False)
        _db_file = self.get('db_file', db_file)
        vasp_db = VaspCalcDb.from_db_file(db_file=_db_file, admin=admin)
        volumes, energies, dos_objs, _calc = get_static_calculations(vasp_db, tag)
        structure = Structure.from_dict(_calc['output']['structure'])

        qha_result = {}
        qha_result['structure'] = structure.as_dict()
        qha_result['formula_pretty'] = structure.composition.reduced_formula
        qha_result['elements'] = sorted([el.name for el in structure.composition.elements])
        qha_result['metadata'] = self.get('metadata', {'tag':tag})
        #qha_result['has_phonon'] = self['phonon']

        poisson = self.get('poisson', 0.363615)
        bp2gru = self.get('bp2gru', 2./3.)

        # phonon properties
        # check if phonon calculations existed
        #always perform phonon calculations when when enough phonon calculations found
        num_phonons = len(list(vasp_db.db['phonon'].find({'$and':[ {'metadata': {'tag':tag}}, {'adopted': True} ]})))       
        #num_phonons = len(list(vasp_db.db['phonon'].find({'$and':[ {'metadata.tag': tag}, {'adopted': True} ]})))       
        qha_result['has_phonon'] = num_phonons >= 5
        #if self['phonon']:
        if qha_result['has_phonon']:
            # get the vibrational properties from the FW spec
            phonon_calculations = list(vasp_db.db['phonon'].find({'$and':[ {'metadata.tag': tag}, {'adopted': True} ]}))
            #phonon_calculations = list(vasp_db.db['phonon'].find({'$and':[ {'metadata': {'tag':tag}}, {'adopted': True} ]}))
            vol_vol = []
            vol_f_vib = []
            vol_s_vib = []
            vol_c_vib = []
            for calc in phonon_calculations:
                if calc['volume'] in vol_vol: continue
                if calc['volume'] not in volumes: continue
                vol_vol.append(calc['volume'])
                vol_f_vib.append(calc['F_vib'])
                vol_s_vib.append(calc['S_vib'])
                vol_c_vib.append(calc['CV_vib'])
            # sort them order of the unit cell volumes
            vol_f_vib = sort_x_by_y(vol_f_vib, vol_vol)
            vol_s_vib = sort_x_by_y(vol_s_vib, vol_vol)
            vol_c_vib = sort_x_by_y(vol_c_vib, vol_vol)
            f_vib = np.vstack(vol_f_vib)

            # by Yi Wang, after a long day debug, finally I fixex the bug below
            # i.e, sometimes, the number of phonon volumes is less than that of static!
            _volumes = []
            _energies = []
            _dos_objs = []
            for iv,vol in enumerate(volumes):
                if vol not in vol_vol:  continue
                _volumes.append(vol)
                _energies.append(energies[iv])
                _dos_objs.append(dos_objs[iv])
            volumes = _volumes
            energies = _energies
            dos_objs = _dos_objs      
               
            qha = Quasiharmonic(energies, volumes, structure, dos_objects=dos_objs, F_vib=f_vib,
                                t_min=self['t_min'], t_max=self['t_max'], t_step=self['t_step'],
                                poisson=poisson, bp2gru=bp2gru)
            qha_result['phonon'] = qha.get_summary_dict()
            qha_result['phonon']['entropies'] = vol_s_vib
            qha_result['phonon']['heat_capacities'] = vol_c_vib
            qha_result['phonon']['temperatures'] = qha_result['phonon']['temperatures'].tolist()

        # calculate the Debye model results no matter what
        qha_debye = Quasiharmonic(energies, volumes, structure, dos_objects=dos_objs, F_vib=None,
                                  t_min=self['t_min'], t_max=self['t_max'], t_step=self['t_step'],
                                  poisson=poisson, bp2gru=bp2gru)

        # fit 0 K EOS for good measure
        eos = Vinet(volumes, energies)
        eos.fit()
        errors = eos.func(volumes) - energies
        sum_square_error = float(np.sum(np.square(errors)))
        eos_res = {}
        eos_res['b0_GPa'] = float(eos.b0_GPa)
        eos_res['b0'] = float(eos.b0)
        eos_res['b1'] = float(eos.b1)
        eos_res['eq_volume'] = float(eos.v0)
        eos_res['eq_energy'] = float(eos.e0)
        eos_res['energies'] = energies
        eos_res['volumes'] = volumes
        eos_res['name'] = 'Vinet'
        eos_res['error'] = {}
        eos_res['error']['difference'] = errors.tolist()  # volume by volume differences
        eos_res['error']['sum_square_error'] = sum_square_error
        qha_result['eos'] = eos_res

        qha_result['debye'] = qha_debye.get_summary_dict()
        qha_result['debye']['poisson'] = poisson
        qha_result['debye']['bp2gru'] = bp2gru
        qha_result['debye']['temperatures'] = qha_result['debye']['temperatures'].tolist()

        qha_result['version_atomate'] = atomate_ver
        qha_result['version_dfttk'] = dfttk_ver
        volumes_false = []
        energies_false = []
        static_falses = vasp_db.collection.find({'$and':[ {'metadata.tag': tag}, {'adopted': False} ]})
        for static_false in static_falses:
            volumes_false.append(static_false['output']['structure']['lattice']['volume'])
            energies_false.append(static_false['output']['energy'])
        qha_result['Volumes_fitting_false'] = volumes_false
        qha_result['Energies_fitting_false'] = energies_false
        print('Volumes_fitting_false : %s' %volumes_false)
        print('Energies_fitting_false: %s' %energies_false)
        print('number of phonon calculations found : %s' %num_phonons)

        """
        # write to JSON for debugging purposes
        import json
        with open('qha_summary.json', 'w') as fp:
            json.dump(qha_result, fp, indent=4)
        """

        if self.get("test_failure", False) : return
        #if self['phonon']:
        if qha_result['has_phonon']:
            vasp_db.db['qha_phonon'].insert_one(qha_result)
            qha_result.pop('phonon')
            vasp_db.db['qha'].insert_one(qha_result)
        else:
            vasp_db.db['qha'].insert_one(qha_result)



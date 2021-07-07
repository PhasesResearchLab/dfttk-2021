from pymatgen.core import Structure
import pymatgen
from pymatgen.io.vasp.inputs import Incar
from fireworks import FWorker, Workflow, LaunchPad
from dfttk import get_wf_gibbs_robust
from dfttk.utils import update_fws_spec
import pytest
import shutil
import os
import time

MODULE_DIR = os.path.dirname(__file__)

# TODO: does not support use_fake_vasp yet. VASP will fail.

DEBUG_MODE = False

POSCAR_STR = """Si2
1.0
3.840198 0.000000 0.000000
1.920099 3.325710 0.000000
0.000000 -2.217138 3.135509
Si
2
direct
0.000000 0.000000 0.000000 Si
0.750000 0.500000 0.750000 Si"""
STRUCT = Structure.from_str(POSCAR_STR, fmt='POSCAR')
TEST_DIR = os.path.join(MODULE_DIR, 'tmp_fw_test_dir')
# LPAD = LaunchPad.from_dict({'host': 'localhost', 'logdir': None, 'name': 'dfttk_unittest', 'password': None, 'port': 27017, 'ssl_ca_file': None, 'strm_lvl': 'DEBUG', 'user_indices': [], 'username': None, 'wf_user_indices': []})

# TODO: enable debug mode by having a launchpad that does not reset
# Can this be done by still having other tests pass?
# Should we only run one test?
# Stop on failure?
@pytest.fixture
def lpad():
    """A LaunchPad object for test instances to use. Always gives a clean (reset) LaunchPad. """
    LPAD.reset(None, require_password=False, max_reset_wo_password=5)
    yield LPAD
    LPAD.connection.close()
    return

@pytest.fixture
def launch_dir():
    test_dir = TEST_DIR + '-' + str(time.time()).split('.')[0]
    os.mkdir(test_dir)
    os.chdir(test_dir)
    yield test_dir
    os.chdir('..')
    shutil.rmtree(test_dir)
    return


@pytest.fixture
def patch_pmg_psp_dir():
    current_psp_dir = pymatgen.SETTINGS.get('PMG_VASP_PSP_DIR')
    if current_psp_dir is None:
        pymatgen.SETTINGS['PMG_VASP_PSP_DIR'] = os.path.join(MODULE_DIR, 'test_potcars')
    yield
    pymatgen.SETTINGS['PMG_VASP_PSP_DIR'] = current_psp_dir


@pytest.fixture(scope='module')
def launch_dir_debug():
    test_dir = TEST_DIR + '-' + str(time.time()).split('.')[0]
    os.mkdir(test_dir)
    os.chdir(test_dir)
    yield test_dir
    os.chdir('..')
    return

if DEBUG_MODE:
    launch_dir = launch_dir_debug


@pytest.fixture
def fworker():
    scratch_dir = os.path.join(MODULE_DIR, 'scratch_dir')
    os.mkdir(scratch_dir)
    yield FWorker(env={"db_file": os.path.join(MODULE_DIR, "db.json"), 'scratch_dir': scratch_dir})
    shutil.rmtree(scratch_dir)


def test_fw_spec_modified_by_powerup():
    wf = get_wf_gibbs_robust(STRUCT, db_file=os.path.join(MODULE_DIR, "db.json"))
    wf = update_fws_spec(wf, {'_preserve_fworker': True})
    assert all([fw.spec['_preserve_fworker'] == True for fw in wf.fws])


def test_gibbs_wf_fireworks_graph():
    """Test that the graph of Fireworks is correct for a Gibbs workflow."""
    wf_phonon = get_wf_gibbs_robust(STRUCT, db_file=os.path.join(MODULE_DIR, "db.json"), phonon=True, \
        phonon_supercell_matrix=[[2,0,0],[0,2,0],[0,0,2]], num_deformations=11)
    assert len(wf_phonon.fws) == 4
    wf_debye = get_wf_gibbs_robust(STRUCT, db_file=os.path.join(MODULE_DIR, "db.json"), num_deformations=11, phonon=False)
    assert len(wf_debye.fws) == 3



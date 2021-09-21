import argparse
from dfttk.scripts.run_dfttk_ext import run_ext_thfind
from pymatgen.core import Structure
import pytest
import shutil
import os
import time
import copy
from os import walk
from atomate.vasp.database import VaspCalcDb
from dfttk.scripts.run_dfttk_ext import ext_EVfind, ext_thfind, ext_thelec


head,tail = os.path.split(__file__)
db_file = os.path.join(head,"db.json")
print("db_file", db_file)
vasp_db = VaspCalcDb.from_db_file(db_file, admin=False)


def _thargs():
    parser = run_ext_thfind(None)
    args, unknown = parser.parse_known_args() 
    return args


@pytest.mark.get_thargs
def test_args():
    args=_thargs()
    for arg in vars(args):
        print(arg,getattr(args, arg))


@pytest.fixture(scope="module", autouse=True)
def failure_tracking_fixture(request):
    tests_failed_before_module = request.session.testsfailed
    yield
    tests_failed_during_module = request.session.testsfailed - tests_failed_before_module


def run_thelec(capsys,tmp_path,compound, phasename, 
    ext_module=ext_thfind, metatag=None, nfiles=5):
    """Test that the postprocess module thfind."""
    thargs = _thargs()
    thargs.within = compound
    thargs.containall = compound
    thargs.debug = False
    thargs.get = True
    thargs.noel = False
    thargs.metatag = metatag
    thargs.plot = "DFT"
    cwd = os.getcwd()
    test_dir = os.path.join(tmp_path,"dfttk-test-" + str(time.time()).split('.')[0])
    os.mkdir(test_dir)
    os.chdir(test_dir)
    #yield test_dir
    #assert cwd == test_dir
    print("\n\n# Testing compound=",compound, "using database", db_file, \
        "if failed, see additional info at",test_dir)

    #vasp_db = VaspCalcDb.from_db_file(db_file, admin=False)
    ext_module(thargs, vasp_db=vasp_db)
    lines, _ = capsys.readouterr()
    line = "Full thermodynamic properties have outputed into: "
    kRec = [rec for rec in lines.split('\n') if rec.startswith(line)]

    assert len(kRec) >= 1
    _, _, filenames = next(walk(phasename))
    assert len(filenames) == nfiles
    _, _, filenames = next(walk(os.path.join(phasename,"figures")))
    os.chdir(cwd)
    assert len(filenames) >= 19
    #if not failure_tracking_fixture(request):
    shutil.rmtree(test_dir)

@pytest.mark.EVfind
def test_EVfind(capsys, tmp_path):
    """Test that the postprocess module thfind."""
    thargs = _thargs()
    thargs.print = False
    thargs.nV = 0
    cwd = os.getcwd()
    test_dir = os.path.join(tmp_path,"dfttk-test-" + str(time.time()).split('.')[0])
    os.mkdir(test_dir)
    os.chdir(test_dir)
    print("\n\n# Testing EVfind, if failed, see additional info at",test_dir)

    ext_EVfind(thargs, vasp_db)
    n_thfind, _ = capsys.readouterr()
    kRec = [rec for rec in n_thfind.split('\n') if rec.startswith("{'tag': '")]
    assert len(kRec) >= 3
    os.chdir(cwd)
    if len(kRec) >= 3:
        shutil.rmtree(test_dir)


@pytest.mark.thfind
def test_thfind(capsys):
    """Test that the postprocess module thfind."""
    thargs = _thargs()
    thargs.nV = 1
    ext_thfind(thargs, vasp_db=vasp_db)
    n_thfind, _ = capsys.readouterr()
    kRec = [rec for rec in n_thfind.split('\n') if rec.startswith("{'tag': '")]
    print ("testing thfind using the test database", db_file)
    assert len(kRec) >= 7


@pytest.mark.MgO
def test_thelec_MgO(capsys,tmp_path):
    compound, phasename ="MgO", "MgO_Fm-3m_225LDA"
    print ("testing thelec for", compound)
    run_thelec(capsys,tmp_path,compound, phasename, nfiles=8)


@pytest.mark.Al
def test_thelec_Al(capsys,tmp_path):
    compound, phasename ="Al", "Al_Fm-3m_225PBE"
    print ("testing thelec for", compound)
    run_thelec(capsys,tmp_path,compound, phasename, 
    ext_module=ext_thelec, metatag='12b9ba49-78f9-44b3-bbb0-b7b6b00e2d61',
    nfiles=5)


@pytest.mark.Al2O3
def test_thelec_Al2O3(capsys,tmp_path):
    compound, phasename = "Al2O3", "Al2O3_R-3c_167LDA"
    print ("testing thelec for", compound)
    run_thelec(capsys,tmp_path,compound, phasename, nfiles=8)


DIR = "/gpfs/scratch/yuw3/v7/block_2021-09-20-21-31-28-177520/launcher_2021-09-20-21-31-49-741147/launcher_2021-09-21-08-53-15-627218"
if not os.path.exists(DIR): DIR = None
@pytest.mark.check_symmetry_magmom
@pytest.mark.skipif(DIR is None, reason="Check magmom at the given path which needs to exist")
def test_check_symmetry_magmom():
    check_symmetry(tol_energy=0.025, tol_strain=0.05, tol_bond=0.10, site_properties=Ture)



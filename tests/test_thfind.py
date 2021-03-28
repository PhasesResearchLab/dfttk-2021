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
expt = os.path.join(head,"ExptData.json")
vasp_db = VaspCalcDb.from_db_file(db_file, admin=False)


class _thargs:
    def __init__(self):
        self.local = ""
        self.within = None
        self.containall = None
        self.excludeall = None
        self.excludeany = None
        self.containany = None
        self.nV = 6
        self.supercellN = 0
        self.findbandgap = False
        self.get = False
        self.check = False
        self.remove = False
        self.pyphon = False
        self.t0 = 0
        self.t1 = 4000
        self.td = 10
        self.xdn = -100
        self.xup = 100
        self.dope = 0
        self.ndosmx = 1001
        self.gaussian = 1000
        self.natom = 1
        self.nT = 257
        self.everyT = 1
        self.outf = "fvib_ele"
        self.noel = False
        self.metatag = None
        self.qhamode = None
        self.phasename = None
        self.jobpath = None
        self.eqmode = 4
        self.elmode = 0
        self.smooth = False
        self.plot = None
        self.renew = False
        self.refresh = False
        self.SGTEfitCp = False
        self.fitF = False
        self.plotonly = False
        self.debug = True
        self.expt = expt
        self.xlim = None
        self.doscar = None
        self.poscar = None
        self.vdos = None


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
    ext_module=ext_thelec, metatag='9e653d55-2766-48db-b2da-05933d31e5ea',
    nfiles=8)


@pytest.mark.Al2O3
def test_thelec_Al2O3(capsys,tmp_path):
    compound, phasename = "Al2O3", "Al2O3_R-3c_167LDA"
    print ("testing thelec for", compound)
    run_thelec(capsys,tmp_path,compound, phasename, nfiles=8)

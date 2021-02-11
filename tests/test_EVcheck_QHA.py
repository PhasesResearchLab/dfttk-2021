#!python
#
import pytest

from dfttk.EVcheck_QHA import *
from pymatgen import Structure
from pymatgen.analysis.eos import EOS
from fireworks import Firework
from atomate.vasp.config import VASP_CMD, DB_FILE

db_file = "db.json"

if db_file == ">>db_file<<":
    #In PengGao's version, some function used the absolute db_file
    from fireworks.fw_config import config_to_dict
    from monty.serialization import loadfn
    db_file = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]

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

def test_extract_accord_index():
    index = [0, 2]
    p_in = [[0, 1], [1, 2], [2, 3], [3, 4]]
    p_out = extract_accord_index(index, p_in)
    assert(p_out == [[0, 1], [2, 3]])

def test_gen_volenergdos():
    result_volume = [-0.01, 0.005]
    result_energy = [-2.0, -4.0]
    result_dos_obj = [1, 3]
    num = [0, 2]
    volumes = [-0.01, -0.005, 0.005, 0.01]
    energies = [-2.0, -3.0, -4.0, -5.0]
    volume, energy = gen_volenergdos(num, volumes, energies)
    assert(all(volume[i] == result_volume[i] for i in range(len(volume))))
    assert(all(energy[i] == result_energy[i] for i in range(len(energy))))
    dos_objs = [1, 2, 3, 4]
    volume, energy, dos_obj = gen_volenergdos(num, volumes, energies, dos_objs)
    assert(all(dos_obj[i] == result_dos_obj[i] for i in range(len(dos_obj))))

def test_check_deformations_in_volumes():
    deformations = [0.9, 0.95, 1.0, 1.05, 1.1]
    volumes = [31, 32, 33, 34, 35]
    orig_vol = 33
    deform_outof_vol = check_deformations_in_volumes(deformations, volumes, orig_vol=orig_vol)
    result = [0.9, 1.1]
    assert(all(deform_outof_vol[i] == result[i] for i in range(len(deform_outof_vol))))

def test_cal_stderr():
    value = [-34.69020102,
        -34.88230787, 
        -34.98749533, 
        -35.02251529, 
        -35.00101837, 
        -34.93416781, 
        -34.83111696, 
        -34.69802042]
    ref = [-34.69037673,
        -34.88126365,
        -34.98844492,
        -35.02444959,
        -34.99974727,
        -34.92873799,
        -34.8383195,
        -34.69550342]
    stderr = cal_stderr(value=value, ref=ref)
    assert(stderr == pytest.approx(1.18844e-5, abs=1e-7))

def test_eosfit_stderr():
    volume = [64.26025658624827,
        66.6402902061661,
        69.02026278463558,
        71.40028395651238,
        73.7802955243384,
        76.16031916837127,
        78.5402791170494,
        80.94268976633485]
    energy = [-34.69037673,
        -34.88126365,
        -34.98844492,
        -35.02444959,
        -34.99974727,
        -34.92873799,
        -34.8383195,
        -34.69550342]
    eos = EOS('vinet')
    eos_fit = eos.fit(volume, energy)
    fit_value = eos_fit.func(volume)
    print(fit_value)
    stderr = eosfit_stderr(eos_fit, volume, energy)
    assert(stderr == pytest.approx(1.18844e-5, abs=1e-7))

def test_EVcheck_QHA():
    wf = Firework(EVcheck_QHA(db_file=db_file,vasp_cmd=VASP_CMD,tag="test",metadata={}))
    #print(wf.as_dict())

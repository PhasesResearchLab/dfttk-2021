#!python
#
#Test for parse_anrl_prototype.py in structure_builders
#

import os
import sys
import pytest
from pymatgen import Structure
from monty.serialization import loadfn
try:
    import dfttk.structure_builders.parse_anrl_prototype as parse_proto
except:
    dfttkhome = os.path.abspath(os.path.join('..'))
    sys.path.append(dfttkhome)
    import dfttk.structure_builders.parse_anrl_prototype as parse_proto
    print(dfttkhome)

def test_totalnumber():
    AFLOW_PROTOTYPE_LIBRARY = loadfn(os.path.join(os.path.dirname(parse_proto.__file__),
                                              "aflow_prototype_db.json"))
    assert(len(AFLOW_PROTOTYPE_LIBRARY) == 590)

def test_multi_replace():
    s = "a,b,c,\alpha,\beta,\gamma"
    s = parse_proto.multi_replace(s, {"\a": "a", "\b": "b", "\g": "g"})
    assert(s == "a,b,c,alpha,beta,gamma")

def test_parse_proto_param():
    poscar = parse_proto.poscar_map(67)
    proto_info = parse_proto.parse_proto_param(poscar)
    #[aflow_proto, sg, pearson, strukturbericht, mineral, param, value, ref]
    assert(proto_info[0] == "AB2_aP12_1_4a_8a")
    assert(proto_info[1] == 1)
    assert(proto_info[2] == "aP12")
    assert(proto_info[4] == "anisotropic Pyrite")
    assert(proto_info[5] == "a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12")
    assert(proto_info[6] == "5.417,1.0,1.0,90.0,90.0,90.0,0.001,0.002,0.003,0.4966,0.0001,0.5036,0.5001,0.502,0.0011,-0.0006,0.5013,0.5038,0.3857,0.3832,0.384,0.1149,0.6114,0.8846,0.8854,0.1157,0.6143,0.6153,0.8865,0.1141,0.6151,0.6132,0.6137,0.8854,0.3818,0.1149,0.1147,0.8856,0.3841,0.3857,0.1161,0.8842")
    assert(proto_info[-1] == "Bayliss, Am. Mineral. 62, 1168-72 (1977)")    

def test_gen_proto_dict():
    poscar = parse_proto.poscar_map(252)
    proto_info = parse_proto.parse_proto_param(poscar)
    struct_dict = Structure.from_str(poscar, fmt="POSCAR").as_dict()
    proto_dict = parse_proto.gen_proto_dict(struct_dict, proto_info)
    assert(proto_dict["tags"]["aflow"] == "A3B4_tI28_141_ad_h")

def test_parse_aflow_proto_single():
    url = "http://www.aflowlib.org/CrystalDatabase/POSCAR/A_cF4_225_a.poscar"
    proto_dict = parse_proto.parse_aflow_proto_single(url)
    assert(proto_dict["tags"]["aflow"] == "A_cF4_225_a")
    url_none = "http://www.aflowlib.org/CrystalDatabase/POSCAR/A_cF4_225_a.poscar2"
    proto_dict = parse_proto.parse_aflow_proto_single(url_none)
    assert(proto_dict is None)
    url_poscar = parse_proto.poscar_map(252)
    proto_dict = parse_proto.parse_aflow_proto_single(url_poscar, fmt="poscar")
    assert(proto_dict["tags"]["aflow"] == "A3B4_tI28_141_ad_h")
    
@pytest.mark.skip
def test_parse_aflow_proto_url():
    #This test will cost about 5-10 minutes
    proto_list = parse_proto.parse_aflow_proto_url(write_json=False)
    assert(len(proto_list) == 590)

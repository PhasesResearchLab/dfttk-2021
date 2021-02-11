# coding:utf-8
# The template for single run of DFTTK

################## COMMOM SETTING #############################
TEMPLATE_STRUCTURE_FILENAME = "POSCAR"
#bool, submit(False) or not(True) the workflows to LaunchPad
DRY_RUN = False
#bool, write out(True) or not(False) the workflow (single workflow) to yaml file
WRITE_OUT_WF = False

################ PARAMETERS FOR WF #############################
#str, the absolute path of db.json file, e.g. /storage/home/mjl6505/atomate/config/db.json
#  If None, it will use the configuration in fireworks 
db_file=None
#list, the MAGMOM of the structure, e.g. [4.0, 4.0, -4.0, -4.0]
magmom = None
#int, the number of initial deformations, e.g. 7
num_deformations = 7
#list/tuple(min, max) or float(-max, max), the maximum amplitude of deformation, e.g. (-0.05, 0.1) means (0.95, 1.1) in volume
deformation_fraction = (-0.1, 0.1)
#float, minimum ratio of Volumes spacing, e.g. 0.03
volume_spacing_min = 0.03
#bool, run phonon(True) or not(False)
phonon=False
#list(3x3), the supercell matrix for phonon, e.g. [[2.0, 0, 0], [0, 2.0, 0], [0, 0, 2.0]]
phonon_supercell_matrix=None
#float, the mimimum of temperature in QHA process, e.g. 5
t_min=5
#float, the maximum of temperature in QHA process, e.g. 2000
t_max=2000
#float, the step of temperature in QHA process, e.g. 5
t_step=5
#float, acceptable value for average RMS, recommend >= 0.005
tolerance = 0.01
#str, the vasp command, if None then find in the FWorker configuration
vasp_cmd=None
#dict, metadata to be included, this parameter is useful for filter the data, e.g. metadata={"phase": "BCC_A2", "tag": "AFM"}
metadata=None
#float, the tolerannce for symmetry, e.g. 0.05
symmetry_tolerance = 0.05
#bool, set True to pass initial VASP running if the results exist in DB, use carefully to keep data consistent.
passinitrun=False
#bool, Whether run isif=2 calculation before isif=4 running
run_isif2=False
#bool, Whether pass isif=4 calculation.
pass_isif4=False
#Set the path already exists for new static calculations; if set as '', will try to get the path from db_file
relax_path=''
#dict, dict of class ModifyIncar with keywords in Workflow name. e.g.
"""
modify_incar_params = { 'Full relax': {'incar_update': {"LAECHG":False,"LCHARG":False,"LWAVE":False}},
                        'PreStatic': {'incar_update': {"LAECHG":False,"LCHARG":False,"LWAVE":False}},
                        'PS2': {'incar_update': {"LAECHG":False,"LCHARG":False,"LWAVE":False}}, 
                        'static': {'incar_update': {"LAECHG":False,"LCHARG":False,"LWAVE":False}},
"""
modify_incar_params={}
#dict, dict of class ModifyKpoints with keywords in Workflow name, similar with modify_incar_params
modify_kpoints_params={}
#bool, print(True) or not(False) some informations, used for debug
verbose=False

###################### DO NOT CHANGE THE FOLLOWING LINES ##############################
from pymatgen import MPRester, Structure
from dfttk.wflows import get_wf_gibbs

structure = Structure.from_file(TEMPLATE_STRUCTURE_FILENAME)

if magmom:
    structure.add_site_property('magmom', magmom)

if not db_file:
    from fireworks.fw_config import config_to_dict
    from monty.serialization import loadfn
    db_file = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]

wf = get_wf_gibbs(structure, num_deformations=num_deformations, deformation_fraction=deformation_fraction, 
            phonon=phonon, phonon_supercell_matrix=phonon_supercell_matrix,  t_min=t_min, t_max=t_max, 
            t_step=t_step, tolerance=tolerance, volume_spacing_min=volume_spacing_min,vasp_cmd=vasp_cmd, 
            db_file=db_file, metadata=metadata, name='EV_QHA', symmetry_tolerance=symmetry_tolerance, 
            run_isif2=run_isif2, pass_isif4=pass_isif4, passinitrun=passinitrun, relax_path=relax_path, 
            modify_incar_params=modify_incar_params, modify_kpoints_params=modify_kpoints_params, 
            verbose=verbose)

if not DRY_RUN:
    from fireworks import LaunchPad
    lpad = LaunchPad.auto_load()
    lpad.add_wf(wf)
if WRITE_OUT_WF:
    wf.to_file("dfttk_wf.yaml")
# -*- coding: utf-8 -*-
# The template for batch run of DFTTK
import argparse
from pymatgen.ext.matproj import MPRester, Structure
from pymatgen.io.vasp.inputs import Potcar, Incar
from dfttk.wflows import get_wf_gibbs, get_wf_EV_bjb, get_wf_gibbs_robust, get_wf_borncharge, get_wf_elastic
from dfttk.utils import recursive_glob
from dfttk.structure_builders.parse_anrl_prototype import multi_replace
from dfttk.scripts.querydb import get_eq_structure_by_metadata
from dfttk.scripts.assign_fworker_name import Customizing_Workflows
import dfttk.scripts.querydb as querydb
from fireworks.fw_config import config_to_dict
from monty.serialization import loadfn, dumpfn
import warnings
import copy
import os
import sys
import shutil
import glob


def get_abspath(path):
    """
    Get the absolute path.

    Parameter
        path: str (path-like)
            An absolute or relative path. It can start from home (~)
    Return
        path: str (path-like)
            The absolute path(start from root "/").
    """
    return os.path.abspath(os.path.expanduser(path))

def creat_folders(folder):
    """
    Create folders if not exist, leave a warning if exists
    """
    if os.path.exists(folder):
        warnings.warn( folder + " exists!")
    else:
        os.makedirs(folder)

def get_structure_file(STR_FOLDER=".", RECURSIVE=False, MATCH_PATTERN="*"):
    """
    Get all file names in STR_FOLDER path [and its sub-folder (RECURSIVE)] with condition (MATCH_PATTEN)

    Parameter
        STR_FOLDER: str (path-like)
            The path containing files, Default: .
        RECURSIVE: bool
            Including the sub-folder (True) or not (Default: False)
        MATCH_PATTERN: str
            The match pattern for the file name
    Return
        STR_FILES: list[str (filename-like)]
            The list of all the filenames matching the conditions
    """
    ## Get the file name of structure file
    STR_FOLDER = get_abspath(STR_FOLDER)
    if not os.path.exists(STR_FOLDER):
        raise FileNotFoundError("The specified folder/file '{}' does not exist.".format(STR_FOLDER))
    if os.path.isfile(STR_FOLDER):
        STR_FILES = [STR_FOLDER]
    else:
        if RECURSIVE:
            STR_FILES = recursive_glob(STR_FOLDER, MATCH_PATTERN)
        else:
            STR_FILES = glob.glob(os.path.join(STR_FOLDER,MATCH_PATTERN))
            for file_i in STR_FILES:
                if os.path.isdir(file_i):
                    STR_FILES.remove(file_i)
    return STR_FILES

def get_user_settings(STR_FILENAME, STR_PATH="./", NEW_SETTING="SETTINGS"):
    """
    Get the filename (without ext) of setting file
    (By default: The setting file should be SETTINGS or start with SETTINGS- or end with -SETTINGS (case insensitive))

    Parameter
        STR_FILENAME: str
            The individual tags for the setting file
        STR_PATH: str
            The path of the setting files
        NEW_SETTING: str
            The str to replace "SETTINGS"
    Return
        user_settings: dict
            User settings, if no setting file, return an empty dict
    """
    user_settings = {}
    SETTING_FILENAMES = ["SETTINGS", "settings",
                         "SETTINGS-" + STR_FILENAME, "settings-" + STR_FILENAME,
                         STR_FILENAME + "-SETTINGS", STR_FILENAME + "-settings"]
    replace_dict = {"SETTINGS": NEW_SETTING.upper(), "settings": NEW_SETTING.lower()}
    SETTING_FILENAMES = [multi_replace(item, replace_dict) for item in SETTING_FILENAMES]

    for SETTING_FILENAME in SETTING_FILENAMES:
        for STR_EXT in ['json', 'yaml']:
            SETTING_FULL_FILENAME = os.path.join(STR_PATH, "{}.{}".format(SETTING_FILENAME, STR_EXT))
            if os.path.exists(SETTING_FULL_FILENAME):
                try:
                    user_settings.update(loadfn(SETTING_FULL_FILENAME))
                except Exception as e:
                    raise TypeError("The file contant or file type is not supported. ref. " +\
                        "http://guide.materialsvirtuallab.org/monty/monty.serialization.html#monty.serialization.loadfn")
    return user_settings

def get_wf_single(structure, WORKFLOW="get_wf_gibbs", settings={}):
    """
    Get a single workflow

    Parameters
        structure: pymatgen.Structure
            The structure
        WORKFLOW: str
            The name of the workflow, now only gibbs energy workflow(get_wf_gibbs) is supported
        settings: dict
            User settings for the workflow
    Return
    """
    ################ PARAMETERS FOR WF #############################
    #str, the absolute path of db.json file, e.g. /storage/home/mjl6505/atomate/config/db.json
    #  If None, it will use the configuration in fireworks
    db_file = settings.get('db_file', None)
    #list, the MAGMOM of the structure, e.g. [4.0, 4.0, -4.0, -4.0]
    magmom = settings.get('magmom', None)
    #int, the number of initial deformations, e.g. 7
    num_deformations = settings.get('num_deformations', 7)
    #list/tuple(min, max) or float(-max, max), the maximum amplitude of deformation, e.g. (-0.15, 0.15) means (0.95, 1.1) in volume
    deformation_fraction = settings.get('deformation_fraction', (-0.15, 0.15))
    #float, minimum ratio of Volumes spacing, e.g. 0.05
    volume_spacing_min = settings.get('volume_spacing_min', 0.05)
    #bool, run phonon(True) or not(False)
    phonon = settings.get('phonon', False)
    #list(3x3), the supercell matrix for phonon, e.g. [[2.0, 0, 0], [0, 2.0, 0], [0, 0, 2.0]]
    phonon_supercell_matrix = settings.get('phonon_supercell_matrix', None)
    phonon_supercell_matrix_min = settings.get('phonon_supercell_matrix_min', None)
    phonon_supercell_matrix_max = settings.get('phonon_supercell_matrix_max', None)
    optimize_sc = settings.get('optimize_sc', False)
    #run phonon always, no matter ISIF=4 passed or not
    force_phonon  = settings.get('force_phonon', False)
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
    #str, the vasp command, if None then find in the FWorker configuration
    vasp_cmd = settings.get('vasp_cmd', None)
    #dict, metadata to be included, this parameter is useful for filter the data, e.g. metadata={"phase": "BCC_A2", "tag": "AFM"}
    metadata = settings.get('metadata', None)
    #It is for RobustOptimizeFW, if run ISIF=4 followed ISIF=7
    isif4 = settings.get('isif4', False)
    #The level for robust optimization
    level = settings.get('level', 1)
    #float, the tolerannce for symmetry, e.g. 0.05
    symmetry_tolerance = settings.get('symmetry_tolerance', 0.05)
    #bool, set True to pass initial VASP running if the results exist in DB, use carefully to keep data consistent.
    passinitrun = settings.get('passinitrun', False)
    #bool, Whether run isif=2 calculation before isif=4 running
    run_isif2 = settings.get('run_isif2', False)
    #bool, Whether pass isif=4 calculation.
    pass_isif4 = settings.get('pass_isif4', False)
    #Set the path already exists for new static calculations; if set as '', will try to get the path from db_file
    relax_path = settings.get('relax_path', '')
    #The symmetry tolerance, including three keys,
    #e.g. override_symmetry_tolerances={'tol_strain': 0.05, 'tol_energy': 0.025, 'tol_bond': 0.10}
    override_symmetry_tolerances = settings.get('override_symmetry_tolerances', None)
    #Global settings for all vasp job, e.g.
    #override_default_vasp_params = {'user_incar_settings': {}, 'user_kpoints_settings': {}, 'user_potcar_functional': str}
    #If some value in 'user_incar_settings' is set to None, it will use vasp's default value
    override_default_vasp_params = settings.get('override_default_vasp_params', {})
    #check if fworker_name is assigned
    powerups = settings.get('powerups', {})
    if len(powerups)>0:
        override_default_vasp_params['user_incar_settings'].update({'powerups':powerups})

    #dict, dict of class ModifyIncar with keywords in Workflow name. e.g.
    """
    modify_incar_params = { 'Full relax': {'incar_update': {"LAECHG":False,"LCHARG":False,"LWAVE":False}},
                            'PreStatic': {'incar_update': {"LAECHG":False,"LCHARG":False,"LWAVE":False}},
                            'PS2': {'incar_update': {"LAECHG":False,"LCHARG":False,"LWAVE":False}},
                            'static': {'incar_update': {"LAECHG":False,"LCHARG":False,"LWAVE":False}},
    """
    modify_incar_params = settings.get('modify_incar_params', {})

    #dict, dict of class ModifyKpoints with keywords in Workflow name, similar with modify_incar_params
    modify_kpoints_params = settings.get('modify_kpoints_params', {})
    #bool, print(True) or not(False) some informations, used for debug
    verbose = settings.get('verbose', False)
    #Save the volume data or not ("chgcar", "aeccar0", "aeccar2", "elfcar", "locpot")
    store_volumetric_data = settings.get('store_volumetric_data', False)

    ## The following settings only work for elastic constants workflow
    strain_states = settings.get('strain_states', None)
    stencils = settings.get('stencils', None)
    analysis = settings.get('analysis', True)
    sym_reduce = settings.get('sym_reduce', False)
    order = settings.get('order', 2)
    conventional = settings.get('conventional', False)

    """
    stencils = settings.get('stencils', [0.01])
    #sym_reduce = settings.get('sym_reduce', True)
    stencils = settings.get('stencils', [-0.01,0.01])
    #conventional = settings.get('conventional', True)
    """

    uis = override_default_vasp_params.get('user_incar_settings', {})
    #Set the default value for phonon_supercell_matrix_min/max
    if isinstance(phonon_supercell_matrix, str) and (phonon_supercell_matrix_min is None):
        if phonon_supercell_matrix.lower().startswith('a'):
            phonon_supercell_matrix_min = 60
            phonon_supercell_matrix_max = 130
        elif phonon_supercell_matrix.lower().startswith('l'):
            phonon_supercell_matrix_min = 8
            phonon_supercell_matrix_max = 12
        elif phonon_supercell_matrix.lower().startswith('v'):
            phonon_supercell_matrix_min = 512
            phonon_supercell_matrix_max = 1728
        else:
            raise ValueError("Unknown parameters for phonon_supercell_matrix({}), support 'atoms', 'lattice' or 'volume' or 3x3 list.".format(phonon_supercell_matrix))

    if magmom:
        structure.add_site_property('magmom', magmom)
    elif 'MAGMOM' in uis:
        magmom = uis['MAGMOM']
        if isinstance(magmom, str):
            magmom = Incar.from_string('MAGMOM={}'.format(magmom)).as_dict()['MAGMOM']
        structure.add_site_property('magmom', magmom)
    if not db_file:
        #from fireworks.fw_config import config_to_dict
        #db_file = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]
        db_file = ">>db_file<<"

    """
    if WORKFLOW == "get_wf_gibbs":
        #Currently, only this workflow is supported
        wf = get_wf_gibbs(structure, num_deformations=num_deformations, deformation_fraction=deformation_fraction,
                    phonon=phonon, phonon_supercell_matrix=phonon_supercell_matrix,  t_min=t_min, t_max=t_max,
                    t_step=t_step, eos_tolerance=eos_tolerance, volume_spacing_min=volume_spacing_min,vasp_cmd=vasp_cmd,
                    db_file=db_file, metadata=metadata, name='EV_QHA', symmetry_tolerance=symmetry_tolerance,
                    run_isif2=run_isif2, pass_isif4=pass_isif4, passinitrun=passinitrun, relax_path=relax_path,
                    modify_incar_params=modify_incar_params, modify_kpoints_params=modify_kpoints_params,
                    verbose=verbose, store_volumetric_data=store_volumetric_data)
    elif WORKFLOW == "eos":
    """
    if WORKFLOW == "eos":
        wf = get_wf_EV_bjb(structure, deformation_fraction=deformation_fraction, store_volumetric_data=store_volumetric_data,
                  num_deformations=num_deformations, override_symmetry_tolerances=override_default_vasp_params, metadata=metadata)
    elif WORKFLOW == "robust" or WORKFLOW == "get_wf_gibbs":
        wf = get_wf_gibbs_robust(structure, num_deformations=num_deformations, deformation_fraction=deformation_fraction,
                 phonon=phonon, phonon_supercell_matrix=phonon_supercell_matrix, t_min=t_min, t_max=t_max, t_step=t_step,
                 eos_tolerance=eos_tolerance, volume_spacing_min=volume_spacing_min, vasp_cmd=">>vasp_cmd<<", db_file=">>db_file<<",
                 isif4=isif4, metadata=metadata, name='EV_QHA', override_symmetry_tolerances=override_symmetry_tolerances,
                 override_default_vasp_params=override_default_vasp_params, modify_incar_params=modify_incar_params,
                 modify_kpoints_params=modify_kpoints_params, verbose=verbose, phonon_supercell_matrix_min=phonon_supercell_matrix_min,
                 phonon_supercell_matrix_max=phonon_supercell_matrix_max, optimize_sc=optimize_sc, level=level,
                 force_phonon=force_phonon, stable_tor=stable_tor, store_volumetric_data=store_volumetric_data)
    elif WORKFLOW == "born":
        wf = get_wf_borncharge(structure=structure, metadata=metadata, db_file=">>db_file<<", isif=2, name="born charge",
                      vasp_cmd=">>vasp_cmd<<", override_default_vasp_params=override_default_vasp_params,
                      modify_incar=modify_incar_params)
    elif WORKFLOW == 'elastic':
            wf = get_wf_elastic(structure=structure, metadata=metadata, vasp_cmd=">>vasp_cmd<<", db_file=">>db_file<<", name="elastic",
                       override_default_vasp_params=override_default_vasp_params, strain_states=strain_states,
                       stencils=stencils, analysis=analysis, sym_reduce=sym_reduce, order=order, conventional=conventional)
    else:
        raise ValueError("Currently, only the gibbs energy workflow is supported.")
    return wf


def run(args):
    """
    Run dfttk
    Currently, only support get_wf_gibbs

    Parameters
        STR_FOLDER = args.STRUCTURE_FOLDER
            folder/file containing structures
        MATCH_PATTERN = args.MATCH_PATTERN
            Match patterns for structure file, e.g. *POSCAR
        RECURSIVE = args.RECURSIVE
            recursive or not
        WORKFLOW = args.WORKFLOW
            workflow, current only get_wf_gibbs
        LAUNCH = args.LAUNCH
            Launch to lpad or not
        MAX_JOB = args.MAX_JOB
            Max job to submit
        SETTINGS = args.SETTINGS
            Settings file
        WRITE_OUT_WF = args.WRITE_OUT_WF
            Write out wf file or not
    """
    STR_FOLDER = args.STRUCTURE_FOLDER  # folder/file containing structures
    MATCH_PATTERN = args.MATCH_PATTERN  # Match patterns for structure file, e.g. *POSCAR
    RECURSIVE = args.RECURSIVE          # recursive or not
    WORKFLOW = args.WORKFLOW            # workflow, current only get_wf_gibbs
    PHONON = args.PHONON                # run phonon
    LAUNCH = args.LAUNCH               # Launch to lpad or not
    MAX_JOB = args.MAX_JOB              # Max job to submit
    SETTINGS = args.SETTINGS            # Settings file
    WRITE_OUT_WF = args.WRITE_OUT_WF    # Write out wf file or not
    TAG = args.TAG                      # Metadata from the command line
    APPEND = args.APPEND                # Append calculations, e.g. appending volumes or phonon or born

    if os.path.exists('db.json'):
        db_file = 'db.json'
    else:
        db_file = None

    ## Initial wfs and metadatas
    wfs = []
    metadatas = {}

    if APPEND:
        if TAG:
            metadatas = {os.path.join(os.path.abspath('./'), 'POSCAR'): {'tag':TAG}}
        elif os.path.exists('METADATAS.yaml'):
            metadatas = loadfn('METADATAS.yaml')
        else:
            raise ValueError('For APPEND model, please provide TAG with -tag or provide METADATAS.yaml file')
        for keyi in metadatas:
            (STR_PATH, STR_FILENAME_WITH_EXT) = os.path.split(keyi)
            (STR_FILENAME, STR_EXT) = os.path.splitext(STR_FILENAME_WITH_EXT)
            user_settings = get_user_settings(STR_FILENAME_WITH_EXT, STR_PATH=STR_PATH, NEW_SETTING=SETTINGS)
            metadata = user_settings.get('metadata', {})
            metadata.update(metadatas[keyi])
            user_settings.update({'metadata': metadata})
            structure = get_eq_structure_by_metadata(metadata=metadata, db_file=db_file)
            if structure is None:
                raise FileNotFoundError('There is no static results under current metadata tag({})'.format(metadata['tag']))
            if PHONON:
                user_settings.update({'phonon': True})
            phonon_supercell_matrix = user_settings.get('phonon_supercell_matrix', None)
            if phonon_supercell_matrix is None:
                user_settings.update({"phonon_supercell_matrix": "atoms"})

            wf = get_wf_single(structure, WORKFLOW=WORKFLOW, settings=user_settings)
            wf = Customizing_Workflows(wf)
            if isinstance(wf, list):
                wfs = wfs + wf
            else:
                wfs.append(wf)

            if WRITE_OUT_WF:
                dfttk_wf_filename = os.path.join(STR_PATH, "dfttk_wf-" + STR_FILENAME_WITH_EXT + ".yaml")
                #dumpfn(wf.to_dict(), dfttk_wf_filename)
                dumpfn(wf, dfttk_wf_filename)
    else:
        if os.path.exists('METADATAS.yaml'):
            metadatas = loadfn('METADATAS.yaml')
        ## Get the file names of files
        STR_FILES = get_structure_file(STR_FOLDER=STR_FOLDER, RECURSIVE=RECURSIVE, MATCH_PATTERN=MATCH_PATTERN)
        ## generat the wf
        for STR_FILE in STR_FILES:
            (STR_PATH, STR_FILENAME_WITH_EXT) = os.path.split(STR_FILE)
            (STR_FILENAME, STR_EXT) = os.path.splitext(STR_FILENAME_WITH_EXT)
            str_filename = STR_FILENAME.lower()
            if (str_filename.endswith("-" + SETTINGS.lower()) or
               str_filename.startswith( SETTINGS.lower() + "-") or
               (str_filename == SETTINGS.lower())):
                print(STR_FILE + " is a setting file, not structure file, and skipped when reading the structure.")
            elif STR_FILE == os.path.abspath(__file__):
                #This is current file
                pass
            else:
                flag_run = False
                try:
                    structure = Structure.from_file(STR_FILE)
                    flag_run = True
                except Exception as e:
                    warnings.warn("The name or the contant of " + STR_FILE + " is not supported by dfttk, and skipped. " + \
                        "Ref. https://pymatgen.org/pymatgen.core.structure.html#pymatgen.core.structure.IStructure.from_file")

                if flag_run:
                    user_settings = get_user_settings(STR_FILENAME_WITH_EXT, STR_PATH=STR_PATH, NEW_SETTING=SETTINGS)
                    metadatai = metadatas.get(STR_FILE, None)
                    if metadatai:
                        user_settings.update({'metadata': metadatai})
                    if PHONON:
                        user_settings.update({'phonon': True})
                    phonon_supercell_matrix = user_settings.get('phonon_supercell_matrix', None)
                    if phonon_supercell_matrix is None:
                        user_settings.update({"phonon_supercell_matrix": "atoms"})

                    wf = get_wf_single(structure, WORKFLOW=WORKFLOW, settings=user_settings)
                    wf = Customizing_Workflows(wf)
                    metadatas[STR_FILE] = wf.as_dict()["metadata"]
                    wfs.append(wf)

                    if WRITE_OUT_WF:
                        dfttk_wf_filename = os.path.join(STR_PATH, "dfttk_wf-" + STR_FILENAME_WITH_EXT + ".yaml")
                        dumpfn(wf.to_dict(), dfttk_wf_filename)

        #Write Out the metadata for POST and continue purpose
        dumpfn(metadatas, "METADATAS.yaml")

    """
    _fws = []
    for wflow in wfs:
        revised_wflow = Customizing_Workflows(wflow,user_settings={})
        _fws.append(revised_wflow)
    fws = _fws
    """

    if LAUNCH:
        from fireworks import LaunchPad
        lpad = LaunchPad.auto_load()

        for wflow in wfs:
            lpad.add_wf(wflow)

        if MAX_JOB:
            # Not False or Empty
            if MAX_JOB == 1:
                os.system("qlaunch singleshot")
            else:
                os.system("qlaunch rapidfire -m " + str(MAX_JOB))

def config(args):
    """
    Config dfttk
    It can be used to config atomate and pymatgen
    """
    import dfttk.scripts.config_dfttk as dfttkconfig
    TEST_CONFIG = args.TEST_CONFIG

    ALL = args.ALL
    PATH_TO_STORE_CONFIG = get_abspath(args.PATH_TO_STORE_CONFIG)

    ATOMATE = args.ATOMATE
    VASP_CMD_FLAG = args.VASP_CMD_FLAG
    CONFIG_FOLDER = args.CONFIG_FOLDER
    QUEUE_SCRIPT = args.QUEUE_SCRIPT
    MACHINE = args.MACHINE
    QUEUE_TYPE = args.QUEUE_TYPE
    NODES= args.NODES
    PPN= args.PPN

    PYMATGEN = args.PYMATGEN
    VASP_PSP_DIR = args.VASP_PSP_DIR
    MAPI_KEY = args.MAPI_KEY
    DEFAULT_FUNCTIONAL = args.DEFAULT_FUNCTIONAL
    ACI = args.ACI
    MACHINES = args.MACHINES
    PMEM = args.PMEM

    if ALL:
        ATOMATE = True
        PYMATGEN = True

    TEST_CONFIG_MAP = {"all": [True, True], "atomate": [True, False],
                       "pymatgen": [False, True], "none": [False, False]}
    [test_atomate, test_pymagen] = TEST_CONFIG_MAP[TEST_CONFIG.lower()]
    if test_atomate or test_pymagen:
        # -t parameter exists in the input
        dfttkconfig.test_config(test_pymagen=True, test_atomate=True)
        exit()

    if PATH_TO_STORE_CONFIG is None:
        PATH_TO_STORE_CONFIG = dfttkconfig.default_path()
    PATH_TO_STORE_CONFIG = get_abspath(PATH_TO_STORE_CONFIG)

    if ATOMATE:
        dfttkconfig.config_atomate(path_to_store_config=PATH_TO_STORE_CONFIG,
            config_folder=CONFIG_FOLDER, machine=MACHINE, machines=MACHINES,
            nodes=NODES, ppn=PPN,pmem=PMEM,
            queue_script=QUEUE_SCRIPT, queue_type=QUEUE_TYPE, vasp_cmd_flag=VASP_CMD_FLAG)

    if PYMATGEN:
        dfttkconfig.config_pymatgen(psp_dir=VASP_PSP_DIR, def_fun=DEFAULT_FUNCTIONAL,
            mapi=MAPI_KEY, path_to_store_psp=os.path.join(PATH_TO_STORE_CONFIG, "vasp_psp"), aci=ACI,
            vasp_cmd=VASP_CMD_FLAG, template=QUEUE_SCRIPT, queue_type=QUEUE_TYPE)

def db_remove(args):
    querydb.remove_data_by_metadata(tag=args.TAG, rem_mode=args.MODE, forcedelete=args.FORCE)

def run_dfttk():
    """
    dfttk command
    """
    from dfttk._version import get_versions
    print("DFTTK version: " + get_versions()["version"])
    print("Copyright \u00a9 Phases Research Lab (https://www.phaseslab.com/)\n")

    parser = argparse.ArgumentParser(description='Run DFTTK jobs.')

    subparsers = parser.add_subparsers()

    #SUB-PROCESS: run
    prun = subparsers.add_parser("run", help="Run dfttk.")
    prun.add_argument("-f", "--structure_folder", dest="STRUCTURE_FOLDER", type=str, default=".",
                      help="The folder/file containing the structure,\n"
                           "Default: '.' ")
    prun.add_argument("-mh", "--match_pattern", dest="MATCH_PATTERN", type=str, default="*",
                      help="The match pattern for structure file, and it should be place in quotes."
                           " e.g. '*POSCAR*'. Default: * -- everything except SETTING files, ref. -s")
    prun.add_argument("-s", "--setting", dest="SETTINGS", type=str, default="SETTINGS",
                      help="Specify the name of SETTINGS files (yaml or json file)\n"
                           "Default: SETTINGS (case insensitive and without ext) \n"
                           "The following filename will be treat as SETTINGS file \n"
                           "\t SETTINGS (global settings in the folder)\n"
                           "\t Start with SETTINGS- (individual settings for struct)\n"
                           "\t End with -SETTINGS (individual settings)")
    prun.add_argument("-r", "--recursive", dest="RECURSIVE", action="store_true",
                      help="Recursive the path.")
    prun.add_argument("-wf", "--workflow", dest="WORKFLOW", type=str, default="robust",
                      help="""Specify the workflow to run.\n
                           Default: robust (run get_wf_gibbs_robust workflow) \n
                           (NOTE: currently, only robust and born are supported.)""")
    prun.add_argument("-ph", "--phonon", dest="PHONON", action="store_true",
                      help="Run phonon. This is equivalent with set phonon=True in SETTINGS file")
    prun.add_argument("-tag", "--metatag", dest="TAG", type=str,
                      help="Specify the tag for continue mode")
    prun.add_argument("-a", "--append", dest="APPEND", action="store_true",
                      help="Append calculation according to metadata, e.g. appending volumes or phonon")
    prun.add_argument("-l", "--launch", dest="LAUNCH", action="store_true",
                      help="Launch the wf to launchpad")
    prun.add_argument("-m", "--max_job", dest="MAX_JOB", nargs="?", type=int, default=0,
                      help="Run the job, only works when -l is specified.\n"
                           "Default: 0 (Not submit to queue) \n"
                           "1: qlaunch singleshot (single job) \n"
                           "N(N>1): qlaunch rapidfire -m N")
    prun.add_argument("-o", "--write_out_wf", dest="WRITE_OUT_WF", action="store_true",
                      help="Write out the workflow")
    prun.set_defaults(func=run)

    #SUB-PROCESS: config
    pconfig = subparsers.add_parser("config", help="Config dfttk.")
    pconfig.add_argument("-all", "--all", dest="ALL", action="store_true",default=True,
                         help="Configure atomate and pymatgen. Default: True")
    pconfig.add_argument("-M", "--machine", dest="MACHINE", type=str,
                         default="aci-roar",
                         help="Computer name to be configured.\n"
                              "Default: aci-roar")
    pconfig.add_argument("-MS", "--machines", dest="MACHINES", type=str,
                         default=None,
                         help="User supplied yaml file containing a list of computers with configuration.\n"
                              "Default: None")
    pconfig.add_argument("-p", "--prefix", dest="PATH_TO_STORE_CONFIG", default=".",
                         help="The folder to store the config files.\n"
                              "Default: . (current folder)")
    pconfig.add_argument("-N", "--nodes", dest="NODES", type=int, default=1,
                         help="Number of nodes. Default: 1")
    pconfig.add_argument("-NP", "--ppn", dest="PPN", type=int, default=16,
                         help="Number of cores per node. Default: 16")
    pconfig.add_argument("-PM", "--pmem", dest="PMEM", type=str, default="8gb",
                         help="RAM required per node. Default: 8gb")
    pconfig.add_argument("-a", "--atomate", dest="ATOMATE", action="store_true",
                         help="Configure atomate.")
    pconfig.add_argument("-c", "--config_folder", dest="CONFIG_FOLDER", default=".",
                         help="The folder containing config files, "
                              "at least contain db.json and my_launchpad.yaml. Default: '.'")
    pconfig.add_argument("-q", "--queue_script", dest="QUEUE_SCRIPT", default="vaspjob.pbs",
                         help="The filename of the script for sumitting vasp job. "
                              "It will search in current folder and sub-folders. Default: vaspjob.pbs")
    pconfig.add_argument("-qt", "--queue_type", dest="QUEUE_TYPE", type=str, default="pbs",
                         help="The type of queue system. Default: pbs")
    pconfig.add_argument("-v", "--vasp_cmd_flag", dest="VASP_CMD_FLAG", type=str, default="vasp_std",
                         help="The flag to distinguish vasp_cmd to othe commands in queue_script. Default: vasp_std")
    pconfig.add_argument("-mp", "--pymatgen", dest="PYMATGEN", action="store_true",
                         help="Configure pymatgen.")
    pconfig.add_argument("-aci", "--aci", dest="ACI", action="store_true",
                         help="Using the pesudopotential on the ACI cluster at PSU.")
    pconfig.add_argument("-psp", "--vasp_psp_dir", dest="VASP_PSP_DIR",
                         help="The path of pseudopotentials.")
    pconfig.add_argument("-mapi", "--mapi_key", dest="MAPI_KEY", type=str,
                         help="The API key of Materials Projects")
    pconfig.add_argument("-df", "--default_functional", dest="DEFAULT_FUNCTIONAL", type=str, default="PBE",
                         choices=sorted(Potcar.FUNCTIONAL_CHOICES),
                         help="The default functional. Default: PBE")
    pconfig.add_argument("-t", "--test_config", dest="TEST_CONFIG", nargs="?", const="all", default="None",
                         choices=["all", "pymatgen", "atomate"],
                         help="Test for configurations. Note: currently only support for pymatgen.")
    pconfig.set_defaults(func=config)

    #SUB-PROCESS: db_romove
    pdbrm = subparsers.add_parser("db_remove", help="Remove data in MongoDb.")
    pdbrm.add_argument('-tag', '--tag', dest='TAG', help='Specify the tag. Default: None')
    pdbrm.add_argument('-m', '--mode', dest='MODE', default='vol', help='Specify the remove mode. Default: vol. '
        '1. vol: all volume except dos and bandstructure. 2. allvol: all volume. 3. all: all data.'
        '4. property: all data except volume data. 5. any single properties or volume data, e.g. chgcar, or dos')
    pdbrm.add_argument('-f', '--force', dest='FORCE', action="store_true", help='Force remove (no query). Default: False.')
    pdbrm.set_defaults(func=db_remove)


    # extension by Yi Wang, finalized on August 4, 2020
    # -----------------------------------
    from dfttk.scripts.run_dfttk_ext import run_ext_thelec
    run_ext_thelec(subparsers)


    args = parser.parse_args()

    try:
        a = getattr(args, "func")
    except AttributeError:
        parser.print_help()
        sys.exit(0)
    args.func(args)

if __name__ == '__main__':
    run_dfttk()

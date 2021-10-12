"""
Custom DFTTK Workflows
"""
from pathlib import Path
import os
#unfortunately atomate 0.9.8 explicitly use the HOME environment which is not compatible
#with Windows, so it is added here
os.environ["HOME"] = str(Path.home())

import numpy as np
import copy
from uuid import uuid4
from copy import deepcopy
from fireworks import Workflow, Firework
from atomate.vasp.config import VASP_CMD, DB_FILE
from dfttk.fworks import OptimizeFW, StaticFW, PhononFW, RobustOptimizeFW, BornChargeFW
from dfttk.ftasks import CheckRelaxScheme
from dfttk.input_sets import PreStaticSet, RelaxSet, ForceConstantsSet, ElasticSet
from dfttk.EVcheck_QHA import EVcheck_QHA, Crosscom_EVcheck_QHA, PreEV_check
from dfttk.utils import check_relax_path, add_modify_incar_by_FWname, add_modify_kpoints_by_FWname, supercell_scaling_by_atom_lat_vol
from dfttk.scripts.querydb import is_property_exist_in_db, get_eq_structure_by_metadata
#from atomate.vasp.workflows.base.elastic import get_wf_elastic_constant
from dfttk.elasticity.elastic import get_wf_elastic_constant
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import sys


from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def scale_lattice_vector(structure, factor: float, axisa=False, axisb=False, axisc=True,):
        """
        Performs a scaling of the lattice vectors so that angles are preserved
        axisx(x=a,b,c) = True:
            x lattice scale.
        Args:
            factor (float): scaling factor.
        """
        struct = copy.deepcopy(structure).as_dict()
        matrix = np.array(struct['lattice']['matrix'])
        if axisa:
            matrix[0] *= factor
        if axisb:
            matrix[1] *= factor
        if axisc:
            matrix[2] *= factor
        struct['lattice']['matrix']=list(matrix)
        return Structure.from_dict(struct)

def _get_deformations(def_frac, num_def):
    if isinstance(def_frac, (list, tuple)):
        return np.linspace(1 + def_frac[0], 1 + def_frac[1], num_def)
    else:
        return np.linspace(1 - def_frac, 1 + def_frac, num_def)


def get_wf_EV_bjb(structure, deformation_fraction=(-0.08, 0.12), store_volumetric_data=False,
                  num_deformations=11, override_symmetry_tolerances=None, metadata=None):
    """
    Perform an E-V curve, robustly relaxating all structures on the curve.

    Parameters
    ----------
    structure : pymatgen.Structure

    deformation_fraction : tuple, optional
        Min and max volume scaling factors
    num_deformations : int, optional
        Number of volumes on the E-V curve
    override_symmetry_tolerances : Dict, optional
        See ``tol_X`` below.
    tol_energy : float, optional
        Symmetry is considered broken if the energy decrease for a relaxation
        step exceeds this value.
    tol_strain : float, optional
        Symmetry is considered broken if the the non-isotropic strain for a
        relaxation step exceeds this value.
    tol_bond : float, optional
        Symmetry is considered broken if the bond length per atom change for a
        relaxation step exceeds this value.

    """
    deformations = _get_deformations(deformation_fraction, num_deformations)*structure.volume

    fws = []
    for defo in deformations:
        struct = deepcopy(structure)
        struct.scale_lattice(defo)
        full_relax_fw = RobustOptimizeFW(struct, isif=5, vasp_cmd=VASP_CMD, db_file=DB_FILE,
                                         store_volumetric_data=store_volumetric_data)
        fws.append(full_relax_fw)
    if metadata is not None and all(x in metadata for x in ('phase_name', 'sublattice_configuration')):
        # create a nicer name for the workflow
        subl_config = ':'.join(','.join(subl_comps) for subl_comps in metadata['sublattice_configuration'])
        wfname = f"{metadata['phase_name']}:{subl_config}:{structure.composition.reduced_formula}"
    else:
        wfname = f"unknown:{structure.composition.reduced_formula}:unknown"
    wf = Workflow(fws, name=wfname, metadata=metadata)
    return wf


def get_wf_singleV(structure, store_volumetric_data=False, metadata=None, override_default_vasp_params=None, settings=None):
    """
    Perform single volume relaxation calculation.

    Parameters
    ----------
    structure : pymatgen.Structure
    """
    metadata = metadata or {}
    tag = metadata.get('tag', '{}'.format(str(uuid4())))
    metadata.update({'tag': tag})
    common_kwargs = {"metadata": metadata, "tag":tag,
        'override_default_vasp_params': override_default_vasp_params}
    
    num_deformations = settings.get('num_deformations', 1)
    #list/tuple(min, max) or float(-max, max), the maximum amplitude of deformation, e.g. (-0.15, 0.15) means (0.95, 1.1) in volume
    deformation_fraction = settings.get('deformation_fraction', (-0.0, +0.0))
    deformation_scheme = settings.get('deformation_scheme', 'volume')
    
    isif = settings.get('run_isif', None)
    if not isif:
        isif=3
        if num_deformations>1: 
            if deformation_scheme=='volume': isif = 4
            else: isif = 2

    if deformation_scheme=='volume':
        dmin = pow(1.0+min(deformation_fraction), 1./3.) - 1.0
        dmax = pow(1.0+max(deformation_fraction), 1./3.) - 1.0
        axisa=True
        axisb=True
        axisc=True
    elif deformation_scheme=='a':
        dmin = deformation_fraction
        dmax = deformation_fraction
        axisa=True
        axisb=False
        axisc=False        
    elif deformation_scheme=='b':
        dmin = deformation_fraction
        dmax = deformation_fraction
        axisa=False
        axisb=True
        axisc=False
    elif deformation_scheme=='c':
        dmin = deformation_fraction
        dmax = deformation_fraction
        axisa=False
        axisb=False
        axisc=True
    elif deformation_scheme=='bc' or deformation_scheme=='cb':
        dmin = pow(1.0+min(deformation_fraction), 2./3.) - 1.0
        dmax = pow(1.0+max(deformation_fraction), 2./3.) - 1.0
        axisa=False
        axisb=True
        axisc=True        
    elif deformation_scheme=='ca' or deformation_scheme=='ac':
        dmin = pow(1.0+min(deformation_fraction), 2./3.) - 1.0
        dmax = pow(1.0+max(deformation_fraction), 2./3.) - 1.0
        axisa=True
        axisb=False
        axisc=True
    elif deformation_scheme=='ab' or deformation_scheme=='ba':
        dmin = pow(1.0+min(deformation_fraction), 2./3.) - 1.0
        dmax = pow(1.0+max(deformation_fraction), 2./3.) - 1.0
        axisa=True
        axisb=True
        axisc=False


    deformations = _get_deformations((dmin,dmax), num_deformations)

    fws = []
    for defo in deformations:
        struct = scale_lattice_vector(structure, defo, axisa=axisa, axisb=axisb, axisc=axisc)
        full_relax_fw = OptimizeFW(struct, isif=isif, vasp_cmd=VASP_CMD, db_file=DB_FILE,
            name='Structure_relax_with_ISIF='+str(isif),
            store_volumetric_data=store_volumetric_data, **common_kwargs)
        fws.append(full_relax_fw)
        static_fw = StaticFW(struct, isif=2, vasp_cmd=VASP_CMD, db_file=DB_FILE, 
            name='Staitc',
            vasp_input_set=None, prev_calc_loc=True, parents=full_relax_fw,
            store_volumetric_data=store_volumetric_data, **common_kwargs)
        fws.append(static_fw)    
    if metadata is not None and all(x in metadata for x in ('phase_name', 'sublattice_configuration')):
        # create a nicer name for the workflow
        subl_config = ':'.join(','.join(subl_comps) for subl_comps in metadata['sublattice_configuration'])
        wfname = f"{metadata['phase_name']}:{subl_config}:{structure.composition.reduced_formula}"
    else:
        wfname = f"unknown:{structure.composition.reduced_formula}:unknown"
    wf = Workflow(fws, name=wfname, metadata=metadata)
    return wf


def get_wf_crosscom(structure, metadata=None, settings=None, 
        new_num_deformations=None, new_deformation_fraction=None):
    """
    Perform cross computer QHA calculation without computer dependent.

    Parameters
    ----------
    structure : pymatgen.Structure
    """
    ################ PARAMETERS FOR WF #############################
    #str, the absolute path of db.json file, e.g. /storage/home/mjl6505/atomate/config/db.json
    #  If None, it will use the configuration in fireworks
    db_file = settings.get('db_file', DB_FILE)
    #str, the vasp command, if None then find in the FWorker configuration
    vasp_cmd = settings.get('vasp_cmd', VASP_CMD)
    #dict, metadata to be included, this parameter is useful for filter the data, e.g. metadata={"phase": "BCC_A2", "tag": "AFM"}
    metadata = settings.get('metadata', None)
    tag = metadata.get('tag', '{}'.format(str(uuid4())))
    metadata.update({'tag': tag})

    #int, the number of initial deformations, e.g. 7
    num_deformations = settings.get('num_deformations', 8)
    #list/tuple(min, max) or float(-max, max), the maximum amplitude of deformation, e.g. (-0.15, 0.15) means (0.95, 1.1) in volume
    deformation_fraction = settings.get('deformation_fraction', (-0.15, 0.20))

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

    if new_deformation_fraction is None:
        new_deformation_fraction = copy.deepcopy(deformation_fraction)
        new_num_deformations = num_deformations

    deformation_scheme = settings.get('deformation_scheme', 'volume')
    single_volume = settings.get('single_volume', False)


    isif = settings.get('run_isif', None)
    if not isif:
        isif=3
        if num_deformations>1: 
            if deformation_scheme=='volume': isif = 4
            else: isif = 2

    t_kwargs = {'t_min': t_min, 't_max': t_max, 't_step': t_step}
    common_kwargs = {'vasp_cmd': vasp_cmd, 'db_file': db_file, "metadata": metadata, "tag": tag,
                     'override_default_vasp_params': override_default_vasp_params}
    vasp_kwargs = {'modify_incar_params': modify_incar_params, 'modify_kpoints_params': modify_kpoints_params}
    eos_kwargs = {'deformations': _deformations, 'vol_spacing': vol_spacing, 'eos_tolerance': eos_tolerance, 'threshold': 14}
    a_kwargs = {"structure":structure, "settings":settings, "eos_kwargs":eos_kwargs,
        "static": True, "phonon":phonon, "phonon_supercell_matrix":phonon_supercell_matrix}

    if deformation_scheme=='volume':
        dmin = pow(1.0+min(new_deformation_fraction), 1./3.) - 1.0
        dmax = pow(1.0+max(new_deformation_fraction), 1./3.) - 1.0
        axisa=True
        axisb=True
        axisc=True
    elif deformation_scheme=='a':
        dmin = new_deformation_fraction
        dmax = new_deformation_fraction
        axisa=True
        axisb=False
        axisc=False        
    elif deformation_scheme=='b':
        dmin = new_deformation_fraction
        dmax = new_deformation_fraction
        axisa=False
        axisb=True
        axisc=False
    elif deformation_scheme=='c':
        dmin = new_deformation_fraction
        dmax = new_deformation_fraction
        axisa=False
        axisb=False
        axisc=True
    elif deformation_scheme=='bc' or deformation_scheme=='cb':
        dmin = pow(1.0+min(new_deformation_fraction), 2./3.) - 1.0
        dmax = pow(1.0+max(new_deformation_fraction), 2./3.) - 1.0
        axisa=False
        axisb=True
        axisc=True        
    elif deformation_scheme=='ca' or deformation_scheme=='ac':
        dmin = pow(1.0+min(new_deformation_fraction), 2./3.) - 1.0
        dmax = pow(1.0+max(new_deformation_fraction), 2./3.) - 1.0
        axisa=True
        axisb=False
        axisc=True
    elif deformation_scheme=='ab' or deformation_scheme=='ba':
        dmin = pow(1.0+min(new_deformation_fraction), 2./3.) - 1.0
        dmax = pow(1.0+max(new_deformation_fraction), 2./3.) - 1.0
        axisa=True
        axisb=True
        axisc=False


    deformations = _get_deformations((dmin,dmax), new_num_deformations)


    fws = []
    for defo in deformations:
        struct = scale_lattice_vector(structure, defo, axisa=axisa, axisb=axisb, axisc=axisc)
        full_relax_fw = OptimizeFW(struct, isif=isif, 
            name='Structure_relax_with_ISIF='+str(isif),
            store_volumetric_data=store_volumetric_data,
            t_kwargs=t_kwargs, a_kwargs=a_kwargs, **common_kwargs)
        fws.append(full_relax_fw)
    if not single_volume:
        check_qha_fw = Firework(Crosscom_EVcheck_QHA(verbose=verbose, stable_tor=stable_tor,
            store_volumetric_data=store_volumetric_data, a_kwargs=a_kwargs,
            **eos_kwargs, **vasp_kwargs, **t_kwargs, **common_kwargs),
            parents=fws, name='{}-EVcheck_QHA'.format(structure.composition.reduced_formula))
        fws.append(check_qha_fw)

    wfname = "{}:{}".format(structure.composition.reduced_formula, 'EV_QHA_crosscom')
    wf = Workflow(fws, name=wfname, metadata=metadata)
    return wf


def vol_in_volumes(vol, volumes):
    for v in volumes:
        ncell = int(v/vol+1.e-12)
        if abs(v-ncell*vol)<1.e-12: return True
    return False


def get_wf_elastic(structure=None, metadata=None, tag=None, vasp_cmd=None, db_file=None, name="elastic",
                   vasp_input_set=None, override_default_vasp_params=None, strain_states=None, stencils=None,
                   analysis=True, sym_reduce=False, order=2, conventional=False, **kwargs):
    '''
    Parameter
    ---------
        structure (Structure): the structure to be calculated.
            if the metedata is exist in db, then get the structure from the databae,
                else using current provided structure.
        db_file (str): path to file containing the database settings.
        vasp_input_set (VaspInputSet): vasp input set to be used.  Defaults to ElasticSet in input_sets
        override_default_vasp_params (dict): change the vasp settings.
            e.g. {'user_incar_settings': {}, 'user_kpoints_settings': {}, 'user_potcar_functional': str}
        tag (str):
        #########The following parameters is taken from atomate directly#########
        strain_states (list of Voigt-notation strains): list of ratios of nonzero elements
            of Voigt-notation strain, e. g. [(1, 0, 0, 0, 0, 0), (0, 1, 0, 0, 0, 0), etc.].
        stencils (list of floats, or list of list of floats): values of strain to multiply
            by for each strain state, i. e. stencil for the perturbation along the strain
            state direction, e. g. [-0.01, -0.005, 0.005, 0.01].  If a list of lists,
            stencils must correspond to each strain state provided.
        conventional (bool): flag to convert input structure to conventional structure,
            defaults to False.
        order (int): order of the tensor expansion to be determined.  Defaults to 2 and
            currently supports up to 3.
        analysis (bool): flag to indicate whether analysis task should be added
            and stresses and strains passed to that task
        sym_reduce (bool): Whether or not to apply symmetry reductions

    '''
    vasp_cmd = vasp_cmd or VASP_CMD
    db_file = db_file or DB_FILE

    override_default_vasp_params = override_default_vasp_params or {}

    metadata = metadata or {}
    tag = metadata.get('tag', '{}'.format(str(uuid4())))
    metadata.update({'tag': tag})

    struct_energy_elasticity = is_property_exist_in_db(metadata=metadata, db_file=db_file)
    if struct_energy_elasticity:
        bandgap = struct_energy_elasticity[2]
        static_setting = struct_energy_elasticity[3]

        volumes_existed_calc = is_property_exist_in_db(metadata=metadata, db_file=db_file, collection='elasticity')
        if not volumes_existed_calc:
            structures = struct_energy_elasticity[0]
            bandgap = struct_energy_elasticity[2]
            volumes = struct_energy_elasticity[4]
        else:
            structures = []
            bandgap = []
            volumes = []
            for i,vol in enumerate(struct_energy_elasticity[4]):
                #if vol in volumes_existed_calc:
                if vol_in_volumes(vol, volumes_existed_calc):
                    print("Elasticity already calculated for volume=", vol)
                    continue
                structures.append(struct_energy_elasticity[0][i])
                bandgap.append(struct_energy_elasticity[2][i])
                volumes.append(struct_energy_elasticity[4][i])

        if len(structures)==0:
            print("Elasticity already calculated for all volumes for", \
            struct_energy_elasticity[0][0].composition.reduced_formula, " with tag:", tag, "\n")
            sys.exit()
        else:
            print("Elasticity will be calculated for volumes=", volumes," for", \
            struct_energy_elasticity[0][0].composition.reduced_formula, " with tag:", tag, "\n")
    else: bandgap=False

    wfs = []
    if bandgap:
        if override_default_vasp_params is None: override_default_vasp_params = {}
        using_incar_new_settings = override_default_vasp_params.get('user_incar_settings', None)
        override_default_vasp_params.update(static_setting)
        if using_incar_new_settings is not None:
            override_default_vasp_params['user_incar_settings'].update(using_incar_new_settings)

        for i,struct in enumerate(structures):
            if conventional:
                _struct = SpacegroupAnalyzer(
                    struct).get_conventional_standard_structure()
            else:
                _struct = struct

            """
            if bandgap[i]==0.0:
                override_default_vasp_params['user_incar_settings'].update({'ISMEAR': 0})
                override_default_vasp_params['user_incar_settings'].update({'SIGMA': 0.05})
            print(vasp_input_set.CONFIG['INCAR'])
            """
            vasp_input_set = vasp_input_set or ElasticSet(structure=_struct, **override_default_vasp_params)

            wf_elastic = get_wf_elastic_constant(struct, metadata, strain_states=strain_states, stencils=stencils,
                                db_file=db_file, conventional=conventional, order=order, vasp_input_set=vasp_input_set,
                                analysis=analysis, sym_reduce=sym_reduce, tag='{}-{}'.format(name, tag),
                                vasp_cmd=vasp_cmd, **kwargs)
            wfs.append(wf_elastic)
    else:
        if structure is None:
                raise ValueError('There is no optimized structure with tag={}, Please provide structure.'.format(tag))
        else:
            wfs = get_wf_elastic_constant(structure, metadata, strain_states=strain_states, stencils=stencils,
                                db_file=db_file, conventional=conventional, order=order, vasp_input_set=vasp_input_set,
                                analysis=analysis, sym_reduce=sym_reduce, tag='{}-{}'.format(name, tag),
                                vasp_cmd=vasp_cmd, **kwargs)
    return wfs



def get_wf_borncharge(structure=None, metadata=None, db_file=None, isif=2, name="born charge",
                      vasp_input_set=None,vasp_cmd=None, override_default_vasp_params=None,
                      tag=None, modify_incar=None, **kwargs):
    '''
    The born charge work flow

    structure or metadata must be given.
        If structure is given, then run borncharge for the structure
        If metadata is given, then it will try to find the static calculations form mongodb,
            then it will run born charge calculaions for all structures, if not exist, raise error
        If both are given, then it will try to find structure from mongodb, if not exist, using the given structure
        If both are not given, raise error

    Parameters
    ----------
        structure: pymatgen.Structure
        metadata: dict
            metadata = {'tag': xxxxxxx}
    Return
    ------
        wf: workflow
            The borncharge workflow
    '''
    vasp_cmd = vasp_cmd or VASP_CMD
    metadata = metadata or {}
    tag = metadata.get('tag', '{}'.format(str(uuid4())))
    metadata.update({'tag': tag})
    struct_energy_bandgap = is_property_exist_in_db(metadata=metadata, db_file=db_file)
    if struct_energy_bandgap:
        bandgap = struct_energy_bandgap[2]
        static_setting = struct_energy_bandgap[3]
        if all(np.array(bandgap) == 0):
            print("\nWARNING! No bandgap found for ", \
            struct_energy_bandgap[0][0].composition.reduced_formula, " with tag:", tag, "\n")

        volumes_existed_calc = is_property_exist_in_db(metadata=metadata, db_file=db_file, collection='borncharge')
        if not volumes_existed_calc:
            structures = struct_energy_bandgap[0]
            bandgap = struct_energy_bandgap[2]
            volumes = struct_energy_bandgap[4]
        else:
            structures = []
            bandgap = []
            volumes = []
            for i,vol in enumerate(struct_energy_bandgap[4]):
                if vol in volumes_existed_calc:
                    print("Born effective charges already calculated for volume=", vol)
                    continue
                structures.append(struct_energy_bandgap[0][i])
                bandgap.append(struct_energy_bandgap[2][i])
                volumes.append(struct_energy_bandgap[4][i])

            #Born charge has been calculated
            #raise ValueError('The borncharge with tag={} has been calculated.'.format(tag))
        if len(structures)==0:
            print("Born effective charges already calculated for all volumes for", \
            struct_energy_bandgap[0][0].composition.reduced_formula, " with tag:", tag, "\n")
            sys.exit()
        else:
            print("Born effective charges will be calculated for volumes=", volumes," for", \
            struct_energy_bandgap[0][0].composition.reduced_formula, " with tag:", tag, "\n")
    else: bandgap=False


    fws = []
    if bandgap:
        if override_default_vasp_params is None: override_default_vasp_params = {}
        using_incar_new_settings = None
        if 'user_incar_settings' in override_default_vasp_params.keys():
            using_incar_new_settings = override_default_vasp_params['user_incar_settings']
        override_default_vasp_params.update(static_setting)
        if using_incar_new_settings is not None:
            override_default_vasp_params['user_incar_settings'].update(using_incar_new_settings)

        #any bandgap > 0
        #if any(np.array(bandgap) > 0):
        if True:
            for i in range(0,len(bandgap)):
                structure = structures[i]
                fw = BornChargeFW(structure, isif=isif, name="{}-{:.3f}".format(name, structure.volume),
                                  vasp_cmd=vasp_cmd, metadata=metadata, modify_incar=modify_incar,
                                  override_default_vasp_params=override_default_vasp_params, tag=tag,
                                  prev_calc_loc=False, db_file=db_file, **kwargs)
                fws.append(fw)
    else:
        if structure is None:
            raise ValueError('You must provide metadata existed in mongodb or structure')
        else:
            fw = BornChargeFW(structure, isif=isif, name="{}-{:.3f}".format(name, structure.volume),
                              vasp_cmd=vasp_cmd, metadata=metadata, modify_incar=modify_incar,
                              override_default_vasp_params=override_default_vasp_params, tag=tag,
                              prev_calc_loc=False, db_file=db_file, **kwargs)
            fws.append(fw)
    if not fws:
        raise ValueError('The system is metal or no static result under given metadata in the mongodb')

    wfname = "{}:{}".format(structure.composition.reduced_formula, name)
    wf = Workflow(fws, name=wfname, metadata=metadata)
    return wf


import subprocess
def supercell_scaling_by_Yphon(structure, supercellsize=64):
    structure.to(filename='t.m.p.POSCAR')
    cmd = "Ycell -SN " + str(supercellsize) +" <t.m.p.POSCAR"
    output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        universal_newlines=True)
    with open ("pMatrix", "r") as f : line = f.readline()
    m = np.array([int(x) for x in line.split(' ') if x!='']).reshape((3, 3))
    os.remove("t.m.p.POSCAR")
    os.remove("pMatrix")
    return m


def get_wf_gibbs_robust(structure, num_deformations=7, deformation_fraction=(-0.1, 0.1), phonon=False, isif4=False,
                        phonon_supercell_matrix=None, override_symmetry_tolerances=None, t_min=5, t_max=2000,
                        t_step=5, eos_tolerance=0.01, volume_spacing_min=0.03, vasp_cmd=None, db_file=None,
                        metadata=None, name='EV_QHA', override_default_vasp_params=None, modify_incar_params={},
                        modify_kpoints_params={}, verbose=False, level=1, phonon_supercell_matrix_min=60,
                        phonon_supercell_matrix_max=120, optimize_sc=False, force_phonon=False, stable_tor=0.01,
                        store_volumetric_data=False):
    """
    E - V
    curve

    workflow
    Parameters
    ------
    structure: pymatgen.Structure
        The initial structure
    num_deformations: int
        The number of deformation
    deformation_fraction: float
        Can be a float (a single value) or a 2-type of a min,max deformation fraction.
        Default is (-0.1, 0.1) leading to volumes of 0.90-1.10. A single value gives plus/minus by default.
    phonon : bool
        Whether to do a phonon calculation. Defaults to False, meaning the Debye model.
    phonon_supercell_matrix : list/str
        3x3 array of the supercell matrix, e.g. [[2,0,0],[0,2,0],[0,0,2]]. Must be specified if phonon is specified.
         if string, choose from 'atoms', 'lattice' or 'volume' (only the first letter works, case insensitive),
            which determine who to determin the matrix
    override_symmetry_tolerances : dict
        The symmetry tolerance. It contains three keys, default: {'tol_energy':0.025, 'tol_strain':0.05, 'tol_bond':0.1}
    t_min : float
        Minimum temperature
    t_step : float
        Temperature step size
    t_max : float
        Maximum temperature (inclusive)
    eos_tolerance: float
        Acceptable value for average RMS, recommend >= 0.005.
    volume_spacing_min: float
        Minimum ratio of Volumes spacing
    vasp_cmd : str
        Command to run VASP. If None (the default) is passed, the command will be looked up in the FWorker.
    db_file : str
        Points to the database JSON file. If None (the default) is passed, the path will be looked up in the FWorker.
    name : str
        Name of the workflow
    metadata : dict
        Metadata to include
    override_default_vasp_params: dict
        Override vasp parameters for all vasp jobs. e.g override_default_vasp_params={'user_incar_settings': {'ISIF': 4}}
    modify_incar_params : dict
        Override vasp settings in firework level
        User can use these params to modify the INCAR set. It is a dict of class ModifyIncar with keywords in Workflow name.
    modify_kpoints_params : dict
        User can use these params to modify the KPOINTS set. It is a dict of class ModifyKpoints with keywords in Workflow name.
        Only 'kpts' supported now.
    phonon_supercell_matrix_min/max: int
        minimum/maximum atoms/lattice/volume(controlled by scale_object, default is using atoms)
    optimize_sc: bool
        Optimize the super cell matrix (True) or not (False)
            If False, then use the closest integer transformation matrix of ideal matrix
    stable_tor: float
        The tolerance for phonon stability (the percentage of negative part frequency)
   """
    vasp_cmd = vasp_cmd or VASP_CMD
    db_file = db_file or DB_FILE

    override_symmetry_tolerances = override_symmetry_tolerances or {'tol_energy':0.025, 'tol_strain':0.05, 'tol_bond':0.10}
    override_default_vasp_params = override_default_vasp_params or {}

    site_properties = deepcopy(structure).site_properties

    metadata = metadata or {}
    tag = metadata.get('tag', '{}'.format(str(uuid4())))
    metadata.update({'tag': tag})

    deformations = _get_deformations(deformation_fraction, num_deformations)
    vol_spacing = max((deformations[-1] - deformations[0]) / (num_deformations - 1), volume_spacing_min)

    common_kwargs = {'vasp_cmd': vasp_cmd, 'db_file': ">>db_file<<", "metadata": metadata, "tag": tag,
                     'override_default_vasp_params': override_default_vasp_params}
    robust_opt_kwargs = {'isif': 7, 'isif4': isif4, 'level': level, 'override_symmetry_tolerances': override_symmetry_tolerances}
    vasp_kwargs = {'modify_incar_params': modify_incar_params, 'modify_kpoints_params': modify_kpoints_params}
    t_kwargs = {'t_min': t_min, 't_max': t_max, 't_step': t_step}
    eos_kwargs = {'deformations': deformations, 'vol_spacing': vol_spacing, 'eos_tolerance': eos_tolerance, 'threshold': 14}

    fws = []

    robust_opt_fw = RobustOptimizeFW(structure, prev_calc_loc=False, name='Full relax', store_volumetric_data=store_volumetric_data,
        **robust_opt_kwargs, **vasp_kwargs, **common_kwargs)
    fws.append(robust_opt_fw)
    check_qha_parent = []

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
        ph_scm_size = np.array(phonon_supercell_matrix).shape
        if not (ph_scm_size[0] == 3 and ph_scm_size[1] == 3):
            raise ValueError('Current phonon_supercell_matrix({}) is not correct.'.format(phonon_supercell_matrix))
        phonon_wf = PhononFW(structure, phonon_supercell_matrix, parents=robust_opt_fw, prev_calc_loc='static',
                             name='structure_{:.3f}-phonon'.format(structure.volume), stable_tor=stable_tor,
                             **t_kwargs, **common_kwargs)
        fws.append(phonon_wf)
        check_qha_parent.append(phonon_wf)

    check_relax_fw = Firework(CheckRelaxScheme(db_file=">>db_file<<", tag=tag), parents=robust_opt_fw,
                              name="{}-CheckRelaxScheme".format(structure.composition.reduced_formula))
    fws.append(check_relax_fw)
    check_qha_parent.append(check_relax_fw)
    check_qha_fw = Firework(EVcheck_QHA(site_properties=site_properties,verbose=verbose, stable_tor=stable_tor,
                                        phonon=phonon, phonon_supercell_matrix=phonon_supercell_matrix, force_phonon=force_phonon,
                                        override_symmetry_tolerances=override_symmetry_tolerances, store_volumetric_data=store_volumetric_data,
                                        **eos_kwargs, **vasp_kwargs, **t_kwargs, **common_kwargs),
                            parents=check_qha_parent, name='{}-EVcheck_QHA'.format(structure.composition.reduced_formula))
    fws.append(check_qha_fw)

    wfname = "{}:{}".format(structure.composition.reduced_formula, name)
    wf = Workflow(fws, name=wfname, metadata=metadata)
    add_modify_incar_by_FWname(wf, modify_incar_params = modify_incar_params)
    add_modify_kpoints_by_FWname(wf, modify_kpoints_params = modify_kpoints_params)

    return wf


def get_wf_gibbs(structure, num_deformations=7, deformation_fraction=(-0.1, 0.1), run_isif2=False,
                 phonon=False, phonon_supercell_matrix=None, pass_isif4=False,
                 t_min=5, t_max=2000, t_step=5, tolerance = 0.01, volume_spacing_min = 0.03,
                 vasp_cmd=None, db_file=None, metadata=None, name='EV_QHA', symmetry_tolerance = 0.05,
                 passinitrun=False, relax_path='', modify_incar_params={},
                 modify_kpoints_params={}, verbose=False, store_volumetric_data=False):
    """
    E - V
    curve

    workflow
    Parameters
    ------
    structure: pymatgen.Structure
    num_deformations: int
    deformation_fraction: float
        Can be a float (a single value) or a 2-type of a min,max deformation fraction.
        Default is (-0.05, 0.1) leading to volumes of 0.95-1.10. A single value gives plus/minus
        by default.
    phonon : bool
        Whether to do a phonon calculation. Defaults to False, meaning the Debye model.
    phonon_supercell_matrix : list
        3x3 array of the supercell matrix, e.g. [[2,0,0],[0,2,0],[0,0,2]]. Must be specified if phonon is specified.
    t_min : float
        Minimum temperature
    t_step : float
        Temperature step size
    t_max : float
        Maximum temperature (inclusive)
    tolerance: float
        Acceptable value for average RMS, recommend >= 0.005.
    volume_spacing_min: float
        Minimum ratio of Volumes spacing
    vasp_cmd : str
        Command to run VASP. If None (the default) is passed, the command will be looked up in the FWorker.
    db_file : str
        Points to the database JSON file. If None (the default) is passed, the path will be looked up in the FWorker.
    name : str
        Name of the workflow
    metadata : dict
        Metadata to include
    passinitrun : bool
        Set True to pass initial VASP running if the results exist in DB, use carefully to keep data consistent.
    relax_path : str
        Set the path already exists for new static calculations; if set as '', will try to get the path from db_file.
    modify_incar_params : dict
        User can use these params to modify the INCAR set. It is a dict of class ModifyIncar with keywords in Workflow name.
    modify_kpoints_params : dict
        User can use these params to modify the KPOINTS set. It is a dict of class ModifyKpoints with keywords in Workflow name.
        Only 'kpts' supported now.
    run_isif2: bool
        Whether run isif=2 calculation before isif=4 running.
    pass_isif4: bool
        Whether pass isif=4 calculation.
   """
    vasp_cmd = vasp_cmd or VASP_CMD
    db_file = db_file or DB_FILE

    if db_file == ">>db_file<<":
        #In PengGao's version, some function used the absolute db_file
        from fireworks.fw_config import config_to_dict
        from monty.serialization import loadfn
        db_file = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]

    site_properties = deepcopy(structure).site_properties

    metadata = metadata or {}
    tag = metadata.get('tag', '{}'.format(str(uuid4())))

    if isinstance(deformation_fraction, (list, tuple)):
        deformations = np.linspace(1+deformation_fraction[0], 1+deformation_fraction[1], num_deformations)
        vol_spacing = max((deformation_fraction[1] - deformation_fraction[0]) / (num_deformations - 0.999999) + 0.001,
                          volume_spacing_min)
    else:
        deformations = np.linspace(1-deformation_fraction, 1+deformation_fraction, num_deformations)
        vol_spacing = max(deformation_fraction / (num_deformations - 0.999999) * 2 + 0.001,
                          volume_spacing_min)

    fws = []
    if 'tag' not in metadata.keys():
        metadata['tag'] = tag
    relax_path, run_isif2, pass_isif4 = check_relax_path(relax_path, db_file, tag, run_isif2, pass_isif4)

    if (relax_path == ''):
        # follow a scheme of
        # 1. Full relax + symmetry check
        # 2. If symmetry check fails, detour to 1. Volume relax, 2. inflection detection
        # 3. Inflection detection
        # 4. Static EV
        # 5. Phonon EV
        # for each FW, we set the structure to the original structure to verify to ourselves that the
        # volume deformed structure is set by input set.

        vis_relax = RelaxSet(structure)
        print('Full relax will be running ...')
        full_relax_fw = OptimizeFW(structure, symmetry_tolerance=symmetry_tolerance, job_type='normal', name='Full relax',
                                   prev_calc_loc=False, vasp_input_set=vis_relax, vasp_cmd=vasp_cmd, db_file=db_file,
                                   metadata=metadata, record_path = True, run_isif2=run_isif2, pass_isif4=pass_isif4,
                                   modify_incar_params=modify_incar_params, modify_kpoints_params = modify_kpoints_params,
                                   store_volumetric_data=store_volumetric_data, spec={'_preserve_fworker': True})
        fws.append(full_relax_fw)
    else:
        full_relax_fw = None

    check_result = Firework(EVcheck_QHA(db_file = db_file, tag = tag, relax_path = relax_path, deformations = deformations, site_properties=site_properties,
                                        tolerance = tolerance, threshold = 14, vol_spacing = vol_spacing, vasp_cmd = vasp_cmd,
                                        metadata = metadata, t_min=t_min, t_max=t_max, t_step=t_step, phonon = phonon, symmetry_tolerance = symmetry_tolerance,
                                        phonon_supercell_matrix = phonon_supercell_matrix, verbose = verbose, run_isif2=run_isif2, pass_isif4=pass_isif4,
                                        modify_incar_params=modify_incar_params, modify_kpoints_params = modify_kpoints_params, store_volumetric_data=store_volumetric_data),
                            parents=full_relax_fw, name='%s-EVcheck_QHA' %structure.composition.reduced_formula)
    fws.append(check_result)

    wfname = "{}:{}".format(structure.composition.reduced_formula, name)
    wf = Workflow(fws, name=wfname, metadata=metadata)
    add_modify_incar_by_FWname(wf, modify_incar_params = modify_incar_params)
    add_modify_kpoints_by_FWname(wf, modify_kpoints_params = modify_kpoints_params)

    return wf


def get_wf_gibbs_SQS(structure, num_deformations=7, deformation_fraction=(-0.1, 0.1),
                 phonon=False, phonon_supercell_matrix=None, run_isif2=False, pass_isif4=False,
                 t_min=5, t_max=2000, t_step=5, tolerance = 0.01, volume_spacing_min = 0.03,
                 vasp_cmd=None, db_file=None, metadata=None, name='EV_QHA', symmetry_tolerance = 0.05,
                 passinitrun=False, relax_path='', modify_incar_params={},
                 modify_kpoints_params={}, verbose=False, store_volumetric_data=False):
    """
    E - V
    curve

    workflow
    Parameters
    ------
    structure: pymatgen.Structure
    num_deformations: int
    deformation_fraction: float
        Can be a float (a single value) or a 2-type of a min,max deformation fraction.
        Default is (-0.05, 0.1) leading to volumes of 0.95-1.10. A single value gives plus/minus
        by default.
    phonon : bool
        Whether to do a phonon calculation. Defaults to False, meaning the Debye model.
    phonon_supercell_matrix : list
        3x3 array of the supercell matrix, e.g. [[2,0,0],[0,2,0],[0,0,2]]. Must be specified if phonon is specified.
    t_min : float
        Minimum temperature
    t_step : float
        Temperature step size
    t_max : float
        Maximum temperature (inclusive)
    tolerance: float
        Acceptable value for average RMS, recommend >= 0.005.
    volume_spacing_min: float
        Minimum ratio of Volumes spacing
    vasp_cmd : str
        Command to run VASP. If None (the default) is passed, the command will be looked up in the FWorker.
    db_file : str
        Points to the database JSON file. If None (the default) is passed, the path will be looked up in the FWorker.
    name : str
        Name of the workflow
    metadata : dict
        Metadata to include
    passinitrun : bool
        Set True to pass initial VASP running if the results exist in DB, use carefully to keep data consistent.
    relax_path : str
        Set the path already exists for new static calculations; if set as '', will try to get the path from db_file.
    modify_incar_params : dict
        User can use these params to modify the INCAR set. It is a dict of class ModifyIncar with keywords in Workflow name.
    modify_kpoints_params : dict
        User can use these params to modify the KPOINTS set. It is a dict of class ModifyKpoints with keywords in Workflow name.
        Only 'kpts' supported now.
    run_isif2: bool
        Whether run isif=2 calculation before isif=4 running, pass to get_gibbs.
    pass_isif4: bool
        Whether pass isif=4 calculation, pass to get_gibbs.
    """
    vasp_cmd = vasp_cmd or VASP_CMD
    db_file = db_file or DB_FILE

    metadata = metadata or {}
    tag = metadata.get('tag', '{}'.format(str(uuid4())))

    if isinstance(deformation_fraction, (list, tuple)):
        deformations = np.linspace(1+deformation_fraction[0], 1+deformation_fraction[1], num_deformations)
        vol_spacing = max((deformation_fraction[1] - deformation_fraction[0]) / (num_deformations - 0.999999) + 0.001,
                          volume_spacing_min)
    else:
        deformations = np.linspace(1-deformation_fraction, 1+deformation_fraction, num_deformations)
        vol_spacing = max(deformation_fraction / (num_deformations - 0.999999) * 2 + 0.001,
                          volume_spacing_min)

    if 'tag' not in metadata.keys():
        metadata['tag'] = tag
    relax_path, run_isif2, pass_isif4 = check_relax_path(relax_path, db_file, tag, run_isif2, pass_isif4)

    fws = []
    prestatic_calcs = []
    if relax_path != '':
        pass
    else:
        for i, deformation in enumerate(np.linspace(0.85, 1.15, 7)):
            structure1 = deepcopy(structure)
            structure1.scale_lattice(deformation * structure.volume)
            vis_PreStatic = PreStaticSet(structure1)
            prestatic = StaticFW(structure=structure1, scale_lattice=deformation, name='VR_%.3f-PreStatic' %deformation,
                               prev_calc_loc=False, vasp_input_set=vis_PreStatic, vasp_cmd=vasp_cmd, db_file=db_file,
                               metadata=metadata, Prestatic=True, store_volumetric_data=store_volumetric_data)

            fws.append(prestatic)
            prestatic_calcs.append(prestatic)

    check_result = Firework(PreEV_check(db_file = db_file, tag = tag, relax_path = relax_path, deformations =deformations, structure = structure,
                                        tolerance = tolerance, threshold = 14, vol_spacing = vol_spacing, vasp_cmd = vasp_cmd, run_isif2=run_isif2,
                                        metadata = metadata, t_min=t_min, t_max=t_max, t_step=t_step, phonon = phonon, symmetry_tolerance = symmetry_tolerance,
                                        phonon_supercell_matrix = phonon_supercell_matrix, verbose = verbose, pass_isif4=pass_isif4,
                                        modify_incar_params=modify_incar_params, modify_kpoints_params = modify_kpoints_params),
                            parents=prestatic_calcs, name='%s-PreEV_check' %structure.composition.reduced_formula)
    fws.append(check_result)

    wfname = "{}:{}".format(structure.composition.reduced_formula, name)
    wf = Workflow(fws, name=wfname, metadata=metadata)
    add_modify_incar_by_FWname(wf, modify_incar_params = modify_incar_params)
    add_modify_kpoints_by_FWname(wf, modify_kpoints_params = modify_kpoints_params)

    return wf

# coding: utf-8

import math
import numpy as np
from uuid import uuid4
from copy import deepcopy
from atomate.vasp.database import VaspCalcDb
from atomate.utils.utils import env_chk
from dfttk.utils import sort_x_by_y, mark_adopted, consistent_check_db, check_relax_path
from itertools import combinations
from pymatgen.analysis.eos import EOS
from fireworks import FiretaskBase, LaunchPad, Workflow, Firework
from fireworks.utilities.fw_utilities import explicit_serialize
from dfttk.input_sets import PreStaticSet, RelaxSet, StaticSet, ForceConstantsSet
from dfttk.fworks import OptimizeFW, StaticFW, PhononFW
from dfttk.ftasks import QHAAnalysis
from dfttk.analysis.quasiharmonic import Quasiharmonic

from atomate.vasp.config import VASP_CMD, DB_FILE


def gen_volenergdos(num, volumes, energies, dos_objs=None):
    """
    Extract volume, energies and dos_obj according to the index (num)
    """
    volume = extract_accord_index(index=num, p_in=volumes)
    energy = extract_accord_index(index=num, p_in=energies)
    if dos_objs:
        dos_obj = extract_accord_index(index=num, p_in=dos_objs)
        return volume, energy, dos_obj
    return volume, energy

def extract_accord_index(index, p_in):
    """
    Extract properties according to index

    Parameter
    ---------
        index : list
            The index
        p_in : list
            The property to be extracted
    Return
    ------
        p_out : list
            The extracted value
    """
    p_out = []
    for i in range(len(index)):
        indexi = index[i]
        if indexi > len(p_in) - 1:
            raise IndexError("The index exceed the max index of given properties")
        p_out.append(p_in[indexi])
    return p_out

def check_deformations_in_volumes(deformations, volumes, orig_vol=None):
    """
    Check if the deformations are in volumes, return the deformations which are not in the volumes

    Parameters
    ----------
        deformations: list
            The deformation, thus the value is always near 1, e.g. [0.9, 0.95, 1.05]
        volumes: list
            The volumes have been applied
        orig_vol: float
            The original volume of the structure, can got by Structure.volume
            default: None which means it equals to the averages of the max and min values of volumes
    Return
    ------
        deformations_not_in_vol: ndarray
    """
    if len(deformations) == 0:
        return(np.array([]))
    elif len(volumes) == 0:
        return(np.array(deformations))
    if orig_vol is None:
        orig_vol = (max(volumes) + min(volumes))/2.
    result = []
    min_vol = min(volumes) / orig_vol 
    max_vol = max(volumes) / orig_vol 
    for deformation in deformations:
        if deformation < min_vol or deformation > max_vol:
            result.append(deformation)
    return(np.array(result))

def init_evcheck_result(**kwargs):
    #Transform the input parameters into dict
    """
    Here it used to initial the dict of evcheck_result
    It has following keys:
        append_run_num, correct, volumes, energies, tolerance, 
        threshold, vol_spacing, error, metadata
    """
    return kwargs

def eosfit_stderr(eos_fit, volume, energy):
    """
    Calculate the error of the eos_fit

    Parameter
    ---------
        eos_fit: pymatgen.analysis.eos.fit class
            The eos_fit class, it have the attribute of func
        volume: list
            The volume
        energy: list
            The energy of first-principles calculation, len(energy)=len(volume)
    Return
    ------
        eos_err: float
            The standard error of the volume-energy pair with eos_fit
    """
    fit_value = eos_fit.func(volume)
    eos_err = cal_stderr(value=fit_value, ref=energy)
    return eos_err

def cal_stderr(value, ref=None):
    """
    Calculate the standard error of a list reference to a *ref* list

    Parameter
    ---------
        value: list
            The list to be calculated
        ref: list
            The reference state, len(ref)=len(value), default: None
            Note if ref=None, it equals to [0., 0., 0., ....]
    Return
    ------
        stderr: float
            The standard error
    """
    stderr = 0.
    n = len(value)
    if ref is None:
        ref = [0. for i in range(n)]
    for i in range(n):
        stderr += math.pow((value[i] - ref[i]), 2)
    stderr = stderr / n
    return stderr

def update_err(temperror, error, verbose, ind, **kwargs):
    """
    Update the error if temperror is less than error

    Parameter
    ---------
        temperror: float
            Current error
        error: float
            Previous error
        verbose: bool
            Print(True) the informations or not(False)
        ind: index list
            previous index
        temp_ind: index list  (**kwargs)
            current index, if it is not provided, it equals to ind
    Return
    ------
        error: float
            updated error
        ind: index list
            if temp_ind is not provided, no such return
    """
    if "temp_ind" not in kwargs:
        temp_ind = ind.copy()
    else:
        temp_ind = kwargs["temp_ind"]
    if verbose:
        print('error = %.4f in %s ' %(temperror, temp_ind))
    if temperror < error:
        error = temperror
        if "temp_ind" in kwargs:
            ind = temp_ind
            return error, ind
    return error


@explicit_serialize
class EVcheck_QHA(FiretaskBase):
    ''' 
    If EVcheck(Energies versus Volumes) meets the tolerance, it will launch QHA;
        otherwise it will append more volumes to VASP calculation and take EVcheck again.
    The maximum appending VASP running times set by run_num;
    
    Important Properties:
    correct: whether result satisfies the tolerance
    points: the selected data index
    error: actual fitting error
    eos_fit: eos fitting

    required_params
        None - but you must specify structure or ensure the structure exists in fw_spec(e.g. run the CheckRelaxScheme)
    optional_paramms
        structure
    '''
    _fw_name = 'EVcheck'
    required_params = []
    optional_params = ['structure', 'tag', 'metadata', 'deformations', 'relax_scheme', 'eos_tolerance', 'threshold', 
                       'del_limited', 'vol_spacing', 't_min', 't_max', 't_step', 'phonon', 'phonon_supercell_matrix', 
                       'verbose', 'modify_incar_params', 'run_num','modify_kpoints_params', 'site_properties', 
                       'override_symmetry_tolerances', 'override_default_vasp_params', 'db_file', 'vasp_cmd',
                       'force_phonon', 'stable_tor', 'store_volumetric_data']

    def run_task(self, fw_spec):
        ''' 
        run_num: maximum number of appending VASP running; this limitation is to avoid always running due to bad settings;
            only for internal usage;

        Important args:
        eos_tolerance: acceptable value for average RMS, recommend >= 0.005;
        threshold: total point number above the value should be reduced, recommend < 16 or much time to run;
        del_limited: maximum deletion ration for large results;
        vol_spacing: the maximum ratio step between two volumes, larger step will be inserted points to calculate;
        '''
        # Get the parameters from the object
        max_run = 10
        db_file = env_chk(self.get('db_file', DB_FILE), fw_spec)
        vasp_cmd = env_chk(self.get('vasp_cmd', VASP_CMD), fw_spec)
        deformations = self.get('deformations', [])
        run_num = self.get('run_num', 0)
        eos_tolerance = self.get('eos_tolerance', 0.005)
        threshold = self.get('threshold', 14)
        del_limited = self.get('del_limited', 0.3)
        vol_spacing = self.get('vol_spacing', 0.03)
        t_min = self.get('t_min', 5) 
        t_max = self.get('t_max', 2000)
        t_step = self.get('t_step', 5)
        phonon = self.get('phonon', False)
        force_phonon = self.get('force_phonon', False)
        phonon_supercell_matrix = self.get('phonon_supercell_matrix', None)
        verbose = self.get('verbose', False)
        modify_incar_params = self.get('modify_incar_params', {})
        modify_kpoints_params = self.get('modify_kpoints_params', {})
        site_properties = self.get('site_properties', None)
        override_default_vasp_params = self.get('override_default_vasp_params', {})
        override_symmetry_tolerances = self.get('override_symmetry_tolerances', {})
        store_volumetric_data = self.get('store_volumetric_data', False)

        stable_tor = self.get('stable_tor', 0.01)
        force_phonon = self.get('force_phonon', False)

        relax_structure = self.get('structure') or fw_spec.get('structure', None)
        relax_scheme = self.get('relax_scheme') or fw_spec.get('relax_scheme', [2])
        relax_phonon = fw_spec.get('relax_phonon', False)

        #Only set phonon=True and ISIF=4 passed, then run phonon
        if not force_phonon:
            phonon = phonon and relax_phonon

        metadata = self.get('metadata', {})
        tag = self.get('tag', metadata.get('tag', None))
        if tag is None:
            tag = str(uuid4())
            metadata['tag'] = tag

        common_kwargs = {'vasp_cmd': vasp_cmd, 'db_file': db_file, "metadata": metadata, "tag": tag,
                         'override_default_vasp_params': override_default_vasp_params,}
        vasp_kwargs = {'modify_incar_params': modify_incar_params, 'modify_kpoints_params': modify_kpoints_params}
        t_kwargs = {'t_min': t_min, 't_max': t_max, 't_step': t_step}
        eos_kwargs = {'vol_spacing': vol_spacing, 'eos_tolerance': eos_tolerance, 'threshold': 14}

        run_num += 1
        
        #Some initial checks
        #TODO: add phonon after RobustOptimizeFW
        if phonon:
            #To check if the consistent of phonon and optimize
            if not consistent_check_db(db_file, tag):
                print('Please check DB, DFTTK running ended!')
                return

        if relax_structure is not None:
            structure = deepcopy(relax_structure)
        else:
            raise ValueError('Not structure in spec, please provide structure as input')

        if site_properties:
            for pkey in site_properties:
                structure.add_site_property(pkey, site_properties[pkey])
        # get original EV curve 
        volumes, energies, dos_objs = self.get_orig_EV(db_file, tag)
        vol_adds = check_deformations_in_volumes(deformations, volumes, structure.volume)
        if (len(vol_adds)) == 0:
            self.check_points(db_file, metadata, eos_tolerance, threshold, del_limited, volumes, energies, verbose)
        else:
            self.correct = True
            self.error = 1e10
        
        EVcheck_result = init_evcheck_result(append_run_num=run_num, correct=self.correct, volumes=volumes, 
                         energies=energies, eos_tolerance=eos_tolerance, threshold=threshold, vol_spacing=vol_spacing, 
                         error=self.error, metadata=metadata)

        if self.correct:
            vol_orig = structure.volume
            if (len(vol_adds)) == 0:
                volume, energy, dos_obj = gen_volenergdos(self.points, volumes, energies, dos_objs)
                vol_adds = self.check_vol_coverage(volume, vol_spacing, vol_orig, run_num, 
                                                   energy, structure, dos_obj, phonon, 
                                                   db_file, tag, t_min, t_max, t_step,
                                                   EVcheck_result)   # Normalized to 1
                EVcheck_result['selected'] = volume
                EVcheck_result['append'] = (vol_adds).tolist()
                # Marked as adopted in db
                mark_adopted(tag, db_file, volume, phonon=phonon)
            lpad = LaunchPad.auto_load()
            fws = []
            if len(vol_adds) > 0:      # VASP calculations need to append
                if run_num < max_run:
                    # Do VASP and check again
                    print('Appending the volumes of : %s to calculate in VASP!' %(vol_adds).tolist())
                    calcs = []
                    #vis_relax = RelaxSet(structure)
                    #vis_static = StaticSet(structure)
                    #isif2 = 5 if 'infdet' in relax_path else 4
                    for vol_add in vol_adds:
                        struct = deepcopy(structure)
                        struct.scale_lattice(structure.volume * vol_add)

                        relax_parents_fw = None
                        for isif_i in relax_scheme:
                            #record_path=record_path
                            relax_fw = OptimizeFW(struct, isif=isif_i, store_volumetric_data=store_volumetric_data,
                                 name="relax_Vol{:.3f}".format(vol_add), vasp_input_set=None, job_type="normal",
                                 override_symmetry_tolerances=override_symmetry_tolerances,
                                 prev_calc_loc=True, parents=relax_parents_fw, db_insert=False, force_gamma=True,
                                 modify_incar={}, **vasp_kwargs, **common_kwargs)
                            relax_parents_fw = deepcopy(relax_fw)
                            fws.append(relax_fw)
                            calcs.append(relax_fw)

                        static_fw = StaticFW(struct, isif=relax_scheme[-1], name='static_Vol{:.3f}'.format(vol_add), 
                                        vasp_input_set=None, prev_calc_loc=True, parents=relax_parents_fw,
                                        store_volumetric_data=store_volumetric_data, **common_kwargs)
                        fws.append(static_fw)
                        calcs.append(static_fw)

                        if phonon:
                            #visphonon = ForceConstantsSet(struct)
                            phonon_fw = PhononFW(struct, phonon_supercell_matrix, vasp_input_set=None, stable_tor=stable_tor,
                                                 name='structure_{:.3f}-phonon'.format(vol_add), prev_calc_loc=True,
                                                 parents=static_fw, **t_kwargs, **common_kwargs)
                            fws.append(phonon_fw)
                            calcs.append(phonon_fw)
                    check_result = Firework(EVcheck_QHA(structure=relax_structure, relax_scheme=relax_scheme, store_volumetric_data=store_volumetric_data,
                                                        run_num=run_num, verbose=verbose, site_properties=site_properties, stable_tor=stable_tor,
                                                        phonon=phonon, phonon_supercell_matrix=phonon_supercell_matrix, force_phonon=force_phonon,
                                                        **eos_kwargs, **vasp_kwargs, **t_kwargs, **common_kwargs), 
                                            parents=calcs, name='{}-EVcheck_QHA'.format(structure.composition.reduced_formula))
                    fws.append(check_result)
                    strname = "{}:{}".format(structure.composition.reduced_formula, 'EV_QHA_Append')
                    wfs = Workflow(fws, name = strname, metadata=metadata)
                    if modify_incar_params != {}:
                        from dfttk.utils import add_modify_incar_by_FWname
                        add_modify_incar_by_FWname(wfs, modify_incar_params = modify_incar_params)
                    if modify_kpoints_params != {}:
                        from dfttk.utils import add_modify_kpoints_by_FWname
                        add_modify_kpoints_by_FWname(wfs, modify_kpoints_params = modify_kpoints_params)
                    lpad.add_wf(wfs)
                else:
                    too_many_run_error()
            else:  # No need to do more VASP calculation, QHA could be running 
                print('Success in Volumes-Energies checking, enter QHA ...')
                debye_fw = Firework(QHAAnalysis(phonon=phonon, t_min=t_min, t_max=t_max, t_step=t_step, db_file=db_file, tag=tag, metadata=metadata), 
                                    name="{}-qha_analysis".format(structure.composition.reduced_formula))
                fws.append(debye_fw)
                '''
                # Debye
                debye_fw = Firework(QHAAnalysis(phonon=False, t_min=t_min, t_max=t_max, t_step=t_step, db_file=db_file, tag=tag, metadata=metadata), 
                                    name="{}-qha_analysis-Debye".format(structure.composition.reduced_formula))
                fws.append(debye_fw)
                if phonon:
                    phonon_supercell_matrix = self.get('phonon_supercell_matrix')
                    # do a Debye run before the phonon, so they can be done in stages.
                    phonon_fw = Firework(QHAAnalysis(phonon=True, t_min=t_min, t_max=t_max, t_step=t_step, db_file=db_file, tag=tag, 
                                                     metadata=metadata), parents=debye_fw, name="{}-qha_analysis-phonon".format(structure.composition.reduced_formula))
                    fws.append(phonon_fw)
                '''
                strname = "{}:{}".format(structure.composition.reduced_formula, 'QHA')
                wfs = Workflow(fws, name = strname, metadata=metadata)
                lpad.add_wf(wfs)
        else:   # failure to meet the tolerance
            if len(volumes) == 0: #self.error == 1e10:   # Bad initial running set
                pass_result_error()
            else:                      # fitting fails
                tol_error()
        import json
        with open('EV_check_summary.json', 'w') as fp:
            json.dump(EVcheck_result, fp, indent=4)  
    
    def get_orig_EV(self, db_file, tag):
        vasp_db = VaspCalcDb.from_db_file(db_file, admin = True)
        energies = []
        volumes = []
        dos_objs = []  # pymatgen.electronic_structure.dos.Dos objects
        if vasp_db.collection.count_documents({'$and':[ {'metadata.tag': tag}, {'adopted': True},
                                            {'output.structure.lattice.volume': {'$exists': True}}]}) <= 4:
            static_calculations = vasp_db.collection.find({'$and':[ {'metadata.tag': tag},
                                                        {'output.structure.lattice.volume': {'$exists': True} }]})
        else:
            static_calculations = vasp_db.collection.find({'$and':[ {'metadata.tag': tag}, {'adopted': True},
                                                        {'output.structure.lattice.volume': {'$exists': True}}]})            
        vol_last = 0
        for calc in static_calculations:
            vol = calc['output']['structure']['lattice']['volume']
            if abs((vol - vol_last) / vol) > 1e-8:
                energies.append(calc['output']['energy'])
                volumes.append(vol)
                dos_objs.append(vasp_db.get_dos(calc['task_id']))
            else:
                energies[-1] = calc['output']['energy']
                volumes[-1] = vol
                dos_objs[-1] = vasp_db.get_dos(calc['task_id'])
            vol_last = vol
        energies = sort_x_by_y(energies, volumes)
        dos_objs = sort_x_by_y(dos_objs, volumes)
        ## should sort firstly?
        volumes = sorted(volumes)
        n = len(volumes) - 1           # Delete duplicated
        while n > 0:
            if abs((volumes[n] - volumes[n - 1]) / volumes[n]) < 1e-8:
                volumes.pop(n)
                energies.pop(n)
                dos_objs.pop(n)
            n -= 1
        print('%s Volumes  = %s' %(len(volumes), volumes))
        print('%s Energies = %s' %(len(energies), energies))
        return(volumes, energies, dos_objs)       
        
    def check_points(self, db_file, metadata, eos_tolerance, threshold, del_limited, volumes, energies, verbose = False):
        """
        Check the existing points if they reached the tolerance, 
        if reached, update self.correct as True
        if not, reduce the point to 4
        """
        self.correct = False
        self.error = 1e11
        error = 1e10
        num = np.arange(len(volumes))
        comb = num
        limit = len(volumes) * del_limited
        
        # For len(num) > threshold case, do a whole number fitting to pass numbers delete if met tolerance
        for i in range(1):     # To avoid the codes after except running, ?
            if (len(num) > threshold):
                try:
                    self.check_fit(volumes, energies)
                except:
                    if verbose:
                        print('Fitting error in: ', comb, '. If you can not achieve QHA result, try to run far negative deformations.')
                    break
                temperror = eosfit_stderr(self.eos_fit, volumes, energies)
                error = update_err(temperror=temperror, error=error, verbose=verbose, ind=comb)
        
        # Decrease the quantity of large results
        while (error > eos_tolerance) and (len(num) > limit) and (len(num) > threshold):
            volume, energy = gen_volenergdos(num, volumes, energies)
            try:
                self.check_fit(volume, energy)
            except:
                if verbose:
                    print('Fetal error in Fitting : ', num)         # Seldom
                break
            fit_value = self.eos_fit.func(volume)
            errors = abs(fit_value - energy)
            num = sort_x_by_y(num, errors)
            errors = sorted(errors)
            for m in range(min(len(errors) - threshold, 1)):
                errors.pop(-1)
                num.pop(-1)
            temperror = cal_stderr(errors)
            error, comb = update_err(temperror=temperror, error=error, verbose=verbose, ind=comb, temp_ind=num)

        # combinations
        len_comb = len(comb)
        if (error > eos_tolerance) and (len_comb <= threshold):
            comb_source = comb
            while (error > eos_tolerance) and (len_comb >= 4):
                print('Combinations in "%s"...' %len_comb)
                combination = combinations(comb_source, len_comb)
                for combs in combination:
                    print(volumes)
                    print(energies)
                    volume, energy = gen_volenergdos(combs, volumes, energies)
                    try:
                        self.check_fit(volume, energy)
                    except:
                        if verbose:
                            print('Fitting error in: ', combs, '. If you can not achieve QHA result, try to run far negative deformations.')
                        continue
                    temperror = eosfit_stderr(self.eos_fit, volume, energy)
                    error, comb = update_err(temperror=temperror, error=error, verbose=verbose, ind=comb, temp_ind=combs)
                len_comb -= 1

        print('Minimum error = %s' %error, comb)
        if error <= eos_tolerance:
            self.correct = True
        comb = list(comb)
        comb.sort()
        self.points = comb
        self.error = error        
    
    def check_vol_coverage(self, volume, vol_spacing, vol_orig, run_num, energy, structure, 
        dos_objects, phonon, db_file, tag, t_min, t_max, t_step, EVcheck_result):
        result = []
        volumer = volume.copy()
        
        # Check minimum spacing
        volumer = [vol_i / vol_orig for vol_i in volumer]
        decimals = 4
        for m in range(len(volumer) - 1):
            vol_space_m = volumer[m + 1] - volumer[m]
            if np.around(vol_space_m, decimals) > np.around(vol_spacing, decimals):
                vol = np.linspace(volumer[m], volumer[m + 1], math.ceil(vol_space_m / vol_spacing) + 1).tolist()
                print('Additional volume({}) is appended for the volume space({}) is larger than specified({})'.format(vol[1:-1], vol_space_m, vol_spacing))
                result.append(vol[1:-1])
        
        # To check (and extend) deformation coverage
        # To make sure that coverage extension smaller than interpolation spacing
        vol_spacing = vol_spacing * 0.98   
        
        qha = Quasiharmonic(energy, volume, structure, dos_objects=dos_objects, F_vib=None,
                            t_min=t_min, t_max=t_max, t_step=t_step, poisson=0.363615, bp2gru=1)
        vol_max = np.nanmax(qha.optimum_volumes)
        vol_min = np.nanmin(qha.optimum_volumes)
        EVcheck_result['debye'] = qha.get_summary_dict()
        EVcheck_result['debye']['temperatures'] = EVcheck_result['debye']['temperatures'].tolist()
        if phonon:
            # get the vibrational properties from the FW spec
            # TODO: add a stable check in Quasiharmonic
            vasp_db = VaspCalcDb.from_db_file(db_file, admin=True)
            phonon_calculations = list(vasp_db.db['phonon'].find({'$and':[ {'metadata.tag': tag}, {'adopted': True} ]}))
            vol_vol = [calc['volume'] for calc in phonon_calculations]  # these are just used for sorting and will be thrown away
            vol_f_vib = [calc['F_vib'] for calc in phonon_calculations]
            # sort them order of the unit cell volumes
            vol_f_vib = sort_x_by_y(vol_f_vib, vol_vol)
            f_vib = np.vstack(vol_f_vib)
            qha_phonon = Quasiharmonic(energy, volume, structure, dos_objects=dos_objects, F_vib=f_vib,
                                t_min=t_min, t_max=t_max, t_step=t_step, poisson=0.363615, bp2gru=1)
            vol_max = max(np.nanmax(qha_phonon.optimum_volumes), vol_max)
            vol_min = min(np.nanmax(qha_phonon.optimum_volumes), vol_min)
            EVcheck_result['phonon'] = qha_phonon.get_summary_dict()
            EVcheck_result['phonon']['temperatures'] = EVcheck_result['phonon']['temperatures'].tolist()
        EVcheck_result['MIN_volume_Evaluated'] = '%.3f' %vol_min
        EVcheck_result['MAX_volume_Evaluated'] = '%.3f' %vol_max
        print('Evaluated MIN volume is %.3f;' %vol_min)
        print('Evaluated MAX volume is %.3f;' %vol_max)
        
        vol_max = vol_max / vol_orig
        vol_min = vol_min / vol_orig
        counter = 1
        # Using max_append for reducing unnecessary calculations because of rough fitting
        max_append = 1 if phonon else 2
        # Over coverage ratio set to 1.01 as following
        if volumer[-1] * 1.01 < vol_max:
            print('The current maximum volume is smaller than maximum volume evaluated by quasiharmonic analysis')
            result.append(volumer[-1] + vol_spacing)
            # counter is set to limit calculation times when exception occurs
            while (counter < max_append) and (result[-1] < vol_max):
                result.append(result[-1] + vol_spacing)
                counter += 1
        counter = 1
        if volumer[0] * 0.99 > vol_min:
            print('The current minimum volume is larger than minimum volume evaluated by quasiharmonic analysis')
            result.append(volumer[0] - vol_spacing)
            while (counter < max_append) and (result[-1] > vol_min):
                result.append(result[-1] - vol_spacing)
                counter += 1
        return(np.array(result))
    
    def check_fit(self, volumes, energies):
        eos = EOS('vinet')
        self.eos_fit = eos.fit(volumes, energies)


@explicit_serialize
class PreEV_check(FiretaskBase):
    ''' 
    If EVcheck(Energies versus Volumes) meets the tolerance, it will launch QHA;
        otherwise it will append more volumes to VASP calculation and take EVcheck again.
    The maximum appending VASP running times set by run_num;
    
    Important Properties:
    correct: whether result satisfies the tolerance
    points: the selected data index
    error: actual fitting error
    eos_fit: eos fitting
    
    '''
    _fw_name = 'PreEV_check'
    required_params = ['db_file', 'tag', 'vasp_cmd', 'metadata']
    optional_params = ['deformations', 'relax_path', 'run_num', 'tolerance', 'threshold', 'del_limited', 'vol_spacing', 't_min',
                       't_max', 't_step', 'phonon', 'phonon_supercell_matrix', 'verbose', 'modify_incar_params', 'structure',
                       'modify_kpoints_params', 'symmetry_tolerance', 'run_isif2', 'pass_isif4', 'site_properties',
                       'store_volumetric_data']
    
    def run_task(self, fw_spec):
        ''' 
        run_num: maximum number of appending VASP running; this limitation is to avoid always running due to bad settings;
            only for internal usage;

        Important args:
        tolerance: acceptable value for average RMS, recommend >= 0.005;
        threshold: total point number above the value should be reduced, recommend < 16 or much time to run;
        del_limited: maximum deletion ration for large results;
        vol_spacing: the maximum ratio step between two volumes, larger step will be inserted points to calculate;
        '''
        max_run = 10
        deformations = self.get('deformations') or []
        db_file = self['db_file']
        tag = self['tag']
        vasp_cmd = self['vasp_cmd']
        metadata = self['metadata']
        relax_path = self['relax_path'] or ''
        structure = self.get('structure') or None
        run_num = self.get('run_num') or 0
        tolerance = self.get('tolerance') or 0.005
        threshold = self.get('threshold') or 14
        del_limited = self.get('del_limited') or 0.3
        vol_spacing = self.get('vol_spacing') or 0.03
        t_min = self.get('t_min') or 5 
        t_max = self.get('t_max') or 2000
        t_step = self.get('t_step') or 5
        phonon = self.get('phonon') or False
        phonon_supercell_matrix = self.get('phonon_supercell_matrix') or None
        verbose = self.get('verbose') or False
        modify_incar_params = self.get('modify_incar_params') or {}
        modify_kpoints_params = self.get('modify_kpoints_params') or {}
        symmetry_tolerance = self.get('symmetry_tolerance') or None
        run_isif2 = self.get('run_isif2') or None
        pass_isif4 = self.get('pass_isif4') or False
        site_properties = self.get('site_properties') or None
        store_volumetric_data = self.get('store_volumetric_data', False)
        run_num += 1
        
        volumes, energies = self.get_orig_EV_structure(db_file, tag)
        self.check_points(db_file, metadata, tolerance, 0.1, del_limited, volumes, energies, verbose)
        
        EVcheck_result = init_evcheck_result(run_num, self.correct, volumes, energies, tolerance, 
                                             threshold, vol_spacing, self.error, metadata)

        structure.scale_lattice(self.minE_value)
        if site_properties:
            for pkey in site_properties:
                structure.add_site_property(pkey, site_properties[pkey])
        vol_orig = structure.volume
        volume, energy = gen_volenergdos(self.points, volumes, energies)
        vol_adds = self.check_vol_coverage(volume, vol_spacing, vol_orig, run_num, 
                                           energy, structure, phonon, 
                                           db_file, tag, t_min, t_max, t_step,
                                           EVcheck_result)   # Normalized to 1
        if self.correct or len(vol_adds) > 0:
            EVcheck_result['sellected'] = volume
            EVcheck_result['minE_value'] = self.minE_value
            EVcheck_result['append'] = (vol_adds).tolist()
            # Marked as adopted in db
            lpad = LaunchPad.auto_load()
            fws = []
            if len(vol_adds) > 0:      # VASP calculations need to append
                if run_num < max_run:
                    # Do VASP and check again
                    print('Appending PreStatic of : %s to calculate in VASP!' %(vol_adds * vol_orig).tolist())

                    fws = []
                    prestatic_calcs = []
                    vis_prestatic = PreStaticSet(structure)
                    for vol_add in vol_adds:
                        prestatic = StaticFW(structure=structure, job_type='normal', name='VR_%.3f-PreStatic' %vol_add, 
                                           prev_calc_loc=False, vasp_input_set=vis_prestatic, vasp_cmd=vasp_cmd, db_file=db_file, 
                                           metadata=metadata, Prestatic=True)
                        fws.append(prestatic)
                        prestatic_calcs.append(prestatic)
                
                    check_result = Firework(PreEV_check(db_file = db_file, tag = tag, relax_path = relax_path, deformations = deformations, run_isif2=run_isif2,
                                                        tolerance = tolerance, threshold = 14, vol_spacing = vol_spacing, vasp_cmd = vasp_cmd, pass_isif4=pass_isif4,
                                                        metadata = metadata, t_min=t_min, t_max=t_max, t_step=t_step, phonon = phonon, symmetry_tolerance = symmetry_tolerance,
                                                        phonon_supercell_matrix = phonon_supercell_matrix, verbose = verbose, site_properties=site_properties,
                                                        modify_incar_params=modify_incar_params, modify_kpoints_params = modify_kpoints_params), 
                                            parents=prestatic_calcs, name='%s-PreEV_check%s' %(structure.composition.reduced_formula, run_num))
                    fws.append(check_result)
                    strname = "{}:{}".format(structure.composition.reduced_formula, 'PreEV_check')
                    wfs = Workflow(fws, name = strname, metadata=metadata)
                    if modify_incar_params != {}:
                        from dfttk.utils import add_modify_incar_by_FWname
                        add_modify_incar_by_FWname(wfs, modify_incar_params = modify_incar_params)
                    if modify_kpoints_params != {}:
                        from dfttk.utils import add_modify_kpoints_by_FWname
                        add_modify_kpoints_by_FWname(wfs, modify_kpoints_params = modify_kpoints_params)
                    lpad.add_wf(wfs)
                else:
                    too_many_run_error()
            else:  # No need to do more VASP calculation, QHA could be running 
                relax_path, run_isif2, pass_isif4 = check_relax_path(relax_path, db_file, tag, run_isif2, pass_isif4)
                if relax_path == '':
                    print('Success in PreStatic calculations, entering Position relax ...')
                    vis_relax = RelaxSet(structure)
                    ps2_relax_fw = OptimizeFW(structure, symmetry_tolerance=symmetry_tolerance, job_type='normal', name='MinE V=%.3f relax' %vol_orig, 
                                               prev_calc_loc=False, vasp_input_set=vis_relax, vasp_cmd=vasp_cmd, db_file=db_file, 
                                               metadata=metadata, record_path = True, modify_incar = {'ISIF': 2}, run_isif2=run_isif2, pass_isif4=pass_isif4, 
                                               modify_incar_params=modify_incar_params, modify_kpoints_params = modify_kpoints_params,
                                               spec={'_preserve_fworker': True}, store_volumetric_data=store_volumetric_data)
                    fws.append(ps2_relax_fw)
                else:
                    print('Initial setting found, enter static claculations ...')
                    ps2_relax_fw = None
                check_result = Firework(EVcheck_QHA(db_file = db_file, tag = tag, relax_path = relax_path, tolerance = tolerance, run_isif2=run_isif2,
                                                    threshold = threshold, vol_spacing = vol_spacing, vasp_cmd = vasp_cmd, run_num = run_num,
                                                    metadata = metadata, t_min = t_min, t_max = t_max, t_step = t_step, phonon = phonon, deformations =deformations,
                                                    phonon_supercell_matrix = phonon_supercell_matrix, symmetry_tolerance = symmetry_tolerance,
                                                    modify_incar_params = modify_incar_params, verbose = verbose, pass_isif4=pass_isif4, 
                                                    modify_kpoints_params = modify_kpoints_params, site_properties=site_properties), 
                                        parents = ps2_relax_fw, name='%s-EVcheck_QHA' %structure.composition.reduced_formula, store_volumetric_data=store_volumetric_data)
                fws.append(check_result)
                strname = "{}:{}".format(structure.composition.reduced_formula, 'prePS2_Relax')
                wfs = Workflow(fws, name = strname, metadata=metadata)
                if modify_incar_params != {}:
                    from dfttk.utils import add_modify_incar_by_FWname
                    add_modify_incar_by_FWname(wfs, modify_incar_params = modify_incar_params)
                if modify_kpoints_params != {}:
                    from dfttk.utils import add_modify_kpoints_by_FWname
                    add_modify_kpoints_by_FWname(wfs, modify_kpoints_params = modify_kpoints_params)
                lpad.add_wf(wfs)
        else:   # failure to meet the tolerance
            if len(volumes) == 0: #self.error == 1e10:   # Bad initial running set
                pass_result_error()
            else:                      # fitting fails
                tol_error()

        import json
        with open('PreStatic_check_summary.json', 'w') as fp:
            json.dump(EVcheck_result, fp)  
    
    def get_orig_EV_structure(self, db_file, tag):
        vasp_db = VaspCalcDb.from_db_file(db_file, admin = True)
        energies = []
        volumes = []
        static_calculations = vasp_db.db["PreStatic"].find({'$and':[ {'metadata.tag': tag},
                                                    {'structure.lattice.volume': {'$exists': True} }]})
        vol_last = 0
        for calc in static_calculations:
            vol = calc['structure']['lattice']['volume']
            if abs((vol - vol_last) / vol) > 1e-8:
                energies.append(calc['energy'])
                volumes.append(vol)
            else:
                energies[-1] = calc['energy']
                volumes[-1] = vol
            vol_last = vol
        self.scale_lattice = calc['scale_lattice']
        # Reset to lattice
        energies = sort_x_by_y(energies, volumes)
        volumes = sorted(volumes)
        n = len(volumes) - 1           # Delete duplicated
        while n > 0:
            if abs((volumes[n] - volumes[n - 1]) / volumes[n]) < 1e-8:
                volumes.pop(n)
                energies.pop(n)
            n -= 1
        print('%s Volumes  = %s' %(len(volumes), volumes))
        print('%s Energies = %s' %(len(energies), energies))
        return(volumes, energies) #, structure)       
        
    def check_points(self, db_file, metadata, tolerance, threshold, del_limited, volumes, energies, verbose = False):
        self.correct = False
        self.error = 1e11
        self.minE_value = -1
        error = 1e10
        num = np.arange(len(volumes))
        comb = num
        limit = len(volumes) * del_limited
        
        # For len(num) > threshold case, do a whole number fitting to pass numbers delete if met tolerance
        for i in range(1):     # To avoid the codes after except running
            if (len(num) > threshold):
                try:
                    self.check_fit(volumes, energies)
                except:
                    if verbose:
                        print('Fitting error in: ', comb)
                    break
                temperror = eosfit_stderr(self.eos_fit, volumes, energies)
                error = update_err(temperror=temperror, error=error, verbose=verbose, ind=comb)
        
        # Decrease the quantity of large results
        while (error > tolerance) and (len(num) > limit) and (len(num) > threshold):
            volume, energy = gen_volenergdos(num, volumes, energies)
            try:
                self.check_fit(volume, energy)
            except:
                if verbose:
                    print('Fetal error in Fitting : ', num)         # Seldom
                break
            fit_value = self.eos_fit.func(volume)
            errors = abs(fit_value - energy)
            num = sort_x_by_y(num, errors)
            errors = sorted(errors)
            for m in range(min(len(errors) - threshold, 1)):
                errors.pop(-1)
                num.pop(-1)
            temperror = cal_stderr(errors)
            if verbose:
                print('Absolutely average offest is: %.4f in %s numbers combination.' %(temperror, len(num)))
            if temperror < error:
                error = temperror
                comb = num
                self.minE_value = self.eos_fit.v0

        # combinations

        len_comb = len(comb)
        if (error > tolerance) and (len_comb <= threshold):
            comb_source = comb
            while (error > tolerance) and (len_comb >= 4):
                print('Combinations in "%s"...' %len_comb)
                combination = combinations(comb_source, len_comb)
                for combs in combination:
                    volume, energy = gen_volenergdos(combs, volumes, energies)
                    try:
                        self.check_fit(volume, energy)
                    except:
                        if verbose:
                            print('Fitting error in: ', comb)
                        continue
                    temperror = eosfit_stderr(self.eos_fit, volume, energy)
                    if verbose:
                        print('error = %.4f in %s ' %(temperror, combs))
                    if temperror < error:
                        error = temperror
                        comb = combs
                        self.minE_value = self.eos_fit.v0
                len_comb -= 1

        print('Minimum error = %s' %error, comb)
        if error <= tolerance:
            self.correct = True
            print('The volume ratio of minmum energy in PreStatic is %.3f.' %self.minE_value)
        comb = list(comb)
        comb.sort()
        self.points = comb
        self.error = error
           
    def check_vol_coverage(self, volume, vol_spacing, vol_orig, run_num, energy, structure, 
        phonon, db_file, tag, t_min, t_max, t_step, EVcheck_result):
        result = []
        energier = np.array(energy)
        min_index = np.argmin(energier)
        vol_min = volume[min_index] * 0.93
        vol_max = volume[min_index] * 1.07
        
        counter = 1
        max_append = 10
        # Over coverage ratio set to 1.01 as following
        if volume[-1] < vol_max:
            result.append(volume[-1] + volume[min_index] * 1.05)
            # counter is set to limit calculation times when exception occurs
            while (counter < max_append) and (result[-1] < vol_max):
                result.append(result[-1] + volume[min_index] * 1.05)
                counter += 1
        counter = 1
        if volume[0] > vol_min:
            result.append(volume[0] - volume[min_index] * 1.05)
            while (counter < max_append) and (result[-1] > vol_min):
                result.append(result[-1] - volume[min_index] * 1.05)
                counter += 1
        return(np.array(result)/vol_orig)

    def check_fit(self, volumes, energies):
        eos = EOS('vinet')
        self.eos_fit = eos.fit(volumes, energies)


def tol_error():
    print('''
    #######################################################################
    #                                                                     #
    #           Can not achieve the tolerance requirement, abort!         #
    #                                                                     #
    #######################################################################''')

def pass_result_error():
    print('''
    #######################################################################
    #                                                                     #
    #  "passinitrun = True" could not set while initial results absent.   #
    #                                                                     #
    #######################################################################''')

def too_many_run_error():
    print('''
    #######################################################################
    #                                                                     #
    #            Too many PreStatic VASP running times, abort!            #
    #                          Please check setting!                      #
    #                                                                     #
    #######################################################################''')

def relax_path_error():
    print('''
    #######################################################################
    #                                                                     #
    #       Cannot find relax path for static calculations, exit!         #
    #               You can modify the tag and run again!                 #
    #                                                                     #
    #######################################################################''')
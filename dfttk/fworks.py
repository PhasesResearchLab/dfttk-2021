import warnings
import copy
import os
import numpy as np
from uuid import uuid4
from copy import deepcopy
from fireworks import Firework, PyTask
from atomate.vasp.firetasks.parse_outputs import VaspToDb
from atomate.vasp.firetasks.write_inputs import WriteVaspFromIOSet, ModifyIncar
from atomate.common.firetasks.glue_tasks import PassCalcLocs
from atomate.vasp.firetasks.glue_tasks import CopyVaspOutputs
from atomate.vasp.firetasks.run_calc import RunVaspCustodian
from dfttk.input_sets import RelaxSet, StaticSet, ForceConstantsSet, ATATIDSet, BornChargeSet
from dfttk.ftasks import WriteVaspFromIOSetPrevStructure, SupercellTransformation, CalculatePhononThermalProperties, \
    CheckSymmetry, CheckRelaxation, ScaleVolumeTransformation, TransmuteStructureFile, WriteATATFromIOSet, RunATATCustodian, RunVaspCustodianNoValidate, \
    Record_relax_running_path, Record_PreStatic_result, CheckSymmetryToDb, PhononStable, BornChargeToDb
from atomate import __version__ as atomate_ver
from dfttk import __version__ as dfttk_ver

"""
from dfttk.run_task_ext import nonscalc, InsertXMLToDb
import dfttk.scripts.user_SETTINGS as user_SETTINGS

from monty.serialization import loadfn, dumpfn
if os.path.exists('SETTINGS.yaml'): #treat settings in 'SETTINGS.yaml' as globally accessible
    user_SETTINGS.user_settings=loadfn('SETTINGS.yaml')
"""
import dfttk.run_task_ext.run_task_ext 

STORE_VOLUMETRIC_DATA = ("chgcar", "aeccar0", "aeccar2", "elfcar", "locpot")

class OptimizeFW(Firework):
    """
    Optimize the given structure.

    Results are not entered into the database.

    Args:
        structure (Structure): Input structure. Note that for prev_calc_loc jobs, the structure
            is only used to set the name of the FW and any structure with the same composition
            can be used.
        name (str): Name for the Firework.
        vasp_input_set (VaspInputSet): input set to use (for jobs w/no parents)
            Defaults to MPStaticSet() if None.
        vasp_cmd (str): Command to run vasp.
        prev_calc_loc (bool or str): If true (default), copies outputs from previous calc. If
            a str value, grabs a previous calculation output by name. If False/None, will create
            new static calculation using the provided structure.
        db_file (str): Path to file specifying db credentials.
        parents (Firework): Parents of this particular Firework. FW or list of FWS.
        db_insert : bool
            Whether to insert the task into the database. Defaults to False.
        **kwargs: Other kwargs that are passed to Firework.__init__.
    """
    def __init__(self, structure, scale_lattice=None, isif=4, override_symmetry_tolerances=None, 
                 name="structure optimization", vasp_input_set=None, job_type="normal", vasp_cmd="vasp", 
                 metadata=None, override_default_vasp_params=None, db_file=None, record_path=False, 
                 prev_calc_loc=True, parents=None, db_insert=False, tag=None,
                 run_isif2=False, pass_isif4=False, force_gamma=True, store_volumetric_data=False,
                 modify_incar=None, modify_incar_params={}, modify_kpoints_params={}, **kwargs):
        metadata = metadata or {}
        tag = tag or metadata.get('tag')
        # generate a tag with a warning
        if tag is None:
            tag = str(uuid4())
            metadata['tag'] = tag
        metadata.update({'tag': tag})

        if isinstance(store_volumetric_data, (list, tuple)):
            store_volumetric_data = store_volumetric_data
        elif isinstance(store_volumetric_data, bool):
            if store_volumetric_data:
                store_volumetric_data = STORE_VOLUMETRIC_DATA
            else:
                store_volumetric_data = ()
        else:
            raise ValueError('The store_volumetric_data should be list or bool')

        override_default_vasp_params = override_default_vasp_params or {}
        override_symmetry_tolerances = override_symmetry_tolerances or {}
        vasp_input_set = vasp_input_set or RelaxSet(structure, isif=isif, force_gamma=force_gamma,
                                                       **override_default_vasp_params)
        site_properties = deepcopy(structure).site_properties

        t = []
        # Avoids delivery (prev_calc_loc == '' (instead by True))
        if type(prev_calc_loc) == str:
            t.append(CopyVaspOutputs(calc_dir=prev_calc_loc, contcar_to_poscar=True))
            t.append(WriteVaspFromIOSetPrevStructure(vasp_input_set=vasp_input_set, site_properties=site_properties))
        elif parents:
            if prev_calc_loc:
                t.append(CopyVaspOutputs(calc_loc=prev_calc_loc, contcar_to_poscar=True))
            t.append(WriteVaspFromIOSetPrevStructure(vasp_input_set=vasp_input_set, site_properties=site_properties))
        else:
        #vasp_input_set = vasp_input_set or RelaxSet(structure)  # ??
            t.append(WriteVaspFromIOSetPrevStructure(structure=structure, vasp_input_set=vasp_input_set, site_properties=site_properties))
        if scale_lattice is not None:
            t.append(ScaleVolumeTransformation(scale_factor=scale_lattice, structure=structure))
        t.append(ModifyIncar(incar_update=">>incar_update<<"))
        if modify_incar != None:
             t.append(ModifyIncar(incar_update=modify_incar))
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, job_type=job_type, gzip_output=False))
        t.append(PassCalcLocs(name=name))
        if record_path:
            t.append(Record_relax_running_path(db_file = ">>db_file<<", metadata = metadata, run_isif2=run_isif2, pass_isif4=pass_isif4))
        if db_insert:
            t.append(VaspToDb(db_file=">>db_file<<", additional_fields={"task_label": name, "metadata": metadata}, store_volumetric_data=store_volumetric_data))
        t.append(CheckSymmetryToDb(db_file=">>db_file<<", tag=tag, override_symmetry_tolerances=override_symmetry_tolerances, site_properties=site_properties))
        super(OptimizeFW, self).__init__(t, parents=parents, name="{}-{}".format(structure.composition.reduced_formula, name), **kwargs)


class RobustOptimizeFW(Firework):
    """
    Optimize the given structure.

    Results are not entered into the database.

    Args:
        structure (Structure): Input structure. Note that for prev_calc_loc jobs, the structure
            is only used to set the name of the FW and any structure with the same composition
            can be used.
        name (str): Name for the Firework.
        vasp_input_set (VaspInputSet): input set to use (for jobs w/no parents)
            Defaults to MPStaticSet() if None.
        vasp_cmd (str): Command to run vasp.
        prev_calc_loc (bool or str): If true (default), copies outputs from previous calc. If
            a str value, grabs a previous calculation output by name. If False/None, will create
            new calculation using the provided structure.
        db_file (str): Path to file specifying db credentials.
        parents (Firework): Parents of this particular Firework. FW or list of FWS.
        db_insert : bool
            Whether to insert the task into the database. Defaults to False.
        \*\*kwargs: Other kwargs that are passed to Firework.__init__.
    """
    def __init__(self, structure, isif=7, name="structure optimization", isif4=False, level=1,
                 override_symmetry_tolerances=None, job_type="normal", vasp_input_set=None,
                 vasp_cmd="vasp", metadata=None, override_default_vasp_params=None, db_file=None,
                 prev_calc_loc=True, parents=None, db_insert=False, tag=None, modify_incar_params={},
                 modify_kpoints_params={}, energy_with_isif={}, store_volumetric_data=False, 
                 **kwargs):
        metadata = metadata or {}
        tag = tag or metadata.get('tag')
        # generate a tag with a warning
        if tag is None:
            tag = str(uuid4())
            metadata['tag'] = tag

        if isinstance(store_volumetric_data, (list, tuple)):
            store_volumetric_data = store_volumetric_data
        elif isinstance(store_volumetric_data, bool):
            if store_volumetric_data:
                store_volumetric_data = STORE_VOLUMETRIC_DATA
            else:
                store_volumetric_data = ()
        else:
            raise ValueError('The store_volumetric_data should be list or bool')

        override_default_vasp_params = override_default_vasp_params or {}
        override_symmetry_tolerances = override_symmetry_tolerances or {}
        vasp_input_set = vasp_input_set or RelaxSet(structure, isif=isif, **override_default_vasp_params)
        site_properties = deepcopy(structure).site_properties

        t = []
        if parents:
            if prev_calc_loc:
                t.append(CopyVaspOutputs(calc_loc=prev_calc_loc, contcar_to_poscar=True))
            t.append(WriteVaspFromIOSetPrevStructure(vasp_input_set=vasp_input_set, site_properties=site_properties))
        else:
            t.append(WriteVaspFromIOSetPrevStructure(structure=structure, vasp_input_set=vasp_input_set, site_properties=site_properties))
        t.append(ModifyIncar(incar_update=">>incar_update<<"))
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, job_type=job_type, gzip_output=False))
        t.append(PassCalcLocs(name=name))
        if db_insert:
            t.append(VaspToDb(db_file=">>db_file<<", additional_fields={"task_label": name, "metadata": metadata}, store_volumetric_data=store_volumetric_data))
        t.append(CheckSymmetryToDb(db_file=">>db_file<<", tag=tag, site_properties=site_properties))

        common_kwargs = {'vasp_cmd': vasp_cmd, 'db_file': ">>db_file<<", "metadata": metadata, "tag": tag,
                         'override_default_vasp_params': override_default_vasp_params}
        static_kwargs = {}
        relax_kwargs = {}
        t.append(CheckRelaxation(db_file=">>db_file<<", metadata=metadata, tag=tag, isif4=isif4, level=level, energy_with_isif=energy_with_isif,
                                 common_kwargs=common_kwargs, relax_kwargs=relax_kwargs, static_kwargs=static_kwargs, site_properties=site_properties,
                                 store_volumetric_data=store_volumetric_data, 
                                 **override_symmetry_tolerances))
        super().__init__(t, parents=parents, name="{}-{}".format(structure.composition.reduced_formula, name), **kwargs)


class StaticFW(Firework):
    """
    Standard static calculation Firework - either from a previous location or from a structure.

    Parameters
    ----------
    structure : pymatgen.Structure
        Input structure. Note that for prev_calc_loc jobs, the structure
        is only used to set the name of the FW and any structure with the same composition
        can be used.
    name : str
        Name for the Firework.
    vasp_input_set : pymategen.io.vasp.inputs.VaspInputSet
        Input set to use. Defaults to StaticSet() if None.
    vasp_cmd : str
        Command to run vasp.
    prev_calc_loc : (bool or str)
        If true (default), copies outputs from previous calc. If
        a str value, grabs a previous calculation output by name. If False/None, will create
        new static calculation using the provided structure.
    db_file : str
        Path to file specifying db credentials.
    parents : Firework
        Parents of this particular Firework. FW or list of FWS.
    **kwargs : dict
        Other kwargs that are passed to Firework.__init__.
    """
    def __init__(self, structure, isif=2, scale_lattice=None, name="static", vasp_input_set=None, 
                 vasp_cmd="vasp", metadata=None, prev_calc_loc=True, Prestatic=False, modify_incar=None, 
                 db_file=None, parents=None, tag=None, override_default_vasp_params=None,
                 store_volumetric_data=False, **kwargs):

        # TODO: @computron - I really don't like how you need to set the structure even for
        # prev_calc_loc jobs. Sometimes it makes appending new FWs to an existing workflow
        # difficult. Maybe think about how to remove this need? -computron
        metadata = metadata or {}
        tag = tag or metadata.get('tag')
        # generate a tag with a warning
        if tag is None:
            tag = str(uuid4())
            metadata['tag'] = tag

        if isinstance(store_volumetric_data, (list, tuple)):
            store_volumetric_data = store_volumetric_data
        elif isinstance(store_volumetric_data, bool):
            if store_volumetric_data:
                store_volumetric_data = STORE_VOLUMETRIC_DATA
            else:
                store_volumetric_data = ()
        else:
            raise ValueError('The store_volumetric_data should be list or bool')

        override_default_vasp_params = override_default_vasp_params or {}
        vasp_input_set = vasp_input_set or StaticSet(structure, isif=isif, **override_default_vasp_params)
        site_properties = deepcopy(structure).site_properties
        # Avoids delivery (prev_calc_loc == '' (instead by True))
        t = []
        if type(prev_calc_loc) == str:
            t.append(CopyVaspOutputs(calc_dir=prev_calc_loc, contcar_to_poscar=True))
            t.append(WriteVaspFromIOSetPrevStructure(vasp_input_set=vasp_input_set, site_properties=site_properties))
        elif parents:
            if prev_calc_loc:
                t.append(CopyVaspOutputs(calc_loc=prev_calc_loc, contcar_to_poscar=True))
            t.append(WriteVaspFromIOSetPrevStructure(vasp_input_set=vasp_input_set, site_properties=site_properties))
        else:
            t.append(WriteVaspFromIOSetPrevStructure(structure=structure, vasp_input_set=vasp_input_set, site_properties=site_properties))
        if (scale_lattice is not None) and not Prestatic:
            t.append(ScaleVolumeTransformation(scale_factor=scale_lattice, structure=structure))
        t.append(ModifyIncar(incar_update=">>incar_update<<"))
        if modify_incar != None:
             t.append(ModifyIncar(incar_update=modify_incar))
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<", gzip_output=False))
        t.append(PassCalcLocs(name=name))
        if Prestatic:
            t.append(Record_PreStatic_result(db_file = ">>db_file<<", metadata = metadata, structure = structure, scale_lattice = scale_lattice))
        else:
            t.append(VaspToDb(db_file=">>db_file<<", parse_dos=True, additional_fields={"task_label": name, "metadata": metadata,
                                "version_atomate": atomate_ver, "version_dfttk": dfttk_ver, "adopted": True, "tag": tag},
                                store_volumetric_data=store_volumetric_data))
            run_task_ext(t,vasp_cmd,">>db_file<<",structure,tag)
            """
            print (user_SETTINGS.user_settings.get('store_raw_vasprunxml', False))

            if user_SETTINGS.user_settings.get('store_raw_vasprunxml', False):
                t.append(nonscalc())
                t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<", gzip_output=False))
                t.append(InsertXMLToDb(db_file=">>db_file<<", structure=structure, 
                    tag=tag, xml="vasprun.xml"))
            """

        t.append(CheckSymmetryToDb(db_file=">>db_file<<", tag=tag, site_properties=site_properties))
        super(StaticFW, self).__init__(t, parents=parents, name="{}-{}".format(
            structure.composition.reduced_formula, name), **kwargs)


class InflectionDetectionFW(Firework):
    """
    Inflection detection Firework. Assumes that there are items in the calc_loc
     of the Firework called 'Full relax' and 'Volume relax'.

    Parameters
    ----------
    structure : pymatgen.Structure
        Input structure. Note that for prev_calc_loc jobs, the structure
        is only used to set the name of the FW and any structure with the same composition
        can be used.
    name : str
        Name for the Firework.
    vasp_input_set : pymategen.io.vasp.inputs.VaspInputSet
        Input set to use. Defaults to StaticSet() if None.
    vasp_cmd : str
        Command to run vasp.
    db_file : str
        Path to file specifying db credentials.
    parents : Firework
        Parents of this particular Firework. FW or list of FWS.
    continuation : bool
        Whether this Firework is continuing from a previous run of the InflectionDetection code.
    **kwargs : dict
        Other kwargs that are passed to Firework.__init__.

    Notes
    -----
    1. Copy vasp outputs from volume and full relax, write to str_beg.out and str_end.out, respectively
    2. Write the vaspid.wrap file
    3. Run ATAT's robustrelax_vasp, detouring a continuation of this Firework if we hit the walltime
    4. Move str_relax.out to CONTCAR

    """
    def __init__(self, structure, name="infdet", input_set=None, metadata=None, prev_calc_loc=True,
                 db_file=None, parents=None, continuation=False, run_isif2=False, pass_isif4=False, **kwargs):
        metadata = metadata or {}
        input_set = input_set or ATATIDSet(structure)

        t = []

        if not continuation:
            # Copy the volume relax CONTCAR to POSCAR and the full relax CONTCAR as CONTCAR. Get the CHGCAR and WAVECAR from the fully relaxed structure
            # There are other ways to do this, but it's important to pay attention
            # to the order so that work is not destoryed because CopyVaspOutputs
            # will always give back a POSCAR (or CONTCAR as POSCAR), KPOINTS, INCAR, POTCAR, OUTCAR, and
            # vasprun.xml.
            # What we do here ensures that
            # 1. We get the the WAVECAR and CHGCAR from the full relax
            # 2. We do not overwrite the structure that we took from the full relax when we copy the volume relax
            t.append(CopyVaspOutputs(calc_loc='Full relax', contcar_to_poscar=False, additional_files=["CONTCAR"]))
            t.append(CopyVaspOutputs(calc_loc='Volume relax', contcar_to_poscar=True))
            # Move the volume relaxed POSCAR to str_beg.out
            t.append(TransmuteStructureFile(input_fname='POSCAR', output_fname='str_beg.out'))
            # Move the fully relaxed CONTCAR to str_end.out
            t.append(TransmuteStructureFile(input_fname='CONTCAR', output_fname='str_end.out'))
            # write the vaspid.wrap file
            t.append(WriteATATFromIOSet(input_set=input_set))
        else:
            # Copy all the files from the previous run.
            files_needed = ['CONTCAR', 'str_beg.out', 'str_end.out', 'str_relax.out', 'epipos.out', 'epidir.out',
                            'epipos.out', 'infdet.log', 'str_current.out']
            t.append(CopyVaspOutputs(calc_loc=prev_calc_loc, contcar_to_poscar=False, additional_files=files_needed))

        # Unfortunately, it seems that PassCalcLocs must happen before
        # running ATAT because it can return a FW action that is dynamic and will
        # skip the remaining Firetasks. We don't really want to do this (I think)
        # because if this fizzles, the calc_locs will still be changed if this is rerun.
        t.append(PassCalcLocs(name=name))
        # Run ATAT's inflection detection
        t.append(RunATATCustodian(continuation=continuation, name=name))
        t.append(Record_relax_running_path(db_file = ">>db_file<<", metadata = metadata, run_isif2=run_isif2, pass_isif4=pass_isif4))
        super(InflectionDetectionFW, self).__init__(t, parents=parents,
                                                    name="{}-{}".format(structure.composition.reduced_formula, name), **kwargs)


class PhononFW(Firework):
    """
    Calculation of phonon thermal properties by direct calculation of force constants.

    Parameters
    ----------
    structure : pymatgen.Structure
        Input structure. Note that for prev_calc_loc jobs, the structure
        is only used to set the name of the FW and any structure with the same composition
        can be used.
    supercell_matrix:
        3x3 array of the supercell matrix, e.g. [[2,0,0],[0,2,0],[0,0,2]].
    name : str
        Name for the Firework.
    vasp_input_set : pymategen.io.vasp.inputs.VaspInputSet
        Input set to use. Defaults to ForceConstantsSet() if None.
    vasp_cmd : str
        Command to run vasp.
    prev_calc_loc : (bool or str)
        If true (default), copies outputs from previous calc. If
        a str value, grabs a previous calculation output by name. If False/None, will create
        new static calculation using the provided structure.
    db_file : str
        Path to file specifying db credentials.
    parents : Firework
        Parents of this particular Firework. FW or list of FWS.
    **kwargs : dict
        Other kwargs that are passed to Firework.__init__.
    """
    def __init__(self, structure, supercell_matrix, t_min=5, t_max=2000, t_step=5,
                 name="phonon", vasp_input_set=None, override_default_vasp_params=None,
                 vasp_cmd="vasp", metadata=None, tag=None, qpoint_mesh=(50, 50, 50),
                 prev_calc_loc=True, db_file=None, parents=None, stable_tor=0.01,
                 **kwargs):

        metadata = metadata or {}
        tag = tag or metadata.get('tag')
        # generate a tag with a warning
        if tag is None:
            tag = str(uuid4())
            warnings.warn('No ``tag`` was passed explicitly or in ``metadata`` to PhononFW. In order to find this Firework later, you should assign one. This was assigned: {}'.format(tag))
            metadata['tag'] = tag

        override_default_vasp_params = override_default_vasp_params or {}
        ncell = int(0.5+np.linalg.det(supercell_matrix))
        tmp = copy.deepcopy(override_default_vasp_params)
        if 'user_incar_settings' in tmp:
          if 'magmom' in tmp['user_incar_settings']:
            mag = tmp['user_incar_settings']['magmom']
            supermag = []
            for site in mag:
                n = str(site).split('*')
                if len(n)==1:
                    supermag.append('{}*{}'.format(ncell,float(n[0])))
                else:
                    supermag.append('{}*{}'.format(ncell*int(n[0]),float(n[1])))
            tmp['user_incar_settings']['magmom']=supermag
            print("phonon setting", tmp)

        vasp_input_set = vasp_input_set or ForceConstantsSet(structure, **tmp)

        supercell_structure = deepcopy(structure)
        supercell_structure.make_supercell(supercell_matrix)
        supercell_site_properties = deepcopy(supercell_structure.site_properties)

        t = []

        # We need to get the POSCAR from the previous run or from the passed Structure
        # so it can be transformed to a supercell in the next step
        if parents:
            if prev_calc_loc:
                t.append(CopyVaspOutputs(calc_loc=prev_calc_loc, contcar_to_poscar=True))
        else:
            # write the input set first, just to get the POSCAR file in the directory
            # the other inputs will get overridden by WriteVaspFromIOSetPrevStructure
            t.append(WriteVaspFromIOSetPrevStructure(structure=structure, vasp_input_set=vasp_input_set, site_properties=site_properties))

        t.append(SupercellTransformation(supercell_matrix=supercell_matrix))
        t.append(WriteVaspFromIOSetPrevStructure(vasp_input_set=vasp_input_set, site_properties=supercell_site_properties))
        t.append(RunVaspCustodianNoValidate(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<", gzip_output=False))
        # we skipped the validation so we can potentially fix the vasprun.xml file.
        # Fix and validate here.
        t.append(PyTask(func='dfttk.vasprun_fix.fix_vasprun', args=['vasprun.xml']))
        t.append(PassCalcLocs(name=name))
        t.append(CalculatePhononThermalProperties(supercell_matrix=supercell_matrix, t_min=t_min, t_max=t_max, t_step=t_step, db_file=">>db_file<<", tag=tag, metadata=metadata))
        t.append(PhononStable(supercell_matrix=supercell_matrix, db_file=">>db_file<<", tag=tag, metadata=metadata, qpoint_mesh=qpoint_mesh, stable_tor=stable_tor))

        super(PhononFW, self).__init__(t, parents=parents, name="{}-{}".format(
            structure.composition.reduced_formula, name), **kwargs)


class BornChargeFW(Firework):
    """
    Standard static calculation Firework - either from a previous location or from a structure.

    Parameters
    ----------
    structure : pymatgen.Structure
        Input structure. Note that for prev_calc_loc jobs, the structure
        is only used to set the name of the FW and any structure with the same composition
        can be used.
    name : str
        Name for the Firework.
    vasp_input_set : pymategen.io.vasp.inputs.VaspInputSet
        Input set to use. Defaults to StaticSet() if None.
    vasp_cmd : str
        Command to run vasp.
    prev_calc_loc : (bool or str)
        If true (default), copies outputs from previous calc. If
        a str value, grabs a previous calculation output by name. If False/None, will create
        new static calculation using the provided structure.
    db_file : str
        Path to file specifying db credentials.
    parents : Firework
        Parents of this particular Firework. FW or list of FWS.
    **kwargs : dict
        Other kwargs that are passed to Firework.__init__.
    """
    def __init__(self, structure, isif=2, scale_lattice=None, name="born charge", vasp_input_set=None,
                 vasp_cmd="vasp", metadata=None, override_default_vasp_params=None, tag=None,
                 prev_calc_loc=True, modify_incar=None, db_file=None, parents=None, **kwargs):

        metadata = metadata or {}
        tag = tag or metadata.get('tag')
        # generate a tag with a warning
        if tag is None:
            tag = str(uuid4())
            metadata['tag'] = tag

        override_default_vasp_params = override_default_vasp_params or {}

        vasp_input_set = vasp_input_set or BornChargeSet(structure, isif=isif, **override_default_vasp_params)
        site_properties = deepcopy(structure).site_properties

        # Avoids delivery (prev_calc_loc == '' (instead by True))
        t = []
        if type(prev_calc_loc) == str:
            t.append(CopyVaspOutputs(calc_dir=prev_calc_loc, contcar_to_poscar=True))
            t.append(WriteVaspFromIOSetPrevStructure(vasp_input_set=vasp_input_set, site_properties=site_properties))
        elif parents:
            if prev_calc_loc:
                t.append(CopyVaspOutputs(calc_loc=prev_calc_loc, contcar_to_poscar=True))
            t.append(WriteVaspFromIOSetPrevStructure(vasp_input_set=vasp_input_set, site_properties=site_properties))
        else:
            t.append(WriteVaspFromIOSetPrevStructure(structure=structure, vasp_input_set=vasp_input_set, site_properties=site_properties))
        if (scale_lattice is not None):
            t.append(ScaleVolumeTransformation(scale_factor=scale_lattice, structure=structure))

        #the following statement may not correct to Born effective charge calculation, so be commented.
        """
        t.append(ModifyIncar(incar_update=">>incar_update<<"))
        if modify_incar != None:
             t.append(ModifyIncar(incar_update=modify_incar))
        """
        #t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<", gzip_output=False))
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=False, gzip_output=False))
        t.append(PassCalcLocs(name=name))
        
        t.append(BornChargeToDb(db_file=">>db_file<<", tag=tag))
        #t.append(CheckSymmetryToDb(db_file=db_file, tag=tag))
        super(BornChargeFW, self).__init__(t, parents=parents, name="{}-{}".format(
            structure.composition.reduced_formula, name), **kwargs)



=============
Settings file
=============

The user can change the settings for dfttk calculations by providing a settings file.


File name instruction
=====================

Globale settings
----------------

By default, the global settings file is named as ``SETTINGS.yaml`` or ``SETTINGS.json``.

The user can change the name by ``-s`` parameter in ``dfttk run``. e.g. if the user run dfttk by ``dfttk run -s SET``, then the files named with ``SET.yaml`` or ``SET.json`` is the global settings file.

Individual settings
-------------------

The user can provide individual settings for some structures whoes settings is different with others. The individual settings file should be named with ``SETTINGS-FILENAME_WITH_EXT.yaml(json)`` or ``FILENAME_WITH_EXT-SETTINGS.yaml(json)``

Note: 
    - ``FILENAME_WITH_EXT`` is the file name of structure with extension.

    - ``SETTINGS`` is case insensitive, both ``settings`` and ``SETTINGS`` are OK.

e.g. In the working folder, there are some structure files as follows:

.. code-block::
       
    Fe3Ni.cif
    POSCAR
    POSCAR-2

Then the following files will be recognized as settings files.

.. code-block::
       
    SETTINGS-Fe3Ni.cif.yaml
    POSCAR-SETTINGS.yaml
    settings-POSCAR-2.yaml

Keywords in the file
====================

Default settings
----------------

The settings file is optional, if the user does not provide the settings file, the default value will be used.

- Common settings

+----------------------------+---------------+-------------------------------------------------------+
|Keywords                    | Default Value |Comments                                               |
+----------------------------+---------------+-------------------------------------------------------+
|magmom                      | No            |The MAGMOM for each atom, e.g. [4, 4, -4, 4].          |
|                            |               |*Note*, the length must agree with the number of atoms |
+----------------------------+---------------+-------------------------------------------------------+
|metadata                    | No            |The metadata of the calculation. If the user provides  |
|                            |               |it, DFTTK will find the existing calculations in the   |
|                            |               |databasse. If not, it will generate by uuid4.          |
+----------------------------+---------------+-------------------------------------------------------+
|isif4                       | False         |If run ISIF=4 following ISIF=7 in RobustOptmizeFW      |
+----------------------------+---------------+-------------------------------------------------------+
|level                       | 1             |Optimize level.                                        |
|                            |               |If run ISIF=2 after last ISIF=4 in RobustOptimizeFW    |
|                            |               | 1 for run and 2 for not                               |
+----------------------------+---------------+-------------------------------------------------------+
|override_symmetry_tolerances| None          |Override the default symmetry tolerance, if None,      |
|                            |               |{'tol_strain':0.05,'tol_energy':0.025, 'tol_bond':0.10}|
+----------------------------+---------------+-------------------------------------------------------+
|override_default_vasp_params| {}            |Override the default vasp settings                     |
|                            |               |The optional keys is 'user_incar_settings',            |
|                            |               |'user_kpoints_settings',                               |
|                            |               |'user_potcar_functional'.                              |
|                            |               |For more details, ref.                                 |
|                            |               |https://pymatgen.org/pymatgen.io.vasp.sets.html        |
+----------------------------+---------------+-------------------------------------------------------+
|modify_incar_params         | {}            |Modify the incar settings in the fireworks level. e.g. |
|                            |               |modify_incar_params =                                  |
|                            |               |{ 'Full relax': {'incar_update': {"LCHARG":False}},    |
|                            |               |'static': {'incar_update': "LAECHG":True}}             |
+----------------------------+---------------+-------------------------------------------------------+
|modify_kpoints_params       | {}            |Modify the kpoints settings in the fireworks level.    |
+----------------------------+---------------+-------------------------------------------------------+
|store_volumetric_data       | False         |Store the volumetric data (True) or not (False)        |
+----------------------------+---------------+-------------------------------------------------------+
|verbose                     | False         |print(True) or not(False) some informations,  for debug|
+----------------------------+---------------+-------------------------------------------------------+
|passinitrun                 | False         |Pass init vasp result.                                 |
|                            |               |**It will be dropped in the future**                   |
+----------------------------+---------------+-------------------------------------------------------+
|run_isif2                   | False         |If run ISIF=2 before ISIF=4 (True) or not (False).     |
|                            |               |**It will be dropped in the future**                   |
+----------------------------+---------------+-------------------------------------------------------+
|pass_isif4                  | False         |Whether pass isif=4 calculation.                       |
|                            |               |**It will be dropped in the future**                   |
+----------------------------+---------------+-------------------------------------------------------+
|symmetry_tolerance          | 0.05          |The tolerannce for symmetry.                           |
|                            |               |**It will be dropped in the future**                   |
+----------------------------+---------------+-------------------------------------------------------+
|relax_path                  |''             |The path of relaxiation.                               |
|                            |               |**It will be dropped in the future**                   |
+----------------------------+---------------+-------------------------------------------------------+

- Phonon settings

+----------------------------+---------------+-------------------------------------------------------+
|Keywords                    | Default Value |Comments                                               |
+----------------------------+---------------+-------------------------------------------------------+
|phonon                      | False         |Run phonon (True) or not(False, Debye model)           |
+----------------------------+---------------+-------------------------------------------------------+
|phonon_supercell_matrix     |atoms          |The supercell matrix for phonon calculations.          |
|                            |               |It can take the following values:                      |
|                            |               |**Matrix**, e.g. [[2, 0, 0], [0, 2, 0], [0, 0, 2]]     |
|                            |               |**String**: atom/lattice/volume(the first letter works)|
|                            |               |Determining the supercell matrix automatically         |
|                            |               |by atoms/lattice/volume ranges                         |
+----------------------------+---------------+-------------------------------------------------------+
|phonon_supercell_matrix_min |60             |The lower boundary for phonon_supercell_matrix(String) |
+----------------------------+---------------+-------------------------------------------------------+
|phonon_supercell_matrix_max |130            |The upper boundary for phonon_supercell_matrix(String) |
+----------------------------+---------------+-------------------------------------------------------+
|force_phonon                |False          |Force run phonon (True) or not(False),                 |
|                            |               |No matter ISIF=4/stable_tor pass or not                |
+----------------------------+---------------+-------------------------------------------------------+
|stable_tor                  |0.01           |Stable torlerance (The percentage of negative dos),    |
|                            |               |If the negative part of DOS is larger than this value, |
|                            |               |DFTTK won't run phonon for this structure.             |
+----------------------------+---------------+-------------------------------------------------------+

- QHA settings

+--------------------+---------------+---------------------------------------------------------+
|Keywords            | Default Value |Comments                                                 |
+--------------------+---------------+---------------------------------------------------------+
|num_deformations    | 7             |The number of deformations/structures                    |
+--------------------+---------------+---------------------------------------------------------+
|deformation_fraction|[-0.1, 0.1]    |The range of deformation, 0.1 means 10%                  |
+--------------------+---------------+---------------------------------------------------------+
|eos_tolerance       | 0.01          |The tolerance for eos fitting. If larger than this value,|
|                    |               |DFTTK will append volumes automatically                  |
+--------------------+---------------+---------------------------------------------------------+
|t_min               | 5             |The mimimum of temperature in QHA process                |
+--------------------+---------------+---------------------------------------------------------+
|t_max               | 2000          |The maximum of temperature in QHA process                |
+--------------------+---------------+---------------------------------------------------------+
|t_step              | 5             |The step of temperature in QHA process                   |
+--------------------+---------------+---------------------------------------------------------+
|volume_spacing_min  | 0.03          |Minimum ratio of Volumes spacing.                        |
|                    |               |This keyword will be dropped in the future               |
+--------------------+---------------+---------------------------------------------------------+


- Elastic settings

+---------------+---------------+---------------------------------------------------------+
|Keywords       | Default Value |Comments                                                 |
+---------------+---------------+---------------------------------------------------------+
|strain_states  |None           |Strain modes, if it is None, it will generated by atomate|
+---------------+---------------+---------------------------------------------------------+
|stencils       |None           |The amplitude of the strain modes/states                 |
+---------------+---------------+---------------------------------------------------------+
|analysis       |True           |Analysis (True) or not (False)                           |
+---------------+---------------+---------------------------------------------------------+
|sym_reduce     |False          |Reduce the strain according to the symmetry or not       |
+---------------+---------------+---------------------------------------------------------+
|order          |2              |The order of the elastic constants                       |
+---------------+---------------+---------------------------------------------------------+
|conventional   |False          |Convert the structure into conventional format or not    |
+---------------+---------------+---------------------------------------------------------+




=======
Recipes
=======


Construct a series of Debye workflows by substituting into a template structure
===============================================================================

   ##############
   # User Input #
   ##############

 .. code-block:: console

   DRY_RUN = True  # Don't submit the workflows
   VERBOSE = True   # Turn on printing of substitutions

   # Filename of the template structure to use, usually a dilute or SQS structure
   TEMPLATE_STRUCTURE_FILENAME = 'structures/FCC_A1.sqs.12Ni-4Fe.POSCAR'

   # sublattice configuration of the template structure.
   # This will be substitued exactly, this does not need to be sorted,
   # however the individual configurations to build should be sorted.
   TEMPLATE_SUBLATTICE_CONFIGURATION = [['Fe', 'Ni']]
   TEMPLATE_SUBLATTICE_OCCUPANCIES = [[0.25, 0.75]]
   SUBLATTICE_SITE_RATIOS = [1.0]

   PHASE_NAME = 'FCC_A1'

   configurations_to_build = [  # list of sublattice configurations in DFTTK format
   [['Cr', 'Ni']],
   [['Fe', 'Ni']],
   # [['V', 'Ni']],  # should not use this because it's out of order, make a second script with the templates flipped
   ]

   # Dictionary of densities for each pure element.
   # Not necessary, using the peridic_table in pymatgen instead
   #DENSITY_DICT = {
   #    'V': 6.313,  # bcc
   #    'Cr': 7.463,  # bcc
   #    'Ni': 9.03,  # fcc
   #    'Ti': 4.58,  # hcp
   #    'Fe': 8.028  # bcc
   #}


##########
# SCRIPT #
##########

   # Should not need to edit below this line.

.. code-block:: console

   from pymatgen import Structure
   from fireworks import LaunchPad
   from dfttk import get_wf_gibbs
   from dfttk.structure_builders.substitutions import substitute_configuration_with_metadata

   workflows = []

   temp_struct = Structure.from_file(TEMPLATE_STRUCTURE_FILENAME)

   for config in configurations_to_build:
       struct, meta = substitute_configuration_with_metadata(temp_struct, TEMPLATE_SUBLATTICE_CONFIGURATION, config, TEMPLATE_SUBLATTICE_OCCUPANCIES, PHASE_NAME, SUBLATTICE_SITE_RATIOS)
       if VERBOSE:
           print("PHASE: {}    CONFIGURATION: {}    OCCUPANCIES: {}    STRUCTURE: {}".format(PHASE_NAME, config, TEMPLATE_SUBLATTICE_OCCUPANCIES, struct.composition.hill_formula))
       workflows.append(get_wf_gibbs(struct, deformation_fraction=(-0.05,0.10), phonon=False, num_deformations=11, t_max=2000, metadata=meta))


   ################################################################################
   # Load everything in the LaunchPad
   ################################################################################

   if VERBOSE:
       print("{} workflows.".format(len(workflows)))

   if DRY_RUN:
       exit()

   if VERBOSE:
       print('Adding workflows to LaunchPad')

   lpad = LaunchPad.auto_load()

   for workflow in workflows:
       lpad.add_wf(workflow)


ESPEI datasets from a QHA database
==================================

The following code snippet will take ESPEI datasets from a QHA database,
optionally writing the (nicely named) files to dict.
The QHA database requires the following metadata schema:

.. code-block:: console

   {
     'metadata': {
       'phase_name': 'FCC_A1',
       'tag': 'ed447049-ad67-4090-ba99-378188d3416b',
       'sublattice': {
         'configuration': [['Cr', 'Ni']],
         'occupancies': [[0.03125, 0.96875]]
        },
     }
   }

.. code-block:: python

   ########
   # EDIT #
   ########

   phase_name = 'BCC_A2'
   configuration_to_find = [['Ni', 'V']]
   sublattice_site_ratios = [1.0]
   db_username = 'BrandonResultsRO'
   db_password = 'piqhg38hap3'
   db_uri = 'mongodb://206.189.190.225:27018'
   WRITE_FILES = True

   temperature_index = 59  # index of 300 K temperature (close to 298 K), found by hand

   refstate_tags = {
       'Fe': '4ac77fce-0e43-4c07-8418-4843a2cd5723',
       'Ni': '0059ee69-4a8f-4e86-9895-9b40cf67dd96',
       'Cr': '8cceb186-2796-4488-ba8c-2380c5278f62',
       'V': 'fba46b6b-1699-419f-b5e1-da9533530701',
   }

   ########################
   ### SCRIPT
   ########################

   from dfttk.analysis.formation_energies import get_formation_energy, get_thermal_props
   from dfttk.espei_compat import make_dataset, dfttk_config_to_espei, dfttk_occupancies_to_espei
   from pymatgen import Structure
   import numpy as np
   from dfttk.utils import recursive_flatten
   import json
   from pymongo import MongoClient

   # create the MongoClient
   cli = MongoClient(db_uri)
   db = cli.results
   db.authenticate(name=db_username, password=db_password)
   coll = db.qha

   # construct the energies, assupmtion of same temperature grid
   # energies are J/mol-atom
   refstate_energies = {}
   for el, tag in refstate_tags.items():
       qha_result = coll.find_one({'metadata.tag': tag})
       refstate_energies[el] = get_thermal_props(qha_result)

   # calculate all the dataset values
   configs     = []
   occupancies = []
   hm_values   = []
   sm_values   = []
   cpm_values  = []

   # we'll change the T to the right temperatures later
   fixed_conds = {'P': 101325, 'T': 0}
   temp_conds = {'P': 101325, 'T': 0}


   for qha_res in coll.find({'metadata.sublattice.configuration': configuration_to_find, 'metadata.phase_name': phase_name}):
       configs.append(qha_res['metadata']['sublattice']['configuration'])
       occupancies.append(qha_res['metadata']['sublattice']['occupancies'])

       tprops = get_thermal_props(qha_res)
       struct = Structure.from_dict(qha_res['structure'])
       hm_form = get_formation_energy(tprops, struct, refstate_energies, 'HM', idx=temperature_index)
       sm_form = get_formation_energy(tprops, struct, refstate_energies, 'SM', idx=temperature_index)
       cpm_form = get_formation_energy(tprops, struct, refstate_energies, 'CPM', thin=10)[:-2]
       fixed_temp = tprops['T'][temperature_index]
       cpm_temps = tprops['T'][::10][:-2]

       hm_values.append(hm_form)
       sm_values.append(sm_form)
       cpm_values.append(cpm_form)

   fixed_conds['T'] = fixed_temp.tolist()
   temp_conds['T'] = cpm_temps.tolist()

   # make the HM, SM, CPM values arrays of the proper shape
   hm_values = np.array([[hm_values]])
   sm_values = np.array([[sm_values]])
   cpm_values = np.array(cpm_values).T[np.newaxis, ...]

   if WRITE_FILES:
       # write JSON files
       comps = [c.upper() for c in sorted(recursive_flatten(configuration_to_find))]
       for prop, vals, conds in [('HM_FORM', hm_values, fixed_conds), ('SM_FORM', sm_values, fixed_conds), ('CPM_FORM', cpm_values, temp_conds)]:
           ds = make_dataset(phase_name, prop, sublattice_site_ratios, configs, conds, vals, occupancies=occupancies, tag=tag)
           with open('{}-{}-{}-DFTTK.json'.format('-'.join(comps), phase_name, prop), 'w') as fp:
               json.dump(ds, fp)


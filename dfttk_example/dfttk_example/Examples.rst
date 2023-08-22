Examples
=========

The Examples folder is designed to keep the data for user to test the DFTTK package. The two subfolder are::

 - dir ``Al/`` - contains all data for test run for Al
 - dir ``ZrSiO4/`` - contains all data for test run for ZrSiO4
 - dir ``ExptData/`` - contains a json file ``ExptData.json`` which saves the experimental thermodynamic data for a collection of materials.

Within the ``Al`` or the ``ZrSiO4`` subfolder, one see two subfloders

 - dir ``input/`` - contain input setup files using ``Al`` as the example on E-V, phonon, and thermodynamic property calculations
 - ``Al_Fm-3m_225PBE/`` - contain outputs by postprocessing data that saved in MongoDB by the above ``Al`` example.

For the data within Al_Fm-3m_225PBE/

 - dir ``Yphon/`` - all data input/output for Yphon, e.x., hessian matrix (superfij.out), calculated phonon dos
 - dir ``figures/`` - plots in png format for most of the thermodynamic properities
 - file ``readme`` - extensive summary of the calculated results in json format
 - file ``fvib_ele`` - tablated data containing the calculated thermodynamic properties
 - file ``fvib_eij`` - tablated data containing the calculated thermal expansion coefficient tensor 
 - file ``record.json`` - SGTE fitting record for heat capacity, Gibbs energy, enthalpy, and entropy at given temperature range

To run the Example for the VASP calculation, say ``Al``, run

.. code-block:: bash

 cd Al/input
 dfttk run -wf robust -f POSCAR.Al -l -m 1

To postprocess calculations after the VASP calculation done, run

.. code-block:: bash

 cd Al
 dfttk thfind -py -td -50 -plot find_or_DFT -eq 4 -smooth -el 1 -renew -get -metatag 0c1887fa-0cb4-4ce2-9559-7f7909ffa11a -expt ../ExptData/ExptData.json
 #note that the key ``0c1887fa-0cb4-4ce2-9559-7f7909ffa11a`` is obtained from the file ``input/METADATAS.yaml`` automatically produced by the VASP calculation step.

The above will produce more thatn 20 figures stored in the folder “Al_Fm-3m_225PBE/figures/” and they
can be view them using the linux command ``display`` to show the figure, for example


.. code-block:: bash

  display Al_Fm-3m_225PBE/figures/LTC.png #to see the linear thermal expansion coefficient

.. image:: _static/Al-LTC.png

.. code-block:: bash

  display Al_Fm-3m_225PBE/figures/Heat_capacities.png #to see the heat capaticity, and so on

.. image:: _static/Al-Heat_capacities.png



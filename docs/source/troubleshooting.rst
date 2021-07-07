***************
Troubleshooting
***************


This document covers how to handle jobs that have fizzled. There are three main sections: common troubleshooting, general troubleshooting workflow to find the root cause of issue and specific issues and their fixes.

Common troubleshooting
======================

1.      In the config stages, one may need properly set up VASP environments, such as

.. code-block:: bash

    module load intel impi vasp

2.      Make the latest `automate <https://atomate.org/>`_ installed;
3.      In ``.bashrc`` or ``.cshrc``, make sure not messed up with your other `atomate` FW config;
4.      Make sure the setup in the ``.bashrc`` file 

.. code-block:: bash

    export FW_CONFIG_FILE=~/dfttk/config/FW_config.yaml

or equivalently in the ``.cshrc`` file

.. code-block:: bash

    setenv FW_CONFIG_FILE /storage/work/y/yuw3/dfttk/config/FW_config.yaml

5.      git push issue for contributors, see https://docs.github.com/en/github/authenticating-to-github/adding-a-new-ssh-key-to-your-github-account
6.      for batch run of the postprocessing modules, make sure compatibilities of non-ascii character by:

.. code-block:: bash

    export LC_ALL='en_US.utf8' #for bsh;
    setenv LC_ALL en_US.utf8 #for csh

7.      For phonon calculations, due to certain reasons (such as temperature range too high), one may not see results in the ``qha_phonon`` or ``qha`` MongoDB collections. In this case, the subcommand ``dfttk thfind`` will try to find results from the ``phonon`` collection and process the data by calling ``Yphon``

8.      When you are interesting in revising the code, if have job running in the system before your changes, the codes in the batch system might not be updated and the results might be not as you assumed. It takes me two days to figure out this problem. The solution is to kill all the dfttk running job and resubmit them.
Following are the steps of adding API key number on DFTTK.

9. How to solve the install warning of 'PMG_MAPI_KEY' is empty.

  1. Go to the materials project website,
  https://materialsproject.org/, under the API section, you will
  easily find you API Keys number.

  2. Go to the .pmgrc.yaml file

  .. code-block:: bash

      vi ~/.pmgrc.yaml

  3. Add your API key number into your .pmgrc.yaml file. For example

  .. code-block:: bash

      PMG_MAPI_KEY: ######(your API key number)


pymatgen 2021 issue
===================

You mag meet numpy version issues using pymatgen, reporting::

    pymatgen 2021.2.16 requires numpy>=1.20.1, but you'll have numpy 1.19.2 which is incompatible.


In such case, please upgrade numpy by::

    pip install numpy --upgrade

conda issues
============

In some cases, such as in the Windows environment, one may meet the error::

    ModuleNotFoundError: No module named 'ruamel' #106

This is due to ``conda`` bug on namespace of ruamel_yaml vs ruamel.yaml. 
 One can resolve this by open the Annaconda Powershell Prompt as adminstrator and reinstall ruamel.yaml by::

    conda install ruamel.yaml


Troubleshooting Workflow
========================


**My job has fizzled!** The following steps can help you get information about how your job. You can imagine it as a decision tree. Check one thing before moving on to the next one.

1. Check that the job ran and has raised an exception with a traceback.

   Run ``lpad get_fws -i <ID> -d more``, replacing ``<ID>`` with the integer id of the Firework.
   Search the output for ``_exception``.
   What you see is the Python exception that was raised when running the Firework.

   *TIP:*: Searching works well when you pipe the output to ``less`` with ``lpad get_fws -i <ID> -d more | less`` and search using ``/``.

   *todo:*: If you don't see a traceback, that means... (this is the first step, but does this actually happen?)


2. Check the traceback is not a common error.

   See the `Common Errors section <CommonErrors>`_


.. _CommonErrors:

Common Errors
=============

Custodian VasprunXMLValidator failed
------------------------------------

In this error, you get a traceback that looks something like:

.. code-block:: python

   Traceback (most recent call last):
     File "/storage/home/bjb54/.conda/envs/wfs/lib/python3.7/site-packages/custodian/custodian.py", line 320, in run
       self._run_job(job_n, job)
     File "/storage/home/bjb54/.conda/envs/wfs/lib/python3.7/site-packages/custodian/custodian.py", line 428, in _run_job
       raise CustodianError(s, True, v)
   custodian.custodian.CustodianError: (CustodianError(...), 'Validation failed: <custodian.vasp.validators.VasprunXMLValidator object at 0x2af45b1d3908>')

   During handling of the above exception, another exception occurred:

   Traceback (most recent call last):
     File "/storage/home/bjb54/.conda/envs/wfs/lib/python3.7/site-packages/fireworks/core/rocket.py", line 262, in run
       m_action = t.run_task(my_spec)
     File "/storage/home/bjb54/.conda/envs/wfs/lib/python3.7/site-packages/atomate/vasp/firetasks/run_calc.py", line 204, in run_task
       c.run()
     File "/storage/home/bjb54/.conda/envs/wfs/lib/python3.7/site-packages/custodian/custodian.py", line 330, in run
       .format(self.total_errors, ex))
   RuntimeError: 0 errors reached: (CustodianError(...), 'Validation failed: <custodian.vasp.validators.VasprunXMLValidator object at 0x2af45b1d3908>'). Exited...


With the key being that Custodian fails to validate the ``vasprun.xml``. After running VASP, Custodian will try to parse the ``vasprun.xml`` file using pymatgen.

There are usually two possible triggers for this failure:

1. VASP failed to run at all (more common) or quit in a way that custodian did not detect (less common)
2. The ``vasprun.xml`` file could not be parsed by pymatgen.

To investigate this, first check that VASP ran (e.g. the ``OUTCAR`` shows that the run completed successfully).
If VASP did not run, find out why and fix that issue.
If VASP did run successfully, it was probably an issue parsing the ``vasprun.xml`` file.
Try parsing the ``vasprun.xml`` file using the ``pymatgen.io.vasp.outputs.Vasprun`` class.
If it throws an error when you try to parse, that's what made Custodian fail and you should fix that.

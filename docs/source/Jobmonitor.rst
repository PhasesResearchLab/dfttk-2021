==============
Job monitoring
==============

The following briefes some common ``automate`` commands for DFT jobs submitted by ``dfttk run``.

- Monitor workflows

To check running status of all submitted dfttk jobs, use

.. code-block:: bash

  lpad get_wflows

To find the running calculations

.. code-block:: bash

  lpad get_fws -s RUNNING 

To find the jobdir for a job with specific fw_id (The number portion of job id reported by ``lpad get_wflows``), use 
lpad get_launchdir  fw_id

- FIZZLED jobs

``FIZZLED`` means failed. The following command can summarize all FIZZLED jobs:

.. code-block:: bash

  lpad get_fws -s FIZZLED; or
  lpad get_wflows -s FIZZLED

One can rerun ``FIZZLED`` by 

.. code-block:: bash

  lpad rerun_fws -s FIZZLED

For more details, see ref. FIZZLED.

To get help for the subcommands, try ``lpad get_fws -h`` or ``lpad get_wflows -s FIZZLED -h``



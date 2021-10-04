=========
Changelog
=========

0.3.4 (2021-10-04)
==================
046cc64 Merge branch '20210826'
20dbd4f mark FM and AFM
27acc3f Update input_sets.py
c14d283 Update fworks.py
c1d1ea6 Update assign_fworker_name.py
a884245 Update EVcheck_QHA.py
dfcb2fc Update utils.py
fbbbdf0 magnetic correction
fd30553 Update test_thfind.py
f367b48 Update pytest.ini
2f7fa79 for vasp6 and born effective charge
5e1b9b6 code reorder
8daff7d for compatibility with vasp 6
8340474 performance improving
10b0a65 Update run_dfttk.py
34f9348 Update run_task_ext.py
414fdf1 update for elasticity
179b763 update elastic code
19f8fe8 Update parse_outputs.py
315750e debug elastic
3eb9511 Update parse_outputs.py
b9d9792 fix bugs on EVcheck and add a pytest therefore
e3ccef3 pull request for merging changes on Debye model (#8)
8d46358 Update troubleshooting.rst
9164935 bug fixes on multiple computer running DFTTK
d1472c7 for qha results repairing
dce5198 accumulative changes after manuscript submission
ba3fc50 Add Yphon for phonon supercell and bug fix
ba6888f Debye-Gruneisen model
ead37f9 pull request before correcting Debye model (#6)
93beaa3 Update querydb.py
bf88aeb Phonon set to 3x3x3
3266f21 for non selfconsistent calculations
7649c21 add store_raw_vasprunxml for not save vasprun.xml by defaul
3e77778 more update on Born
06793fb config update
ede041c Update thermal_electronic.py
927c1eb docs update

0.3.0 (2020-12-10)
==================

(Contributor: @YiWang, @mxf469, @bocklund)
- Change List:
 - formally release dfttk version 0.3.0
 - Revised dfttk/script/run_dfttk.py
 - made workflow "get_wf_gibbs" run on "robust" with get_wf_gibbs_robust changed test module from get_wf_gibbs into get_wf_gibbs_robust
 - Add elasticity calculation model
 - Add Born effective charge calculations

0.2.2 (2020-08-18)
================

(Contributor: `@YiWang`_, `@mxf469`_)

- Change List:
 - added codes in the dfttk/scripts directory:
  - run_dfttk_ext.py
   - handle the argumetns for the thelec and thfind modules

 - added python code in the dfttk directory:
  - pyfind.py
   - database search engine

  - pythelec.py
   - for compatibiliy with Yphon
   - generating the majority of thermodynamic properties, such as thermal expansion coefficient, Seebech coefficients, Lorenz number etc

  - pyphon.py for
   - calculate the phonon contributions to the various thermodynamic properties

 - added python code in the dfttk/analysis directory:
  - database
   - for plot phonon dispersions for all crystalline systems

  - ywutils.py
   - general utils code

  - ywplot.py
   - for plots of ~20 different phonon and thermodynamic properties in the png format 

* made Yphon compatibile with phonopy
 - added codes in the CRO-soc directory:
   - phonopy2yphon, phonopy2yphon.py
    - convert the phonopy force constant matrix in hdf5 format into superfij.out format used by Yphon

 - changed codes:
  - in the dfttk/scripts directory:
   - run_dfttk.py
    - added the following lines aimed to handle the argumetns for the thelec and thfind modules

    # extension by Yi Wang, finalized on August 4, 2020
    # -----------------------------------
    from dfttk.scripts.run_dfttk_ext import run_ext_thelec
    run_ext_thelec(subparsers)

  - in the dfttk/analysis directory:
   - debye.py is renamed as debye_ext.py
    - to include the vibrational entropy (S_vib) and heat capacity (C_vib) into the "qha" MongoDB collection

   - quasiharmonic.py:
    - copy the S_vib and C_vib from the "phonon" collection into the "qha_phonon" MongoDB collection

0.2 (2020-03-30)
================

New features

(Contributor: `@bocklund`_ , @Peng_Gao, `@hitliaomq`_ )

* The relax scheme is optimized. (from ``ISIF=3`` to ``ISIF=2`` followed by ``ISIF=4``) (@Peng_Gao)
* Change the static workflow to dynamic workflow. (``EVcheck_QHA.py`` increase the data points atomately if the fitting of initial points is incorrect) (@Peng_Gao)
* Support run dfttk by command. (Add ``dfttk run [options]``) (`@hitliaomq`_)
* Support configrate dfttk automately. (Add ``dfttk config [options]``) (`@hitliaomq`_)
* Documents' enhance. (`@hitliaomq`_)
* Bug fix. (Including `#8`_ ) (`@bocklund`_, @Peng_Gao, `@hitliaomq`_)

.. _`#8`: https://github.com/PhasesResearchLab/dfttk/issues/8

0.1 (2018-08-28)
================

Initial release. Includes

(Contributor: `@bocklund`_, `@mxf469`_)

* Gibbs workflow for stable structures
* Analysis code and libraries for calculation quasiharmonic Gibbs energies with 0K, vibrational and thermal electronic contributions
* Useful utilities for interfacing with structure, calculations and the Materials Project

.. _`@bocklund`: https://github.com/bocklund
.. _`@mxf469`: https://github.com/mxf469
.. _`@hitliaomq`: https://github.com/hitliaomq
.. _`@YiWang`: https://github.com/yiwang62

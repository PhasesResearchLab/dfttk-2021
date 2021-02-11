.. title::FAQ

FAQ
===


What do the references to DFTTK-style or ESPEI-style configurations or occupancies mean?
----------------------------------------------------------------------------------------

`ESPEI <http://espei.org>`_ is separate software for fitting thermodynamic
databases within the CALPHAD method. The CALPHAD method requires descriptions of
the formation and mixing energies of phases described by sublattices, as in the
compound energy formalism.

DFTTK is an excellent tool for calculating the temperature and composition
dependent thermodynamic properties of sublattice configurations within the
compound energy formalism.

In both ESPEI and DFTTK, there are three key aspects to describing a phase in the
CEF: sublattice site ratios, sublattice configurations, and sublattice occupancies.
Sublattice site ratios are a list of integer or fraction ratios of the sublattice.
There are no constraints on the order or values. Some valid ones would be:

.. code-block::

   [1]  # 1 sublattice
   [3.0, 2.0]  # 2 sublattices with 3 occurrences of species in sublattice 1 and 2 in sublattice 2
   [1, 0.5]  # 2 sublattices with 1 occurrence of species in sublattice 1 and 1/2 in sublattice 2
   [16, 11]  # 2 sublattices with 16 occurrences of species in sublattice 1 and 11 in sublattice 2

DFTTK describes a sublattice configuration as a list of sublattices,
where each sublattice is also a list of the components occupying that sublattice.
Some examples are:

.. code-block::

   [[Fe, Ni]]  # 1 sublattice with mixing of Fe and Ni
   [[Mg], [Cu]]  # 2 sublattices with no mixing in each sublattice.
   [[Cr, Fe], [Cr], [Cr, Fe]]  # 3 sublattices with mixing in the first and last sublattice

The ordering of sublattices themselves have no distinct meaning. As a convention,
DFTTK sorts the components in each sublattice by component name, e.g., Cr
(starts with C) always comes before Fe (starts with F)). ESPEI requires
each sublattice to be ordered in this way.

Sublattice occupancies or occupations describe how a sublattice is filled with
each particular element. The shape of these corresponds exactly to the configuration
that the occupancy list describes. The sorting of the occupancies therefore is
also in the same order of the configuration.

DFTTK includes several tools to help build ESPEI-compatible sublattice descriptions.
The only difference between DFTTK and ESPEI conigurations and occupancies is that
a singly occuped sublattice, e.g. [[Fe]] in DFTTK, is not wrapped in the sublattice
parenthesis for both the configuration and occupancies. Some examples of ESPEI configuration
and occupancy lists are given below:

.. code-block::

   [[FE, NI]], [[0.5, 0.5]]  # this is the same as DFTTK, there are no single occupations
   [MG, CU], [1.0, 1.0]  # each sublattice is singly-occupied and is not wrapped in the parenthesis
   [[CR, FE], CR, [CR, FE]], [[0.5, 0.5], 1, [0.3, 0.7]] # the second sublattice is singly occupied and is not wrapped

As a final note, ESPEI format typically uses all uppercase, consistent with the
CALPHAD community, while DFTTK usually uses standard element capitalization. The
conversion tools take this into account.

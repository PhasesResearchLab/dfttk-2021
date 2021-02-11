==============================
Dfttk postprocessing mechanism
==============================

Flow controls
-------------

Depending on the given condition, postprocessing can take data from

    | The ``qha-phonon`` collection
    | The ``qha`` collection (Debye model)
    | The ``phonon`` collection (in case of qha-phonon calculations failed)

Then, invoke ``Yphon`` to recalculate the phonon properties based on the force constants obtained from the database 

Scheme to find equilibrium volume
---------------------------------

    | If the data quality is excellent, use central symmetric 7-point difference
    | If the data quality is very good, fit the free energies use 4-parameter Birch-Murnaghan 
    | If the data quality is good (may use as default), fit the 0-K total energies using 4-parameter Birch-Murnaghan 
    | Fit the finite T part of the free energies by UnivariateSpline

Scheme to calculate derivative
------------------------------

 1. 7\-point symmetrical central difference

  .. math::

    Deriv= \frac{1}{3}\sum_{i=1}^{3}{\frac{f(X_{N+i})-f(X_{N-i})}{X_{N+i}-X_{N-i}}}


 2. Birch-Mannhan Euqations of state fitting

  .. math::

    Deriv=-\frac{2}{3}bx^{-\frac{5}{3}}-\frac{4}{3}cx^{-\frac{7}{3}}-\frac{6}{3}dx^{-\frac{9}{3}}

Scheme for LTC
--------------

 1. By entropy derivative

  .. math::

    \alpha =\frac{1}{B_{T}}\frac{\partial S}{\partial V}
    
where :math:`B_{T}` is isothermal bulk modulus

 2. By internal energy derivative via Birch-Mannhan fitting

  .. math::

    \alpha =\frac{1}{B_{T} T}(\frac{\partial U}{\partial V}+P)

where :math:`P` is pressure


 3. When the data quality is fair (~20% cases)
 
  | Fit the 0-K total energies using 4-parameter Birch-Murnaghan; and 
  | fit the finite T part of the free energies by linear function f=a+b*V



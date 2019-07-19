.. index::
   single: Program; FFPT
   single: FFPT
   single: Properties; Finite-field PT
   single: Properties; FFPT

.. _TUT\:sec\:ffpt:

:program:`FFPT` --- A Finite Field Perturbation Program
========================================================

Many molecular properties of wave functions can be computed using the
:program:`FFPT` program module in |molcas|. It adds the requested operator to
the integrals computed by the :program:`seward` module. This must be done
before the |molcas| module calculating the required wave function is requested
so the :program:`FFPT` module is best run directly after the :program:`seward`
module.

The :kword:`TITLe` keyword behaviors in a similar fashion to other
|molcas| modules.
The sample input below contains the :program:`FFPT` input
requesting that the dipole moment operator be added to the integrals
using the :kword:`DIPOle` keyword.
The size and direction is specified using the :kword:`COMP` keyword
which accepts free format input. We can compute the dipole of the
molecule by numerical determination of the gradient of the energy
curve determined for several values of the dipole operator. From the second
derivative we can obtain the polarizability component.

Sample input requesting the FFPT module to
include a dipole moment operator in the integral file: ::

  &FFPT
  Title= Finite Perturbation with a dipole in the x negative of strength 0.1 au
  FFPT
  Dipole
   Comp
   X -0.1

:program:`ffpt` Output
----------------------

The :program:`ffpt` section of the output is short and self
explanatory. The :file:`ONEINT` file is updated with the requested
operator.

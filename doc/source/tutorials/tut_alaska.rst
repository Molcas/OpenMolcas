.. index::
   single: Program; Alaska
   single: Alaska

:program:`ALASKA` --- A Program for Integral Derivatives
=======================================================

:program:`ALASKA` computes the first derivatives of the one- and two-electron
integrals with respect to the nuclear displacements. The derivatives are contracted
with the one- and two-electron densities to form the molecular gradients, which
will be used by the program :program:`SLAPAF`. At present the :program:`ALASKA`
module computes SCF/DFT and MCSCF gradients analytically, the rest are computed
numerically. The :program:`ALASKA` module is automatically invoked when needed if
the user has not explicitly requested the module to be executed. We postpone the
discussion about :program:`ALASKA` to section :ref:`TUT:sec:structure`.

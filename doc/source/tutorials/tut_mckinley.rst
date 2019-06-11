.. index::
   single: Program; McKinley
   single: McKinley

.. _TUT\:sec\:mckinley:

:program:`MCKINLEY` --- A Program for Integral Second Derivatives
=================================================================

:program:`MCKINLEY` computes the analytic second derivatives of the one- and two-electron
integrals with respect to the nuclear positions at the SCF and CASSCF level of theory.
The differentiated integrals
can be used by program :program:`MCLR` to performs response calculations on
single and multiconfigurational SCF wave functions. One of the basic uses
of :program:`MCKINLEY` and :program:`MCLR` is to compute analytical Hessians
(vibrational frequencies, IR intensities, etc).
Note that :program:`MCKINLEY` for a normal frequency calculations will automatically
start the :program:`MCLR` module!
For all other methods a numerical procedure is automatically
invoked by :program:`MCKINLEY` to compute the vibrational frequencies.

:program:`MCKINLEY` --- Basic and Most Common Keywords
------------------------------------------------------

.. class:: keywordlist

:kword:`PERTurbation`
  Suboptions Geometry (for geometry optimizations) or Hessian (full Hessian)

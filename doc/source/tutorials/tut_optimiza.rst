.. index::
   single: Optimization
   single: Geometry

.. _TUT\:sec\:structure:

:program:`ALASKA` and :program:`SLAPAF` --- A Molecular Structure Optimization
==============================================================================

One of the most powerful functions of *ab initio* calculations is geometry
predictions. The minimum energy structure of a molecule for a given method and
basis set is instructive especially when experiment is unable to determine the
actual geometry. |molcas| performs a geometry optimization with analytical
gradients at the SCF or RASSCF level of calculation, and with numerical
gradients at the CASPT2 level.

In order to perform geometry optimization an input file must contain
a loop, which includes several calls: calculation of integrals (:program:`SEWARD`),
calculation of energy (:program:`SCF`, :program:`RASSCF`, :program:`CASPT2`), calculation of gradients
(:program:`ALASKA`), and calculation of the new geometry (:program:`SLAPAF`).

This is an example of such input ::


  &GATEWAY
   coord= file.xyz
   basis= ANO-S-MB
  >> EXPORT MOLCAS_MAXITER=25
  >> Do While <<
  &SEWARD
  &SCF
  &SLAPAF
  >> EndDo <<

The initial coordinates will be taken from xyz file :file:`file.xyz`, and the geometry
will be optimized at the SCF level in this case. After the wave function calculation,
calculation of gradients is required, although code :program:`ALASKA` is automatically
called by |molcas|. :program:`SLAPAF` in this case required the calculation of an
energy minimum (no input). Other options are transition states (:kword:`TS`), minimum energy
paths (:kword:`MEP-search`), etc
The loop will be terminated if the geometry
converges, or maximum number of iterations (:variable:`MOLCAS_MAXITER`) will be reached (the
default value is 50).

There are several EMIL commands
(see :numref:`UG:sec:emil_commands`) which can be
used to control geometry optimization. For example, it is possible to execute
some |molcas| modules only once: ::

  >> IF ( ITER = 1 )
  * this part of the input will be executed only during the first iteration
  >> ENDIF

Program :program:`SLAPAF` is tailored to use analytical or numerical gradients produced
by :program:`ALASKA` to relax the geometry of a molecule towards an energy
minimum (default, no input required then) or a transition state. The program is also used for
finding inter state crossings (ISC), conical intersections (CI),
to compute reaction paths, intrinsic reaction coordinate (IRC) paths, etc.

.. Examples as how to use the :program:`SLAPAF` code is displayed in following :numref:`TUT:sec:structure`.

:program:`SLAPAF` --- Basic and Most Common Keywords
----------------------------------------------------

.. class:: keywordlist

:kword:`TS`
  Computing a transition state

:kword:`FindTS`
  Computing a transition state with a constraint

:kword:`MEP-search`
  Computing a steepest-descent minimum reaction path

:kword:`ITER`
  Number of iterations

:kword:`INTErnal`
  Definition of the internal coordinates

  .. :kword:`IRC`
     Intrinsic reaction coordinate analysis of a TS

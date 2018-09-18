.. index::
   single: Program; MCLR
   single: MCLR

.. _TUT\:sec\:mclr:

:program:`MCLR` --- A Program for Linear Response Calculations
==============================================================

:program:`MCLR` computes response calculations on single and multiconfigurational
SCF wave functions. One of the basic uses of :program:`MCKINLEY` and :program:`MCLR`
is to compute analytical Hessians (vibrational frequencies, IR intensities, etc).
:program:`MCLR` can also calculate the Lagrangian multipliers for
a MCSCF state included in a state average optimization and construct the effective
densities required for analytical gradients of such a state.
The use of keyword :kword:`RLXRoot` in the :program:`RASSCF` program is required.
In both cases the explicit request of executing the :program:`MCLR` module is not
required and will be automatic.
We postpone further
discussion about :program:`MCLR` to section :ref:`TUT:sec:structure`.

It follows an example of how to optimize an excited state from a previous
State-Average (SA) CASSCF calculation.

.. extractfile:: tutorials/MCLR.acrolein.input

  &GATEWAY
  Title= acrolein minimum optimization in excited state 2
  Coord=$MOLCAS/Coord/Acrolein.xyz
  Basis= sto-3g
  Group=NoSym
  >>> Do while
  &SEWARD
  &RASSCF
  Title= acrolein
  Spin= 1; nActEl= 6 0 0; Inactive= 12; Ras2= 5
  CiRoot= 3 3 1
  Rlxroot= 2
  &SLAPAF
  >>> EndDo

The root selected for optimization has been selected here with the keyword
:kword:`Rlxroot` in :program:`RASSCF`, but it is also possible to select it
with keyword :kword:`SALA` in :program:`MCLR`.

Now if follows an example as how to compute the analytical hessian for the lowest
state of each symmetry in a CASSCF calculation (SCF, DFT, and RASSCF analytical
Hessians are also available).

.. extractfile:: tutorials/MCLR.benzoquinone.input

  &GATEWAY
  Title=p-benzoquinone anion. Casscf optimized geometry.
  Coord = $MOLCAS/Coord/benzoquinone.xyz
  Basis= sto-3g
  Group= X Y Z
  &SEWARD
  &RASSCF
  TITLE=p-benzoquinone anion. 2B3u state.
  SYMMETRY=2; SPIN=2; NACTEL=9 0 0
  INACTIVE=8  0  5  0  7  0  4  0
  RAS2    =0  3  0  1  0  3  0  1

  &MCKINLEY; Perturbation=Hessian

The :program:`MCLR` is automatically called after :program:`MCKINLEY`
and it is not needed in the input.

:program:`MCLR` program --- Basic and Most Common Keywords
----------------------------------------------------------

.. class:: keywordlist

:kword:`SALA`
  Root to relax in geometry optimizations

:kword:`ITER`
  Number of iterations

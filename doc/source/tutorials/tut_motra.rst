.. index::
   single: Program; MOTRA
   single: MOTRA
   single: Integrals; Integral transformation

.. _TUT\:sec\:motra:

:program:`MOTRA` --- An Integral Transformation Program
=======================================================

Integrals saved by the :program:`SEWARD` module
are stored in the Atomic Orbital (AO) basis. Some programs have their own
procedures to transform the integrals into the Molecular Orbital (MO) basis.
The |molcas| :program:`MOTRA` module performs this task for
Configuration Interaction (CI), Coupled- and Modified Coupled-Pair (CPF and
MCPF, respectively) and Coupled-Cluster (CC) calculations.

The sample input below contains the :program:`motra` input
information for our continuing water calculation. We firstly specify that the
:program:`RASSCF` module interface file will be the source of the
orbitals using the keyword :kword:`JOBIph`. The keyword
:kword:`FROZen` is used to specify the number of orbitals in each
symmetry which will not be correlated in
subsequent calculations. This can also be performed in the corresponding
:program:`MRCI`, :program:`CPF` or CC programs
but is more efficient to freeze them here.
Virtual orbitals can be deleted using the :kword:`DELEte` keyword. ::

  &MOTRA
  JobIph
  Frozen= 1 0 0 0

:program:`motra` Output
-----------------------

The :program:`motra` section of the output is short and self
explanatory. The integral files produced by :program:`SEWARD`, :file:`ONEINT`
and :file:`ORDINT`, are used as input by the
:program:`MOTRA` module which produces the transformed symbolic files
:file:`TRAONE` and :file:`TRAINT`, respectively. In our case, the files
are called :file:`water.TraOne` and :file:`water.TraInt`, respectively.


The :program:`motra` module also requires input orbitals.
If the :kword:`LUMOrb` keyword is specified the orbitals are taken
from the :file:`INPORB` file which can be any formated orbital
file such as :file:`water.ScfOrb` or :file:`water.RasOrb`. The
:kword:`JOBIph` keyword causes the :program:`MOTRA` module to
read the required orbitals from the :file:`JOBIPH` file.

:program:`MOTRA` --- Basic and Most Common Keywords
---------------------------------------------------

.. class:: keywordlist

:kword:`FROZEN`
  By symmetry: non-correlated orbitals (default: core)

:kword:`RFPErt`
  Previous reaction field introduced as a perturbation

:kword:`LUMORB`
  Input orbital file as ASCII (INPORB)

:kword:`JOBIPH`
  Input orbital file as binary (JOBOLD)

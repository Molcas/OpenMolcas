.. index::
   single: Program; MRCI
   single: MRCI
   single: CI
   single: ACPF

.. _TUT\:sec\:mrci:

:program:`MRCI` --- A Configuration Interaction Program
=======================================================

Multi Reference Single and Doubles Configuration Interaction (MR-SDCI)
wave functions are produced by the :program:`MRCI` program module in
the |molcas| codes.
The :kword:`SDCI` keyword requests an
ordinary Multi Reference Single and Doubles Configuration Interaction
calculation. This is the default and is mutually exclusive with the
:kword:`ACPF` keyword which requests an Average Coupled Pair Function
calculation. The final keyword, :kword:`ROOT`, specifies the number
of the CI root the calculation should compute. The second CI root is
the first excited state and since the :program:`GUGA` module has computed the
coupling coefficients for a triplet state, the :program:`MRCI` module will
converge to the first excited triplet state.

:program:`mrci` Output
----------------------

The :program:`mrci` section of the output lists the number of each type
of orbital in each symmetry including pre-frozen orbitals that were
frozen by the :program:`guga` module. There is a list of the
reference configurations with the inactive orbitals included. An empty
orbital is listed as "``0``" and a doubly occupied as "``3``". The
spin of a singly occupied orbital by "``1``" (spin up) or "``2``"
(spin down). The total
number of configuration state functions (CSFs) is listed below the reference
configurations.

Sample input requesting the the MRCI module to calculate the first
excited MRCI energy for neutral triplet water in :math:`C_{2v}` symmetry with six
electrons in the active space: ::

  &MRCI
  Title= MR-SDCI of 2nd CI root of C2v Water
  SDCI; Root= 2

A listing of the possible CI roots is followed by the CI iteration and
convergence information. The Davidson and ACPF corrections are included
along with the important CSFs in the CI wave function. The molecular
orbitals are listed near the end of the output.

There are four input files to the :program:`MRCI` module; :file:`CIGUGA`
from :program:`GUGA`, :file:`TRAONE` and :file:`TRAINT` from
:program:`MOTRA` and :file:`ONEINT` from :program:`SEWARD`. The orbitals
are saved in :file:`CIORBnn` where :file:`nn` is the number of the CI root.

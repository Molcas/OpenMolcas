.. index::
   single: Program; CASPT2
   single: CASPT2

.. _TUT\:sec\:caspt2:

:program:`CASPT2` --- A Many Body Perturbation Program
======================================================

Dynamic correlation energy of a molecular system can be calculated using
the :program:`CASPT2` program module in |molcas|. A :program:`CASPT2`
calculation gives a second order perturbation estimate of the full CI energy
using the CASSCF wave function of the system.
The program can also perform Multi-State CASPT2 calculations (MS-CASPT2) in
which different CASPT2 states are coupled using an effective Hamiltonian
computed to second order in perturbation theory. This is necessary in cases
where different CASSCF wave functions are strongly dependent on dynamical
correlation effects. The wave function have to be obtained in a previous
State-Average CASSCF calculation.

.. index::
   single: CASPT2; Input

A sample input is given in :numref:`block:caspt2_input`. The
:kword:`FROZen` keyword specifies the number of orbitals of each
symmetry which will not be included in the correlation. We have
chosen the :program:`RASSCF` :kword:`INACtive` orbitals to be frozen
for this calculation (the default is to freeze all core orbitals, so the input
is strictly not needed). The remaining two keywords, :kword:`CONVergence` and
:kword:`MAXIter`, are included with there default values. The
:kword:`MULTistate` keyword is included for clarity even if not needed in this single
state calculation. A single line follows indicating the number of
simultaneously treated CASPT2 roots and the number of the roots in the previous
SA-CASSCF calculation.

.. index::
   single: CASPT2; Output

:program:`CASPT2` Output
------------------------

In :numref:`TUT:sec:pt2out` the meaning and significance of most of the
features used and printed by the :program:`CASPT2` program are explained in the
context of an actual example. We suggest a careful reading of that section
because understanding the results of a CASPT2 calculation is important for
the analysis of problems like intruder states, large coefficients, convergence,
etc.

.. code-block:: none
   :caption: Sample input requesting the :program:`CASPT2` module to calculate the CASPT2
             energy of a water molecule in :math:`C_{2v}` symmetry with one frozen orbital.
   :name: block:caspt2_input

   &CASPT2
   Frozen= 1 0 0 0
   Multistate= 1 1
   MaxIter= 40

The output of the :program:`CASPT2` program begins with the title
from the input as well as the title from the :program:`SEWARD` input.
It also contains the cartesian coordinates of the molecule and the
CASSCF wave function and orbital specifications. This is followed by
details about the type of Fock and :math:`H_0` operator used and, eventually,
the value of the level-shift parameter employed. It is possible then
to obtain, by input specifications, the quasi-canonical orbitals in
which the wave function will be represented. The following CI vector
and occupation number analysis will be performed using the
quasi-canonical orbitals.

Two important sections follow. First a detailed report on small energy
denominators, large components, and large energy contributions which will
inform about the reliability of the calculation
(see :numref:`TUT:sec:pt2out`)
and finally the :program:`CASPT2` property section
including the natural orbitals obtained
as defined in the output and a number of approximated molecular properties.

If the :kword:`MULTistate` option is used, the program will perform one CASPT2
calculation for each one of the selected roots, and finally the complete
effective Hamiltonian containing the selected states will be solved to obtain
the final MS-CASPT2 energies and PM-CASSCF wave functions :cite:`Finley:98b`.

The :program:`CASPT2` module needs the integral files in :file:`$WorkDir` and the
:file:`RUNFILE` file from the and the :file:`JOBIPH` file from the
:program:`RASSCF` module. The orbitals are saved in the :file:`PT2ORB` file.
The new PM-CASSCF wave functions generated in a MS-CASPT2 calculation
is saved in the :file:`JOBMIX` file.

:program:`CASPT2` --- Basic and Most Common Keywords
----------------------------------------------------

.. class:: keywordlist

:kword:`MULTistate`
  Multi-State CASPT2 calculation: number of roots and roots (Ex. 3 1 2 3)

:kword:`IMAG`
  Value for the imaginary shift for the zero order Hamiltonian

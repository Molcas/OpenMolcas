.. index::
   single: Program; CCSDT
   single: Program; CCSD
   single: Program; CCT3
   single: Program; CCSORT
   single: CCSD(T)
   single: CCSD

.. _TUT\:sec\:ccsdt:

:program:`CCSDT` --- A Set of Coupled-Cluster Programs
======================================================

The |molcas| program :program:`CCSDT`
computes Coupled-Cluster Singles Doubles, CCSD, and Coupled-Cluster Singles
Doubles and Non-iterative Triples Correction CCSD(T) wave functions
for restricted single reference
both closed- and open-shell systems.

In addition to the :file:`ONEINT` and :file:`ORDINT` integral files
(in non-Cholesky calculations),
the :program:`CCSDT` code requires the :file:`JOBIPH` file containing the
reference wave function (remember that it is not possible to
compute open-shell systems with the :program:`SCF` program) and
the transformed two-electron integrals produced by the :program:`MOTRA`
module and stored in the :file:`TRAINT` file.


Previously to execute the :program:`CCSDT` module, wave functions
and integrals have to be prepared. First, a RASSCF calculation has
to be run in such a way that the resulting wave function has one
single reference. In closed-shell situations this means to include
all the orbitals as inactive and set the number of active electrons to zero.
Keyword :kword:`OUTOrbitals` followed by the specification :kword:`CANOnical`
must be used in
the :program:`RASSCF` input to activate the construction of canonical
orbitals and the calculation of the CI-vectors on the basis of
the canonical orbitals.
After that the :program:`MOTRA` module has to
be run to transform the two-electron integrals using the molecular
orbitals provided by the :program:`RASSCF` module.
The files :file:`JOBIPH` or :file:`RASORB` from the
:program:`RASSCF` calculation can be used directly by :program:`MOTRA`
using the keywords :kword:`JOBIph` or :kword:`LUMOrb` in the :program:`MOTRA` input.
Frozen or
deleted orbitals can be introduced in the transformation step
by the proper options in the :program:`MOTRA` input.

:program:`CCSDT` Outputs
------------------------

The section of the |molcas| output corresponding to the CC program
is self explanatory. The default output simply contains
the wave function specifications from the previous RASSCF calculation,
the orbital specifications, the diagonal Fock matrix elements and orbital
energies, the technical description of the calculation, the iterations leading to the CCSD energy,
and the five largest amplitudes of each type, which will help to evaluate
the calculation. If triples excitations have been required the description
of the employed method (from the three available) to compute perturbatively
the triple excited contributions to the CC energy, the value of the
correction, and the energy decomposition into spin parts will be available.

Example of a CCSD(T) calculation
--------------------------------

:numref:`block:ccsdt_input` contains the input files required by the
:program:`seward`, :program:`scf`,
:program:`rasscf`, :program:`motra` and :program:`ccsdt`
programs to compute the ground state of the :math:`\ce{HF^+}` cation.
molecule, which is a doublet of :math:`\Sigma^+` symmetry. A more detailed
description of the different options included in the input of the
programs can be found in the CCSDT section of the user's guide.
This example describes how to calculate CCSD(T) energy for :math:`\ce{HF^+}` cation.
This cation can be safely represented by the single determinant as a reference
function, so one can assume that CCSD(T) method will be suitable for its
description.

The calculation can be divided into few steps:

#. Run :program:`SEWARD` to generate AO integrals.

#. Calculate the HF molecule at the one electron level using :program:`SCF` to
   prepare an estimate of MO for the :program:`RASSCF` run.

#. Calculate :math:`\ce{HF^+}` cation by subtracting one electron from the orbital with
   the first symmetry. There is only one electron in one active orbital
   so only one configuration is created. Hence, we obtain a simple single
   determinant ROHF reference.

#. Perform MO transformation exploiting :program:`MOTRA` using MO coefficients
   from the :program:`RASSCF` run.

#. Perform the Coupled Cluster calculation using :program:`CCSDT` program. First,
   the data produced by the programs :program:`RASSCF` and :program:`MOTRA` need
   to be reorganized, then the CCSD calculation follows, with the chosen spin
   adaptation being T2 DDVV. Finally, the noniterative triple excitation contribution
   calculation is following, where the CCSD amplitudes are used.

This is an open shell case, so it is suitable to choose CCSD(T) method
as it is defined by Watts *et al.* :cite:`t3_watts`.
Since CCSD amplitudes produced by previous :program:`CCSD` run are partly
spin adapted and denominators are produced from the corresponding diagonal
Fock matrix elements,
final energy is sometimes referred as SA1 CCSD(T)\ :math:`_d` (see
:cite:`t3_neo`).

.. A suitable shell script to run these calculations can be found at the end of
   section :ref:`UG:sec:cct3` of the user's guide.

.. extractcode-block:: none
   :filename: tutorials/CCSDT.HF.input
   :caption: Sample input containing the files required by the :program:`SEWARD`, :program:`SCF`,
             :program:`RASSCF`, :program:`MOTRA`, :program:`CCSORT`, :program:`CCSD`, and
             :program:`CCT3` programs to compute the ground state of the :math:`\ce{HF^+}` cation.
   :name: block:ccsdt_input

   &SEWARD &END
   Title= HF molecule
   Symmetry
   X Y
   Basis set
   F.ANO-S-VDZ
   F      0.00000   0.00000   1.73300
   End of basis
   Basis set
   H.ANO-S-VDZ
   H      0.00000   0.00000   0.00000
   End of basis
   End of input
   &SCF
   &RASSCF
   Title= HF(+) cation
   OUTOrbitals= Canonical
   Symmetry= 1; Spin= 2
   nActEl= 1 0 0; Inactive= 2 1 1 0; Ras2= 1 0 0 0
   LumOrb; OUTOrbitals= Canonical
   &MOTRA; JobIph; Frozen= 1 0 0 0
   &CCSDT
   Iterations= 50; Shift= 0.2,0.2; Accuracy= 1.0d-7
   Denominators= 2; Extrapolation= 5,4
   Adaptation= 1; Triples= 3; T3Denominators= 0

:program:`RASSCF` calculates the HF ionized state by removing one electron
from the orbital in the first symmetry.
Do not forget to use keyword
:kword:`CANONICAL`.
In the :program:`CCSDT` run, the number of iterations is limited to 50.
Denominators will be formed using orbital energies. (This corresponds to the
chosen spin adaptation.) Orbitals will be shifted by 0.2 au,
what will accelerate the convergence. However, final energy will not be
affected by the chosen type of denominators and orbital shifts. Required
accuracy is 10\ :math:`^{-7}` au. for the energy. T2 DDVV class of CCSD amplitudes will
be spin adapted.
To accelerate the convergence,
DIIS procedure is exploited. It will start after 5th iteration and
the last four iterations will be taken into account in each extrapolation step.

In the triples step the CCSD(T) procedure as defined
by Watts *et al.* :cite:`t3_watts` will be performed.
Corresponding denominators will be produced using diagonal Fock matrix elements.

:program:`CCSDT` --- Basic and Most Common Keywords
---------------------------------------------------

.. class:: keywordlist

:kword:`CCSD`
  Coupled-cluster singles and doubles method

:kword:`CCT`
  CCSD plus a non iterative triples (T) calculation

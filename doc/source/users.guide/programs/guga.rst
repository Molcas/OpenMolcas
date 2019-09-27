.. index::
   single: Program; GUGA
   single: GUGA

.. _UG\:sec\:guga:

:program:`Guga`
===============

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. xmldoc:: <MODULE NAME="GUGA">
            %%Description:
            <HELP>
            The GUGA program generates coupling coefficients
            used by the MRCI and the CPF programs. The program was written
            by P. E. M. Siegbahn, and has since been slightly modified to fit MOLCAS.
            </HELP>

The :program:`GUGA` program generates :index:`coupling coefficients <single: Coupling coefficients; GUGA>`
used in the :program:`MRCI` and the :program:`CPF` programs
in :index:`Direct CI` calculations :cite:`Roos:72`.
These coupling coefficients are evaluated by the :index:`Graphical Unitary
Group Approach` :cite:`Shavitt:77,Shavitt:78,Siegbahn:80`,
for wavefunctions with at most two electrons excited from a set of
reference configurations. The program was written by :index:`P. E. M. Siegbahn <single: Siegbahn, P. E. M>`,
Institute of Physics, Stockholm University, Sweden.
Only the :program:`MRCI` program can use several reference
configurations. The reference configurations can be specified as a
list, where the occupation numbers are given for each active orbital
(see below) in each reference configuration, or as a :index:`Full CI`
the space defined by the active orbitals. In the :program:`GUGA`, :program:`MRCI`
and :program:`CPF` programs, the orbitals are classified as follows:
Frozen, Inactive, Active, Secondary, and Deleted orbitals. Within each
symmetry type, they follow this order. For the :program:`GUGA` program,
only the inactive and active orbitals are relevant.

* **Inactive:** :index:`Inactive orbitals <single: GUGA; Inactive>` are doubly occupied
  in all reference configurations, but excitations out of this orbital
  space are allowed in the final CI wavefunction, i.e., they are
  correlated but have two electrons in all *reference* configurations.
  Since only single and double excitations are allowed, there can be no
  more than two holes in the active orbitals.
  Using keyword :index:`NoCorr <single: GUGA; NoCorr>` (See input description) a subset of the
  inactive orbitals can be selected, and at most a single hole
  is then allowed in the selected set. This allows the core-polarization
  part of core-valence correlation, while preventing large but usually
  inaccurate double-excitation core correlation.

* **Active:** :index:`Active orbitals <single: GUGA; Active>` are those which may have
  different occupation in different reference configurations.
  Using keyword :index:`OneOcc <single: GUGA; OneOcc>` (See input description) a restriction may be
  imposed on some selection of active orbitals, so that the selected
  orbitals are always singly occupied. This may be useful for transition
  metal compounds or for deep inner holes.

.. index::
   pair: Dependencies; GUGA

.. _UG\:sec\:guga_dependencies:

Dependencies
------------

The :program:`GUGA` program does not depend on any other program for its
execution.

.. index::
   pair: Files; GUGA

.. _UG\:sec\:guga_files:

Files
-----

Input files
...........

The :program:`GUGA` program does not need any input files apart from the file of
input keywords.

Output files
............

.. class:: filelist

:file:`CIGUGA`
  This file contains the coupling coefficients that are needed in
  subsequent CI calculations. For information about how these
  coefficients are structured you are referred to the source
  code :cite:`Siegbahn:80`. The theoretical background for the
  coefficient can be found in Refs :cite:`Shavitt:77,Shavitt:78,Siegbahn:80` and
  references therein.

.. index::
   pair: Input; GUGA

.. _UG\:sec\:guga_input:

Input
-----

This section describes the input to the
:program:`GUGA` program in the |molcas| program system, with the program name: ::

  &GUGA

.. index::
   pair: Keywords; GUGA

Keywords
........

Formally, there are no compulsory keyword. Obviously, some
input must be given for a meaningful calculation.

.. class:: keywordlist

:kword:`TITLe`
  The line following this keyword is treated as title line

  .. xmldoc:: <KEYWORD MODULE="GUGA" NAME="TITLE" APPEAR="Title" KIND="STRING" LEVEL="BASIC">
              %%Keyword: Title <basic>
              <HELP>
              The line following this keyword is treated as title line
              </HELP>
              </KEYWORD>

:kword:`SPIN`
  The spin degeneracy number, i.e. 2S+1. The value is read from the
  line following the keyword, in free format. The default value is
  1, meaning a singlet wave function.

  .. xmldoc:: <KEYWORD MODULE="GUGA" NAME="SPIN" APPEAR="Spin (2S+1)" KIND="INT" LEVEL="BASIC">
              %%Keyword: Spin <basic>
              <HELP>
              Spin degeneracy number (multiplicity), 2S+1. Default 1=Singlet.
              </HELP>
              </KEYWORD>

:kword:`ELECtrons`
  The number of electrons to be correlated in the CI of CPF calculation.
  The value is read from the line following the keyword, in free format.
  Note that this number should include the nr of electrons in inactive
  orbitals. An alternative input specification is NACTEL.
  Default: Twice nr of inactive orbitals.

  .. xmldoc:: <KEYWORD MODULE="GUGA" NAME="ELECTRONS" APPEAR="Nr of electrons." KIND="INT" LEVEL="BASIC">
              %%Keyword: Electrons <basic>
              <HELP>
              Number of electrons to be correlated.
              </HELP>
              </KEYWORD>

:kword:`NACTel`
  The number of electrons in active orbitals in the reference configurations.
  The value is read from the line following the keyword, in free format.
  Note that this number includes only the of electrons in active
  orbitals. An alternative input specification is ELECTRONS.
  Default: Zero.

  .. xmldoc:: <KEYWORD MODULE="GUGA" NAME="NACTEL" APPEAR="Number of active electrons." KIND="INT" LEVEL="BASIC">
              %%Keyword: NACTEL <basic>
              <HELP>
              Number of active electrons in the reference CI (if multireference).
              </HELP>
              </KEYWORD>

:kword:`INACtive`
  The number of inactive orbitals, i.e. orbitals that have
  occupation numbers of 2 in all reference configurations. Specified for
  each of the symmetries. The values are read from the line
  following the keyword, in free format.

  .. xmldoc:: <GROUP MODULE="GUGA" NAME="ORBITALS" APPEAR="Orbitals" KIND="BOX">
              <HELP>
              Various orbital spaces.
              </HELP>

  .. xmldoc:: <KEYWORD MODULE="GUGA" NAME="INACTIVE" APPEAR="Inactive orbitals" KIND="INTS_LOOKUP" SIZE="NSYM" LEVEL="BASIC">
              <HELP>
              Number of inactive orbitals for each irrep.
              </HELP>
              %%Keyword: Inactive <basic>
              List which tells, for each symmetry species, how many orbitals
              to keep fully occupied always. Default is 0 in all symmetries.
              </KEYWORD>

:kword:`ACTIve`
  The number of active orbitals, i.e. orbitals that have varying
  occupation numbers in the reference configurations. Specified for each
  of the symmetries. The values are read from the line following
  the keyword, in free format.

  At least one of the :kword:`Inactive` or :kword:`Active` keywords must
  be present for a meaningful calculation. If one of them is left out,
  the default is 0 in all symmetries.

  .. xmldoc:: <KEYWORD MODULE="GUGA" NAME="ACTIVE" APPEAR="Active orbitals" KIND="INTS_LOOKUP" SIZE="NSYM" LEVEL="BASIC">
              <HELP>
              Number of active orbitals for each irrep.
              </HELP>
              %%Keyword: Active <basic>
              List which tells, for each symmetry species, how many orbitals
              that are active. Default is 0 in all symmetries.
              </KEYWORD>

:kword:`ONEOcc`
  Specify a number of active orbitals per symmetry that are required to have occupation
  number one in all configurations. These orbitals are the first active orbitals.
  The input is read from the line after the keyword, in free format.

  .. xmldoc:: <KEYWORD MODULE="GUGA" NAME="ONEOCC" APPEAR="Singly occupied orbitals" KIND="INTS_LOOKUP" SIZE="NSYM" LEVEL="ADVANCED">
              <HELP>
              Number of always open active orbitals.
              </HELP>
              %%Keyword: OneOcc <advanced>
              List which tells, for each symmetry species, how many orbitals
              that are required to be singly occupied always. Default is 0 in all symmetries.
              </KEYWORD>

:kword:`NOCOrr`
  Specify the number of inactive orbitals per symmetry out of which at most one electron
  (total) is excited. These orbitals are the first inactive orbitals.
  The input is read from the line after the keyword, in free format.

  .. xmldoc:: <KEYWORD MODULE="GUGA" NAME="NOCORR" APPEAR="Always non-empty orbitals" KIND="INTS_LOOKUP" SIZE="NSYM" LEVEL="ADVANCED">
              <HELP>
              Number of always non-empty active orbitals.
              </HELP>
              %%Keyword: NoCorr <advanced>
              List which tells, for each symmetry species, how many orbitals
              that are not allowed to be empty. Default is 0 in all symmetries.
              </KEYWORD>

  .. xmldoc:: </GROUP>

:kword:`REFErence`
  Specify selected reference configurations. The additional input
  that is required usually spans more than one line. The first line
  after the keyword contains the number of reference configurations, and
  the total number of active orbitals, and these two numbers are
  read by free format. Thereafter the input has one line per
  reference configuration, specifying the occupation number for each of
  the active orbitals, read by 80I1 format. Note that
  :kword:`Reference` and :kword:`CIall` are mutually exclusive.

  .. xmldoc:: <GROUP MODULE="GUGA" NAME="REF_SPACE" APPEAR="Reference space" KIND="BOX" LEVEL="BASIC">

  .. xmldoc:: <KEYWORD MODULE="GUGA" NAME="REFERENCE" APPEAR="Reference occupations" KIND="CUSTOM" LEVEL="BASIC">
              <HELP>
              A single string like '22010' for occupations.
              </HELP>
              %%Keyword: Reference <basic>
              One way of specifying the reference space -- see manual.
              One of the two keywords REFERENCE and CIALL should be chosen.
              </KEYWORD>

:kword:`CIALl`
  Use a Full CI within the subspace of the active orbitals as
  reference configurations. The symmetry of the wavefunction must be
  specified. The value is read from the line following the keyword, in
  free format. Note that
  :kword:`CIall` and :kword:`Reference` are mutually exclusive.
  One of these two alternatives must be chosen for a meaningful calculation.

  .. xmldoc:: <KEYWORD MODULE="GUGA" NAME="CIALL" APPEAR="Full reference" KIND="INT" LEVEL="BASIC">
              <HELP>
              Use a full reference.
              </HELP>
              %%Keyword: CIAll <basic>
              Use a Full CI space as reference -- see manual.
              One of the two keywords REFERENCE and CIALL should be chosen.
              </KEYWORD>

  .. xmldoc:: </GROUP>

:kword:`FIRSt`
  Perform a first order calculation, i.e. only single excitations
  from the reference space. No additional input is required.

  .. xmldoc:: <KEYWORD MODULE="GUGA" NAME="FIRST" APPEAR="First-order" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: First-order <advanced>
              <HELP>
              Allow only single excitations from the reference space.
              </HELP>
              </KEYWORD>

:kword:`NONInteracting space`
  By default, those double excitations from inactive
  to virtual orbitals are excluded, where the inactive and virtual electrons
  would couple to a resulting triplet.
  With the NonInteracting Space option, such 'non-interacting' configurations
  are included as well.

  .. xmldoc:: <KEYWORD MODULE="GUGA" NAME="NONINTERACT" APPEAR="Non-interacting space" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: NonInteracting <advanced>
              <HELP>
              Include triplet-coupled double excitations from inactive to virtual orbitals.
              </HELP>
              </KEYWORD>

  .. xmldoc:: <KEYWORD MODULE="GUGA" NAME="INTERACT" KIND="SINGLE" LEVEL="UNDOCUMENTED" />

:kword:`PRINt`
  Printlevel of the program. Default printlevel (0) produces very
  little output. Printlevel 5 gives some information that may be of
  interest. The value is read from the line following the keyword, in free
  format.

  .. xmldoc:: <KEYWORD MODULE="GUGA" NAME="PRINT" APPEAR="Print level" KIND="INT" LEVEL="ADVANCED">
              <HELP>
              Enter print level, from 0 (default) up to 5.
              </HELP>
              </KEYWORD>
              %%Keyword: PrintLevel <advanced>
              Requested print level. Default 0. 5 is reasonable.

Input example
.............

::

  &GUGA
  Title
   Water molecule. 2OH correlated.
  Electrons =     4
  Spin      =     1
  Active    =     2    2    0    0
  Interacting space
  Reference
      3    4
    2020 ; 0220 ; 2002

.. xmldoc:: </MODULE>

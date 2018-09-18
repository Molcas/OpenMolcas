.. index::
   single: Program; GUGADRT
   single: GUGADRT

.. _UG\:sec\:gugadrt:

:program:`GUGADRT`
==================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. xmldoc:: <MODULE NAME="GUGADRT">
            %%Description:
            <HELP>
            The GUGADRT program generates distict row table
            used by the MRCI. The program was written
            by Yubin Wang and Bingbing Suo,
            and has since been slightly modified to fit MOLCAS.
            </HELP>

The :program:`GUGADRT` program generates distict row table (DRT)
used in the :program:`GUGACI`
in :index:`Direct CI` calculations :cite:`Roos:72`.
Only DRT in active space are generated because
the hole-particle symmetry is used in :program:`GUGACI` :cite:`YBWang:1,BSuo:1`.
These DRT are used to evaluated the coupling coefficients
by the :index:`Graphical Unitary
Group Approach` :cite:`Shavitt:77,Shavitt:78,Siegbahn:80`,
for wavefunctions with at most two electrons excited from a set of
reference configurations. The reference configurations can be specified as a
list, where the occupation numbers are given for each active orbital
(see below) in each reference configuration, or as a :index:`Full CI`
within
the space defined by the active orbitals. In the :program:`GUGADRT` and :program:`GUGACI`
the orbitals are classified as follows:
Frozen, Inactive, Active, Secondary, and Deleted orbitals. Within each
symmetry type, they follow this order. For the :program:`GUGADRT` program,
only the active orbitals are relevant.

* **Inactive:** :index:`Inactive orbitals <single: GUGADRT; Inactive>` are doubly occupied
  in all reference configurations, but excitations out of this orbital
  space are allowed in the final CI wavefunction, i.e., they are
  correlated but have two electrons in all *reference* configurations.
  Since only single and double excitations are allowed, there can be no
  more than two holes in the active orbitals.

* **Active:** :index:`Active orbitals <single: GUGADRT; Active>` are those which may have
  different occupation in different reference configurations.

.. index::
   pair: Dependencies; GUGADRT

.. _UG\:sec\:gugadrt_dependencies:

Dependencies
------------

.. index::
   pair: Files; GUGADRT

.. _UG\:sec\:gugadrt_files:

Files
-----

Input files
...........

.. class:: filelist

:file:`TRAONE`
  Transformed one-electron integrals from :program:`MOTRA`. Orbital information
  such as frozen, deleted orbitals will be read from this file.

Output files
............

.. class:: filelist

:file:`GUGADRT`
  This file contains the DRT that is needed in
  subsequent CI calculations.

.. index::
   pair: Input; GUGADRT

.. _UG\:sec\:gugadrt_input:

Input
-----

This section describes the input to the
:program:`GUGADRT` program in the |molcas| program system, with the program name: ::

  &GUGADRT

The first four characters of the keywords are
decoded and the rest are ignored.

.. index::
   pair: Keywords; GUGADRT

Keywords
........

Formally, there are no compulsory keyword. Obviously, some
input must be given for a meaningful calculation.

.. class:: keywordlist

:kword:`TITLe`
  The lines following this keyword are treated as title lines, until
  another keyword is encountered.

  .. xmldoc:: <KEYWORD MODULE="GUGADRT" NAME="TITLE" APPEAR="Title" KIND="STRING" LEVEL="BASIC" >
              %%Keyword: Title <basic>
              <HELP>
              The lines following this keyword are treated as title lines, until
              another keyword is encountered.
              Enter at most ten lines of title.
              </HELP>
              </KEYWORD>

:kword:`SPIN`
  The spin degeneracy number, i.e. 2S+1. The value is read from the
  line following the keyword, in free format. The default value is
  1, meaning a singlet wave function.

  .. xmldoc:: <KEYWORD MODULE="GUGADRT" NAME="SPIN" APPEAR="Spin (2S+1)" KIND="INT" LEVEL="BASIC" >
              %%Keyword: SPIN <basic>
              <HELP>
              Enter spin multiplicity, 2S+1. Default 1=Singlet.
              </HELP>
              </KEYWORD>

:kword:`ELECtrons`
  The number of electrons to be correlated in the CI calculation.
  The value is read from the line following the keyword, in free format.
  Note that this number should include the nr of electrons in inactive
  orbitals. An alternative input specification is NACTEL.
  Default: Twice nr of inactive orbitals.

  .. xmldoc:: <KEYWORD MODULE="GUGADRT" NAME="ELECTRONS" APPEAR="Number of electrons" KIND="INT" LEVEL="BASIC" >
              %%Keyword: Electrons <basic>
              <HELP>
              Enter number of electrons to be correlated.
              </HELP>
              </KEYWORD>

:kword:`NACTel`
  The number of electrons in active orbitals in the reference configurations.
  The value is read from the line following the keyword, in free format.
  Note that this number includes only the of electrons in active
  orbitals. An alternative input specification is ELECTRONS.
  Default: Zero.

  .. xmldoc:: <KEYWORD MODULE="GUGADRT" NAME="NACTEL" APPEAR="Number of active electrons" KIND="INT" LEVEL="BASIC" >
              %%Keyword: NACTEL <basic>
              <HELP>
              Number of active electrons (if multireference).
              </HELP>
              </KEYWORD>

:kword:`INACtive`
  The number of inactive orbitals, i.e. orbitals that have
  occupation numbers of 2 in all reference configurations. Specified for
  each of the symmetries. The values are read from the line
  following the keyword, in free format.

  .. xmldoc:: <KEYWORD MODULE="GUGADRT" NAME="INACTIVE" APPEAR="Inactive orbitals" KIND="INTS_LOOKUP" SIZE="NSYM" LEVEL="BASIC" >
              %%Keyword: Inactive <basic>
              <HELP>
              Number of inactive orbitals for each irrep.
              Enter list which tells, for each symmetry species, how many orbitals
              to keep fully occupied always. Default is 0 in all symmetries.
              </HELP>
              </KEYWORD>

:kword:`ACTIve`
  The number of active orbitals, i.e. orbitals that have varying
  occupation numbers in the reference configurations. Specified for each
  of the symmetries. The values are read from the line following
  the keyword, in free format.

  At least one of the :kword:`Inactive` or :kword:`Active` keywords must
  be present for a meaningful calculation. If one of them is left out,
  the default is 0 in all symmetries.

  .. xmldoc:: <KEYWORD MODULE="GUGADRT" NAME="ACTIVE" APPEAR="Active orbitals" KIND="INTS_LOOKUP" SIZE="NSYM" LEVEL="BASIC" >
              %%Keyword: Active <basic>
              <HELP>
              Number of active orbitals for each irrep.
              Enter list which tells, for each symmetry species, how many orbitals
              that are active. Default is 0 in all symmetries.
              </HELP>
              </KEYWORD>

:kword:`REFErence`
  Specify selected reference configurations. The additional input
  that is required usually spans more than one line. The first line
  after the keyword contains the number of reference configurations, and
  the total number of active orbitals, and these two numbers are
  read by free format. Thereafter the input has one line per
  reference configuration, specifying the occupation number for each of
  the active orbitals, read by 80I1 format. Note that
  :kword:`Reference` and :kword:`CIall` are mutually exclusive.

  .. xmldoc:: <KEYWORD MODULE="GUGADRT" NAME="REFERENCE" APPEAR="Reference occupations" KIND="STRINGS" LEVEL="BASIC" >
              %%Keyword: REFERENCE <basic>
              <HELP>
              One way of specifying the reference space -- see manual.
              </HELP>
              One of the two keywords REFERENCE and CIALL should be chosen.
              </KEYWORD>

:kword:`SYMMetry`
  Specify the selected symmetry type (the irrep) of the wave function
  as a number between 1 and 8 (see SEWARD). Default is 1, which
  always denote the totally symmetric irrep.

  .. xmldoc:: <KEYWORD MODULE="GUGADRT" NAME="CIALL" APPEAR="Full reference" KIND="SINGLE" LEVEL="BASIC" >
              <HELP>
              Use a full reference.
              </HELP>
              %%Keyword: CIAll <basic>
              Use a Full CI space as reference -- see manual.
              One of the two keywords REFERENCE and CIALL should be chosen.
              </KEYWORD>

:kword:`CIALl`
  Use a Full CI within the subspace of the active orbitals as
  reference configurations. The symmetry of the wavefunction must be
  specified. The value is read from the line following the keyword, in
  free format. Note that
  :kword:`CIall` and :kword:`Reference` are mutually exclusive.
  One of these two alternatives must be chosen for a meaningful calculation.

  .. xmldoc:: <KEYWORD MODULE="GUGADRT" NAME="SYMMetry" APPEAR="Symmetry of state" KIND="INT" LEVEL="BASIC" >
              %%Keyword: Symmetry <basic>
              <HELP>
              Specify symmetry of the state to be calculated.
              Default value is 1.
              </HELP>
              </KEYWORD>

:kword:`PRINt`
  Printlevel of the program. Default printlevel (0) produces very
  little output. Printlevel 5 gives some information that may be of
  interest. The value is read from the line following the keyword, in free
  format.

  .. xmldoc:: <KEYWORD MODULE="GUGADRT" NAME="PRINT" APPEAR="Print level" KIND="INT" LEVEL="ADVANCED" >
              %%Keyword: PrintLevel <advanced>
              <HELP>
              Enter print level, from 0 (default) up to 5.
              </HELP>
              </KEYWORD>

Input example
.............

::

  &GUGADRT
  Title     =  CH2 molecule.
  Electrons =  8
  Spin      =  1
  Inactive  =  1    0    0    0
  Active    =  2    2    2    0
  Symmetry  =  1
  Ciall

.. xmldoc:: </MODULE>

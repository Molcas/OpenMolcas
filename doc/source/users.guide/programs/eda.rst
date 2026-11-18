.. index::
   single: Program; EDA
   single: EDA
   single: CMOEDA

.. _UG\:sec\:eda:

:program:`EDA`
==============

.. xmldoc:: <MODULE NAME="EDA" APPEAR="EDA">
            %%Description:
            <HELP>
            Closed-shell canonical molecular orbital energy decomposition
            analysis (CMOEDA).
            </HELP>

Description
-----------

The :program:`EDA` module performs the closed-shell canonical molecular
orbital energy decomposition analysis (CMOEDA), following the
energy-decomposition scheme of Su and Li :cite:`Su2009`.

The current implementation supports closed-shell fragments and the
reference methods AUTO, RHF, and MBPT2.

Input
-----

The input for the :program:`EDA` module begins with

::

  &EDA

The molecular system is divided into fragments by specifying the number
of fragments, the number of atoms in each fragment, their charges and
spin multiplicities, and the electronic-structure reference used for
the decomposition.

Keywords
--------

NFRAG
.....

.. xmldoc:: <KEYWORD MODULE="EDA" NAME="NFRAG" APPEAR="Number of fragments" KIND="INT" LEVEL="BASIC" MIN_VALUE="1">
            %%Keyword: NFRAG <basic>
            <HELP>
            Number of molecular fragments participating in the decomposition.
            </HELP>
            </KEYWORD>

Specifies the number of molecular fragments participating in the
decomposition. The value must be a positive integer.

NATOM
.....

.. xmldoc:: <KEYWORD MODULE="EDA" NAME="NATOM" APPEAR="Atoms per fragment" KIND="INTS_LOOKUP" SIZE="NFRAG" LEVEL="BASIC" REQUIRE="NFRAG">
            %%Keyword: NATOM <basic>
            <HELP>
            Number of atoms belonging to each fragment. One integer must be
            supplied for every fragment.
            </HELP>
            </KEYWORD>

Specifies the number of atoms belonging to each fragment. One integer
must be supplied for every fragment, so the number of values must equal
the value specified by ``NFRAG``.

CHARGE
......

.. xmldoc:: <KEYWORD MODULE="EDA" NAME="CHARGE" APPEAR="Fragment charges" KIND="INTS_LOOKUP" SIZE="NFRAG" LEVEL="BASIC" REQUIRE="NFRAG">
            %%Keyword: CHARGE <basic>
            <HELP>
            Charge of each fragment. One integer must be supplied for every
            fragment.
            </HELP>
            </KEYWORD>

Specifies the charge of each molecular fragment. One integer must be
supplied for every fragment, in the same fragment order used for
``NATOM``.

MULT
....

.. xmldoc:: <KEYWORD MODULE="EDA" NAME="MULT" APPEAR="Fragment multiplicities" KIND="INTS_LOOKUP" SIZE="NFRAG" LEVEL="BASIC" REQUIRE="NFRAG">
            %%Keyword: MULT <basic>
            <HELP>
            Spin multiplicity of each fragment. The current implementation
            requires closed-shell singlet fragments with MULT=1.
            </HELP>
            </KEYWORD>

Specifies the spin multiplicity of each molecular fragment. One value
must be supplied for every fragment. The current implementation supports
only closed-shell singlet fragments, and therefore every multiplicity
must be ``1``.

REFERENCE
.........

.. xmldoc:: <KEYWORD MODULE="EDA" NAME="REFERENCE" APPEAR="Reference method" KIND="CHOICE" LIST="AUTO:Automatic reference selection,RHF:Closed-shell RHF,MBPT2:Closed-shell MBPT2" LEVEL="BASIC">
            %%Keyword: REFERENCE <basic>
            <HELP>
            Selects the electronic-structure reference. Accepted values are
            AUTO, RHF, and MBPT2.
            </HELP>
            </KEYWORD>

Selects the electronic-structure reference used for the energy
decomposition. The accepted values are:

``AUTO``
  Select the reference method automatically.

``RHF``
  Use a closed-shell restricted Hartree--Fock reference.

``MBPT2``
  Use a closed-shell second-order Moller--Plesset reference.

FBASIS
......

.. xmldoc:: <KEYWORD MODULE="EDA" NAME="FBASIS" APPEAR="Fragment basis convention" KIND="CHOICE" LIST="OWN:Own fragment basis,SUPER:Supermolecule basis" LEVEL="BASIC" DEFAULT_VALUE="OWN">
            %%Keyword: FBASIS <basic>
            <HELP>
            Selects the fragment-basis convention. Accepted values are OWN
            and SUPER. OWN is the currently validated mode.
            </HELP>
            </KEYWORD>

Selects the basis convention used for the fragment calculations. The
accepted values are:

``OWN``
  Use the basis functions belonging to each individual fragment. This
  is the currently validated mode and the default.

``SUPER``
  Use the full supermolecule basis for each fragment.

DEBUG
.....

.. xmldoc:: <KEYWORD MODULE="EDA" NAME="DEBUG" APPEAR="Debug output" KIND="SINGLE" LEVEL="ADVANCED">
            %%Keyword: DEBUG <advanced>
            <HELP>
            Print additional diagnostic information.
            </HELP>
            </KEYWORD>

Requests additional diagnostic information from the :program:`EDA`
module. This keyword is intended primarily for debugging and
development.

Example
-------

For a system divided into two closed-shell fragments, an
:program:`EDA` input can, e.g., contain

::

  &EDA
  NFRAG
    2
  NATOM
    2 1
  CHARGE
    0 0
  MULT
    1 1
  REFERENCE
    RHF
  FBASIS
    OWN

In this example the system is divided into two neutral singlet
fragments containing two and one atoms, respectively. The
decomposition uses an RHF reference with each fragment represented in
its own basis.

.. xmldoc:: </MODULE>

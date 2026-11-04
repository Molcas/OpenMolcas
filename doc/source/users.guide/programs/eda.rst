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
energy-decomposition scheme of Su and Li [J. Chem. Phys. 131, 014102
(2009), DOI: 10.1063/1.3159673].

The current implementation supports closed-shell fragments and the
reference methods AUTO, RHF, and MBPT2.

Input
-----

The input for the module begins with

::

  &EDA

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

NATOM
.....

.. xmldoc:: <KEYWORD MODULE="EDA" NAME="NATOM" APPEAR="Atoms per fragment" KIND="INTS_LOOKUP" SIZE="NFRAG" LEVEL="BASIC" REQUIRE="NFRAG">
            %%Keyword: NATOM <basic>
            <HELP>
            Number of atoms belonging to each fragment. One integer must be
            supplied for every fragment.
            </HELP>
            </KEYWORD>

CHARGE
......

.. xmldoc:: <KEYWORD MODULE="EDA" NAME="CHARGE" APPEAR="Fragment charges" KIND="INTS_LOOKUP" SIZE="NFRAG" LEVEL="BASIC" REQUIRE="NFRAG">
            %%Keyword: CHARGE <basic>
            <HELP>
            Charge of each fragment. One integer must be supplied for every
            fragment.
            </HELP>
            </KEYWORD>

MULT
....

.. xmldoc:: <KEYWORD MODULE="EDA" NAME="MULT" APPEAR="Fragment multiplicities" KIND="INTS_LOOKUP" SIZE="NFRAG" LEVEL="BASIC" REQUIRE="NFRAG">
            %%Keyword: MULT <basic>
            <HELP>
            Spin multiplicity of each fragment. The current implementation
            requires closed-shell singlet fragments with MULT=1.
            </HELP>
            </KEYWORD>

REFERENCE
.........

.. xmldoc:: <KEYWORD MODULE="EDA" NAME="REFERENCE" APPEAR="Reference method" KIND="CHOICE" LIST="AUTO:Automatic reference selection,RHF:Closed-shell RHF,MBPT2:Closed-shell MBPT2" LEVEL="BASIC">
            %%Keyword: REFERENCE <basic>
            <HELP>
            Selects the electronic-structure reference. Accepted values are
            AUTO, RHF, and MBPT2.
            </HELP>
            </KEYWORD>

FBASIS
......

.. xmldoc:: <KEYWORD MODULE="EDA" NAME="FBASIS" APPEAR="Fragment basis convention" KIND="CHOICE" LIST="OWN:Own fragment basis,SUPER:Supermolecule basis" LEVEL="BASIC" DEFAULT_VALUE="OWN">
            %%Keyword: FBASIS <basic>
            <HELP>
            Selects the fragment-basis convention. Accepted values are OWN
            and SUPER. OWN is the currently validated mode.
            </HELP>
            </KEYWORD>

DEBUG
.....

.. xmldoc:: <KEYWORD MODULE="EDA" NAME="DEBUG" APPEAR="Debug output" KIND="SINGLE" LEVEL="ADVANCED">
            %%Keyword: DEBUG <advanced>
            <HELP>
            Print additional diagnostic information.
            </HELP>
            </KEYWORD>

.. xmldoc:: </MODULE>


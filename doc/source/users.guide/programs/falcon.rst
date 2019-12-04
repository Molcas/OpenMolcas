.. index::
   single: Program; FALCON
   single: FALCON

.. _UG\:sec\:falcon:

:program:`falcon` |extramark|
=============================

.. warning::

   This program is not available in |openmolcas|

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. _UG\:sec\:falcon_description:

Description
-----------

.. xmldoc:: <MODULE NAME="FALCON">
            %%Description:
            <HELP>
            FALCON calculates total energy and orbitals of the large system
            based on the fragment method.
            </HELP>

.. compound::

  :program:`FALCON` calculates total energy of the large system based on
  the fragment approach.
  Total energy of the whole system is calculated from total energies of
  fragments as follows,

  .. math:: E^{\text{whole}}=\sum C_i^{\text{fragment}} E_i^{\text{fragment}},

  where :math:`E_i^{\text{fragment}}` is the total energy of fragment :math:`i`, and
  :math:`C_i^{\text{fragment}}` is its coefficient.

.. compound::

  In addition to the total energy, :program:`FALCON` can calculate orbitals
  of the whole system.
  Fock matrix and overlap matrix of the whole system are calculated from
  ones of fragments using following equations,

  .. math:: \mat{F}^{\text{whole}}=\sum C_i^{\text{fragment}} \mat{F}_i^{\text{fragment}},

  and

  .. math:: \mat{S}^{\text{whole}}=\sum C_i^{\text{fragment}} \mat{S}_i^{\text{fragment}},

  where
  :math:`F_i` and :math:`S_i` are the Fock matrix and overlap matrix,
  respectively, of fragment :math:`i`.

.. compound::

  Then

  .. math:: \mat{F}\mat{C}=\mat{S}\mat{C}\mat{\varepsilon}

  is solved to obtain the orbitals, :math:`\mat{C}`, and orbitals energies,
  :math:`\mat{\varepsilon}`.

.. index::
   pair: Input; FALCON

.. _UG\:sec\:falcon_inpfalcon:

Input
-----

Below follows a description of the input to :program:`FALCON`.

The input for each module is preceded by its name like: ::

  &FALCON

Argument(s) to a keyword, either individual or composed by several entries,
can be placed in a separated line or in the same line separated by a semicolon.
If in the same line, the first argument requires an equal sign after the
name of the keyword.

Keywords
........

.. class:: keywordlist

:kword:`TITLe`
  One-line title.

  .. xmldoc:: <KEYWORD MODULE="FALCON" NAME="TITLE" KIND="STRING" LEVEL="BASIC">
              %%Keyword: Title <basic>
              <HELP>
              One line title.
              </HELP>
              </KEYWORD>

:kword:`FRAGment`
  Takes one, two or three argument(s).
  The first value (integer) defines the fragment number,
  the second value (real) determines coefficient,
  and the third value (integer) is the fragment number that is equivalent
  to this fragment when translational symmetry is used.
  A default for the second value is 1.0 where the first and third values have
  no default.
  Other keyword(s) specific to this fragment must follow this keyword.

  .. xmldoc:: <KEYWORD MODULE="FALCON" NAME="FRAGMENT" KIND="INT" LEVEL="BASIC">
              %%Keyword: Fragment <basic> GUI:keyword
              <HELP>
              Takes one, two or three argument(s).
              The first value defines the fragment number, the second value determines coefficient,
              and the third value is the fragment number that is equivalent to this fragment
              when translational symmetry is used.
              Other keyword(s) specific to this fragment must follow this keyword.
              </HELP>
              </KEYWORD>

:kword:`OPERator`
  A real value following this keyword represents a coefficient, :math:`C_i^{\text{fragment}}`,
  of fragment :math:`i` (current fragment), where :math:`i` is a value specified by FRAGMENT keyword.
  This keyword is equivalent with the second value of keyword, FRAGMENT.

  .. xmldoc:: <KEYWORD MODULE="FALCON" NAME="OPERATOR" KIND="REAL" LEVEL="BASIC">
              %%Keyword: Operator <basic> GUI:keyword
              <HELP>
              An integer following this keyword specifies the number of fragments.
              </HELP>
              </KEYWORD>

:kword:`EQUIvalence`
  An integer, :math:`j`, following this keyword declares that current fragment
  is translationally equivalent with fragment :math:`j`, and information provided for
  fragment :math:`j` are tranfered to current fragment.
  This keyword is equivalent with the third value of keyword, FRAGMENT.

  .. xmldoc:: <KEYWORD MODULE="FALCON" NAME="EQUIVALENCE" KIND="INT" LEVEL="BASIC">
              %%Keyword: Equivalence <basic> GUI:keyword
              <HELP>
              An integer following this keyword specifies the number of fragments.
              </HELP>
              </KEYWORD>

:kword:`TRANslate`
  Three real numbers following this keyword specifies the translational
  vector by which the current fragment is translated to give new coordinate.
  A unit of either bohr or angstrom can follow. The default unit is angstrom.
  This keyword takes effect only when the equivalent fragment is specified.

  .. xmldoc:: <KEYWORD MODULE="FALCON" NAME="TRANSLATE" KIND="REALS" SIZE="3" LEVEL="BASIC">
              %%Keyword: Translate <basic> GUI:keyword
              <HELP>
              An integer following this keyword specifies the number of fragments.
              </HELP>
              </KEYWORD>

:kword:`RUNFile`
  Following this keyword specifies the name of RunFile file for the
  corresponding fragment.

  .. xmldoc:: <KEYWORD MODULE="FALCON" NAME="RUNFILE" KIND="STRING" LEVEL="BASIC">
              %%Keyword: RunFile <basic> GUI:keyword
              <HELP>
              Following this keyword specifies the name of RunFile file for the
              corresponding fragment.
              </HELP>
              </KEYWORD>

:kword:`ONEInt`
  Following this keyword specifies the name of OneInt file for the
  corresponding fragment.

  .. xmldoc:: <KEYWORD MODULE="FALCON" NAME="ONEINT" KIND="STRING" LEVEL="BASIC">
              %%Keyword: OneInt <basic> GUI:keyword
              <HELP>
              Following this keyword specifies the name of OneInt file for the
              corresponding fragment.
              </HELP>
              </KEYWORD>

:kword:`NFRAgment`
  An integer following this keyword specifies the number of fragments.
  If this keyword is not given, the largest fragment number given by
  FRAGMENT keyword is set to be the number of fragment.

  .. xmldoc:: <KEYWORD MODULE="FALCON" NAME="NFRAGMENT" KIND="INT" LEVEL="BASIC">
              %%Keyword: nFragment <basic> GUI:keyword
              <HELP>
              An integer following this keyword specifies the number of fragments.
              </HELP>
              </KEYWORD>

:kword:`NIRRep`
  An integer following this keyword specifies the number of irreducible
  representation of point group symmetry.

  .. xmldoc:: <KEYWORD MODULE="FALCON" NAME="NIRREP" KIND="INT" LEVEL="BASIC">
              %%Keyword: nIrrep <basic> GUI:keyword
              <HELP>
              An integer following this keyword specifies the number of irreducible
              representation of point group symmetry.
              </HELP>
              </KEYWORD>

:kword:`OCCUpation`
  A list of integer(s) following this keyword specifies the number of
  occupied orbitals in each symmetry representation in the unfragmented
  system.

  .. xmldoc:: <KEYWORD MODULE="FALCON" NAME="OCCUPATION" KIND="INTS_LOOKUP" SIZE="NSYM" LEVEL="BASIC">
              %%Keyword: Occupation <basic> GUI:keyword
              <HELP>
              A list of integer(s) following this keyword specifies the number of
              occupied orbitals in each symmetry representation.
              </HELP>
              </KEYWORD>

:kword:`DISTance`
  A real number following this keyword specifies the distance
  of two atoms that are equivalent to each other, followed by a unit that
  is eather angstrom or bohr.
  Default is angstrom.

  .. xmldoc:: <KEYWORD MODULE="FALCON" NAME="DISTANCE" KIND="REAL" LEVEL="BASIC">
              %%Keyword: Distance <basic> GUI:keyword
              <HELP>
              A real number following this keyword specifies the distance
              of two atoms that are equivalent to each other.
              </HELP>
              </KEYWORD>

:kword:`NEAR`
  A real number following this keyword specifies the distance
  of two atoms within which atoms are considered to be too close each other.
  An unit that is eather angstrom or bohr can follow.
  Default is angstrom.

  .. xmldoc:: <KEYWORD MODULE="FALCON" NAME="NEAR" KIND="REAL" LEVEL="BASIC">
              %%Keyword: Near <basic> GUI:keyword
              <HELP>
              A real number following this keyword specifies the distance
              of two atoms within which atoms are considered to be too close each other.
              </HELP>
              </KEYWORD>

:kword:`PRINt`
  An integer following this keyword specifies the format of orbital print out.

  .. xmldoc:: <KEYWORD MODULE="FALCON" NAME="PRINT" KIND="INT" LEVEL="BASIC">
              %%Keyword: Print <basic> GUI:keyword
              <HELP>
              An integer following this keyword specifies the format of orbital
              print out.
              </HELP>
              </KEYWORD>

:kword:`ORBEne`
  A real number follwing this keyword stands for the threshold for orbital print
  out.
  The orbitals with orbital energy below this value are print out.

  .. xmldoc:: <KEYWORD MODULE="FALCON" NAME="ORBENE" KIND="REAL" LEVEL="BASIC">
              %%Keyword: OrbEne <basic> GUI:keyword
              <HELP>
              An integer following this keyword specifies the format of orbital
              print out.
              </HELP>
              </KEYWORD>

:kword:`ORBOcc`
  A real number follwing this keyword stands for the threshold for orbital print
  out.
  The orbitals with occupation number above this value are print out.

  .. xmldoc:: <KEYWORD MODULE="FALCON" NAME="ORBOCC" KIND="REAL" LEVEL="BASIC">
              %%Keyword: OrbOcc <basic> GUI:keyword
              <HELP>
              An integer following this keyword specifies the format of orbital
              print out.
              </HELP>
              </KEYWORD>

Input examples
..............

.. compound::

  Below shows an example of input file for the three fragment system of which
  energy, :math:`E^{\text{whole}}`, is written as

  .. math:: E^{\text{whole}}= E_1^{\text{fragment}} + E_2^{\text{fragment}} - E_3^{\text{fragment}},

  by fragment energies, :math:`E_1^{\text{fragment}}`, :math:`E_2^{\text{fragment}}`, and :math:`E_3^{\text{fragment}}`. ::

    &FALCON
    Fragment=1,  1.0
    Fragment=2,  1.0
    Fragment=3, -1.0

  which can be simplified as, ::

    &FALCON
    Fragment=3, -1.0

The next example is a two fragment system in which fragment 1 and fragment 2
are equivalent except for their positons.
When their difference in position is described by a vector, (1.0, 1.0, -1.0),
a translational symmetry can be used and the input becomes as follows, ::

  &FALCON
  Fragment=2, 1.0, 1
  Translate=1.0, 1.0, -1.0

.. compound::

  If the total energy of the whole system is given by the sum of total energies
  of three fragment,

  .. math:: E^{\text{whole}}= E_1^{\text{fragment}} + E_2^{\text{fragment}} + E_3^{\text{fragment}},

  input is simplly as follows, ::

    &FALCON
    nFragment=3

.. xmldoc:: </MODULE>

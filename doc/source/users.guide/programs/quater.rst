.. index::
   single: Program; Quater
   single: Quater

.. _UG\:sec\:quater:

:program:`quater`
=================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. xmldoc:: <MODULE NAME="QUATER">
            %%Description:
            <HELP>
            This program aligns two molecules or frames of reference
            </HELP>

.. index::
   pair: Dependencies; Quater

.. _UG\:sec\:quater_dependencies:

Dependencies
------------

The :program:`QUATER` is free-standing and does not depend on any
other program.

.. index::
   pair: Files; Quater

.. _UG\:sec\:quater_files:

Files
-----

Input files
...........

The calculation of vibrational wave functions and spectroscopic
constants uses no input files (except for the standard input).

.. index::
   pair: Input; Quater

.. _UG\:sec\:quater_input:

Input
-----

This section describes the input to the :program:`QUATER` program in the
|molcas| program system. The program name is ::

  &QUATER

.. index::
   pair: Keywords; Quater

Keywords
........

.. class:: keywordlist

:kword:`NOROtation`
  No rotation is performed by the program.
  Only the rotation matrix is printed out.

  .. xmldoc:: <KEYWORD MODULE="QUATER" NAME="NOROTATION" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: NOROtation <basic>
              <HELP>
              No rotation is performed by the program.
              Only the rotation matrix is printed out.
              </HELP>
              </KEYWORD>

:kword:`NOTRanslation`
  No translation is performed by the program.

  .. xmldoc:: <KEYWORD MODULE="QUATER" NAME="NOTRANSLATION" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: NOTRanslation <basic>
              <HELP>
              No translation is performed by the program.
              </HELP>
              </KEYWORD>

:kword:`DEBUg`
  Turn on DEBUG printout

  .. xmldoc:: <KEYWORD MODULE="QUATER" NAME="DEBUG" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: DEBUg <basic>
              <HELP>
              Turn on DEBUG printout
              </HELP>
              </KEYWORD>

:kword:`AXIS`
  Define the old frame of reference

  .. xmldoc:: <KEYWORD MODULE="QUATER" NAME="AXIS" KIND="REALS" SIZE="6" LEVEL="BASIC">
              %%Keyword: AXIS <basic>
              <HELP>
              Define the old frame of reference
              </HELP>
              </KEYWORD>

:kword:`NEWAxis`
  Define the new frame of reference

  .. xmldoc:: <KEYWORD MODULE="QUATER" NAME="NEWAXIS" KIND="REALS" SIZE="6" LEVEL="BASIC">
              %%Keyword: NEWAxis <basic>
              <HELP>
              Define the new frame of reference
              </HELP>
              </KEYWORD>

:kword:`GEO1`
  Define the first geometry

  .. xmldoc:: <KEYWORD MODULE="QUATER" NAME="GEO1" KIND="CUSTOM" LEVEL="BASIC">
              %%Keyword: GEO1 <basic>
              <HELP>
              Define the first geometry
              </HELP>
              </KEYWORD>

:kword:`GEO2`
  Define the second geometry

  .. xmldoc:: <KEYWORD MODULE="QUATER" NAME="GEO2" KIND="CUSTOM" LEVEL="BASIC">
              %%Keyword: GEO2 <basic>
              <HELP>
              Define the second geometry
              </HELP>
              </KEYWORD>

:kword:`XYZ1`
  Define the origin and two axes for the orientation of the first geometry by
  the index of three atoms of this geometry.

  .. xmldoc:: <KEYWORD MODULE="QUATER" NAME="XYZ1" KIND="INTS" SIZE="3" LEVEL="BASIC">
              %%Keyword: XYZ1 <basic>
              <HELP>
              Define the origin and two axes for the orientation of the first geometry by
              the index of three atoms of this geometry.
              </HELP>
              </KEYWORD>

:kword:`XYZ2`
  Define the origin and two axes for the orientation of the second geometry by
  the index of three atoms of this geometry.

  .. xmldoc:: <KEYWORD MODULE="QUATER" NAME="XYZ2" KIND="INTS" SIZE="3" LEVEL="BASIC">
              %%Keyword: XYZ2 <basic>
              <HELP>
              Define the origin and two axes for the orientation of the second geometry by
              the index of three atoms of this geometry.
              </HELP>
              </KEYWORD>

:kword:`END`
  End of input

  .. xmldoc:: <KEYWORD MODULE="QUATER" NAME="END" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: END <basic>
              <HELP>
              End of input
              </HELP>
              </KEYWORD>

:program:`QUATER` will perform a vib-rot analysis and compute
spectroscopic constants.

Input example
.............

::

  &QUATER

  GEO1
     19
  titre
   C     0.000000     0.000000     0.000000
   O     0.000000     0.000000     1.400000
   H     0.895670     0.000000     1.716663
   C    -0.683537    -1.183920    -0.483333
   H    -0.513360     0.889165    -0.363000
   C    -0.683537    -1.183920    -1.933333
   H    -0.170177    -2.073085    -0.120333
   H    -1.710256    -1.183920    -0.120333
   C     0.683537    -1.183920    -2.416667
   H    -1.196896    -2.073085    -2.296333
   H    -1.196896    -0.294755    -2.296333
   C     1.367073     0.000000    -1.933333
   H     1.196896    -2.073085    -2.053667
   H     0.683537    -1.183920    -3.505667
   C     1.367073     0.000000    -0.483333
   H     2.393792     0.000000    -2.296333
   H     0.853714     0.889165    -2.296333
   H     1.880433    -0.889165    -0.120333
   H     1.880433     0.889165    -0.120333
  END
  GEO2
     23
  titre
   C     0.000000     0.000000     0.000000
   H     0.000000     0.000000     1.089000
   C     1.367075     0.000000    -0.483328
   H    -0.334267    -0.970782    -0.363000
   C     1.367081     0.000000    -1.933328
   H     1.880433     0.889165    -0.120326
   H     1.880433    -0.889165    -0.120326
   C     0.683546     1.183920    -2.416664
   H     2.393801     0.000000    -2.296324
   H     0.853722    -0.889165    -2.296330
   C    -0.683529     1.183920    -1.933336
   H     1.196904     2.073085    -2.053662
   O     0.683551     1.183920    -3.816664
   C    -0.683535     1.183920    -0.483336
   H    -1.196887     2.073085    -2.296338
   H    -1.196887     0.294755    -2.296338
   O    -0.023570     2.327015    -0.016667
   H    -1.710255     1.183920    -0.120340
   H     0.237132     1.957142    -4.132332
   C    -0.023576     2.327015     1.383333
   H     0.489783     3.216180     1.746335
   H    -1.050296     2.327015     1.746329
   H     0.489783     1.437850     1.746335
  END
  XYZ1
  15 12 9
  XYZ2
  11 14 1
  END

This input will perform the alignment of the second geometry (GEO2) on the first one (GEO1).
Atom number 11 (C11) of the second geometry will be moved to the position of atom number
15 of the first geometry (C15).
The vector C11 C14 in GEO1 will be aligned with the vector C15 C12 of GEO1.
Finally the plane 11 14 1 of GEO1 will be aligned with the plane 15 12 9 of GEO2.

.. xmldoc:: </MODULE>

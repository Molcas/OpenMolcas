************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine find_aniso_format(old_aniso_format)

      Implicit None
      Logical :: old_aniso_format
      Character(Len=280) :: line
      Integer :: LINENR, Input

      Input=5
      old_aniso_format=.false.

C=========== End of default settings====================================
      REWIND(Input)
50    READ(Input,'(A280)',End=998) LINE
      Call NORMAL(LINE)
      If(LINE(1:11).ne.'&POLY_ANISO') Go To 50
      LINENR=0
100   READ(Input,'(A280)',End=998) line
      LINENR=LINENR+1
      Call NORMAL(LINE)
      If (LINE(1:1).eq.'*') Go To 100
      If (LINE.eq.' ') Go To 100
      If( LINE(1:4).ne.'OLDA') Go To 100
      If((LINE(1:4).eq.'END ').OR.(LINE(1:4).eq.'    ' )) Go To 200

      If (line(1:4).eq.'OLDA') Then

          old_aniso_format=.true.

          LINENR=LINENR+1
          Go To 100
      End If

200   Continue

      Go To 190
C------ errors ------------------------------
998   continue
      Write(6,*)' READIN: Unexpected End of input file.'

190   Continue
      Write(6,*) 'find_aniso_format::  old_aniso_format=',
     &  old_aniso_format
      Return
      End

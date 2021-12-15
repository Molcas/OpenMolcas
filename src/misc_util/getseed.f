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
      Subroutine getSeed(iseed)
      Implicit None
      Integer,Intent(Out) :: iseed
      Integer   i
      External  datimx
      Character Line*72
      Integer*8 hours,minutes,seconds,days
      Character*100, External :: Get_ProgName
*
*     Externally defined seed
      Call getenvf('MOLCAS_RANDOM_SEED',Line)
      If (Line.ne.' ') Then
        Read(Line,*) iseed
        Return
      End If
*
*     Somewhat reproducible if inside verification
      Call getenvf('MOLCAS_TEST',Line)
      If (Line.ne.' ') Then
        Call getenvf('MOLCAS_ITER',Line)
        Read(Line,*) iseed
        Call getenvf('MOLCAS_PRINT',Line)
        Do i=1,Len_Trim(Line)
          iseed = iseed+iChar(Line(i:i))
        End Do
        Line=Trim(Get_ProgName())
        Do i=1,Len_Trim(Line)
          iseed = iseed+iChar(Line(i:i))
        End Do
        Return
      End If
*
*     Default: based on time and project name
      Call datimx(Line)
      Read(Line,'(8x,i2,1x,i2,1x,i2,1x,i2)') days,hours,minutes,seconds
      iseed=Int(((days*24+hours)*60+minutes)*60+seconds,Kind(iseed))
      Call getenvf('Project',Line)
      Do i=1,Len_Trim(Line)
        iseed = iseed+iChar(Line(i:i))
      End Do
*
      Return
      End

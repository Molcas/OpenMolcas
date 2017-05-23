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
      Integer Function nToken(Line)
      Character*(*) Line
      Logical On
*
      nLine=Len(Line)
      nToken = 0
      On = .True.
*
      iLine = 1
 99   Continue
      If (iLine.eq.nLine) Return
      If (Line(iLine:iLine).eq.' ') Then
         On = .True.
      Else
         If (On) nToken = nToken + 1
         On = .False.
      End If
      iLine = iLine + 1
      Go To 99
*
c      Write (6,*) 'Error in nToken!'
c      Call Abend()
      End

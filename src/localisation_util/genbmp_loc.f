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
      SubRoutine GenBMp_Loc(X,nRow,nCol,FilNam,Color)
      Implicit None
      Integer nRow, nCol
      Real*8  X(nRow,nCol)
      Character*(*) FilNam
      Character*1   Color

      Character*10 SecNam
      Parameter (SecNam = 'GenBMp_Loc')

      Integer  isFreeUnit
      External isFreeUnit

      Character*80 Txt
      Integer Lunit, irc, nStp
      Real*8  StpSiz

      Lunit = isFreeUnit(11)
      Call Molcas_Open(Lunit,FilNam)
      irc = 0
      nStp = -1
      StpSiz = -1.0d0
      Call GenBMp(irc,X,nRow,nCol,Lunit,nStp,StpSiz,Color)
      If (irc .ne. 0) Then
         Write(Txt,'(A,I4)') 'GenBMp returned',irc
         Call SysAbendMsg(SecNam,'GenBMp failed!',Txt)
      End If
      Close(Lunit,Status='Keep')

      End

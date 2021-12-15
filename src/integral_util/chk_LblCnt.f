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
      Subroutine Chk_LblCnt(Lbl,nList)
      use Center_Info
      Character*(*) Lbl
      Character(Len=72) Warning
*
      Do iList = 1, nList
         If (Lbl.Eq.dc(iList)%LblCnt) then
            Write (Warning,'(A,A)') 'ChkLbl: Duplicate label;'//
     &                 ' Lbl=',Lbl
            Call WarningMessage(2,Warning)
            Call Quit_OnUserError()
         End If
      End Do
*
      Return
      End

!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      Subroutine FoundAtomicNumber(LuWr,Symb,NAT,iErr)
      Implicit Integer (i-n)
      Character*6 Symb
#include "periodic_table.fh"
      If ((Symb(1:1).le.'z').and.(Symb(1:1).ge.'a'))                    &
     &     Symb(1:1)=char(ichar(Symb(1:1))-32)
      If ((Symb(2:2).le.'Z').and.(Symb(2:2).ge.'A'))                    &
     &     Symb(2:2)=char(ichar(Symb(2:2))+32)
      iErr = 1
! --- Z atoms are ghost atoms that will disappear in seward
      If (Symb(1:1).EQ.'Z' .AND.                                        &
     &    Symb(2:2).NE.'n' .AND. Symb(2:2).NE.'r') then
        iErr = 0
        NAT = -1
        Return
      EndIf
! --- X atoms are dummy atoms that will appear in seward
      If (Symb(1:1).EQ.'X' .AND. Symb(2:2).NE.'e') then
        iErr = 0
        NAT =  0
        Return
      EndIf

      Do i = 1, Num_Elem
         If (Symb(1:2).EQ.AdjustL(PTab(i))) then
          iErr = 0
          NAT = i
          Return
        EndIf
      EndDo
      Do i = 1, Num_Elem
        If (Symb(1:1)//' '.EQ.AdjustL(PTab(i))) then
          iErr = 0
          NAT = i
          Return
        EndIf
      EndDo
      Write(LuWr,*) '   [FoundAtomicNumber]: Wrong atomic symbol !'
      Return
      End

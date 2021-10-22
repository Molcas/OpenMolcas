!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Francesco Aquilante                                    *
!               2021, Roland Lindh                                     *
!***********************************************************************
      SUBROUTINE swap_tosqrt(irc,iLoc,nRS,JSYM,XLT,Xab)
      use ChoArr, only: iRS2F
      use Data_Structures, only: NDSBA_Type

      Implicit Real*8 (a-h,o-z)
      Integer  irc, iLoc, JSYM
      Type (NDSBA_Type) XLT
      Real*8 Xab(nRS)

      Integer, External:: cho_isao

#include "cholesky.fh"
#include "choorb.fh"

      Integer i, j, MulD2h
!                                                                      *
!***********************************************************************
!                                                                      *
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
!                                                                      *
!***********************************************************************
!                                                                      *
      If (JSYM.ne.1) then
!      ! NON TOTAL-SYMMETRIC

         Do jRab=1,nnBstR(jSym,iLoc)

            kRab = iiBstr(jSym,iLoc) + jRab ! already in 1st red set

            iag   = iRS2F(1,kRab)  !global address
            ibg   = iRS2F(2,kRab)

            iSyma = cho_isao(iag)  !symmetry block
            iSymb = MulD2h(jSym,iSyma) ! sym(a) .gt. sym(b)

            ias   = iag - ibas(iSyma)
            ibs   = ibg - ibas(iSymb)

            XLT%SB(iSyma,iSymb)%A2(ias,ibs) = sqrt(abs(Xab(kRab)))

         End Do  ! jRab loop


      ElseIf (JSYM.eq.1) then

         Do jRab=1,nnBstR(jSym,iLoc)

            kRab = iiBstr(jSym,iLoc) + jRab ! already in 1st red set

            iag   = iRS2F(1,kRab)  !global address
            ibg   = iRS2F(2,kRab)

            iSyma = cho_isao(iag)  ! sym(a)=sym(b)

            ias   = iag - ibas(iSyma)  !address within that symm block
            ibs   = ibg - ibas(iSyma)

            XLT%SB(iSyma,iSyma)%A2(ias,ibs) = sqrt(abs(Xab(kRab)))
            XLT%SB(iSyma,iSyma)%A2(ibs,ias) = sqrt(abs(Xab(kRab)))

         End Do  ! jRab loop

      EndIf

      irc = 0

      Return
      End

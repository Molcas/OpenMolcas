************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Francesco Aquilante                                    *
************************************************************************
      SUBROUTINE play_rassi_sto(irc,iLoc,JSYM,ISLT,ISSQ,
     &                        ipXLT,ipXab,mode)
      use ChoArr, only: iRS2F
      use ChoSwp, only: IndRed
      Implicit Real*8 (a-h,o-z)
      Integer ISLT(8),ISSQ(8,8),cho_isao
      External cho_isao
      Integer ipXLT,ipXab
      Character*6 mode

#include "cholesky.fh"
#include "choorb.fh"
#include "WrkSpc.fh"

************************************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
******
      iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
************************************************************************



      If (mode.eq.'toreds'.and.JSYM.eq.1) then ! TOTAL SYMMETRIC

         Do jRab=1,nnBstR(jSym,iLoc)

            kRab = iiBstr(jSym,iLoc) + jRab
            iRab = IndRed(kRab,iLoc)

            iag   = iRS2F(1,iRab)  !global address
            ibg   = iRS2F(2,iRab)

            iSyma = cho_isao(iag)  !symmetry block; Sym(b)=Sym(a)

            ias   = iag - ibas(iSyma)  !address within that symm block
            ibs   = ibg - ibas(iSyma)
            iab   = iTri(ias,ibs)

            kfrom = ipXLT + isLT(iSyma) + iab - 1

            Work(ipXab+jRab-1) = Work(kfrom)

         End Do  ! jRab loop


      ElseIf (mode.eq.'tofull'.and.JSYM.eq.1) then  ! TOTAL SYMMETRIC

         Do jRab=1,nnBstR(jSym,iLoc)

            kRab = iiBstr(jSym,iLoc) + jRab
            iRab = IndRed(kRab,iLoc)

            iag   = iRS2F(1,iRab)  !global address
            ibg   = iRS2F(2,iRab)

            iSyma = cho_isao(iag)  !symmetry block; Sym(b)=Sym(a)

            ias   = iag - ibas(iSyma)  !address within that symm block
            ibs   = ibg - ibas(iSyma)
            iab   = iTri(ias,ibs)

            kto = ipXLT + isLT(iSyma) + iab - 1

            Work(kto) = Work(kto)
     &                + Work(ipXab+jRab-1)

         End Do  ! jRab loop


      ElseIf (mode.eq.'tosqrt'.and.JSYM.ne.1) then
c      ! NON TOTAL-SYMMETRIC

         Do jRab=1,nnBstR(jSym,iLoc)

            kRab = iiBstr(jSym,iLoc) + jRab ! already in 1st red set

            iag   = iRS2F(1,kRab)  !global address
            ibg   = iRS2F(2,kRab)

            iSyma = cho_isao(iag)  !symmetry block
            iSymb = MulD2h(jSym,iSyma) ! sym(a) .gt. sym(b)

            ias   = iag - ibas(iSyma)
            ibs   = ibg - ibas(iSymb)

            iab   = nBas(iSyma)*(ibs-1) + ias

            kto = ipXLT - 1 + isSQ(iSyma,iSymb) + iab

            Work(kto) = sqrt(abs(Work(ipXab+kRab-1)))

         End Do  ! jRab loop


      ElseIf (mode.eq.'tosqrt'.and.JSYM.eq.1) then

         Do jRab=1,nnBstR(jSym,iLoc)

            kRab = iiBstr(jSym,iLoc) + jRab ! already in 1st red set

            iag   = iRS2F(1,kRab)  !global address
            ibg   = iRS2F(2,kRab)

            iSyma = cho_isao(iag)  ! sym(a)=sym(b)

            ias   = iag - ibas(iSyma)  !address within that symm block
            ibs   = ibg - ibas(iSyma)

            iab   = nBas(iSyma)*(ibs-1) + ias
            iba   = nBas(iSyma)*(ias-1) + ibs

            kto = ipXLT - 1 + isSQ(iSyma,iSyma)

            Work(kto+iab) = sqrt(abs(Work(ipXab+kRab-1)))

            Work(kto+iba) = sqrt(abs(Work(ipXab+kRab-1)))

         End Do  ! jRab loop


      Else

         write(6,*)'Wrong input parameters. JSYM,mode = ',JSYM,mode
         irc = 66
         Call abend()

      EndIf

      irc = 0

      Return
      End

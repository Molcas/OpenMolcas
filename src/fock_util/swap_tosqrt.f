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
*               2021, Roland Lindh                                     *
************************************************************************
      SUBROUTINE swap_tosqrt(irc,iLoc,nRS,JSYM,ipXLT,Xab)
      use ChoArr, only: iRS2F
      use ChoSwp, only: IndRed
      Implicit Real*8 (a-h,o-z)
      Integer  irc, iLoc, JSYM
      Integer ipXLT
      Real*8 Xab(nRS)

      Integer  ISSQ(8,8)
      Integer, External:: cho_isao

#include "cholesky.fh"
#include "choorb.fh"
#include "WrkSpc.fh"

      Integer i, j, MulD2h
*                                                                      *
************************************************************************
*                                                                      *
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
*                                                                      *
************************************************************************
*                                                                      *
      nnBSQ=0
      DO LSYM=1,NSYM
         DO KSYM=LSYM,NSYM
            ISSQ(KSYM,LSYM) = nnBSQ
            ISSQ(LSYM,KSYM) = nnBSQ ! symmetrization
            nnBSQ = nnBSQ + nBas(kSym)*nBas(lSym)
         END DO
      END DO

      If (JSYM.ne.1) then
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

            Work(kto) = sqrt(abs(Xab(kRab)))

         End Do  ! jRab loop


      ElseIf (JSYM.eq.1) then

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

            Work(kto+iab) = sqrt(abs(Xab(kRab)))

            Work(kto+iba) = sqrt(abs(Xab(kRab)))

         End Do  ! jRab loop

      EndIf

      irc = 0

      Return
      End

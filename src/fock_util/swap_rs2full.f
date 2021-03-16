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
      SUBROUTINE play_sto(irc,iLoc,nDen,JSYM,ipXLT,ipXab,mode,add)
      use ChoArr, only: iRS2F
      use ChoSwp, only: IndRed
      Implicit Real*8 (a-h,o-z)
      Integer  irc, iLoc, nDen, JSYM
      Integer ipXLT(nDen),ipXab(jDen)
      Logical add
      Character*6 mode

      Integer  ISLT(8),ISSQ(8,8)
      Integer, External:: cho_isao

#include "cholesky.fh"
#include "choorb.fh"
#include "WrkSpc.fh"

      Integer i, j, MulD2h, iTri
*                                                                      *
************************************************************************
*                                                                      *
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
      iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
*                                                                      *
************************************************************************
*                                                                      *

      ISLT(1)=0
      DO ISYM=2,NSYM
         ISLT(ISYM) = ISLT(ISYM-1)
     &              + NBAS(ISYM-1)*(NBAS(ISYM-1)+1)/2
      END DO

      nnBSQ=0
      DO LSYM=1,NSYM
         DO KSYM=LSYM,NSYM
            ISSQ(KSYM,LSYM) = nnBSQ
            ISSQ(LSYM,KSYM) = nnBSQ ! symmetrization
            nnBSQ = nnBSQ + nBas(kSym)*nBas(lSym)
         END DO
      END DO

      If (mode.eq.'toreds'.and.JSYM.eq.1) then ! TOTAL SYMMETRIC

         Do jRab=1,nnBstR(jSym,iLoc)

            kRab = iiBstr(jSym,iLoc) + jRab
            iRab = IndRed(kRab,iLoc)

            iag   = iRS2F(1,iRab)  !global address
            ibg   = iRS2F(2,iRab)

            iSyma = cho_isao(iag)  !symmetry block; Sym(b)=Sym(a)

            ias   = iag - ibas(iSyma)
c           !address within that symm block
            ibs   = ibg - ibas(iSyma)
            iab   = iTri(ias,ibs)

            Do jDen=1,nDen

               kfrom = ipXLT(jDen) + isLT(iSyma) + iab - 1

               Work(ipXab(jDen)+jRab-1) =  Work(kfrom)

            End Do

         End Do  ! jRab loop


      ElseIf (mode.eq.'tofull'.and.JSYM.eq.1) then
c      ! TOTAL SYMMETRIC

         If (.NOT.add) Then
            nTot = ISLT(NSYM) + NBAS(NSYM)*(NBAS(NSYM)+1)/2
            Do jDen = 1, nDen
               Call FZero(Work(ipXLT(jDen)),nTot)
            End Do
         End If

         Do jRab=1,nnBstR(jSym,iLoc)

            kRab = iiBstr(jSym,iLoc) + jRab
            iRab = IndRed(kRab,iLoc)

            iag   = iRS2F(1,iRab)  !global address
            ibg   = iRS2F(2,iRab)

            iSyma = cho_isao(iag)  !symmetry block; Sym(b)=Sym(a)

            ias   = iag - ibas(iSyma)  !address within that symm block
            ibs   = ibg - ibas(iSyma)
            iab   = iTri(ias,ibs)

            Do jDen=1,nDen

               kto = ipXLT(jDen) + isLT(iSyma) + iab - 1

               Work(kto) = Work(kto) + Work(ipXab(jDen)+jRab-1)

            End Do

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

            Do jDen=1,nDen

               kto = ipXLT(jDen) - 1 + isSQ(iSyma,iSymb) + iab

               Work(kto) = sqrt(abs(Work(ipXab(jDen)+kRab-1)))

            End Do

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

            Do jDen=1,nDen

               kto = ipXLT(jDen) - 1 + isSQ(iSyma,iSyma)

               Work(kto+iab) = sqrt(abs(Work(ipXab(jDen)+kRab-1)))

               Work(kto+iba) = sqrt(abs(Work(ipXab(jDen)+kRab-1)))

            End Do

         End Do  ! jRab loop


      Else

         write(6,*)'Wrong input parameters. JSYM,mode = ',JSYM,mode
         irc = 66
         Call abend()

      EndIf

      irc = 0

      Return
      End

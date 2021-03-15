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
      SUBROUTINE swap_rs2full(irc,iLoc,nDen,ipXLT,ipXab,mode,add)
      use ChoArr, only: iRS2F
      use ChoSwp, only: IndRed
      Implicit Real*8 (a-h,o-z)
      Integer  nDen
      Integer ipXLT(nDen),ipXab(nDen)
      Logical add
      Character*6 mode

      Integer, External :: cho_isao
      Integer  ISLT(8)
      Integer, Parameter:: JSYM=1

#include "cholesky.fh"
#include "choorb.fh"
#include "WrkSpc.fh"

************************************************************************
      iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
************************************************************************
      ISLT(1)=0
      DO ISYM=2,NSYM
         ISLT(ISYM) = ISLT(ISYM-1)
     &              + NBAS(ISYM-1)*(NBAS(ISYM-1)+1)/2
      END DO

      If (mode.eq.'toreds') then ! TOTAL SYMMETRIC

         If (.NOT.add) Then
            nTot = nnBstR(jSym,iLoc)
            Call FZero(Work(ipXab(1)),nDen*nTot)
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

               kfrom = ipXLT(jDen) + isLT(iSyma) + iab - 1

               Work(ipXab(jDen)+jRab-1) = Work(ipXab(jDen)+jRab-1)
     &                                  + Work(kfrom)

            End Do

         End Do  ! jRab loop

      ElseIf (mode.eq.'tofull') then  ! TOTAL SYMMETRIC

         If (.NOT.add) Then
            nTot = ISLT(NSYM) + NBAS(NSYM)*(NBAS(NSYM)+1)/2
            Call FZero(Work(ipXLT(1)),nDen*nTot)
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

               Work(kto) = Work(kto)
     &                   + Work(ipXab(jDen)+jRab-1)

            End Do

         End Do  ! jRab loop


      Else

         write(6,*)'Wrong input parameters. mode = ',mode
         irc = 66
         Call abend()

      EndIf

      irc = 0

      Return
      End

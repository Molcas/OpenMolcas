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
      SUBROUTINE swap_rs2full(irc,iLoc,nRS,nDen,JSYM,XLT,Xab,mode,add)
      use ChoArr, only: iRS2F
      use ChoSwp, only: IndRed
      use Data_Structures, only: DSBA_Type
      Implicit Real*8 (a-h,o-z)
      Integer  irc, iLoc, nDen, JSYM
      Type (DSBA_Type) XLT(nDen)
      Real*8 Xab(nRS,nDen)
      Logical add
      Character*6 mode

      Integer, External:: cho_isao

#include "real.fh"
#include "cholesky.fh"
#include "choorb.fh"

      Integer i, j, iTri
*                                                                      *
************************************************************************
*                                                                      *
      iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
*                                                                      *
************************************************************************
*                                                                      *
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

               Xab(jRab,jDen) =  XLT(jDen)%SB(iSyma)%A1(iab)

            End Do

         End Do  ! jRab loop


      ElseIf (mode.eq.'tofull'.and.JSYM.eq.1) then
c      ! TOTAL SYMMETRIC

         If (.NOT.add) Then
            Do jDen = 1, nDen
               XLT(jDen)%A0(:)=Zero
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

               XLT(jDen)%SB(iSyma)%A1(iab) = XLT(jDen)%SB(iSyma)%A1(iab)
     &                                     + Xab(jRab,jDen)

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

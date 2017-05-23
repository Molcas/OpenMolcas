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
* Copyright (C) 2005, Giovanni Ghigo                                   *
************************************************************************
      Subroutine MkExSB11(iAddSB,LenSB,
     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
************************************************************************
* Author :  Giovanni Ghigo                                             *
*           Lund University, Sweden & Torino University, Italy         *
*           February 2005                                              *
*----------------------------------------------------------------------*
* Purpuse:  Generation of the SubBlock(1,1) (p,q inactive) of the      *
*           two-electron integral matrix for each i,j occupied couple. *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
#include "rasdim.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "cho_tra.fh"
      Logical SameLx

*   - SubBlock 1 1
      LenSB = nIsh(iSymA) * nIsh(iSymB)
      Call GetMem('SB','Allo','Real',iAddSB,LenSB)

*     Build Lx
      Call GetMem('Lx','Allo','Real',iAddLx0,nIsh(iSymA)*numV)
      LxType=0
      iIx=0
      SameLx=.False.
      Call MkL1(iSymA,iSymI,iI,numV, LxType,iIx, iAddLx0,SameLx)

*     Build Ly
      Call GetMem('Ly','Allo','Real',iAddLy0,nIsh(iSymB)*numV)
      If(iSymA.EQ.iSymB) SameLx=.True.
      Call MkL1(iSymB,iSymJ,iJ,numV, LxType,iIx, iAddLy0,SameLx)

*     Generate the SubBlock (Ly*Lx)
      If (.NOT.SameLx) then
        Call DGEMM_('N','T',nIsh(iSymB),nIsh(iSymA),numV,1.0d0,
     &    Work(iAddLy0),nIsh(iSymB), Work(iAddLx0),nIsh(iSymA),
     &                     0.0d0,Work(iAddSB),nIsh(iSymB) )
      else
        Call DGEMM_('N','T',nIsh(iSymA),nIsh(iSymA),numV,1.0d0,
     &    Work(iAddLx0),nIsh(iSymA), Work(iAddLx0),nIsh(iSymA),
     &                     0.0d0,Work(iAddSB),nIsh(iSymA) )
      EndIf

      Call GetMem('Ly','Free','Real',iAddLy0,nIsh(iSymB)*numV)
      Call GetMem('Lx','Free','Real',iAddLx0,nIsh(iSymA)*numV)

      Return
      End

      Subroutine MkCouSB11(iAddSB,LenSB,
     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
************************************************************************
* Author :  Giovanni Ghigo                                             *
*           Lund University, Sweden & Torino University, Italy         *
*           July 2005                                                  *
*----------------------------------------------------------------------*
* Purpuse:  Generation of the SubBlock(1,1) (p,q inactive) of the      *
*           two-electron integral matrix for each i,j occupied couple. *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
#include "rasdim.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "cho_tra.fh"

*   - SubBlock 1 1
      LenSB = nIsh(iSymA) * nIsh(iSymB)
      Call GetMem('SB','Allo','Real',iAddSB,LenSB)

*     Define Lab
      iAddAB = iMemTCVX(1,iSymA,iSymB,1)
      LenAB  = LenSB
CGG   ------------------------------------------------------------------
c      If(IfTest) then
c      Write(6,*)'     MkCouSB11: TCVA(',iSymA,',',iSymB,')'
c      Write(6,'(8F10.6)')(Work(iAddAB+k),k=0,LenAB*numV-1)
c      Call XFlush(6)
c      EndIf
CGG   ------------------------------------------------------------------

*     Build Lij
      LenLij = numV
      Call GetMem('Lij','Allo','Real',iAddLij,LenLij)
      Call MkLij(iSymI,iSymJ,iI,iJ,numV, iAddLij)

*     Generate the SubBlock
      Call DGEMM_('N','N',LenAB,1,numV,1.0d0,
     &    Work(iAddAB),LenAB, Work(iAddLij),LenLij,
     &                0.0d0,Work(iAddSB),LenSB )

      Call GetMem('Lij','Free','Real',iAddLij,LenLij)

      Return
      End

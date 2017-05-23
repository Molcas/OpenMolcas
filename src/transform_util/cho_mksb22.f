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
      Subroutine MkExSB22(iAddSB,LenSB,
     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
************************************************************************
* Author :  Giovanni Ghigo                                             *
*           Lund University, Sweden & Torino University, Italy         *
*           February 2005                                              *
*----------------------------------------------------------------------*
* Purpuse:  Generation of the SubBlock(2,2) (p,q active) of the        *
*           two-electron integral matrix for each i,j occupied couple. *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
#include "rasdim.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "cho_tra.fh"
      Logical SameLx

*   - SubBlock 2 2
      LenSB = nAsh(iSymA) * nAsh(iSymB)
      Call GetMem('SB','Allo','Real',iAddSB,LenSB)

*     Build Lx
      Call GetMem('Lx','Allo','Real',iAddLx0,nAsh(iSymA)*numV)
      LxType=0
      iIx=0
      SameLx=.False.
      Call MkL2(iSymA,iSymI,iI,numV, LxType,iIx, iAddLx0,SameLx)

*     Build Ly
      Call GetMem('Ly','Allo','Real',iAddLy0,nAsh(iSymB)*numV)
      If(iSymA.EQ.iSymB) SameLx=.True.
      Call MkL2(iSymB,iSymJ,iJ,numV, LxType,iIx, iAddLy0,SameLx)

*     Generate the SubBlock
      If (.NOT.SameLx) then
        Call DGEMM_('N','T',nAsh(iSymB),nAsh(iSymA),numV,1.0d0,
     &    Work(iAddLy0),nAsh(iSymB), Work(iAddLx0),nAsh(iSymA),
     &                     0.0d0,Work(iAddSB),nAsh(iSymB) )
      else
        Call DGEMM_('N','T',nAsh(iSymA),nAsh(iSymA),numV,1.0d0,
     &    Work(iAddLx0),nAsh(iSymA), Work(iAddLx0),nAsh(iSymA),
     &                     0.0d0,Work(iAddSB),nAsh(iSymA) )
      EndIf

      Call GetMem('Ly','Free','Real',iAddLy0,nAsh(iSymB)*numV)
      Call GetMem('Lx','Free','Real',iAddLx0,nAsh(iSymA)*numV)

      Return
      End

      Subroutine MkCouSB22(iAddSB,LenSB,
     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
************************************************************************
* Author :  Giovanni Ghigo                                             *
*           Lund University, Sweden & Torino University, Italy         *
*           July 2005                                                  *
*----------------------------------------------------------------------*
* Purpuse:  Generation of the SubBlock(2,2) (p,q active) of the        *
*           two-electron integral matrix for each i,j occupied couple. *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
#include "rasdim.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "cho_tra.fh"

*   - SubBlock 2 2
      LenSB = nAsh(iSymA) * nAsh(iSymB)
      Call GetMem('SB','Allo','Real',iAddSB,LenSB)

*     Define Lab
      iAddAB = iMemTCVX(4,iSymA,iSymB,1)
      LenAB  = LenSB
CGG   ------------------------------------------------------------------
c      If(IfTest) then
c      Write(6,*)'     MkCouSB22: TCVD(',iSymA,',',iSymB,')'
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

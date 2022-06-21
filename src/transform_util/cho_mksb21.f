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
      Subroutine MkExSB21(AddSB,iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV,
     &                    AddSBt)
************************************************************************
* Author :  Giovanni Ghigo                                             *
*           Lund University, Sweden & Torino University, Italy         *
*           February 2005                                              *
*----------------------------------------------------------------------*
* Purpuse:  Generation of the SubBlock(2,1) (p,q active,inactive) of   *
*           two-electron integral matrix for each i,j occupied couple. *
************************************************************************
      use Cho_Tra
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
      Real*8, Allocatable:: AddSB(:)
      Integer iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV
      Real*8 AddSBt(*)
#include "rasdim.fh"
#include "stdalloc.fh"
#include "SysDef.fh"
      Logical SameLx

      Real*8, Allocatable:: Lx0(:), Ly0(:)

*   - SubBlock 2 1
      LenSB = nAsh(iSymA) * nIsh(iSymB)
      Call mma_allocate(AddSB,LenSB,Label='AddSB')
      If (iSymA.EQ.iSymB .and. iSymI.EQ.iSymJ .and. iI.EQ.iJ) then
*       SB 2,1 = (SB 1,2)+
        Call Trnsps(nAsh(iSymB),nIsh(iSymA),AddSBt,AddSB)
        Return
      EndIf

*     Build Lx
      Call mma_allocate(Lx0,nAsh(iSymA)*numV,Label='Lx0')
      LxType=0
      iIx=0
      SameLx=.False.
      Call MkL2(iSymA,iSymI,iI,numV, LxType,iIx, Lx0,SameLx)

*     Build Ly
      Call mma_allocate(Ly0,nIsh(iSymB)*numV,Label='Ly0')
      Call MkL1(iSymB,iSymJ,iJ,numV, LxType,iIx, Ly0,SameLx)

*     Generate the SubBlock
      Call DGEMM_('N','T',nIsh(iSymB),nAsh(iSymA),numV,
     &            1.0d0,Ly0,nIsh(iSymB),
     &                  Lx0,nAsh(iSymA),
     &            0.0d0,AddSB,nIsh(iSymB) )

      Call mma_deallocate(Ly0)
      Call mma_deallocate(Lx0)

      Return
      End


      Subroutine MkCouSB21(AddSB,iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
************************************************************************
* Author :  Giovanni Ghigo                                             *
*           Lund University, Sweden & Torino University, Italy         *
*           July 2005                                                  *
*----------------------------------------------------------------------*
* Purpuse:  Generation of the SubBlock(2,1) (p active, q inactive) of  *
*           two-electron integral matrix for each i,j occupied couple. *
************************************************************************
      use Cho_Tra
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
      Real*8, Allocatable:: AddSB(:)
#include "rasdim.fh"
#include "stdalloc.fh"
#include "SysDef.fh"

      Real*8, Allocatable:: Lij(:)
      Real*8, Allocatable:: AddSBt(:)

*   - SubBlock 2 1
      LenSB = nAsh(iSymA) * nIsh(iSymB)
      Call mma_allocate(AddSB,LenSB,Label='AddSB')

      Call mma_allocate(AddSBt,LenSB,Label='AddSBt')

*     Define Lab
      LenAB  = LenSB
CGG   ------------------------------------------------------------------
c      If(IfTest) then
c      Write(6,*)'     MkCouSB21: TCVB(',iSymA,',',iSymB,')'
c      Write(6,'(8F10.6)') TCVX(2,iSymA,iSymB)%A(:,:)
c      Call XFlush(6)
c      EndIf
CGG   ------------------------------------------------------------------

*     Build Lij
      Call mma_allocate(Lij,NumV,Label='Lij')
      Call MkLij(iSymI,iSymJ,iI,iJ,numV,Lij)

*     Generate the SubBlock
      Call DGEMM_('N','N',LenAB,1,numV,
     &            1.0d0,TCVX(2,iSymA,iSymB)%A,LenAB,
     &                  Lij,NumV,
     &            0.0d0,AddSBt,LenSB )
      Call Trnsps(nAsh(iSymA),nIsh(iSymB),AddSBt,AddSB)

      Call mma_deallocate(Lij)
      Call mma_deallocate(AddSBt)

      Return
      End

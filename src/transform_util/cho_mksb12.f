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
      Subroutine MkExSB12(AddSB,iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
************************************************************************
* Author :  Giovanni Ghigo                                             *
*           Lund University, Sweden & Torino University, Italy         *
*           February 2005                                              *
*----------------------------------------------------------------------*
* Purpuse:  Generation of the SubBlock(1,2) (p,q inactive,active) of   *
*           two-electron integral matrix for each i,j occupied couple. *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
      Real*8, Allocatable:: AddSB(:)
      Integer iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV
#include "rasdim.fh"
#include "stdalloc.fh"
#include "SysDef.fh"
#include "cho_tra.fh"
      Logical SameLx

      Real*8, Allocatable:: Lx0(:), Ly0(:)

*   - SubBlock 1 2
      LenSB = nIsh(iSymA) * nAsh(iSymB)
      Call mma_allocate(AddSB,LenSB,Label='AddSB')

*     Build Lx
      Call mma_allocate(Lx0,nIsh(iSymA)*numV,Label='Lx0')
      LxType=0
      iIx=0
      SameLx=.False.
      Call MkL1(iSymA,iSymI,iI,numV, LxType,iIx, Lx0,SameLx)

*     Build Ly
      Call mma_allocate(Ly0,nAsh(iSymB)*numV,Label='Ly0')
      Call MkL2(iSymB,iSymJ,iJ,numV, LxType,iIx, Ly0,SameLx)

*     Generate the SubBlock
      Call DGEMM_('N','T',nAsh(iSymB),nIsh(iSymA),numV,
     &            1.0d0,Ly0,nAsh(iSymB),
     &                  Lx0,nIsh(iSymA),
     &            0.0d0,AddSB,nAsh(iSymB) )

      Call mma_deallocate(Ly0)
      Call mma_deallocate(Lx0)

      Return
      End

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
      Subroutine MkExSB23(AddSB,iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
************************************************************************
* Author :  Giovanni Ghigo                                             *
*           Lund University, Sweden & Torino University, Italy         *
*           February 2005                                              *
*----------------------------------------------------------------------*
* Purpuse:  Generation of the SubBlock(2,3) (p,q active,secondary) of  *
*           two-electron integral matrix for each i,j occupied couple. *
************************************************************************
      use Cho_Tra
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
      Real*8, Allocatable:: AddSB(:)
      Integer iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV
#include "rasdim.fh"
#include "stdalloc.fh"
#include "SysDef.fh"
      Logical SameLx

      Real*8, Allocatable:: Lx0(:), Ly0(:)

*   - SubBlock 2 3
      LenSB = nAsh(iSymA) * nSsh(iSymB)
      Call mma_allocate(AddSB,LenSB,Label='AddSB')

*     Build Lx
      Call mma_allocate(Lx0,nAsh(iSymA)*numV,Label='Lx0')
      LxType=0
      iIx=0
      SameLx=.False.
      Call MkL2(iSymA,iSymI,iI,numV, LxType,iIx, Lx0,SameLx)

*     Build Ly
      Call mma_allocate(Ly0,nSsh(iSymB)*numV,Label='Ly0')
      Call MkL3(iSymB,iSymJ,iJ,numV, LxType,iIx, Ly0,SameLx)

*     Generate the SubBlock
      Call DGEMM_('N','T',nSsh(iSymB),nAsh(iSymA),numV,
     &            1.0d0,Ly0,nSsh(iSymB),
     &                  Lx0,nAsh(iSymA),
     &            0.0d0,AddSB,nSsh(iSymB) )

      Call mma_deallocate(Ly0)
      Call mma_deallocate(Lx0)

      Return
      End

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
      Subroutine ChoMP2_GenE(iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV,
     &                       AddEx,LenE)
************************************************************************
* Author :  Giovanni Ghigo                                             *
*           Lund University, Sweden & Torino University, Italy         *
*           Jannuary-February 2005                                     *
* Modified for Cholesky-MP2 May 2005                                   *
************************************************************************
      use Cho_Tra
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
      Integer iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV,LenE
      Real*8 AddEx(LenE)
#include "rasdim.fh"
#include "stdalloc.fh"
#include "SysDef.fh"
      Logical SameLx

      Real*8, Allocatable:: Lx0(:), Ly0(:)

*   - SubBlock 3 3

*     Build Lx
      Call mma_allocate(Lx0,nSsh(iSymA)*numV,Label='Lx0')
      LxType=0
      iIx=0
      SameLx=.False.
      Call MkL3(iSymA,iSymI,iI,numV, LxType,iIx, Lx0,SameLx)

*     Build Ly
      Call mma_allocate(Ly0,nSsh(iSymB)*numV,Label='Ly0')
      If(iSymA.EQ.iSymB) SameLx=.True.
      Call MkL3(iSymB,iSymJ,iJ,numV, LxType,iIx, Ly0,SameLx)

*     Generate the SubBlock
      If (.NOT.SameLx) then
        Call DGEMM_('N','T',nSsh(iSymB),nSsh(iSymA),numV,
     &               1.0d0,Ly0,nSsh(iSymB),
     &                     Lx0,nSsh(iSymA),
     &               1.0d0,AddEx,nSsh(iSymB) )
      else
        Call DGEMM_('N','T',nSsh(iSymA),nSsh(iSymA),numV,
     &              1.0d0,Lx0,nSsh(iSymA),
     &                    Lx0,nSsh(iSymA),
     &              1.0d0,AddEx,nSsh(iSymA) )
      EndIf

      Call mma_deallocate(Ly0)
      Call mma_deallocate(Lx0)

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(LenE)
      End

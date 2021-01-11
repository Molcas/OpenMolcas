************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine FIXIC(nFix,SS,mInt,B,NDIM,F,Label,u)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
      Real*8 SS(mInt), B(nDim*mInt), F(nDim), u(nDim)
      Character(LEN=8) Label(mInt)
      Real*8, Allocatable:: uInv(:,:), uB(:,:)
*
*
*     write out the internal coordinates which will be fixed
*
      WRITE (6,*)
      WRITE (6,*)
     &      ' Following internal coordinates are fixed'
      WRITE (6,*)
*
*     loop over all internal coordinates to be fixed
*
      Do I = mInt-nFix+1,mInt
         WRITE (6,'(A,A,E10.3,A)')
     &             Label(i),' with a gradient of ',SS(I),
     &             ' is frozen and the gradient is annihilated'
         SS(i) = Zero
      End Do
*
*     now transform remaining internal coordinates back to cartesian ba
*                          -1 +
*                    fx = u  B  fq
*
      Call mma_allocate(uInv,nDim,nDim,Label='uInv')
      uInv(:,:)=Zero
      Do i = 1, nDim
         uInv(i,i)=One/u(i)
      End Do
*     Call RecPrt('uInv',' ',uInv,nDim,nDim)

      Call mma_allocate(uB,mInt,nDim,Label='uB')
      uB(:,:)=Zero

      Call DGEMM_('N','N',
     &            nDim,mInt,nDim,
     &            One,uInv,nDim,
     &            B,nDim,
     &            Zero,uB,nDim)
*     Call RecPrt('uInvB',' ',uB,nDim,mInt)
      Call DGEMM_('N','N',
     &            nDim,1,mInt,
     &            One,uB,nDim,
     &            SS,mInt,
     &            Zero,F,nDim)
*     Call RecPrt('F',' ',F,mInt,1)

      Call mma_deallocate(uB)
      Call mma_deallocate(uInv)
*
      Return
      End

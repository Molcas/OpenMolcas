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
      Subroutine TRPGen(nDim,nAtom,Coor,mTR,CofM,TRVec)
      use Slapaf_Info, only: Degen, Smmtrc
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
#include "print.fh"
      Real*8 Coor(3,nAtom), TRVec(3*nAtom*6)
      Logical CofM
      Logical, Save:: g12K=.True.
      Real*8, Allocatable:: TR(:), Scrt(:), G(:), EVal(:), EVec(:),
     &                      U(:)
*
      iRout=135
      iPrint = nPrint(iRout)
*
      Call mma_allocate(TR,18*nAtom,Label='TR')
*
*-----Compute the symmetric translations and rotations
*
*     B    (nTR x nDim)
*      tr
*
      Call TRMake(TR,Coor,nAtom,nTR,Degen,nDim,CofM)
*
      TRVec(1:nTR*nDim) = TR(1:nTR*nDim)
*
      Call mma_allocate(Scrt,(3*nAtom)*nTR,Label='Scrt')
      Call mma_allocate(G,nTR**2,Label='G')
      Call mma_allocate(EVal,nTR*(nTR+1)/2,Label='EVal')
      Call mma_allocate(EVec,nTR**2,Label='EVec')
*
*-----Eliminate redundancy and produce an orthogonal representation.
*
*        -1/2
*     K g        (nTR x mTR)
*
      Call mma_allocate(U,nDim,Label='U')
      U(:) = One
      i=0
      Do iAtom = 1, nAtom
      Do ixyz = 1, 3
         If (Smmtrc(ixyz,iAtom)) Then
            i = i + 1
            Call DScal_(nTR,Sqrt(Degen(ixyz,iAtom)),
     &                  TRVec((i-1)*nTR+1),1)
         End If
      End Do
      End Do
*
      Thr_ElRed=1.0D-12
      Call ElRed(TRVec,nTR,nDim,G,EVal,EVec,mTR,U,Scrt,g12K,Thr_ElRed)
*
      If (mTR.gt.0) Then
         TRVec(1:3*nAtom*nTR) = Zero
         Call DGEMM_('T','N',
     &               nDim,mTR,nTR,
     &               1.0d0,TR,nTR,
     &                     EVec,nTR,
     &               0.0d0,TRVec,nDim)
      End If
*
      Call mma_deallocate(U)
      Call mma_deallocate(EVec)
      Call mma_deallocate(EVal)
      Call mma_deallocate(G)
      Call mma_deallocate(Scrt)
      Call mma_deallocate(TR)
*
      Return
      End

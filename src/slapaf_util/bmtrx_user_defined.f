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
      Subroutine BMtrx_User_Defined(nsAtom,Coor,nDim,nIter,mTR,nQQ)
      use Slapaf_Info, only: Gx, qInt, dqInt, KtB, BMx, Degen, Smmtrc,
     &                       Lbl
      use Slapaf_Parameters, only: iInt, nFix, nBVec, Analytic_Hessian,
     &                             MaxItr, iOptC, BSet, HSet, lOld,
     &                             Numerical
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
#include "stdalloc.fh"
      Real*8 Coor(3,nsAtom)
      Logical Proc_dB
      Real*8, Allocatable:: Degen2(:)
*                                                                      *
************************************************************************
*                                                                      *
*#define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
*.... Section for user defined internal coordinates
*
      Call Rd_UDIC(iInt,nFix,nRowH) ! nRowH is not used!
      nQQ=iInt+nFix
*
      If (Allocated(qInt).and.SIZE(qInt,1)/=nQQ) Then
         Call mma_deallocate(qInt)
         Call mma_deallocate(dqInt)
      End If
      If (.NOT.Allocated(qInt)) Then
         Call mma_allocate(qInt,nQQ,MaxItr,Label='qInt')
         Call mma_allocate(dqInt,nQQ,MaxItr,Label='dqInt')
         qInt(:,:) = Zero
         dqInt(:,:) = Zero
      End If
      Call mma_allocate(BMx,3*nsAtom,nQQ,Label='BMx')
      BMx(:,:)=Zero

*
*-----Compute the B matrix in symmetry distinct basis and the
*     internal coordinates.
*
*     iOptC(256) = constrained optimization
      Proc_dB=HSet.and..Not.lOld.and.
     &           (Analytic_Hessian.or.Numerical.or.
     &            iAnd(iOptC,256).eq.256)
*     Compute and store dBQQ in the reference structure
      If (Proc_dB) Then
*        Not implimented, sorry
      End If
*
      Call DefInt(nBVec,BMx,nQQ,nsAtom,qInt(:,nIter),Lbl,Coor,nDim-mTR)
*                                                                      *
************************************************************************
*                                                                      *
*     Compute the gradient
*
      If (BSet) Call Force(nFix,Gx(:,:,nIter),nsAtom,nQQ,BMx,
     &                     nIter,dqInt,Lbl,Degen)
*                                                                      *
************************************************************************
*                                                                      *
      If (HSet.and..NOT.lOld.and.BSet) Then
         Call mma_allocate(KtB,nDim,nQQ,Label='KtB')
*
         Call mma_allocate(Degen2,nDim,Label='Degen2')
         i=0
         Do ix = 1, 3*nsAtom
            iAtom = (ix+2)/3
            ixyz = ix - (iAtom-1)*3
            If (Smmtrc(ixyz,iAtom)) Then
               i = i + 1
               Degen2(i) = Degen(ixyz,iAtom)
            End If
         End Do
*
         Do j = 1, nQQ
            i = 0
            Do ix = 1, 3*nsAtom
               iAtom = (ix+2)/3
               ixyz = ix - (iAtom-1)*3
               If (Smmtrc(ixyz,iAtom)) Then
                  i = i + 1
                  KtB(i,j) = BMx(ix,j)
               End If
            End Do
         End Do
*
         Do iInter = 1, nQQ
            Do iDim = 1, nDim
*              KtB(iDim,iInter) = KtB(iDim,iInter) / Sqrt(Degen2(iDim))
               KtB(iDim,iInter) = KtB(iDim,iInter) / Degen2(iDim)
            End Do
         End Do
         Call mma_deallocate(Degen2)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End

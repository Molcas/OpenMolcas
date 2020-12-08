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
      Subroutine BMtrx_User_Defined(
     &                 nLines,nBVec,ipBMx,nAtom,nInter,
     &                 Lbl,Coor,nDim,
     &                 Name,Smmtrc,
     &                 Degen,BSet,HSet,nIter,
     &                 nStab,jStab,Numerical,
     &                 Analytic_Hessian,
     &                 iOptC,mxdc,lOld,
     &                 nFix,mTR,nQQ,Redundant,MaxItr)
      use Slapaf_Info, only: Gx, dMass, qInt, dqInt, KtB
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
      Real*8 Coor(3,nAtom), Degen(3*nAtom)
      Character Lbl(nInter)*8, Name(nAtom)*(LENIN)
      Integer   nStab(nAtom), jStab(0:7,nAtom)
      Logical Smmtrc(3*nAtom), BSet, HSet, Redundant,
     &        Numerical, Analytic_Hessian, lOld, Proc_dB
      Real*8, Allocatable:: Degen2(:)
      Character(LEN=8), Allocatable:: Lab(:)
      Real*8, Allocatable:: Mult(:), BVec(:), Val(:)
*                                                                      *
************************************************************************
*                                                                      *
*.... Section for user defined internal coordinates
*
      Call Rd_UDIC(nLines,iInt,nFix,nRowH)
      nQQ=iInt+nFix
*
      If (.NOT.Allocated(qInt)) Then
         Call mma_allocate(qInt,nQQ,MaxItr,Label='qInt')
         Call mma_allocate(dqInt,nQQ,MaxItr,Label='dqInt')
         qInt(:,:) = Zero
         dqInt(:,:) = Zero
      End If
      Call Allocate_Work(ipBmx,3*nAtom*nQQ)
      Call FZero(Work(ipBMx),3*nAtom*nQQ)

      Call mma_allocate(Mult,nBVec,Label='Mult')
      Call mma_allocate(BVec,nBVec*3*nAtom,Label='BVec')
      BVec(:)=Zero
      Call mma_allocate(Val,nBVec,Label='Val')
      Val(:)=Zero
      Call mma_allocate(Lab,nBVec,Label='Lab')
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
      Call DefInt(BVec,nBVec,Lab,Work(ipBMx),nQQ,
     &            nAtom,nLines,Val,qInt(:,nIter),Lbl,Name,
     &            Coor,dMass,jStab,nStab,mxdc,Mult,
     &            nDim-mTR,Redundant)
*
      Call mma_deallocate(Lab)
      Call mma_deallocate(Val)
      Call mma_deallocate(BVec)
      Call mma_deallocate(Mult)
*                                                                      *
************************************************************************
*                                                                      *
*     Compute the gradient
*
      If (BSet) Call Force(nFix,Gx(:,:,nIter),nAtom,nQQ,Work(ipBMx),
     &                     Name,nIter,dqInt,Lbl,Degen)
*                                                                      *
************************************************************************
*                                                                      *
      If (HSet.and..NOT.lOld.and.BSet) Then
         Call mma_allocate(KtB,nDim,nQQ,Label='KtB')
*
         Call mma_allocate(Degen2,nDim,Label='Degen2')
         i=0
         Do ix = 1, 3*nAtom
            If (Smmtrc(ix)) Then
               i = i + 1
               Degen2(i) = Degen(ix)
            End If
         End Do
*
         Do j = 1, nQQ
            i = 0
            Do ix = 1, 3*nAtom
               If (Smmtrc(ix)) Then
                  i = i + 1
                  ixj= (j-1)*3*nAtom + ix - 1 + ipBmx
                  KtB(i,j) = Work(ixj)
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

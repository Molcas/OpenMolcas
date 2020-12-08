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
      Subroutine BMtrx_Cartesian(
     &                 ipBMx,nAtom,nInter,nDim,
     &                 Name,Smmtrc,Degen,BSet,HSet,
     &                 nIter,mTtAtm,
     &                 PrQ,lOld,mTR,TRVec,EVal,Hss_x,
     &                 nQQ,Redundant,MaxItr,nWndw)
      use Slapaf_Info, only: Cx, Gx, qInt, dqInt, KtB
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "print.fh"
      Real*8 Degen(3*nAtom), TRVec(nDim,mTR)
      Character Name(nAtom)*(LENIN)
      Logical Smmtrc(3*nAtom), BSet, HSet, Redundant, PrQ, lOld
      Real*8 Eval(3*mTtAtm*(3*mTtAtm+1)/2)
      Real*8 Hss_x((3*mTtAtm)**2)
      Real*8, Allocatable:: EVec(:), Hi(:,:), iHi(:), Degen2(:)
      Integer, Allocatable:: Ind(:)

*                                                                      *
************************************************************************
*                                                                      *
      iRout=133
      iPrint=nPrint(iRout)
*
*-----Recompute the B matrix once each macroIteration, this is
*     not done if a numerical Hessian is computed.
*                                                                      *
************************************************************************
*                                                                      *
*     R E D U N D A N T  C A R T E S I A N  C O O R D S
*
      If (Redundant) Then
         nQQ = nDim
         If (.NOT.Allocated(qInt)) Then
            Call mma_allocate(qInt,nQQ,MaxItr,Label='qInt')
            Call mma_allocate(dqInt,nQQ,MaxItr,Label='dqInt')
            qInt(:,:) = Zero
            dqInt(:,:) = Zero
         End If
         Call mma_allocate(EVec,nDim**2,Label='EVec')
         EVec(:)=Zero
         call dcopy_(nDim,[One],0,EVec,nDim+1)
*                                                                      *
************************************************************************
*                                                                      *
*------- Move over the eigenvectors putting to BMx
*
         Call Allocate_Work(ipBMx,(3*nAtom)*nQQ)
         Call FZero(Work(ipBMx), (3*nAtom)*nQQ)
         ipFrom = 1
         Call BPut(EVec(ipFrom),nDim,Work(ipBMx),3*nAtom,Smmtrc,
     &             nQQ,Degen)
         If (iPrint.ge.19) Call RecPrt('In Bmtrx: B',' ',Work(ipBMx),
     &                                 3*nAtom,nQQ)
*                                                                      *
************************************************************************
*                                                                      *
         If (PrQ.and.nAtom.le.5)
     &      Call List2('Cartesian Redundant',
     &                 Name,Work(ipBMx),nAtom,nQQ,Smmtrc)
*                                                                      *
************************************************************************
*                                                                      *
*------- Project the model Hessian with respect to rotations and
*        translations. The eigenvalues are shifted to large positive
*        eigenvalues to effectively remove any displacements in the
*        rotational and translations directions and to make sure that
*        the matrix is not singular.
*
         Call mma_allocate(Ind,nDim,Label='Ind')
         iInd=0
         Do i = 1, 3*nAtom
            If (Smmtrc(i)) Then
               iInd=iInd+1
               Ind(iInd)=i
            End If
         End Do
*
*        Compute H|i>
*
         Call mma_allocate(Hi,nDim,mTR,Label='Hi')
         Call mma_allocate(iHi,mTR,Label='iHi')
         Hi(:,:)=Zero
*
         Do j = 1, mTR
            Do i = 1, nDim
               Temp = 0.0D0
               Do k = 1, nDim
                  kx = Ind(k)
                  ik=(i-1)*nDim+k
                  Temp = Temp
     &                 + Hss_X(ik) * Sqrt(Degen(kx))
     &                 * TRVec(k,j)
               End Do
               Hi(i,j) = Temp
            End Do
         End Do
*        Call RecPrt('Hi',' ',Hi,nDim,mTR)
         Do iTR = 1, mTR
            iHi(iTR) = DDot_(nDim,TRVec(1,iTR),1,Hi(:,iTR),1)
         End Do
*        Call RecPrt('iHi',' ',iHi,mTR,1)
*
         Do i = 1, nDim
            ix = Ind(i)
            Do j = 1, i
               jx = Ind(j)
               ij = (j-1)*nDim + i
               ji = (i-1)*nDim + j
               Temp = Half*(Hss_x(ij)+Hss_x(ji))
*define UNIT_MM
#ifndef UNIT_MM
*
*              Here we shift the eigenvectors corresponding to
*              translations and rotations to a large positive values.
*
               Do iTR = 1, mTR
                  Omega = 1.0D+5
                  Hii   = iHi(iTR)
                  Temp = Temp
     &                 + Sqrt(Degen(ix)) * (
     &                 - TRVec(i,iTR) * Hi(j,iTR)
     &                 - Hi(i,iTR) * TRVec(j,iTR)
     &                 + TRVec(i,iTR) * (Omega+Hii) * TRVec(j,iTR)
     &                                     )* Sqrt(Degen(jx))
               End Do
#endif
*
               Hss_X(ij)=Temp
               Hss_X(ji)=Temp
            End Do
         End Do
         Call mma_deallocate(iHi)
         Call mma_deallocate(Hi)
         Call mma_deallocate(Ind)
*
*        Clean up the gradient wrt translational and rotational
*        component.
*
*        |g> = |g> - Sum(TR) |i><i|g>
*
         If (BSet) Then
*
*           Call RecPrt('Gx',' ',Gx(:,:,nIter),1,3*nAtom)
*           Call RecPrt('TRVec',' ',TRVec,nDim,mTR)
*
            Do iTR = 1, mTR
*
*              <i|g>
*
               Temp=0.0D0
               iInd=0
               i=0
               Do iAtom = 1, nAtom
               Do j = 1, 3
                  i = i + 1
                  If (Smmtrc(i)) Then
                     iInd=iInd+1
                     Temp = Temp + Degen(i)*Gx(j,iAtom,nIter)
     &                                     *TRVec(iInd,iTR)
                  End If
               End Do
               End Do
*
               iInd=0
               i=0
               Do iAtom = 1, nAtom
               Do j = 1, 3
                  i = i + 1
                  If (Smmtrc(i)) Then
                     iInd=iInd+1
                     Gx(j,iAtom,nIter) = Gx(j,iAtom,nIter)
     &                                 - TRVec(iInd,iTR)*Temp
                  End If
               End Do
               End Do
*
            End Do
*
*           Call RecPrt('Gx',' ',Gx(:,:,nIter),1,3*nAtom)
*
         End If
*                                                                      *
************************************************************************
*                                                                      *
      Else
*                                                                      *
************************************************************************
*                                                                      *
*     N O N - R E D U N D A N T  C A R T E S I A N  C O O R D S
*
         nQQ=nInter
         If (.NOT.Allocated(qInt)) Then
            Call mma_allocate(qInt,nQQ,MaxItr,Label='qInt')
            Call mma_allocate(dqInt,nQQ,MaxItr,Label='dqInt')
            qInt(:,:) = Zero
            dqInt(:,:) = Zero
         End If
*
*------- Project the model Hessian with respect to rotations and
*        translations. The eigenvalues are shifted to negative
*        eigenvalues.
*
         Call mma_allocate(Ind,nDim,Label='Ind')
         iInd=0
         Do i = 1, 3*nAtom
            If (Smmtrc(i)) Then
               iInd=iInd+1
               Ind(iInd)=i
            End If
         End Do
*
*        Compute H|i>
*
         Call mma_allocate(Hi,nDim,mTR,Label='Hi')
         Call mma_allocate(iHi,mTR,Label='iHi')
         Hi(:,:)=Zero
*
         Do j = 1, mTR
            Do i = 1, nDim
               Temp = 0.0D0
               Do k = 1, nDim
                  kx = Ind(k)
                  ik=(i-1)*nDim+k
                  Temp = Temp
     &                 + Hss_X(ik) * Sqrt(Degen(kx))
     &                 * TRVec(k,j)
               End Do
               Hi(i,j) = Temp
            End Do
         End Do
*        Call RecPrt('Hi',' ',Hi,nDim,mTR)
         Do iTR = 1, mTR
            iHi(iTR) = DDot_(nDim,TRVec(1,iTR),1,Hi(:,iTR),1)
         End Do
*        Call RecPrt('iHi',' ',iHi,mTR,1)
*
         Do i = 1, nDim
            ix = Ind(i)
            Do j = 1, i
               jx = Ind(j)
               ijTri=i*(i-1)/2 + j
               ij = (j-1)*nDim + i
               ji = (i-1)*nDim + j
               EVal(ijTri) = Half*(Hss_X(ij)+Hss_X(ji))
*
*              Here we shift the eigenvectors corresponding to tran-
*              lations and rotations down to negative faked eigen-
*              values.
*
               Do iTR = 1, mTR
                  Omega = -DBLE(iTR)
                  Hii   = iHi(iTR)
                  Eval(ijTri) = Eval(ijTri)
     &                 + Sqrt(Degen(ix)) * (
     &                 - TRVec(i,iTR) * Hi(j,iTR)
     &                 - Hi(i,iTR) * TRVec(j,iTR)
     &                 + TRVec(i,iTR) * (Omega+Hii) * TRVec(j,iTR)
     &                                     )* Sqrt(Degen(jx))
               End Do
*
               Hss_X(ij)=EVal(ijTri)
               Hss_X(ji)=EVal(ijTri)
            End Do
         End Do
         Call mma_deallocate(iHi)
         Call mma_deallocate(Hi)
         Call mma_deallocate(Ind)
#ifdef _DEBUGPRINT_
         Call TriPrt(' The Projected Model Hessian','(5G20.10)',
     &               EVal,nDim)
         Call RecPrt(' The Projected Model Hessian','(5G20.10)',
     &               Hss_x,nDim,nDim)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*        Compute the eigen vectors for the Cartesian Hessian
*
         Call mma_allocate(EVec,(3*mTtAtm)**2,Label='EVec')
         Call Hess_Vec(mTtAtm,EVal,EVec,nAtom,nDim)
*                                                                      *
************************************************************************
*                                                                      *
*------- Move over the eigenvectors putting to BMx
*
         Call Allocate_Work(ipBMx,(3*nAtom)**2)
         Call FZero(Work(ipBMx), (3*nAtom)**2)
         ipFrom = 1 + mTR*nDim
         Call BPut(EVec(ipFrom),nDim,Work(ipBMx),3*nAtom,Smmtrc,
     &             nQQ,Degen)
         If (iPrint.ge.19) Call RecPrt('In Bmtrx: B',' ',Work(ipBMx),
     &                                 3*nAtom,nQQ)
*                                                                      *
************************************************************************
*                                                                      *
         If (PrQ.and.nAtom.le.5)
     &      Call List2('Cartesian Approximate Normal Modes',
     &                  Name,Work(ipBMx),nAtom,nQQ,Smmtrc)
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      If (HSet.and..NOT.lOld) Then
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
         call dcopy_(nDim*nQQ,EVec(ipFrom),1,KtB,1)
         Do iInter = 1, nQQ
            Do iDim = 1, nDim
               KtB(iDim,iInter) = KtB(iDim,iInter) / Sqrt(Degen2(iDim))
            End Do
         End Do
         Call mma_deallocate(Degen2)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deallocate(EVec)
*                                                                      *
************************************************************************
*                                                                      *
*
*---- Compute the value and gradient vectors in the new basis.
*
      Call ValANM(nAtom,nQQ,nIter,Work(ipBmx),Degen,qInt,Cx,'Values',
     &            nWndw)
      If (BSet) Call ValANM(nAtom,nQQ,nIter,Work(ipBMx),Degen,
     &                      dqInt,Gx,'Gradients',nWndw)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End

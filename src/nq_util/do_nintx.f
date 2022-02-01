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
************************************************************************
*                                                                      *
      Subroutine Do_NIntX(AOInt,mGrid,TabAO1,TabAO2,nBfn,nD,mAO,nFn)
*                                                                      *
************************************************************************
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "stdalloc.fh"
#include "nq_info.fh"
      Real*8 AOInt(nBfn,nBfn,nD), TabAO1(nFn,mGrid,nBfn,nD),
     &       TabAO2(mAO,mGrid,nBfn)
      Real*8, Allocatable:: A1(:,:), A2(:,:), A_tri(:)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*#define _ANALYSIS_
#ifdef _ANALYSIS_
      Write (6,*)
      Write (6,*)  ' Analysing TabAO1'
      Thr=1.0D-14
      Do iD = 1, nD
      Do iFn = 1, nFn
         lBfn = 0
         Total=0.0D0
         Do iBfn = 1, nBfn
            lGrid = 0
            Do iGrid = 1, mGrid
               If (Abs(TabAO1(iFn,iGrid,iBfn,iD))<Thr) lGrid = lGrid + 1
            End Do
            If (lGrid==mGrid) lBfn=lBfn+1
            Total = Total + DBLE(lGrid)/DBLE(mGrid)
         End Do
         Total = Total / DBLE(nBfn)
         Write (6,*) 'Sparsity analysis, iD, iFn', iD, iFn
         Write (6,*) ' Total sparcity in %:', 1.0D2 * Total
         Write (6,*) ' Complete Bfn sparcity in %:',
     &                1.0D2*DBLE(lBfn)/DBLE(nBfn)
         Write (6,*)
      End Do
      End Do
      Write (6,*)
      Write (6,*)  ' Analysing TabAO2'
      Do iAO = 1, mAO
         lBfn = 0
         Total=0.0D0
         Do iBfn = 1, nBfn
            lGrid = 0
            Do iGrid = 1, mGrid
               If (Abs(TabAO2(iAO,iGrid,iBfn))<Thr) lGrid = lGrid + 1
            End Do
            If (lGrid==mGrid) lBfn=lBfn+1
            Total = Total + DBLE(lGrid)/DBLE(mGrid)
         End Do
         Total = Total / DBLE(nBfn)
         Write (6,*) 'Sparsity analysis, iAO', iAO
         Write (6,*) ' Total sparcity in %:', 1.0D2 * Total
         Write (6,*) ' Complete Bfn sparcity in %:',
     &                1.0D2*DBLE(lBfn)/DBLE(nBfn)
         Write (6,*)
      End Do
#endif
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Call mma_allocate(A1,mGrid,nBfn,Label='A1')
      Call mma_allocate(A2,mGrid,nBfn,Label='A2')

*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Select Case (Functional_type)
*                                                                      *
      Case (LDA_type)
      Do iD = 1, nD
      Call DCopy_(mGrid*nBfn,TabAO1(1,1,1,iD),nFn,A1,1)
      Call DCopy_(mGrid*nBfn,TabAO2(1,1,1)   ,mAO,A2,1)
      Call DGEMM_('T','N',nBfn,nBfn,mGrid,
     &             One,A1,mGrid,
     &                 A2,mGrid,
     &             Zero,AOInt(1,1,iD),nBfn)
      End Do
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Case (GGA_type)
*                                                                      *
      Do iD = 1, nD
      Call DCopy_(mGrid*nBfn,TabAO1(1,1,1,iD),nFn,A1,1)
      Call DCopy_(mGrid*nBfn,TabAO2(1,1,1)   ,mAO,A2,1)
      Call DGEMM_('T','N',nBfn,nBfn,mGrid,
     &             One,A1,mGrid,
     &                 A2,mGrid,
     &             Zero,AOInt(1,1,iD),nBfn)
      Do iBfn = 1, nBfn
         Do jBfn = 1, iBfn
            AOInt_Sym = AOInt(iBfn,jBfn,iD) +  AOInt(jBfn,iBfn,iD)
            AOInt(iBfn,jBfn,iD) = AOInt_Sym
            AOInt(jBfn,iBfn,iD) = AOInt_Sym
         End Do
      End Do
      End Do
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Case (meta_GGA_type1,meta_GGA_type2)
*                                                                      *
      Do iD = 1, nD
      Call DCopy_(mGrid*nBfn,TabAO1(1,1,1,iD),nFn,A1,1)
      Call DCopy_(mGrid*nBfn,TabAO2(1,1,1)   ,mAO,A2,1)
      Call DGEMM_('T','N',nBfn,nBfn,mGrid,
     &             One,A1,mGrid,
     &                 A2,mGrid,
     &             Zero,AOInt(1,1,iD),nBfn)
      Do iBfn = 1, nBfn
         Do jBfn = 1, iBfn
            AOInt_Sym = AOInt(iBfn,jBfn,iD) +  AOInt(jBfn,iBfn,iD)
            AOInt(iBfn,jBfn,iD) = AOInt_Sym
            AOInt(jBfn,iBfn,iD) = AOInt_Sym
         End Do
      End Do
      End Do
      Call mma_deallocate(A1)
      Call mma_deallocate(A2)

      Call mma_allocate(A1,3*mGrid,nBfn,Label='A1')
      Call mma_allocate(A2,3*mGrid,nBfn,Label='A2')

      Call mma_Allocate(A_tri,nBfn*(nBfn+1)/2,Label='A_tri')
      Do iD = 1, nD
      Call DCopy_(mGrid*nBfn,TabAO1(2,1,1,iD),nFn,A1(1,1),3)
      Call DCopy_(mGrid*nBfn,TabAO2(2,1,1)   ,mAO,A2(1,1),3)
      Call DCopy_(mGrid*nBfn,TabAO1(3,1,1,iD),nFn,A1(2,1),3)
      Call DCopy_(mGrid*nBfn,TabAO2(3,1,1)   ,mAO,A2(2,1),3)
      Call DCopy_(mGrid*nBfn,TabAO1(4,1,1,iD),nFn,A1(3,1),3)
      Call DCopy_(mGrid*nBfn,TabAO2(4,1,1)   ,mAO,A2(3,1),3)
      Call DGEMM_Tri('T','N',nBfn,nBfn,3*mGrid,
     &               One,A1,3*mGrid,
     &                   A2,3*mGrid,
     &               Zero,A_Tri,nBfn)
      ijBfn = 0
      Do iBfn = 1, nBfn
         Do jBfn = 1, iBfn-1
            ijBfn = ijBfn + 1
            AOInt_Sym = A_tri(ijBfn)
            AOInt(iBfn,jBfn,iD) = AOInt(iBfn,jBfn,iD) + AOInt_Sym
            AOInt(jBfn,iBfn,iD) = AOInt(jBfn,iBfn,iD) + AOInt_Sym
         End Do
         ijBfn = ijBfn + 1
         AOInt_Sym = A_tri(ijBfn)
         AOInt(iBfn,iBfn,iD) = AOInt(iBfn,iBfn,iD) + AOInt_Sym
      End Do
      End Do
      Call mma_deAllocate(A_tri)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Case default
*                                                                      *
         Write (6,*) 'DFT_Int: Illegal functional type!'
         Write (6,*) Functional_type
         Call Abend()
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      End Select
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Call mma_deallocate(A1)
      Call mma_deallocate(A2)

*#define _ANALYSIS_
#ifdef _ANALYSIS_
      Write (6,*)
      Write (6,*)  ' Analysing AOInt'
      Thr=1.0D-14
      Do iD = 1, nD
         lBfn = 0
         Do iBfn = 1, nBfn
            Do jBfn = 1, nBfn
               If (Abs(AOInt(iBfn,jBfn,iD))<Thr) lBfn = lBfn + 1
            End Do
         End Do
         Total = DBLE(lBfn)/DBLE(nBfn**2)
         Write (6,*) 'Sparsity analysis, iD', iD
         Write (6,*) ' Total parcity in %:', 1.0D2 * Total
      End Do
#endif
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      End

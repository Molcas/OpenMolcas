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
#define _ALTERNATIVE_CODE_
#ifdef _ALTERNATIVE_CODE_
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
      Real*8, Allocatable:: A1(:,:), A2(:,:)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Call mma_allocate(A1,mGrid,nBfn,Label='A1')
      Call mma_allocate(A2,mGrid,nBfn,Label='A2')
      If (Functional_type.eq.LDA_type) Then
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Do iD = 1, nD
      Call DCopy_(mGrid*nBfn,TabAO1(1,1,1,iD),nFn,A1,1)
      Call DCopy_(mGrid*nBfn,TabAO2(1,1,1)   ,mAO,A2,1)
      Call DGEMM_('T','N',nBfn,nBfn,mGrid,
     &             One,A1,mGrid,
     &                 A2,mGrid,
     &             Zero,AOInt(1,1,iD),nBfn)
      End Do
      Call mma_deallocate(A1)
      Call mma_deallocate(A2)
      Return
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Else If (Functional_type.eq.GGA_type) Then
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Do iD = 1, nD
      Call DCopy_(mGrid*nBfn,TabAO1(1,1,1,iD),nFn,A1,1)
      Call DCopy_(mGrid*nBfn,TabAO2(1,1,1)   ,mAO,A2,1)
      Call DGEMM_('T','N',nBfn,nBfn,mGrid,
     &             One,A1,mGrid,
     &                 A2,mGrid,
     &             Zero,AOInt(1,1,iD),nBfn)
      Call DCopy_(mGrid*nBfn,TabAO1(2,1,1,iD),nFn,A1,1)
      Call DCopy_(mGrid*nBfn,TabAO2(2,1,1)   ,mAO,A2,1)
      Call DGEMM_('T','N',nBfn,nBfn,mGrid,
     &             One,A1,mGrid,
     &                 A2,mGrid,
     &             One ,AOInt(1,1,iD),nBfn)
      Call DCopy_(mGrid*nBfn,TabAO1(3,1,1,iD),nFn,A1,1)
      Call DCopy_(mGrid*nBfn,TabAO2(3,1,1)   ,mAO,A2,1)
      Call DGEMM_('T','N',nBfn,nBfn,mGrid,
     &             One,A1,mGrid,
     &                 A2,mGrid,
     &             One ,AOInt(1,1,iD),nBfn)
      Call DCopy_(mGrid*nBfn,TabAO1(4,1,1,iD),nFn,A1,1)
      Call DCopy_(mGrid*nBfn,TabAO2(4,1,1)   ,mAO,A2,1)
      Call DGEMM_('T','N',nBfn,nBfn,mGrid,
     &             One,A1,mGrid,
     &                 A2,mGrid,
     &             One ,AOInt(1,1,iD),nBfn)
      End Do
      Call mma_deallocate(A1)
      Call mma_deallocate(A2)
      Return
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Else If (Functional_type.eq.meta_GGA_type2) Then
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Do iD = 1, nD
      Call DCopy_(mGrid*nBfn,TabAO1(1,1,1,iD),nFn,A1,1)
      Call DCopy_(mGrid*nBfn,TabAO2(1,1,1)   ,mAO,A2,1)
      Call DGEMM_('T','N',nBfn,nBfn,mGrid,
     &             One,A1,mGrid,
     &                 A2,mGrid,
     &             Zero,AOInt(1,1,iD),nBfn)
      Call DCopy_(mGrid*nBfn,TabAO1(2,1,1,iD),nFn,A1,1)
      Call DCopy_(mGrid*nBfn,TabAO2(2,1,1)   ,mAO,A2,1)
      Call DGEMM_('T','N',nBfn,nBfn,mGrid,
     &             One,A1,mGrid,
     &                 A2,mGrid,
     &             One ,AOInt(1,1,iD),nBfn)
      Call DCopy_(mGrid*nBfn,TabAO1(3,1,1,iD),nFn,A1,1)
      Call DCopy_(mGrid*nBfn,TabAO2(3,1,1)   ,mAO,A2,1)
      Call DGEMM_('T','N',nBfn,nBfn,mGrid,
     &             One,A1,mGrid,
     &                 A2,mGrid,
     &             One ,AOInt(1,1,iD),nBfn)
      Call DCopy_(mGrid*nBfn,TabAO1(4,1,1,iD),nFn,A1,1)
      Call DCopy_(mGrid*nBfn,TabAO2(4,1,1)   ,mAO,A2,1)
      Call DGEMM_('T','N',nBfn,nBfn,mGrid,
     &             One,A1,mGrid,
     &                 A2,mGrid,
     &             One ,AOInt(1,1,iD),nBfn)
      Call DCopy_(mGrid*nBfn,TabAO1(5,1,1,iD),nFn,A1,1)
      Call DCopy_(mGrid*nBfn,TabAO2(5,1,1)   ,mAO,A2,1)
      Call DGEMM_('T','N',nBfn,nBfn,mGrid,
     &             One,A1,mGrid,
     &                 A2,mGrid,
     &             One ,AOInt(1,1,iD),nBfn)
      Call DCopy_(mGrid*nBfn,TabAO2(8,1,1)   ,mAO,A2,1)
      Call DGEMM_('T','N',nBfn,nBfn,mGrid,
     &             One,A1,mGrid,
     &                 A2,mGrid,
     &             One ,AOInt(1,1,iD),nBfn)
      Call DCopy_(mGrid*nBfn,TabAO2(10,1,1)   ,mAO,A2,1)
      Call DGEMM_('T','N',nBfn,nBfn,mGrid,
     &             One,A1,mGrid,
     &                 A2,mGrid,
     &             One ,AOInt(1,1,iD),nBfn)
      End Do
      Call mma_deallocate(A1)
      Call mma_deallocate(A2)
      Return
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Else If (Functional_type.eq.meta_GGA_type1) Then
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Do iD = 1, nD
      Call DCopy_(mGrid*nBfn,TabAO1(1,1,1,iD),nFn,A1,1)
      Call DCopy_(mGrid*nBfn,TabAO2(1,1,1)   ,mAO,A2,1)
      Call DGEMM_('T','N',nBfn,nBfn,mGrid,
     &             One,A1,mGrid,
     &                 A2,mGrid,
     &             Zero,AOInt(1,1,iD),nBfn)
      Call DCopy_(mGrid*nBfn,TabAO1(2,1,1,iD),nFn,A1,1)
      Call DCopy_(mGrid*nBfn,TabAO2(2,1,1)   ,mAO,A2,1)
      Call DGEMM_('T','N',nBfn,nBfn,mGrid,
     &             One,A1,mGrid,
     &                 A2,mGrid,
     &             One ,AOInt(1,1,iD),nBfn)
      Call DCopy_(mGrid*nBfn,TabAO1(3,1,1,iD),nFn,A1,1)
      Call DCopy_(mGrid*nBfn,TabAO2(3,1,1)   ,mAO,A2,1)
      Call DGEMM_('T','N',nBfn,nBfn,mGrid,
     &             One,A1,mGrid,
     &                 A2,mGrid,
     &             One ,AOInt(1,1,iD),nBfn)
      Call DCopy_(mGrid*nBfn,TabAO1(4,1,1,iD),nFn,A1,1)
      Call DCopy_(mGrid*nBfn,TabAO2(4,1,1)   ,mAO,A2,1)
      Call DGEMM_('T','N',nBfn,nBfn,mGrid,
     &             One,A1,mGrid,
     &                 A2,mGrid,
     &             One ,AOInt(1,1,iD),nBfn)
      End Do
      Call mma_deallocate(A1)
      Call mma_deallocate(A2)
      Return
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Else
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
         Write (6,*) 'DFT_Int: Illegal functional type!'
         Call Abend()
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      End
#else
      Subroutine Do_NIntX(AOInt,mGrid,TabAO1,TabAO2,nBfn,nD,mAO,nFn)
*                                                                      *
************************************************************************
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_info.fh"
      Real*8 AOInt(nBfn,nBfn,nD), TabAO1(nFn,mGrid,nBfn,nD),
     &       TabAO2(mAO,mGrid,nBfn)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      If (Functional_type.eq.LDA_type) Then
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Do iD = 1, nD
      Do iCB = 1, nBfn
*                                                                      *
************************************************************************
*                                                                      *
         Do jCB = 1, iCB
*                                                                      *
************************************************************************
*                                                                      *
            ToAdd1=Zero
            Do iGrid = 1, mGrid
               ToAdd1 = ToAdd1
     &                + TabAO1(1,iGrid,iCB,iD)*TabAO2(1,iGrid,jCB)
            End Do
            AOInt(iCB,jCB,iD) = ToAdd1
            AOInt(jCB,iCB,iD) = ToAdd1
*                                                                      *
************************************************************************
*                                                                      *
         End Do
      End Do
      End Do
      Return
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Else If (Functional_type.eq.GGA_type) Then
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Do iD = 1, nD
      Do iCB = 1, nBfn
*                                                                      *
************************************************************************
*                                                                      *
*
         Do jCB = 1, iCB
*                                                                      *
************************************************************************
*                                                                      *
            ToAdd1=Zero
*
            Do iGrid = 1, mGrid
               ToAdd1= ToAdd1
     &               + TabAO1(1,iGrid,iCB,iD) * TabAO2(1,iGrid,jCB)
     &               + TabAO1(2,iGrid,iCB,iD) * TabAO2(2,iGrid,jCB)
     &               + TabAO1(3,iGrid,iCB,iD) * TabAO2(3,iGrid,jCB)
     &               + TabAO1(4,iGrid,iCB,iD) * TabAO2(4,iGrid,jCB)
            End Do
            AOInt(iCB,jCB,iD) = ToAdd1
            AOInt(jCB,iCB,iD) = ToAdd1
*                                                                      *
************************************************************************
*                                                                      *
         End Do
      End Do
      End Do
      Return
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Else If (Functional_type.eq.meta_GGA_type2) Then
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Do iD = 1, nD
      Do iCB = 1, nBfn
*                                                                      *
************************************************************************
*                                                                      *
*
         Do jCB = 1, iCB
*                                                                      *
************************************************************************
*                                                                      *
            ToAdd1=Zero
*
            Do iGrid = 1, mGrid
               ToAdd1= ToAdd1
     &               + TabAO1(1,iGrid,iCB,iD) * TabAO2(1,iGrid,jCB)
     &               + TabAO1(2,iGrid,iCB,iD) * TabAO2(2,iGrid,jCB)
     &               + TabAO1(3,iGrid,iCB,iD) * TabAO2(3,iGrid,jCB)
     &               + TabAO1(4,iGrid,iCB,iD) * TabAO2(4,iGrid,jCB)
     &               + TabAO1(5,iGrid,iCB,iD) *(TabAO2(5,iGrid,jCB)
     &                                         +TabAO2(8,iGrid,jCB)
     &                                        +TabAO2(10,iGrid,jCB))
            End Do
            AOInt(iCB,jCB,iD) = ToAdd1
            AOInt(jCB,iCB,iD) = ToAdd1
*                                                                      *
************************************************************************
*                                                                      *
         End Do
      End Do
      End Do
      Return
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Else If (Functional_type.eq.meta_GGA_type1) Then
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Do iD = 1, nD
      Do iCB = 1, nBfn
*                                                                      *
************************************************************************
*                                                                      *
*
         Do jCB = 1, iCB
*                                                                      *
************************************************************************
*                                                                      *
            ToAdd1=Zero
*
            Do iGrid = 1, mGrid
               ToAdd1= ToAdd1
     &               + TabAO1(1,iGrid,iCB,iD) * TabAO2(1,iGrid,jCB)
     &               + TabAO1(2,iGrid,iCB,iD) * TabAO2(2,iGrid,jCB)
     &               + TabAO1(3,iGrid,iCB,iD) * TabAO2(3,iGrid,jCB)
     &               + TabAO1(4,iGrid,iCB,iD) * TabAO2(4,iGrid,jCB)
            End Do
            AOInt(iCB,jCB,iD) = ToAdd1
            AOInt(jCB,iCB,iD) = ToAdd1
*                                                                      *
************************************************************************
*                                                                      *
         End Do
      End Do
      End Do
      Return
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Else
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
         Write (6,*) 'DFT_Int: Illegal functional type!'
         Call Abend()
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      End
#endif

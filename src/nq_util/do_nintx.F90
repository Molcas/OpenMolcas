!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1991,2021,2022, Roland Lindh                           *
!***********************************************************************
!***********************************************************************
!                                                                      *
      Subroutine Do_NIntX()
!                                                                      *
!***********************************************************************
!***********************************************************************
      use nq_Grid, only: Grid_AO, TabAO
      use nq_Grid, only: AOInt => Dens_AO
      use nq_Grid, only: iBfn_Index
      use nq_Info
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "stdalloc.fh"
      Real*8, Allocatable:: A1(:,:), A2(:,:), A_tri(:)
      Real*8, Allocatable:: A3(:,:,:), A4(:,:,:)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
      mGrid=SIZE(TabAO,2)
      nBfn =SIZE(iBfn_Index,2)
      nD   =SIZE(Grid_AO,4)

!#define _ANALYSIS_
#ifdef _ANALYSIS_
      mAO  =SIZE(TabAO,1)
      nFn  =SIZE(Grid_AO,1)
      Write (6,*)
      Write (6,*)  ' Analysing Grid_AO'
      Thr=1.0D-14
      Do iD = 1, nD
      Do iFn = 1, nFn
         lBfn = 0
         Total=0.0D0
         Do iBfn = 1, nBfn
            lGrid = 0
            Do iGrid = 1, mGrid
               If (Abs(Grid_AO(iFn,iGrid,iBfn,iD))<Thr) lGrid = lGrid + 1
            End Do
            If (lGrid==mGrid) lBfn=lBfn+1
            Total = Total + DBLE(lGrid)/DBLE(mGrid)
         End Do
         Total = Total / DBLE(nBfn)
         Write (6,*) 'Sparsity analysis, iD, iFn', iD, iFn
         Write (6,*) ' Total sparcity in %:', 1.0D2 * Total
         Write (6,*) ' Complete Bfn sparcity in %:',                    &
     &                1.0D2*DBLE(lBfn)/DBLE(nBfn)
         Write (6,*)
      End Do
      End Do
      Write (6,*)
      Write (6,*)  ' Analysing TabAO'
      Do iAO = 1, mAO
         lBfn = 0
         Total=0.0D0
         Do iBfn = 1, nBfn
            lGrid = 0
            Do iGrid = 1, mGrid
               If (Abs(TabAO(iAO,iGrid,iBfn))<Thr) lGrid = lGrid + 1
            End Do
            If (lGrid==mGrid) lBfn=lBfn+1
            Total = Total + DBLE(lGrid)/DBLE(mGrid)
         End Do
         Total = Total / DBLE(nBfn)
         Write (6,*) 'Sparsity analysis, iAO', iAO
         Write (6,*) ' Total sparcity in %:', 1.0D2 * Total
         Write (6,*) ' Complete Bfn sparcity in %:',                    &
     &                1.0D2*DBLE(lBfn)/DBLE(nBfn)
         Write (6,*)
      End Do
#endif
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
      Call mma_allocate(A1,mGrid,nBfn,Label='A1')
      Call mma_allocate(A2,mGrid,nBfn,Label='A2')
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
      Select Case (Functional_type)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
      Case (LDA_type)
!                                                                      *
      Call mma_Allocate(A_tri,nBfn*(nBfn+1)/2,Label='A_tri')
      AOInt(:,:,:)=Zero
      A2(1:mGrid,1:nBfn)=TabAO(1,1:mGrid,1:nBfn)
      Do iD = 1, nD
      A1(1:mGrid,1:nBfn)=Grid_AO(1,1:mGrid,1:nBfn,iD)
      Call DGEMM_Tri('T','N',nBfn,nBfn,mGrid,                           &
     &              One,A1,mGrid,                                       &
     &                  A2,mGrid,                                       &
     &              Zero,A_Tri,nBfn)
      Call Sym_Dist()
      End Do
      Call mma_deAllocate(A_tri)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
      Case (GGA_type)
!                                                                      *
      A2(1:mGrid,1:nBfn)=TabAO(1,1:mGrid,1:nBfn)
      Do iD = 1, nD
      A1(1:mGrid,1:nBfn)=Grid_AO(1,1:mGrid,1:nBfn,iD)
      Call DGEMM_('T','N',nBfn,nBfn,mGrid,                              &
     &             One,A1,mGrid,                                        &
     &                 A2,mGrid,                                        &
     &             Zero,AOInt(1,1,iD),nBfn)
      Call Symmetrize()
      End Do
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
      Case (meta_GGA_type1,meta_GGA_type2)
!                                                                      *
      A2(1:mGrid,1:nBfn)=TabAO(1,1:mGrid,1:nBfn)
      Do iD = 1, nD
      A1(1:mGrid,1:nBfn)=Grid_AO(1,1:mGrid,1:nBfn,iD)
      Call DGEMM_('T','N',nBfn,nBfn,mGrid,                              &
     &             One,A1,mGrid,                                        &
     &                 A2,mGrid,                                        &
     &             Zero,AOInt(1,1,iD),nBfn)
      Call Symmetrize()
      End Do

      Call mma_allocate(A3,3,mGrid,nBfn,Label='A1')
      Call mma_allocate(A4,3,mGrid,nBfn,Label='A2')

      Call mma_Allocate(A_tri,nBfn*(nBfn+1)/2,Label='A_tri')
      A4(1:3,1:mGrid,1:nBfn)=TabAO(2:4,1:mGrid,1:nBfn)
      Do iD = 1, nD
      A3(1:3,1:mGrid,1:nBfn)=Grid_AO(2:4,1:mGrid,1:nBfn,iD)
      Call DGEMM_Tri('T','N',nBfn,nBfn,3*mGrid,                         &
     &               One,A3,3*mGrid,                                    &
     &                   A4,3*mGrid,                                    &
     &               Zero,A_Tri,nBfn)
      Call Sym_Dist()
      End Do
      Call mma_deallocate(A3)
      Call mma_deallocate(A4)
      Call mma_deAllocate(A_tri)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
      Case default
!                                                                      *
         Write (6,*) 'DFT_Int: Illegal functional type!'
         Write (6,*) Functional_type
         Call Abend()
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
      End Select
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
      Call mma_deallocate(A1)
      Call mma_deallocate(A2)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
!#define _ANALYSIS_
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
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
      Contains
        Subroutine Sym_Dist()
           Integer iBfn, jBfn, ijBfn
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
        End Subroutine Sym_Dist
        Subroutine Symmetrize()
           Integer iBfn, jBfn
           Do iBfn = 1, nBfn
              Do jBfn = 1, iBfn
                 AOInt_Sym = AOInt(iBfn,jBfn,iD) +  AOInt(jBfn,iBfn,iD)
                 AOInt(iBfn,jBfn,iD) = AOInt_Sym
                 AOInt(jBfn,iBfn,iD) = AOInt_Sym
              End Do
           End Do
        End Subroutine Symmetrize
      End

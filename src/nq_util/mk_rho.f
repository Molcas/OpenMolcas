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
* Copyright (C) 2000,2021, Roland Lindh                                *
************************************************************************
      Subroutine Mk_Rho(list_s,nlist_s,Fact,mdc,list_bas,Index,nIndex,
     &                  Do_Grad)
************************************************************************
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN.  2000                                   *
************************************************************************
      use iSD_data
      use k2_arrays, only: DeDe, ipDijS
      use nq_grid, only: Rho, TabAO, Dens_AO, Grid_AO, TabAO_Short
      use nq_grid, only: GradRho, Tau, Lapl, kAO
      use nq_Grid, only: dRho_dR, iBfn_Index
      use nq_Grid, only: List_G
#ifdef _DEBUGPRINT_
      use nq_grid, only: nRho
#endif
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
#include "debug.fh"
#include "nq_info.fh"
#include "nsd.fh"
#include "setup.fh"
#include "stdalloc.fh"
      Integer Index(nIndex)
      Real*8 Fact(mdc**2)
      Integer ipD(2)
      Integer list_s(2,nlist_s), list_bas(2,nlist_s)
      Integer, Parameter :: Index_d2(3,3)=
     &    Reshape([5,6,7, 6,8,9, 7,9,10],[3,3])
      Integer, Parameter :: Index_d3(3,3) =
     &    Reshape([11,14,16, 12,17,19, 13,18,20],[3,3])
      Logical Do_Grad
      Integer, Allocatable:: Ind_Grd(:,:)
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
*
      nD = SIZE(Dens_AO,3)
      nAO = SIZE(Dens_AO,1)
      Dens_AO(:,:,:)=Zero
      mAO = SIZE(TabAO,1)
      mGrid = SIZE(TabAO,2)

*#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Write (6,*) 'mAO=',mAO
      Write (6,*) 'mGrid=',mGrid
      Write (6,*) 'nlist_s=',nlist_s
      Call RecPrt('Rho: TabAO',' ',TabAO,mAO*mGrid,nAO)
#endif
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Generate the one-particle density matrix, D(mu,nu)               *
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*
      If (Do_Grad) Then
         Call mma_Allocate(Ind_Grd,3,nAO,Label='Ind_Grd')
         Ind_Grd(:,:)=0
      End If

      nBfn=SIZE(iBfn_Index,2)
      If (nBfn/=nAO) Then
         Write (6,*) 'mk_Rho: internal error!'
         Call Abend()
      End If
      Factor = DBLE(2/nD)
      Do iBfn = 1, nBfn
         ilist_s=iBfn_Index(2,iBfn)
         i1     =iBfn_Index(3,iBfn)
         i2     =iBfn_Index(4,iBfn)
         iSkal = list_s(1,ilist_s)
         kDCRE = list_s(2,ilist_s)
         iCmp  = iSD( 2,iSkal)
         iBas  = iSD( 3,iSkal)
         mdci  = iSD(10,iSkal)
         iShell= iSD(11,iSkal)
         iBas_Eff=list_bas(1,ilist_s)
         index_i =list_bas(2,ilist_s)
         nFunc_i=iBas*iCmp

         i_R=(i1-1)*iBas_Eff+i2
         iCB = Index(index_i-1+i_R)

         If (Do_Grad) Ind_Grd(:,iBfn)=List_g(:,ilist_s)

         Do jBfn = 1, iBfn
            jlist_s=iBfn_Index(2,jBfn)
            j1     =iBfn_Index(3,jBfn)
            j2     =iBfn_Index(4,jBfn)
            jSkal = list_s(1,jlist_s)
            kDCRR = list_s(2,jlist_s)
            jCmp  = iSD( 2,jSkal)
            jBas  = iSD( 3,jSkal)
            mdcj  = iSD(10,jSkal)
            jShell= iSD(11,jSkal)
            jBas_Eff=list_bas(1,jlist_s)
            index_j =list_bas(2,jlist_s)
            nFunc_j=jBas*jCmp

            j_R=(j1-1)*jBas_Eff+j2
            jCB = Index(index_j-1+j_R)

            ijS=iTri(iShell,jShell)
            ip_Tmp=ipDijs
            Call Dens_Info(ijS,ipDij,ipDSij,mDCRij,ipDDij,ip_Tmp,nD)
            ij = (mdcj-1)*mdc + mdci

            iER=iEOr(kDCRE,kDCRR)
            lDCRER=NrOpr(iER)

            mDij=nFunc_i*nFunc_j
            ip_D_a=ipDij+lDCRER*mDij
            ip_D_b=ip_D_a
            If (nD.ne.1) ip_D_b=ipDSij+lDCRER*mDij
            ipD(1)=ip_D_a
            ipD(2)=ip_D_b

            ij_D = (jCB-1)*nFunc_i + iCB - 1
            Do iD = 1, nD
               DAij =DeDe(ipD(iD)+ij_D)*Fact(ij)*Factor
               Dens_AO(iBfn,jBfn,iD) = DAij
               Dens_AO(jBfn,iBfn,iD) = DAij
            End Do

         End Do

      End Do
*                                                                      *
************************************************************************
*                                                                      *
*#define _ANALYSIS_
#ifdef _ANALYSIS_
      Thr=1.0D-15
      Write (6,*)
      Write (6,*)  ' Sparsity analysis of D(i,j)'
      Write (6,*)  ' Threshold: ',Thr
      Write (6,*)  ' Grid size: ',mGrid
      Write (6,*)  ' Dimension: ',n,' x ',n
      n=SIZE(Dens_AO,1)
      n2 = n**2
      Do iD = 1, nD
      m=0
      Do i = 1, n
         Do j = 1, n
            If (Abs(Dens_AO(i,j,iD))<Thr) m=m+1
         End Do
      End Do
      Write (6,*) 'Total Sparsity in %', 1.0D2*DBLE(m)/DBLE(n2)
      k=0
      Do i = 1, n
         m = 0
         Do j = 1, n
            If (Abs(Dens_AO(i,j,iD))<Thr) m=m+1
         End Do
         If (m==n) k=k+1
      End Do
      Write (6,*) 'Column Sparsity in %', 1.0D2*DBLE(k)/DBLE(n)
      k=0
      Do j = 1, n
         m = 0
         Do i = 1, n
            If (Abs(Dens_AO(i,j,iD))<Thr) m=m+1
         End Do
         If (m==n) k=k+1
      End Do
      Write (6,*) 'Row Sparsity in %', 1.0D2*DBLE(k)/DBLE(n)
      End Do
#endif
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Construct: Sum_i D_ij TabAO(i,iGrid,iAO)
*          D_ij is the one-electron density
*          TabAO(i,iGrid,iAO) are values with respect to the ith AO
*          i=1 is the value of the AO
*          i=2-4 are the values of the first order derivatives
*          i=5-10 are the values of the second order derivatives
*          i=11-20 are the values of the third order derivatives
*
*     During a gradient calculation the size of the fast index of
*     TabAO is larger than that of Grid_AO. In those cases we copy
*     the part of TabAO which we need to TabAO_Short before we make the
*     contraction with the 1-particle density matrix.
*
      If (Do_Grad) Then
         TabAO_Short(1:kAO,1:mGrid,:) = TabAO(1:kAO,1:mGrid,:)
         Call DGEMM_('N','N',kAO*mGrid,nAO*nD,nAO,
     &               One,TabAO_Short,kAO*mGrid,
     &                   Dens_AO,nAO,
     &               Zero,Grid_AO,kAO*mGrid)
      Else
         Call DGEMM_('N','N',kAO*mGrid,nAO*nD,nAO,
     &               One,TabAO,mAO*mGrid,
     &                   Dens_AO,nAO,
     &               Zero,Grid_AO,kAO*mGrid)
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      If (Allocated(dRho_dR))  dRho_dR(:,:,:)=Zero

      Select Case(Functional_Type)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Case (LDA_Type)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
         Rho(:,1:mGrid)=Zero
         Do iD = 1, nD
            Do iAO = 1, nAO

               Do iGrid = 1, mGrid
                  Rho(iD,iGrid) = Rho(iD,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
               End Do

               If (Do_Grad) Then
*
*------------- Loop over cartesian components
*
               Do iCar = 1, 3

                  Ind_xyz=Ind_Grd(iCar,iAO)
                  j = iCar + 1

                  If (Ind_xyz/=0) Then
                     Do iGrid = 1, mGrid
*
*                       Cartesian derivative of the density.
*
                        dRho_dR(iD,iGrid,Ind_xyz)
     &                                = dRho_dR(iD,iGrid,Ind_xyz)
     &                                + Two * Grid_AO(1,iGrid,iAO,iD)
     &                                * TabAO(j,iGrid,iAO)
                     End Do
                  End If

               End Do

               End If
            End Do
         End Do
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Case (GGA_Type)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
         Rho(:,1:mGrid)=Zero
         GradRho(:,1:mGrid)=Zero
         Do iD = 1, nD
            ix = (iD-1)*3 + 1
            iy = (iD-1)*3 + 2
            iz = (iD-1)*3 + 3
            Do iAO = 1, nAO

               Do iGrid = 1, mGrid
                  Rho(iD,iGrid) = Rho(iD,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
                  GradRho(ix,iGrid)=GradRho(ix,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(2,iGrid,iAO)
     &                          + Grid_AO(2,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
                  GradRho(iy,iGrid)=GradRho(iy,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(3,iGrid,iAO)
     &                          + Grid_AO(3,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
                  GradRho(iz,iGrid)=GradRho(iz,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(4,iGrid,iAO)
     &                          + Grid_AO(4,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
               End Do

               If (Do_Grad) Then
*
*------------- Loop over cartesian components
*
               Do iCar = 1, 3

                  Ind_xyz=Ind_Grd(iCar,iAO)! index of  nuclear gradient

                  j = iCar + 1             ! index derivative of AO

                  iDx = nD + (iD-1)*3 + 1  ! index of grad rho component
                  iDy = iDx + 1
                  iDz = iDy + 1

                  idjx = Index_d2(1,iCar)
                  idjy = Index_d2(2,iCar)
                  idjz = Index_d2(3,iCar)
                  If (Ind_xyz/=0) Then
                     Do iGrid = 1, mGrid
*
*                       Cartesian derivative of rho
*
                        dRho_dR(iD,iGrid,Ind_xyz)
     &                             = dRho_dR(iD,iGrid,Ind_xyz)
     &                        + Two * Grid_AO(1,iGrid,iAO,iD)
     &                                * TabAO(j,iGrid,iAO)
*
*                       Cartesian derivatives of grad rho
*
                        dRho_dR(iDx,iGrid,Ind_xyz)
     &                              = dRho_dR(iDx,iGrid,Ind_xyz)
     &                         + Two * TabAO(idjx,iGrid,iAO)
     &                                * Grid_AO(1,iGrid,iAO,iD)
     &                            + Two * TabAO(j,iGrid,iAO)
     &                                * Grid_AO(2,iGrid,iAO,iD)
                        dRho_dR(iDy,iGrid,Ind_xyz)
     &                              = dRho_dR(iDy,iGrid,Ind_xyz)
     &                         + Two * TabAO(idjy,iGrid,iAO)
     &                                * Grid_AO(1,iGrid,iAO,iD)
     &                            + Two * TabAO(j,iGrid,iAO)
     &                                * Grid_AO(3,iGrid,iAO,iD)
                        dRho_dR(iDz,iGrid,Ind_xyz)
     &                              = dRho_dR(iDz,iGrid,Ind_xyz)
     &                         + Two * TabAO(idjz,iGrid,iAO)
     &                                * Grid_AO(1,iGrid,iAO,iD)
     &                            + Two * TabAO(j,iGrid,iAO)
     &                                * Grid_AO(4,iGrid,iAO,iD)
                     End Do
                  End If

               End Do
               End If

            End Do
         End Do
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Case (meta_GGA_Type1)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
         Rho(:,1:mGrid)=Zero
         GradRho(:,1:mGrid)=Zero
         Tau(:,1:mGrid)=Zero
         Do iD = 1, nD
            ix = (iD-1)*3 + 1
            iy = (iD-1)*3 + 2
            iz = (iD-1)*3 + 3
            Do iAO = 1, nAO

               Do iGrid = 1, mGrid
                  Rho(iD,iGrid) = Rho(iD,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
                  GradRho(ix,iGrid)=GradRho(ix,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(2,iGrid,iAO)
     &                          + Grid_AO(2,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
                  GradRho(iy,iGrid)=GradRho(iy,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(3,iGrid,iAO)
     &                          + Grid_AO(3,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
                  GradRho(iz,iGrid)=GradRho(iz,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(4,iGrid,iAO)
     &                          + Grid_AO(4,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
                  Tau(iD,iGrid) = Tau(iD,iGrid)
     &                          + Grid_AO(2,iGrid,iAO,iD)
     &                          * TabAO(2,iGrid,iAO)
     &                          + Grid_AO(3,iGrid,iAO,iD)
     &                          * TabAO(3,iGrid,iAO)
     &                          + Grid_AO(4,iGrid,iAO,iD)
     &                          * TabAO(4,iGrid,iAO)
               End Do

               If (Do_Grad) Then
*
*------------- Loop over cartesian components
*
               Do iCar = 1, 3

                  Ind_xyz=Ind_Grd(iCar,iAO)! index of  nuclear gradient

                  j = iCar + 1             ! index derivative of AO

                  iDx = nD + (iD-1)*3 + 1  ! index of grad rho component
                  iDy = iDx + 1
                  iDz = iDy + 1

                  iT  = nD*4 + iD      ! index of tau component

                  idjx = Index_d2(1,iCar)
                  idjy = Index_d2(2,iCar)
                  idjz = Index_d2(3,iCar)
                  If (Ind_xyz/=0) Then
                     Do iGrid = 1, mGrid
*
*                       Cartesian derivative of rho
*
                        dRho_dR(iD,iGrid,Ind_xyz)
     &                             = dRho_dR(iD,iGrid,Ind_xyz)
     &                        + Two * Grid_AO(1,iGrid,iAO,iD)
     &                                * TabAO(j,iGrid,iAO)
*
*                       Cartesian derivatives of grad rho
*
                        dRho_dR(iDx,iGrid,Ind_xyz)
     &                              = dRho_dR(iDx,iGrid,Ind_xyz)
     &                         + Two * TabAO(idjx,iGrid,iAO)
     &                                * Grid_AO(1,iGrid,iAO,iD)
     &                            + Two * TabAO(j,iGrid,iAO)
     &                                * Grid_AO(2,iGrid,iAO,iD)
                        dRho_dR(iDy,iGrid,Ind_xyz)
     &                              = dRho_dR(iDy,iGrid,Ind_xyz)
     &                         + Two * TabAO(idjy,iGrid,iAO)
     &                                * Grid_AO(1,iGrid,iAO,iD)
     &                            + Two * TabAO(j,iGrid,iAO)
     &                                * Grid_AO(3,iGrid,iAO,iD)
                        dRho_dR(iDz,iGrid,Ind_xyz)
     &                              = dRho_dR(iDz,iGrid,Ind_xyz)
     &                         + Two * TabAO(idjz,iGrid,iAO)
     &                                * Grid_AO(1,iGrid,iAO,iD)
     &                            + Two * TabAO(j,iGrid,iAO)
     &                                * Grid_AO(4,iGrid,iAO,iD)
*
*                       Cartesian derivatives of tau
*
                        dRho_dR(iT,iGrid,Ind_xyz)
     &                             = dRho_dR(iT,iGrid,Ind_xyz)
     &                       + Four* TabAO(idjx,iGrid,iAO)
     &                              * Grid_AO(2,iGrid,iAO,iD)
     &                       + Four* TabAO(idjy,iGrid,iAO)
     &                              * Grid_AO(3,iGrid,iAO,iD)
     &                       + Four* TabAO(idjz,iGrid,iAO)
     &                              * Grid_AO(4,iGrid,iAO,iD)
                     End Do
                  End If

               End Do

               End If
            End Do
         End Do
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Case (meta_GGA_Type2)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
         Rho(:,1:mGrid)=Zero
         GradRho(:,1:mGrid)=Zero
         Tau(:,1:mGrid)=Zero
         Lapl(:,1:mGrid)=Zero
         Do iD = 1, nD
            ix = (iD-1)*3 + 1
            iy = (iD-1)*3 + 2
            iz = (iD-1)*3 + 3
            Do iAO = 1, nAO

               Do iGrid = 1, mGrid
                  Rho(iD,iGrid) = Rho(iD,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
                  GradRho(ix,iGrid)=GradRho(ix,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(2,iGrid,iAO)
     &                          + Grid_AO(2,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
                  GradRho(iy,iGrid)=GradRho(iy,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(3,iGrid,iAO)
     &                          + Grid_AO(3,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
                  GradRho(iz,iGrid)=GradRho(iz,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(4,iGrid,iAO)
     &                          + Grid_AO(4,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
                  Tau(iD,iGrid) = Tau(iD,iGrid)
     &                          + Grid_AO(2,iGrid,iAO,iD)
     &                          * TabAO(2,iGrid,iAO)
     &                          + Grid_AO(3,iGrid,iAO,iD)
     &                          * TabAO(3,iGrid,iAO)
     &                          + Grid_AO(4,iGrid,iAO,iD)
     &                          * TabAO(4,iGrid,iAO)
                  Lapl(iD,iGrid) = Lapl(iD,iGrid)
     &                          + TabAO( 1,iGrid,iAO) *
     &                          ( Grid_AO( 5,iGrid,iAO,iD)
     &                          + Grid_AO( 8,iGrid,iAO,iD)
     &                          + Grid_AO(10,iGrid,iAO,iD) )
     &                          + Two *
     &                          ( Grid_AO(2,iGrid,iAO,iD)
     &                          * TabAO(2,iGrid,iAO)
     &                          + Grid_AO(3,iGrid,iAO,iD)
     &                          * TabAO(3,iGrid,iAO)
     &                          + Grid_AO(4,iGrid,iAO,iD)
     &                          * TabAO(4,iGrid,iAO) )
     &                          +(TabAO( 5,iGrid,iAO)
     &                          + TabAO( 8,iGrid,iAO)
     &                          + TabAO(10,iGrid,iAO) )
     &                          * Grid_AO(1,iGrid,iAO,iD)
               End Do

               If (Do_Grad) Then
*
*------------- Loop over cartesian components
*
               Do iCar = 1, 3

                  Ind_xyz=Ind_Grd(iCar,iAO)! index of  nuclear gradient

                  j = iCar + 1             ! index derivative of AO

                  iDx = nD + (iD-1)*3 + 1  ! index of grad rho component
                  iDy = iDx + 1
                  iDz = iDy + 1

                  iT  = nD*4 + iD      ! index of tau component

                  iL  = nD*5 + iD      ! index if laplacian component

                  idjx = Index_d2(1,iCar)
                  idjy = Index_d2(2,iCar)
                  idjz = Index_d2(3,iCar)

                  idjx2 = Index_d3(1,iCar)
                  idjy2 = Index_d3(2,iCar)
                  idjz2 = Index_d3(3,iCar)
                  idx2  = Index_d2(1,1)
                  idy2  = Index_d2(2,2)
                  idz2  = Index_d2(3,3)
                  If (Ind_xyz/=0) Then
                     Do iGrid = 1, mGrid
*
*                       Cartesian derivative of rho
*
                        dRho_dR(iD,iGrid,Ind_xyz)
     &                             = dRho_dR(iD,iGrid,Ind_xyz)
     &                        + Two * Grid_AO(1,iGrid,iAO,iD)
     &                                * TabAO(j,iGrid,iAO)
*
*                       Cartesian derivatives of grad rho
*
                        dRho_dR(iDx,iGrid,Ind_xyz)
     &                              = dRho_dR(iDx,iGrid,Ind_xyz)
     &                         + Two * TabAO(idjx,iGrid,iAO)
     &                                * Grid_AO(1,iGrid,iAO,iD)
     &                            + Two * TabAO(j,iGrid,iAO)
     &                                * Grid_AO(2,iGrid,iAO,iD)
                        dRho_dR(iDy,iGrid,Ind_xyz)
     &                              = dRho_dR(iDy,iGrid,Ind_xyz)
     &                         + Two * TabAO(idjy,iGrid,iAO)
     &                                * Grid_AO(1,iGrid,iAO,iD)
     &                            + Two * TabAO(j,iGrid,iAO)
     &                                * Grid_AO(3,iGrid,iAO,iD)
                        dRho_dR(iDz,iGrid,Ind_xyz)
     &                              = dRho_dR(iDz,iGrid,Ind_xyz)
     &                         + Two * TabAO(idjz,iGrid,iAO)
     &                                * Grid_AO(1,iGrid,iAO,iD)
     &                            + Two * TabAO(j,iGrid,iAO)
     &                                * Grid_AO(4,iGrid,iAO,iD)
*
*                       Cartesian derivatives of tau
*
                        dRho_dR(iT,iGrid,Ind_xyz)
     &                             = dRho_dR(iT,iGrid,Ind_xyz)
     &                       + Four* TabAO(idjx,iGrid,iAO)
     &                              * Grid_AO(2,iGrid,iAO,iD)
     &                       + Four* TabAO(idjy,iGrid,iAO)
     &                              * Grid_AO(3,iGrid,iAO,iD)
     &                       + Four* TabAO(idjz,iGrid,iAO)
     &                              * Grid_AO(4,iGrid,iAO,iD)
*
*                       Cartesian derivatives of the laplacian
*
                        dRho_dR(iL,iGrid,Ind_xyz)
     &                             = dRho_dR(iL,iGrid,Ind_xyz)

     &                             + Two * Grid_AO(1,iGrid,iAO,iD)
     &                             *   ( TabAO(idjx2,iGrid,iAO)
     &                                  +TabAO(idjy2,iGrid,iAO)
     &                                  +TabAO(idjz2,iGrid,iAO))

     &                             + Two *(Grid_AO(idx2,iGrid,iAO,iD)
     &                                    +Grid_AO(idy2,iGrid,iAO,iD)
     &                                    +Grid_AO(idz2,iGrid,iAO,iD))
     &                                    *TabAO(j,iGrid,iAO)

     &                             + Four*(Grid_AO(2,iGrid,iAO,iD)
     &                                    *TabAO(idjx,iGrid,iAO)
     &                                    +Grid_AO(3,iGrid,iAO,iD)
     &                                    *TabAO(idjy,iGrid,iAO)
     &                                    +Grid_AO(4,iGrid,iAO,iD)
     &                                    *TabAO(idjz,iGrid,iAO))

                     End Do
                  End If

               End Do

             End If

            End Do
         End Do
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Case default
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
         Call abend()
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      End Select
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
#ifdef _ANALYSIS_
      Write (6,*)
      Write (6,*) 'Rho Sparsity analysis'
      n=0
      Do iGrid = 1, mGrid
         tmp = Zero
         Do iD = 1, nD
            tmp = tmp + Rho(iD,iGrid)
         End Do
         If (tmp<Thr) n=n+1
      End Do
      Write (6,*) 'Rho Sparsity in %: ',1.0D2*DBLE(n)/DBLE(mGrid)
#endif
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
       If (Allocated(Tau)) Tau(:,1:mGrid)=Half*Tau(:,1:mGrid)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      Do iD = 1, nD
         Call RecPrt('Dens_AO',' ',Dens_AO(:,:,iD),nAO,nAO)
         Call RecPrt('Grid_AO',' ',Grid_AO(:,:,:,iD),mAO*mGrid,
     &                                               nAO)
      End Do
      If (Do_Grad) Then
         nGrad_Eff = SIZE(dRho_dR,3)
         Call RecPrt('dRho_dR_LDA: dRho_dR',' ',dRho_dR,
     &                       SIZE(dRho_dR,1)*mGrid,nGrad_Eff)
      End If
#endif
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      If (Allocated(Ind_grd))   Call mma_deAllocate(Ind_Grd)
      Return
      End

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
! Copyright (C) 2000,2021, Roland Lindh                                *
!***********************************************************************

subroutine Mk_Rho(list_s,nlist_s,Fact,mdc,list_bas,Indx,nIndex,Do_Grad)
!***********************************************************************
!      Author:Roland Lindh, Department of Chemical Physics, University *
!             of Lund, SWEDEN.  2000                                   *
!***********************************************************************

use iSD_data, only: iSD
use k2_arrays, only: DeDe, ipDijS
use nq_Grid, only: Dens_AO, dRho_dR, GradRho, Grid_AO, iBfn_Index, kAO, Lapl, List_G, Rho, TabAO, TabAO_Short, Tau
use nq_Info, only: Functional_type, GGA_Type, LDA_Type, meta_GGA_Type1, meta_GGA_Type2
use Index_Functions, only: iTri
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Four, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nlist_s, list_s(2,nlist_s), mdc, list_bas(2,nlist_s), nIndex, Indx(nIndex)
real(kind=wp), intent(in) :: Fact(mdc,mdc)
logical(kind=iwp), intent(in) :: Do_Grad
integer(kind=iwp) :: i1, i2, i_R, iAO, iBas, iBas_Eff, iBfn, iCar, iCB, iCmp, iD, idjx, idjx2, idjy, idjy2, idjz, idjz2, iDx, &
                     idx2, iDy, idy2, iDz, idz2, iER, iGrid, ij_D, ijS, iL, ilist_s, Ind_xyz, index_i, index_j, ip_D_a, ip_D_b, &
                     ip_Tmp, ipD(2), ipDDij, ipDij, ipDSij, iShell, iSkal, iT, ix, iy, iz, j, j1, j2, j_R, jBas, jBas_Eff, jBfn, &
                     jCB, jCmp, jlist_s, jShell, jSkal, kDCRE, kDCRR, lDCRER, mAO, mdci, mdcj, mDCRij, mDij, mGrid, nAO, nBfn, nD, &
                     nFunc_i, nFunc_j
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: nGrad_Eff
#endif
real(kind=wp) :: DAij, Factor
integer(kind=iwp), parameter :: Index_d2(3,3) = reshape([5,6,7,6,8,9,7,9,10],[3,3]), &
                                Index_d3(3,3) = reshape([11,14,16,12,17,19,13,18,20],[3,3])
integer(kind=iwp), allocatable :: Ind_Grd(:,:)
integer(kind=iwp), external :: NrOpr

!                                                                      *
!***********************************************************************
!                                                                      *
nD = size(Dens_AO,3)
nAO = size(Dens_AO,1)
Dens_AO(:,:,:) = Zero
mAO = size(TabAO,1)
mGrid = size(TabAO,2)

!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
write(u6,*) 'mAO=',mAO
write(u6,*) 'mGrid=',mGrid
write(u6,*) 'nlist_s=',nlist_s
call RecPrt('Rho: TabAO',' ',TabAO,mAO*mGrid,nAO)
#endif
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
!     Generate the one-particle density matrix, D(mu,nu)               *
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *

if (Do_Grad) then
  call mma_Allocate(Ind_Grd,3,nAO,Label='Ind_Grd')
  Ind_Grd(:,:) = 0
end if

nBfn = size(iBfn_Index,2)
if (nBfn /= nAO) then
  write(u6,*) 'mk_Rho: internal error!'
  call Abend()
end if
Factor = real(2/nD,kind=wp)
do iBfn=1,nBfn
  ilist_s = iBfn_Index(2,iBfn)
  i1 = iBfn_Index(3,iBfn)
  i2 = iBfn_Index(4,iBfn)
  iSkal = list_s(1,ilist_s)
  kDCRE = list_s(2,ilist_s)
  iCmp = iSD(2,iSkal)
  iBas = iSD(3,iSkal)
  mdci = iSD(10,iSkal)
  iShell = iSD(11,iSkal)
  iBas_Eff = list_bas(1,ilist_s)
  index_i = list_bas(2,ilist_s)
  nFunc_i = iBas*iCmp

  i_R = (i1-1)*iBas_Eff+i2
  iCB = Indx(index_i-1+i_R)

  if (Do_Grad) Ind_Grd(:,iBfn) = List_g(:,ilist_s)

  do jBfn=1,iBfn
    jlist_s = iBfn_Index(2,jBfn)
    j1 = iBfn_Index(3,jBfn)
    j2 = iBfn_Index(4,jBfn)
    jSkal = list_s(1,jlist_s)
    kDCRR = list_s(2,jlist_s)
    jCmp = iSD(2,jSkal)
    jBas = iSD(3,jSkal)
    mdcj = iSD(10,jSkal)
    jShell = iSD(11,jSkal)
    jBas_Eff = list_bas(1,jlist_s)
    index_j = list_bas(2,jlist_s)
    nFunc_j = jBas*jCmp

    j_R = (j1-1)*jBas_Eff+j2
    jCB = Indx(index_j-1+j_R)

    ijS = iTri(iShell,jShell)
    ip_Tmp = ipDijs
    call Dens_Info(ijS,ipDij,ipDSij,mDCRij,ipDDij,ip_Tmp,nD)

    iER = ieor(kDCRE,kDCRR)
    lDCRER = NrOpr(iER)

    mDij = nFunc_i*nFunc_j
    ip_D_a = ipDij+lDCRER*mDij
    ip_D_b = ip_D_a
    if (nD /= 1) ip_D_b = ipDSij+lDCRER*mDij
    ipD(1) = ip_D_a
    ipD(2) = ip_D_b

    ij_D = (jCB-1)*nFunc_i+iCB-1
    do iD=1,nD
      DAij = DeDe(ipD(iD)+ij_D)*Fact(mdci,mdcj)*Factor
      Dens_AO(iBfn,jBfn,iD) = DAij
      Dens_AO(jBfn,iBfn,iD) = DAij
    end do

  end do

end do
!                                                                      *
!***********************************************************************
!                                                                      *
!#define _ANALYSIS_
#ifdef _ANALYSIS_
Thr = 1.0e-15_wp
write(u6,*)
write(u6,*) ' Sparsity analysis of D(i,j)'
write(u6,*) ' Threshold: ',Thr
write(u6,*) ' Grid size: ',mGrid
write(u6,*) ' Dimension: ',n,' x ',n
n = size(Dens_AO,1)
n2 = n**2
do iD=1,nD
  m = 0
  do i=1,n
    do j=1,n
      if (abs(Dens_AO(i,j,iD)) < Thr) m = m+1
    end do
  end do
  write(u6,*) 'Total Sparsity in %',100.0_wp*real(m,kind=wp)/real(n2,kind=wp)
  k = 0
  do i=1,n
    m = 0
    do j=1,n
      if (abs(Dens_AO(i,j,iD)) < Thr) m = m+1
    end do
    if (m == n) k = k+1
  end do
  write(u6,*) 'Column Sparsity in %',100.0_wp*real(k,kind=wp)/real(n,kind=wp)
  k = 0
  do j=1,n
    m = 0
    do i=1,n
      if (abs(Dens_AO(i,j,iD)) < Thr) m = m+1
    end do
    if (m == n) k = k+1
  end do
  write(u6,*) 'Row Sparsity in %',100.0_wp*real(k,kind=wp)/real(n,kind=wp)
end do
#endif
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Construct: Sum_i D_ij TabAO(i,iGrid,iAO)
!      D_ij is the one-electron density
!      TabAO(i,iGrid,iAO) are values with respect to the ith AO
!      i=1 is the value of the AO
!      i=2-4 are the values of the first order derivatives
!      i=5-10 are the values of the second order derivatives
!      i=11-20 are the values of the third order derivatives
!
! During a gradient calculation the size of the fast index of
! TabAO is larger than that of Grid_AO. In those cases we copy
! the part of TabAO which we need to TabAO_Short before we make the
! contraction with the 1-particle density matrix.

if (Do_Grad) then
  TabAO_Short(1:kAO,1:mGrid,:) = TabAO(1:kAO,1:mGrid,:)
  call DGEMM_('N','N',kAO*mGrid,nAO*nD,nAO,One,TabAO_Short,kAO*mGrid,Dens_AO,nAO,Zero,Grid_AO,kAO*mGrid)
else
  call DGEMM_('N','N',kAO*mGrid,nAO*nD,nAO,One,TabAO,mAO*mGrid,Dens_AO,nAO,Zero,Grid_AO,kAO*mGrid)
end if
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
if (allocated(dRho_dR)) dRho_dR(:,:,:) = Zero

select case (Functional_Type)
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  case (LDA_Type)
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
    Rho(:,1:mGrid) = Zero
    do iD=1,nD
      do iAO=1,nAO

        do iGrid=1,mGrid
          Rho(iD,iGrid) = Rho(iD,iGrid)+Grid_AO(1,iGrid,iAO,iD)*TabAO(1,iGrid,iAO)
        end do

        if (Do_Grad) then

          ! Loop over cartesian components

          do iCar=1,3

            Ind_xyz = Ind_Grd(iCar,iAO)
            j = iCar+1

            if (Ind_xyz /= 0) then
              do iGrid=1,mGrid

                ! Cartesian derivative of the density.

                dRho_dR(iD,iGrid,Ind_xyz) = dRho_dR(iD,iGrid,Ind_xyz)+Two*Grid_AO(1,iGrid,iAO,iD)*TabAO(j,iGrid,iAO)
              end do
            end if

          end do

        end if
      end do
    end do
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
  case (GGA_Type)
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
    Rho(:,1:mGrid) = Zero
    GradRho(:,1:mGrid) = Zero
    do iD=1,nD
      ix = (iD-1)*3+1
      iy = (iD-1)*3+2
      iz = (iD-1)*3+3
      do iAO=1,nAO

        do iGrid=1,mGrid
          Rho(iD,iGrid) = Rho(iD,iGrid)+Grid_AO(1,iGrid,iAO,iD)*TabAO(1,iGrid,iAO)
          GradRho(ix,iGrid) = GradRho(ix,iGrid)+Grid_AO(1,iGrid,iAO,iD)*TabAO(2,iGrid,iAO)+ &
                              Grid_AO(2,iGrid,iAO,iD)*TabAO(1,iGrid,iAO)
          GradRho(iy,iGrid) = GradRho(iy,iGrid)+Grid_AO(1,iGrid,iAO,iD)*TabAO(3,iGrid,iAO)+ &
                              Grid_AO(3,iGrid,iAO,iD)*TabAO(1,iGrid,iAO)
          GradRho(iz,iGrid) = GradRho(iz,iGrid)+Grid_AO(1,iGrid,iAO,iD)*TabAO(4,iGrid,iAO)+ &
                              Grid_AO(4,iGrid,iAO,iD)*TabAO(1,iGrid,iAO)
        end do

        if (Do_Grad) then

          ! Loop over cartesian components

          do iCar=1,3

            Ind_xyz = Ind_Grd(iCar,iAO) ! index of  nuclear gradient

            j = iCar+1 ! index derivative of AO

            iDx = nD+(iD-1)*3+1 ! index of grad rho component
            iDy = iDx+1
            iDz = iDy+1

            idjx = Index_d2(1,iCar)
            idjy = Index_d2(2,iCar)
            idjz = Index_d2(3,iCar)
            if (Ind_xyz /= 0) then
              do iGrid=1,mGrid

                ! Cartesian derivative of rho

                dRho_dR(iD,iGrid,Ind_xyz) = dRho_dR(iD,iGrid,Ind_xyz)+Two*Grid_AO(1,iGrid,iAO,iD)*TabAO(j,iGrid,iAO)

                ! Cartesian derivatives of grad rho

                dRho_dR(iDx,iGrid,Ind_xyz) = dRho_dR(iDx,iGrid,Ind_xyz)+Two*TabAO(idjx,iGrid,iAO)*Grid_AO(1,iGrid,iAO,iD)+ &
                                             Two*TabAO(j,iGrid,iAO)*Grid_AO(2,iGrid,iAO,iD)
                dRho_dR(iDy,iGrid,Ind_xyz) = dRho_dR(iDy,iGrid,Ind_xyz)+Two*TabAO(idjy,iGrid,iAO)*Grid_AO(1,iGrid,iAO,iD)+ &
                                             Two*TabAO(j,iGrid,iAO)*Grid_AO(3,iGrid,iAO,iD)
                dRho_dR(iDz,iGrid,Ind_xyz) = dRho_dR(iDz,iGrid,Ind_xyz)+Two*TabAO(idjz,iGrid,iAO)*Grid_AO(1,iGrid,iAO,iD)+ &
                                             Two*TabAO(j,iGrid,iAO)*Grid_AO(4,iGrid,iAO,iD)
              end do
            end if

          end do
        end if

      end do
    end do
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
  case (meta_GGA_Type1)
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
    Rho(:,1:mGrid) = Zero
    GradRho(:,1:mGrid) = Zero
    Tau(:,1:mGrid) = Zero
    do iD=1,nD
      ix = (iD-1)*3+1
      iy = (iD-1)*3+2
      iz = (iD-1)*3+3
      do iAO=1,nAO

        do iGrid=1,mGrid
          Rho(iD,iGrid) = Rho(iD,iGrid)+Grid_AO(1,iGrid,iAO,iD)*TabAO(1,iGrid,iAO)
          GradRho(ix,iGrid) = GradRho(ix,iGrid)+Grid_AO(1,iGrid,iAO,iD)*TabAO(2,iGrid,iAO)+ &
                              Grid_AO(2,iGrid,iAO,iD)*TabAO(1,iGrid,iAO)
          GradRho(iy,iGrid) = GradRho(iy,iGrid)+Grid_AO(1,iGrid,iAO,iD)*TabAO(3,iGrid,iAO)+ &
                              Grid_AO(3,iGrid,iAO,iD)*TabAO(1,iGrid,iAO)
          GradRho(iz,iGrid) = GradRho(iz,iGrid)+Grid_AO(1,iGrid,iAO,iD)*TabAO(4,iGrid,iAO)+ &
                              Grid_AO(4,iGrid,iAO,iD)*TabAO(1,iGrid,iAO)
          Tau(iD,iGrid) = Tau(iD,iGrid)+Grid_AO(2,iGrid,iAO,iD)*TabAO(2,iGrid,iAO)+Grid_AO(3,iGrid,iAO,iD)*TabAO(3,iGrid,iAO)+ &
                          Grid_AO(4,iGrid,iAO,iD)*TabAO(4,iGrid,iAO)
        end do

        if (Do_Grad) then

          ! Loop over cartesian components

          do iCar=1,3

            Ind_xyz = Ind_Grd(iCar,iAO) ! index of  nuclear gradient

            j = iCar+1 ! index derivative of AO

            iDx = nD+(iD-1)*3+1 ! index of grad rho component
            iDy = iDx+1
            iDz = iDy+1

            iT = nD*4+iD ! index of tau component

            idjx = Index_d2(1,iCar)
            idjy = Index_d2(2,iCar)
            idjz = Index_d2(3,iCar)
            if (Ind_xyz /= 0) then
              do iGrid=1,mGrid

                ! Cartesian derivative of rho

                dRho_dR(iD,iGrid,Ind_xyz) = dRho_dR(iD,iGrid,Ind_xyz)+Two*Grid_AO(1,iGrid,iAO,iD)*TabAO(j,iGrid,iAO)

                ! Cartesian derivatives of grad rho

                dRho_dR(iDx,iGrid,Ind_xyz) = dRho_dR(iDx,iGrid,Ind_xyz)+Two*TabAO(idjx,iGrid,iAO)*Grid_AO(1,iGrid,iAO,iD)+ &
                                             Two*TabAO(j,iGrid,iAO)*Grid_AO(2,iGrid,iAO,iD)
                dRho_dR(iDy,iGrid,Ind_xyz) = dRho_dR(iDy,iGrid,Ind_xyz)+Two*TabAO(idjy,iGrid,iAO)*Grid_AO(1,iGrid,iAO,iD)+ &
                                             Two*TabAO(j,iGrid,iAO)*Grid_AO(3,iGrid,iAO,iD)
                dRho_dR(iDz,iGrid,Ind_xyz) = dRho_dR(iDz,iGrid,Ind_xyz)+Two*TabAO(idjz,iGrid,iAO)*Grid_AO(1,iGrid,iAO,iD)+ &
                                             Two*TabAO(j,iGrid,iAO)*Grid_AO(4,iGrid,iAO,iD)

                ! Cartesian derivatives of tau

                dRho_dR(iT,iGrid,Ind_xyz) = dRho_dR(iT,iGrid,Ind_xyz)+Four*TabAO(idjx,iGrid,iAO)*Grid_AO(2,iGrid,iAO,iD)+ &
                                            Four*TabAO(idjy,iGrid,iAO)*Grid_AO(3,iGrid,iAO,iD)+ &
                                            Four*TabAO(idjz,iGrid,iAO)*Grid_AO(4,iGrid,iAO,iD)
              end do
            end if

          end do

        end if
      end do
    end do
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
  case (meta_GGA_Type2)
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
    Rho(:,1:mGrid) = Zero
    GradRho(:,1:mGrid) = Zero
    Tau(:,1:mGrid) = Zero
    Lapl(:,1:mGrid) = Zero
    do iD=1,nD
      ix = (iD-1)*3+1
      iy = (iD-1)*3+2
      iz = (iD-1)*3+3
      do iAO=1,nAO

        do iGrid=1,mGrid
          Rho(iD,iGrid) = Rho(iD,iGrid)+Grid_AO(1,iGrid,iAO,iD)*TabAO(1,iGrid,iAO)
          GradRho(ix,iGrid) = GradRho(ix,iGrid)+Grid_AO(1,iGrid,iAO,iD)*TabAO(2,iGrid,iAO)+ &
                              Grid_AO(2,iGrid,iAO,iD)*TabAO(1,iGrid,iAO)
          GradRho(iy,iGrid) = GradRho(iy,iGrid)+Grid_AO(1,iGrid,iAO,iD)*TabAO(3,iGrid,iAO)+ &
                              Grid_AO(3,iGrid,iAO,iD)*TabAO(1,iGrid,iAO)
          GradRho(iz,iGrid) = GradRho(iz,iGrid)+Grid_AO(1,iGrid,iAO,iD)*TabAO(4,iGrid,iAO)+ &
                              Grid_AO(4,iGrid,iAO,iD)*TabAO(1,iGrid,iAO)
          Tau(iD,iGrid) = Tau(iD,iGrid)+Grid_AO(2,iGrid,iAO,iD)*TabAO(2,iGrid,iAO)+Grid_AO(3,iGrid,iAO,iD)*TabAO(3,iGrid,iAO)+ &
                          Grid_AO(4,iGrid,iAO,iD)*TabAO(4,iGrid,iAO)
          Lapl(iD,iGrid) = Lapl(iD,iGrid)+TabAO(1,iGrid,iAO)*(Grid_AO(5,iGrid,iAO,iD)+Grid_AO(8,iGrid,iAO,iD)+ &
                           Grid_AO(10,iGrid,iAO,iD))+Two*(Grid_AO(2,iGrid,iAO,iD)*TabAO(2,iGrid,iAO)+ &
                           Grid_AO(3,iGrid,iAO,iD)*TabAO(3,iGrid,iAO)+Grid_AO(4,iGrid,iAO,iD)*TabAO(4,iGrid,iAO))+ &
                           Grid_AO(1,iGrid,iAO,iD)*(TabAO(5,iGrid,iAO)+TabAO(8,iGrid,iAO)+TabAO(10,iGrid,iAO))
        end do

        if (Do_Grad) then

          ! Loop over cartesian components

          do iCar=1,3

            Ind_xyz = Ind_Grd(iCar,iAO) ! index of  nuclear gradient

            j = iCar+1 ! index derivative of AO

            iDx = nD+(iD-1)*3+1 ! index of grad rho component
            iDy = iDx+1
            iDz = iDy+1

            iT = nD*4+iD ! index of tau component

            iL = nD*5+iD ! index of laplacian component

            idjx = Index_d2(1,iCar)
            idjy = Index_d2(2,iCar)
            idjz = Index_d2(3,iCar)

            idjx2 = Index_d3(1,iCar)
            idjy2 = Index_d3(2,iCar)
            idjz2 = Index_d3(3,iCar)
            idx2 = Index_d2(1,1)
            idy2 = Index_d2(2,2)
            idz2 = Index_d2(3,3)
            if (Ind_xyz /= 0) then
              do iGrid=1,mGrid

                ! Cartesian derivative of rho

                dRho_dR(iD,iGrid,Ind_xyz) = dRho_dR(iD,iGrid,Ind_xyz)+Two*Grid_AO(1,iGrid,iAO,iD)*TabAO(j,iGrid,iAO)

                ! Cartesian derivatives of grad rho

                dRho_dR(iDx,iGrid,Ind_xyz) = dRho_dR(iDx,iGrid,Ind_xyz)+Two*TabAO(idjx,iGrid,iAO)*Grid_AO(1,iGrid,iAO,iD)+ &
                                             Two*TabAO(j,iGrid,iAO)*Grid_AO(2,iGrid,iAO,iD)
                dRho_dR(iDy,iGrid,Ind_xyz) = dRho_dR(iDy,iGrid,Ind_xyz)+Two*TabAO(idjy,iGrid,iAO)*Grid_AO(1,iGrid,iAO,iD)+ &
                                             Two*TabAO(j,iGrid,iAO)*Grid_AO(3,iGrid,iAO,iD)
                dRho_dR(iDz,iGrid,Ind_xyz) = dRho_dR(iDz,iGrid,Ind_xyz)+Two*TabAO(idjz,iGrid,iAO)*Grid_AO(1,iGrid,iAO,iD)+ &
                                             Two*TabAO(j,iGrid,iAO)*Grid_AO(4,iGrid,iAO,iD)

                ! Cartesian derivatives of tau

                dRho_dR(iT,iGrid,Ind_xyz) = dRho_dR(iT,iGrid,Ind_xyz)+Four*TabAO(idjx,iGrid,iAO)*Grid_AO(2,iGrid,iAO,iD)+ &
                                            Four*TabAO(idjy,iGrid,iAO)*Grid_AO(3,iGrid,iAO,iD)+ &
                                            Four*TabAO(idjz,iGrid,iAO)*Grid_AO(4,iGrid,iAO,iD)

                ! Cartesian derivatives of the laplacian

                dRho_dR(iL,iGrid,Ind_xyz) = dRho_dR(iL,iGrid,Ind_xyz)+Two*Grid_AO(1,iGrid,iAO,iD)*(TabAO(idjx2,iGrid,iAO)+ &
                                            TabAO(idjy2,iGrid,iAO)+TabAO(idjz2,iGrid,iAO))+Two*(Grid_AO(idx2,iGrid,iAO,iD)+ &
                                            Grid_AO(idy2,iGrid,iAO,iD)+Grid_AO(idz2,iGrid,iAO,iD))*TabAO(j,iGrid,iAO)+ &
                                            Four*Grid_AO(2,iGrid,iAO,iD)*TabAO(idjx,iGrid,iAO)+ &
                                            Four*Grid_AO(3,iGrid,iAO,iD)*TabAO(idjy,iGrid,iAO)+ &
                                            Four*Grid_AO(4,iGrid,iAO,iD)*TabAO(idjz,iGrid,iAO)

              end do
            end if

          end do

        end if

      end do
    end do
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
  case default
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
    call abend()
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
end select
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
#ifdef _ANALYSIS_
write(u6,*)
write(u6,*) 'Rho Sparsity analysis'
n = 0
do iGrid=1,mGrid
  tmp = Zero
  do iD=1,nD
    tmp = tmp+Rho(iD,iGrid)
  end do
  if (tmp < Thr) n = n+1
end do
write(u6,*) 'Rho Sparsity in %: ',100.0_wp*real(n,kind=wp)/real(mGrid,kind=wp)
#endif
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
!      Scale Tau to compy with the Libxc definition.
!
if (allocated(Tau)) Tau(:,1:mGrid) = Half*Tau(:,1:mGrid)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
do iD=1,nD
  call RecPrt('Dens_AO',' ',Dens_AO(:,:,iD),nAO,nAO)
  call RecPrt('Grid_AO',' ',Grid_AO(:,:,:,iD),mAO*mGrid,nAO)
end do
if (Do_Grad) then
  nGrad_Eff = size(dRho_dR,3)
  call RecPrt('dRho_dR_LDA: dRho_dR',' ',dRho_dR,size(dRho_dR,1)*mGrid,nGrad_Eff)
end if
#endif
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
if (allocated(Ind_grd)) call mma_deAllocate(Ind_Grd)

return

end subroutine Mk_Rho

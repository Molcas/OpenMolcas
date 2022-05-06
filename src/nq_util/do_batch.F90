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
!               2021, Jie Bao                                          *
!***********************************************************************

subroutine Do_Batch(Kernel,Func,mGrid,list_s,nlist_s,List_Exp,List_Bas,Index,nIndex,FckInt,nFckDim,nFckInt,mAO,nD,nP2_ontop,Do_Mo, &
                    TabMO,TabSO,nMOs,Do_Grad,Grad,nGrad,ndRho_dR,nGrad_Eff,iNQ,EG_OT,nTmpPUVX,PDFTPot1,PDFTFocI,PDFTFocA)
!***********************************************************************
!      Author:Roland Lindh, Department of Chemical Physics, University *
!             of Lund, SWEDEN. November 2000                           *
!***********************************************************************

use iSD_data
use SOAO_Info, only: iAOtSO
use Real_Spherical
use Basis_Info
use Center_Info
use Phase_Info
use KSDFT_Info
use nq_Grid, only: Grid, Weights, Rho, nRho
use nq_Grid, only: GradRho, Sigma
use nq_Grid, only: l_CASDFT, TabAO, TabAO_Pack, dRho_dR
use nq_Grid, only: F_xc, F_xca, F_xcb, kAO, Grid_AO
use nq_Grid, only: Fact, Angular, Mem
use nq_Grid, only: D1UnZip, P2UnZip
use nq_Grid, only: Dens_AO, iBfn_Index
use nq_pdft
use nq_MO, only: CMO, D1MO, P2_ontop
use Grid_On_Disk
use nq_Info

implicit real*8(A-H,O-Z)
external Kernel
#include "SysDef.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "debug.fh"
#include "nsd.fh"
#include "setup.fh"
#include "pamint.fh"
integer list_s(2,nlist_s), List_Exp(nlist_s), index(nIndex), List_Bas(2,nlist_s)
real*8 A(3), RA(3), Grad(nGrad), FckInt(nFckInt,nFckDim), TabMO(mAO,mGrid,nMOs), TabSO(mAO,mGrid,nMOs), PDFTPot1(nPot1), &
       PDFTFocI(nPot1), PDFTFocA(nPot1)
logical Do_Grad, Do_Mo
logical l_tanhr
real*8 P2_ontop_d(nP2_ontop,nGrad_Eff,mGrid)
real*8, dimension(:), allocatable :: P2MOCube, P2MOCubex, P2MOCubey, P2MOCubez, MOs, MOx, MOy, MOz
! MOs,MOx,MOy and MOz are for active MOs.
! MOas is for all MOs.
integer nPMO3p
real*8 EG_OT(nTmpPUVX)
real*8, allocatable :: RhoI(:,:), RhoA(:,:)
real*8, allocatable :: TabAO_Tmp(:)
integer :: TabAO_Size(2)
integer, allocatable :: Tmp_Index(:,:)
! Statement function
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

!                                                                      *
!***********************************************************************
!                                                                      *
nCMO = size(CMO)
l_tanhr = .false.

if (l_casdft) then
  call PDFTMemAlloc(mGrid,nOrbt)
  mRho = nP2_ontop
  call mma_allocate(RhoI,mRho,mGrid,Label='RhoI')
  call mma_allocate(RhoA,mRho,mGrid,Label='RhoA')
  RhoI(:,:) = Zero
  RhoA(:,:) = Zero
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Set up an indexation translation between the running index of
! the AOIntegrals and the actual basis function index

nBfn = 0
do iList_s=1,nList_s
  iBas_Eff = List_Bas(1,ilist_s)
  iSkal = list_s(1,ilist_s)
  iCmp = iSD(2,iSkal)
  nBfn = nBfn+iBas_Eff*iCmp
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Evaluate the AOs on the grid points.                                 *
!                                                                      *
!***********************************************************************
!                                                                      *
TabAO(:,:,:) = Zero
TabAO_Size(:) = 0
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute AO's or retrive from disk

if ((NQ_Direct == Off) .and. ((Grid_Status == Use_Old) .and. (.not. Do_Grad) .and. (Functional_Type == Old_Functional_Type))) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Retrieve (and unpack) the AOs from disc

  call iDaFile(Lu_Grid,2,TabAO_Size,2,iDisk_Grid)
  if (TabAO_Size(1) == 0) then
    call Terminate()
    return
  end if
  nBfn = TabAO_Size(1)
  call mma_Allocate(iBfn_Index,6,nBfn,Label='iBfn_Index')
  call iDaFile(Lu_Grid,2,iBfn_Index,size(iBfn_Index),iDisk_Grid)

  nByte = TabAO_Size(2)
  if (Packing == On) then
    mTabAO = (nByte+RtoB-1)/RtoB
  else
    mTabAO = nByte
  end if
  call dDaFile(Lu_Grid,2,TabAO,mTabAO,iDisk_Grid)

  if (Packing == On) then
    nData = size(TabAO)
    call mma_Allocate(TabAO_Tmp,nData,Label='TabAO_Tmp')
    nByte = TabAO_Size(2)
    call UpkR8(0,nData,nByte,TabAO_Pack,TabAO_Tmp)
    TabAO_Pack(:) = TabAO_Tmp(:)
    call mma_deAllocate(TabAO_Tmp)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Compute the AO's

  call mma_Allocate(iBfn_Index,6,nBfn,Label='iBfn_Index')
  iBfn_Index(:,:) = 0

  ! Generate the values of the AOs on the grid

  Mem(:) = Zero
  ipxyz = 1

! define _ANALYSIS_
# ifdef _ANALYSIS_
  Thr = T_Y
  write(6,*)
  write(6,*) ' Sparsity analysis of AO blocks'
  mlist_s = 0
# endif
  !iOff = 1
  iBfn = 0
  iBfn_s = 0
  iBfn_e = 0
  do ilist_s=1,nlist_s
    ish = list_s(1,ilist_s)

    iShll = iSD(0,iSh)
    iAng = iSD(1,iSh)
    iCmp = iSD(2,iSh)
    iBas = iSD(3,iSh)
    iBas_Eff = List_Bas(1,ilist_s)
    iPrim = iSD(5,iSh)
    iPrim_Eff = List_Exp(ilist_s)
    iAO = iSD(7,iSh)
    mdci = iSD(10,iSh)
    iShll = iSD(0,iSh)
    iCnttp = iSD(13,iSh)
    iCnt = iSD(14,iSh)
    A(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)

    ! Set up the unsifted version of iBfn_Index

    iAdd = iBas-iBas_Eff
    iBfn_s = iBfn+1
    do i1=1,iCmp
      iSO1 = iAOtSO(iAO+i1,0) ! just used when nIrrep=1
      do i2=1,iBas_Eff
        IndAO1 = i2+iAdd
        Indi = iSO1+IndAO1-1

        iBfn = iBfn+1
        iBfn_Index(1,iBfn) = Indi
        iBfn_Index(2,iBfn) = ilist_s
        iBfn_Index(3,iBfn) = i1
        iBfn_Index(4,iBfn) = i2
        iBfn_Index(5,iBfn) = mdci
        iBfn_Index(6,iBfn) = IndAO1
      end do
    end do
    iBfn_e = iBfn

    nDrv = mRad-1
    nForm = 0
    do iDrv=0,nDrv
      nForm = nForm+nElem(iDrv)
    end do
    nTerm = 2**nDrv
    nxyz = mGrid*3*(iAng+mRad)
    !nRadial = iBas_Eff*mGrid*mRad
    ipRadial = ipxyz+nxyz

    iR = list_s(2,ilist_s)

    ipx = iPhase(1,iR)
    ipy = iPhase(2,iR)
    ipz = iPhase(3,iR)
    px = dble(iPhase(1,iR))
    py = dble(iPhase(2,iR))
    pz = dble(iPhase(3,iR))
    RA(1) = px*A(1)
    RA(2) = py*A(2)
    RA(3) = pz*A(3)

    ! Evaluate AOs at RA

    call AOEval(iAng,mGrid,Grid,Mem(ipxyz),RA,Shells(iShll)%Transf,RSph(ipSph(iAng)),nElem(iAng),iCmp,Angular,nTerm,nForm,T_Y, &
                mRad,iPrim,iPrim_Eff,Shells(iShll)%Exp,Mem(ipRadial),iBas_Eff,Shells(iShll)%pCff(1,iBas-iBas_Eff+1), &
                TabAO(:,:,iBfn_s:),mAO,px,py,pz,ipx,ipy,ipz)
#   ifdef _ANALYSIS_
    ix = iDAMax_(mAO*mGrid*iBas_Eff*iCmp,TabAO_Pack(iOff),1)
    TMax = abs(TabAO_Pack(iOff-1+ix))
    if (TMax < Thr) then
      mlist_s = mlist_s+1
      write(6,*) ' ilist_s: ',ilist_s
      write(6,*) ' TMax:    ',TMax
    end if
#   endif

    ! At this time eliminate individual basis functions which have
    ! an insignificant contribution to any of the grid points we
    ! are processing at this stage.

    Thr = T_Y
    iSkip = 0
    kBfn = iBfn_s-1
    do jBfn=iBfn_s,iBfn_e
      call Spectre(SMax)
      if (SMax < Thr) then
        iSkip = iSkip+1
      else
        kBfn = kBfn+1
        if (kBfn /= jBfn) then
          TabAO(:,:,kBfn) = TabAO(:,:,jBfn)
          iBfn_Index(:,kBfn) = iBfn_Index(:,jBfn)
        end if
      end if
    end do
    iBfn = kBfn

    !iOff = iBfn*mAO*mGrid+1

    !                                                                  *
    !*******************************************************************
    !                                                                  *
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! reduced the size of the table to be exactly that of the
  ! number of functions that have non-zero contributions.
  if (iBfn /= nBfn) then
    if (iBfn == 0) then
      TabAO_Size(:) = 0
      call iDaFile(Lu_Grid,1,TabAO_Size,2,iDisk_Grid)
      call Terminate()
      return
    end if
    call mma_allocate(Tmp_Index,6,iBfn,Label='Tmp_Index')
    Tmp_Index(:,1:iBfn) = iBfn_Index(:,1:iBfn)
    call mma_deallocate(iBfn_Index)
    nBfn = iBfn
    call mma_Allocate(iBfn_Index,6,nBfn,Label='iBfn_Index')
    iBfn_Index(:,:) = Tmp_Index(:,:)
    call mma_deallocate(Tmp_Index)
  end if
  TabAO_Size(1) = nBfn

# ifdef _ANALYSIS_
  write(6,*) ' % AO blocks that can be eliminated: ',1.0d2*dble(mlist_s)/dble(nlist_s)
  write(6,*)
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_Allocate(Dens_AO,nBfn,nBfn,nD,Label='Dens_AO')
call mma_Allocate(Grid_AO,kAO,mGrid,nBfn,nD,Label='Grid_AO')
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Evaluate some MOs on the grid                                        *
!                                                                      *
!***********************************************************************
!                                                                      *
if (Do_MO) then

  ! First, symmatry adapt the AOs
  TabSO(:,:,:) = Zero
  jlist_s = 0
  call mk_SOs(TabSO,mAO,mGrid,nMOs,List_s,List_Bas,nList_s,jlist_s)

  ! Second, transform SOs to MOs
  TabMO(:,:,:) = Zero
  call mk_MOs(TabSO,mAO,mGrid,TabMO,nMOs,CMO,nCMO)
end if
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Compute Rho, Grad Rho, Tau, Laplacian, and the Sigma vectors.
! In case of gradient calculations compute Cartesian derivatives
! of Rho, Grad Rho, Tau, and the Laplacian.

call Mk_Rho(list_s,nlist_s,Fact,ndc,list_bas,Index,nIndex,Do_Grad)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
if (l_casdft) then
  Dens_t1 = Dens_t1+Comp_d(Weights,mGrid,Rho,nRho,nD,0)
  Dens_a1 = Dens_a1+Comp_d(Weights,mGrid,Rho,nRho,nD,1)
  Dens_b1 = Dens_b1+Comp_d(Weights,mGrid,Rho,nRho,nD,2)

  nPMO3p = 1
  if (lft .and. lGGA) nPMO3p = mGrid*NASHT

  call mma_allocate(P2MOCube,mGrid*NASHT)
  call mma_allocate(P2MOCubex,nPMO3p)
  call mma_allocate(P2MOCubey,nPMO3p)
  call mma_allocate(P2MOCubez,nPMO3p)
  call mma_allocate(MOs,mGrid*NASHT)
  call mma_allocate(MOx,mGrid*NASHT)
  call mma_allocate(MOy,mGrid*NASHT)
  call mma_allocate(MOz,mGrid*NASHT)

  call CalcP2MOCube(P2MOCube,P2MOCubex,P2MOCubey,P2MOCubez,nPMO3p,MOs,MOx,MOy,MOz,TabMO,P2Unzip,mAO,mGrid,nMOs,Do_Grad)
  call Fzero(P2_ontop,nP2_ontop*mGrid)

  if (.not. Do_Grad) then !regular MO-based run
    call Do_PI2(D1MO,size(D1MO),TabMO,mAO,mGrid,nMOs,P2_ontop,nP2_ontop,RhoI,RhoA,mRho,Do_Grad,P2MOCube,MOs,MOx,MOy,MOz)
  else !AO-based run for gradients
    !nP2_ontop_d = nP2_ontop*mGrid*nGrad_Eff
    P2_ontop_d(:,:,:) = 0
    call Do_Pi2grad(mAO,mGrid,P2_ontop,nP2_ontop,nGrad_Eff,list_s,nlist_s,list_bas,D1MO,size(D1MO),TabMO,P2_ontop_d,RhoI,RhoA, &
                    mRho,nMOs,CMO,nCMO,TabSO,lft,P2MOCube,P2MOCubex,P2MOCubey,P2MOCubez,nPMO3p,MOs,MOx,MOy,MOz)
  end if

  call TranslateDens(P2_OnTop,dRho_dr,P2_OnTop_d,l_tanhr,nRho,mGrid,nP2_OnTop,ndRho_dR,nGrad_Eff,Do_Grad)

  call mma_deallocate(P2MOCube)
  call mma_deallocate(P2MOCubex)
  call mma_deallocate(P2MOCubey)
  call mma_deallocate(P2MOCubez)
  call mma_deallocate(MOs)
  call mma_deallocate(MOx)
  call mma_deallocate(MOy)
  call mma_deallocate(MOz)

  ! Integrate out the number of electrons
  Dens_t2 = Dens_t2+Comp_d(Weights,mGrid,Rho,nRho,nD,0)
  Dens_a2 = Dens_a2+Comp_d(Weights,mGrid,Rho,nRho,nD,1)
  Dens_b2 = Dens_b2+Comp_d(Weights,mGrid,Rho,nRho,nD,2)

end if
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
if (allocated(Sigma)) then
  if (size(Sigma,1) == 1) then
    do iGrid=1,mGrid
      Sigma(1,iGrid) = GradRho(1,iGrid)**2+GradRho(2,iGrid)**2+GradRho(3,iGrid)**2
    end do
  else
    do iGrid=1,mGrid
      Sigma(1,iGrid) = GradRho(1,iGrid)**2+GradRho(2,iGrid)**2+GradRho(3,iGrid)**2
      Sigma(2,iGrid) = GradRho(1,iGrid)*GradRho(4,iGrid)+GradRho(2,iGrid)*GradRho(5,iGrid)+GradRho(3,iGrid)*GradRho(6,iGrid)
      Sigma(3,iGrid) = GradRho(4,iGrid)**2+GradRho(5,iGrid)**2+GradRho(6,iGrid)**2
    end do
  end if
end if
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Integrate out the number of electrons, |grad|, and tau

Dens_I = Dens_I+Compute_Rho(Weights,mGrid,nD)
select case (Functional_type)
  case (LDA_Type)
  case (GGA_type)
    Grad_I = Grad_I+Compute_Grad(Weights,mGrid,nD)
  case (meta_GGA_type1,meta_GGA_type2)
    Grad_I = Grad_I+Compute_Grad(Weights,mGrid,nD)
    Tau_I = Tau_I+Compute_Tau(Weights,mGrid,nD)
end select
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Evaluate the functional on the grid                                  *
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! evaluate the energy density, the derivative of the functional with
! respect to rho and grad rho.

call Kernel(mGrid,nD)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Integrate the energy of the functional

Func = Func+DDot_(mGrid,Weights,1,F_xc,1)
if (l_casdft) then
  Funcaa = Funcaa+DDot_(mGrid,Weights,1,F_xca,1)
  Funcbb = Funcbb+DDot_(mGrid,Weights,1,F_xcb,1)
  Funccc = Func-Funcaa-Funcbb
end if
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
if (Do_Grad) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Compute the DFT contribution to the gradient                       *
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call DFT_Grad(Grad,nGrad,nD,Grid,mGrid,dRho_dR,ndRho_dR,nGrad_Eff,Weights,iNQ)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (l_casdft) then
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! For MC-PDFT optionally compute stuff for the CP-MC-PDFT          *
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    if (do_pdftPot) then
      call mma_allocate(MOs,mGrid*NASHT)
      call TransferMO(MOas,TabMO,mAO,mGrid,nMOs,1)
      if (lft .and. lGGA) then
        call TransferMO(MOax,TabMO,mAO,mGrid,nMOs,2)
        call TransferMO(MOay,TabMO,mAO,mGrid,nMOs,3)
        call TransferMO(MOaz,TabMO,mAO,mGrid,nMOs,4)
      end if
      call TransActMO(MOs,TabMO,mAO,mGrid,nMOs)
      call Calc_Pot1(PDFTPot1,TabMO,mAO,mGrid,nMOs,P2_ontop,nP2_ontop,MOas)
      call Calc_Pot2(EG_OT,mGrid,P2_ontop,nP2_ontop)
      call PDFTFock(PDFTFocI,PDFTFocA,D1Unzip,mGrid,MOs)
      call mma_deallocate(MOs)
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  else
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Compute the DFT contribution to the Fock matrix                  *
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    call DFT_Int(list_s,nlist_s,FckInt,nFckInt,nD,Fact,ndc)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! AOs on the grid are (packed and) written to disk.

if ((NQ_Direct == Off) .and. ((Grid_Status == Regenerate) .and. (.not. Do_Grad))) then

  TabAO_Size(1) = nBfn
  if (Packing == On) then

    ! Pack before they are put on disc

    nData = mAO*mGrid*nBfn
    call mma_Allocate(TabAO_Tmp,nData,Label='TabAO_Tmp')
    TabAO_Tmp(1:nData) = TabAO_Pack(1:nData)
    call PkR8(0,nData,nByte,TabAO_Tmp,TabAO_Pack)
    mData = (nByte+RtoB-1)/RtoB
    if (mData > nData) then
      call WarningMessage(2,'mData > nData')
      write(6,*) 'nData=',nData
      write(6,*) 'nData=',nData
      call Abend()
    end if
    TabAO_Size(2) = nByte
    call mma_deAllocate(TabAO_Tmp)
  else
    mData = mAO*mGrid*nBfn
    TabAO_Size(2) = mData
  end if

  call iDaFile(Lu_Grid,1,TabAO_Size,2,iDisk_Grid)
  call iDaFile(Lu_Grid,1,iBfn_Index,size(iBfn_Index),iDisk_Grid)
  mTabAO = mData
  call dDaFile(Lu_Grid,1,TabAO,mTabAO,iDisk_Grid)

end if
!                                                                      *
!***********************************************************************
!                                                                      *
call Terminate()

return

contains

subroutine Terminate()

  if (l_casdft) call PDFTMemDeAlloc()

  if (allocated(RhoI)) then
    call mma_deallocate(RhoI)
    call mma_deallocate(RhoA)
  end if
  if (allocated(iBfn_Index)) call mma_deAllocate(iBfn_Index)
  if (allocated(Grid_AO)) call mma_deAllocate(Grid_AO)
  if (allocated(Dens_AO)) call mma_deAllocate(Dens_AO)

end subroutine Terminate

subroutine Spectre(SMax)

  integer iGrid, iAO
  real*8 SMax

  SMax = Zero
  do iGrid=1,mGrid
    do iAO=1,mAO
      SMax = max(SMax,abs(Weights(iGrid)*TabAO(iAO,iGrid,jBfn)))
    end do
  end do
  !if (SMax < Thr) Write(6,*) SMax,TMax

end subroutine Spectre

end subroutine Do_Batch

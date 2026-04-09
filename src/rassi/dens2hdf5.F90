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
! Copyright (C) 2019, Ignacio Fdez. Galvan                             *
!***********************************************************************

module Dens2HDF5

use Symmetry_Info, only: nIrrep, MUL
use rassi_data, only: NBASF
use Cntrl, only: NSTATE
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), allocatable :: IdxState(:,:)
real(kind=wp), parameter :: Thrs = 1.0e-14_wp

public :: UpdateIdx
#ifdef _HDF5_
public :: StoreDens
#endif

contains

!***********************************************************************
!  UpdateIdx
!
!> @brief
!>   Update index of TDMs to save to HDF5
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Update the table of indices of states for the TDMs that will be saved
!> on the HDF5 file. If SUBSets keyword is not used, all TDMs are saved,
!> otherwise only those for the selected transitions. Since the TDMs are
!> saved in SF basis, for a SO calculation we need to figure out which
!> SF TDMs contribute to the desired transitions.
!>
!> @param[in] IndexE  SF states sorted by energy
!> @param[in] nSS     number of SO states
!> @param[in] USOR    SO coefficients in SF basis (real part)
!> @param[in] USOI    SO coefficients in SF basis (imaginary part)
!> @param[in] MapSt   map of SF states expanded by multiplicity
!***********************************************************************
subroutine UpdateIdx(IndexE,nSS,USOR,USOI,MapSt)

  use Cntrl, only: LOOPDIVIDE, LOOPMAX, REDUCELOOP

  integer(kind=iwp), intent(in) :: IndexE(nState), nSS
  real(kind=wp), intent(in), optional :: USOR(nSS,nSS), USOI(nSS,nSS)
  integer(kind=iwp), intent(in), optional :: MapSt(nSS)
  integer(kind=iwp) :: i, i_, iEnd, iSS, iState, j, j_, jEnd, jSS, jStart, jState
  real(kind=wp) :: f1, f2

  if (.not. allocated(IdxState)) then
    call mma_Allocate(IdxState,nState,nState,Label='IdxState')
    IdxState(:,:) = 0
  end if
  jEnd = nState
  if (nSS > 0) jEnd = nSS
  if (ReduceLoop) then
    iEnd = LoopDivide
    jStart = LoopDivide+1
    if (LoopMax > 0) jEnd = min(jEnd,LoopDivide+LoopMax)
  else
    iEnd = nState
    if (nSS > 0) iEnd = nSS
    jStart = 1
  end if
  ! list states for/between which we want to store the density matrices
  if (nSS > 0) then
    do iSS=1,nSS
      do jSS=1,nSS
        if (jSS /= iSS) then
          if (jSS < jStart) cycle
          if (iSS > iEnd) cycle
          if (jSS > jEnd) cycle
        end if
        do i=1,nSS
          iState = MapSt(i)
          do j=1,i
            jState = MapSt(j)
            i_ = max(jState,iState)
            j_ = min(jState,iState)
            if (IdxState(i_,j_) > 0) cycle
            f1 = (USOR(i,iSS)*USOR(j,jSS)+USOI(i,iSS)*USOI(j,jSS))**2+(USOR(i,iSS)*USOI(j,jSS)-USOI(i,iSS)*USOR(j,jSS))**2
            f2 = (USOR(j,iSS)*USOR(i,jSS)+USOI(j,iSS)*USOI(i,jSS))**2+(USOR(j,iSS)*USOI(i,jSS)-USOI(j,iSS)*USOR(i,jSS))**2
            ! this should be Thrs**2, but let's be looser with SO states
            if (max(f1,f2) >= Thrs) IdxState(i_,j_) = 2
          end do
        end do
      end do
    end do
  else
    do i=1,nState
      iState = IndexE(i)
      IdxState(iState,iState) = 1
      if (i > iEnd) cycle
      do j=max(i+1,jStart),nState
        if (j > jEnd) cycle
        jState = IndexE(j)
        i_ = max(jState,iState)
        j_ = min(jState,iState)
        IdxState(i_,j_) = 1
      end do
    end do
  end if

end subroutine UpdateIdx

!***********************************************************************
!  StoreDens
!
!> @brief
!>   Store the density matrices to HDF5
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Save the selected state and transition density matrices to the HDF5
!> file. The TDMs were computed in input state basis, but will be stored
!> in SF eigen states, so a transformation is needed.
!>
!> @param[in] EigVec  coefficients of the SF eigen states in the input state basis
!***********************************************************************
#ifdef _HDF5_
subroutine StoreDens(EigVec)

  use mh5, only: mh5_put_dset
  use rassi_aux, only: iDisk_TDM
  use rassi_global_arrays, only: JbNum
  use rassi_data, only: NTDMZZ
  use Cntrl, only: IFSO, IRREP, LuTDM
  use RASSIWfn, only: wfn_SFS_TDM, wfn_SFS_TSDM, wfn_SFS_WETDM
  use Constants, only: Zero

  real(kind=wp), intent(in) :: EigVec(nState,nState)
  integer(kind=iwp) :: iDisk, iEmpty, iGo, iOpt, iState, iSy12, iSym1, iSym2, Job1, Job2, jState, k, l, nThisTDMZZ
  logical(kind=iwp) :: isZero(3)
  real(kind=wp) :: f1, f2
  real(kind=wp), allocatable :: TDMIJ(:), TDMZZ(:), TSDMIJ(:), TSDMZZ(:), WDMIJ(:), WDMZZ(:)

  if (.not. allocated(IdxState)) return
  ! Transform TDMs to SF eigenstates
  call mma_Allocate(TDMZZ,nTDMZZ,Label='TDMZZ')
  call mma_Allocate(TSDMZZ,nTDMZZ,Label='TSDMZZ')
  call mma_Allocate(WDMZZ,nTDMZZ,Label='WDMZZ')
  call mma_Allocate(TDMIJ,nTDMZZ,Label='TDMIJ')
  call mma_Allocate(TSDMIJ,nTDMZZ,Label='TSDMIJ')
  call mma_Allocate(WDMIJ,nTDMZZ,Label='WDMIJ')
  do iState=1,nState
    do jState=1,iState
      if (IdxState(iState,jState) == 0) cycle
      TDMIJ(:) = Zero
      TSDMIJ(:) = Zero
      WDMIJ(:) = Zero
      isZero = [.true.,.true.,.true.]
      do k=1,nState
        Job1 = JbNum(k)
        iSym1 = Irrep(Job1)
        do l=1,k
          f1 = EigVec(k,iState)*EigVec(l,jState)
          f2 = EigVec(l,iState)*EigVec(k,jState)
          if (max(abs(f1),abs(f2)) < Thrs) cycle
          Job2 = JbNum(l)
          iSym2 = Irrep(Job2)
          iSy12 = Mul(iSym1,iSym2)
          iDisk = iDisk_TDM(k,l,1)
          iEmpty = iDisk_TDM(k,l,2)
          iOpt = 2
          iGo = 3
          if (IfSO) iGo = iGo+4
          call dens2file(TDMZZ,TSDMZZ,WDMZZ,nTDMZZ,LuTDM,iDisk,iEmpty,iOpt,iGo,k,l)
          if (btest(iEmpty,0)) then
            isZero(1) = .false.
            if (abs(f1) >= Thrs) TDMIJ(:) = TDMIJ(:)+f1*TDMZZ(:)
            if ((k /= l) .and. (abs(f2) >= Thrs)) then
              call Transpose_TDM(TDMZZ,iSy12)
              TDMIJ(:) = TDMIJ(:)+f2*TDMZZ(:)
            end if
          end if
          if (btest(iEmpty,1)) then
            isZero(2) = .false.
            if (abs(f1) >= Thrs) TSDMIJ(:) = TSDMIJ(:)+f1*TSDMZZ(:)
            if ((k /= l) .and. (abs(f2) >= Thrs)) then
              call Transpose_TDM(TSDMZZ,iSy12)
              TSDMIJ(:) = TSDMIJ(:)+f2*TSDMZZ(:)
            end if
          end if
          if (IFSO .and. btest(iEmpty,2)) then
            isZero(3) = .false.
            if (abs(f1) >= Thrs) WDMIJ(:) = WDMIJ(:)+f1*WDMZZ(:)
            if ((k /= l) .and. (abs(f2) >= Thrs)) then
              call Transpose_TDM(WDMZZ,iSy12)
              WDMIJ(:) = WDMIJ(:)+f2*WDMZZ(:)
            end if
          end if
        end do
      end do
      if (all(isZero)) cycle
      nThisTDMZZ = 0
      do iSym1=1,nIrrep
        iSym2 = Mul(iSy12,iSym1)
        nThisTDMZZ = nThisTDMZZ+NBASF(iSym1)*NBASF(iSym2)
      end do
      if (.not. isZero(1)) call mh5_put_dset(wfn_sfs_tdm,TDMIJ,[nThisTDMZZ,1,1],[0,iState-1,jState-1])
      if (.not. isZero(2)) call mh5_put_dset(wfn_sfs_tsdm,TSDMIJ,[nThisTDMZZ,1,1],[0,iState-1,jState-1])
      if (IFSO .and. (.not. isZero(3))) call mh5_put_dset(wfn_sfs_wetdm,WDMIJ,[nThisTDMZZ,1,1],[0,iState-1,jState-1])
    end do
  end do
  call mma_deAllocate(IdxState)
  call mma_deAllocate(TDMZZ)
  call mma_deAllocate(TSDMZZ)
  call mma_deAllocate(WDMZZ)
  call mma_deAllocate(TDMIJ)
  call mma_deAllocate(TSDMIJ)
  call mma_deAllocate(WDMIJ)

end subroutine StoreDens

!***********************************************************************
!  Transpose_TDM
!
!> @brief
!>   Transpose a transition density matrix in place.
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Transpose a transition density matrix, stored in symmetry blocks,
!> replacing the original matrix. The matrices contain only the symmetry
!> blocks that match the total symmetry of the transition.
!>
!> @param[in,out] TDM       Transition density matrix
!> @param[in]     Symmetry  Symmetry of the transition
!***********************************************************************
subroutine Transpose_TDM(TDM,Symmetry)

  real(kind=wp), intent(inout) :: TDM(*)
  integer(kind=iwp), intent(in) :: Symmetry
  integer(kind=iwp) :: i, iBlock(0:8), iSym1, iSym2, j, nTot
  real(kind=wp), allocatable :: Tmp(:)

  ! Compute the location of all the stored symmetry blocks
  nTot = 0
  iBlock(0) = 0
  do iSym1=1,nIrrep
    iSym2 = Mul(Symmetry,iSym1)
    nTot = nTot+nBasF(iSym1)*nBasF(iSym2)
    iBlock(iSym1) = nTot
  end do
  ! Make a copy so we can transpose in place
  call mma_Allocate(Tmp,nTot,Label='Tmp')
  Tmp(:) = TDM(1:nTot)
  ! Transpose symmetry block (a,b) onto symmetry block (b,a)
  do iSym1=1,nIrrep
    iSym2 = Mul(Symmetry,iSym1)
    do i=1,nBasF(iSym2)
      do j=1,nBasF(iSym1)
        TDM(iBlock(iSym2-1)+(j-1)*nBasF(iSym2)+i) = Tmp(iBlock(iSym1-1)+(i-1)*nBasF(iSym1)+j)
      end do
    end do
  end do
  call mma_deAllocate(Tmp)

end subroutine Transpose_TDM
#endif

end module Dens2HDF5

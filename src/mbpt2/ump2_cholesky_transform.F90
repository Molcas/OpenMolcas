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
! Copyright (C) 2026, Okan Koeksal                                     *
!***********************************************************************

module UMP2_Cholesky_Transform

use Definitions, only: wp, iwp
use MBPT2_Global, only: MBPT2_All_Finite
use Constants, only: Zero, One
use stdalloc, only: mma_allocate, mma_deallocate
use linalg_mod, only: mult
use UMP2_Global, only: nSym, nBas, nOrb, nOccA, nOccB, &
                       nFro, nOA, nOB, nVA, nVB, CAlpha, CBeta

implicit none
private

public :: UMP2_Transform_AO_Batch

contains

subroutine UMP2_Transform_AO_Batch(LAO,La,Lb,ierr,Message)

  ! Transform a common batch of fully unpacked AO Cholesky vectors.
  !
  ! Input:
  !   LAO(mu,nu,P)
  !
  ! Outputs:
  !   La(a,i,P): alpha virtual, alpha occupied, vector
  !   Lb(b,j,P): beta  virtual, beta  occupied, vector
  !
  ! Required integral convention:
  !
  !   (mu nu | rho sigma)
  !       = sum_P LAO(mu,nu,P)*LAO(rho,sigma,P)
  !
  ! LAO must contain actual symmetric AO matrix elements.
  ! Reduced-pair storage and any normalization factors must already
  ! have been resolved by the OpenMolcas vector reader.
  !
  ! The same vector batch and ordering are used for both spins.
  ! Frozen occupied orbitals are excluded from the transformation.
  !
  ! This transformation is intentionally rank-local. In MPI-v2 each
  ! rank transforms only its locally owned AO Cholesky vectors into its
  ! assigned global P slice. The driver subsequently uses GADGOp to
  ! assemble the complete transformed vector set on every rank before
  ! UMP2_Energy distributes excitation pairs across ranks.
  !
  ! The caller allocates La and Lb before calling this routine.
  ! A module interface is required because the arrays are assumed-shape.

  real(kind=wp), intent(in) :: LAO(:,:,:)
  real(kind=wp), intent(out) :: La(:,:,:), Lb(:,:,:)

  integer(kind=iwp), intent(out) :: ierr
  character(len=*), intent(out) :: Message

  integer(kind=iwp) :: nVec, iVec, nWork
  real(kind=wp) :: Scale, Asymmetry
  real(kind=wp), parameter :: SymTol = 1.0e-12_wp
  real(kind=wp), allocatable :: Work(:,:)

  ierr = 0
  Message = ''

  La(:,:,:) = Zero
  Lb(:,:,:) = Zero

  if (nSym /= 1) then
    call Fail('UMP2 transformation requires C1 symmetry.')
    return
  end if

  if ((nBas <= 0) .or. (nOrb <= 0) .or. (nOrb > nBas)) then
    call Fail('Invalid basis or orbital dimensions.')
    return
  end if

  if ((nOccA < 0) .or. (nOccA > nOrb) .or. &
      (nOccB < 0) .or. (nOccB > nOrb)) then
    call Fail('Invalid spin occupation counts.')
    return
  end if

  if ((nFro < 0) .or. (nFro > min(nOccA,nOccB))) then
    call Fail('Invalid frozen-core count.')
    return
  end if

  if ((nOA /= nOccA-nFro) .or. (nOB /= nOccB-nFro) .or. &
      (nVA /= nOrb-nOccA) .or. (nVB /= nOrb-nOccB)) then
    call Fail('UMP2 orbital spaces have not been set up consistently.')
    return
  end if

  ! Allocation must be checked separately from array dimensions.

  if (.not. allocated(CAlpha)) then
    call Fail('Alpha MO coefficients are missing.')
    return
  end if

  if (.not. allocated(CBeta)) then
    call Fail('Beta MO coefficients are missing.')
    return
  end if

  if ((size(CAlpha,1) /= nBas) .or. &
      (size(CAlpha,2) /= nOrb)) then
    call Fail('Incorrect alpha coefficient dimensions.')
    return
  end if

  if ((size(CBeta,1) /= nBas) .or. &
      (size(CBeta,2) /= nOrb)) then
    call Fail('Incorrect beta coefficient dimensions.')
    return
  end if

  if ((size(LAO,1) /= nBas) .or. &
      (size(LAO,2) /= nBas)) then
    call Fail('AO Cholesky matrices must have dimensions nBas by nBas.')
    return
  end if

  nVec = size(LAO,3)

  if ((size(La,1) /= nVA) .or. &
      (size(La,2) /= nOA) .or. &
      (size(La,3) /= nVec)) then
    call Fail('Incorrect alpha transformed-vector dimensions.')
    return
  end if

  if ((size(Lb,1) /= nVB) .or. &
      (size(Lb,2) /= nOB) .or. &
      (size(Lb,3) /= nVec)) then
    call Fail('Incorrect beta transformed-vector dimensions.')
    return
  end if

  if (.not. MBPT2_All_Finite(CAlpha)) then
    call Fail('Nonfinite alpha MO coefficients.')
    return
  end if

  if (.not. MBPT2_All_Finite(CBeta)) then
    call Fail('Nonfinite beta MO coefficients.')
    return
  end if

  if (.not. MBPT2_All_Finite(LAO)) then
    call Fail('Nonfinite AO Cholesky vector elements.')
    return
  end if

  ! Check the full-matrix contract before transforming anything.
  ! Do not silently symmetrize incorrectly unpacked vectors.

  do iVec=1,nVec

    Scale = max(One,maxval(abs(LAO(:,:,iVec))))
    Asymmetry = maxval(abs(LAO(:,:,iVec)- &
                          transpose(LAO(:,:,iVec))))

    if (Asymmetry > SymTol*Scale) then
      call Fail('AO Cholesky vector is not a symmetric full matrix.')
      return
    end if

  end do

  ! An empty batch has no transformation work.
  ! This does not certify that the complete calculation has no vectors.

  if (nVec == 0) return

  nWork = max(1,nOA,nOB)
  call mma_allocate(Work,nBas,nWork,label='UMP2 AO MO work')

  do iVec=1,nVec

    ! Alpha:
    !
    ! Work(mu,i) = sum_nu LAO(mu,nu,P)*CAlpha(nu,i)
    !
    ! Dimensions:
    !
    !   LAO                 : nBas x nBas
    !   CAlpha(occupied)    : nBas x nOA
    !   Work                : nBas x nOA
    !
    ! La(a,i,P) = sum_mu CAlpha(mu,a)*Work(mu,i)
    !
    ! Dimensions:
    !
    !   transpose(CAlpha virtual) : nVA x nBas
    !   Work                      : nBas x nOA
    !   La                        : nVA x nOA

    if ((nOA > 0) .and. (nVA > 0)) then

      call mult(LAO(:,:,iVec), &
                CAlpha(:,nFro+1:nOccA), &
                Work(:,1:nOA))

      call mult(CAlpha(:,nOccA+1:nOrb), &
                Work(:,1:nOA), &
                La(:,:,iVec), &
                transpA=.true.)

    end if

    ! Beta:
    !
    ! The occupied/virtual boundary is independent of alpha.
    !
    ! Work(mu,j) = sum_nu LAO(mu,nu,P)*CBeta(nu,j)
    !
    ! Lb(b,j,P) = sum_mu CBeta(mu,b)*Work(mu,j)

    if ((nOB > 0) .and. (nVB > 0)) then

      call mult(LAO(:,:,iVec), &
                CBeta(:,nFro+1:nOccB), &
                Work(:,1:nOB))

      call mult(CBeta(:,nOccB+1:nOrb), &
                Work(:,1:nOB), &
                Lb(:,:,iVec), &
                transpA=.true.)

    end if

  end do

  call mma_deallocate(Work)

  if ((.not. MBPT2_All_Finite(La)) .or. &
      (.not. MBPT2_All_Finite(Lb))) then
    call Fail('Nonfinite transformed Cholesky vector elements.')
    return
  end if

contains

  subroutine Fail(Text)

    character(len=*), intent(in) :: Text

    ierr = 1
    Message = Text

    ! Never leave partially transformed output on a detected error.

    La(:,:,:) = Zero
    Lb(:,:,:) = Zero

  end subroutine Fail

end subroutine UMP2_Transform_AO_Batch

end module UMP2_Cholesky_Transform
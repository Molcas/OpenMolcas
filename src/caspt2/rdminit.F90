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
! Copyright (C) 2019, Stefano Battaglia                                *
!***********************************************************************

subroutine rdminit()

use PrintLevel, only: DEBUG
#ifdef _DMRG_
use iso_c_binding, only: c_int
use qcmaquis_interface, only: qcmaquis_interface_set_state
use caspt2_module, only: DMRG, mState
#endif
use caspt2_global, only: CMO, CMO_Internal, DMIX, DREF, DWGT, iPrGlb, LUONEM, NCMO
use caspt2_module, only: iAd1m, ISCF, nConf, nState
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: I, J, iDisk
real(kind=wp) :: Wij
real(kind=wp), allocatable :: CI(:)

if (IPRGLB >= DEBUG) write(u6,*) ' Entered rdminit.'

! Get CASSCF MO coefficients
call mma_allocate(CMO_Internal,NCMO,Label='CMO_Internal')
CMO => CMO_Internal
IDISK = IAD1M(1)
call ddafile(LUONEM,2,CMO,NCMO,IDISK)

! Allocate memory for CI vector
call mma_allocate(CI,Nconf,Label='CI')

! Initialize array of 1-RDMs with zeros
DMIX(:,:) = Zero

! Start long loop over all states and compute the weighted density
! of each state using the weights in DWGT
do I=1,Nstate

  if (ISCF /= 0) then
    ! Then we still need the "CI array": It is used in subroutine calls
    CI(1) = One
  else
    ! Get the CI array
    call loadCI(CI,I)
  end if

  ! compute 1-particle active density matrix GAMMA1
# ifdef _DMRG_
  if (DMRG) then
    ! set state number here because in poly1 we have no reference
    ! to which state we are computing
    if (iPrGlb >= DEBUG) write(u6,*) 'STINI setting DMRG state number to ',mstate(i)-1
    ! Convert to the root number despite having
    ! set only the checkpoint file paths for the desired state(s)
    call qcmaquis_interface_set_state(int(mstate(i)-1,c_int))
  end if
# endif

  ! Compute 1-particle active density matrix GAMMA1
  call POLY1(CI,nConf)
  ! Restructure GAMMA1 as DREF array, but keep it in DMIX
  call GETDREF(DREF,size(DREF))

  ! Loop over states to compute the contribution of state I to states J
  do J=1,Nstate
    ! Retrieve the weight of the contribution of state I to the density
    ! of state J
    wij = DWGT(I,J)
    ! Multiply density of state I with weight wij and add it to whatever
    ! is already in DMIX (contributions of other states already computed)
    ! and store it in DMIX
    DMIX(:,J) = DMIX(:,J)+wij*DREF(:)
    !call daxpy_(SIZE(DREF),wij,DREF,1,DMIX(:,J),1)
  end do

  ! End of long loop over states
end do

! Deallocate everything
call mma_deallocate(CMO_Internal)
nullify(CMO)
call mma_deallocate(CI)

end subroutine rdminit

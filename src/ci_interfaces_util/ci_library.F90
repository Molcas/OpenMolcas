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
! Copyright (C) 2026, Roland Lindh                                     *
!***********************************************************************

!> The CI library of the CI interface, the module the programs use: the
!> generic CI operations CI_Lib_* of the core (CI_Operations, see there),
!> with the checks of MOLCAS_CI_VERIFY (CI_Lib_Checks) around those that
!> are verified. The checks apply to the CI state set up by
!> CI_Lib_Initialize from the CI input (CI_Lib_Input_State); a state the
!> program set up itself gets the core operations unchanged.
module CI_Library

use CI_Operations, Raw_Transform => CI_Lib_Transform

implicit none
private :: Raw_Transform

public :: CI_Lib_Transform

contains

!> The CI vector CI (basis of the selected library, symmetry iSym of state
!> iState) for the active orbitals rotated by TU (per symmetry the square
!> nAsh x nAsh blocks): the transformation of the core (CI_Lib_Transform of
!> CI_Operations). Under MOLCAS_CI_VERIFY=SIGMA compared with Lucia (SGUGA,
!> Faroald) and with SGUGA (the other libraries).
subroutine CI_Lib_Transform(iState,iSym,nTU,TU,n,CI)

  use CI_Lib_Checks, only: CI_Basis_To_File, CI_Basis_To_Lucia, CI_File_To_Basis, CI_Lucia_To_Basis, Faro_CI, SG_Avail, &
                           SG_CI, Tim_Off, Tim_On, Tim_TraPsi, TraPsi_Faro_Verified, TraPsi_Verified
  use CI_Verify_Options, only: CI_Verify
  use Lucia_Interface, only: Lucia_Util
  use sguga, only: sg_traci
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: One
  use Definitions, only: wp, iwp

  integer(kind=iwp), intent(in) :: iState, iSym, nTU, n
  real(kind=wp), intent(in) :: TU(nTU)
  real(kind=wp), intent(inout) :: CI(n)
  real(kind=wp) :: Dev, TimT(2)
  logical(kind=iwp) :: Ref_Lucia, Verify_SG
  real(kind=wp), allocatable :: CI_Ref(:), CI_SGUGA(:), ctemp(:)

  call Tim_On(TimT)
  Ref_Lucia = CI_Verify('SIGMA') .and. CI_Lib_Input_State(iState) .and. (SG_CI() .or. Faro_CI())
  Verify_SG = CI_Verify('SIGMA') .and. CI_Lib_Input_State(iState) .and. SG_Avail() .and. (.not. SG_CI())
  if (Ref_Lucia) then
    ! the reference of Lucia, in the Lucia basis
    call mma_allocate(CI_Ref,n,Label='CI_Lucia')
    call CI_Basis_To_Lucia(iState,iSym,n,CI,CI_Ref)
    call Lucia_Util('TRACI',CI_Vector=CI_Ref(:),TU=TU(:),nTU=nTU)
  end if
  if (Verify_SG) then
    ! the reference of SGUGA, in the order of the file
    call mma_allocate(CI_SGUGA,n,Label='CI_SGUGA')
    call CI_Basis_To_File(iState,iSym,n,CI,CI_SGUGA)
    call SG_TRACI(iState,iSym,TU,nTU,CI_SGUGA,n)
  end if

  call Raw_Transform(iState,iSym,nTU,TU,n,CI)

  if (Ref_Lucia) then
    call mma_allocate(ctemp,n,Label='ctemp')
    call CI_Lucia_To_Basis(iState,iSym,n,CI_Ref,ctemp)
    Dev = maxval(abs(CI(:)-ctemp(:)))
    if (SG_CI()) then
      call Report_Dev('Mk_Tra_Psi (SGCI): SGUGA CI transformation vs Lucia',Dev,1.0e-10_wp,TraPsi_Verified)
    else
      call Report_Dev('Mk_Tra_Psi (FARO): Faroald CI transformation vs Lucia',Dev,1.0e-10_wp,TraPsi_Faro_Verified)
    end if
    call mma_deallocate(ctemp)
    call mma_deallocate(CI_Ref)
  end if
  if (Verify_SG) then
    call mma_allocate(ctemp,n,Label='ctemp')
    call CI_File_To_Basis(iState,iSym,n,CI_SGUGA,ctemp)
    call Report_Dev('Mk_Tra_Psi: SGUGA CI transformation vs the CI library',maxval(abs(CI(:)-ctemp(:))), &
                    1.0e-10_wp*max(One,maxval(abs(CI(:)))),TraPsi_Verified)
    call mma_deallocate(ctemp)
    call mma_deallocate(CI_SGUGA)
  end if
  call Tim_Off(Tim_TraPsi,TimT)

end subroutine CI_Lib_Transform

end module CI_Library

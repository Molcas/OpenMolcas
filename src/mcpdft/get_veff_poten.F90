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
! Copyright (C) 2024, Matthew R. Hennefarth                            *
!               2026, Bhavnesh Jangid                                  *
!***********************************************************************
!  get_veff_poten
!
!> @brief get the effective two-electron potential in AO basis
!>
!> @details
!>   The Fock builder generates J - 0.5*K when ExFac=1 for the spin-summed
!>   RASSCF densities.  This wraps that path and returns the combined
!>   effective potential directly.
!>
!> @param[in]  cmo      MO coefficients
!> @param[in]  dm1_core core density in AO basis
!> @param[in]  dm1_cas  active space density in AO basis
!> @param[out] veff     effective potential in AO basis
!***********************************************************************

subroutine get_veff_poten(cmo,dm1_core,dm1_cas,veff)

use PrintLevel, only: DEBUG, INSANE
use rasscf_global, only: lsquare, NACPR2, nfint
use general_data, only: ntot1, ntot2
use mcpdft_output, only: iprglb
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One, Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: cmo(ntot2), dm1_core(ntot2), dm1_cas(ntot2)
real(kind=wp), intent(out) :: veff(ntot1)
integer(kind=iwp) :: print_level
real(kind=wp), allocatable :: focka(:), focki(:), puvx(:), tuvx(:)

select case (iPrGlb)
  case (DEBUG)
    print_level = 5
  case (INSANE)
    print_level = 10
  case default
    print_level = 0
end select

call mma_allocate(puvx,nfint,label='puvx')
call mma_allocate(tuvx,nacpr2,label='tuvx')
call mma_allocate(focki,ntot1,label='focki')
call mma_allocate(focka,ntot1,label='focka')

puvx(:) = Zero
tuvx(:) = Zero
focki(:) = Zero
focka(:) = Zero
call tractl2(cmo,puvx,tuvx,dm1_core,focki,dm1_cas,focka,print_level,lsquare,One)

veff(:ntot1) = focki(:)+focka(:)

call mma_deallocate(focki)
call mma_deallocate(focka)
call mma_deallocate(tuvx)
call mma_deallocate(puvx)

end subroutine get_veff_poten

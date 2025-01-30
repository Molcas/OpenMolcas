!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2024, Matthew R. Hennefarth                            *
!***********************************************************************

!> @brief get coulomb integrals in AO basis
!>
!> @details
!>   Coulomb integrals generated from \p dm1_core and \p dm1_cas. Currently implementation uses the Fock matrices with ExFac set to
!    zero. As such, it just wraps the TRACTL2 call.
!>
!> @author Matthew R. Hennefarth
!>
!> @param[in] cmo MO coefficients
!> @param[in] dm1_core core density in AO basis
!> @param[in] dm1_cas active space density in AO basis
!> @param[in] coul coulomb integrals in AO basis
subroutine get_coulomb(cmo,dm1_core,dm1_cas,coul)
  use definitions,only:iwp,wp
  use constants,only:zero
  use printlevel,only:debug,insane
  use stdalloc,only:mma_allocate,mma_deallocate
  use rasscf_global,only:lsquare,NACPR2,nfint
  use general_data,only:ntot1
  use mcpdft_output,only:iprglb
  implicit none

  real(kind=wp),intent(in) :: cmo(*),dm1_core(*),dm1_cas(*)
  real(kind=wp),intent(out) :: coul(*)

  integer(kind=iwp) :: print_level
  real(kind=wp),parameter :: exfac = zero
  real(kind=wp),allocatable :: puvx(:),tuvx(:),focki(:),focka(:)

  select case(iPrGlb)
  case(debug)
    print_level = 5
  case(insane)
    print_level = 10
  case default
    print_level = 0
  endselect

  call mma_allocate(puvx,nfint,label='puvx')
  call mma_allocate(tuvx,nacpr2,label='tuvx')
  call mma_allocate(focki,ntot1,label='focki')
  call mma_allocate(focka,ntot1,label='focka')

  ! Quite honestly, this is overkill. It generates the inactive fock, active fock, and transformed integrals puvx/tuvx etc. This
  ! is completely not necessary for what we need.
  call tractl2(cmo,puvx,tuvx,dm1_core,focki,dm1_cas,focka,print_level,lsquare,exfac)

  coul(:ntot1) = focki(:)+focka(:)

  call mma_deallocate(puvx)
  call mma_deallocate(tuvx)
  call mma_deallocate(focki)
  call mma_deallocate(focka)
endsubroutine get_coulomb

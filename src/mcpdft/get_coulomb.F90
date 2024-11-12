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
!    zero.
!>
!> @author Matthew R. Hennefarth
!>
!> @param[in] cmo MO coefficients
!> @param[in] dm1_core core density in AO basis
!> @param[in] dm1_cas active space density in AO basis
!> @param[in] coul coulomb integrals in AO basis
subroutine get_coulomb(cmo,dm1_core,dm1_cas,coul)
  use definitions,only:wp
  use constants,only:zero
  use stdalloc,only:mma_allocate,mma_deallocate
  use rasscf_global,only:exfac,ipr,lsquare,NACPR2,nfint
  use general_data,only:ntot1
  implicit none

  real(kind=wp),intent(in) :: cmo(*),dm1_core(*),dm1_cas(*)
  real(kind=wp),intent(out) :: coul(*)

  real(kind=wp),allocatable :: puvx(:),tuvx(:),focki(:),focka(:)

  call mma_allocate(puvx,nfint,label='puvx')
  call mma_allocate(tuvx,nacpr2,label='tuvx')
  call mma_allocate(focki,ntot1,label='focki')
  call mma_allocate(focka,ntot1,label='focka')

  call tractl2(cmo,puvx,tuvx,dm1_core,focki,dm1_cas,focka,ipr,lsquare,exfac)

  coul(:ntot1) = focki(:)+focka(:)

  call mma_deallocate(puvx)
  call mma_deallocate(tuvx)
  call mma_deallocate(focki)
  call mma_deallocate(focka)
endsubroutine get_coulomb

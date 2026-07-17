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
! Copyright (C) 1996, Markus P. Fuelscher                              *
!***********************************************************************

!***********************************************************************
! History:                                                             *
! Matthew R. Hennefarth May 31, 2023, upgraded to F90                  *
!***********************************************************************

subroutine cas_mo_to_ao(CMO, CAS_MO, CAS_AO)
  ! Convert CAS_MO from active space-MO to AO basis. Essentially
  ! performs CMOA.transpose() * CAS_MO * CMOA where CMOA is the active
  ! space only part of the MO coefficients
  !
  ! Args:
  !   CMO: ndarray (ntot2)
  !     MO coefficients
  !   CAS_MO: ndarray (nacpar)
  !     Matrix definied in MO active space only
  !
  ! Returns:
  !   CAS_AO: ndarray (ntot2)
  !     Matrix converted into AO basis
  use definitions, only: wp
  use stdalloc, only: mma_allocate, mma_deallocate
  implicit none

  ! For ntot2, nsym, nbas, nash, nish, nfro
#include "rasdim.fh"
#include "general.fh"
  ! For nacpar
#include "rasscf.fh"
  real(kind=wp), parameter :: zero = 0.0D0
  real(kind=wp), dimension(ntot2), intent(in) :: CMO
  real(kind=wp), dimension(NACPAR), intent(in) :: CAS_MO
  real(kind=wp), dimension(ntot2), intent(out) :: CAS_AO

  real(kind=wp), allocatable :: tmp1(:), tmp2(:)
  integer :: iSym, nBasis, nActive, nInactive, nFrozen
  integer, dimension(2) :: offset
  offset = 1

  do iSym=1, nSym
    nBasis = nbas(iSym)
    nActive = nash(iSym)
    nInactive = nIsh(iSym)
    nFrozen = nFro(iSym)
    call dCopy_(nBasis*nBasis, [zero], 0, CAS_AO(offset(2)), 1)
    if(nActive /= 0) then
      call mma_allocate(tmp1, nActive*nActive, label='Scr1')
      call mma_allocate(tmp2, nActive*nBasis, label='Scr2')
      ! This is needed since CAS_MO is stored as upper triangle
      ! and we need iTmp1 to be a full square matrix.
      call square(CAS_MO(offset(1)), tmp1, 1, nactive, nactive)
      ! Does essentially C_a * D1A * C_a converting from MO -> AO basis
      call dgemm_('n','t', nbasis, nactive, nactive, 1.0D0, &
              CMO(offset(2)+(nFrozen+nInactive)*nBasis), nbasis, tmp1, &
              nActive, zero, tmp2, nbasis)
      call dgemm_('n', 't', nbasis, nbasis, nactive, 1.0d0, tmp2, &
              nbasis, CMO(offset(2) + (nFrozen+nInactive)*nBasis), &
              nbasis, zero, CAS_AO(offset(2)), nbasis)

      call mma_deallocate(tmp1)
      call mma_deallocate(tmp2)
    end if
    offset(1) = offset(1) + (nActive**2+nactive)/2
    offset(2) = offset(2) + nbasis**2
  end do

end subroutine cas_mo_to_ao
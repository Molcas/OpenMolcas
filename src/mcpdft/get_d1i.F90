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

subroutine get_d1i_mcpdft(cmo, d1i)
  ! Computes the one-body density in the ao basis from the frozen and
  ! inactive orbitals.
  !
  ! Args:
  !   CMO: ndarray (ntot2)
  !     MO coefficients
  !
  ! Returns:
  !   D1I: ndarray (ntot2)
  !     Frozen/Inactive one-body reduced density matrix in AO basis
  use definitions, only: wp
  implicit none
#include "rasdim.fh"
#include "general.fh"

  real(kind=wp), dimension(ntot2) :: cmo
  real(kind=wp), dimension(ntot2) :: d1i

  real(kind=wp), parameter :: zero=0.0d0
  real(kind=wp), parameter :: two=2.0d0

  integer :: isym, nb, nfi, nbsq, ista

  ista = 1
  do isym=1,nsym
    nb = nbas(isym)
    nbsq = nb**2
    nfi = nfro(isym) + nish(isym)
    if(nb > 0) then
      call dcopy_(nbsq, [zero], 0, d1i(ista), 1)
      if(nfi > 0) then
        ! D1I = 2* C \dot C within the correct symmetry
        call DGEMM_('n', 't', nb, nb, nfi, two, &
             cmo(ista), nb, cmo(ista), nb, zero, d1i(ista), nb)
      end if
      ista = ista + nbsq
    end if
  end do
end subroutine get_d1i_mcpdft
!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

! The diagonal one-electron integrals of the active orbitals, in Lucia's
! (type) order of the orbitals, in which CIDIA4 addresses them (IREOTS:
! type -> symmetry order)
subroutine GT1DIA_MCLR(H1DIA)

use MCLR_Data, only: FIMO, ipCM, IREOTS, NACOB
use Orbital_Spaces, only: nAsh, nIsh, nSym
use input_mclr, only: nOrb
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: H1DIA(*)
integer(kind=iwp) :: i, iAsh, ii, iS
real(kind=wp), allocatable :: HSym(:)

call mma_allocate(HSym,NACOB,Label='HSym')
i = 1
do iS=1,nSym
  ii = ipCM(iS)+nOrb(iS)*(nIsh(iS)-1)+nIsh(iS)-1
  do iAsh=1,nAsh(iS)
    ii = ii+nOrb(iS)+1
    HSym(i) = FIMO(ii)
    i = i+1
  end do
end do
do i=1,NACOB
  H1DIA(i) = HSym(IREOTS(i))
end do
call mma_deallocate(HSym)

end subroutine GT1DIA_MCLR

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

subroutine GT1DIA_MCLR(H1DIA)

use MCLR_Data, only: FIMO, ipCM
use input_mclr, only: nAsh, nIsh, nOrb, nSym
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: H1DIA(*)
integer(kind=iwp) :: i, iAsh, ii, iS

i = 1
do iS=1,nSym
  ii = ipCM(iS)+nOrb(iS)*(nIsh(iS)-1)+nIsh(iS)-1
  do iAsh=1,nAsh(iS)
    ii = ii+nOrb(iS)+1
    H1DIA(i) = FIMO(ii)
    i = i+1
  end do
end do

end subroutine GT1DIA_MCLR

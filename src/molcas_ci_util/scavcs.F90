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

subroutine SCAVCS(VECO,VECI,IDX,NDIM)
! Scatter vector with sign change
!
! vecO(abs(idx(i))) = veci(i)*sign(idx(i))

use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NDIM, IDX(NDIM)
real(kind=wp), intent(_OUT_) :: VECO(*)
real(kind=wp), intent(in) :: VECI(NDIM)
integer(kind=iwp) :: I

do I=1,NDIM
  VECO(abs(IDX(I))) = VECI(I)*real(sign(1,IDX(I)),kind=wp)
end do

return

end subroutine SCAVCS

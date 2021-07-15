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

subroutine GATVCS(VECO,VECI,IDX,NDIM)
! Gather vector allowing for sign change
!
! VECO(I) = VECI(IDX(I))

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NDIM, IDX(NDIM)
real(kind=wp), intent(out) :: VECO(NDIM)
real(kind=wp), intent(in) :: VECI(*)
integer(kind=iwp) :: I

do I=1,NDIM
  VECO(I) = VECI(abs(IDX(I)))*real(sign(1,IDX(I)),kind=wp)
end do

return

end subroutine GATVCS

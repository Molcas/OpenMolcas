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

subroutine VVTOV(VECIN1,VECIN2,VECUT,NDIM)
! VECUT(I) = VECIN1(I) * VECIN2(I)

use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: VECIN1(*), VECIN2(*), VECUT(*)
integer(kind=iwp) :: NDIM
integer(kind=iwp) :: I

do I=1,NDIM
  VECUT(I) = VECIN1(I)*VECIN2(I)
end do

end subroutine VVTOV

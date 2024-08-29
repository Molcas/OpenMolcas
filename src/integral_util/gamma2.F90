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
! Copyright (C) 1992, Roland Lindh                                     *
!***********************************************************************

function Gamma2(m,T)
!***********************************************************************
!                                                                      *
! Object: to compute the auxiliary function in the high argument       *
!         approximation.                                               *
!                                                                      *
!***********************************************************************

use Constants, only: Two, Half, Pi
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: Gamma2
integer(kind=iwp), intent(in) :: m
real(kind=wp) :: T
integer(kind=iwp) :: i

Gamma2 = Half*sqrt(Pi/T)
do i=1,m
  Gamma2 = (real(2*i-1,kind=wp)/(Two*T))*Gamma2
end do

return

end function Gamma2

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

function G(Arg)

use Constants, only: One, Half, Quart, Pi
use Definitions, only: wp

implicit none
real(kind=wp) :: G
real(kind=wp), intent(in) :: Arg
real(kind=wp) :: Arg_, rG

g = -1.0e3_wp

Arg_ = real(int(Arg),kind=wp)
if (abs(Arg-Arg_) < Quart) then
  ! Integer argument
  G = One
  rG = One
else
  ! Fractional argument
  G = sqrt(Pi)
  rG = Half
end if

do
  if (abs(rG-Arg) < Quart) exit
  G = rG*G
  rG = rG+One
end do

return

end function G

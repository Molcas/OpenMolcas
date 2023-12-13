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

function GAM(M)

use Constants, only: One, Half, Pi
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: GAM
integer(kind=iwp), intent(in) :: M
integer(kind=iwp) :: I, MA
real(kind=wp) :: G

if (mod(M,2) == 1) then
  MA = (M+1)/2
  G = One
  if (MA /= 1) then
    do I=2,MA
      G = G*real(I-1,kind=wp)
    end do
  end if
else
  MA = M
  G = sqrt(Pi)
  if (MA /= 0) then
    do I=1,MA,2
      G = G*Half*real(I,kind=wp)
    end do
  end if
end if
GAM = G

return

end function GAM

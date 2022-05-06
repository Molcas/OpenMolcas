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

implicit real*8(a-h,o-z)
#include "real.fh"
real*8 G

g = -1000.0d0

Arg_ = dble(int(Arg))
if (abs(Arg-Arg_) < Half/Two) then
  ! Integer argument
  G = One
  rG = One
else
  ! Fractional argument
  G = sqrt(Pi)
  rG = Half
end if

99 continue
if (abs(rG-Arg) < Half/Two) goto 666
G = rG*G
rG = rG+One
Go To 99
666 continue

return

end function G

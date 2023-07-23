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
! Copyright (C) 2007, Ten-no Research Group                            *
!               2012, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine LVCLR(A,INCA,N)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: INCA, N
real(kind=wp) :: A(*)
integer(kind=iwp) :: L, LA

! ----- ZERO OUT VECTOR -A-, USING INCREMENT -INCA- -----

if (INCA /= 1) GO TO 100
do L=1,N
  A(L) = Zero
end do

return

100 continue
LA = 1-INCA
do L=1,N
  LA = LA+INCA
  A(LA) = Zero
end do

return

end subroutine LVCLR

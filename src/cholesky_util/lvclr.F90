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

implicit real*8(A-H,O-Z)
real*8 A(*)
parameter(ZERO=0.0D+00)

! ----- ZERO OUT VECTOR -A-, USING INCREMENT -INCA- -----

if (INCA /= 1) GO TO 100
do L=1,N
  A(L) = ZERO
end do

return

100 continue
LA = 1-INCA
do L=1,N
  LA = LA+INCA
  A(LA) = ZERO
end do

return

end subroutine LVCLR

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

subroutine ZEROMA(W,I1,I2)

implicit none
integer I1, I2, I
real*8 ZERO, W
parameter(ZERO=0.d0)
dimension W(*)

if (I2 < I1) return
do I=I1,I2
  W(I) = ZERO
end do

return

end subroutine ZEROMA

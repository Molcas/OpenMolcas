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

function CHO_IDOT(N,IVEC1,INC1,IVEC2,INC2)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: CHO_IDOT
integer(kind=iwp) :: N, IVEC1(*), INC1, IVEC2(*), INC2
integer(kind=iwp) :: I, I1, I2, ISUM

CHO_IDOT = 0
if (N < 1) return
ISUM = 0
if ((INC1 == 1) .and. (INC2 == 1)) then
  do I=1,N
    ISUM = ISUM+IVEC1(I)*IVEC2(I)
  end do
else
  I1 = 1
  I2 = 1
  if (INC1 < 0) I1 = (-N+1)*INC1+1
  if (INC2 < 0) I2 = (-N+1)*INC2+1
  do I=1,N
    ISUM = ISUM+IVEC1(I1)*IVEC2(I2)
    I1 = I1+INC1
    I2 = I2+INC2
  end do
end if
CHO_IDOT = ISUM

end function CHO_IDOT

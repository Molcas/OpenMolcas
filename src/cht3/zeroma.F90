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

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: W(*)
integer(kind=iwp) :: I1, I2
integer(kind=iwp) :: I

if (I2 < I1) return
do I=I1,I2
  W(I) = Zero
end do

return

end subroutine ZEROMA

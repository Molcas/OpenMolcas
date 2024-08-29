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

subroutine bino(lmax)

use welcom, only: binom
use Constants, only: Zero, One
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lMax
integer(kind=iwp) :: i, j

binom(:,:) = Zero
binom(0,0) = One
do i=1,lmax
  do j=0,i
    binom(i,j) = binom(i-1,j-1)+binom(i-1,j)
  end do
end do

return

end subroutine bino

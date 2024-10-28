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

subroutine ylmnor(lmax)

use welcom, only: aNorm, FiInt, TetInt
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lmax
integer(kind=iwp) :: i, lm2, j, k
real(kind=wp) :: Tal

do i=0,lmax
  lm2 = i/2
  do j=0,lm2
    do k=0,j
      anorm(i,j,k) = fiint(j-k,k)*tetint(i,j)
    end do
  end do
end do
do i=0,lmax
  tal = One/anorm(i,0,0)
  lm2 = i/2
  do j=0,lm2
    anorm(i,j,0:j) = anorm(i,j,0:j)*tal
  end do
end do

return

end subroutine ylmnor

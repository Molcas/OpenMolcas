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

subroutine fiin(lmax)

use Constants, only: Zero, One, Two, Pi
use welcom, only: fiInt, binom
use Definitions, only: wp

implicit none
integer lmax
integer i, j, k, l, iexp
real*8 a, tal, al

fiint(0,0) = Pi*Two
do i=0,lmax
  do j=0,lmax-i
    fiint(i,j) = Zero
    do k=0,j
      a = binom(j,k)
      iexp = i+k
      tal = pi*Two*a*(-One)**k
      if (iexp /= 0) then
        do l=1,iexp
          al = Two*real(l,kind=wp)
          tal = tal*(al-One)/al
        end do
      end if
      fiint(i,j) = fiint(i,j)+tal
    end do
  end do
end do

return

end subroutine fiin

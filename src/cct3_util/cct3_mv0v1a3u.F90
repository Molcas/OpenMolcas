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

subroutine cct3_mv0v1a3u(rowa,cola,ddx,ddy,nopi,nopj,incx,incy,a,x,y)
! Y(iy) = Y(iy) + A * X(ix)

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: rowa, cola, ddx, ddy, nopi, nopj, incx, incy
real(kind=wp) :: a(rowa,cola), x(ddx), y(ddy)
#include "t31.fh"
integer(kind=iwp) :: i, ix, iy, j

if (mhkey == 1) then
  ! ESSL
  !call dgemx(nopi,nopj,One,a,rowa,x,incx,y,incy)
  call dgemv_('N',nopi,nopj,One,a,rowa,x,incx,One,y,incy)

else
  ! Fortran matrix handling

  if ((incx == 1) .and. (incy == 1)) then

    ! Inc's = 1

    do j=1,nopj
      do i=1,nopi

        y(i) = y(i)+a(i,j)*x(j)
      end do
    end do

  else

    ! Other type inc's

    ix = 1
    do j=1,nopj
      iy = 1
      do i=1,nopi
        y(iy) = a(i,j)*x(ix)+y(iy)
        iy = iy+incy
      end do
      ix = ix+incx
    end do

  end if

end if

return

end subroutine cct3_mv0v1a3u

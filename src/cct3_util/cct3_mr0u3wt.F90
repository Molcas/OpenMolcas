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
       subroutine cct3_mr0u3wt                                          &
     & (ddx, ddy,                                                       &
     & nop, incx, incy,                                                 &
     & x, y,                                                            &
     & scal)
!
!     scalar = sum (x(ix) * y(iy))
!
#include "t31.fh"
       integer ddx, ddy
       integer nop, incx, incy
       real*8  x(1:ddx), y(1:ddy)
       real*8  scal
       real*8 ddot_
!
!      help vraiables
!
       integer i, ix, iy
!
       if (mhkey.eq.1) then
!      ESSL
       scal=ddot_(nop,x,incx,y,incy)
!
       else
!      Fortran matrix handling
!
!      return for no operations
!
       scal = 0.0d0
       if (nop.le.0) return
!
       if (incx.eq.1.and.incy.eq.1) then
!
!      inc's = 1
!
       do 10 i = 1, nop
       scal = scal + x(i)*y(i)
 10    continue
!
       else
!
!      other type increments
!
       ix = 1
       iy = 1
       if (incx.lt.0) ix = 1-(nop-1)*incx
       if (incy.lt.0) iy = 1-(nop-1)*incy
!
       do 20 i = 1,nop
       scal = scal + x(ix)*y(iy)
       ix = ix + incx
       iy = iy + incy
 20    continue
!
       end if
!
       end if
!
       return
       end

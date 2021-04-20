************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
       subroutine mr0u3wt
     & (ddx, ddy,
     & nop, incx, incy,
     & x, y,
     & scal)
C
C     scalar = sum (x(ix) * y(iy))
C
#include "ccsd1.fh"
       integer ddx, ddy
       integer nop, incx, incy
       real*8  x(1:ddx), y(1:ddy)
       real*8  scal
       real*8 ddot_
C
C      help vraiables
C
       integer i, ix, iy
C
       if (mhkey.eq.1) then
c      ESSL
       scal=ddot_(nop,x,incx,y,incy)
C
       else
c      Fortran matrix handling
C
C      return for no operations
C
       scal = 0.0d0
       if (nop.le.0) return
C
       if (incx.eq.1.and.incy.eq.1) then
C
C      inc's = 1
C
       do 10 i = 1, nop
       scal = scal + x(i)*y(i)
 10    continue
C
       else
C
c      other type increments
C
       ix = 1
       iy = 1
       if (incx.lt.0) ix = 1-(nop-1)*incx
       if (incy.lt.0) iy = 1-(nop-1)*incy
C
       do 20 i = 1,nop
       scal = scal + x(ix)*y(iy)
       ix = ix + incx
       iy = iy + incy
 20    continue
C
       end if
C
       end if
c
       return
       end

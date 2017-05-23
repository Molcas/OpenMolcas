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
      Subroutine ICOPY(n,x,incX,y,incY)
************************************************************************
*                                                                      *
*     Copy vector X into vector Y                                      *
*                                                                      *
*     calling arguments:                                               *
*     N       : Integer, input.                                        *
*               Number of input elements.                              *
*     X       : Array of Integer, input.                               *
*               Input vector, X.                                       *
*     incX    : Integer, input.                                        *
*               Stride of vector X.                                    *
*     Y       : Array of Integer, output.                              *
*               Output vector, Y.                                      *
*     incY    : Integer, input.                                        *
*               Stride of vector Y.                                    *
*                                                                      *
************************************************************************
c
      Implicit Integer (A-Z)
c
c     copies integer vector, x, to integer vector, y.
c     uses unrolled loops for increments equal to one.
c
      integer x(*),y(*)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        y(iy) = x(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        y(i) = x(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        y(i) = x(i)
        y(i + 1) = x(i + 1)
        y(i + 2) = x(i + 2)
        y(i + 3) = x(i + 3)
        y(i + 4) = x(i + 4)
        y(i + 5) = x(i + 5)
        y(i + 6) = x(i + 6)
   50 continue
      return
      end

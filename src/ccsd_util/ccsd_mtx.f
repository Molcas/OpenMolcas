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
c
c     mv0zero
c     mc0c1a3b
c     mv0v1a3u
c     mr0u3wt
c     mc0c1at3b
c
c     ---------------
c
       subroutine mv0zero (dd,length,mat)
C
c      mat = 0
c
#include "ccsd1.fh"
       integer dd
       integer length
       real*8  mat(1:dd)
C
C      help variables
C
       integer init
       real*8  zero
       data    zero/0.0d0/
C
C
       if (mhkey.eq.1) then
c      ESSL

       call dcopy_(length,zero,0,mat,1)
c
       else
c      Fortran matrix handling
C
       do 10 init=1,length
       mat(init) = zero
 10    continue
C
       end if
c
       return
       end
C
C-------------------------------------
c
       subroutine mc0c1a3b
     & (rowa,cola,rowb,colb,rowc,colc,
     & row,sum,col,a,b,c)
C
C     C = C + A*B
C
#include "ccsd1.fh"
       integer rowa,cola,rowb,colb,rowc,colc
       integer row,sum,col
       real*8  A(1:rowa,1:cola)
       real*8  B(1:rowb,1:colb)
       real*8  C(1:rowc,1:colc)
C
C      help variables
C
       integer i,j,k
c
       if (mhkey.eq.1) then
c      ESSL
       call DGEMM_('N','N',row,col,sum,1.0d0,a,rowa,b,rowb,
     & 1.0d0,c,rowc)
C
       else
c      Fortran matrix handling
C
       do 50 j=1, col
       do 40 k=1, sum
       do 30 i=1, row
       c(i,j) = c(i,j) + a(i,k)*b(k,j)
 30    continue
 40    continue
 50    continue
C
       end if
c
       return
       end
C
C-------------------------------------
C
       subroutine mv0v1a3u
     & (rowa,cola,ddx,ddy,
     & nopi,nopj,incx,incy,
     & a,x,y)
C
C     Y(iy) = Y(iy) + A * X(ix)
C
#include "ccsd1.fh"
       integer rowa,cola
       integer ddx, ddy
       integer nopi, nopj, incx, incy
       real*8  a(1:rowa,1:cola)
       real*8  x(1:ddx), y(1:ddy)
C
C      help variables
C
       integer i, j, ix, iy
c
       if (mhkey.eq.1) then
c      ESSL
c      call  dgemx(nopi,nopj,1.0d0,a,rowa,x,incx,y,incy)
       call dgemv_('N',nopi,nopj,1.0d0,a,rowa,x,incx,1.0d0,y,incy)
C
       else
c      Fortran matrix handling
C
C
       if (incx.eq.1.and.incy.eq.1) then
C
C      Inc's = 1
C
       do 30 j = 1, nopj
       do 20 i = 1, nopi
C
       y(i) = y(i) + a(i,j)*x(j)
 20    continue
 30    continue
C
       else
C
C      Other type inc's
c
       ix = 1
       do 60 j = 1, nopj
       iy = 1
       do 50 i = 1, nopi
       y(iy) = a(i,j)*x(ix) + y(iy)
       iy = iy + incy
 50    continue
       ix = ix + incx
 60    continue
C
       end if
C
       end if
c
       return
       end
C
C-------------------------------------
c
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
C
C---------------------------------------
C
       subroutine mc0c1at3b
     & (rowa,cola,rowb,colb,rowc,colc,
     & row,sum,col,
     & a,b,c)
C
C     C= A(T)*B
C
#include "ccsd1.fh"
       integer rowa,cola,rowb,colb,rowc,colc
       integer row,sum,col
       real*8  a(1:rowa,1:cola)
       real*8  b(1:rowb,1:colb)
       real*8  c(1:rowc,1:colc)
C
C      help variables
C
       integer i,j,k
C
       if (mhkey.eq.1) then
c      ESSL
       call DGEMM_('T','N',row,col,sum,1.0d0,a,rowa,b,rowb,
     & 1.0d0,c,rowc)
C
       else
c      Fotran
C
       do 30 j=1, col
       do 20 i=1, row
       do 10 k=1, sum
       c(i,j) = c(i,j) + a(k,i)*b(k,j)
 10    continue
 20    continue
 30    continue
C
       end if
c
       return
       end
C

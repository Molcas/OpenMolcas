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
c     mc0c2a3b
c     mv0v1a3u
c     mr0u3wt
c     mc0c1at3b
c     mv0u
c     mv0v1u
c     mv0sv
c
c     ---------------
c
       subroutine mc0c2a3b
     & (rowa,cola,rowb,colb,rowc,colc,
     & row,sum,col,a,b,c)
C
C     C = C + A*B
C
       implicit none
#include "chcc1.fh"
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
       call DGEMM_('N','N',row,col,sum,-1.0d0,a,rowa,b,rowb,
     & 1.0d0,c,rowc)
C
       else
c      Fortran matrix handling
C
       do 50 j=1, col
       do 40 k=1, sum
       do 30 i=1, row
       c(i,j) = c(i,j) - a(i,k)*b(k,j)
 30    continue
 40    continue
 50    continue
C
       end if
c
       return
       end
C
C---------------------------------------
C
       subroutine mv0u (length,X,incX,Y,incY)
C
C      Y = X
C
       implicit none
#include "chcc1.fh"
c
       integer length,incx,incy
       real*8  X(*)
       real*8  Y(*)
C
C      help variables
C
       integer i
C
        if (mhkey.eq.1) then
c      ESSL
       call dcopy_(length,X(1),incx,Y(1),incY)
C
       else
c      Fotran
C
       if (incx*incy.eq.1) then
c
         do 10 i=1, length
         y(i) = x(i)
 10      continue
       else
         do 20 i=0, length-1
         y(1+i*incy) = x(1+i*incx)
 20      continue
       end if
C
       end if
c
       return
       end
C
C---------------------------------------
C
       subroutine mv0v1u (length,X,incX,Y,incY,f)
C
C      Y = Y + f . X
C
       implicit none
#include "chcc1.fh"
c
       integer length,incx,incy
       real*8  f
       real*8  X(*)
       real*8  Y(*)
C
C      help variables
C
       integer i
C
       if (mhkey.eq.1) then
c      ESSL
       call daxpy_(length,f,X(1),incx,Y(1),incY)
C
       else
c      Fotran
C
       if (incx*incy.eq.1) then
c
         do 10 i=1, length
         y(i) = y(i)+f*x(i)
 10      continue
       else
         do 20 i=0, length-1
         y(1+i*incy) = y(1+i*incy)+f*x(1+i*incx)
 20      continue
       end if
C
       end if
c
       return
       end
c
c     ---------------
c
       subroutine mv0sv (dd,length,mat,f)
C
c      A = A .f
c
#include "chcc1.fh"
       integer dd
       integer length
       real*8  mat(1:dd)
       real*8 f
C
C      help variables
C
       integer init
C
C
       if (mhkey.eq.1) then
c      ESSL
c
       do  5 init=1,length
       mat(init) = mat(init)*f
  5    continue
c
       else
c      Fortran matrix handling
C
       do 10 init=1,length
       mat(init) = mat(init)*f
 10    continue
c
       end if
c
c
        return
        end

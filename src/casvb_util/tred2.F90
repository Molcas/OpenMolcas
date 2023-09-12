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
!  ********************************************************
!  ** Public-domain library routines used by casvb only. **
!  ********************************************************
!  **********************
!  ** EISPACK ROUTINES **
!  **********************
      subroutine casvb_tred2(nm,n,a,d,e,z)
!
      integer i,j,k,l,n,ii,nm,jp1
      real*8 a(nm,n),d(n),e(n),z(nm,n)
      real*8 f,g,h,hh,scale
!
!     this subroutine is a translation of the algol procedure tred2,
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine reduces a real symmetric matrix to a
!     symmetric tridiagonal matrix using and accumulating
!     orthogonal similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        a contains the real symmetric input matrix.  only the
!          lower triangle of the matrix need be supplied.
!
!     on output
!
!        d contains the diagonal elements of the tridiagonal matrix.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e(1) is set to zero.
!
!        z contains the orthogonal transformation matrix
!          produced in the reduction.
!
!        a and z may coincide.  if distinct, a is unaltered.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      do 100 i = 1, n
!
         do 80 j = i, n
         z(j,i) = a(j,i)
   80    continue
!
         d(i) = a(n,i)
  100 continue
!
      if (n .eq. 1) go to 510
!     .......... for i=n step -1 until 2 do -- ..........
      do 300 ii = 2, n
         i = n + 2 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 2) go to 130
!     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
         scale = scale + abs(d(k))
  120    continue
!
         if (scale .ne. 0.0d0) go to 140
  130    e(i) = d(l)
!
         do 135 j = 1, l
            d(j) = z(l,j)
            z(i,j) = 0.0d0
            z(j,i) = 0.0d0
  135    continue
!
         go to 290
!
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
!
         f = d(l)
         g = -sign(sqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
!     .......... form a*u ..........
         do 170 j = 1, l
         e(j) = 0.0d0
  170    continue
!
         do 240 j = 1, l
            f = d(j)
            z(j,i) = f
            g = e(j) + z(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
!
            do 200 k = jp1, l
               g = g + z(k,j) * d(k)
               e(k) = e(k) + z(k,j) * f
  200       continue
!
  220       e(j) = g
  240    continue
!     .......... form p ..........
         f = 0.0d0
!
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
!
         hh = f / (h + h)
!     .......... form q ..........
         do 250 j = 1, l
         e(j) = e(j) - hh * d(j)
  250    continue
!     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
!
            do 260 k = j, l
            z(k,j) = z(k,j) - f * e(k) - g * d(k)
  260       continue
!
            d(j) = z(l,j)
            z(i,j) = 0.0d0
  280    continue
!
  290    d(i) = h
  300 continue
!     .......... accumulation of transformation matrices ..........
      do 500 i = 2, n
         l = i - 1
         z(n,l) = z(l,l)
         z(l,l) = 1.0d0
         h = d(i)
         if (h .eq. 0.0d0) go to 380
!
         do 330 k = 1, l
         d(k) = z(k,i) / h
  330    continue
!
         do 360 j = 1, l
            g = 0.0d0
!
            do 340 k = 1, l
            g = g + z(k,i) * z(k,j)
  340       continue
!
            do 365 k = 1, l
               z(k,j) = z(k,j) - g * d(k)
  365       continue
  360    continue
!
  380    do 400 k = 1, l
         z(k,i) = 0.0d0
  400    continue
!
  500 continue
!
  510 do 520 i = 1, n
         d(i) = z(n,i)
         z(n,i) = 0.0d0
  520 continue
!
      z(n,n) = 1.0d0
      e(1) = 0.0d0
      return
      end

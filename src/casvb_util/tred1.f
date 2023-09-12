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
      subroutine casvb_tred1(nm,n,a,d,e,e2)
!
      integer i,j,k,l,n,ii,nm,jp1
      real*8 a(nm,n),d(n),e(n),e2(n)
      real*8 f,g,h,scale
!
!     this subroutine is a translation of the algol procedure tred1,
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine reduces a real symmetric matrix
!     to a symmetric tridiagonal matrix using
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
!        a contains information about the orthogonal trans-
!          formations used in the reduction in its strict lower
!          triangle.  the full upper triangle of a is unaltered.
!
!        d contains the diagonal elements of the tridiagonal matrix.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e(1) is set to zero.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2 may coincide with e if the squares are not needed.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      do 100 i = 1, n
         d(i) = a(n,i)
         a(n,i) = a(i,i)
  100 continue
!     .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 1) go to 130
!     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
         scale = scale + abs(d(k))
  120    continue
!
         if (scale .ne. 0.0d0) go to 140
!
         do 125 j = 1, l
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = 0.0d0
  125    continue
!
  130    e(i) = 0.0d0
         e2(i) = 0.0d0
         go to 300
!
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
!
         e2(i) = scale * scale * h
         f = d(l)
         g = -sign(sqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
         if (l .eq. 1) go to 285
!     .......... form a*u ..........
         do 170 j = 1, l
         e(j) = 0.0d0
  170    continue
!
         do 240 j = 1, l
            f = d(j)
            g = e(j) + a(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
!
            do 200 k = jp1, l
               g = g + a(k,j) * d(k)
               e(k) = e(k) + a(k,j) * f
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
         h = f / (h + h)
!     .......... form q ..........
         do 250 j = 1, l
         e(j) = e(j) - h * d(j)
  250    continue
!     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
!
            do 260 k = j, l
            a(k,j) = a(k,j) - f * e(k) - g * d(k)
  260       continue
!
  280    continue
!
  285    do 290 j = 1, l
            f = d(j)
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = f * scale
  290    continue
!
  300 continue
!
      return
      end

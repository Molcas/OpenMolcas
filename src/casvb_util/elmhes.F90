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
      subroutine elmhes(nm,n,low,igh,a,int)
!
      integer i,j,m,n,la,nm,igh,kp1,low,mm1,mp1
      real*8 a(nm,n)
      real*8 x,y
      integer int(igh)
!
!     this subroutine is a translation of the algol procedure elmhes,
!     num. math. 12, 349-368(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
!
!     given a real general matrix, this subroutine
!     reduces a submatrix situated in rows and columns
!     low through igh to upper hessenberg form by
!     stabilized elementary similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          subroutine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n.
!
!        a contains the input matrix.
!
!     on output
!
!        a contains the hessenberg matrix.  the multipliers
!          which were used in the reduction are stored in the
!          remaining triangle under the hessenberg matrix.
!
!        int contains information on the rows and columns
!          interchanged in the reduction.
!          only elements low through igh are used.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
!
      do 180 m = kp1, la
         mm1 = m - 1
         x = 0.0d0
         i = m
!
         do 100 j = m, igh
            if (abs(a(j,mm1)) .le. abs(x)) go to 100
            x = a(j,mm1)
            i = j
  100    continue
!
         int(m) = i
         if (i .eq. m) go to 130
!     .......... interchange rows and columns of a ..........
         do 110 j = mm1, n
            y = a(i,j)
            a(i,j) = a(m,j)
            a(m,j) = y
  110    continue
!
         do 120 j = 1, igh
            y = a(j,i)
            a(j,i) = a(j,m)
            a(j,m) = y
  120    continue
!     .......... end interchange ..........
  130    if (x .eq. 0.0d0) go to 180
         mp1 = m + 1
!
         do 160 i = mp1, igh
            y = a(i,mm1)
            if (y .eq. 0.0d0) go to 160
            y = y / x
            a(i,mm1) = y
!
            do 140 j = m, n
            a(i,j) = a(i,j) - y * a(m,j)
  140       continue
!
            do 150 j = 1, igh
            a(j,m) = a(j,m) + y * a(j,i)
  150       continue
!
  160    continue
!
  180 continue
!
  200 return
      end

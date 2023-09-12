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
      subroutine balbak(nm,n,low,igh,scale,m,z)
!
      integer i,j,k,m,n,ii,nm,igh,low
      REAL*8 scale(n),z(nm,m)
      REAL*8 s
!
!     this subroutine is a translation of the algol procedure balbak,
!     num. math. 13, 293-304(1969) by parlett and reinsch.
!     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
!
!     this subroutine forms the eigenvectors of a real general
!     matrix by back transforming those of the corresponding
!     balanced matrix determined by  balanc.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by  balanc.
!
!        scale contains information determining the permutations
!          and scaling factors used by  balanc.
!
!        m is the number of columns of z to be back transformed.
!
!        z contains the real and imaginary parts of the eigen-
!          vectors to be back transformed in its first m columns.
!
!     on output
!
!        z contains the real and imaginary parts of the
!          transformed eigenvectors in its first m columns.
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (m .eq. 0) go to 200
      if (igh .eq. low) go to 120
!
      do 110 i = low, igh
         s = scale(i)
!     .......... left hand eigenvectors are back transformed
!                if the foregoing statement is replaced by
!                s=1.0d0/scale(i). ..........
         do 100 j = 1, m
         z(i,j) = z(i,j) * s
  100    continue
!
  110 continue
!     ......... for i=low-1 step -1 until 1,
!               igh+1 step 1 until n do -- ..........
  120 do 140 ii = 1, n
         i = ii
         if (i .ge. low .and. i .le. igh) go to 140
         if (i .lt. low) i = low - ii
         k = Int(scale(i))
         if (k .eq. i) go to 140
!
         do 130 j = 1, m
            s = z(i,j)
            z(i,j) = z(k,j)
            z(k,j) = s
  130    continue
!
  140 continue
!
  200 return
      end

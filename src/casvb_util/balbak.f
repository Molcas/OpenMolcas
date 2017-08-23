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
c  ********************************************************
c  ** Public-domain library routines used by casvb only. **
c  ********************************************************
c  **********************
c  ** EISPACK ROUTINES **
c  **********************
      subroutine balbak(nm,n,low,igh,scale,m,z)
c
      integer i,j,k,m,n,ii,nm,igh,low
      REAL*8 scale(n),z(nm,m)
      REAL*8 s
c
c     this subroutine is a translation of the algol procedure balbak,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
c
c     this subroutine forms the eigenvectors of a real general
c     matrix by back transforming those of the corresponding
c     balanced matrix determined by  balanc.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by  balanc.
c
c        scale contains information determining the permutations
c          and scaling factors used by  balanc.
c
c        m is the number of columns of z to be back transformed.
c
c        z contains the real and imaginary parts of the eigen-
c          vectors to be back transformed in its first m columns.
c
c     on output
c
c        z contains the real and imaginary parts of the
c          transformed eigenvectors in its first m columns.
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (m .eq. 0) go to 200
      if (igh .eq. low) go to 120
c
      do 110 i = low, igh
         s = scale(i)
c     .......... left hand eigenvectors are back transformed
c                if the foregoing statement is replaced by
c                s=1.0d0/scale(i). ..........
         do 100 j = 1, m
  100    z(i,j) = z(i,j) * s
c
  110 continue
c     ......... for i=low-1 step -1 until 1,
c               igh+1 step 1 until n do -- ..........
  120 do 140 ii = 1, n
         i = ii
         if (i .ge. low .and. i .le. igh) go to 140
         if (i .lt. low) i = low - ii
         k = Int(scale(i))
         if (k .eq. i) go to 140
c
         do 130 j = 1, m
            s = z(i,j)
            z(i,j) = z(k,j)
            z(k,j) = s
  130    continue
c
  140 continue
c
  200 return
      end

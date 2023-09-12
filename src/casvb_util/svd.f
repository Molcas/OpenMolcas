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
      subroutine svd(nm,m,n,a,w,matu,u,matv,v,ierr,rv1)
!  TESTING FOR EQUALITY OR ZERO NOW DONE USING THRESHOLD, TT, 110100.
!
      integer i,j,k,l,m,n,ii,i1,kk,k1,ll,l1,mn,nm,its,ierr
      real*8 a(nm,n),w(n),u(nm,n),v(nm,n),rv1(n)
      real*8 c,f,g,h,s,x,y,z,tst1,tst2,scale,casvb_pythag,thresh
      logical matu,matv
      save thresh
      data thresh/1d-16/
!
!     this subroutine is a translation of the algol procedure svd,
!     num. math. 14, 403-420(1970) by golub and reinsch.
!     handbook for auto. comp., vol ii-linear algebra, 134-151(1971).
!
!     this subroutine determines the singular value decomposition
!          t
!     a=usv  of a real m by n rectangular matrix.  householder
!     bidiagonalization and a variant of the qr algorithm are used.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.  note that nm must be at least
!          as large as the maximum of m and n.
!
!        m is the number of rows of a (and u).
!
!        n is the number of columns of a (and u) and the order of v.
!
!        a contains the rectangular input matrix to be decomposed.
!
!        matu should be set to .true. if the u matrix in the
!          decomposition is desired, and to .false. otherwise.
!
!        matv should be set to .true. if the v matrix in the
!          decomposition is desired, and to .false. otherwise.
!
!     on output
!
!        a is unaltered (unless overwritten by u or v).
!
!        w contains the n (non-negative) singular values of a (the
!          diagonal elements of s).  they are unordered.  if an
!          error exit is made, the singular values should be correct
!          for indices ierr+1,ierr+2,...,n.
!
!        u contains the matrix u (orthogonal column vectors) of the
!          decomposition if matu has been set to .true.  otherwise
!          u is used as a temporary array.  u may coincide with a.
!          if an error exit is made, the columns of u corresponding
!          to indices of correct singular values should be correct.
!
!        v contains the matrix v (orthogonal) of the decomposition if
!          matv has been set to .true.  otherwise v is not referenced.
!          v may also coincide with a if u is not needed.  if an error
!          exit is made, the columns of v corresponding to indices of
!          correct singular values should be correct.
!
!        ierr is set to
!          zero       for normal return,
!          k          if the k-th singular value has not been
!                     determined after 30 iterations.
!
!        rv1 is a temporary storage array.
!
!     calls pythag for  sqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!

! The next two initializations are to appease a compiler
      f = 0.0d0
      h = 0.0d0

      ierr = 0
      l  = 0 ! dummy initialize
      l1 = 0 ! dummy initialize
!
      do 100 i = 1, m
!
         do 101 j = 1, n
            u(i,j) = a(i,j)
  101    continue
  100 continue
!     .......... householder reduction to bidiagonal form ..........
      g = 0.0d0
      scale = 0.0d0
      x = 0.0d0
!
      do 300 i = 1, n
         l = i + 1
         rv1(i) = scale * g
         g = 0.0d0
         s = 0.0d0
         scale = 0.0d0
         if (i .gt. m) go to 210
!
         do 120 k = i, m
         scale = scale + abs(u(k,i))
  120    continue
!
         if (abs(scale) .lt. thresh) go to 210
!
         do 130 k = i, m
            u(k,i) = u(k,i) / scale
            s = s + u(k,i)**2
  130    continue
!
         f = u(i,i)
         g = -sign(sqrt(s),f)
         h = f * g - s
         u(i,i) = f - g
         if (i .eq. n) go to 190
!
         do 150 j = l, n
            s = 0.0d0
!
            do 140 k = i, m
            s = s + u(k,i) * u(k,j)
  140       continue
!
            f = s / h
!
            do 160 k = i, m
               u(k,j) = u(k,j) + f * u(k,i)
  160       continue
  150    continue
!
  190    do 200 k = i, m
         u(k,i) = scale * u(k,i)
  200    continue
!
  210    w(i) = scale * g
         g = 0.0d0
         s = 0.0d0
         scale = 0.0d0
         if (i .gt. m .or. i .eq. n) go to 290
!
         do 220 k = l, n
         scale = scale + abs(u(i,k))
  220    continue
!
         if (abs(scale) .lt. thresh) go to 290
!
         do 230 k = l, n
            u(i,k) = u(i,k) / scale
            s = s + u(i,k)**2
  230    continue
!
         f = u(i,l)
         g = -sign(sqrt(s),f)
         h = f * g - s
         u(i,l) = f - g
!
         do 240 k = l, n
         rv1(k) = u(i,k) / h
  240    continue
!
         if (i .eq. m) go to 270
!
         do 260 j = l, m
            s = 0.0d0
!
            do 250 k = l, n
            s = s + u(j,k) * u(i,k)
  250       continue
!
            do 255 k = l, n
               u(j,k) = u(j,k) + s * rv1(k)
  255       continue
  260    continue
!
  270    do 280 k = l, n
         u(i,k) = scale * u(i,k)
  280    continue
!
  290    x = max(x,abs(w(i))+abs(rv1(i)))
  300 continue
!     .......... accumulation of right-hand transformations ..........
      if (.not. matv) go to 410
!     .......... for i=n step -1 until 1 do -- ..........
      do 400 ii = 1, n
         i = n + 1 - ii
         if (i .eq. n) go to 390
         if (abs(g) .lt. thresh) go to 360
!
         do 320 j = l, n
!     .......... double division avoids possible underflow ..........
         v(j,i) = (u(i,j) / u(i,l)) / g
  320    continue
!
         do 350 j = l, n
            s = 0.0d0
!
            do 340 k = l, n
            s = s + u(i,k) * v(k,j)
  340       continue
!
            do 355 k = l, n
               v(k,j) = v(k,j) + s * v(k,i)
  355       continue
  350    continue
!
  360    do 380 j = l, n
            v(i,j) = 0.0d0
            v(j,i) = 0.0d0
  380    continue
!
  390    v(i,i) = 1.0d0
         g = rv1(i)
         l = i
  400 continue
!     .......... accumulation of left-hand transformations ..........
  410 if (.not. matu) go to 510
!     ..........for i=min(m,n) step -1 until 1 do -- ..........
      mn = n
      if (m .lt. n) mn = m
!
      do 500 ii = 1, mn
         i = mn + 1 - ii
         l = i + 1
         g = w(i)
         if (i .eq. n) go to 430
!
         do 420 j = l, n
         u(i,j) = 0.0d0
  420    continue
!
  430    if (abs(g) .lt. thresh) go to 475
         if (i .eq. mn) go to 460
!
         do 450 j = l, n
            s = 0.0d0
!
            do 440 k = l, m
            s = s + u(k,i) * u(k,j)
  440       continue
!     .......... double division avoids possible underflow ..........
            f = (s / u(i,i)) / g
!
            do 455 k = i, m
               u(k,j) = u(k,j) + f * u(k,i)
  455       continue
  450    continue
!
  460    do 470 j = i, m
         u(j,i) = u(j,i) / g
  470    continue
!
         go to 490
!
  475    do 480 j = i, m
         u(j,i) = 0.0d0
  480    continue
!
  490    u(i,i) = u(i,i) + 1.0d0
  500 continue
!     .......... diagonalization of the bidiagonal form ..........
  510 tst1 = x
!     .......... for k=n step -1 until 1 do -- ..........
      do 700 kk = 1, n
         k1 = n - kk
         k = k1 + 1
         its = 0
!     .......... test for splitting.
!                for l=k step -1 until 1 do -- ..........
  520    do 530 ll = 1, k
            l1 = k - ll
            l = l1 + 1
            tst2 = tst1 + abs(rv1(l))
            if (abs(tst2-tst1) .lt. thresh) go to 565
!     .......... rv1(1) is always zero, so there is no exit
!                through the bottom of the loop ..........
            tst2 = tst1 + abs(w(l1))
            if (abs(tst2-tst1) .lt. thresh) go to 540
  530    continue
!     .......... cancellation of rv1(l) if l greater than 1 ..........
  540    c = 0.0d0
         s = 1.0d0
!
         do 560 i = l, k
            f = s * rv1(i)
            rv1(i) = c * rv1(i)
            tst2 = tst1 + abs(f)
            if (abs(tst2-tst1) .lt. thresh) go to 565
            g = w(i)
            h = casvb_pythag(f,g)
            w(i) = h
            c = g / h
            s = -f / h
            if (.not. matu) go to 560
!
            do 550 j = 1, m
               y = u(j,l1)
               z = u(j,i)
               u(j,l1) = y * c + z * s
               u(j,i) = -y * s + z * c
  550       continue
!
  560    continue
!     .......... test for convergence ..........
  565    z = w(k)
         if (l .eq. k) go to 650
!     .......... shift from bottom 2 by 2 minor ..........
         if (its .eq. 30) go to 1000
         its = its + 1
         x = w(l)
         y = w(k1)
         g = rv1(k1)
         h = rv1(k)
         f = 0.5d0 * (((g + z) / h) * ((g - z) / y) + y / h - h / y)
         g = casvb_pythag(f,1.0d0)
         f = x - (z / x) * z + (h / x) * (y / (f + sign(g,f)) - h)
!     .......... next qr transformation ..........
         c = 1.0d0
         s = 1.0d0
!
         do 600 i1 = l, k1
            i = i1 + 1
            g = rv1(i)
            y = w(i)
            h = s * g
            g = c * g
            z = casvb_pythag(f,h)
            rv1(i1) = z
            c = f / z
            s = h / z
            f = x * c + g * s
            g = -x * s + g * c
            h = y * s
            y = y * c
            if (.not. matv) go to 575
!
            do 570 j = 1, n
               x = v(j,i1)
               z = v(j,i)
               v(j,i1) = x * c + z * s
               v(j,i) = -x * s + z * c
  570       continue
!
  575       z = casvb_pythag(f,h)
            w(i1) = z
!     .......... rotation can be arbitrary if z is zero ..........
            if (abs(z) .lt. thresh) go to 580
            c = f / z
            s = h / z
  580       f = c * g + s * y
            x = -s * g + c * y
            if (.not. matu) go to 600
!
            do 590 j = 1, m
               y = u(j,i1)
               z = u(j,i)
               u(j,i1) = y * c + z * s
               u(j,i) = -y * s + z * c
  590       continue
!
  600    continue
!
         rv1(l) = 0.0d0
         rv1(k) = f
         w(k) = x
         go to 520
!     .......... convergence ..........
  650    if (z .ge. 0.0d0) go to 700
!     .......... w(k) is made non-negative ..........
         w(k) = -z
         if (.not. matv) go to 700
!
         do 690 j = 1, n
         v(j,k) = -v(j,k)
  690    continue
!
  700 continue
!
      go to 1001
!     .......... set error -- no convergence to a
!                singular value after 30 iterations ..........
 1000 ierr = k
 1001 return
      end

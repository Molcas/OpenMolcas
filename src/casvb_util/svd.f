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
      subroutine svd(nm,m,n,a,w,matu,u,matv,v,ierr,rv1)
c  TESTING FOR EQUALITY OR ZERO NOW DONE USING THRESHOLD, TT, 110100.
c
      integer i,j,k,l,m,n,ii,i1,kk,k1,ll,l1,mn,nm,its,ierr
      real*8 a(nm,n),w(n),u(nm,n),v(nm,n),rv1(n)
      real*8 c,f,g,h,s,x,y,z,tst1,tst2,scale,casvb_pythag,thresh
      logical matu,matv
      save thresh
      data thresh/1d-16/
c
c     this subroutine is a translation of the algol procedure svd,
c     num. math. 14, 403-420(1970) by golub and reinsch.
c     handbook for auto. comp., vol ii-linear algebra, 134-151(1971).
c
c     this subroutine determines the singular value decomposition
c          t
c     a=usv  of a real m by n rectangular matrix.  householder
c     bidiagonalization and a variant of the qr algorithm are used.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.  note that nm must be at least
c          as large as the maximum of m and n.
c
c        m is the number of rows of a (and u).
c
c        n is the number of columns of a (and u) and the order of v.
c
c        a contains the rectangular input matrix to be decomposed.
c
c        matu should be set to .true. if the u matrix in the
c          decomposition is desired, and to .false. otherwise.
c
c        matv should be set to .true. if the v matrix in the
c          decomposition is desired, and to .false. otherwise.
c
c     on output
c
c        a is unaltered (unless overwritten by u or v).
c
c        w contains the n (non-negative) singular values of a (the
c          diagonal elements of s).  they are unordered.  if an
c          error exit is made, the singular values should be correct
c          for indices ierr+1,ierr+2,...,n.
c
c        u contains the matrix u (orthogonal column vectors) of the
c          decomposition if matu has been set to .true.  otherwise
c          u is used as a temporary array.  u may coincide with a.
c          if an error exit is made, the columns of u corresponding
c          to indices of correct singular values should be correct.
c
c        v contains the matrix v (orthogonal) of the decomposition if
c          matv has been set to .true.  otherwise v is not referenced.
c          v may also coincide with a if u is not needed.  if an error
c          exit is made, the columns of v corresponding to indices of
c          correct singular values should be correct.
c
c        ierr is set to
c          zero       for normal return,
c          k          if the k-th singular value has not been
c                     determined after 30 iterations.
c
c        rv1 is a temporary storage array.
c
c     calls pythag for  sqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c

c The next two initializations are to appease a compiler
      f = 0.0d0
      h = 0.0d0

      ierr = 0
      l  = 0 ! dummy initialize
      l1 = 0 ! dummy initialize
c
      do 100 i = 1, m
c
         do 100 j = 1, n
            u(i,j) = a(i,j)
  100 continue
c     .......... householder reduction to bidiagonal form ..........
      g = 0.0d0
      scale = 0.0d0
      x = 0.0d0
c
      do 300 i = 1, n
         l = i + 1
         rv1(i) = scale * g
         g = 0.0d0
         s = 0.0d0
         scale = 0.0d0
         if (i .gt. m) go to 210
c
         do 120 k = i, m
  120    scale = scale + abs(u(k,i))
c
         if (abs(scale) .lt. thresh) go to 210
c
         do 130 k = i, m
            u(k,i) = u(k,i) / scale
            s = s + u(k,i)**2
  130    continue
c
         f = u(i,i)
         g = -sign(sqrt(s),f)
         h = f * g - s
         u(i,i) = f - g
         if (i .eq. n) go to 190
c
         do 150 j = l, n
            s = 0.0d0
c
            do 140 k = i, m
  140       s = s + u(k,i) * u(k,j)
c
            f = s / h
c
            do 150 k = i, m
               u(k,j) = u(k,j) + f * u(k,i)
  150    continue
c
  190    do 200 k = i, m
  200    u(k,i) = scale * u(k,i)
c
  210    w(i) = scale * g
         g = 0.0d0
         s = 0.0d0
         scale = 0.0d0
         if (i .gt. m .or. i .eq. n) go to 290
c
         do 220 k = l, n
  220    scale = scale + abs(u(i,k))
c
         if (abs(scale) .lt. thresh) go to 290
c
         do 230 k = l, n
            u(i,k) = u(i,k) / scale
            s = s + u(i,k)**2
  230    continue
c
         f = u(i,l)
         g = -sign(sqrt(s),f)
         h = f * g - s
         u(i,l) = f - g
c
         do 240 k = l, n
  240    rv1(k) = u(i,k) / h
c
         if (i .eq. m) go to 270
c
         do 260 j = l, m
            s = 0.0d0
c
            do 250 k = l, n
  250       s = s + u(j,k) * u(i,k)
c
            do 260 k = l, n
               u(j,k) = u(j,k) + s * rv1(k)
  260    continue
c
  270    do 280 k = l, n
  280    u(i,k) = scale * u(i,k)
c
  290    x = max(x,abs(w(i))+abs(rv1(i)))
  300 continue
c     .......... accumulation of right-hand transformations ..........
      if (.not. matv) go to 410
c     .......... for i=n step -1 until 1 do -- ..........
      do 400 ii = 1, n
         i = n + 1 - ii
         if (i .eq. n) go to 390
         if (abs(g) .lt. thresh) go to 360
c
         do 320 j = l, n
c     .......... double division avoids possible underflow ..........
  320    v(j,i) = (u(i,j) / u(i,l)) / g
c
         do 350 j = l, n
            s = 0.0d0
c
            do 340 k = l, n
  340       s = s + u(i,k) * v(k,j)
c
            do 350 k = l, n
               v(k,j) = v(k,j) + s * v(k,i)
  350    continue
c
  360    do 380 j = l, n
            v(i,j) = 0.0d0
            v(j,i) = 0.0d0
  380    continue
c
  390    v(i,i) = 1.0d0
         g = rv1(i)
         l = i
  400 continue
c     .......... accumulation of left-hand transformations ..........
  410 if (.not. matu) go to 510
c     ..........for i=min(m,n) step -1 until 1 do -- ..........
      mn = n
      if (m .lt. n) mn = m
c
      do 500 ii = 1, mn
         i = mn + 1 - ii
         l = i + 1
         g = w(i)
         if (i .eq. n) go to 430
c
         do 420 j = l, n
  420    u(i,j) = 0.0d0
c
  430    if (abs(g) .lt. thresh) go to 475
         if (i .eq. mn) go to 460
c
         do 450 j = l, n
            s = 0.0d0
c
            do 440 k = l, m
  440       s = s + u(k,i) * u(k,j)
c     .......... double division avoids possible underflow ..........
            f = (s / u(i,i)) / g
c
            do 450 k = i, m
               u(k,j) = u(k,j) + f * u(k,i)
  450    continue
c
  460    do 470 j = i, m
  470    u(j,i) = u(j,i) / g
c
         go to 490
c
  475    do 480 j = i, m
  480    u(j,i) = 0.0d0
c
  490    u(i,i) = u(i,i) + 1.0d0
  500 continue
c     .......... diagonalization of the bidiagonal form ..........
  510 tst1 = x
c     .......... for k=n step -1 until 1 do -- ..........
      do 700 kk = 1, n
         k1 = n - kk
         k = k1 + 1
         its = 0
c     .......... test for splitting.
c                for l=k step -1 until 1 do -- ..........
  520    do 530 ll = 1, k
            l1 = k - ll
            l = l1 + 1
            tst2 = tst1 + abs(rv1(l))
            if (abs(tst2-tst1) .lt. thresh) go to 565
c     .......... rv1(1) is always zero, so there is no exit
c                through the bottom of the loop ..........
            tst2 = tst1 + abs(w(l1))
            if (abs(tst2-tst1) .lt. thresh) go to 540
  530    continue
c     .......... cancellation of rv1(l) if l greater than 1 ..........
  540    c = 0.0d0
         s = 1.0d0
c
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
c
            do 550 j = 1, m
               y = u(j,l1)
               z = u(j,i)
               u(j,l1) = y * c + z * s
               u(j,i) = -y * s + z * c
  550       continue
c
  560    continue
c     .......... test for convergence ..........
  565    z = w(k)
         if (l .eq. k) go to 650
c     .......... shift from bottom 2 by 2 minor ..........
         if (its .eq. 30) go to 1000
         its = its + 1
         x = w(l)
         y = w(k1)
         g = rv1(k1)
         h = rv1(k)
         f = 0.5d0 * (((g + z) / h) * ((g - z) / y) + y / h - h / y)
         g = casvb_pythag(f,1.0d0)
         f = x - (z / x) * z + (h / x) * (y / (f + sign(g,f)) - h)
c     .......... next qr transformation ..........
         c = 1.0d0
         s = 1.0d0
c
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
c
            do 570 j = 1, n
               x = v(j,i1)
               z = v(j,i)
               v(j,i1) = x * c + z * s
               v(j,i) = -x * s + z * c
  570       continue
c
  575       z = casvb_pythag(f,h)
            w(i1) = z
c     .......... rotation can be arbitrary if z is zero ..........
            if (abs(z) .lt. thresh) go to 580
            c = f / z
            s = h / z
  580       f = c * g + s * y
            x = -s * g + c * y
            if (.not. matu) go to 600
c
            do 590 j = 1, m
               y = u(j,i1)
               z = u(j,i)
               u(j,i1) = y * c + z * s
               u(j,i) = -y * s + z * c
  590       continue
c
  600    continue
c
         rv1(l) = 0.0d0
         rv1(k) = f
         w(k) = x
         go to 520
c     .......... convergence ..........
  650    if (z .ge. 0.0d0) go to 700
c     .......... w(k) is made non-negative ..........
         w(k) = -z
         if (.not. matv) go to 700
c
         do 690 j = 1, n
  690    v(j,k) = -v(j,k)
c
  700 continue
c
      go to 1001
c     .......... set error -- no convergence to a
c                singular value after 30 iterations ..........
 1000 ierr = k
 1001 return
      end

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
      subroutine hqr2(nm,n,low,igh,h,wr,wi,z,ierr)
!
      integer i,j,k,l,m,n,en,ii,jj,ll,mm,na,nm,nn,                      &
     &        igh,itn,its,low,mp2,enm2,ierr
      real*8 h(nm,n),wr(n),wi(n),z(nm,n)
      real*8 p,q,r,s,t,w,x,y,ra,sa,vi,vr,zz,norm,tst1,tst2
      logical notlas
!
!     this subroutine is a translation of the algol procedure hqr2,
!     num. math. 16, 181-204(1970) by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a real upper hessenberg matrix by the qr method.  the
!     eigenvectors of a real general matrix can also be found
!     if  elmhes  and  eltran  or  orthes  and  ortran  have
!     been used to reduce this general matrix to hessenberg form
!     and to accumulate the similarity transformations.
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
!        h contains the upper hessenberg matrix.
!
!        z contains the transformation matrix produced by  eltran
!          after the reduction by  elmhes, or by  ortran  after the
!          reduction by  orthes, if performed.  if the eigenvectors
!          of the hessenberg matrix are desired, z must contain the
!          identity matrix.
!
!     on output
!
!        h has been destroyed.
!
!        wr and wi contain the real and imaginary parts,
!          respectively, of the eigenvalues.  the eigenvalues
!          are unordered except that complex conjugate pairs
!          of values appear consecutively with the eigenvalue
!          having the positive imaginary part first.  if an
!          error exit is made, the eigenvalues should be correct
!          for indices ierr+1,...,n.
!
!        z contains the real and imaginary parts of the eigenvectors.
!          if the i-th eigenvalue is real, the i-th column of z
!          contains its eigenvector.  if the i-th eigenvalue is complex
!          with positive imaginary part, the i-th and (i+1)-th
!          columns of z contain the real and imaginary parts of its
!          eigenvector.  the eigenvectors are unnormalized.  if an
!          error exit is made, none of the eigenvectors has been found.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the limit of 30*n iterations is exhausted
!                     while the j-th eigenvalue is being sought.
!
!     calls cdiv for complex division.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      w = 0.0d0
      zz = 0.0d0
      ierr = 0
      norm = 0.0d0
      k = 1
      l = 0     ! dummy initialize
      m = 0     ! dummy initialize
      p = 0.0D0 ! dummy initialize
      q = 0.0D0 ! dummy initialize
      r = 0.0D0 ! dummy initialize
      s = 0.0D0 ! dummy initialize
!     .......... store roots isolated by balanc
!                and compute matrix norm ..........
      do 50 i = 1, n
!
         do 40 j = k, n
         norm = norm + abs(h(i,j))
   40    continue
!
         k = i
         if (i .ge. low .and. i .le. igh) go to 50
         wr(i) = h(i,i)
         wi(i) = 0.0d0
   50 continue
!
      en = igh
      t = 0.0d0
      itn = 30*n
!     .......... search for next eigenvalues ..........
   60 if (en .lt. low) go to 340
      its = 0
      na = en - 1
      enm2 = na - 1
!     .......... look for single small sub-diagonal element
!                for l=en step -1 until low do -- ..........
   70 do 80 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 100
         s = abs(h(l-1,l-1)) + abs(h(l,l))
         if (s .eq. 0.0d0) s = norm
         tst1 = s
         tst2 = tst1 + abs(h(l,l-1))
         if (tst2 .eq. tst1) go to 100
   80 continue
!     .......... form shift ..........
  100 x = h(en,en)
      if (l .eq. en) go to 270
      y = h(na,na)
      w = h(en,na) * h(na,en)
      if (l .eq. na) go to 280
      if (itn .eq. 0) go to 1000
      if (its .ne. 10 .and. its .ne. 20) go to 130
!     .......... form exceptional shift ..........
      t = t + x
!
      do 120 i = low, en
      h(i,i) = h(i,i) - x
  120 continue
!
      s = abs(h(en,na)) + abs(h(na,enm2))
      x = 0.75d0 * s
      y = x
      w = -0.4375d0 * s * s
  130 its = its + 1
      itn = itn - 1
!     .......... look for two consecutive small
!                sub-diagonal elements.
!                for m=en-2 step -1 until l do -- ..........
      do 140 mm = l, enm2
         m = enm2 + l - mm
         zz = h(m,m)
         r = x - zz
         s = y - zz
         p = (r * s - w) / h(m+1,m) + h(m,m+1)
         q = h(m+1,m+1) - zz - r - s
         r = h(m+2,m+1)
         s = abs(p) + abs(q) + abs(r)
         p = p / s
         q = q / s
         r = r / s
         if (m .eq. l) go to 150
         tst1 = abs(p)*(abs(h(m-1,m-1)) + abs(zz) + abs(h(m+1,m+1)))
         tst2 = tst1 + abs(h(m,m-1))*(abs(q) + abs(r))
         if (tst2 .eq. tst1) go to 150
  140 continue
!
  150 mp2 = m + 2
!
      do 160 i = mp2, en
         h(i,i-2) = 0.0d0
         if (i .eq. mp2) go to 160
         h(i,i-3) = 0.0d0
  160 continue
!     .......... double qr step involving rows l to en and
!                columns m to en ..........
      do 260 k = m, na
         notlas = k .ne. na
         if (k .eq. m) go to 170
         p = h(k,k-1)
         q = h(k+1,k-1)
         r = 0.0d0
         if (notlas) r = h(k+2,k-1)
         x = abs(p) + abs(q) + abs(r)
         if (x .eq. 0.0d0) go to 260
         p = p / x
         q = q / x
         r = r / x
  170    s = sign(sqrt(p*p+q*q+r*r),p)
         if (k .eq. m) go to 180
         h(k,k-1) = -s * x
         go to 190
  180    if (l .ne. m) h(k,k-1) = -h(k,k-1)
  190    p = p + s
         x = p / s
         y = q / s
         zz = r / s
         q = q / p
         r = r / p
         if (notlas) go to 225
!     .......... row modification ..........
         do 200 j = k, n
            p = h(k,j) + q * h(k+1,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
  200    continue
!
         j = min(en,k+3)
!     .......... column modification ..........
         do 210 i = 1, j
            p = x * h(i,k) + y * h(i,k+1)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
  210    continue
!     .......... accumulate transformations ..........
         do 220 i = low, igh
            p = x * z(i,k) + y * z(i,k+1)
            z(i,k) = z(i,k) - p
            z(i,k+1) = z(i,k+1) - p * q
  220    continue
         go to 255
  225    continue
!     .......... row modification ..........
         do 230 j = k, n
            p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
            h(k+2,j) = h(k+2,j) - p * zz
  230    continue
!
         j = min(en,k+3)
!     .......... column modification ..........
         do 240 i = 1, j
            p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
            h(i,k+2) = h(i,k+2) - p * r
  240    continue
!     .......... accumulate transformations ..........
         do 250 i = low, igh
            p = x * z(i,k) + y * z(i,k+1) + zz * z(i,k+2)
            z(i,k) = z(i,k) - p
            z(i,k+1) = z(i,k+1) - p * q
            z(i,k+2) = z(i,k+2) - p * r
  250    continue
  255    continue
!
  260 continue
!
      go to 70
!     .......... one root found ..........
  270 h(en,en) = x + t
      wr(en) = h(en,en)
      wi(en) = 0.0d0
      en = na
      go to 60
!     .......... two roots found ..........
  280 p = (y - x) / 2.0d0
      q = p * p + w
      zz = sqrt(abs(q))
      h(en,en) = x + t
      x = h(en,en)
      h(na,na) = y + t
      if (q .lt. 0.0d0) go to 320
!     .......... real pair ..........
      zz = p + sign(zz,p)
      wr(na) = x + zz
      wr(en) = wr(na)
      if (zz .ne. 0.0d0) wr(en) = x - w / zz
      wi(na) = 0.0d0
      wi(en) = 0.0d0
      x = h(en,na)
      s = abs(x) + abs(zz)
      p = x / s
      q = zz / s
      r = sqrt(p*p+q*q)
      p = p / r
      q = q / r
!     .......... row modification ..........
      do 290 j = na, n
         zz = h(na,j)
         h(na,j) = q * zz + p * h(en,j)
         h(en,j) = q * h(en,j) - p * zz
  290 continue
!     .......... column modification ..........
      do 300 i = 1, en
         zz = h(i,na)
         h(i,na) = q * zz + p * h(i,en)
         h(i,en) = q * h(i,en) - p * zz
  300 continue
!     .......... accumulate transformations ..........
      do 310 i = low, igh
         zz = z(i,na)
         z(i,na) = q * zz + p * z(i,en)
         z(i,en) = q * z(i,en) - p * zz
  310 continue
!
      go to 330
!     .......... complex pair ..........
  320 wr(na) = x + p
      wr(en) = x + p
      wi(na) = zz
      wi(en) = -zz
  330 en = enm2
      go to 60
!     .......... all roots found.  backsubstitute to find
!                vectors of upper triangular form ..........
  340 if (norm .eq. 0.0d0) go to 1001
!     .......... for en=n step -1 until 1 do -- ..........
      do 800 nn = 1, n
         en = n + 1 - nn
         p = wr(en)
         q = wi(en)
         na = en - 1
!        if (q) 710, 600, 800
         if (q.lt.0.0d0) Goto 710
         if (q.eq.0.0d0) Goto 600
         if (q.gt.0.0d0) Goto 800
!     .......... real vector ..........
  600    m = en
         h(en,en) = 1.0d0
         if (na .eq. 0) go to 800
!     .......... for i=en-1 step -1 until 1 do -- ..........
         do 700 ii = 1, na
            i = en - ii
            w = h(i,i) - p
            r = 0.0d0
!
            do 610 j = m, en
            r = r + h(i,j) * h(j,en)
  610       continue
!
            if (wi(i) .ge. 0.0d0) go to 630
            zz = w
            s = r
            go to 700
  630       m = i
            if (wi(i) .ne. 0.0d0) go to 640
            t = w
            if (t .ne. 0.0d0) go to 635
               tst1 = norm
               t = tst1
  632          t = 0.01d0 * t
               tst2 = norm + t
               if (tst2 .gt. tst1) go to 632
  635       h(i,en) = -r / t
            go to 680
!     .......... solve real equations ..........
  640       x = h(i,i+1)
            y = h(i+1,i)
            q = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i)
            t = (x * s - zz * r) / q
            h(i,en) = t
            if (abs(x) .le. abs(zz)) go to 650
            h(i+1,en) = (-r - w * t) / x
            go to 680
  650       h(i+1,en) = (-s - y * t) / zz
!
!     .......... overflow control ..........
  680       t = abs(h(i,en))
            if (t .eq. 0.0d0) go to 700
            tst1 = t
            tst2 = tst1 + 1.0d0/tst1
            if (tst2 .gt. tst1) go to 700
            do 690 j = i, en
               h(j,en) = h(j,en)/t
  690       continue
!
  700    continue
!     .......... end real vector ..........
         go to 800
!     .......... complex vector ..........
  710    m = na
!     .......... last vector component chosen imaginary so that
!                eigenvector matrix is triangular ..........
         if (abs(h(en,na)) .le. abs(h(na,en))) go to 720
         h(na,na) = q / h(en,na)
         h(na,en) = -(h(en,en) - p) / h(en,na)
         go to 730
  720    call cdiv(0.0d0,-h(na,en),h(na,na)-p,q,h(na,na),h(na,en))
  730    h(en,na) = 0.0d0
         h(en,en) = 1.0d0
         enm2 = na - 1
         if (enm2 .eq. 0) go to 800
!     .......... for i=en-2 step -1 until 1 do -- ..........
         do 795 ii = 1, enm2
            i = na - ii
            w = h(i,i) - p
            ra = 0.0d0
            sa = 0.0d0
!
            do 760 j = m, en
               ra = ra + h(i,j) * h(j,na)
               sa = sa + h(i,j) * h(j,en)
  760       continue
!
            if (wi(i) .ge. 0.0d0) go to 770
            zz = w
            r = ra
            s = sa
            go to 795
  770       m = i
            if (wi(i) .ne. 0.0d0) go to 780
            call cdiv(-ra,-sa,w,q,h(i,na),h(i,en))
            go to 790
!     .......... solve complex equations ..........
  780       x = h(i,i+1)
            y = h(i+1,i)
            vr = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i) - q * q
            vi = (wr(i) - p) * 2.0d0 * q
            if (vr .ne. 0.0d0 .or. vi .ne. 0.0d0) go to 784
               tst1 = norm * (abs(w) + abs(q) + abs(x)                  &
     &                      + abs(y) + abs(zz))
               vr = tst1
  783          vr = 0.01d0 * vr
               tst2 = tst1 + vr
               if (tst2 .gt. tst1) go to 783
  784       call cdiv(x*r-zz*ra+q*sa,x*s-zz*sa-q*ra,vr,vi,              &
     &                h(i,na),h(i,en))
            if (abs(x) .le. abs(zz) + abs(q)) go to 785
            h(i+1,na) = (-ra - w * h(i,na) + q * h(i,en)) / x
            h(i+1,en) = (-sa - w * h(i,en) - q * h(i,na)) / x
            go to 790
  785       call cdiv(-r-y*h(i,na),-s-y*h(i,en),zz,q,                   &
     &                h(i+1,na),h(i+1,en))
!
!     .......... overflow control ..........
  790       t = max(abs(h(i,na)), abs(h(i,en)))
            if (t .eq. 0.0d0) go to 795
            tst1 = t
            tst2 = tst1 + 1.0d0/tst1
            if (tst2 .gt. tst1) go to 795
            do 792 j = i, en
               h(j,na) = h(j,na)/t
               h(j,en) = h(j,en)/t
  792       continue
!
  795    continue
!     .......... end complex vector ..........
  800 continue
!     .......... end back substitution.
!                vectors of isolated roots ..........
      do 840 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 840
!
         do 820 j = i, n
         z(i,j) = h(i,j)
  820    continue
!
  840 continue
!     .......... multiply by transformation matrix to give
!                vectors of original full matrix.
!                for j=n step -1 until low do -- ..........
      do 880 jj = low, n
         j = n + low - jj
         m = min(j,igh)
!
         do 881 i = low, igh
            zz = 0.0d0
!
            do 860 k = low, m
            zz = zz + z(i,k) * h(k,j)
  860       continue
!
            z(i,j) = zz
  881    continue
  880 continue
!
      go to 1001
!     .......... set error -- all eigenvalues have not
!                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end

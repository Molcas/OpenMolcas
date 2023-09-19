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
!********************************************************
!** Public-domain library routines used by casvb only. **
!********************************************************
!**********************
!** EISPACK ROUTINES **
!**********************

subroutine hqr2(nm,n,low,igh,h,wr,wi,z,ierr)
! this subroutine is a translation of the algol procedure hqr2,
! num. math. 16, 181-204(1970) by peters and wilkinson.
! handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
!
! this subroutine finds the eigenvalues and eigenvectors
! of a real upper hessenberg matrix by the qr method.  the
! eigenvectors of a real general matrix can also be found
! if  elmhes  and  eltran  or  orthes  and  ortran  have
! been used to reduce this general matrix to hessenberg form
! and to accumulate the similarity transformations.
!
! on input
!
!    nm must be set to the row dimension of two-dimensional
!      array parameters as declared in the calling program
!      dimension statement.
!
!    n is the order of the matrix.
!
!    low and igh are integers determined by the balancing
!      subroutine  balanc.  if  balanc  has not been used,
!      set low=1, igh=n.
!
!    h contains the upper hessenberg matrix.
!
!    z contains the transformation matrix produced by  eltran
!      after the reduction by  elmhes, or by  ortran  after the
!      reduction by  orthes, if performed.  if the eigenvectors
!      of the hessenberg matrix are desired, z must contain the
!      identity matrix.
!
! on output
!
!    h has been destroyed.
!
!    wr and wi contain the real and imaginary parts,
!      respectively, of the eigenvalues.  the eigenvalues
!      are unordered except that complex conjugate pairs
!      of values appear consecutively with the eigenvalue
!      having the positive imaginary part first.  if an
!      error exit is made, the eigenvalues should be correct
!      for indices ierr+1,...,n.
!
!    z contains the real and imaginary parts of the eigenvectors.
!      if the i-th eigenvalue is real, the i-th column of z
!      contains its eigenvector.  if the i-th eigenvalue is complex
!      with positive imaginary part, the i-th and (i+1)-th
!      columns of z contain the real and imaginary parts of its
!      eigenvector.  the eigenvectors are unnormalized.  if an
!      error exit is made, none of the eigenvectors has been found.
!
!    ierr is set to
!      zero       for normal return,
!      j          if the limit of 30*n iterations is exhausted
!                 while the j-th eigenvalue is being sought.
!
! calls cdiv for complex division.
!
! questions and comments should be directed to burton s. garbow,
! mathematics and computer science div, argonne national laboratory
!
! this version dated august 1983.
!
! Updated to Fortran 90+ (Sep. 2023)
! ----------------------------------------------------------------------

integer i, j, k, l, m, n, en, ii, jj, ll, mm, na, nm, nn, igh, itn, its, low, mp2, enm2, ierr, nroot
real*8 h(nm,n), wr(n), wi(n), z(nm,n)
real*8 p, q, r, s, t, w, x, y, ra, sa, vi, vr, zz, norm, tst1, tst2
logical notlas

w = 0.0d0
zz = 0.0d0
ierr = 0
norm = 0.0d0
k = 1
l = 0     ! dummy initialize
m = 0     ! dummy initialize
p = 0.0d0 ! dummy initialize
q = 0.0d0 ! dummy initialize
r = 0.0d0 ! dummy initialize
s = 0.0d0 ! dummy initialize
nroot = 0 ! dummy initialize
! .......... store roots isolated by balanc and compute matrix norm ..........
do i=1,n

  do j=k,n
    norm = norm+abs(h(i,j))
  end do

  k = i
  if ((i >= low) .and. (i <= igh)) cycle
  wr(i) = h(i,i)
  wi(i) = 0.0d0
end do

en = igh
t = 0.0d0
itn = 30*n
do
  ! .......... search for next eigenvalues ..........
  if (en < low) exit
  its = 0
  na = en-1
  enm2 = na-1
  do
    ! .......... look for single small sub-diagonal element for l=en step -1 until low do -- ..........
    do ll=low,en
      l = en+low-ll
      if (l == low) exit
      s = abs(h(l-1,l-1))+abs(h(l,l))
      if (s == 0.0d0) s = norm
      tst1 = s
      tst2 = tst1+abs(h(l,l-1))
      if (tst2 == tst1) exit
    end do
    ! .......... form shift ..........
    x = h(en,en)
    if (l == en) then
      nroot = 1
      exit
    end if
    y = h(na,na)
    w = h(en,na)*h(na,en)
    if (l == na) then
      nroot = 2
      exit
    end if
    if (itn == 0) then
      ! .......... set error -- all eigenvalues have not converged after 30*n iterations ..........
      ierr = en
      return
    end if
    if ((its == 10) .or. (its == 20)) then
      ! .......... form exceptional shift ..........
      t = t+x

      do i=low,en
        h(i,i) = h(i,i)-x
      end do

      s = abs(h(en,na))+abs(h(na,enm2))
      x = 0.75d0*s
      y = x
      w = -0.4375d0*s*s
    end if
    its = its+1
    itn = itn-1
    ! .......... look for two consecutive small sub-diagonal elements.
    !            for m=en-2 step -1 until l do -- ..........
    do mm=l,enm2
      m = enm2+l-mm
      zz = h(m,m)
      r = x-zz
      s = y-zz
      p = (r*s-w)/h(m+1,m)+h(m,m+1)
      q = h(m+1,m+1)-zz-r-s
      r = h(m+2,m+1)
      s = abs(p)+abs(q)+abs(r)
      p = p/s
      q = q/s
      r = r/s
      if (m == l) exit
      tst1 = abs(p)*(abs(h(m-1,m-1))+abs(zz)+abs(h(m+1,m+1)))
      tst2 = tst1+abs(h(m,m-1))*(abs(q)+abs(r))
      if (tst2 == tst1) exit
    end do

    mp2 = m+2

    do i=mp2,en
      h(i,i-2) = 0.0d0
      if (i == mp2) cycle
      h(i,i-3) = 0.0d0
    end do
    ! .......... double qr step involving rows l to en and columns m to en ..........
    do k=m,na
      notlas = k /= na
      if (k /= m) then
        p = h(k,k-1)
        q = h(k+1,k-1)
        r = 0.0d0
        if (notlas) r = h(k+2,k-1)
        x = abs(p)+abs(q)+abs(r)
        if (x == 0.0d0) cycle
        p = p/x
        q = q/x
        r = r/x
      end if
      s = sign(sqrt(p*p+q*q+r*r),p)
      if (k == m) then
        if (l /= m) h(k,k-1) = -h(k,k-1)
      else
        h(k,k-1) = -s*x
      end if
      p = p+s
      x = p/s
      y = q/s
      zz = r/s
      q = q/p
      r = r/p
      if (notlas) then
        ! .......... row modification ..........
        do j=k,n
          p = h(k,j)+q*h(k+1,j)+r*h(k+2,j)
          h(k,j) = h(k,j)-p*x
          h(k+1,j) = h(k+1,j)-p*y
          h(k+2,j) = h(k+2,j)-p*zz
        end do

        j = min(en,k+3)
        ! .......... column modification ..........
        do i=1,j
          p = x*h(i,k)+y*h(i,k+1)+zz*h(i,k+2)
          h(i,k) = h(i,k)-p
          h(i,k+1) = h(i,k+1)-p*q
          h(i,k+2) = h(i,k+2)-p*r
        end do
        ! .......... accumulate transformations ..........
        do i=low,igh
          p = x*z(i,k)+y*z(i,k+1)+zz*z(i,k+2)
          z(i,k) = z(i,k)-p
          z(i,k+1) = z(i,k+1)-p*q
          z(i,k+2) = z(i,k+2)-p*r
        end do
      else
        ! .......... row modification ..........
        do j=k,n
          p = h(k,j)+q*h(k+1,j)
          h(k,j) = h(k,j)-p*x
          h(k+1,j) = h(k+1,j)-p*y
        end do

        j = min(en,k+3)
        ! .......... column modification ..........
        do i=1,j
          p = x*h(i,k)+y*h(i,k+1)
          h(i,k) = h(i,k)-p
          h(i,k+1) = h(i,k+1)-p*q
        end do
        ! .......... accumulate transformations ..........
        do i=low,igh
          p = x*z(i,k)+y*z(i,k+1)
          z(i,k) = z(i,k)-p
          z(i,k+1) = z(i,k+1)-p*q
        end do
      end if

    end do
  end do

  if (nroot == 1) then
    ! .......... one root found ..........
    h(en,en) = x+t
    wr(en) = h(en,en)
    wi(en) = 0.0d0
    en = na
  else if (nroot == 2) then
    ! .......... two roots found ..........
    p = (y-x)/2.0d0
    q = p*p+w
    zz = sqrt(abs(q))
    h(en,en) = x+t
    x = h(en,en)
    h(na,na) = y+t
    if (q < 0.0d0) then
      ! .......... complex pair ..........
      wr(na) = x+p
      wr(en) = x+p
      wi(na) = zz
      wi(en) = -zz
    else
      ! .......... real pair ..........
      zz = p+sign(zz,p)
      wr(na) = x+zz
      wr(en) = wr(na)
      if (zz /= 0.0d0) wr(en) = x-w/zz
      wi(na) = 0.0d0
      wi(en) = 0.0d0
      x = h(en,na)
      s = abs(x)+abs(zz)
      p = x/s
      q = zz/s
      r = sqrt(p*p+q*q)
      p = p/r
      q = q/r
      ! .......... row modification ..........
      do j=na,n
        zz = h(na,j)
        h(na,j) = q*zz+p*h(en,j)
        h(en,j) = q*h(en,j)-p*zz
      end do
      ! .......... column modification ..........
      do i=1,en
        zz = h(i,na)
        h(i,na) = q*zz+p*h(i,en)
        h(i,en) = q*h(i,en)-p*zz
      end do
      ! .......... accumulate transformations ..........
      do i=low,igh
        zz = z(i,na)
        z(i,na) = q*zz+p*z(i,en)
        z(i,en) = q*z(i,en)-p*zz
      end do
    end if

    en = enm2
  end if
end do
! .......... all roots found.  backsubstitute to find vectors of upper triangular form ..........
if (norm == 0.0d0) return
! .......... for en=n step -1 until 1 do -- ..........
do nn=1,n
  en = n+1-nn
  p = wr(en)
  q = wi(en)
  na = en-1
  if (q == 0.0d0) then
    ! .......... real vector ..........
    m = en
    h(en,en) = 1.0d0
    if (na == 0) cycle
    ! .......... for i=en-1 step -1 until 1 do -- ..........
    do ii=1,na
      i = en-ii
      w = h(i,i)-p
      r = 0.0d0

      do j=m,en
        r = r+h(i,j)*h(j,en)
      end do

      if (wi(i) < 0.0d0) then
        zz = w
        s = r
        cycle
      end if
      m = i
      if (wi(i) /= 0.0d0) then
        ! .......... solve real equations ..........
        x = h(i,i+1)
        y = h(i+1,i)
        q = (wr(i)-p)*(wr(i)-p)+wi(i)*wi(i)
        t = (x*s-zz*r)/q
        h(i,en) = t
        if (abs(x) <= abs(zz)) then
          h(i+1,en) = (-s-y*t)/zz
        else
          h(i+1,en) = (-r-w*t)/x
        end if
      else
        t = w
        if (t == 0.0d0) then
          tst1 = norm
          t = tst1
          do
            t = 0.01d0*t
            tst2 = norm+t
            if (tst2 <= tst1) exit
          end do
        end if
        h(i,en) = -r/t
      end if

      ! .......... overflow control ..........
      t = abs(h(i,en))
      if (t == 0.0d0) cycle
      tst1 = t
      tst2 = tst1+1.0d0/tst1
      if (tst2 > tst1) cycle
      do j=i,en
        h(j,en) = h(j,en)/t
      end do

    end do
    ! .......... end real vector ..........
  else if (q < 0.0d0) then
    ! .......... complex vector ..........
    m = na
    ! .......... last vector component chosen imaginary so that eigenvector matrix is triangular ..........
    if (abs(h(en,na)) <= abs(h(na,en))) then
      call cdiv(0.0d0,-h(na,en),h(na,na)-p,q,h(na,na),h(na,en))
    else
      h(na,na) = q/h(en,na)
      h(na,en) = -(h(en,en)-p)/h(en,na)
    end if
    h(en,na) = 0.0d0
    h(en,en) = 1.0d0
    enm2 = na-1
    if (enm2 == 0) cycle
    ! .......... for i=en-2 step -1 until 1 do -- ..........
    do ii=1,enm2
      i = na-ii
      w = h(i,i)-p
      ra = 0.0d0
      sa = 0.0d0

      do j=m,en
        ra = ra+h(i,j)*h(j,na)
        sa = sa+h(i,j)*h(j,en)
      end do

      if (wi(i) < 0.0d0) then
        zz = w
        r = ra
        s = sa
      else
        m = i
        if (wi(i) /= 0.0d0) then
          ! .......... solve complex equations ..........
          x = h(i,i+1)
          y = h(i+1,i)
          vr = (wr(i)-p)*(wr(i)-p)+wi(i)*wi(i)-q*q
          vi = (wr(i)-p)*2.0d0*q
          if ((vr == 0.0d0) .and. (vi == 0.0d0)) then
            tst1 = norm*(abs(w)+abs(q)+abs(x)+abs(y)+abs(zz))
            vr = tst1
            do
              vr = 0.01d0*vr
              tst2 = tst1+vr
              if (tst2 <= tst1) exit
            end do
          end if
          call cdiv(x*r-zz*ra+q*sa,x*s-zz*sa-q*ra,vr,vi,h(i,na),h(i,en))
          if (abs(x) <= abs(zz)+abs(q)) then
            call cdiv(-r-y*h(i,na),-s-y*h(i,en),zz,q,h(i+1,na),h(i+1,en))
          else
            h(i+1,na) = (-ra-w*h(i,na)+q*h(i,en))/x
            h(i+1,en) = (-sa-w*h(i,en)-q*h(i,na))/x
          end if
        else
          call cdiv(-ra,-sa,w,q,h(i,na),h(i,en))
        end if

        ! .......... overflow control ..........
        t = max(abs(h(i,na)),abs(h(i,en)))
        if (t == 0.0d0) cycle
        tst1 = t
        tst2 = tst1+1.0d0/tst1
        if (tst2 > tst1) cycle
        do j=i,en
          h(j,na) = h(j,na)/t
          h(j,en) = h(j,en)/t
        end do
      end if

    end do
    ! .......... end complex vector ..........
  end if
end do
! .......... end back substitution. vectors of isolated roots ..........
do i=1,n
  if ((i >= low) .and. (i <= igh)) cycle

  do j=i,n
    z(i,j) = h(i,j)
  end do

end do
! .......... multiply by transformation matrix to give vectors of original full matrix.
!            for j=n step -1 until low do -- ..........
do jj=low,n
  j = n+low-jj
  m = min(j,igh)

  do i=low,igh
    zz = 0.0d0

    do k=low,m
      zz = zz+z(i,k)*h(k,j)
    end do

    z(i,j) = zz
  end do
end do

return

end subroutine hqr2

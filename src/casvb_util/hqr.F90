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

subroutine hqr(nm,n,low,igh,h,wr,wi,ierr)
! this subroutine is a translation of the algol procedure hqr,
! num. math. 14, 219-231(1970) by martin, peters, and wilkinson.
! handbook for auto. comp., vol.ii-linear algebra, 359-371(1971).
!
! this subroutine finds the eigenvalues of a real
! upper hessenberg matrix by the qr method.
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
!    h contains the upper hessenberg matrix.  information about
!      the transformations used in the reduction to hessenberg
!      form by  elmhes  or  orthes, if performed, is stored
!      in the remaining triangle under the hessenberg matrix.
!
! on output
!
!    h has been destroyed.  therefore, it must be saved
!      before calling  hqr  if subsequent calculation and
!      back transformation of eigenvectors is to be performed.
!
!    wr and wi contain the real and imaginary parts,
!      respectively, of the eigenvalues.  the eigenvalues
!      are unordered except that complex conjugate pairs
!      of values appear consecutively with the eigenvalue
!      having the positive imaginary part first.  if an
!      error exit is made, the eigenvalues should be correct
!      for indices ierr+1,...,n.
!
!    ierr is set to
!      zero       for normal return,
!      j          if the limit of 30*n iterations is exhausted
!                 while the j-th eigenvalue is being sought.
!
! questions and comments should be directed to burton s. garbow,
! mathematics and computer science div, argonne national laboratory
!
! this version dated september 1989.
!
! ----------------------------------------------------------------------
!
!  RESTORED CORRECT INDICES OF LOOPS (200,210,230,240). (9/29/89 BSG)

integer i, j, k, l, m, n, en, ll, mm, na, nm, igh, itn, its, low, mp2, enm2, ierr
real*8 h(nm,n), wr(n), wi(n)
real*8 p, q, r, s, t, w, x, y, zz, norm, tst1, tst2
logical notlas

ierr = 0
norm = 0.0d0
k = 1
l = 0     ! dummy initialize
m = 0     ! dummy initialize
p = 0.0d0 ! dummy initialize
q = 0.0d0 ! dummy initialize
r = 0.0d0 ! dummy initialize
! .......... store roots isolated by balanc and compute matrix norm ..........
do i=1,n

  do j=k,n
    norm = norm+abs(h(i,j))
  end do

  k = i
  if ((i >= low) .and. (i <= igh)) go to 50
  wr(i) = h(i,i)
  wi(i) = 0.0d0
50 continue
end do

en = igh
t = 0.0d0
itn = 30*n
! .......... search for next eigenvalues ..........
60 if (en < low) go to 1001
its = 0
na = en-1
enm2 = na-1
! .......... look for single small sub-diagonal element for l=en step -1 until low do -- ..........
70 do ll=low,en
  l = en+low-ll
  if (l == low) go to 100
  s = abs(h(l-1,l-1))+abs(h(l,l))
  if (s == 0.0d0) s = norm
  tst1 = s
  tst2 = tst1+abs(h(l,l-1))
  if (tst2 == tst1) go to 100
end do
! .......... form shift ..........
100 x = h(en,en)
if (l == en) go to 270
y = h(na,na)
w = h(en,na)*h(na,en)
if (l == na) go to 280
if (itn == 0) go to 1000
if ((its /= 10) .and. (its /= 20)) go to 130
! .......... form exceptional shift ..........
t = t+x

do i=low,en
  h(i,i) = h(i,i)-x
end do
!
s = abs(h(en,na))+abs(h(na,enm2))
x = 0.75d0*s
y = x
w = -0.4375d0*s*s
130 its = its+1
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
  if (m == l) go to 150
  tst1 = abs(p)*(abs(h(m-1,m-1))+abs(zz)+abs(h(m+1,m+1)))
  tst2 = tst1+abs(h(m,m-1))*(abs(q)+abs(r))
  if (tst2 == tst1) go to 150
end do

150 mp2 = m+2

do i=mp2,en
  h(i,i-2) = 0.0d0
  if (i == mp2) go to 160
  h(i,i-3) = 0.0d0
160 continue
end do
! .......... double qr step involving rows l to en and columns m to en ..........
do k=m,na
  notlas = k /= na
  if (k == m) go to 170
  p = h(k,k-1)
  q = h(k+1,k-1)
  r = 0.0d0
  if (notlas) r = h(k+2,k-1)
  x = abs(p)+abs(q)+abs(r)
  if (x == 0.0d0) go to 260
  p = p/x
  q = q/x
  r = r/x
170 s = sign(sqrt(p*p+q*q+r*r),p)
  if (k == m) go to 180
  h(k,k-1) = -s*x
  go to 190
180 if (l /= m) h(k,k-1) = -h(k,k-1)
190 p = p+s
  x = p/s
  y = q/s
  zz = r/s
  q = q/p
  r = r/p
  if (notlas) go to 225
  ! .......... row modification ..........
  do j=k,EN
    p = h(k,j)+q*h(k+1,j)
    h(k,j) = h(k,j)-p*x
    h(k+1,j) = h(k+1,j)-p*y
  end do

  j = min(en,k+3)
  ! .......... column modification ..........
  do i=L,j
    p = x*h(i,k)+y*h(i,k+1)
    h(i,k) = h(i,k)-p
    h(i,k+1) = h(i,k+1)-p*q
  end do
  go to 255
225 continue
  ! .......... row modification ..........
  do j=k,EN
    p = h(k,j)+q*h(k+1,j)+r*h(k+2,j)
    h(k,j) = h(k,j)-p*x
    h(k+1,j) = h(k+1,j)-p*y
    h(k+2,j) = h(k+2,j)-p*zz
  end do

  j = min(en,k+3)
  ! .......... column modification ..........
  do i=L,j
    p = x*h(i,k)+y*h(i,k+1)+zz*h(i,k+2)
    h(i,k) = h(i,k)-p
    h(i,k+1) = h(i,k+1)-p*q
    h(i,k+2) = h(i,k+2)-p*r
  end do
255 continue

260 continue
end do

go to 70
! .......... one root found ..........
270 wr(en) = x+t
wi(en) = 0.0d0
en = na
go to 60
! .......... two roots found ..........
280 p = (y-x)/2.0d0
q = p*p+w
zz = sqrt(abs(q))
x = x+t
if (q < 0.0d0) go to 320
! .......... real pair ..........
zz = p+sign(zz,p)
wr(na) = x+zz
wr(en) = wr(na)
if (zz /= 0.0d0) wr(en) = x-w/zz
wi(na) = 0.0d0
wi(en) = 0.0d0
go to 330
! .......... complex pair ..........
320 wr(na) = x+p
wr(en) = x+p
wi(na) = zz
wi(en) = -zz
330 en = enm2
go to 60
! .......... set error -- all eigenvalues have not converged after 30*n iterations ..........
1000 ierr = en

1001 return

end subroutine hqr

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
! Updated to Fortran 90+ (Sep. 2023)
! ----------------------------------------------------------------------
!
!  RESTORED CORRECT INDICES OF LOOPS (200,210,230,240). (9/29/89 BSG)

use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nm, n, low, igh
real(kind=wp), intent(inout) :: h(nm,n)
real(kind=wp), intent(out) :: wr(n), wi(n)
integer(kind=iwp), intent(out) :: ierr
integer(kind=iwp) :: en, enm2, i, itn, its, j, k, l, ll, m, mm, mp2, na, nroot
real(kind=wp) :: norm, p, q, r, s, t, tst1, tst2, w, x, y, zz
logical(kind=iwp) :: notlas

ierr = 0
norm = Zero
k = 1
l = 0     ! dummy initialize
m = 0     ! dummy initialize
p = Zero  ! dummy initialize
q = Zero  ! dummy initialize
r = Zero  ! dummy initialize
nroot = 0 ! dummy initialize
! .......... store roots isolated by balanc and compute matrix norm ..........
do i=1,n

  do j=k,n
    norm = norm+abs(h(i,j))
  end do

  k = i
  if ((i >= low) .and. (i <= igh)) cycle
  wr(i) = h(i,i)
  wi(i) = Zero
end do

en = igh
t = Zero
itn = 30*n
do
  ! .......... search for next eigenvalues ..........
  if (en < low) return
  its = 0
  na = en-1
  enm2 = na-1
  do
    ! .......... look for single small sub-diagonal element for l=en step -1 until low do -- ..........
    do ll=low,en
      l = en+low-ll
      if (l == low) exit
      s = abs(h(l-1,l-1))+abs(h(l,l))
      if (s == Zero) s = norm
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
      x = 0.75_wp*s
      y = x
      w = -0.4375_wp*s*s
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
      h(i,i-2) = Zero
      if (i == mp2) cycle
      h(i,i-3) = Zero
    end do
    ! .......... double qr step involving rows l to en and columns m to en ..........
    do k=m,na
      notlas = k /= na
      if (k /= m) then
        p = h(k,k-1)
        q = h(k+1,k-1)
        r = Zero
        if (notlas) r = h(k+2,k-1)
        x = abs(p)+abs(q)+abs(r)
        if (x == Zero) cycle
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
      else
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
      end if

    end do
  end do

  if (nroot == 1) then
    ! .......... one root found ..........
    wr(en) = x+t
    wi(en) = Zero
    en = na
  else if (nroot == 2) then
    ! .......... two roots found ..........
    p = (y-x)*Half
    q = p*p+w
    zz = sqrt(abs(q))
    x = x+t
    if (q < Zero) then
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
      if (zz /= Zero) wr(en) = x-w/zz
      wi(na) = Zero
      wi(en) = Zero
    end if
    en = enm2
  end if
end do

return

end subroutine hqr

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

subroutine svd(nm,m,n,a,w,matu,u,matv,v,ierr,rv1)
! this subroutine is a translation of the algol procedure svd,
! num. math. 14, 403-420(1970) by golub and reinsch.
! handbook for auto. comp., vol ii-linear algebra, 134-151(1971).
!
! this subroutine determines the singular value decomposition
!      t
! a=usv  of a real m by n rectangular matrix.  householder
! bidiagonalization and a variant of the qr algorithm are used.
!
! on input
!
!    nm must be set to the row dimension of two-dimensional
!      array parameters as declared in the calling program
!      dimension statement.  note that nm must be at least
!      as large as the maximum of m and n.
!
!    m is the number of rows of a (and u).
!
!    n is the number of columns of a (and u) and the order of v.
!
!    a contains the rectangular input matrix to be decomposed.
!
!    matu should be set to .true. if the u matrix in the
!      decomposition is desired, and to .false. otherwise.
!
!    matv should be set to .true. if the v matrix in the
!      decomposition is desired, and to .false. otherwise.
!
! on output
!
!    a is unaltered (unless overwritten by u or v).
!
!    w contains the n (non-negative) singular values of a (the
!      diagonal elements of s).  they are unordered.  if an
!      error exit is made, the singular values should be correct
!      for indices ierr+1,ierr+2,...,n.
!
!    u contains the matrix u (orthogonal column vectors) of the
!      decomposition if matu has been set to .true.  otherwise
!      u is used as a temporary array.  u may coincide with a.
!      if an error exit is made, the columns of u corresponding
!      to indices of correct singular values should be correct.
!
!    v contains the matrix v (orthogonal) of the decomposition if
!      matv has been set to .true.  otherwise v is not referenced.
!      v may also coincide with a if u is not needed.  if an error
!      exit is made, the columns of v corresponding to indices of
!      correct singular values should be correct.
!
!    ierr is set to
!      zero       for normal return,
!      k          if the k-th singular value has not been
!                 determined after 30 iterations.
!
!    rv1 is a temporary storage array.
!
! calls pythag for  sqrt(a*a + b*b) .
!
! questions and comments should be directed to burton s. garbow,
! mathematics and computer science div, argonne national laboratory
!
! this version dated august 1983.
!
! Updated to Fortran 90+ (Sep. 2023)
! ----------------------------------------------------------------------
!
!  TESTING FOR EQUALITY OR ZERO NOW DONE USING THRESHOLD, TT, 110100.

use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nm, m, n
real(kind=wp), intent(in) :: a(nm,n)
real(kind=wp), intent(out) :: w(n), u(nm,n), rv1(n)
real(kind=wp), intent(inout) :: v(nm,n)
logical(kind=iwp), intent(in) :: matu, matv
integer(kind=iwp), intent(out) :: ierr
integer(kind=iwp) :: i, i1, ii, its, j, k, k1, kk, l, l1, ll, mn
real(kind=wp) :: c, f, g, h, pythag, s, scl, tst1, tst2, x, y, z
logical(kind=iwp) :: skip
real(kind=wp), parameter :: thresh = 1.0e-16_wp

! The next two initializations are to appease a compiler
f = Zero
h = Zero

ierr = 0
l = 0 ! dummy initialize
l1 = 0 ! dummy initialize

u(1:m,1:n) = a(1:m,1:n)
! .......... householder reduction to bidiagonal form ..........
g = Zero
scl = Zero
x = Zero

do i=1,n
  l = i+1
  rv1(i) = scl*g
  g = Zero
  s = Zero
  scl = Zero

  if (i <= m) then
    scl = scl+sum(abs(u(i:m,i)))

    if (abs(scl) >= thresh) then
      u(i:m,i) = u(i:m,i)/scl
      s = s+sum(u(i:m,i)**2)

      f = u(i,i)
      g = -sign(sqrt(s),f)
      h = f*g-s
      u(i,i) = f-g

      if (i /= n) then
        do j=l,n

          s = sum(u(i:m,i)*u(i:m,j))

          f = s/h

          u(i:m,j) = u(i:m,j)+f*u(i:m,i)
        end do
      end if

      u(i:m,i) = scl*u(i:m,i)
    end if
  end if

  w(i) = scl*g
  g = Zero
  s = Zero
  scl = Zero

  if ((i <= m) .and. (i /= n)) then
    scl = scl+sum(abs(u(i,l:)))

    if (abs(scl) >= thresh) then
      u(i,l:) = u(i,l:)/scl
      s = s+sum(u(i,l:)**2)

      f = u(i,l)
      g = -sign(sqrt(s),f)
      h = f*g-s
      u(i,l) = f-g

      rv1(l:) = u(i,l:)/h

      if (i /= m) then
        do j=l,m
          s = sum(u(j,l:)*u(i,l:))

          u(j,l:) = u(j,l:)+s*rv1(l:)
        end do
      end if

      u(i,l:) = scl*u(i,l:)
    end if
  end if

  x = max(x,abs(w(i))+abs(rv1(i)))
end do
! .......... accumulation of right-hand transformations ..........
if (matv) then
  ! .......... for i=n step -1 until 1 do -- ..........
  do ii=1,n
    i = n+1-ii

    if (i /= n) then
      if (abs(g) >= thresh) then
        ! .......... double division avoids possible underflow ..........
        v(l:n,i) = (u(i,l:n)/u(i,l))/g

        do j=l,n
          s = sum(u(i,l:)*v(l:n,j))

          v(l:n,j) = v(l:n,j)+s*v(l:n,i)
        end do
      end if

      v(i,l:n) = Zero
      v(l:n,i) = Zero
    end if

    v(i,i) = One
    g = rv1(i)
    l = i
  end do
end if
! .......... accumulation of left-hand transformations ..........
if (matu) then
  ! ..........for i=min(m,n) step -1 until 1 do -- ..........
  mn = n
  if (m < n) mn = m

  do ii=1,mn
    i = mn+1-ii
    l = i+1
    g = w(i)

    if (i /= n) u(i,l:) = Zero

    if (abs(g) < thresh) then
      u(i:m,i) = Zero
    else
      if (i /= mn) then
        do j=l,n
          s = sum(u(l:m,i)*u(l:m,j))
          ! .......... double division avoids possible underflow ..........
          f = (s/u(i,i))/g

          u(i:m,j) = u(i:m,j)+f*u(i:m,i)
        end do
      end if

      u(i:m,i) = u(i:m,i)/g
    end if

    u(i,i) = u(i,i)+One
  end do
end if
! .......... diagonalization of the bidiagonal form ..........
tst1 = x
! .......... for k=n step -1 until 1 do -- ..........
do kk=1,n
  k1 = n-kk
  k = k1+1
  its = 0

  do
    ! .......... test for splitting. for l=k step -1 until 1 do -- ..........
    skip = .false.
    do ll=1,k
      l1 = k-ll
      l = l1+1
      tst2 = tst1+abs(rv1(l))
      if (abs(tst2-tst1) < thresh) then
        skip = .true.
        exit
      end if
      ! .......... rv1(1) is always zero, so there is no exit through the bottom of the loop ..........
      tst2 = tst1+abs(w(l1))
      if (abs(tst2-tst1) < thresh) exit
    end do
    if (.not. skip) then
      ! .......... cancellation of rv1(l) if l greater than 1 ..........
      c = Zero
      s = One

      do i=l,k
        f = s*rv1(i)
        rv1(i) = c*rv1(i)
        tst2 = tst1+abs(f)
        if (abs(tst2-tst1) < thresh) exit
        g = w(i)
        h = pythag(f,g)
        w(i) = h
        c = g/h
        s = -f/h

        if (matu) then
          do j=1,m
            y = u(j,l1)
            z = u(j,i)
            u(j,l1) = y*c+z*s
            u(j,i) = -y*s+z*c
          end do
        end if

      end do
    end if
    ! .......... test for convergence ..........
    z = w(k)
    if (l == k) then
      ! .......... convergence ..........
      if (z >= Zero) exit
      ! .......... w(k) is made non-negative ..........
      w(k) = -z

      if (matv) v(1:n,k) = -v(1:n,k)
      exit
    else
      ! .......... shift from bottom 2 by 2 minor ..........
      if (its == 30) then
        ! .......... set error -- no convergence to a singular value after 30 iterations ..........
        ierr = k
        return
      end if
      its = its+1
      x = w(l)
      y = w(k1)
      g = rv1(k1)
      h = rv1(k)
      f = Half*(((g+z)/h)*((g-z)/y)+y/h-h/y)
      g = pythag(f,One)
      f = x-(z/x)*z+(h/x)*(y/(f+sign(g,f))-h)
      ! .......... next qr transformation ..........
      c = One
      s = One

      do i1=l,k1
        i = i1+1
        g = rv1(i)
        y = w(i)
        h = s*g
        g = c*g
        z = pythag(f,h)
        rv1(i1) = z
        c = f/z
        s = h/z
        f = x*c+g*s
        g = -x*s+g*c
        h = y*s
        y = y*c

        if (matv) then
          do j=1,n
            x = v(j,i1)
            z = v(j,i)
            v(j,i1) = x*c+z*s
            v(j,i) = -x*s+z*c
          end do
        end if

        z = pythag(f,h)
        w(i1) = z
        ! .......... rotation can be arbitrary if z is zero ..........
        if (abs(z) >= thresh) then
          c = f/z
          s = h/z
        end if
        f = c*g+s*y
        x = -s*g+c*y

        if (matu) then
          do j=1,m
            y = u(j,i1)
            z = u(j,i)
            u(j,i1) = y*c+z*s
            u(j,i) = -y*s+z*c
          end do
        end if

      end do

      rv1(l) = Zero
      rv1(k) = f
      w(k) = x
    end if

  end do

end do

return

end subroutine svd

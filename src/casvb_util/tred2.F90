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

subroutine tred2(nm,n,a,d,e,z)
! this subroutine is a translation of the algol procedure tred2,
! num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
! handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
! this subroutine reduces a real symmetric matrix to a
! symmetric tridiagonal matrix using and accumulating
! orthogonal similarity transformations.
!
! on input
!
!    nm must be set to the row dimension of two-dimensional
!      array parameters as declared in the calling program
!      dimension statement.
!
!    n is the order of the matrix.
!
!    a contains the real symmetric input matrix.  only the
!      lower triangle of the matrix need be supplied.
!
! on output
!
!    d contains the diagonal elements of the tridiagonal matrix.
!
!    e contains the subdiagonal elements of the tridiagonal
!      matrix in its last n-1 positions.  e(1) is set to zero.
!
!    z contains the orthogonal transformation matrix
!      produced in the reduction.
!
!    a and z may coincide.  if distinct, a is unaltered.
!
! questions and comments should be directed to burton s. garbow,
! mathematics and computer science div, argonne national laboratory
!
! this version dated august 1983.
!
! Updated to Fortran 90+ (Sep. 2023)
! ----------------------------------------------------------------------

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nm, n
real(kind=wp) :: a(nm,n), d(n), e(n), z(nm,n)
integer(kind=iwp) :: i, ii, j, jp1, k, l
real(kind=wp) :: f, g, h, hh, scl

do i=1,n

  do j=i,n
    z(j,i) = a(j,i)
  end do

  d(i) = a(n,i)
end do

! .......... for i=n step -1 until 2 do -- ..........
do ii=2,n
  i = n+2-ii
  l = i-1
  h = Zero
  scl = Zero
  if (l >= 2) then
    ! .......... scale row (algol tol then not needed) ..........
    do k=1,l
      scl = scl+abs(d(k))
    end do
  end if

  if (scl == Zero) then
    e(i) = d(l)

    do j=1,l
      d(j) = z(l,j)
      z(i,j) = Zero
      z(j,i) = Zero
    end do

  else

    do k=1,l
      d(k) = d(k)/scl
      h = h+d(k)*d(k)
    end do

    f = d(l)
    g = -sign(sqrt(h),f)
    e(i) = scl*g
    h = h-f*g
    d(l) = f-g
    ! .......... form a*u ..........
    do j=1,l
      e(j) = Zero
    end do

    do j=1,l
      f = d(j)
      z(j,i) = f
      g = e(j)+z(j,j)*f
      jp1 = j+1

      do k=jp1,l
        g = g+z(k,j)*d(k)
        e(k) = e(k)+z(k,j)*f
      end do

      e(j) = g
    end do
    ! .......... form p ..........
    f = Zero

    do j=1,l
      e(j) = e(j)/h
      f = f+e(j)*d(j)
    end do

    hh = f/(h+h)
    ! .......... form q ..........
    do j=1,l
      e(j) = e(j)-hh*d(j)
    end do
    ! .......... form reduced a ..........
    do j=1,l
      f = d(j)
      g = e(j)

      do k=j,l
        z(k,j) = z(k,j)-f*e(k)-g*d(k)
      end do

      d(j) = z(l,j)
      z(i,j) = Zero
    end do

  end if

  d(i) = h
end do
! .......... accumulation of transformation matrices ..........
do i=2,n
  l = i-1
  z(n,l) = z(l,l)
  z(l,l) = One
  h = d(i)

  if (h /= Zero) then
    do k=1,l
      d(k) = z(k,i)/h
    end do

    do j=1,l
      g = Zero

      do k=1,l
        g = g+z(k,i)*z(k,j)
      end do

      do k=1,l
        z(k,j) = z(k,j)-g*d(k)
      end do
    end do
  end if

  do k=1,l
    z(k,i) = Zero
  end do

end do

do i=1,n
  d(i) = z(n,i)
  z(n,i) = Zero
end do

z(n,n) = One
e(1) = Zero

return

end subroutine tred2

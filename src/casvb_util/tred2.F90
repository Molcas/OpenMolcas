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
integer(kind=iwp), intent(in) :: nm, n
real(kind=wp), intent(inout) :: a(nm,n)
real(kind=wp), intent(out) :: d(n), e(n), z(nm,n)
integer(kind=iwp) :: i, ii, j, jp1, l
real(kind=wp) :: f, g, h, hh, scl

do i=1,n

  z(i:n,i) = a(i:n,i)

  d(i) = a(n,i)
end do

! .......... for i=n step -1 until 2 do -- ..........
do ii=2,n
  i = n+2-ii
  l = i-1
  h = Zero
  ! .......... scale row (algol tol then not needed) ..........
  scl = sum(abs(d(1:l)))

  if (scl == Zero) then
    e(i) = d(l)

    d(1:l) = z(l,1:l)
    z(i,1:l) = Zero
    z(1:l,i) = Zero

  else

    d(1:l) = d(1:l)/scl
    h = h+sum(d(1:l)**2)

    f = d(l)
    g = -sign(sqrt(h),f)
    e(i) = scl*g
    h = h-f*g
    d(l) = f-g
    ! .......... form a*u ..........
    e(1:l) = Zero

    do j=1,l
      f = d(j)
      z(j,i) = f
      g = e(j)+z(j,j)*f
      jp1 = j+1

      g = g+sum(z(jp1:l,j)*d(jp1:l))
      e(jp1:l) = e(jp1:l)+z(jp1:l,j)*f

      e(j) = g
    end do
    ! .......... form p ..........
    e(1:l) = e(1:l)/h
    f = sum(e(1:l)*d(1:l))

    hh = f/(h+h)
    ! .......... form q ..........
    e(1:l) = e(1:l)-hh*d(1:l)
    ! .......... form reduced a ..........
    do j=1,l
      f = d(j)
      g = e(j)

      z(j:l,j) = z(j:l,j)-f*e(j:l)-g*d(j:l)

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
    d(1:l) = z(1:l,i)/h

    do j=1,l
      g = Zero

      g = sum(z(1:l,i)*z(1:l,j))

      z(1:l,j) = z(1:l,j)-g*d(1:l)
    end do
  end if

  z(1:l,i) = Zero

end do

d(:) = z(n,:)
z(n,:) = Zero

z(n,n) = One
e(1) = Zero

return

end subroutine tred2

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

subroutine tred1(nm,n,a,d,e,e2)
! this subroutine is a translation of the algol procedure tred1,
! num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
! handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
! this subroutine reduces a real symmetric matrix
! to a symmetric tridiagonal matrix using
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
!    a contains information about the orthogonal trans-
!      formations used in the reduction in its strict lower
!      triangle.  the full upper triangle of a is unaltered.
!
!    d contains the diagonal elements of the tridiagonal matrix.
!
!    e contains the subdiagonal elements of the tridiagonal
!      matrix in its last n-1 positions.  e(1) is set to zero.
!
!    e2 contains the squares of the corresponding elements of e.
!      e2 may coincide with e if the squares are not needed.
!
! questions and comments should be directed to burton s. garbow,
! mathematics and computer science div, argonne national laboratory
!
! this version dated august 1983.
!
! Updated to Fortran 90+ (Sep. 2023)
! ----------------------------------------------------------------------

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nm, n
real(kind=wp) :: a(nm,n), d(n), e(n), e2(n)
integer(kind=iwp) :: i, ii, j, jp1, k, l
real(kind=wp) :: f, g, h, scl

do i=1,n
  d(i) = a(n,i)
  a(n,i) = a(i,i)
end do
! .......... for i=n step -1 until 1 do -- ..........
do ii=1,n
  i = n+1-ii
  l = i-1
  h = Zero
  scl = Zero
  ! .......... scale row (algol tol then not needed) ..........
  do k=1,l
    scl = scl+abs(d(k))
  end do

  if (scl == Zero) then

    do j=1,l
      d(j) = a(l,j)
      a(l,j) = a(i,j)
      a(i,j) = Zero
    end do

    e(i) = Zero
    e2(i) = Zero

  else

    do k=1,l
      d(k) = d(k)/scl
      h = h+d(k)*d(k)
    end do

    e2(i) = scl*scl*h
    f = d(l)
    g = -sign(sqrt(h),f)
    e(i) = scl*g
    h = h-f*g
    d(l) = f-g
    if (l /= 1) then
      ! .......... form a*u ..........
      do j=1,l
        e(j) = Zero
      end do

      do j=1,l
        f = d(j)
        g = e(j)+a(j,j)*f
        jp1 = j+1

        do k=jp1,l
          g = g+a(k,j)*d(k)
          e(k) = e(k)+a(k,j)*f
        end do

        e(j) = g
      end do
      ! .......... form p ..........
      f = Zero

      do j=1,l
        e(j) = e(j)/h
        f = f+e(j)*d(j)
      end do

      h = f/(h+h)
      ! .......... form q ..........
      do j=1,l
        e(j) = e(j)-h*d(j)
      end do
      ! .......... form reduced a ..........
      do j=1,l
        f = d(j)
        g = e(j)

        do k=j,l
          a(k,j) = a(k,j)-f*e(k)-g*d(k)
        end do

      end do
    end if

    do j=1,l
      f = d(j)
      d(j) = a(l,j)
      a(l,j) = a(i,j)
      a(i,j) = f*scl
    end do

  end if

end do

return

end subroutine tred1

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
! ----------------------------------------------------------------------

integer i, j, k, l, n, ii, nm, jp1
real*8 a(nm,n), d(n), e(n), z(nm,n)
real*8 f, g, h, hh, scale

do i=1,n

  do j=i,n
    z(j,i) = a(j,i)
  end do

  d(i) = a(n,i)
end do

if (n == 1) go to 510
! .......... for i=n step -1 until 2 do -- ..........
do ii=2,n
  i = n+2-ii
  l = i-1
  h = 0.0d0
  scale = 0.0d0
  if (l < 2) go to 130
  ! .......... scale row (algol tol then not needed) ..........
  do k=1,l
    scale = scale+abs(d(k))
  end do

  if (scale /= 0.0d0) go to 140
130 e(i) = d(l)

  do j=1,l
    d(j) = z(l,j)
    z(i,j) = 0.0d0
    z(j,i) = 0.0d0
  end do

  go to 290

140 do k=1,l
    d(k) = d(k)/scale
    h = h+d(k)*d(k)
  end do

  f = d(l)
  g = -sign(sqrt(h),f)
  e(i) = scale*g
  h = h-f*g
  d(l) = f-g
  ! .......... form a*u ..........
  do j=1,l
    e(j) = 0.0d0
  end do

  do j=1,l
    f = d(j)
    z(j,i) = f
    g = e(j)+z(j,j)*f
    jp1 = j+1
    if (l < jp1) go to 220

    do k=jp1,l
      g = g+z(k,j)*d(k)
      e(k) = e(k)+z(k,j)*f
    end do

220 e(j) = g
  end do
  ! .......... form p ..........
  f = 0.0d0

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
    z(i,j) = 0.0d0
  end do

290 d(i) = h
end do
! .......... accumulation of transformation matrices ..........
do i=2,n
  l = i-1
  z(n,l) = z(l,l)
  z(l,l) = 1.0d0
  h = d(i)
  if (h == 0.0d0) go to 380

  do k=1,l
    d(k) = z(k,i)/h
  end do

  do j=1,l
    g = 0.0d0

    do k=1,l
      g = g+z(k,i)*z(k,j)
    end do

    do k=1,l
      z(k,j) = z(k,j)-g*d(k)
    end do
  end do

380 do k=1,l
    z(k,i) = 0.0d0
  end do

end do

510 do i=1,n
  d(i) = z(n,i)
  z(n,i) = 0.0d0
end do

z(n,n) = 1.0d0
e(1) = 0.0d0

return

end subroutine tred2

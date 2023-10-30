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

subroutine elmhes(nm,n,low,igh,a,intx)
! this subroutine is a translation of the algol procedure elmhes,
! num. math. 12, 349-368(1968) by martin and wilkinson.
! handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
!
! given a real general matrix, this subroutine
! reduces a submatrix situated in rows and columns
! low through igh to upper hessenberg form by
! stabilized elementary similarity transformations.
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
!    a contains the input matrix.
!
! on output
!
!    a contains the hessenberg matrix.  the multipliers
!      which were used in the reduction are stored in the
!      remaining triangle under the hessenberg matrix.
!
!    intx contains information on the rows and columns
!      interchanged in the reduction.
!      only elements low through igh are used.
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
integer(kind=iwp), intent(in) :: nm, n, low, igh
real(kind=wp), intent(inout) :: a(nm,n)
integer(kind=iwp), intent(out) :: intx(igh)
integer(kind=iwp) :: i, j, kp1, la, m, mm1, mp1
real(kind=wp) :: x, y

la = igh-1
kp1 = low+1

do m=kp1,la
  mm1 = m-1
  x = Zero
  i = m

  do j=m,igh
    if (abs(a(j,mm1)) <= abs(x)) cycle
    x = a(j,mm1)
    i = j
  end do

  intx(m) = i
  if (i /= m) then
    ! .......... interchange rows and columns of a ..........
    do j=mm1,n
      y = a(i,j)
      a(i,j) = a(m,j)
      a(m,j) = y
    end do

    do j=1,igh
      y = a(j,i)
      a(j,i) = a(j,m)
      a(j,m) = y
    end do
    ! .......... end interchange ..........
  end if
  if (x == Zero) cycle
  mp1 = m+1

  do i=mp1,igh
    y = a(i,mm1)
    if (y == Zero) cycle
    y = y/x
    a(i,mm1) = y

    a(i,m:) = a(i,m:)-y*a(m,m:)

    a(1:igh,m) = a(1:igh,m)+y*a(1:igh,i)

  end do

end do

return

end subroutine elmhes

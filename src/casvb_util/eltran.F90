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

subroutine eltran(nm,n,low,igh,a,intx,z)
! this subroutine is a translation of the algol procedure elmtrans,
! num. math. 16, 181-204(1970) by peters and wilkinson.
! handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
!
! this subroutine accumulates the stabilized elementary
! similarity transformations used in the reduction of a
! real general matrix to upper hessenberg form by  elmhes.
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
!    a contains the multipliers which were used in the
!      reduction by  elmhes  in its lower triangle
!      below the subdiagonal.
!
!    intx contains information on the rows and columns
!      interchanged in the reduction by  elmhes.
!      only elements low through igh are used.
!
! on output
!
!    z contains the transformation matrix produced in the
!      reduction by  elmhes.
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
integer(kind=iwp), intent(in) :: nm, n, low, igh, intx(igh)
real(kind=wp), intent(in) :: a(nm,igh)
real(kind=wp), intent(out) :: z(nm,n)
integer(kind=iwp) :: i, j, kl, mm, mp, mp1

! .......... initialize z to identity matrix ..........
z(1:n,:) = Zero
do j=1,n
  z(j,j) = One
end do

kl = igh-low-1
! .......... for mp=igh-1 step -1 until low+1 do -- ..........
do mm=1,kl
  mp = igh-mm
  mp1 = mp+1

  z(mp1:igh,mp) = a(mp1:igh,mp-1)

  i = intx(mp)
  if (i == mp) cycle

  z(mp,mp:igh) = z(i,mp:igh)
  z(i,mp) = One
  z(i,mp+1:igh) = Zero

end do

return

end subroutine eltran

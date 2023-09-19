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

subroutine eltran(nm,n,low,igh,a,int,z)
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
!    int contains information on the rows and columns
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

integer i, j, n, kl, mm, mp, nm, igh, low, mp1
real*8 a(nm,igh), z(nm,n)
integer int(igh)

! .......... initialize z to identity matrix ..........
do j=1,n

  do i=1,n
    z(i,j) = 0.0d0
  end do

  z(j,j) = 1.0d0
end do

kl = igh-low-1
! .......... for mp=igh-1 step -1 until low+1 do -- ..........
do mm=1,kl
  mp = igh-mm
  mp1 = mp+1

  do i=mp1,igh
    z(i,mp) = a(i,mp-1)
  end do

  i = int(mp)
  if (i == mp) cycle

  do j=mp,igh
    z(mp,j) = z(i,j)
    z(i,j) = 0.0d0
  end do

  z(i,mp) = 1.0d0
end do

return

end subroutine eltran

!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1992, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine SetUpA(nZeta,A,Pxyz)
!***********************************************************************
!                                                                      *
! Object: to set up the transformation matrix from the local coordinate*
!         system which has the z-axis going through P to the global    *
!         coordinate system. Formula by P.-A. Malmqvist.               *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!              October '92.                                            *
!***********************************************************************

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta
real(kind=wp), intent(out) :: A(nZeta,3,3)
real(kind=wp), intent(in) :: Pxyz(nZeta,3)
integer(kind=iwp) :: iZeta
real(kind=wp) :: r, sgn, x, y, z

#ifdef _DEBUGPRINT_
call RecPrt(' In SetupA: Pxyz',' ',Pxyz,nZeta,3)
#endif

! Set up the transformation matrix

do iZeta=1,nZeta
  x = Pxyz(iZeta,1)
  y = Pxyz(iZeta,2)
  z = Pxyz(iZeta,3)
  r = sqrt(x**2+y**2+z**2)
  sgn = One
  if (z < Zero) then
    x = -x
    y = -y
    z = -z
    sgn = -sgn
  end if
  if (r == Zero) then
    A(iZeta,1,1) = One
    A(iZeta,2,1) = Zero
    A(iZeta,1,2) = A(iZeta,2,1)
    A(iZeta,3,1) = Zero
    A(iZeta,1,3) = A(iZeta,3,1)
    A(iZeta,2,2) = One
    A(iZeta,2,3) = Zero
    A(iZeta,3,2) = A(iZeta,2,3)
    A(iZeta,3,3) = One
  else
    A(iZeta,1,1) = sgn*(One-x**2/(r*(r+z)))
    A(iZeta,2,1) = sgn*(-x*y/(r*(r+z)))
    A(iZeta,1,2) = A(iZeta,2,1)
    A(iZeta,3,1) = sgn*(-x/r)
    A(iZeta,1,3) = A(iZeta,3,1)
    A(iZeta,2,2) = sgn*(One-y**2/(r*(r+z)))
    A(iZeta,2,3) = sgn*(-y/r)
    A(iZeta,3,2) = A(iZeta,2,3)
    A(iZeta,3,3) = sgn*(-z/r)
  end if
end do

#ifdef _DEBUGPRINT_
call RecPrt(' The transformation matrix',' ',A,nZeta,9)
#endif

end subroutine SetUpA

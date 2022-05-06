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

subroutine Angular_Grid()
!***********************************************************************
!                                                                      *
!     Computes data useful for the angular quadrature.                 *
!                                                                      *
!***********************************************************************

use nq_Structure, only: Info_Ang
use nq_Info

implicit real*8(a-h,o-z)
#include "itmax.fh"
#include "real.fh"
#include "debug.fh"
logical Check
! Statement function
Check(i,j) = iand(i,2**(j-1)) /= 0

!                                                                      *
!***********************************************************************
!                                                                      *
nAngularGrids = 0
if (Check(iOpt_Angular,3)) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Generate angular grid a la Lebedev

  call Lebedev_Grid(L_Quad)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else if (Check(iOpt_Angular,1)) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Generate angular grid a la Lobatto

  call Lobatto_Grid(L_Quad)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Generate angular grid from Gauss and Gauss-Legendre quadrature

  call GGL_Grid(L_Quad)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (Debug) then
  do iSet=1,nAngularGrids
    nGP = Info_Ang(iSet)%nPoints
    l = Info_Ang(iSet)%L_eff
    write(6,*) 'l=',l
    call RecPrt('Angular grid',' ',Info_Ang(iSet)%R,4,nGP)
  end do
end if

return

end subroutine Angular_Grid

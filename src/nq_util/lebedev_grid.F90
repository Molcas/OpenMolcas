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

subroutine Lebedev_Grid(L_Max)
!***********************************************************************
!                                                                      *
!     Computes datas useful for the angular quadrature.                *
!                                                                      *
!***********************************************************************

use nq_Structure, only: Info_Ang
use nq_Info

implicit real*8(a-h,o-z)
#include "real.fh"
parameter(nSet=11)
integer Lebedev_order(nSet)
data Lebedev_order/5,7,11,17,23,29,35,41,47,53,59/
!                                                                      *
!***********************************************************************
!                                                                      *
interface
  subroutine Do_GGL(L_Eff,nPoints,R)
    implicit none
    integer L_Eff, nPoints
    real*8, allocatable :: R(:,:)
  end subroutine Do_GGL
  subroutine Do_Lebedev(L_Eff,nPoints,R)
    implicit none
    integer L_Eff, nPoints
    real*8, allocatable :: R(:,:)
  end subroutine Do_Lebedev
end interface

!                                                                      *
!***********************************************************************
!                                                                      *
! Use the GGL grids to minimize the number of grid points.

if (L_Max < 3) return
nAngularGrids = nAngularGrids+1
Info_Ang(nAngularGrids)%L_Eff = 3
call Do_GGL(3,Info_Ang(nAngularGrids)%nPoints,Info_Ang(nAngularGrids)%R)
!                                                                      *
!***********************************************************************
!                                                                      *
! Generate angular grid a la Lebedev

do iSet=1,nSet
  if (Lebedev_order(iSet) <= L_Max) then
    nAngularGrids = nAngularGrids+1
    L_Eff = Lebedev_order(iSet)

    Info_Ang(nAngularGrids)%L_Eff = L_Eff
    call Do_Lebedev(L_Eff,Info_Ang(nAngularGrids)%nPoints,Info_Ang(nAngularGrids)%R)
  else

    return

  end if
end do
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine Lebedev_Grid

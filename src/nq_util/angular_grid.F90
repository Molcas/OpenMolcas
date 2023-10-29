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

use nq_Info, only: iOpt_Angular, L_Quad, nAngularGrids
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
use nq_Structure, only: Info_Ang
use Definitions, only: iwp, u6
#endif

implicit none
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: iSet, l, nGP
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
nAngularGrids = 0
if (btest(iOpt_Angular,2)) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Generate angular grid a la Lebedev

  call Lebedev_Grid(L_Quad)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else if (btest(iOpt_Angular,0)) then
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
#ifdef _DEBUGPRINT_
do iSet=1,nAngularGrids
  nGP = Info_Ang(iSet)%nPoints
  l = Info_Ang(iSet)%L_eff
  write(u6,*) 'l=',l
  call RecPrt('Angular grid',' ',Info_Ang(iSet)%R,4,nGP)
end do
#endif

return

end subroutine Angular_Grid

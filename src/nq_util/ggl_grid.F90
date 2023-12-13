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

subroutine GGL_Grid(L_Max)
!***********************************************************************
!                                                                      *
!     Computes data useful for the angular quadrature.                 *
!                                                                      *
!***********************************************************************

use do_grid, only: Do_GGL
use nq_Structure, only: Info_Ang
use nq_Info, only: nAngularGrids
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: L_Max
integer(kind=iwp) :: L_Eff

!                                                                      *
!***********************************************************************
!                                                                      *
! Generate angular grid from Gauss and Gauss-Legendre quadrature
!
!-- Theta (polar angle): 0 <= theta <= pi
!   Gauss-Legendre Quadrature (L_Quad+1)/2 points
!-- Phi (azimuthal angle): 0 <= phi <= 2*pi
!   Gauss-Quadrature (L_Quad+1) points

do L_Eff=1,L_Max
  nAngularGrids = nAngularGrids+1

  Info_Ang(nAngularGrids)%L_eff = L_eff
  call Do_GGL(L_Eff,Info_Ang(nAngularGrids)%nPoints,Info_Ang(nAngularGrids)%R)

end do ! L_Eff
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine GGL_Grid

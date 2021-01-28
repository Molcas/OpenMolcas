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
! Copyright (C) Yannick Carissan                                       *
!***********************************************************************
!  SetVectTrans
!
!> @brief
!>   Computes the translation vector from geometry 2 to geometry 1
!> @author Y. Carissan
!>
!> @details
!> Computes the translation vector from geometry 2 to geometry 1.
!>
!> @param[in]  nat1  Number of atoms in \p geom1
!> @param[in]  Geom1 First Geometry to be considered, xyz coordinates
!> @param[in]  XYZ1  atom index array in the first geometry
!> @param[in]  nat2  Number of atoms in \p geom2
!> @param[in]  Geom2 Second Geometry to be considered, xyz coordinates
!> @param[in]  XYZ2  atom index array in the second geometry
!> @param[out] V     Output vector
!***********************************************************************

subroutine SetVectTrans(nat1,Geom1,XYZ1,nat2,Geom2,XYZ2,V)

use Quater_globals, only: debug
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nat1, nat2, XYZ1(3), XYZ2(3)
real(kind=wp), intent(in) :: Geom1(3,nat1), Geom2(3,nat2)
real(kind=wp), intent(out) :: V(3)
real(kind=wp) :: O1(3), O2(3)

O1(:) = Geom1(:,XYZ1(1))
O2(:) = Geom2(:,XYZ2(1))
V(:) = O1(:) - O2(:)

if (debug) then
  call RecPrt('Vtrans',' ',V,3,1)
end if

return

end subroutine SetVectTrans

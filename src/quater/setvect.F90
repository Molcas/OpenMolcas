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
!  setvect
!
!> @brief
!>   Computes the vectors \f$ V_1 = \mathit{XYZ}(2)-\mathit{XYZ}(1) \f$ and
!>   \f$ V_2 = \mathit{XYZ}(3)-\mathit{XYZ}(1) \f$ using the geometry \p geom
!> @author Y. Carissan
!>
!> @details
!> Computes the vectors \f$ V_1 = \mathit{XYZ}(2)-\mathit{XYZ}(1) \f$ and
!> \f$ V_2 = \mathit{XYZ}(3)-\mathit{XYZ}(1) \f$ using the geometry \p geom.
!>
!> @param[in]  natoms Number of atoms
!> @param[in]  Geom   Geometry to be considered, xyz coordinates
!> @param[in]  XYZ    atom index array
!> @param[out] V1     output vector
!> @param[out] V2     output vector
!***********************************************************************

subroutine setvect(natoms,Geom,XYZ,V1,V2)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: natoms, XYZ(3)
real(kind=wp), intent(in) :: Geom(3,natoms)
real(kind=wp), intent(out) :: V1(3), V2(3)
real(kind=wp) :: O(3), A(3), B(3)

O(:) = Geom(:,XYZ(1))
A(:) = Geom(:,XYZ(2))
B(:) = Geom(:,XYZ(3))

V1(:) = A(:)-O(:)
V2(:) = B(:)-O(:)

return

end subroutine setvect

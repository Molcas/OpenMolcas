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
!  translategeom
!
!> @brief
!>   Performs the translation of \p Geom with the \p Vtrans vector and stores the result in \p Geom
!> @author Y. Carissan
!>
!> @details
!> Performs the translation of \p Geom with the \p Vtrans vector
!> and stores the result in \p Geom.
!>
!> @param[in]     Vtrans Input translation vector
!> @param[in]     nat    Number of atoms
!> @param[in,out] Geom   Input and output geometry, xyz coordinates
!***********************************************************************

subroutine translategeom(Vtrans,nat,Geom)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nat
real(kind=wp), intent(in) :: Vtrans(3)
real(kind=wp), intent(inout) :: Geom(3,nat)
integer(kind=iwp) :: iat

do iat=1,nat
  Geom(:,iat) = Geom(:,iat) + Vtrans(:)
end do

return

end subroutine translategeom

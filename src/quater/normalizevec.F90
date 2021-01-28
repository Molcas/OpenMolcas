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
!  NormalizeVec
!
!> @brief
!>   Normalize the vector \p V
!> @author Y. Carissan
!>
!> @details
!> Normalize the vector \p V.
!>
!> @param[in,out] V Vector to be normalized
!***********************************************************************

subroutine normalizeVec(V)

use Definitions, only: wp, r8

implicit none
real(kind=wp), intent(inout) :: V(3)
real(kind=wp) :: norm
real(kind=r8), external :: dnrm2_

norm = dnrm2_(3,V,1)
V(:) = V(:)/norm

return

end subroutine normalizeVec

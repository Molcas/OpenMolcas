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
!  GetKandC
!
!> @brief
!>   Utility routine for quaternion resolution
!> @author Y. Carissan
!>
!> @details
!> Performs the following operation:
!>
!> \f[ K = (V_1-U_1) \times (V_2-U_2) \\
!>     C = K . (U_1 \times U_2) \f]
!>
!> @param[in]  U1 Input vector \f$ U_1 \f$
!> @param[in]  U2 Input vector \f$ U_2 \f$
!> @param[in]  V1 Input vector \f$ V_1 \f$
!> @param[in]  V2 Input vector \f$ V_2 \f$
!> @param[out] K  Output vector \f$ K \f$
!> @param[out] C  Output value \f$ C \f$
!***********************************************************************

subroutine GetKandC(U1,U2,V1,V2,K,C)

use Quater_globals, only: debug
use Definitions, only: wp, r8, u6

implicit none
real(kind=wp), intent(in) :: U1(3), U2(3), V1(3), V2(3)
real(kind=wp), intent(out) :: K(3), C
real(kind=wp) :: T1(3), T2(3)
real(kind=r8), external :: ddot_

T1(:) = V1(:) - U1(:)     ! T1 = V1-U1
T2(:) = V2(:) - U2(:)     ! T2 = V2-U2

call cross(T1,T2,K)       ! K = T1 x T2 = (V1-U1) x (V2-U2)
call cross(U1,U2,T1)      ! T1 = U1 x U2
C = ddot_(3,K,1,T1,1)     ! C = K . (U1 x U2)

if (debug) then
  call RecPrt('K',' ',K,3,1)
  write(u6,*) 'C',C
end if

return ! K and C

end subroutine GetKandC

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
! Copyright (C) 2000, Roland Lindh                                     *
!***********************************************************************

subroutine Overlap(mGrid,iSpin)
!***********************************************************************
!      Author:Roland Lindh, Department of Chemical Physics, University *
!             of Lund, SWEDEN. November 2000                           *
!***********************************************************************

use nq_Grid, only: F_xc, Rho, vRho
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: mGrid, iSpin
integer(kind=iwp) :: iGrid
real(kind=wp) :: d_alpha, d_beta, DTot, Rho_Min
real(kind=wp), parameter :: T_x = 1.0e-20_wp

!                                                                      *
!***********************************************************************
!                                                                      *

vRho(:,:) = Zero
Rho_Min = T_X*1.0e-2_wp
if (iSpin == 1) then
  ! iSpin=1
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  do iGrid=1,mGrid

    d_alpha = Rho(1,iGrid)
    DTot = Two*d_alpha
    if (DTot < T_X) cycle

    ! Accumulate contributions to the integrated density

    F_xc(iGrid) = F_xc(iGrid)+Dtot

    vRho(1,iGrid) = One

  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else
  ! iSpin=/=1
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  do iGrid=1,mGrid

    d_alpha = max(Rho_min,Rho(1,iGrid))
    d_beta = max(Rho_min,Rho(2,iGrid))
    DTot = d_alpha+d_beta
    if (DTot < T_X) cycle

    ! Accumulate contributions to the integrated density

    F_xc(iGrid) = F_xc(iGrid)+Dtot

    vRho(1,iGrid) = One
    vRho(2,iGrid) = One

  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine Overlap

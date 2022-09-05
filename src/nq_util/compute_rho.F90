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

function Compute_Rho(Weights,mGrid,iSpin)
!***********************************************************************
!      Author:Roland Lindh, Department of Chemical Physics, University *
!             of Lund, SWEDEN. November 2000                           *
!***********************************************************************

use nq_Grid, only: Rho
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: Compute_Rho
integer(kind=iwp), intent(in) :: mGrid, iSpin
real(kind=wp), intent(in) :: Weights(mGrid)
integer(kind=iwp) :: iGrid
real(kind=wp) :: d_alpha, d_beta, DTot

!                                                                      *
!***********************************************************************
!                                                                      *
Compute_Rho = Zero
if (iSpin == 1) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! iSpin=1

  do iGrid=1,mGrid

    d_alpha = Rho(1,iGrid)
    DTot = d_alpha

    ! Accumulate contributions to the integrated density

    Compute_Rho = Compute_Rho+DTot*Weights(iGrid)

  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! iSpin=/=1

  do iGrid=1,mGrid

    d_alpha = Rho(1,iGrid)
    d_beta = Rho(2,iGrid)
    DTot = d_alpha+d_beta

    ! Accumulate contributions to the integrated density

    Compute_Rho = Compute_Rho+DTot*Weights(iGrid)

  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *

return

end function Compute_Rho

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

function Compute_Grad(Weights,mGrid,iSpin)
!***********************************************************************
!      Author:Roland Lindh, Department of Chemical Physics, University *
!             of Lund, SWEDEN. November 2000                           *
!***********************************************************************

use nq_Grid, only: Sigma
use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: Compute_Grad
integer(kind=iwp), intent(in) :: mGrid, iSpin
real(kind=wp), intent(in) :: Weights(mGrid)
integer(kind=iwp) :: iGrid
real(kind=wp) :: Gmma

!                                                                      *
!***********************************************************************
!                                                                      *
Compute_Grad = Zero

if (iSpin == 1) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! iSpin == 1

  do iGrid=1,mGrid

    Gmma = sqrt(Sigma(1,iGrid))

    ! Accumulate contributions to the integrated Tau

    Compute_Grad = Compute_Grad+Two*Gmma*Weights(iGrid)

  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! iSpin /= 1

  do iGrid=1,mGrid

    Gmma = sqrt(Sigma(1,iGrid)+Two*Sigma(2,iGrid)+Sigma(3,iGrid))

    ! Accumulate contributions to the integrated density

    Compute_Grad = Compute_Grad+Gmma*Weights(iGrid)

  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *

return

end function Compute_Grad

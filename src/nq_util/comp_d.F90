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
! Copyright (C) Giovanni Li Manni                                      *
!***********************************************************************

function Comp_d(Weights,mGrid,Rho,nRho,iSpin,iSwitch)
!***********************************************************************
!                                                                      *
! Object: integrate densities (alpha, beta, total, gradients....)      *
!         the object integrated is dictaded by iSwitch value:          *
!         iSwitch = 0  .... total density                              *
!         iSwitch = 1  .... alpha density                              *
!         iSwitch = 2  .... beta density                               *
!                                                                      *
! Author: G. Li Manni... taking Sir R. Lindh as model                  *
!***********************************************************************

use Constants, only: Zero, Two, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: Comp_d
integer(kind=iwp), intent(in) :: mGrid, nRho, iSpin, iSwitch
real(kind=wp), intent(in) :: Weights(mGrid), Rho(nRho,mGrid)
integer(kind=iwp) :: iGrid
real(kind=wp) :: d_alpha, d_beta, DTot

Comp_d = Zero
if (iSpin == 1) then
  !*********************************************************************
  ! iSpin == 1
  !*********************************************************************
  do iGrid=1,mGrid
    d_alpha = Half*Rho(1,iGrid)
    if (iSwitch == 1) then
      DTot = d_alpha
    else if (iSwitch == 2) then
      DTot = d_alpha
    else !if (iSwitch == 0) then
      DTot = Two*d_alpha
    end if
    Comp_d = Comp_d+DTot*Weights(iGrid)
  end do
else
  !*********************************************************************
  ! iSpin /= 1
  !*********************************************************************
  do iGrid=1,mGrid
    d_alpha = Rho(1,iGrid)
    d_beta = Rho(2,iGrid)
    if (iSwitch == 1) then
      DTot = d_alpha
    else if (iSwitch == 2) then
      DTot = d_beta
    else !if (iSwitch == 0) then
      DTot = d_alpha+d_beta
    end if
    Comp_d = Comp_d+DTot*Weights(iGrid)
  end do
end if
!***********************************************************************

return

end function Comp_d

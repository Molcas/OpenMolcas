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

subroutine Do_NucAtt_(mGrid,iSpin,Grid,RA,ZA,mCenter)
!***********************************************************************
!      Author:Roland Lindh, Department of Chemical Physics, University *
!             of Lund, SWEDEN. November 2000                           *
!***********************************************************************

use nq_Grid, only: F_xc, Rho, vRho
use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: mGrid, iSpin, mCenter
real(kind=wp), intent(in) :: Grid(3,mGrid), RA(3,mCenter), ZA(mCenter)
integer(kind=iwp) :: i, iGrid
real(kind=wp) :: Attr, d_alpha, d_beta, DTot, Fact, x, y, z

!                                                                      *
!***********************************************************************
!                                                                      *

vRho(:,:) = Zero
if (iSpin == 1) then
  ! iSpin=1
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  do iGrid=1,mGrid

    d_alpha = Rho(1,iGrid)
    DTot = Two*d_alpha

    ! Accumulate contributions to the nuclear attraction energy

    Attr = Zero
    do i=1,mCenter
      x = Grid(1,iGrid)-RA(1,i)
      y = Grid(2,iGrid)-RA(2,i)
      z = Grid(3,iGrid)-RA(3,i)
      Fact = ZA(i)/sqrt(x**2+y**2+z**2)
      Attr = Attr+Fact
    end do
    F_xc(iGrid) = F_xc(iGrid)-Attr*DTot

    vRho(1,iGrid) = -Attr

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

    d_alpha = Rho(1,iGrid)
    d_beta = Rho(2,iGrid)
    DTot = d_alpha+d_beta

    ! Accumulate contributions to the nuclear attraction energy

    Attr = Zero
    do i=1,mCenter
      x = Grid(1,iGrid)-RA(1,i)
      y = Grid(2,iGrid)-RA(2,i)
      z = Grid(3,iGrid)-RA(3,i)
      Fact = ZA(i)/sqrt(x**2+y**2+z**2)
      Attr = Attr+Fact
    end do
    F_xc(iGrid) = F_xc(iGrid)-Attr*DTot

    vRho(1,iGrid) = -Attr
    vRho(2,iGrid) = -Attr

  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine Do_NucAtt_

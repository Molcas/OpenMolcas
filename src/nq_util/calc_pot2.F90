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
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************

! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Dec. 08, 2021, created this file.               *
! ****************************************************************
subroutine Calc_Pot2(Pot2,mGrid,Pi,nPi)

use nq_Grid, only: Weights
use nq_pdft, only: d2RdRhodPi, d2ZdR2, dEdPi, dEdPiMO, dEdPix, dEdPiy, dEdPiz, dF_dRhoamb, dF_dRhoxamb, dF_dRhoyamb, dF_dRhozamb, &
                   dRdPi, dZdR, GdEdPiMO, GradRdFdRho, GradRhodFdRho, lft, lGGA, MOas, MOax, MOay, MOaz, Pass1, Pass2, Pass3, RhoAB
use nq_Info, only: nPot2
use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: Pot2(nPot2)
integer(kind=iwp), intent(in) :: mGrid, nPi
real(kind=wp), intent(in) :: Pi(nPi,mGrid)
integer(kind=iwp) :: iGrid
real(kind=wp) :: ftggaterm, ggaterm, predEdPip
real(kind=wp), parameter :: ThrsPi = 1.0e-30_wp

if (lGGA .and. lft) then
  dEdPix(:) = Zero
  dEdPiy(:) = Zero
  dEdPiz(:) = Zero
  GdEdPiMO(:,:) = Zero
end if

do iGrid=1,mGrid
  if (Pass1(iGrid) .and. (Pi(1,iGrid) > ThrsPi)) then
    if (Pass2(iGrid) .or. Pass3(iGrid)) then
      if (lGGA) then
        ggaterm = GradRhodFdRho(iGrid)
        if (lft) then
          ftggaterm = (d2ZdR2(iGrid)*dRdPi(iGrid)*GradRdFdRho(iGrid)+ &
                       d2RdRhodPi(iGrid)*dZdR(iGrid)*GradRhodFdRho(iGrid))*RhoAB(iGrid)
          predEdPip = RhoAB(iGrid)*dZdR(iGrid)*dRdPi(iGrid)*Weights(iGrid)
          dEdPix(iGrid) = predEdPip*dF_dRhoxamb(iGrid)
          dEdPiy(iGrid) = predEdPip*dF_dRhoyamb(iGrid)
          dEdPiz(iGrid) = predEdPip*dF_dRhozamb(iGrid)
        else
          ftggaterm = Zero
        end if
      else
        ggaterm = Zero
        ftggaterm = Zero
      end if
      dEdPi(iGrid) = Weights(iGrid)*(dZdR(iGrid)*dRdPi(iGrid)*(RhoAB(iGrid)*dF_dRhoamb(iGrid)+ggaterm)+ftggaterm)
    else
      dEdPi(iGrid) = Zero
    end if
  else
    dEdPi(iGrid) = Zero
  end if
end do

dEdPi(:) = Half*dEdPi
if (lGGA .and. lft) then
  dEdPix(:) = Half*dEdPix
  dEdPiy(:) = Half*dEdPiy
  dEdPiz(:) = Half*dEdPiz
end if

do iGrid=1,mGrid
  dEdPiMO(iGrid,:) = MOas(iGrid,:)*dEdPi(iGrid)
end do

if (lft .and. lGGA) then
  do iGrid=1,mGrid
    GdEdPiMO(iGrid,:) = GdEdPiMO(iGrid,:)+dEdPix(iGrid)*MOax(iGrid,:)+dEdPiy(iGrid)*MOay(iGrid,:)+dEdPiz(iGrid)*MOaz(iGrid,:)
  end do
  dEdPiMO(:,:) = dEdPiMO+GdEdPiMO
end if

! dEdPiMO is practically (Phi_p*dEdPi+Phi_p'*dEdPi')
! The subroutine below calculates
! (Phi_p*dEdPi+Phi_p'*dEdPi')*Phi_u*Phi_v*Phi_x
call Calc_Pot2_Inner(Pot2,mGrid,dEdPiMO,MOas,MOas,MOas,.false.)

if (lft .and. lGGA) then
  call Calc_Pot2_Inner(Pot2,mGrid,MOas,MOas,MOas,GdEdPiMO,.true.)
end if

return

end subroutine Calc_Pot2

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
use nq_pdft
use nq_Info

#include "stdalloc.fh"
! Input
integer mGrid, nPi
real*8, dimension(nPi,mGrid) :: Pi
! Output
real*8, dimension(nPot2) :: Pot2
! Internal
integer iGrid, nGOrb
real*8 ThrsPi, ggaterm, ftggaterm, predEdPip

ThrsPi = 1.0d-30

if (lGGA .and. lft) then
  call FZero(dEdPix,mGrid)
  call FZero(dEdPiy,mGrid)
  call FZero(dEdPiz,mGrid)
  call FZero(GdEdPiMO,mGrid*nOrbt)
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
          ftggaterm = 0.0d0
        end if
      else
        ggaterm = 0.0d0
        ftggaterm = 0.0d0
      end if
      dEdPi(iGrid) = Weights(iGrid)*(dZdR(iGrid)*dRdPi(iGrid)*(RhoAB(iGrid)*dF_dRhoamb(iGrid)+ggaterm)+ftggaterm)
    else
      dEdPi(iGrid) = 0.0d0
    end if
  else
    dEdPi(iGrid) = 0.0d0
  end if
end do
nGOrb = mGrid*nOrbt

call DSCal_(mGrid,0.5d0,dEdPi,1)
if (lGGA .and. lft) then
  call DSCal_(mGrid,0.5d0,dEdPix,1)
  call DSCal_(mGrid,0.5d0,dEdPiy,1)
  call DSCal_(mGrid,0.5d0,dEdPiz,1)
end if

call DCopy_(nGOrb,MOas,1,dEdPiMO,1)

do iGrid=1,mGrid
  call DScal_(nOrbt,dEdPi(iGrid),dEdPiMO(iGrid),mGrid)
end do

if (lft .and. lGGA) then
  do iGrid=1,mGrid
    call DAXpY_(nOrbt,dEdPix(iGrid),MOax(iGrid),mGrid,GdEdPiMO(iGrid),mGrid)
    call DAXpY_(nOrbt,dEdPiy(iGrid),MOay(iGrid),mGrid,GdEdPiMO(iGrid),mGrid)
    call DAXpY_(nOrbt,dEdPiz(iGrid),MOaz(iGrid),mGrid,GdEdPiMO(iGrid),mGrid)
  end do
  call DAXpY_(nGOrb,1.0d0,GdEdPiMO,1,dEdPiMO,1)
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

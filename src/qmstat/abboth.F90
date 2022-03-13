!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

! Routine for the case where both centres are diffuse. Since these
! formulas are pretty nasty and apparently with little general
! structure, each type of interaction is hard-coded.
subroutine ABBoth(iLA,iLB,dMulA,Tau,dKappa,Rho,RhoA,RhoB,Rinv,lTooSmall,Colle)

implicit real*8(a-h,o-z)
parameter(MxMltp=2)
dimension dMulA((MxMltp+1)*(MxMltp+2)/2), Colle(3)
logical lTooSmall
#include "warnings.h"

! To calculate the interaction Sigma is the product of both multipoles
! in A and in B but since we need potential, field and field gradient
! for the QM system whe do not multiply for multipoles in B, but we
! have to take into account to move the result for the original
! coordinate system in QmStat.

do i=1,3
  Colle(i) = 0.0d0
end do

if ((iLA == 0) .and. (iLB == 0)) then
  ! s-s interaction. There is only sigma-components, hence simple.

  Sigma = dMulA(1)
  if (lTooSmall) then
    Ex = exp((-2.0d0)*Rho)
    Colle(1) = Sigma*CoulT0_1(Rho,Rinv,Ex)
  else
    ExA = exp((-2.0d0)*RhoA)
    ExB = exp((-2.0d0)*RhoB)
    Colle(1) = Sigma*CoulTN_1(RhoA,RhoB,dKappa,Rinv,ExA,ExB)
  end if

else if ((iLA == 1) .and. (iLB == 0)) then
  ! s-p interaction. Only the z-component of the dipole interacts
  ! through a sigma-interaction with the s-distribution. Observe
  ! that in the case that iLA > iLB, then the formulas by Roothan
  ! has to be reversed, i.e. RhoA and RhoB change place and
  ! Tau and Kappa changes sign.

  Sigma = dMulA(3)
  if (lTooSmall) then
    Ex = exp((-2.0d0)*Rho)
    Colle(1) = Sigma*CoulT0_2(Rho,Rinv,Ex)
  else
    ExA = exp((-2.0d0)*RhoA)
    ExB = exp((-2.0d0)*RhoB)
    Colle(1) = Sigma*CoulTN_2(Rho,-Tau,RhoB,RhoA,-dKappa,Rinv,ExB,ExA)
  end if
else if ((iLA == 0) .and. (iLB == 1)) then
  Sigma = dMulA(1)
  if (lTooSmall) then
    Ex = exp((-2.0d0)*Rho)
    Colle(1) = Sigma*CoulT0_2(Rho,Rinv,Ex)
  else
    ExA = exp((-2.0d0)*RhoA)
    ExB = exp((-2.0d0)*RhoB)
    Colle(1) = Sigma*CoulTN_2(Rho,Tau,RhoA,RhoB,dKappa,Rinv,ExA,ExB)
  end if

else if ((iLA == 1) .and. (iLB == 1)) then
  ! p-p interaction. The z-z combination gives a sigma-interaction,
  ! and the x-x and y-y combinations give pi-interactions.

  ! The sigma-component.

  Sigma = dMulA(3)
  if (lTooSmall) then
    Ex = exp((-2.0d0)*Rho)
    Colle(1) = Sigma*CoulT0_4(Rho,Rinv,Ex)
  else
    ExA = exp((-2.0d0)*RhoA)
    ExB = exp((-2.0d0)*RhoB)
    Colle(1) = Sigma*CoulTN_4(Rho,Tau,RhoA,RhoB,dKappa,Rinv,ExA,ExB)
  end if

  ! The two pi-components.

  Pi1 = dMulA(1)
  Pi2 = dMulA(2)
  if (lTooSmall) then
    Ex = exp((-2.0d0)*Rho)
    Width = CoulT0_5(Rho,Rinv,Ex)
    Colle(2) = Pi1*Width
    Colle(3) = Pi2*Width
  else
    ExA = exp((-2.0d0)*RhoA)
    ExB = exp((-2.0d0)*RhoB)
    Width = CoulTN_5(Rho,Tau,RhoA,RhoB,dKappa,Rinv,ExA,ExB)
    Colle(2) = Pi1*Width
    Colle(3) = Pi2*Width
  end if

else
  ! Higher angular momentum interactions.

  write(6,*) 'Too high angular momentum'
  write(6,*) 'at least you start to implement.'
  call Quit(_RC_IO_ERROR_READ_)
end if

return

end subroutine ABBoth

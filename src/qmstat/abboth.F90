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
subroutine ABBoth(iLA,iLB,dMulA,dKappa,Rho,RhoA,RhoB,Rinv,lTooSmall,Colle)

use Index_Functions, only: nTri_Elem1
use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6

! Maximum multipole implemented
#define _MxM_ 2

implicit none
integer(kind=iwp), intent(in) :: iLA, iLB
real(kind=wp), intent(in) :: dMulA(nTri_Elem1(_MxM_)), dKappa, Rho, RhoA, RhoB, Rinv
logical(kind=iwp), intent(in) :: lTooSmall
real(kind=wp), intent(out) :: Colle(3)
real(kind=wp) :: Ex, ExA, ExB, Pi1, Pi2, Sigma, Width
real(kind=wp), external :: CoulT0_1, CoulT0_2, CoulT0_4, CoulT0_5, CoulTN_1, CoulTN_2, CoulTN_4, CoulTN_5
#include "warnings.h"

! To calculate the interaction Sigma is the product of both multipoles
! in A and in B but since we need potential, field and field gradient
! for the QM system whe do not multiply for multipoles in B, but we
! have to take into account to move the result for the original
! coordinate system in QmStat.

Colle(:) = Zero

if ((iLA == 0) .and. (iLB == 0)) then
  ! s-s interaction. There is only sigma-components, hence simple.

  Sigma = dMulA(1)
  if (lTooSmall) then
    Ex = exp(-Two*Rho)
    Colle(1) = Sigma*CoulT0_1(Rho,Rinv,Ex)
  else
    ExA = exp(-Two*RhoA)
    ExB = exp(-Two*RhoB)
    Colle(1) = Sigma*CoulTN_1(RhoA,RhoB,dKappa,Rinv,ExA,ExB)
  end if

else if ((iLA == 1) .and. (iLB == 0)) then
  ! s-p interaction. Only the z-component of the dipole interacts
  ! through a sigma-interaction with the s-distribution. Observe
  ! that in the case that iLA > iLB, then the formulas by Roothan
  ! has to be reversed, i.e. RhoA and RhoB change place and
  ! Kappa changes sign.

  Sigma = dMulA(3)
  if (lTooSmall) then
    Ex = exp(-Two*Rho)
    Colle(1) = Sigma*CoulT0_2(Rho,Rinv,Ex)
  else
    ExA = exp(-Two*RhoA)
    ExB = exp(-Two*RhoB)
    Colle(1) = Sigma*CoulTN_2(RhoB,RhoA,-dKappa,Rinv,ExB,ExA)
  end if
else if ((iLA == 0) .and. (iLB == 1)) then
  Sigma = dMulA(1)
  if (lTooSmall) then
    Ex = exp(-Two*Rho)
    Colle(1) = Sigma*CoulT0_2(Rho,Rinv,Ex)
  else
    ExA = exp(-Two*RhoA)
    ExB = exp(-Two*RhoB)
    Colle(1) = Sigma*CoulTN_2(RhoA,RhoB,dKappa,Rinv,ExA,ExB)
  end if

else if ((iLA == 1) .and. (iLB == 1)) then
  ! p-p interaction. The z-z combination gives a sigma-interaction,
  ! and the x-x and y-y combinations give pi-interactions.

  ! The sigma-component.

  Sigma = dMulA(3)
  if (lTooSmall) then
    Ex = exp(-Two*Rho)
    Colle(1) = Sigma*CoulT0_4(Rho,Rinv,Ex)
  else
    ExA = exp(-Two*RhoA)
    ExB = exp(-Two*RhoB)
    Colle(1) = Sigma*CoulTN_4(RhoA,RhoB,dKappa,Rinv,ExA,ExB)
  end if

  ! The two pi-components.

  Pi1 = dMulA(1)
  Pi2 = dMulA(2)
  if (lTooSmall) then
    Ex = exp(-Two*Rho)
    Width = CoulT0_5(Rho,Rinv,Ex)
    Colle(2) = Pi1*Width
    Colle(3) = Pi2*Width
  else
    ExA = exp(-Two*RhoA)
    ExB = exp(-Two*RhoB)
    Width = CoulTN_5(RhoA,RhoB,dKappa,Rinv,ExA,ExB)
    Colle(2) = Pi1*Width
    Colle(3) = Pi2*Width
  end if

else
  ! Higher angular momentum interactions.

  write(u6,*) 'Too high angular momentum'
  write(u6,*) 'at least you start to implement.'
  call Quit(_RC_IO_ERROR_READ_)
end if

return

end subroutine ABBoth

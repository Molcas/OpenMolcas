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
! Copyright (C) 2020, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************

subroutine ThetaOpt2(R,theta,deltaQ,SPair,NP,GD,Vee,G)

use rasscf_global, only: lRoots, NAC
use Constants, only: Zero

implicit none
integer NP
real*8, dimension(NP) :: theta
real*8 Change, deltaQ
integer, dimension(NP,2) :: SPair
real*8, dimension(lroots,lroots) :: R
real*8, dimension(lroots) :: Vee
real*8, dimension(LRoots*(LRoots+1)/2,NAC,NAC) :: GD
real*8, dimension(NAC,NAC,NAC,NAC) :: G
integer IP, I, J
#include "warnings.h"

deltaQ = Zero
do IP=1,NP
  I = SPair(IP,1)
  J = SPair(IP,2)
  call OptOneAngle2(theta(iP),change,R,GD,I,J,Vee,G)
  deltaQ = deltaQ+change
end do

do IP=NP-1,1,-1
  I = SPair(IP,1)
  J = SPair(IP,2)
  call OptOneAngle2(theta(iP),change,R,GD,I,J,Vee,G)
  deltaQ = deltaQ+change
end do

end subroutine ThetaOpt2

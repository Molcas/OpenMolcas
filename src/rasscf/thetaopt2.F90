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
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: NP
real(kind=wp) :: R(lRoots,lRoots), theta(NP), deltaQ, GD(lRoots*(lRoots+1)/2,NAC,NAC), Vee(lRoots), G(NAC,NAC,NAC,NAC)
integer(kind=iwp) :: SPair(NP,2)
integer(kind=iwp) :: I, IP, J
real(kind=wp) :: Change

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

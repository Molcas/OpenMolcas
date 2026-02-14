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

subroutine ThetaOpt(FRot,theta,SumVee,StatePair,NPairs,DDg)

use rasscf_global, only: lRoots
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: NPairs
real(kind=wp) :: FRot(lRoots,lRoots), theta(NPairs), SumVee, DDG(lRoots,lRoots,lRoots,lRoots)
integer(kind=iwp) :: StatePair(NPairs,2)
integer(kind=iwp) :: IPair, IState, JState

do IPair=1,NPairs
  IState = StatePair(IPair,1)
  JState = StatePair(IPair,2)
  call OptOneAngle(theta(iPair),SumVee,FRot,DDg,IState,JState,lRoots)
end do
do IPair=NPairs-1,1,-1
  IState = StatePair(IPair,1)
  JState = StatePair(IPair,2)
  call OptOneAngle(theta(iPair),SumVee,FRot,DDg,IState,JState,lRoots)
end do

end subroutine ThetaOpt

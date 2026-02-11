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
      Subroutine ThetaOpt(FRot,theta,SumVee,StatePair,NPairs,DDg)
      use rasscf_global, only: lRoots
      Implicit None

#include "warnings.h"
      INTEGER NPairs
      Real*8 SumVee
      Real*8,DIMENSION(lRoots,lRoots,lRoots,lRoots)::DDG
      Real*8,DIMENSION(lroots,lroots)::FRot
      INTEGER,DIMENSION(NPairs,2)::StatePair
      Real*8,DIMENSION(NPairs)::theta

      INTEGER IPair,IState,JState
!      Real*8,DIMENSION(NPairs)::thetanew

      DO IPair=1,NPairs
       IState=StatePair(IPair,1)
       JState=StatePair(IPair,2)
       CALL                                                             &
     & OptOneAngle(theta(iPair),SumVee,FRot,DDg,IState,JState,lRoots)
      END DO
      DO IPair=NPairs-1,1,-1
       IState=StatePair(IPair,1)
       JState=StatePair(IPair,2)
       CALL                                                             &
     & OptOneAngle(theta(iPair),SumVee,FRot,DDg,IState,JState,lRoots)
      END DO
      END SUBROUTINE ThetaOpt

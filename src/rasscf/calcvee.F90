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
      Subroutine CalcVee(Vee,RMat,DDg)
      use rasscf_global, only: lRoots
      Implicit None


#include "warnings.h"
      Real*8,DIMENSION(lRoots,lRoots,lRoots,lRoots)::DDG
      Real*8,DIMENSION(lroots,lroots)::RMat
      Real*8,DIMENSION(lroots)::Vee
      INTEGER IState,iJ,iK,iL,iM
      DO IState=1,lRoots
       Vee(IState)=0.0d0
       Do iJ=1,lRoots
        Do iK=1,lRoots
         Do iL=1,lRoots
          Do iM=1,lRoots
          Vee(Istate)=Vee(IState)+RMat(IState,iJ)*RMat(IState,iK)*      &
     &RMat(IState,iL)*RMat(IState,iM)*DDG(iJ,iK,iL,iM)
          End Do
         End Do
        End Do
       End Do
       Vee(IState)=Vee(IState)/2
!       write(6,'(A,I2,A,F10.6)')'The classic coulomb energy for state ',
!     & IState,' is ',Vee(IState)
      END DO
      END SUBROUTINE CalcVee

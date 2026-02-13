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

subroutine CalcVee(Vee,RMat,DDg)

use rasscf_global, only: lRoots
use Constants, only: Zero

implicit none
real*8, dimension(lRoots,lRoots,lRoots,lRoots) :: DDG
real*8, dimension(lroots,lroots) :: RMat
real*8, dimension(lroots) :: Vee
integer IState, iJ, iK, iL, iM
#include "warnings.h"

do IState=1,lRoots
  Vee(IState) = Zero
  do iJ=1,lRoots
    do iK=1,lRoots
      do iL=1,lRoots
        do iM=1,lRoots
          Vee(Istate) = Vee(IState)+RMat(IState,iJ)*RMat(IState,iK)*RMat(IState,iL)*RMat(IState,iM)*DDG(iJ,iK,iL,iM)
        end do
      end do
    end do
  end do
  Vee(IState) = Vee(IState)/2
  !write(u6,'(A,I2,A,F10.6)') 'The classic coulomb energy for state ',IState,' is ',Vee(IState)
end do

end subroutine CalcVee

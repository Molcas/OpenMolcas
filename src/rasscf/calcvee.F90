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
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: Vee(lRoots), RMat(lRoots,lRoots), DDG(lRoots,lRoots,lRoots,lRoots)
integer(kind=iwp) :: iJ, iK, iL, iM, iState

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

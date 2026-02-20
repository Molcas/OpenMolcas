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
use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: Vee(lRoots)
real(kind=wp), intent(in) :: RMat(lRoots,lRoots), DDG(lRoots,lRoots,lRoots,lRoots)
integer(kind=iwp) :: iK, iL, iM, iState

Vee(:) = Zero
do IState=1,lRoots
  do iL=1,lRoots
    do iM=1,lRoots
      do iK=1,lRoots
        Vee(Istate) = Vee(IState)+RMat(IState,iL)*RMat(IState,iM)*RMat(IState,iK)*sum(RMat(IState,:)*DDG(:,iK,iL,iM))
      end do
    end do
  end do
  !write(u6,'(A,I2,A,F10.6)') 'The classic coulomb energy for state ',IState,' is ',Half*Vee(IState)
end do
Vee(:) = Half*Vee(:)

end subroutine CalcVee

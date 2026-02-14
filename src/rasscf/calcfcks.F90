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

subroutine CalcFckS(FckO,GDMat,FckS)

use rasscf_global, only: lRoots, nAc
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: FckO(NAC,NAC), GDMat(lRoots*(lRoots+1)/2,NAC,NAC), FckS(lRoots,lRoots)
integer(kind=iwp) :: iOrb, IState, jOrb, JState

FckS(:,:) = Zero

do IState=1,lRoots
  do JState=1,IState
    do IOrb=1,NAC
      do JOrb=1,NAC
        FckS(IState,JState) = FckS(IState,JState)+FckO(IOrb,JOrb)*GDMat(IState*(IState-1)/2+JState,IOrb,JOrb)
      end do
    end do
    FckS(JState,IState) = FckS(IState,JState)
  end do
end do

!call PrintMat('XMS_Mat','test',FckS,LRoots,LRoots,0,4,'N')

end subroutine CalcFckS

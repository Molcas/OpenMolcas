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

use Index_Functions, only: iTri, nTri_Elem
use rasscf_global, only: lRoots, nAc
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: FckO(NAC,NAC), GDMat(nTri_Elem(lRoots),NAC,NAC), FckS(lRoots,lRoots)
integer(kind=iwp) :: IState, JState

FckS(:,:) = Zero

do IState=1,lRoots
  do JState=1,IState
    FckS(IState,JState) = FckS(IState,JState)+sum(FckO(:,:)*GDMat(iTri(IState,JState),:,:))
    FckS(JState,IState) = FckS(IState,JState)
  end do
end do

!call PrintMat('XMS_Mat','test',FckS,LRoots,LRoots,0,4,'N')

end subroutine CalcFckS

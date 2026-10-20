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

subroutine NStateOpt2(RotMat,GDMat,Gtuvx)

use Index_Functions, only: nTri_Elem
use CMS, only: CMSNotConverged
use rasscf_global, only: CMSThreshold, iCMSIterMax, iCMSIterMin, lRoots, NAC
use PrintLevel, only: USUAL
use output_ras, only: IPRLOC
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, deg2rad
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: RotMat(lRoots,lRoots), GDMat(nTri_Elem(lRoots),NAC,NAC)
real(kind=wp), intent(in) :: Gtuvx(NAC,NAC,NAC,NAC)
integer(kind=iwp) :: ICMSIter, IPair, iPrLev, IState, JState, NPairs
real(kind=wp) :: VeeSumChange, VeeSumNew, VeeSumOld
logical(kind=iwp) :: Converged
integer(kind=iwp), allocatable :: StatePair(:,:)
real(kind=wp), allocatable :: FRot(:,:), theta(:), Vee(:)

IPRLEV = IPRLOC(6)

NPairs = nTri_Elem(lRoots-1)
call mma_allocate(StatePair,NPairs,2)
call mma_allocate(theta,NPairs)
call mma_allocate(Vee,lRoots)
call mma_allocate(FRot,lRoots,lRoots)
IPair = 0
do IState=1,lRoots
  do JState=1,IState-1
    IPair = IPair+1
    StatePair(IPair,1) = IState
    StatePair(IPair,2) = JState
  end do
end do
Converged = .false.
FRot(:,:) = RotMat(:,:)
call RotGDMat(FRot,GDMat)
call CalcVee2(Vee,GDMat,Gtuvx)
VeeSumOld = sum(Vee(:))
ICMSIter = 0
!write(u6,'(6X,I4,8X,F16.8,8X,ES16.4E3)') ICMSIter,VeeSumOld,0.0d0
do while (.not. Converged)
  theta(:) = Zero
  ICMSIter = ICMSIter+1
  call ThetaOpt2(FRot,theta,VeeSumChange,StatePair,NPairs,GDMat,Vee,Gtuvx)
  VeeSumNew = VeeSumOld+VeeSumChange
  if (IPRLEV >= USUAL) then
    if (lRoots > 2) then
      write(u6,'(6X,I4,8X,F16.8,8X,ES16.4E3)') ICMSIter,VeeSumNew,VeeSumChange
      !call RecPrt(' ',' ',Vee,lRoots,1)
      !write(u6,*) sum(Vee(:))
    else
      write(u6,'(6X,I4,8X,F6.1,9X,F16.8,5X,ES16.4E3)') ICMSIter,asin(FRot(2,1))/deg2rad,VeeSumNew,VeeSumChange
      !call RecPrt(' ',' ',Vee,lRoots,1)
      !write(u6,*) sum(Vee(:))
    end if
  end if
  if (abs(VeeSumChange) < CMSThreshold) then
    if (ICMSIter >= ICMSIterMin) then
      Converged = .true.
      if (IPRLEV >= USUAL) write(u6,'(4X,A)') 'CONVERGENCE REACHED'
    end if
  else
    if (ICMSIter >= ICMSIterMax) then
      Converged = .true.
      CMSNotConverged = .true.
      write(u6,'(4X,A)') 'NOT CONVERGED AFTER MAX NUMBER OF CYCLES'
      write(u6,'(4X,A)') 'TEMPORARY ROTATION MATRIX SAVED'
    end if
  end if
  !Converged = .true.
  VeeSumOld = VeeSumNew
end do
if (IPRLEV >= USUAL) write(u6,*) repeat('=',71)

RotMat(:,:) = FRot(:,:)
call mma_deallocate(StatePair)
call mma_deallocate(theta)
call mma_deallocate(Vee)
call mma_deallocate(FRot)

end subroutine NStateOpt2

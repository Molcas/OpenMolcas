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

use CMS, only: CMSNotConverged
use rasscf_global, only: lRoots, NAC, CMSThreshold, iCMSIterMax, iCMSIterMin
use PrintLevel, only: USUAL
use output_ras, only: IPRLOC
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
real*8, dimension(LRoots*(LRoots+1)/2,NAC,NAC) :: GDMat
real*8, dimension(lRoots,lRoots) :: RotMat
real*8, dimension(NAC,NAC,NAC,NAC) :: Gtuvx
integer IState, JState, NPairs, IPair, ICMSIter
real*8 VeeSumOld, VeeSumNew, Threshold, VeeSumChange
integer, dimension(:,:), allocatable :: StatePair
real*8, dimension(:), allocatable :: theta
real*8, dimension(:), allocatable :: Vee
real*8, dimension(:,:), allocatable :: FRot
logical Converged
real*8, external :: SumArray
integer iPrLev
#include "warnings.h"

IPRLEV = IPRLOC(6)

call mma_allocate(StatePair,LRoots*(LRoots-1)/2,2)
call mma_allocate(theta,LRoots*(LRoots-1)/2)
call mma_allocate(Vee,LRoots)
call mma_allocate(FRot,lRoots,lRoots)
Threshold = CMSThreshold
NPairs = lRoots*(lRoots-1)/2
IPair = 0
do IState=1,lRoots
  do JState=1,IState-1
    IPair = IPair+1
    StatePair(IPair,1) = IState
    StatePair(IPair,2) = JState
  end do
end do
Converged = .false.
call Copy2DMat(FRot,RotMat,lRoots,lRoots)
call RotGDMat(FRot,GDMat)
call CalcVee2(Vee,GDMat,Gtuvx)
VeeSumOld = SumArray(Vee,lRoots)
ICMSIter = 0
!write(6,'(6X,I4,8X,F16.8,8X,ES16.4E3)') ICMSIter,VeeSumOld,0.0d0
do while (.not. Converged)
  do IPair=1,NPairs
    theta(IPair) = 0.0d0
  end do
  ICMSIter = ICMSIter+1
  call ThetaOpt2(FRot,theta,VeeSumChange,StatePair,NPairs,GDMat,Vee,Gtuvx)
  VeeSumNew = VeeSumOld+VeeSumChange
  if (IPRLEV >= USUAL) then
    if (lRoots > 2) then
      write(6,'(6X,I4,8X,F16.8,8X,ES16.4E3)') ICMSIter,VeeSumNew,VeeSumChange
      !call RecPrt(' ',' ',Vee,lRoots,1)
      !write(6,*) SumArray(Vee,lRoots)
    else
      write(6,'(6X,I4,8X,F6.1,9X,F16.8,5X,ES16.4E3)') ICMSIter,asin(FRot(2,1))/atan(1.0d0)*45.0d0,VeeSumNew,VeeSumChange
      !call RecPrt(' ',' ',Vee,lRoots,1)
      !write(6,*) SumArray(Vee,lRoots)
    end if
  end if
  if (abs(VeeSumChange) < Threshold) then
    if (ICMSIter >= ICMSIterMin) then
      Converged = .true.
      if (IPRLEV >= USUAL) write(6,'(4X,A)') 'CONVERGENCE REACHED'
    end if
  else
    if (ICMSIter >= ICMSIterMax) then
      Converged = .true.
      CMSNotConverged = .true.
      write(6,'(4X,A)') 'NOT CONVERGED AFTER MAX NUMBER OF CYCLES'
      write(6,'(4X,A)') 'TEMPORARY ROTATION MATRIX SAVED'
    end if
  end if
  !Converged = .true.
  VeeSumOld = VeeSumNew
end do
if (IPRLEV >= USUAL) write(6,*) repeat('=',71)

call Copy2DMat(RotMat,FRot,lRoots,lRoots)
call mma_deallocate(StatePair)
call mma_deallocate(theta)
call mma_deallocate(Vee)
call mma_deallocate(FRot)

end subroutine NStateOpt2

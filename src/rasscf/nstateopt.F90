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

subroutine NStateOpt(RotMat,DDg)

use CMS, only: CMSNotConverged
use rasscf_global, only: lRoots, CMSThreshold, iCMSIterMax, iCMSIterMin
use PrintLevel, only: USUAL
use output_ras, only: IPRLOC
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, deg2rad
use Definitions, only: u6

implicit none
real*8, dimension(lRoots,lRoots,lRoots,lRoots) :: DDG
real*8, dimension(lroots,lroots) :: RotMat
integer IState, JState, NPairs, IPair, ICMSIter
real*8 VeeSumOld, VeeSumNew, Threshold
integer, dimension(:,:), allocatable :: StatePair
real*8, dimension(:), allocatable :: theta
real*8, dimension(:,:), allocatable :: FRot
logical Converged
real*8 CalcNSumVee
external CalcNSumVee
integer iPrLev
#include "warnings.h"

IPRLEV = IPRLOC(6)

call mma_allocate(StatePair,LRoots*(LRoots-1)/2,2)
call mma_allocate(theta,LRoots*(LRoots-1)/2)
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
VeeSumOld = CalcNSumVee(RotMat,DDg)
ICMSIter = 0
do while (.not. Converged)
  do IPair=1,NPairs
    theta(IPair) = Zero
  end do
  ICMSIter = ICMSIter+1
  call ThetaOpt(FRot,theta,VeeSumNew,StatePair,NPairs,DDg)
  if (IPRLEV >= USUAL) then
    if (lRoots > 2) then
      write(u6,'(6X,I4,8X,F16.8,8X,ES16.4E3)') ICMSIter,VeeSumNew,VeeSumNew-VeeSumOld
    else
      write(u6,'(6X,I4,8X,F6.1,9X,F16.8,5X,ES16.4E3)') ICMSIter,asin(FRot(2,1))/deg2rad,VeeSumNew,VeeSumNew-VeeSumOld
    end if
  end if
  if (abs(VeeSumNew-VeeSumOld) < Threshold) then
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
  VeeSumOld = VeeSumNew
end do
if (IPRLEV >= USUAL) write(u6,*) repeat('=',71)
call Copy2DMat(RotMat,FRot,lRoots,lRoots)
call mma_deallocate(StatePair)
call mma_deallocate(theta)
call mma_deallocate(FRot)

end subroutine NStateOpt

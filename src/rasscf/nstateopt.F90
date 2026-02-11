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
      Subroutine NStateOpt(RotMat,DDg)
      use stdalloc, only : mma_allocate, mma_deallocate
      use CMS, only: CMSNotConverged
      use rasscf_global, only: lRoots, CMSThreshold, iCMSIterMax,       &
     &                         iCMSIterMin
      use PrintLevel, only: USUAL
      use output_ras, only: LF,IPRLOC
      Implicit None


#include "warnings.h"
      Real*8,DIMENSION(lRoots,lRoots,lRoots,lRoots)::DDG
      Real*8,DIMENSION(lroots,lroots)::RotMat

      INTEGER IState,JState,NPairs,IPair,ICMSIter
      Real*8 VeeSumOld,VeeSumNew,Threshold
      INTEGER,DIMENSION(:,:),Allocatable::StatePair
      Real*8,DIMENSION(:),Allocatable::theta
      Real*8,DIMENSION(:,:),Allocatable::FRot
      Logical Converged
      Real*8 CalcNSumVee
      External CalcNSumVee
      Integer iPrLev

      IPRLEV=IPRLOC(6)

      CALL mma_allocate(StatePair,LRoots*(LRoots-1)/2,2)
      CALL mma_allocate(theta,LRoots*(LRoots-1)/2)
      CALL mma_allocate(FRot,lRoots,lRoots)
      Threshold=CMSThreshold
      NPairs=lRoots*(lRoots-1)/2
      IPair=0
      DO IState=1,lRoots
       Do JState=1,IState-1
        IPair=IPair+1
        StatePair(IPair,1)=IState
        StatePair(IPair,2)=JState
       End Do
      END DO
      Converged=.false.
      CALL Copy2DMat(FRot,RotMat,lRoots,lRoots)
      VeeSumOld=CalcNSumVee(RotMat,DDg)
      ICMSIter=0
      DO WHILE(.not.Converged)
       Do IPair=1,NPairs
        theta(IPair)=0.0d0
       End Do
       ICMSIter=ICMSIter+1
       CALL ThetaOpt(FRot,theta,VeeSumNew,StatePair,NPairs,DDg)
       IF(IPRLEV.ge.USUAL) THEN
       IF(lRoots.gt.2) THEN
       write(LF,'(6X,I4,8X,F16.8,8X,ES16.4E3)')                         &
     & ICMSIter,VeeSumNew,VeeSumNew-VeeSumOld
       ELSE
       write(LF,'(6X,I4,8X,F6.1,9X,F16.8,5X,ES16.4E3)')                 &
     & ICMSIter,asin(FRot(2,1))/atan(1.0d0)*45.0d0,VeeSumNew            &
     & ,VeeSumNew-VeeSumOld
       END IF
       END IF
       IF(ABS(VeeSumNew-VeeSumOld).lt.Threshold) THEN
        If(ICMSIter.ge.ICMSIterMin) Then
         Converged=.true.
         IF(IPRLEV.ge.USUAL) write(6,'(4X,A)')'CONVERGENCE REACHED'
        End If
       ELSE
        if(ICMSIter.ge.ICMSIterMax) then
         Converged=.true.
         CMSNotConverged=.true.
         write(LF,'(4X,A)')'NOT CONVERGED AFTER MAX NUMBER OF CYCLES'
         write(LF,'(4X,A)')'TEMPORARY ROTATION MATRIX SAVED'
        end if
       END IF
       VeeSumOld=VeeSumNew
      END DO
      IF(IPRLEV.ge.USUAL) write(6,*) repeat('=',71)
      CALL Copy2DMat(RotMat,FRot,lRoots,lRoots)
      CALL mma_deallocate(StatePair)
      CALL mma_deallocate(theta)
      CALL mma_deallocate(FRot)
      END SUBROUTINE NStateOpt

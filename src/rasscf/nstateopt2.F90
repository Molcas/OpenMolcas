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
      Subroutine NStateOpt2(RotMat,GDMat,Gtuvx)
      use stdalloc, only : mma_allocate, mma_deallocate
      use CMS, only: CMSNotConverged
      use rasscf_global, only: lRoots, NAC, CMSThreshold, iCMSIterMax,  &
     &                         iCMSIterMin
      use PrintLevel, only: USUAL
      use output_ras, only: IPRLOC
      Implicit None


#include "warnings.h"

      Real*8,DIMENSION(LRoots*(LRoots+1)/2,NAC,NAC)::GDMat
      Real*8,DIMENSION(lRoots,lRoots)::RotMat
      Real*8,DIMENSION(NAC,NAC,NAC,NAC)::Gtuvx

      INTEGER IState,JState,NPairs,IPair,ICMSIter
      Real*8 VeeSumOld,VeeSumNew,Threshold,VeeSumChange
      INTEGER,DIMENSION(:,:),Allocatable::StatePair
      Real*8,DIMENSION(:),Allocatable::theta
      Real*8,DIMENSION(:),Allocatable::Vee
      Real*8,DIMENSION(:,:),Allocatable::FRot
      Logical Converged
      Real*8, External :: SumArray
      Integer iPrLev

      IPRLEV=IPRLOC(6)

      CALL mma_allocate(StatePair,LRoots*(LRoots-1)/2,2)
      CALL mma_allocate(theta,LRoots*(LRoots-1)/2)
      CALL mma_allocate(Vee,LRoots)
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
      CALL RotGDMat(FRot,GDMat)
      CALL CalcVee2(Vee,GDMat,Gtuvx)
      VeeSumOld=SumArray(Vee,lRoots)
      ICMSIter=0
!        write(6,'(6X,I4,8X,F16.8,8X,ES16.4E3)')
!     &  ICMSIter,VeeSumOld,0.0d0
      DO WHILE(.not.Converged)
       Do IPair=1,NPairs
        theta(IPair)=0.0d0
       End Do
       ICMSIter=ICMSIter+1
       CALL ThetaOpt2                                                   &
     & (FRot,theta,VeeSumChange,StatePair,NPairs,GDMat,Vee,Gtuvx)
       VeeSumNew=VeeSumOld+VeeSumChange
       IF(IPRLEV.ge.USUAL) THEN
       IF(lRoots.gt.2) THEN
        write(6,'(6X,I4,8X,F16.8,8X,ES16.4E3)')                         &
     &  ICMSIter,VeeSumNew,VeeSumChange
!        CALL RecPrt(' ',' ',Vee,lRoots,1)
!        write(6,*) SumArray(Vee,lRoots)
       ELSE
       write(6,'(6X,I4,8X,F6.1,9X,F16.8,5X,ES16.4E3)')                  &
     & ICMSIter,asin(FRot(2,1))/atan(1.0d0)*45.0d0,VeeSumNew            &
     & ,VeeSumChange
!       CALL RecPrt(' ',' ',Vee,lRoots,1)
!       write(6,*) SumArray(Vee,lRoots)
       END IF
       END IF
       IF(ABS(VeeSumChange).lt.Threshold) THEN
        If(ICMSIter.ge.ICMSIterMin) Then
         Converged=.true.
         IF(IPRLEV.ge.USUAL) write(6,'(4X,A)')'CONVERGENCE REACHED'
        End If
       ELSE
        if(ICMSIter.ge.ICMSIterMax) then
         Converged=.true.
         CMSNotConverged=.true.
         write(6,'(4X,A)')'NOT CONVERGED AFTER MAX NUMBER OF CYCLES'
         write(6,'(4X,A)')'TEMPORARY ROTATION MATRIX SAVED'
        end if
       END IF
!         Converged=.true.
       VeeSumOld=VeeSumNew
      END DO
      IF(IPRLEV.ge.USUAL) write(6,*) repeat('=',71)

      CALL Copy2DMat(RotMat,FRot,lRoots,lRoots)
      CALL mma_deallocate(StatePair)
      CALL mma_deallocate(theta)
      CALL mma_deallocate(Vee)
      CALL mma_deallocate(FRot)
      END SUBROUTINE NStateOpt2

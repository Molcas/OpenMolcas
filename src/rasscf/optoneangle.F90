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
      SUBROUTINE OptOneAngle(Angle,SumVee,RotMat,DDg,I1,I2,lRoots)
      use stdalloc, only : mma_allocate, mma_deallocate
      Implicit None
      real*8 Angle,SumVee
      INTEGER I1,I2,lRoots
      Real*8,DIMENSION(lRoots,lRoots)::RotMat
      Real*8,DIMENSION(lRoots,lRoots,lRoots,lRoots)::DDG

      Logical Converged
      INTEGER Iter,IterMax,IA,IMax
      Real*8 Threshold,StepSize,SumOld
      Real*8,DIMENSION(:),Allocatable::Angles,Sums
      Real*8,DIMENSION(:),Allocatable::ScanA,ScanS
      Real*8,DIMENSION(:,:),Allocatable::RTmp

      Real*8, External ::CalcNSumVee
      INTEGER, External ::RMax

      CALL mma_allocate(Angles,4)
      CALL mma_allocate(Sums,4)
      CALL mma_allocate(ScanA,31)
      CALL mma_allocate(ScanS,31)
      CALL mma_allocate(RTmp,lRoots,lRoots)

      Converged=.false.
      stepsize=dble(atan(1.0d0))/15
      Threshold=1.0d-8

!       write(6,'(A,2(I2,2X))')
!     &'scanning rotation angles for ',I1,I2
      Angles(2)=0.0d0
      DO Iter=1,31
       ScanA(Iter)=(Iter-16)*stepsize*2
       CALL Copy2DMat(RTmp,RotMat,lRoots,lRoots)
       CALL CMSMatRot(RTmp,ScanA(Iter),I1,I2,lRoots)
       ScanS(Iter)=CalcNSumVee(RTmp,DDg)
!       IF(I2.eq.1) write(6,*) Iter,ScanA(Iter),ScanS(Iter)
      END DO

      IMax=RMax(ScanS,31)

      Iter=0
      IterMax=100
      SumOld=ScanS(IMax)
      Angles(2)=ScanA(IMax)
      DO WHILE(.not.Converged)
       Iter=Iter+1
       Angles(1)=Angles(2)-stepsize
       Angles(3)=Angles(2)+stepsize
       Do iA=1,3
        CALL Copy2DMat(RTmp,RotMat,lRoots,lRoots)
        CALL CMSMatRot(RTmp,Angles(iA),I1,I2,lRoots)
        Sums(iA)=CalcNSumVee(RTmp,DDg)
       End Do
       CALL CMSFitTrigonometric(Angles,Sums)
       CALL Copy2DMat(RTmp,RotMat,lRoots,lRoots)
       CALL CMSMatRot(RTmp,Angles(4),I1,I2,lRoots)
       Sums(4)=CalcNSumVee(RTmp,DDg)
       IF(ABS(Sums(4)-SumOld).lt.Threshold) THEN
        Converged=.true.
        Angle=Angles(4)
        CALL CMSMatRot(RotMat,Angle,I1,I2,lRoots)
        SumVee=CalcNSumVee(RotMat,DDg)
!        write(6,'(A,I3,A)')
!     &'Convergence reached after ',Iter,' micro cycles'
       ELSE
        If(Iter.eq.IterMax) Then
         Converged=.true.
        write(6,'(A,I3,A)')                                             &
     &'No convergence reached after ',Iter,' micro cycles'
        Else
         Angles(2)=Angles(4)
         SumOld=Sums(4)
        End If
       END IF
      END DO
      CALL mma_deallocate(Angles)
      CALL mma_deallocate(Sums)
      CALL mma_deallocate(ScanA)
      CALL mma_deallocate(ScanS)
      CALL mma_deallocate(RTmp)

      END SUBROUTINE OptOneAngle

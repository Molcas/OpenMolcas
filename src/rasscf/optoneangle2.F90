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
      SubRoutine OptOneAngle2(ang,change,R,GD,I1,I2,Vee,G)
      use stdalloc, only : mma_allocate, mma_deallocate
      use rasscf_global, only: lRoots, NAC
      Implicit None


#include "warnings.h"
      Real*8 ang,change
      Integer I1,I2
      Real*8,DIMENSION(lRoots,lRoots)::R
      Real*8,DIMENSION(LRoots*(LRoots+1)/2,NAC,NAC)::GD
      Real*8,DIMENSION(lRoots)::Vee
      Real*8,DIMENSION(NAC,NAC,NAC,NAC)::G

      Logical Converged
      INTEGER Itera,Itermax,IA,IMax
      Real*8 Threshold,StepSize,SumOld,Vee1,Vee2,SumOld2
      Real*8,DIMENSION(:),Allocatable::Angles,Sums
      Real*8,DIMENSION(:),Allocatable::ScanA,ScanS

      INTEGER, External :: RMax

      CALL mma_allocate(Angles,4)
      CALL mma_allocate(Sums,4)
      CALL mma_allocate(ScanA,31)
      CALL mma_allocate(ScanS,31)

      Converged=.false.
      stepsize=dble(atan(1.0d0))/15
      Threshold=1.0d-8

      Angles(2)=0.0d0
      DO Itera=1,31
       ScanA(Itera)=(Itera-16)*stepsize*2
       CALL                                                             &
     &SumVeeNew(ScanS(Itera),ScanA(Itera),GD,I1,I2,G,Vee1,Vee2,.false.)
!       IF(I2.eq.1) write(6,*) Iter,ScanA(Iter),ScanS(Iter)
      END DO

      IMax=RMax(ScanS,21)

      Itera=0
      IterMax=100
      SumOld=Vee(I1)+Vee(I2)
      SumOld2=SumOld
      Angles(2)=ScanA(IMax)
      DO WHILE(.not.Converged)
       Itera=Itera+1
       Angles(1)=Angles(2)-stepsize
       Angles(3)=Angles(2)+stepsize
       Do iA=1,3
       CALL SumVeeNew(Sums(iA),Angles(iA),GD,I1,I2,G,Vee1,Vee2,.false.)
       End Do
       CALL CMSFitTrigonometric(Angles,Sums)
       CALL SumVeeNew(Sums(4),Angles(4),GD,I1,I2,G,Vee1,Vee2,.false.)
       change=Sums(4)-SumOld
       IF(ABS(change).lt.Threshold) THEN
        Converged=.true.
        Ang=Angles(4)
        Vee(I1)=Vee1
        Vee(I2)=Vee2
        CALL SumVeeNew(Sums(4),Ang,GD,I1,I2,G,Vee1,Vee2,.true.)
        CALL CMSMatRot(R,Ang,I1,I2,lRoots)
       ELSE
        If(Itera.eq.IterMax) Then
         Converged=.true.
        write(6,'(A,I3,A)')                                             &
     &'No convergence reached after ',Itera,' micro cycles'
        Else
         Angles(2)=Angles(4)
         SumOld=Sums(4)
        End If
       END IF
      END DO
      change=Vee(I1)+Vee(I2)-SumOld2
      CALL mma_deallocate(Angles)
      CALL mma_deallocate(Sums)
      CALL mma_deallocate(ScanA)
      CALL mma_deallocate(ScanS)
      End Subroutine OptOneAngle2

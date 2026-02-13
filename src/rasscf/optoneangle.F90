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

subroutine OptOneAngle(Angle,SumVee,RotMat,DDg,I1,I2,lRoots)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Three, deg2rad
use Definitions, only: wp, u6

implicit none
real*8 Angle, SumVee
integer I1, I2, lRoots
real*8, dimension(lRoots,lRoots) :: RotMat
real*8, dimension(lRoots,lRoots,lRoots,lRoots) :: DDG
logical Converged
integer Iter, IterMax, IA, IMax
real*8 Threshold, StepSize, SumOld
real*8, dimension(:), allocatable :: Angles, Sums
real*8, dimension(:), allocatable :: ScanA, ScanS
real*8, dimension(:,:), allocatable :: RTmp
real*8, external :: CalcNSumVee
integer, external :: RMax

call mma_allocate(Angles,4)
call mma_allocate(Sums,4)
call mma_allocate(ScanA,31)
call mma_allocate(ScanS,31)
call mma_allocate(RTmp,lRoots,lRoots)

Converged = .false.
stepsize = Three*deg2rad
Threshold = 1.0e-8_wp

!write(u6,'(A,2(I2,2X))') 'scanning rotation angles for ',I1,I2
Angles(2) = Zero
do Iter=1,31
  ScanA(Iter) = (Iter-16)*stepsize*2
  call Copy2DMat(RTmp,RotMat,lRoots,lRoots)
  call CMSMatRot(RTmp,ScanA(Iter),I1,I2,lRoots)
  ScanS(Iter) = CalcNSumVee(RTmp,DDg)
  !if (I2 == 1) write(u6,*) Iter,ScanA(Iter),ScanS(Iter)
end do

IMax = RMax(ScanS,31)

Iter = 0
IterMax = 100
SumOld = ScanS(IMax)
Angles(2) = ScanA(IMax)
do while (.not. Converged)
  Iter = Iter+1
  Angles(1) = Angles(2)-stepsize
  Angles(3) = Angles(2)+stepsize
  do iA=1,3
    call Copy2DMat(RTmp,RotMat,lRoots,lRoots)
    call CMSMatRot(RTmp,Angles(iA),I1,I2,lRoots)
    Sums(iA) = CalcNSumVee(RTmp,DDg)
  end do
  call CMSFitTrigonometric(Angles,Sums)
  call Copy2DMat(RTmp,RotMat,lRoots,lRoots)
  call CMSMatRot(RTmp,Angles(4),I1,I2,lRoots)
  Sums(4) = CalcNSumVee(RTmp,DDg)
  if (abs(Sums(4)-SumOld) < Threshold) then
    Converged = .true.
    Angle = Angles(4)
    call CMSMatRot(RotMat,Angle,I1,I2,lRoots)
    SumVee = CalcNSumVee(RotMat,DDg)
    !write(u6,'(A,I3,A)') 'Convergence reached after ',Iter,' micro cycles'
  else
    if (Iter == IterMax) then
      Converged = .true.
      write(u6,'(A,I3,A)') 'No convergence reached after ',Iter,' micro cycles'
    else
      Angles(2) = Angles(4)
      SumOld = Sums(4)
    end if
  end if
end do
call mma_deallocate(Angles)
call mma_deallocate(Sums)
call mma_deallocate(ScanA)
call mma_deallocate(ScanS)
call mma_deallocate(RTmp)

end subroutine OptOneAngle

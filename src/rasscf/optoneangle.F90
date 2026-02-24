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
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(out) :: Angle, SumVee
integer(kind=iwp), intent(in) :: I1, I2, lRoots
real(kind=wp), intent(inout) :: RotMat(lRoots,lRoots)
real(kind=wp), intent(in) :: DDg(lRoots,lRoots,lRoots,lRoots)
integer(kind=iwp) :: IA, IMax, Iter, IterMax
real(kind=wp) :: StepSize, SumOld
logical(kind=iwp) :: Converged
real(kind=wp) :: Angles(4), ScanA(31), ScanS(31), Sums(4)
real(kind=wp), allocatable :: RTmp(:,:)
real(kind=wp), parameter :: Threshold = 1.0e-8_wp
real(kind=wp), external :: CalcNSumVee

call mma_allocate(RTmp,lRoots,lRoots)

Converged = .false.
stepsize = Three*deg2rad

!write(u6,'(A,2(I2,2X))') 'scanning rotation angles for ',I1,I2
Angles(2) = Zero
do Iter=1,31
  ScanA(Iter) = (Iter-16)*stepsize*2
  RTmp(:,:) = RotMat(:,:)
  call CMSMatRot(RTmp,ScanA(Iter),I1,I2,lRoots)
  ScanS(Iter) = CalcNSumVee(RTmp,DDg)
  !if (I2 == 1) write(u6,*) Iter,ScanA(Iter),ScanS(Iter)
end do

IMax = sum(maxloc(ScanS(:)))

Iter = 0
IterMax = 100
SumOld = ScanS(IMax)
Angles(2) = ScanA(IMax)
do while (.not. Converged)
  Iter = Iter+1
  Angles(1) = Angles(2)-stepsize
  Angles(3) = Angles(2)+stepsize
  do iA=1,3
    RTmp(:,:) = RotMat(:,:)
    call CMSMatRot(RTmp,Angles(iA),I1,I2,lRoots)
    Sums(iA) = CalcNSumVee(RTmp,DDg)
  end do
  call CMSFitTrigonometric(Angles,Sums)
  RTmp(:,:) = RotMat(:,:)
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
call mma_deallocate(RTmp)

end subroutine OptOneAngle

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

subroutine OptOneAngle2(ang,change,R,GD,I1,I2,Vee,G)

use rasscf_global, only: lRoots, NAC
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Three, deg2rad
use Definitions, only: wp, u6

implicit none
real*8 ang, change
integer I1, I2
real*8, dimension(lRoots,lRoots) :: R
real*8, dimension(LRoots*(LRoots+1)/2,NAC,NAC) :: GD
real*8, dimension(lRoots) :: Vee
real*8, dimension(NAC,NAC,NAC,NAC) :: G
logical Converged
integer Itera, Itermax, IA, IMax
real*8 Threshold, StepSize, SumOld, Vee1, Vee2, SumOld2
real*8, dimension(:), allocatable :: Angles, Sums
real*8, dimension(:), allocatable :: ScanA, ScanS
integer, external :: RMax
#include "warnings.h"

call mma_allocate(Angles,4)
call mma_allocate(Sums,4)
call mma_allocate(ScanA,31)
call mma_allocate(ScanS,31)

Converged = .false.
stepsize = Three*deg2rad
Threshold = 1.0e-8_wp

Angles(2) = Zero
do Itera=1,31
  ScanA(Itera) = (Itera-16)*stepsize*2
  call SumVeeNew(ScanS(Itera),ScanA(Itera),GD,I1,I2,G,Vee1,Vee2,.false.)
  !if (I2 == 1) write(u6,*) Iter,ScanA(Iter),ScanS(Iter)
end do

IMax = RMax(ScanS,21)

Itera = 0
IterMax = 100
SumOld = Vee(I1)+Vee(I2)
SumOld2 = SumOld
Angles(2) = ScanA(IMax)
do while (.not. Converged)
  Itera = Itera+1
  Angles(1) = Angles(2)-stepsize
  Angles(3) = Angles(2)+stepsize
  do iA=1,3
    call SumVeeNew(Sums(iA),Angles(iA),GD,I1,I2,G,Vee1,Vee2,.false.)
  end do
  call CMSFitTrigonometric(Angles,Sums)
  call SumVeeNew(Sums(4),Angles(4),GD,I1,I2,G,Vee1,Vee2,.false.)
  change = Sums(4)-SumOld
  if (abs(change) < Threshold) then
    Converged = .true.
    Ang = Angles(4)
    Vee(I1) = Vee1
    Vee(I2) = Vee2
    call SumVeeNew(Sums(4),Ang,GD,I1,I2,G,Vee1,Vee2,.true.)
    call CMSMatRot(R,Ang,I1,I2,lRoots)
  else
    if (Itera == IterMax) then
      Converged = .true.
      write(u6,'(A,I3,A)') 'No convergence reached after ',Itera,' micro cycles'
    else
      Angles(2) = Angles(4)
      SumOld = Sums(4)
    end if
  end if
end do
change = Vee(I1)+Vee(I2)-SumOld2
call mma_deallocate(Angles)
call mma_deallocate(Sums)
call mma_deallocate(ScanA)
call mma_deallocate(ScanS)

end subroutine OptOneAngle2

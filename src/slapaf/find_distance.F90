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
! Copyright (C) 2020, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine Find_Distance(Ref,Point,Dir,Fact,Dist,nAtom,BadConstraint)

use Slapaf_Info, only: RefGeo
use Slapaf_Parameters, only: MEP_Type

implicit none
#include "real.fh"
#include "stdalloc.fh"
integer, intent(In) :: nAtom
real*8, intent(In) :: Ref(3,nAtom), Dir(3,nAtom), Fact, Dist
real*8, intent(Out) :: Point(3,nAtom)
logical, intent(Out) :: BadConstraint
real*8, allocatable :: OldRef(:,:), Dummy(:), Not_Allocated(:,:)
real*8 :: R, CurFact, PrevR, Correct
real*8, parameter :: Thr = 1.0d-6
integer :: nCoor, i
real*8 rDum(1,1,1,1)
!                                                                      *
!***********************************************************************
!                                                                      *
interface
  subroutine SphInt(xyz,nCent,OfRef,RR0,Bf,l_Write,Label,dBf,ldB)
    integer nCent
    real*8 xyz(3,nCent)
    real*8, allocatable, target :: OfRef(:,:)
    real*8 RR0
    real*8 Bf(3,nCent)
    logical l_Write
    character(len=8) Label
    real*8 dBf(3,nCent,3,nCent)
    logical ldB
  end subroutine SphInt
end interface

!                                                                      *
!***********************************************************************
!                                                                      *
nCoor = 3*nAtom
call mma_allocate(OldRef,3,nAtom,Label='OldRef')
call mma_allocate(Dummy,nCoor,Label='Dummy')
OldRef(:,:) = RefGeo(:,:)
RefGeo(:,:) = Ref(:,:)

R = Zero
CurFact = Zero
Correct = Fact
i = 0
do while (abs(One-R/Dist) > Thr)

  ! Add the scaled direction vector
  CurFact = CurFact+Correct
  Point(:,:) = Ref(:,:)+CurFact*Dir(:,:)

  ! Align and measure distance
  PrevR = R
  call Align(Point(:,:),Ref(:,:),nAtom)
  if (MEP_Type == 'SPHERE') then
    call SphInt(Point,nAtom,Not_Allocated,R,Dummy,.false.,'dummy   ',rDum,.false.)
  else if (MEP_Type == 'TRANSVERSE') then
    call Transverse(Point,nAtom,R,Dummy,.false.,'dummy   ',rDum,.false.)
  end if

  ! Stop if too many iterations or if the constraint is moving
  ! in the wrong direction
  i = i+1
  if ((i > 5) .or. (Correct*(R-PrevR) < Zero)) exit
  Correct = (One-R/Dist)*Fact
end do
BadConstraint = (abs(One-R/Dist) > Thr)

RefGeo(:,:) = OldRef(:,:)
call mma_deallocate(OldRef)
call mma_deallocate(Dummy)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Find_Distance

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

use Slapaf_Info, only: MEP_Type, RefGeo
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtom
real(kind=wp), intent(in) :: Ref(3,nAtom), Dir(3,nAtom), Fact, Dist
real(kind=wp), intent(out) :: Point(3,nAtom)
logical(kind=iwp), intent(out) :: BadConstraint
integer(kind=iwp) :: i, nCoor
real(kind=wp) :: Correct, CurFact, PrevR, R, rDum(1,1,1,1)
real(kind=wp), allocatable :: Dummy(:), OldRef(:,:)
real(kind=wp), parameter :: Thr = 1.0e-6_wp

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
    call SphInt(Point,nAtom,Point,.false.,R,Dummy,.false.,'dummy   ',rDum,.false.)
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

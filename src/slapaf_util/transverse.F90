!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Transverse(xyz,nCent,HDist,Bf,l_Write,Label,dBf,ldB)

use Slapaf_Info, only: GradRef, R12, RefGeo, Weights
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nCent
real(kind=wp), intent(in) :: xyz(3,nCent)
real(kind=wp), intent(out) :: HDist, Bf(3,nCent)
logical(kind=iwp), intent(in) :: l_Write, ldB
character(len=8), intent(in) :: Label
real(kind=wp), intent(inout) :: dBf(3,nCent,3,nCent)
integer(kind=iwp) :: i, iCent, nTrans
real(kind=wp) :: f, RR_R12, SqInvTWeight, TWeight, xWeight
logical(kind=iwp) :: lTrans
real(kind=wp), allocatable, target :: TV(:,:)
real(kind=wp), pointer :: r12_p(:,:)
integer(kind=iwp), external :: iDeg

!                                                                      *
!***********************************************************************
!                                                                      *
! The reference direction (normal to the hyperplane) is taken from
! the input (GradRef) if allocated, or from the stored
! transverse direction if there is one, or from the R-P vector
! (R12) otherwise

if (allocated(GradRef)) then
  lTrans = .false.
  r12_p => GradRef
  !write(u6,*), 'Using Reference Gradient'
else
  call qpg_dArray('Transverse',lTrans,nTrans)
  if (lTrans) then
    call mma_allocate(TV,3,nCent,Label='TV')
    call Get_dArray('Transverse',TV,3*nCent)
    r12_p => TV
    !write(u6,*), 'Using stored Transverse'
  else
    r12_p => R12
    !write(u6,*), 'Using R-P Vector'
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
!call RecPrt('Ref',' ',RefGeo,3,nCent)
!call RecPrt('R12',' ',r12_p,3,nCent)

! Length of the direction vector (weighted)

RR_R12 = Zero
TWeight = Zero
do iCent=1,nCent
  xWeight = real(iDeg(xyz(1,iCent)),kind=wp)*Weights(iCent)
  TWeight = TWeight+xWeight
  do i=1,3
    RR_R12 = RR_R12+xWeight*(r12_p(i,iCent))**2
  end do
end do
RR_R12 = sqrt(RR_R12)

! The distance scaling will be 1/Sqrt(TWeight)

SqInvTWeight = One/sqrt(TWeight)

! Dot product between the x-x0 vector and the direction (weighted)

f = Zero
do iCent=1,nCent
  xWeight = real(iDeg(xyz(1,iCent)),kind=wp)*Weights(iCent)
  do i=1,3
    f = f+xWeight*(xyz(i,iCent)-RefGeo(i,iCent))*r12_p(i,iCent)
  end do
end do

! The distance to the plane is the dot product / direction length

if (RR_R12 == Zero) then
  HDist = Zero
else
  HDist = f/RR_R12*SqInvTWeight
end if
!write(u6,*) 'f, RR_R12=',f,RR_R12

if (l_Write) write(u6,'(2A,F18.8,A)') Label,' : Hyperplane distance =',HDist,' au (weighted/sqrt(total weight)'

!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the WDC B-matrix
! If the direction is null, the derivative is not defined

Bf(:,:) = Zero
if (RR_R12 > Zero) then

  ! The derivative is simply the unit direction vector
  ! (with weighting and scaling accounted for)

  do iCent=1,nCent
    xWeight = real(iDeg(xyz(1,iCent)),kind=wp)*Weights(iCent)
    do i=1,3
      Bf(:,iCent) = xWeight*r12_p(:,iCent)/RR_R12*SqInvTWeight
    end do
  end do
end if
!call RecPrt('Bf',' ',Bf,3,nCent)
!                                                                      *
!***********************************************************************
!                                                                      *
! The second derivative is null, as the derivative is constant

if (ldB) dBf(:,:,:,:) = Zero

if (lTrans) call mma_deallocate(TV)
nullify(r12_p)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Transverse

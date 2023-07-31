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

subroutine SphInt(xyz,nCent,OfRef,lOR,RR0,Bf,l_Write,Label,dBf,ldB)

use Slapaf_Info, only: RefGeo, Weights
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nCent
real(kind=wp), intent(in) :: xyz(3,nCent), OfRef(3,nCent)
real(kind=wp), intent(out) :: RR0, Bf(3,nCent)
logical(kind=iwp), intent(in) :: lOR, l_Write, ldB
character(len=8), intent(in) :: Label
real(kind=wp), intent(inout) :: dBf(3,nCent,3,nCent)
integer(kind=iwp) :: iCar, iCent, ixyz, jCent, jxyz
real(kind=wp) :: Fact, RR0_unscaled, SqInvTWeight, temp, Tweight, xWeight, yWeight
real(kind=wp), allocatable :: xyz_0(:,:)
integer(kind=iwp), external :: iDeg

!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the radius of the hypersphere

call mma_allocate(xyz_0,3,nCent,Label='xyz_0')
!call RecPrt('SphInt: xyz',' ',xyz,3,nCent)
if (lOR) then
  xyz_0(:,:) = xyz(:,:)-OfRef(:,:)
  !call RecPrt('Ref:',' ',OfRef,3,nCent)
else
  xyz_0(:,:) = xyz(:,:)-RefGeo(:,:)
  !call RecPrt('Ref:',' ',RefGeo,3,nCent)
end if
RR0 = Zero
TWeight = Zero
do iCent=1,nCent
  Fact = real(iDeg(xyz(1,iCent)),kind=wp)
  xWeight = Fact*Weights(iCent)
  TWeight = TWeight+xWeight
  !write(u6,*) 'xWeight=',xWeight
  do ixyz=1,3
    !write(u6,*) xyz_0(ixyz,iCent)
    RR0 = RR0+xWeight*xyz_0(ixyz,iCent)**2
  end do
end do
RR0_unscaled = sqrt(RR0)

! RR0_unscaled is the real (weighthed) distance,
! RR0 is scaled by 1/Sqrt(TWeight) and so are the derivatives

SqInvTWeight = One/sqrt(TWeight)
RR0 = RR0_unscaled*SqInvTWeight

if (l_Write) write(u6,'(2A,F18.8,A)') Label,' : Radius of h-sphere= ',RR0,' au (weighted/sqrt(total weight))'
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the WDC B-matrix

!FIXME: revise the symmetry
do iCent=1,nCent
  Fact = real(iDeg(xyz(1,iCent)),kind=wp)
  xWeight = Fact*Weights(iCent)
  do iCar=1,3
    if (RR0_unscaled /= Zero) then
      Bf(iCar,iCent) = xWeight*xyz_0(iCar,iCent)/RR0_unscaled*SqInvTWeight
    else

      ! If we are standing on the reference point the gradient
      ! is not well defined.

      Bf(iCar,iCent) = Zero
    end if
  end do
end do
!call RecPrt('Bf',' ',Bf,3,nCent)
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the cartesian derivative of the B-Matrix.

!FIXME: revise the symmetry
if (ldB) then
  dBf(:,:,:,:) = Zero
  if (RR0 /= Zero) then
    do iCent=1,nCent
      Fact = real(iDeg(xyz(1,iCent)),kind=wp)
      xWeight = Fact*Weights(iCent)
      do ixyz=1,3
        do jCent=1,nCent
          Fact = real(iDeg(xyz(1,jCent)),kind=wp)
          yWeight = Fact*Weights(jCent)
          do jxyz=1,3
            temp = Zero
            if ((ixyz == jxyz) .and. (iCent == jCent)) temp = RR0_unscaled
            temp = temp-yWeight*xyz_0(ixyz,iCent)*xyz_0(jxyz,jCent)/RR0_unscaled
            temp = (xWeight*temp)/RR0_unscaled**2
            dBf(ixyz,iCent,jxyz,jCent) = temp*SqInvTWeight
          end do
        end do
      end do
    end do
  end if
  !call RecPrt('dBf',' ',dBf,3*nCent,3*nCent)

end if

call mma_deallocate(xyz_0)

return

end subroutine SphInt

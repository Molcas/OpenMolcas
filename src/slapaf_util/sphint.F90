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

subroutine SphInt(xyz,nCent,OfRef,RR0,Bf,l_Write,Label,dBf,ldB)

use Slapaf_Info, only: Weights, RefGeo

implicit real*8(a-h,o-z)
#include "real.fh"
real*8 Bf(3,nCent), xyz(3,nCent), dBf(3,nCent,3,nCent)
real*8, allocatable, target :: OfRef(:,:)
logical l_Write, ldB
character*8 Label
real*8, pointer :: xyz0(:,:)

!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the radius of the hypersphere

if (.not. allocated(OfRef)) then
  xyz0 => RefGeo(1:3,1:nCent)
else
  xyz0 => OfRef(1:3,1:nCent)
end if
!call RecPrt('SphInt: xyz',' ',xyz,3,nCent)
!call RecPrt('Ref: xyz0',' ',OfRef,3,nCent)
RR0 = Zero
TWeight = Zero
do iCent=1,nCent
  Fact = dble(iDeg(xyz(1,iCent)))
  xWeight = Fact*Weights(iCent)
  TWeight = TWeight+xWeight
  !write(6,*) 'xWeight=',xWeight
  do ixyz=1,3
    !write(6,*) xyz(ixyz,iCent),xyz0(ixyz,iCent)
    temp = xyz(ixyz,iCent)-xyz0(ixyz,iCent)
    RR0 = RR0+xWeight*temp**2
  end do
end do
RR0_unscaled = sqrt(RR0)

! RR0_unscaled is the real (weighthed) distance,
! RR0 is scaled by 1/Sqrt(TWeight) and so are the derivatives

SqInvTWeight = One/sqrt(TWeight)
RR0 = RR0_unscaled*SqInvTWeight

if (l_Write) then
  write(6,'(2A,F18.8,A)') Label,' : Radius of h-sphere= ',RR0,' au (weighted/sqrt(total weight))'
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the WDC B-matrix

!FIXME: revise the symmetry
do iCent=1,nCent
  Fact = dble(iDeg(xyz(1,iCent)))
  xWeight = Fact*Weights(iCent)
  do iCar=1,3
    temp = xyz(iCar,iCent)-xyz0(iCar,iCent)
    if (RR0_unscaled /= Zero) then
      temp = xyz(iCar,iCent)-xyz0(iCar,iCent)
      Bf(iCar,iCent) = xWeight*temp/RR0_unscaled*SqInvTWeight
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
  call FZero(dBf,(3*nCent)**2)
  if (RR0 == Zero) Go To 99
  do iCent=1,nCent
    Fact = dble(iDeg(xyz(1,iCent)))
    xWeight = Fact*Weights(iCent)
    do ixyz=1,3
      tempi = xyz(ixyz,iCent)-xyz0(ixyz,iCent)
      do jCent=1,nCent
        Fact = dble(iDeg(xyz(1,jCent)))
        yWeight = Fact*Weights(jCent)
        do jxyz=1,3
          tempj = xyz(jxyz,jCent)-xyz0(jxyz,jCent)
          temp = Zero
          if ((ixyz == jxyz) .and. (iCent == jCent)) temp = RR0_unscaled
          temp = temp-yWeight*tempi*tempj/RR0_unscaled
          temp = (xWeight*temp)/RR0_unscaled**2
          dBf(ixyz,iCent,jxyz,jCent) = temp*SqInvTWeight
        end do
      end do
    end do
  end do
99 continue
!call RecPrt('dBf',' ',dBf,3*nCent,3*nCent)

end if

xyz0 => null()

return

end subroutine SphInt

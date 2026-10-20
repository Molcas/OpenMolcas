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

subroutine SumVeeNew(SV,A,GD,I1,I2,G,V1,V2,Update)

use Index_Functions, only: iTri, nTri_Elem
use rasscf_global, only: lRoots, NAC
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: SV, GD(nTri_Elem(lRoots),NAC,NAC), V1, V2
real(kind=wp), intent(in) :: A, G(NAC,NAC,NAC,NAC)
integer(kind=iwp), intent(in) :: I1, I2
logical(kind=iwp), intent(in) :: Update
integer(kind=iwp) :: i11, i12, I1J, i22, I2J, J, u, v, x
real(kind=wp), allocatable :: D11(:,:), D1J(:,:,:), D22(:,:), D2J(:,:,:)

i11 = nTri_Elem(I1)
i22 = nTri_Elem(I2)
i12 = iTri(I1,I2)
if (Update) then
  call mma_allocate(D1J,lRoots,NAC,NAC)
  call mma_allocate(D2J,lRoots,NAC,NAC)
  ! calculating
  do J=1,I2-1      !(J<I2<I1)
    I1J = iTri(I1,J)
    I2J = iTri(I2,J)
    D2J(J,:,:) = cos(A)*GD(I2J,:,:)+sin(A)*GD(I1J,:,:)
    D1J(J,:,:) = -sin(A)*GD(I2J,:,:)+cos(A)*GD(I1J,:,:)
  end do
  J = I2           !(J=I2<I1)
  do u=1,NAC
    D2J(i2,:,u) = GD(i11,:,u)*sin(A)**2+GD(i22,:,u)*cos(A)**2+cos(A)*sin(A)*(GD(i12,u,:)+GD(i12,:,u))
    D1J(i2,:,u) = cos(A)*sin(A)*(GD(i11,:,u)-GD(i22,:,u))+GD(i12,:,u)*cos(A)**2-GD(i12,u,:)*sin(A)**2
  end do
  do J=I2+1,I1-1   !(I2<J<I1)
    I1J = iTri(I1,J)
    I2J = iTri(J,I2)
    do u=1,NAC
      D2J(J,:,u) = cos(A)*GD(I2J,u,:)+sin(A)*GD(I1J,:,u)
      D1J(J,:,u) = -sin(A)*GD(I2J,u,:)+cos(A)*GD(I1J,:,u)
    end do
  end do
  J = I1           !(I2<J=I1)
  do u=1,NAC
    D1J(I1,:,u) = GD(i11,:,u)*cos(A)**2+GD(i22,:,u)*sin(A)**2-cos(A)*sin(A)*(GD(i12,:,u)+GD(i12,u,:))
  end do

  do J=I1+1,lRoots !(I2<I1<J)
    I1J = iTri(J,I1)
    I2J = iTri(J,I2)
    do u=1,NAC
      D2J(J,:,u) = cos(A)*GD(I2J,u,:)+sin(A)*GD(I1J,u,:)
      D1J(J,:,u) = -sin(A)*GD(I2J,u,:)+cos(A)*GD(I1J,u,:)
    end do
  end do
  ! updating
  do J=1,I2-1      !(J<I2<I1)
    I1J = iTri(I1,J)
    I2J = iTri(I2,J)
    GD(I2J,:,:) = D2J(J,:,:)
    GD(I1J,:,:) = D1J(J,:,:)
  end do
  J = I2           !(J=I2<I1)
  GD(I22,:,:) = D2J(I2,:,:)
  GD(I12,:,:) = D1J(I2,:,:)
  do J=I2+1,I1-1   !(I2<J<I1)
    I1J = iTri(I1,J)
    I2J = iTri(J,I2)
    do u=1,NAC
      GD(I2J,:,u) = D2J(J,u,:)
      GD(I1J,:,u) = D1J(J,:,u)
    end do
  end do
  J = I1           !(I2<J=I1)
  GD(i11,:,:) = D1J(I1,:,:)

  do J=I1+1,lRoots !(I2<I1<J)
    I1J = iTri(J,I1)
    I2J = iTri(J,I2)
    do u=1,NAC
      GD(I2J,:,u) = D2J(J,u,:)
      GD(I1J,:,u) = D1J(J,u,:)
    end do
  end do
  call mma_deallocate(D1J)
  call mma_deallocate(D2J)
else
  call mma_allocate(D11,NAC,NAC)
  call mma_allocate(D22,NAC,NAC)
  do u=1,NAC
    D11(:,u) = GD(i11,:,u)*cos(A)**2+GD(i22,:,u)*sin(A)**2-cos(A)*sin(A)*(GD(i12,:,u)+GD(i12,u,:))
    D22(:,u) = GD(i11,:,u)*sin(A)**2+GD(i22,:,u)*cos(A)**2+cos(A)*sin(A)*(GD(i12,u,:)+GD(i12,:,u))
  end do
  V1 = Zero
  V2 = Zero
  do v=1,NAC
    do x=1,NAC
      V1 = V1+sum(D11(:,:)*D11(v,x)*G(:,:,v,x))
      V2 = V2+sum(D22(:,:)*D22(v,x)*G(:,:,v,x))
    end do
  end do
  V1 = Half*V1
  V2 = Half*V2
  SV = V1+V2
  call mma_deallocate(D11)
  call mma_deallocate(D22)
end if

end subroutine SumVeeNew

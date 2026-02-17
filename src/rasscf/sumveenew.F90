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
real(kind=wp) :: SV, A, GD(nTri_Elem(lRoots),NAC,NAC), G(NAC,NAC,NAC,NAC), V1, V2
integer(kind=iwp) :: I1, I2
logical(kind=iwp) :: Update
integer(kind=iwp) :: i11, i12, I1J, i22, I2J, J, t, u, v, x
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
    do t=1,NAC
      do u=1,NAC
        D2J(J,t,u) = cos(A)*GD(I2J,t,u)+sin(A)*GD(I1J,t,u)
        D1J(J,t,u) = -sin(A)*GD(I2J,t,u)+cos(A)*GD(I1J,t,u)
      end do
    end do
  end do
  J = I2           !(J=I2<I1)
  do t=1,NAC
    do u=1,NAC
      D2J(i2,t,u) = GD(i11,t,u)*sin(A)**2+GD(i22,t,u)*cos(A)**2+cos(A)*sin(A)*(GD(i12,u,t)+GD(i12,t,u))
      D1J(i2,t,u) = cos(A)*sin(A)*(GD(i11,t,u)-GD(i22,t,u))+GD(i12,t,u)*cos(A)**2-GD(i12,u,t)*sin(A)**2
    end do
  end do
  do J=I2+1,I1-1   !(I2<J<I1)
    I1J = iTri(I1,J)
    I2J = iTri(J,I2)
    do t=1,NAC
      do u=1,NAC
        D2J(J,t,u) = cos(A)*GD(I2J,u,t)+sin(A)*GD(I1J,t,u)
        D1J(J,t,u) = -sin(A)*GD(I2J,u,t)+cos(A)*GD(I1J,t,u)
      end do
    end do
  end do
  J = I1           !(I2<J=I1)
  do t=1,NAC
    do u=1,NAC
      D1J(i1,t,u) = GD(i11,t,u)*cos(A)**2+GD(i22,t,u)*sin(A)**2-cos(A)*sin(A)*(GD(i12,t,u)+GD(i12,u,t))
    end do
  end do

  do J=I1+1,lRoots !(I2<I1<J)
    I1J = iTri(J,I1)
    I2J = iTri(J,I2)
    do t=1,NAC
      do u=1,NAC
        D2J(J,t,u) = cos(A)*GD(I2J,u,t)+sin(A)*GD(I1J,u,t)
        D1J(J,t,u) = -sin(A)*GD(I2J,u,t)+cos(A)*GD(I1J,u,t)
      end do
    end do
  end do
  ! updating
  do J=1,I2-1      !(J<I2<I1)
    I1J = iTri(I1,J)
    I2J = iTri(I2,J)
    do t=1,NAC
      do u=1,NAC
        GD(I2J,t,u) = D2J(J,t,u)
        GD(I1J,t,u) = D1J(J,t,u)
      end do
    end do
  end do
  J = I2           !(J=I2<I1)
  do t=1,NAC
    do u=1,NAC
      GD(I22,t,u) = D2J(I2,t,u)
      GD(I12,t,u) = D1J(I2,t,u)
    end do
  end do
  do J=I2+1,I1-1   !(I2<J<I1)
    I1J = iTri(I1,J)
    I2J = iTri(J,I2)
    do t=1,NAC
      do u=1,NAC
        GD(I2J,t,u) = D2J(J,u,t)
        GD(I1J,t,u) = D1J(J,t,u)
      end do
    end do
  end do
  J = I1           !(I2<J=I1)
  do t=1,NAC
    do u=1,NAC
      GD(i11,t,u) = D1J(I1,t,u)
    end do
  end do

  do J=I1+1,lRoots !(I2<I1<J)
    I1J = iTri(J,I1)
    I2J = iTri(J,I2)
    do t=1,NAC
      do u=1,NAC
        GD(I2J,t,u) = D2J(J,u,t)
        GD(I1J,t,u) = D1J(J,u,t)
      end do
    end do
  end do
  call mma_deallocate(D1J)
  call mma_deallocate(D2J)
else
  call mma_allocate(D11,NAC,NAC)
  call mma_allocate(D22,NAC,NAC)
  V1 = Zero
  V2 = V1
  do t=1,NAC
    do u=1,NAC
      D11(t,u) = GD(i11,t,u)*cos(A)**2+GD(i22,t,u)*sin(A)**2-cos(A)*sin(A)*(GD(i12,t,u)+GD(i12,u,t))
      D22(t,u) = GD(i11,t,u)*sin(A)**2+GD(i22,t,u)*cos(A)**2+cos(A)*sin(A)*(GD(i12,u,t)+GD(i12,t,u))
    end do
  end do
  do t=1,NAC
    do u=1,NAC
      do v=1,NAC
        do x=1,NAC
          V1 = V1+D11(t,u)*D11(v,x)*G(t,u,v,x)
          V2 = V2+D22(t,u)*D22(v,x)*G(t,u,v,x)
        end do
      end do
    end do
  end do
  V1 = Half*V1
  V2 = Half*V2
  SV = V1+V2
  call mma_deallocate(D11)
  call mma_deallocate(D22)
end if

end subroutine SumVeeNew

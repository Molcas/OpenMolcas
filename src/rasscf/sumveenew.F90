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

use rasscf_global, only: lRoots, NAC
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
real*8 SV, A, V1, V2
integer I1, I2
real*8, dimension(LRoots*(LRoots+1)/2,NAC,NAC) :: GD
real*8, dimension(NAC,NAC,NAC,NAC) :: G
real*8, dimension(:,:), allocatable :: D11, D22
real*8, dimension(:,:,:), allocatable :: D1J, D2J
logical Update
integer t, u, v, x, i11, i22, i12
integer J, I1J, I2J
#include "warnings.h"

if (Update) then
  call mma_allocate(D1J,lRoots,NAC,NAC)
  call mma_allocate(D2J,lRoots,NAC,NAC)
  ! calculating
  do J=1,I2-1      !(J<I2<I1)
    I1J = (I1-1)*I1/2+J
    I2J = (I2-1)*I2/2+J
    do t=1,NAC
      do u=1,NAC
        D2J(J,t,u) = cos(A)*GD(I2J,t,u)+sin(A)*GD(I1J,t,u)
        D1J(J,t,u) = -sin(A)*GD(I2J,t,u)+cos(A)*GD(I1J,t,u)
      end do
    end do
  end do
  J = I2           !(J=I2<I1)
  i11 = (I1+1)*I1/2
  i22 = (I2+1)*I2/2
  i12 = (I1-1)*I1/2+I2
  do t=1,NAC
    do u=1,NAC
      D2J(i2,t,u) = GD(i11,t,u)*sin(A)**2+GD(i22,t,u)*cos(A)**2+cos(A)*sin(A)*(GD(i12,u,t)+GD(i12,t,u))
      D1J(i2,t,u) = cos(A)*sin(A)*(GD(i11,t,u)-GD(i22,t,u))+GD(i12,t,u)*cos(A)**2-GD(i12,u,t)*sin(A)**2
    end do
  end do
  do J=I2+1,I1-1   !(I2<J<I1)
    I1J = (I1-1)*I1/2+J
    I2J = (J-1)*J/2+I2
    do t=1,NAC
      do u=1,NAC
        D2J(J,t,u) = cos(A)*GD(I2J,u,t)+sin(A)*GD(I1J,t,u)
        D1J(J,t,u) = -sin(A)*GD(I2J,u,t)+cos(A)*GD(I1J,t,u)
      end do
    end do
  end do
  J = I1           !(I2<J=I1)
  i11 = (I1+1)*I1/2
  i22 = (I2+1)*I2/2
  i12 = (I1-1)*I1/2+I2
  do t=1,NAC
    do u=1,NAC
      D1J(i1,t,u) = GD(i11,t,u)*cos(A)**2+GD(i22,t,u)*sin(A)**2-cos(A)*sin(A)*(GD(i12,t,u)+GD(i12,u,t))
    end do
  end do

  do J=I1+1,lRoots !(I2<I1<J)
    I1J = (J-1)*J/2+I1
    I2J = (J-1)*J/2+I2
    do t=1,NAC
      do u=1,NAC
        D2J(J,t,u) = cos(A)*GD(I2J,u,t)+sin(A)*GD(I1J,u,t)
        D1J(J,t,u) = -sin(A)*GD(I2J,u,t)+cos(A)*GD(I1J,u,t)
      end do
    end do
  end do
!      updating
  do J=1,I2-1      !(J<I2<I1)
    I1J = (I1-1)*I1/2+J
    I2J = (I2-1)*I2/2+J
    do t=1,NAC
      do u=1,NAC
        GD(I2J,t,u) = D2J(J,t,u)
        GD(I1J,t,u) = D1J(J,t,u)
      end do
    end do
  end do
  J = I2           !(J=I2<I1)
!      i11=(I1+1)*I1/2
  i22 = (I2+1)*I2/2
  i12 = (I1-1)*I1/2+I2
  do t=1,NAC
    do u=1,NAC
      GD(I22,t,u) = D2J(I2,t,u)
      GD(I12,t,u) = D1J(I2,t,u)
    end do
  end do
  do J=I2+1,I1-1   !(I2<J<I1)
    I1J = (I1-1)*I1/2+J
    I2J = (J-1)*J/2+I2
    do t=1,NAC
      do u=1,NAC
        GD(I2J,t,u) = D2J(J,u,t)
        GD(I1J,t,u) = D1J(J,t,u)
      end do
    end do
  end do
  J = I1           !(I2<J=I1)
  i11 = (I1+1)*I1/2
  do t=1,NAC
    do u=1,NAC
      GD(i11,t,u) = D1J(I1,t,u)
    end do
  end do

  do J=I1+1,lRoots !(I2<I1<J)
    I1J = (J-1)*J/2+I1
    I2J = (J-1)*J/2+I2
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
  i11 = (I1+1)*I1/2
  i22 = (I2+1)*I2/2
  i12 = (I1-1)*I1/2+I2
  V1 = 0.0d0
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
  V1 = V1/2.0d0
  V2 = V2/2.0d0
  SV = V1+V2
  call mma_deallocate(D11)
  call mma_deallocate(D22)
end if

end subroutine SumVeeNew

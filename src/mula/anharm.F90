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
! Copyright (C) 1996, Niclas Forsberg                                  *
!***********************************************************************

subroutine Anharm(eigenVec,harmfreq,D3,D4,Gprime,Gdbleprime,x,nOsc)
!  Purpose:
!    Calculate the anharmonicity constants.
!    This routine assumes that the curvilinear coordinates are such
!    that they coincide with dimensionless normal coordinates for
!    small displacements.
!
!  Input:
!    eigenVec   : Real*8 two dimensional array - eigenvectors.
!    harmfreq   : Real*8 array - harmonic frequencies.
!    D3         : Real*8 three dimensional array - cubic force
!                 constants.
!    D4         : Real*8 four dimensional array - quartic force
!                 constants.
!    Gprime     : Real*8 three dimensional array - first
!                 derivatives of the inverse mass tensor.
!    Gdbleprime : Real*8 four dimensional array - second
!                 derivatives of the inverse mass tensor.
!
!  Output:
!    x          : Real*8 two dimensional array - anharmonicity
!                 constants.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

use mula_global, only: ngdim
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Three, Four, Five, Six, Eight, Nine, Half, Quart
use Definitions, only: wp

implicit real*8(a-h,o-z)
real*8 eigenVec(nosc,nosc)
real*8 harmfreq(nosc)
real*8 D3(ngdim,ngdim,ngdim)
real*8 D4(ngdim,ngdim,ngdim,ngdim)
real*8 Gprime(ngdim,ngdim,ngdim)
real*8 Gdbleprime(ngdim,ngdim,ngdim,ngdim)
real*8 x(nosc,nosc)
real*8, allocatable :: C(:,:), T3(:,:,:), T4(:,:,:,:), Temp(:,:), V3(:,:,:), V4(:,:,:,:)

NumInt = nOsc

call mma_allocate(C,nOsc,nOsc,label='C')
call mma_allocate(Temp,nOsc,nOsc,label='Temp')
call mma_allocate(V3,nOsc,nOsc,nOsc,label='V3')
call mma_allocate(T3,nOsc,nOsc,nOsc,label='T3')
call mma_allocate(V4,nOsc,nOsc,nOsc,nOsc,label='V4')
call mma_allocate(T4,nOsc,nOsc,nOsc,nOsc,label='T4')

! Calculate the eigenvector matrix, C, in dimensionless normal coordinates.
Temp(:,:) = Zero
do i=1,NumInt
  Temp(i,i) = One/sqrt(harmfreq(i))
end do
call DGEMM_('n','n',NumInt,NumInt,NumInt,One,eigenVec,NumInt,Temp,NumInt,Zero,C,NumInt)

call mma_deallocate(Temp)

! Transform cubic force constants to dimensionless normal coordinates
do i=1,NumInt
  do j=1,NumInt
    do k=1,NumInt
      sum1 = Zero
      sum2 = Zero
      do i1=1,NumInt
        do j1=1,NumInt
          do k1=1,NumInt
            coef = C(i1,i)*C(j1,j)*C(k1,k)
            sum1 = sum1+coef*D3(i1,j1,k1)
            sum2 = sum2+coef*Gprime(i1,j1,k1)
          end do
        end do
      end do
      V3(i,j,k) = sum1
      T3(i,j,k) = sum2
    end do
  end do
end do

! Transform quartic force constants to dimensionless normal coordinates
do i=1,NumInt
  do j=1,NumInt
    do k=1,NumInt
      do l=1,NumInt
        sum1 = Zero
        sum2 = Zero
        do i1=1,NumInt
          do j1=1,NumInt
            do k1=1,NumInt
              do l1=1,NumInt
                coef = C(i1,i)*C(j1,j)*C(k1,k)*C(l1,l)
                sum1 = sum1+coef*D4(i1,j1,k1,l1)
                sum2 = sum2+coef*Gdbleprime(i1,j1,k1,l1)
              end do
            end do
          end do
        end do
        V4(i,j,k,l) = sum1
        T4(i,j,k,l) = sum2
      end do
    end do
  end do
end do

call mma_deallocate(C)

! Calculate diagonal anharmonicity constants.
do i=1,NumInt
  x(i,i) = V4(i,i,i,i)/16.0_wp
  tmp = One/(48.0_wp*harmfreq(i))
  x(i,i) = x(i,i)-tmp*(V3(i,i,i)*(Five*V3(i,i,i)+Six*T3(i,i,i))+Nine*T3(i,i,i)**2)
  do j=1,NumInt
    if (j /= i) then
      tmp = -One/(16.0_wp*harmfreq(j)*(Four*harmfreq(i)**2-harmfreq(j)**2))
      x(i,i) = x(i,i)+tmp*(Eight*harmfreq(i)**2-Three*harmfreq(j)**2)*(V3(i,i,j)**2+T3(i,i,j)**2)
      x(i,i) = x(i,i)+tmp*(Eight*harmfreq(i)**2-harmfreq(j)**2)*(Two*V3(i,i,j)*T3(i,i,j))
      x(i,i) = x(i,i)-tmp*Eight*harmfreq(i)*harmfreq(j)*T3(i,j,i)*(V3(i,i,j)-T3(i,i,j))
      x(i,i) = x(i,i)-tmp*Four*harmfreq(j)**2*T3(i,j,i)**2
    end if
  end do
end do

! Calculate off-diagonal anharmonicity constants.
do i=1,NumInt
  do j=1,NumInt
    if (i /= j) then
      x(i,j) = V4(i,i,j,j)*Quart
      x(i,j) = x(i,j)+T4(i,i,j,j)*Half
      x(i,j) = x(i,j)-(V3(i,i,i)+T3(i,i,i))*(V3(i,j,j)+T3(j,j,i))/(Four*harmfreq(i))
      x(i,j) = x(i,j)-(V3(i,i,j)+T3(i,i,j))*(V3(j,j,j)+T3(j,j,j))/(Four*harmfreq(j))
      do k=1,NumInt
        if ((k /= i) .and. (k /= j)) then
          x(i,j) = x(i,j)-(V3(i,i,k)+T3(i,i,k))*(V3(j,j,k)+T3(j,j,k))/(Two*harmfreq(k))
        end if
      end do
      tmp = Half*harmfreq(i)/(Four*harmfreq(i)**2-harmfreq(j)**2)
      x(i,j) = x(i,j)-tmp*((V3(i,i,j)-T3(i,i,j))**2+Two*(V3(j,j,i)-T3(j,j,i))*T3(i,j,j)+Four*T3(i,j,i)**2)
      tmp = Half*harmfreq(j)/(Four*harmfreq(j)**2-harmfreq(i)**2)
      x(i,j) = x(i,j)-tmp*((V3(j,j,i)-T3(j,j,i))**2+Two*(V3(i,i,j)-T3(i,i,j))*T3(i,j,i)+Four*T3(i,j,j)**2)
      do k=1,NumInt
        if ((k /= i) .and. (k /= j)) then
          delta = One
          sum = harmfreq(i)+harmfreq(j)+harmfreq(k)
          delta = delta*sum
          sum = harmfreq(i)-harmfreq(j)-harmfreq(k)
          delta = delta*sum
          sum = -harmfreq(i)+harmfreq(j)-harmfreq(k)
          delta = delta*sum
          sum = -harmfreq(i)-harmfreq(j)+harmfreq(k)
          delta = delta*sum
          deltaInv = One/delta
          x(i,j) = x(i,j)+deltaInv*Half*harmfreq(k)*(harmfreq(i)**2+harmfreq(j)**2-harmfreq(k)**2)* &
                   (V3(i,j,k)**2+T3(i,j,k)**2+T3(j,k,i)**2+T3(k,i,j)**2)
          x(i,j) = x(i,j)+deltaInv*Two*(harmfreq(i)*harmfreq(j)*harmfreq(k))*(V3(i,j,k)*T3(i,j,k)-T3(j,k,i)*T3(k,i,j))
          x(i,j) = x(i,j)-deltaInv*harmfreq(i)*(harmfreq(k)**2+harmfreq(j)**2-harmfreq(i)**2)* &
                   (V3(i,j,k)*T3(k,i,j)-T3(i,j,k)*T3(j,k,i))
          x(i,j) = x(i,j)-deltaInv*harmfreq(j)*(harmfreq(k)**2+harmfreq(i)**2-harmfreq(j)**2)* &
                   (V3(i,j,k)*T3(j,k,i)-T3(i,j,k)*T3(k,i,j))
        end if
      end do
    end if
  end do
end do

call mma_deallocate(V3)
call mma_deallocate(T3)
call mma_deallocate(V4)
call mma_deallocate(T4)

end subroutine Anharm

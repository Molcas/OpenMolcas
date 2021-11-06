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

subroutine Anharm(eigenVec,harmfreq,D3,D4,Gprime,Gdbleprime,x,max_term,nOsc,C,Temp,V3,T3,V4,T4)
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
!    max_term   : Integer - highest power of term in polynomial fit.
!
!  Output:
!    x          : Real*8 two dimensional array - anharmonicity
!                 constants.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

implicit real*8(a-h,o-z)
#include "Constants_mula.fh"
#include "dims.fh"
real*8 eigenVec(nosc,nosc)
real*8 harmfreq(nosc)
real*8 D3(ngdim,ngdim,ngdim)
real*8 D4(ngdim,ngdim,ngdim,ngdim)
real*8 Gprime(ngdim,ngdim,ngdim)
real*8 Gdbleprime(ngdim,ngdim,ngdim,ngdim)
real*8 x(nosc,nosc)
real*8 C(nOsc,nOsc)
real*8 V3(nOsc,nOsc,nOsc)
real*8 V4(nOsc,nOsc,nOsc,nOsc)
real*8 T3(nOsc,nOsc,nOsc)
real*8 T4(nOsc,nOsc,nOsc,nOsc)
real*8 Temp(nOsc,nOsc)

NumInt = nOsc

! Calculate the eigenvector matrix, C, in dimensionless normal coordinates.
call dcopy_(NumInt**2,[0.0d0],0,Temp,1)
do i=1,NumInt
  Temp(i,i) = 1.0d0/sqrt(harmfreq(i))
end do
call DGEMM_('n','n',NumInt,NumInt,NumInt,1.0d0,eigenVec,NumInt,Temp,NumInt,0.0d0,C,NumInt)

! Transform cubic force constants to dimensionless normal coordinates
do i=1,NumInt
  do j=1,NumInt
    do k=1,NumInt
      sum1 = 0.0d0
      sum2 = 0.0d0
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
        sum1 = 0.0d0
        sum2 = 0.0d0
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

! Calculate diagonal anharmonicity constants.
do i=1,NumInt
  x(i,i) = V4(i,i,i,i)/16.0d0
  tmp = 1.0d0/(48.0d0*harmfreq(i))
  x(i,i) = x(i,i)-tmp*(V3(i,i,i)*(5.0d0*V3(i,i,i)+6.0d0*T3(i,i,i))+9.0d0*T3(i,i,i)**2)
  do j=1,NumInt
    if (j /= i) then
      tmp = -1.0d0/(16.0d0*harmfreq(j)*(4.0d0*harmfreq(i)**2-harmfreq(j)**2))
      x(i,i) = x(i,i)+tmp*(8.0d0*harmfreq(i)**2-3.0d0*harmfreq(j)**2)*(V3(i,i,j)**2+T3(i,i,j)**2)
      x(i,i) = x(i,i)+tmp*(8.0d0*harmfreq(i)**2-harmfreq(j)**2)*(2.0d0*V3(i,i,j)*T3(i,i,j))
      x(i,i) = x(i,i)-tmp*8.0d0*harmfreq(i)*harmfreq(j)*T3(i,j,i)*(V3(i,i,j)-T3(i,i,j))
      x(i,i) = x(i,i)-tmp*4.0d0*harmfreq(j)**2*T3(i,j,i)**2
    end if
  end do
end do

! Calculate off-diagonal anharmonicity constants.
do i=1,NumInt
  do j=1,NumInt
    if (i /= j) then
      x(i,j) = V4(i,i,j,j)/4.0d0
      x(i,j) = x(i,j)+T4(i,i,j,j)/2.0d0
      x(i,j) = x(i,j)-(V3(i,i,i)+T3(i,i,i))*(V3(i,j,j)+T3(j,j,i))/(4.0d0*harmfreq(i))
      x(i,j) = x(i,j)-(V3(i,i,j)+T3(i,i,j))*(V3(j,j,j)+T3(j,j,j))/(4.0d0*harmfreq(j))
      do k=1,NumInt
        if ((k /= i) .and. (k /= j)) then
          x(i,j) = x(i,j)-(V3(i,i,k)+T3(i,i,k))*(V3(j,j,k)+T3(j,j,k))/(2.0d0*harmfreq(k))
        end if
      end do
      tmp = 0.5d0*harmfreq(i)/(4.0d0*harmfreq(i)**2-harmfreq(j)**2)
      x(i,j) = x(i,j)-tmp*((V3(i,i,j)-T3(i,i,j))**2+2.0d0*(V3(j,j,i)-T3(j,j,i))*T3(i,j,j)+4.0d0*T3(i,j,i)**2)
      tmp = 0.5d0*harmfreq(j)/(4.0d0*harmfreq(j)**2-harmfreq(i)**2)
      x(i,j) = x(i,j)-tmp*((V3(j,j,i)-T3(j,j,i))**2+2.0d0*(V3(i,i,j)-T3(i,i,j))*T3(i,j,i)+4.0d0*T3(i,j,j)**2)
      do k=1,NumInt
        if ((k /= i) .and. (k /= j)) then
          delta = 1.0d0
          sum = harmfreq(i)+harmfreq(j)+harmfreq(k)
          delta = delta*sum
          sum = harmfreq(i)-harmfreq(j)-harmfreq(k)
          delta = delta*sum
          sum = -harmfreq(i)+harmfreq(j)-harmfreq(k)
          delta = delta*sum
          sum = -harmfreq(i)-harmfreq(j)+harmfreq(k)
          delta = delta*sum
          deltaInv = 1.0d0/delta
          x(i,j) = x(i,j)+deltaInv*0.5d0*harmfreq(k)*(harmfreq(i)**2+harmfreq(j)**2-harmfreq(k)**2)* &
                   (V3(i,j,k)**2+T3(i,j,k)**2+T3(j,k,i)**2+T3(k,i,j)**2)
          x(i,j) = x(i,j)+deltaInv*2.0d0*(harmfreq(i)*harmfreq(j)*harmfreq(k))*(V3(i,j,k)*T3(i,j,k)-T3(j,k,i)*T3(k,i,j))
          x(i,j) = x(i,j)-deltaInv*harmfreq(i)*(harmfreq(k)**2+harmfreq(j)**2-harmfreq(i)**2)* &
                   (V3(i,j,k)*T3(k,i,j)-T3(i,j,k)*T3(j,k,i))
          x(i,j) = x(i,j)-deltaInv*harmfreq(j)*(harmfreq(k)**2+harmfreq(i)**2-harmfreq(j)**2)* &
                   (V3(i,j,k)*T3(j,k,i)-T3(i,j,k)*T3(k,i,j))
        end if
      end do
    end if
  end do
end do

! Avoid unused argument warnings
if (.false.) call Unused_integer(max_term)

end subroutine Anharm

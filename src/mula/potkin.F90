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
! Copyright (C) 1996,1999, Niclas Forsberg                             *
!               1996,1999, Anders Bernhardsson                         *
!***********************************************************************

!module PotKin

!  Contains:
!    PotEnergy      (A,nMat,energy,grad,Hess,D3,D4,max_term)
!    KinEnergy      (A,nMat,G,Gprime,Gdbleprime,max_term)
!
!  Written by:
!    Niclas Forsberg & Anders Bernhardsson,
!    Dept. of Theoretical Chemistry, Lund University, 1996.
!    Dept. of Theoretical Chemistry, Lund University, 1999.

!vv private

!contains

subroutine KinEnergy(A,nMat,iCre,iAnn,G,Gprime,Gdbleprime,max_term,C,W,alpha1,alpha2,beta,r_diff,max_Ord,nOsc,nOscOld)
! Out : A
! In  :
!      Real G,Gprime,Gdbleprime,W,C,alpha1,alpha2,beta,r_diff
!      integer nMat,iCre,iAnn,max_term
!
!  Purpose:
!    Calculate matrix elements of kinetic energy terms.
!
!  Written by:
!    Niclas Forsberg,Anders Bernhardsson
!    Dept. of Theoretical Chemistry, Lund University, 1996.
!    Dept. of Theoretical Chemistry, Lund University, 1999.

use mula_global, only: ndim1, ndim2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Three, Six, Half, OneHalf
use Definitions, only: wp

!use TabMod
implicit real*8(a-h,o-z)
parameter(Thrs=1.0e-15_wp)
real*8 A(0:max_Ord,0:max_Ord)
integer nMat(0:ndim1,ndim2)
integer iAnn(0:ndim1,ndim2)
integer iCre(0:ndim1,ndim2)
real*8 rdx(4)
real*8 alpha1(nosc,nosc)
real*8 alpha2(nosc,nosc)
real*8 beta(nosc,nosc)
real*8 C(nosc,nosc)
real*8 W(nosc,nosc)
real*8 G(nosc,nosc)
real*8 r_diff(noscold)
real*8 Gprime(nosc,nosc,nosc)
real*8 Gdbleprime(nosc,nosc,nosc,nosc)
real*8, allocatable :: G_2(:,:), r_temp(:), T1(:), T2(:), Temp(:,:), Temp1(:,:,:), Temp2(:,:,:), Temp3(:,:,:,:), Temp4(:,:,:,:), &
                       Tempa(:,:), Tempb(:,:)

! Initialize.
call mma_allocate(Temp,nOsc,nOsc,label='Temp')
call mma_allocate(Tempa,nOsc,nOsc,label='Tempa')
call mma_allocate(Tempb,nOsc,nOsc,label='Tempb')
call mma_allocate(Temp1,nOsc,nOsc,nOsc,label='Temp1')
call mma_allocate(Temp2,nOsc,nOsc,nOsc,label='Temp2')
call mma_allocate(Temp3,nOsc,nOsc,nOsc,nOsc,label='Temp3')
call mma_allocate(Temp4,nOsc,nOsc,nOsc,nOsc,label='Temp4')
call mma_allocate(r_temp,nOsc,label='r_temp')
call mma_allocate(G_2,nOsc,nOsc,label='G_2')
call mma_allocate(T1,nOsc,label='T1')
call mma_allocate(T2,nOsc,label='T2')
nOscSqr = nOsc**2
r_norm = Dnrm2_(nOsc,r_diff,1)
ran = One
call Dgesub(alpha1,nOsc,'N',alpha2,nOsc,'N',Tempa,nOsc,nOsc,nOsc)
alpha_norm = Dnrm2_(nOscSqr,Tempa,1)
if ((r_norm > Thrs) .or. (alpha_norm > Thrs)) then
  call Dgesub(alpha1,nOsc,'N',alpha2,nOsc,'N',Tempa,nOsc,nOsc,nOsc)
  call DGEMM_('N','N',nOsc,nOsc,nOsc,One,Tempa,nOsc,W,nOsc,Zero,Tempb,nOsc)
  call DGEMM_('N','N',nOsc,1,nOsc,-Two,beta,nOsc,r_diff,nOsc,Zero,r_temp,nOsc)
  rdx(1) = One
  rdx(2) = -One
  call DGEMM_('T','N',nosc,nosc,nosc,One,G,nosc,Tempb,nosc,Zero,Temp,nosc)
  call DGEMM_('T','T',nosc,nosc,nosc,-ran,Temp,nosc,C,nosc,Zero,G_2,nosc)
  call Mul2(nmat,A,iCre,iAnn,G_2,max_ord,nosc,rdx)
  rdx(1) = -One
  rdx(2) = One
  call DGEMM_('T','T',nosc,nosc,nosc,One,G,nosc,C,nosc,Zero,Temp,nosc)
  call DGEMM_('T','N',nosc,nosc,nosc,-ran,Temp,nosc,Tempb,nosc,Zero,G_2,nosc)
  call Mul2(nmat,A,iCre,iAnn,G_2,max_ord,nosc,rdx)
  rdx(1) = One
  rdx(2) = One
  call DGEMM_('T','N',nosc,nosc,nosc,One,G,nosc,Tempb,nosc,Zero,Temp,nosc)
  call DGEMM_('T','N',nosc,nosc,nosc,-One,Temp,nosc,Tempb,nosc,Zero,G_2,nosc)
  call Mul2(nmat,A,iCre,iAnn,G_2,max_ord,nosc,rdx)
  call dGeMV_('N',nosc,nosc,One,G,nosc,r_temp,1,Zero,T1,1)
  call dGeMV_('T',nosc,nosc,One,G,nosc,r_temp,1,One,T1,1)
  call dGeMV_('N',nosc,nosc,-Half*ran,C,nosc,T1,1,Zero,T2,1)
  rdx(1) = -One
  call Mul1(nmat,A,iCre,iAnn,T2,max_ord,nosc,rdx)
  call dGeMV_('N',nosc,nosc,One,G,nosc,r_temp,1,Zero,T1,1)
  call dGeMV_('T',nosc,nosc,One,G,nosc,r_temp,1,One,T1,1)
  call dGeMV_('T',nosc,nosc,-Half,tempb,nosc,T1,1,Zero,T2,1)
  rdx(1) = One
  call Mul1(nmat,A,iCre,iAnn,T2,max_ord,nosc,rdx)
  call dGeMV_('N',nosc,nosc,One,G,nosc,r_temp,1,Zero,T1,1)
  r = Ddot_(nosc,T1,1,r_temp,1)
  do i=0,max_Ord
    A(i,i) = A(i,i)-Half*r
  end do
  if (max_term > 2) then
    call DGEMM_('T','N',nosc**2,nosc,nosc,One,Gprime,nosc,tempb,nosc,Zero,Temp1,nosc**2)
    call DGEMM_('T','T',nosc**2,nosc,nosc,One,Temp1,nosc,C,nosc,Zero,Temp2,nosc**2)
    call DGEMM_('T','N',nosc**2,nosc,nosc,One,Temp2,nosc,W,nosc,Zero,Temp1,nosc**2)
    rdx(1) = One
    rdx(2) = One
    rdx(3) = -One
    do i=1,nosc
      do j=1,nosc
        do k=1,nosc
          Temp2(i,k,j) = -Three*ran*Temp1(i,j,k)
        end do
      end do
    end do
    call Mul3(nmat,A,iCre,iAnn,Temp2,max_ord,nosc,rdx)
    call DGEMM_('T','T',nosc**2,nosc,nosc,One,Gprime,nosc,C,nosc,Zero,Temp1,nosc**2)
    call DGEMM_('T','N',nosc**2,nosc,nosc,One,Temp1,nosc,Tempb,nosc,Zero,Temp2,nosc**2)
    call DGEMM_('T','N',nosc**2,nosc,nosc,One,Temp2,nosc,W,nosc,Zero,Temp1,nosc**2)
    rdx(1) = -One
    rdx(2) = One
    rdx(3) = One
    do i=1,nosc
      do j=1,nosc
        do k=1,nosc
          Temp2(i,k,j) = -Three*ran*Temp1(i,j,k)
        end do
      end do
    end do
    call Mul3(nmat,A,iCre,iAnn,Temp2,max_ord,nosc,rdx)
    call DGEMM_('T','N',nosc**2,nosc,nosc,One,Gprime,nosc,Tempb,nosc,Zero,Temp1,nosc**2)
    call DGEMM_('T','N',nosc**2,nosc,nosc,One,Temp1,nosc,Tempb,nosc,Zero,Temp2,nosc**2)
    call DGEMM_('T','N',nosc**2,nosc,nosc,One,Temp2,nosc,W,nosc,Zero,Temp1,nosc**2)
    rdx(1) = One
    rdx(2) = One
    rdx(3) = One
    do i=1,nosc
      do j=1,nosc
        do k=1,nosc
          Temp2(i,k,j) = -Three*Temp1(i,j,k)
        end do
      end do
    end do
    call Mul3(nmat,A,iCre,iAnn,Temp2,max_ord,nosc,rdx)

    Temp(:,:) = Zero
    do i=1,nosc
      do j=1,nosc
        do k=1,nosc
          Temp(k,j) = Temp(k,j)-ran*Gprime(i,j,k)*r_temp(i)
        end do
      end do
    end do
    rdx(1) = One
    rdx(2) = -One
    call DGEMM_('T','N',nosc,nosc,nosc,One,Temp,nosc,W,nosc,Zero,G_2,nosc)
    call DGEMM_('T','T',nosc,nosc,nosc,One,G_2,nosc,C,nosc,Zero,Temp,nosc)
    call Mul2(nmat,A,iCre,iAnn,Temp,max_ord,nosc,rdx)

    Temp(:,:) = Zero
    do i=1,nosc
      do j=1,nosc
        do k=1,nosc
          Temp(k,j) = Temp(k,j)-Gprime(i,j,k)*r_temp(i)
        end do
      end do
    end do
    rdx(1) = One
    rdx(2) = One
    call DGEMM_('T','N',nosc,nosc,nosc,One,Temp,nosc,W,nosc,Zero,G_2,nosc)
    call DGEMM_('T','N',nosc,nosc,nosc,One,G_2,nosc,Tempb,nosc,Zero,Temp,nosc)
    call Mul2(nmat,A,iCre,iAnn,Temp,max_ord,nosc,rdx)

    Temp(:,:) = Zero
    do i=1,nosc
      do j=1,nosc
        do k=1,nosc
          Temp(i,k) = Temp(i,k)-ran*Gprime(i,j,k)*r_temp(j)
        end do
      end do
    end do
    rdx(1) = -One
    rdx(2) = One
    call DGEMM_('T','T',nosc,nosc,nosc,One,Temp,nosc,C,nosc,Zero,G_2,nosc)
    call DGEMM_('T','N',nosc,nosc,nosc,One,G_2,nosc,W,nosc,Zero,Temp,nosc)
    call Mul2(nmat,A,iCre,iAnn,Temp,max_ord,nosc,rdx)

    Temp(:,:) = Zero
    do i=1,nosc
      do j=1,nosc
        do k=1,nosc
          Temp(i,k) = Temp(i,k)-Gprime(i,j,k)*r_temp(j) !!
        end do
      end do
    end do
    rdx(1) = One
    rdx(2) = One
    call DGEMM_('T','n',nosc,nosc,nosc,One,Temp,nosc,Tempb,nosc,Zero,G_2,nosc)
    call DGEMM_('T','n',nosc,nosc,nosc,One,G_2,nosc,W,nosc,Zero,Temp,nosc)
    call Mul2(nmat,A,iCre,iAnn,Temp,max_ord,nosc,rdx)

    t1(:) = Zero
    do i=1,nosc
      do j=1,nosc
        do k=1,nosc
          T1(k) = -Half*Gprime(i,j,k)*r_temp(i)*r_temp(j)+t1(k)
        end do
      end do
    end do
    call DGEMM_('T','n',1,nosc,nosc,One,T1,nosc,W,nosc,Zero,T2,1)
    rdx(1) = One
    call Mul1(nmat,A,iCre,iAnn,T2,max_ord,nosc,rdx)
  end if
  if (max_term > 3) then
    call DGEMM_('T','N',nosc**3,nosc,nosc,One,Gdbleprime,nosc,tempb,nosc,Zero,Temp3,nosc**3)
    call DGEMM_('T','T',nosc**3,nosc,nosc,One,Temp3,nosc,C,nosc,Zero,Temp4,nosc**3)
    call DGEMM_('T','N',nosc**3,nosc,nosc,One,Temp4,nosc,W,nosc,Zero,Temp3,nosc**3)
    call DGEMM_('T','N',nosc**3,nosc,nosc,One,Temp3,nosc,W,nosc,Zero,Temp4,nosc**3)
    rdx(1) = One
    rdx(2) = One
    rdx(3) = One
    rdx(4) = -One
    do i=1,nosc
      do j=1,nosc
        do k=1,nosc
          do l=1,nosc
            Temp3(i,k,l,j) = -ran*Six*Temp4(i,j,k,l)
          end do
        end do
      end do
    end do
    call Mul4(nmat,A,iCre,iAnn,Temp3,max_ord,nosc,rdx)
    call DGEMM_('T','T',nosc**3,nosc,nosc,One,Gdbleprime,nosc,C,nosc,Zero,Temp3,nosc**3)
    call DGEMM_('T','N',nosc**3,nosc,nosc,One,Temp3,nosc,Tempb,nosc,Zero,Temp4,nosc**3)
    call DGEMM_('T','N',nosc**3,nosc,nosc,One,Temp4,nosc,W,nosc,Zero,Temp3,nosc**3)
    call DGEMM_('T','N',nosc**3,nosc,nosc,One,Temp3,nosc,W,nosc,Zero,Temp4,nosc**3)
    rdx(1) = -One
    rdx(2) = One
    rdx(3) = One
    rdx(4) = One
    do i=1,nosc
      do j=1,nosc
        do k=1,nosc
          do l=1,nosc
            Temp3(i,k,l,j) = -ran*Six*Temp4(i,j,k,l)
          end do
        end do
      end do
    end do
    call Mul4(nmat,A,iCre,iAnn,Temp3,max_ord,nosc,rdx)
    call DGEMM_('T','n',nosc**3,nosc,nosc,One,Gdbleprime,nosc,Tempb,nosc,Zero,Temp3,nosc**3)
    call DGEMM_('T','N',nosc**3,nosc,nosc,One,Temp3,nosc,Tempb,nosc,Zero,Temp4,nosc**3)
    call DGEMM_('T','N',nosc**3,nosc,nosc,One,Temp4,nosc,W,nosc,Zero,Temp3,nosc**3)
    call DGEMM_('T','N',nosc**3,nosc,nosc,One,Temp3,nosc,W,nosc,Zero,Temp4,nosc**3)
    rdx(1) = One
    rdx(2) = One
    rdx(3) = One
    rdx(4) = One
    do i=1,nosc
      do j=1,nosc
        do k=1,nosc
          do l=1,nosc
            Temp3(i,k,l,j) = -Six*Temp4(i,j,k,l)
          end do
        end do
      end do
    end do
    call Mul4(nmat,A,iCre,iAnn,Temp3,max_ord,nosc,rdx)

    Temp1(:,:,:) = Zero
    do i=1,nosc
      do j=1,nosc
        do k=1,nosc
          do l=1,nosc
            Temp1(k,l,j) = Temp1(k,l,j)-ran*OneHalf*Gdbleprime(i,j,k,l)*r_temp(i)
          end do
        end do
      end do
    end do
    rdx(1) = One
    rdx(2) = One
    rdx(3) = -One
    call DGEMM_('T','n',nosc**2,nosc,nosc,One,Temp1,nosc,W,nosc,Zero,Temp2,nosc**2)
    call DGEMM_('T','n',nosc**2,nosc,nosc,One,Temp2,nosc,W,nosc,Zero,Temp1,nosc**2)
    call DGEMM_('T','T',nosc**2,nosc,nosc,One,Temp1,nosc,C,nosc,Zero,Temp2,nosc**2)
    call Mul3(nmat,A,iCre,iAnn,Temp2,max_ord,nosc,rdx)
    Temp1(:,:,:) = Zero
    do i=1,nosc
      do j=1,nosc
        do k=1,nosc
          do l=1,nosc
            Temp1(k,l,j) = Temp1(k,l,j)-OneHalf*Gdbleprime(i,j,k,l)*r_temp(i)
          end do
        end do
      end do
    end do
    rdx(1) = One
    rdx(2) = One
    rdx(3) = One
    call DGEMM_('T','n',nosc**2,nosc,nosc,One,Temp1,nosc,W,nosc,Zero,Temp2,nosc**2)
    call DGEMM_('T','n',nosc**2,nosc,nosc,One,Temp2,nosc,W,nosc,Zero,Temp1,nosc**2)
    call DGEMM_('T','n',nosc**2,nosc,nosc,One,Temp1,nosc,Tempb,nosc,Zero,Temp2,nosc**2)
    call Mul3(nmat,A,iCre,iAnn,Temp2,max_ord,nosc,rdx)

    Temp1(:,:,:) = Zero
    do i=1,nosc
      do j=1,nosc
        do k=1,nosc
          do l=1,nosc
            Temp1(i,k,l) = Temp1(i,k,l)-ran*OneHalf*Gdbleprime(i,j,k,l)*r_temp(j)
          end do
        end do
      end do
    end do
    rdx(1) = -One
    rdx(2) = One
    rdx(3) = One
    call DGEMM_('T','T',nosc**2,nosc,nosc,One,Temp1,nosc,C,nosc,Zero,Temp2,nosc**2)
    call DGEMM_('T','N',nosc**2,nosc,nosc,One,Temp2,nosc,W,nosc,Zero,Temp1,nosc**2)
    call DGEMM_('T','N',nosc**2,nosc,nosc,One,Temp1,nosc,W,nosc,Zero,Temp2,nosc**2)
    call Mul3(nmat,A,iCre,iAnn,Temp2,max_ord,nosc,rdx)

    Temp1(:,:,:) = Zero
    do i=1,nosc
      do j=1,nosc
        do k=1,nosc
          do l=1,nosc
            Temp1(i,k,l) = Temp1(i,k,l)-OneHalf*Gdbleprime(i,j,k,l)*r_temp(j)
          end do
        end do
      end do
    end do
    rdx(1) = One
    rdx(2) = One
    rdx(3) = One
    call DGEMM_('T','n',nosc**2,nosc,nosc,One,Temp1,nosc,Tempb,nosc,Zero,Temp2,nosc**2)
    call DGEMM_('T','N',nosc**2,nosc,nosc,One,Temp2,nosc,W,nosc,Zero,Temp1,nosc**2)
    call DGEMM_('T','N',nosc**2,nosc,nosc,One,Temp1,nosc,W,nosc,Zero,Temp2,nosc**2)
    call Mul3(nmat,A,iCre,iAnn,Temp2,max_ord,nosc,rdx)

    rdx(1) = One
    rdx(2) = One
    Temp(:,:) = Zero
    do i=1,nosc
      do j=1,nosc
        do k=1,nosc
          do l=1,nosc
            Temp(k,l) = Temp(k,l)-Half*Gdbleprime(i,j,k,l)*r_temp(i)*r_temp(j)
          end do
        end do
      end do
    end do
    call DGEMM_('T','N',nosc,nosc,nosc,One,Temp,nosc,W,nosc,Zero,G_2,nosc)
    call DGEMM_('T','N',nosc,nosc,nosc,One,G_2,nosc,W,nosc,Zero,Temp,nosc)
    call Mul2(nmat,A,iCre,iAnn,Temp,max_ord,nosc,rdx)
  end if
end if

call mma_deallocate(Tempa)
call mma_deallocate(Tempb)
call mma_deallocate(r_temp)
call mma_deallocate(T1)
call mma_deallocate(T2)

! If higher terms than quadratic are used in the polynomial
! fit of the potential surface, then we have to use a
! Taylor expansion of the inverse mass tensor.
rdx(1) = -One
rdx(2) = -One
call DGEMM_('T','T',nosc,nosc,nosc,One,G,nosc,C,nosc,Zero,Temp,nosc)
call DGEMM_('T','T',nosc,nosc,nosc,-One,Temp,nosc,C,nosc,Zero,G_2,nosc)
call Mul2(nmat,A,iCre,iAnn,G_2,max_ord,nosc,rdx)
if (max_term > 2) then
  call DGEMM_('T','T',nosc**2,nosc,nosc,One,Gprime,nosc,C,nosc,Zero,Temp1,nosc**2)
  call DGEMM_('T','T',nosc**2,nosc,nosc,One,Temp1,nosc,C,nosc,Zero,Temp2,nosc**2)
  call DGEMM_('T','N',nosc**2,nosc,nosc,One,Temp2,nosc,W,nosc,Zero,Temp1,nosc**2)
  rdx(1) = -One
  rdx(2) = One
  rdx(3) = -One
  do i=1,nosc
    do j=1,nosc
      do k=1,nosc
        Temp2(i,k,j) = -Three*Temp1(i,j,k)
      end do
    end do
  end do
  call Mul3(nmat,A,iCre,iAnn,Temp2,max_ord,nosc,rdx)
end if
if (max_term > 3) then
  call DGEMM_('T','T',nosc**3,nosc,nosc,One,Gdbleprime,nosc,C,nosc,Zero,Temp3,nosc**3)
  call DGEMM_('T','T',nosc**3,nosc,nosc,One,Temp3,nosc,C,nosc,Zero,Temp4,nosc**3)
  call DGEMM_('T','N',nosc**3,nosc,nosc,One,Temp4,nosc,W,nosc,Zero,Temp3,nosc**3)
  call DGEMM_('T','N',nosc**3,nosc,nosc,One,Temp3,nosc,W,nosc,Zero,Temp4,nosc**3)

  rdx(1) = -One
  rdx(2) = One
  rdx(3) = One
  rdx(4) = -One
  do i=1,nosc
    do j=1,nosc
      do k=1,nosc
        do l=1,nosc
          Temp3(i,k,l,j) = -Six*Temp4(i,j,k,l)
        end do
      end do
    end do
  end do
  call Mul4(nmat,A,iCre,iAnn,Temp3,max_ord,nosc,rdx)
end if

call mma_deallocate(Temp)
call mma_deallocate(Temp1)
call mma_deallocate(Temp2)
call mma_deallocate(Temp3)
call mma_deallocate(Temp4)
call mma_deallocate(G_2)

end subroutine KinEnergy

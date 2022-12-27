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

!contains

subroutine PotEnergy(A,nMat,iCre,iAnn,energy,grad,Hess,D3,D4,max_term,W,max_Ord,nOsc,nOscOld)
!  Purpose:
!    Calculate matrix elements of potential energy terms.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

use mula_global, only: ndim1, ndim2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nMat(0:ndim1,ndim2), iAnn(0:ndim1,ndim2), iCre(0:ndim1,ndim2), max_term, max_Ord, nOsc, nOscOld
real(kind=wp), intent(out) :: A(0:max_Ord,0:max_Ord)
real(kind=wp), intent(in) :: energy, grad(nOscOld), Hess(nOscOld,nOscOld), D3(nOscOld,nOscOld,nOscOld), &
                             D4(nOscOld,nOscOld,nOscOld,nOscOld), W(nOscOld,nOsc)
integer(kind=iwp) :: i
real(kind=wp) :: rdx(4)
real(kind=wp), allocatable :: D3_2(:,:,:), D4_2(:,:,:,:), grad_2(:), Hess_2(:,:), Temp(:)

! Zeroth order term.
do i=0,max_Ord
  A(i,i) = Energy
end do
rdx(1) = One
rdx(2) = One
rdx(3) = One
rdx(4) = One
call mma_allocate(Temp,nOscOld**4,label='Temp')

! First order terms.
if (max_term > 0) then
  call mma_allocate(grad_2,nOsc,label='grad_2')
  call DGEMM_('T','N',1,nOsc,nOscOld,One,grad,nOscOld,W,nOscOld,Zero,grad_2,1)
  call Mul1(nMat,A,icre,iann,grad_2,max_Ord,nOsc,rdx)
  call mma_deallocate(grad_2)
end if

! Second order terms.
if (max_term > 1) then
  call mma_allocate(Hess_2,nOsc,nOsc,label='Hess_2')
  call DGEMM_('T','N',nOscOld,nOsc,nOscOld,One,Hess,nOscOld,W,nOscOld,Zero,Temp,nOscOld)
  call DGEMM_('T','N',nOsc,nOsc,nOscOld,One,Temp,nOscOld,W,nOscOld,Zero,Hess_2,nOsc)
  call Mul2(nMat,A,icre,iann,Hess_2,max_Ord,nOsc,rdx)
  call mma_deallocate(Hess_2)
end if

! Third order terms.
if (max_term > 2) then
  call mma_allocate(D3_2,nOsc,nOsc,nOsc,label='D3_2')
  call DGEMM_('T','N',nOscOld**2,nOsc,nOscOld,One,D3,nOscOld,W,nOscOld,Zero,Temp,nOscOld**2)
  call DGEMM_('T','N',nOsc*nOscOld,nOsc,nOscOld,One,Temp,nOscOld,W,nOscOld,Zero,D3_2,nOsc*nOscOld)
  call DGEMM_('T','N',nOsc**2,nOsc,nOscOld,One,D3_2,nOscOld,W,nOscOld,Zero,Temp,nOsc**2)
  call dcopy_(nOsc**3,Temp,1,D3_2,1)
  call Mul3(nMat,A,icre,iann,D3_2,max_Ord,nOsc,rdx)
  call mma_deallocate(D3_2)
end if

! Fourth order terms.
if (max_term > 3) then
  call mma_allocate(D4_2,nOsc,nOsc,nOsc,nOsc,label='D4_2')
  call DGEMM_('T','N',nOscOld**3,nOsc,nOscOld,One,D4,nOscOld,W,nOscOld,Zero,Temp,nOscOld**3)
  call DGEMM_('T','N',nOsc*nOscOld**2,nOsc,nOscOld,One,Temp,nOscOld,W,nOscOld,Zero,D4_2,nOsc*nOscOld**2)
  call DGEMM_('T','N',nOsc**2*nOscOld,nOsc,nOscOld,One,D4_2,nOscOld,W,nOscOld,Zero,Temp,nOsc**2*nOscOld)
  call DGEMM_('T','N',nOsc**3,nOsc,nOscOld,One,Temp,nOscOld,W,nOscOld,Zero,D4_2,nOsc**3)
  call Mul4(nMat,A,icre,iann,D4_2,max_Ord,nOsc,rdx)
  call mma_deallocate(D4_2)
end if

call mma_deallocate(Temp)

end subroutine PotEnergy

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
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nMat(0:ndim1,ndim2), iCre(0:ndim1,ndim2), iAnn(0:ndim1,ndim2), max_term, max_Ord, nOsc, nOscOld
real(kind=wp), intent(out) :: A(0:max_Ord,0:max_Ord)
real(kind=wp), intent(in) :: G(nOsc,nOsc), Gprime(nOsc,nOsc,nOsc), Gdbleprime(nOsc,nOsc,nOsc,nOsc), C(nOsc,nOsc), W(nOsc,nOsc), &
                             alpha1(nOsc,nOsc), alpha2(nOsc,nOsc), beta(nOsc,nOsc), r_diff(nOscOld)
integer(kind=iwp) :: i, j, k, l, nOscSqr
real(kind=wp) :: alpha_norm, r, r_norm, ran, rdx(4)
real(kind=wp), allocatable :: G_2(:,:), r_temp(:), T1(:), T2(:), Temp(:,:), Temp1(:,:,:), Temp2(:,:,:), Temp3(:,:,:,:), &
                              Temp4(:,:,:,:), Tempa(:,:), Tempb(:,:)
real(kind=wp), parameter :: Thrs = 1.0e-15_wp
real(kind=wp), external :: Ddot_, Dnrm2_

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

!end module PotKin

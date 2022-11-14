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
! Copyright (C) 1995, Niclas Forsberg                                  *
!***********************************************************************

!module FCMod

!  Contains:
!    FCval     (C1,W1,det1,r01,C2,W2,det2,r02,FC,max_mOrd,max_nOrd,max_nOrd2,max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,nInc,mDec,
!               nDec,C,W,det0,L,U,FC00,alpha1,alpha2,beta)
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

!contains

subroutine FCval(C1,W1,det1,r01,C2,W2,det2,r02,FC,max_mOrd,max_nOrd,max_nOrd2,max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,nInc, &
                 mDec,nDec,C,W,det0,L,U,FC00,alpha1,alpha2,beta,nOsc)
!  Purpose:
!    Calculate multidimensional Franck Condon factors.
!
!  Input:
!    W1,W2      : Real two dimensional arrays - eigenvectors scaled by the square root of the eigenvalues.
!    C1,C2      : Real two dimensional arrays - inverses of W1 and W2.
!    det1,det2  : Real variables - determinants of C1 and C2.
!    r01,r02    : Real arrays - coordinates of the two oscillators.
!    max_mOrd,
!    max_nOrd,
!    max_nOrd2
!    max_mInc,
!    max_nInc,
!    max_nInc2  : Integer variables
!    mMat,nMat,
!    mInc,nInc,
!    mDec,nDec  : Two dimensional integer arrays
!
!  Output:
!    W          : Real two dimensional array - eigenvectors of intermediate oscillator scaled by the square root of the eigenvalues.
!    C          : Real two dimensional array - inverse of W.
!    det0       : Real variable - determinant of C.
!    L,U        : Real two dimensional arrays
!    FC00       : Real variable - zero-zero overlap.
!    FC         : Real
!    alpha1     : Real two dimensional array - 0.5*C1(T)*C1.
!    alpha2     : Real two dimensional array - 0.5*C2(T)*C2.
!    beta       : Real two dimensional array - 0.5*alpha1*alpha^(-1)*alpha2
!
!  Calls:
!    Dool_MULA
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

use mula_global, only: mdim1, mdim2, ndim1, ndim2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: max_mOrd, max_nOrd, max_nOrd2, max_mInc, max_nInc, max_nInc2, mMat(0:mdim1,mdim2), &
                                 nMat(0:ndim1,ndim2), mInc(0:mdim1,mdim2), nInc(0:ndim1,ndim2), mDec(0:mdim1,mdim2), &
                                 nDec(0:ndim1,ndim2), nOsc
real(kind=wp), intent(in) :: C1(nOsc,nOsc), W1(nOsc,nOsc), det1, r01(nOsc), C2(nOsc,nOsc), W2(nOsc,nOsc), det2, r02(nOsc), &
                             C(nOsc,nOsc), W(nOsc,nOsc), det0
real(kind=wp), intent(out) :: FC(0:max_mord,0:max_nord), L(0:max_mord,0:max_ninc2), U(0:max_nord,0:max_nord2), FC00, &
                              alpha1(nOsc,nOsc), alpha2(nOsc,nOsc), beta(nOsc,nOsc)
integer(kind=iwp) :: i, iOrd, j, jOrd, kOrd, kOsc, kOsc_start, lOsc, nMaxMat, nTabDim
real(kind=wp) :: const, det, FC00_exp, rsum
real(kind=wp), allocatable :: A1(:,:), A1B1T(:,:), A2(:,:), A2B2T(:,:), alpha(:,:), B1(:,:), B2(:,:), d1(:), d2(:), r_temp1(:), &
                              r_temp2(:), sqr(:), temp(:,:), temp1(:,:), temp2(:,:)
real(kind=wp), external :: Ddot_

! Initialize.
nMaxMat = max(max_mord+1,max_nord+1)
nTabDim = max(nMaxMat,8)
call mma_allocate(temp,nOsc,nOsc,label='temp')
call mma_allocate(temp1,nOsc,nOsc,label='temp1')
call mma_allocate(temp2,nOsc,nOsc,label='temp2')

! Setup sqr table.
call mma_allocate(sqr,[0,nTabDim+1],label='sqr')
do i=0,nTabDim+1
  sqr(i) = sqrt(real(i,kind=wp))
end do

! Calculate alpha1, alpha2 and alpha.
call mma_allocate(alpha,nOsc,nOsc,label='alpha')
call DGEMM_('T','N',nOsc,nOsc,nOsc,Half,C1,nOsc,C1,nOsc,Zero,alpha1,nOsc)
call DGEMM_('T','N',nOsc,nOsc,nOsc,Half,C2,nOsc,C2,nOsc,Zero,alpha2,nOsc)
temp(:,:) = alpha1+alpha2
alpha(:,:) = Half*temp

!call xxDgemul(C,nOsc,'T',C,nOsc,'N',alpha,nOsc,nOsc,nOsc,nOsc)
!alpha(:,:) = Half*alpha

! Calculate C using a Cholesky factorization of 2*alpha.
!call Cholesky(temp,C)

! Calculate W.
!call unitmat(W,nOsc)
!temp(:,:) = C
!call Dool_MULA(temp,W,det0)

! Calculate r00.
call mma_allocate(r_temp1,nOsc,label='r_temp1')
call mma_allocate(r_temp2,nOsc,label='r_temp2')

! Calculate beta.
do i=1,nOsc
  do j=1,nOsc
    temp(j,i) = C1(i,j)
  end do
end do
temp1(:,:) = alpha1
temp(:,:) = Two*alpha

call Dool_MULA(temp,nOsc,nOsc,temp1,nOsc,nOsc,det)
call DGEMM_('N','N',nOsc,nOsc,nOsc,One,alpha2,nOsc,temp1,nOsc,Zero,beta,nOsc)

! Calculate FC00.
!r_temp1(:) = r02-r01
r_temp1(:) = r01-r02

call DGEMM_('N','N',nOsc,1,nOsc,One,beta,nOsc,r_temp1,nOsc,Zero,r_temp2,nOsc)
FC00_exp = Ddot_(nOsc,r_temp1,1,r_temp2,1)
FC00 = (sqrt(det1)*sqrt(det2)/det0)*exp(-FC00_exp)

! Calculate A, B and d matrices.
call mma_allocate(A1,nOsc,nOsc,label='A1')
call mma_allocate(B1,nOsc,nOsc,label='B1')
call mma_allocate(A2,nOsc,nOsc,label='A2')
call mma_allocate(B2,nOsc,nOsc,label='B2')
call mma_allocate(d1,nOsc,label='d1')
call mma_allocate(d2,nOsc,label='d2')
call DGEMM_('N','N',nOsc,nOsc,nOsc,One,C1,nOsc,W,nOsc,Zero,A1,nOsc)
call DGEMM_('T','T',nOsc,nOsc,nOsc,One,W1,nOsc,C,nOsc,Zero,temp,nOsc)
B1(:,:) = A1-temp

const = sqr(8)
call DGEMM_('T','N',nOsc,1,nOsc,const,W1,nOsc,r_temp2,nOsc,Zero,d1,nOsc)
call DGEMM_('N','N',nOsc,nOsc,nOsc,One,C2,nOsc,W,nOsc,Zero,A2,nOsc)
call DGEMM_('T','T',nOsc,nOsc,nOsc,One,W2,nOsc,C,nOsc,Zero,temp,nOsc)
B2(:,:) = A2-temp

const = -sqr(8)
call DGEMM_('T','N',nOsc,1,nOsc,const,W2,nOsc,r_temp2,nOsc,Zero,d2,nOsc)

! Calculate A1B1T and A2B2T.
call mma_allocate(A1B1T,nOsc,nOsc,label='A1B1T')
call mma_allocate(A2B2T,nOsc,nOsc,label='A2B2T')

call DGEMM_('N','T',nOsc,nOsc,nOsc,One,A1,nOsc,B1,nOsc,Zero,A1B1T,nOsc)
call DGEMM_('N','T',nOsc,nOsc,nOsc,One,A2,nOsc,B2,nOsc,Zero,A2B2T,nOsc)

call mma_deallocate(temp)
call mma_deallocate(temp1)
call mma_deallocate(temp2)
call mma_deallocate(alpha)
call mma_deallocate(r_temp1)
call mma_deallocate(r_temp2)

! Initialize L matrix.
L(:,:) = Zero
L(0,0) = One

! If max_mOrd > 0 then set up L(m,0).
if (max_mOrd > 0) then
  do kOsc=1,nOsc
    L(mInc(0,kOsc),0) = d1(kOsc)
  end do
  if (max_mInc > 0) then
    do iOrd=1,max_mInc
      kOsc_start = nOsc
      do while ((mMat(iOrd,kOsc_start) == 0) .and. (kOsc_start > 1))
        kOsc_start = kOsc_start-1
      end do
      do kOsc=kOsc_start,nOsc
        do lOsc=1,nOsc
          if (mMat(iOrd,lOsc) > 0) then
            L(mInc(iOrd,kOsc),0) = L(mInc(iOrd,kOsc),0)+sqr(mMat(iOrd,lOsc))*A1B1T(kOsc,lOsc)*L(mDec(iOrd,lOsc),0)
          end if
        end do
        L(mInc(iOrd,kOsc),0) = (L(mInc(iOrd,kOsc),0)+d1(kOsc)*L(iOrd,0))/sqr(mMat(mInc(iOrd,kOsc),kOsc))
      end do
    end do
  end if

  ! Use recursion formula to obtain the rest of L.
  write(u6,*) 'FCVAL Test prints:'
  write(u6,*) 'In this loop, indices jOrd and iOrd will vary.'
  write(u6,*) 'They are used to address L(iOrd,jOrd).'
  write(u6,'(1x,a,2i8)') '(1) jOrd goes from 1, up to max_mord=',max_mord
  write(u6,'(1x,a,2i8)') '(2) iOrd goes from 1, up to max_mord=',max_mord
  !vv stop
  !do jOrd=1,max_nInc2
  do jOrd=1,max_mOrd
    write(u6,'(1x,a,2i8)') ' jOrd=',jOrd
    lOsc = nOsc
    do while ((mMat(jOrd,lOsc) == 0) .and. (lOsc > 1))
      lOsc = lOsc-1
    end do
    do iOrd=0,max_mOrd
      do kOsc=1,nOsc
        if (mMat(iOrd,kOsc) > 0) then
          L(iOrd,jOrd) = L(iOrd,jOrd)+sqr(mMat(iOrd,kOsc))/sqr(mMat(jOrd,lOsc))*A1(kOsc,lOsc)*L(mDec(iOrd,kOsc),mDec(jOrd,lOsc))
        end if
      end do
    end do
  end do
end if

! Initialize U matrix.
U(:,:) = Zero
U(0,0) = One

! If max_nOrd > 0 then set up U(n,0).
if (max_nOrd > 0) then
  do kOsc=1,nOsc
    U(nInc(0,kOsc),0) = d2(kOsc)
  end do
  if (max_nInc > 0) then
    do iOrd=1,max_nInc
      kOsc_start = nOsc
      do while ((nMat(iOrd,kOsc_start) == 0) .and. (kOsc_start > 1))
        kOsc_start = kOsc_start-1
      end do
      do kOsc=kOsc_start,nOsc
        do lOsc=1,nOsc
          if (nMat(iOrd,lOsc) > 0) then

            U(nInc(iOrd,kOsc),0) = U(nInc(iOrd,kOsc),0)+sqr(nMat(iOrd,lOsc))*A2B2T(kOsc,lOsc)*U(nDec(iOrd,lOsc),0)
          end if
        end do

        U(nInc(iOrd,kOsc),0) = (U(nInc(iOrd,kOsc),0)+d2(kOsc)*U(iOrd,0))/sqr(nMat(nInc(iOrd,kOsc),kOsc))
      end do
    end do
  end if

  ! Use recursion formula to obtain the rest of U.
  do jOrd=1,max_nOrd2
    lOsc = nOsc
    do while ((nMat(jOrd,lOsc) == 0) .and. (lOsc > 1))
      lOsc = lOsc-1
    end do
    do iOrd=0,max_nOrd
      do kOsc=1,nOsc
        if (nMat(iOrd,kOsc) > 0) then
          U(iOrd,jOrd) = U(iOrd,jOrd)+sqr(nMat(iOrd,kOsc))/sqr(nMat(jOrd,lOsc))*A2(kOsc,lOsc)*U(nDec(iOrd,kOsc),nDec(jOrd,lOsc))
        end if
      end do
    end do
  end do
end if

call mma_deallocate(A1)
call mma_deallocate(B1)
call mma_deallocate(A2)
call mma_deallocate(B2)
call mma_deallocate(A1B1T)
call mma_deallocate(A2B2T)
call mma_deallocate(d1)
call mma_deallocate(d2)

! Calculate Franck-Condon factors.
FC(:,:) = Zero
do jOrd=0,max_nOrd
  do iOrd=0,max_mOrd
    rsum = Zero
    do kOrd=0,min(max_mOrd,max_nOrd)
      rsum = rsum+L(iOrd,kOrd)*U(jOrd,kOrd)
    end do
    FC(iOrd,jOrd) = FC00*rsum
  end do
end do

call mma_deallocate(sqr)

end subroutine FCval

!end module FCMod

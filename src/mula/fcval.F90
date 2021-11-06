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
!    FCval     (C1,W1,det1,r01,C2,W2,det2,r02,FC,max_mOrd,max_nOrd,
!               max_nOrd2,max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,
!               nInc,mDec,nDec,C,W,det0,r00,L,U,FC00,alpha1,alpha2,beta)
!
!  Uses:
!    TabMod
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

!contains

subroutine FCval(C1,W1,det1,r01,C2,W2,det2,r02,FC,max_mOrd,max_nOrd,max_nOrd2,max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,nInc, &
                 mDec,nDec,C,W,det0,r00,L,U,FC00,alpha1,alpha2,beta,nOsc,nnsiz)
!  Purpose:
!    Calculate multidimensional Franck Condon factors.
!
!  Input:
!    W1,W2      : Real*8 two dimensional arrays - eigenvectors
!                 scaled by the square root of the eigenvalues.
!    C1,C2      : Real*8 two dimensional arrays - inverses
!                 of W1 and W2.
!    det1,det2  : Real*8 variables - determinants of C1 and C2.
!    r01,r02    : Real*8 arrays - coordinates of the two
!                 oscillators.
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
!    W          : Real*8 two dimensional array - eigenvectors
!                 of intermediate oscillator scaled by the square root
!                 of the eigenvalues.
!    C          : Real*8 two dimensional array - inverse
!                 of W.
!    det0       : Real*8 variable - determinant of C.
!    r00        : Real*8 array - coordinates of intermediate
!                 oscillator.
!    L,U        : Real*8 two dimensional arrays
!    FC00       : Real*8 variable - zero-zero overlap.
!    FC         : Real*8
!    alpha1     : Real*8 two dimensional array -
!                 0.5*C1(T)*C1.
!    alpha2     : Real*8 two dimensional array -
!                 0.5*C2(T)*C2.
!    beta       : Real*8 two dimensional array -
!                 0.5*alpha1*alpha^(-1)*alpha2
!
!  Calls:
!    Cholesky    (LinAlg)
!    Dcopy       (Essl)
!    DGEMM_      (Essl)
!    Dool_MULA   (LinAlg)
!    Dscal       (Essl)
!
!  Uses:
!    Linalg
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

implicit real*8(a-h,o-z)
#include "dims.fh"
real*8 FC(0:max_mord,0:max_nord)
real*8 C1(nosc,nosc), C2(nosc,nosc), W1(nosc,nosc)
real*8 W2(nosc,nosc), C(nosc,nosc), W(nosc,nosc)
real*8 L(0:max_mord,0:max_ninc2)
real*8 U(0:max_nord,0:max_nord2)
real*8 alpha1(nosc,nosc), alpha2(nosc,nosc), beta(nosc,nosc)
real*8 r00(nosc), r01(nosc), r02(nosc)
integer mMat(0:mdim1,mdim2), mInc(0:mdim1,mdim2),mDec(0:mdim1,mdim2)
integer nMat(0:ndim1,ndim2), nInc(0:ndim1,ndim2),nDec(0:ndim1,ndim2)
#include "WrkSpc.fh"

! Initialize.
nMaxMat = max(max_mord+1,max_nord+1)
nTabDim = max(nMaxMat,8)
nOscSqr = nOsc**2
call GetMem('temp','Allo','Real',iptemp,nOscSqr)
call GetMem('temp1','Allo','Real',iptemp1,nOscSqr)
call GetMem('temp2','Allo','Real',iptemp2,nOscSqr)

! Setup sqr table.
n = nTabDim+1
call GetMem('sqr','Allo','Real',ipsqr,n+1)
do i=0,nTabDim+1
  Work(ipsqr+i) = sqrt(dble(i))
end do

! Calculate alpha1, alpha2 and alpha.
call GetMem('alpha','Allo','Real',ipalpha,nOscSqr)
call DGEMM_('T','N',nOsc,nOsc,nOsc,1.0d0,C1,nOsc,C1,nOsc,0.0d0,alpha1,nOsc)
call dscal_(nOscSqr,0.5d0,alpha1,1)
call DGEMM_('T','N',nOsc,nOsc,nOsc,1.0d0,C2,nOsc,C2,nOsc,0.0d0,alpha2,nOsc)
call dscal_(nOscSqr,0.5d0,alpha2,1)
!temp = alpha1+alpha2
call dcopy_(nOscSqr,alpha1,1,Work(iptemp),1)
call Daxpy_(nOscSqr,1.0d0,alpha2,1,Work(iptemp),1)
!alpha = 0.5d0*temp
call dcopy_(nOscSqr,[0.0d0],0,Work(ipalpha),1)
call Daxpy_(nOscSqr,0.5d0,Work(iptemp),1,Work(ipalpha),1)

!call xxDgemul(C,nOsc,'T',C,nOsc,'N',alpha,nOsc,nOsc,nOsc,nOsc)
!call dscal_(nOscSqr,0.5d0,alpha,1)

! Calculate C using a Cholesky factorization of 2*alpha.
!call Cholesky(temp,C)

! Calculate W.
!call dcopy_(nOscSqr,[0.0d0],0,W,1)
!call dcopy_(nOsc,[1.0d0],0,W,nOsc+1)
!temp = C
!call Dool_MULA(temp,W,det0)

! Calculate r00.
call GetMem('r_temp1','Allo','Real',ipr_temp1,nOsc)
call GetMem('r_temp2','Allo','Real',ipr_temp2,nOsc)

! Calculate beta.
do i=1,nOsc
  do j=1,nOsc
    Work(iptemp1+j+nOsc*(i-1)-1) = C1(i,j)
  end do
end do
!temp1 = alpha1
call dcopy_(nOscSqr,alpha1,1,Work(iptemp1),1)
!temp = 2.0d0*alpha
call dcopy_(nOscSqr,[0.0d0],0,Work(iptemp),1)
call Daxpy_(nOscSqr,2.0d0,Work(ipalpha),1,Work(iptemp),1)

call Dool_MULA(Work(iptemp),nOsc,nOsc,Work(iptemp1),nOsc,nOsc,det)
call DGEMM_('N','N',nOsc,nOsc,nOsc,1.0d0,alpha2,nOsc,Work(iptemp1),nOsc,0.0d0,beta,nOsc)

! Calculate FC00.
!r_temp1 = r02-r01
call dcopy_(nOsc,r01,1,Work(ipr_temp1),1)
call Daxpy_(nOsc,-1.0d0,r02,1,Work(ipr_temp1),1)

call DGEMM_('N','N',nOsc,1,nOsc,1.0d0,beta,nOsc,Work(ipr_temp1),nOsc,0.0d0,Work(ipr_temp2),nOsc)
FC00_exp = Ddot_(nOsc,Work(ipr_temp1),1,Work(ipr_temp2),1)
FC00 = (sqrt(det1)*sqrt(det2)/det0)*exp(-FC00_exp)

! Calculate A, B and d matrices.
call GetMem('A1','Allo','Real',ipA1,nOscSqr)
call GetMem('B1','Allo','Real',ipB1,nOscSqr)
call GetMem('A2','Allo','Real',ipA2,nOscSqr)
call GetMem('B2','Allo','Real',ipB2,nOscSqr)
call GetMem('d1','Allo','Real',ipd1,nOsc)
call GetMem('d2','Allo','Real',ipd2,nOsc)
call DGEMM_('N','N',nOsc,nOsc,nOsc,1.0d0,C1,nOsc,W,nOsc,0.0d0,Work(ipA1),nOsc)
call DGEMM_('T','T',nOsc,nOsc,nOsc,1.0d0,W1,nOsc,C,nOsc,0.0d0,Work(iptemp),nOsc)
!B1 = A1-temp
call dcopy_(nOscSqr,Work(ipA1),1,Work(ipB1),1)
call Daxpy_(nOscSqr,-1.0d0,Work(iptemp),1,Work(ipB1),1)

call DGEMM_('T','N',nOsc,1,nOsc,1.0d0,W1,nOsc,Work(ipr_temp2),nOsc,0.0d0,Work(ipd1),nOsc)
const = Work(ipsqr+8)
call dscal_(nOsc,const,Work(ipd1),1)
call DGEMM_('N','N',nOsc,nOsc,nOsc,1.0d0,C2,nOsc,W,nOsc,0.0d0,Work(ipA2),nOsc)
call DGEMM_('T','T',nOsc,nOsc,nOsc,1.0d0,W2,nOsc,C,nOsc,0.0d0,Work(iptemp),nOsc)
!B2 = A2-temp
call dcopy_(nOscSqr,Work(ipA2),1,Work(ipB2),1)
call Daxpy_(nOscSqr,-1.0d0,Work(iptemp),1,Work(ipB2),1)

call DGEMM_('T','N',nOsc,1,nOsc,1.0d0,W2,nOsc,Work(ipr_temp2),nOsc,0.0d0,Work(ipd2),nOsc)
const = -Work(ipsqr+8)
call dscal_(nOsc,const,Work(ipd2),1)

! Calculate A1B1T and A2B2T.
call GetMem('A1B1T','Allo','Real',ipA1B1T,nOscSqr)
call GetMem('A2B2T','Allo','Real',ipA2B2T,nOscSqr)

call DGEMM_('N','T',nOsc,nOsc,nOsc,1.0d0,Work(ipA1),nOsc,Work(ipB1),nOsc,0.0d0,Work(ipA1B1T),nOsc)
call DGEMM_('N','T',nOsc,nOsc,nOsc,1.0d0,Work(ipA2),nOsc,Work(ipB2),nOsc,0.0d0,Work(ipA2B2T),nOsc)

call GetMem('temp','Free','Real',iptemp,nOscSqr)
call GetMem('temp1','Free','Real',iptemp1,nOscSqr)
call GetMem('temp2','Free','Real',iptemp2,nOscSqr)
call GetMem('alpha','Free','Real',ipalpha,nOscSqr)
call GetMem('r_temp1','Free','Real',ipr_temp1,nOsc)
call GetMem('r_temp2','Free','Real',ipr_temp2,nOsc)

! Initialize L matrix.
call dcopy_((max_mord+1)*(max_ninc2+1),[0.0d0],0,L,1)
!L = 0.0d0
L(0,0) = 1.0d0

! If max_mOrd > 0 then set up L(m,0).
if (max_mOrd > 0) then
  do kOsc=1,nOsc
    L(mInc(0,kOsc),0) = Work(ipd1+kOsc-1)
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
            L(mInc(iOrd,kOsc),0) = L(mInc(iOrd,kOsc),0)+Work(ipsqr+mMat(iOrd,lOsc))*Work(ipA1B1T+kOsc+nOsc*(lOsc-1)-1)* &
                                   L(mDec(iOrd,lOsc),0)
          end if
        end do
        L(mInc(iOrd,kOsc),0) = (L(mInc(iOrd,kOsc),0)+Work(ipd1+kOsc-1)*L(iOrd,0))/Work(ipsqr+mMat(mInc(iOrd,kOsc),kOsc))
      end do
    end do
  end if

  ! Use recursion formula to obtain the rest of L.
  write(6,*) 'FCVAL Test prints:'
  write(6,*) 'In this loop, indices jOrd and iOrd will vary.'
  write(6,*) 'They are used to address L(iOrd,jOrd).'
  write(6,'(1x,a,2i8)') '(1) jOrd goes from 1, up to max_mord=',max_mord
  write(6,'(1x,a,2i8)') '(2) iOrd goes from 1, up to max_mord=',max_mord
  !vv stop
  !do jOrd=1,max_nInc2
  do jOrd=1,max_mOrd
    write(6,'(1x,a,2i8)') ' jOrd=',jOrd
    lOsc = nOsc
    do while ((mMat(jOrd,lOsc) == 0) .and. (lOsc > 1))
      lOsc = lOsc-1
    end do
    do iOrd=0,max_mOrd
      do kOsc=1,nOsc
        if (mMat(iOrd,kOsc) > 0) then
          L(iOrd,jOrd) = L(iOrd,jOrd)+(Work(ipsqr+mMat(iOrd,kOsc))/Work(ipsqr+mMat(jOrd,lOsc)))*Work(ipA1+kOsc+nOsc*(lOsc-1)-1)* &
                         L(mDec(iOrd,kOsc),mDec(jOrd,lOsc))
        end if
      end do
    end do
  end do
end if

! Initialize U matrix.
!U = 0.0d0
call dcopy_((max_nord+1)*(max_nord2+1),[0.0d0],0,U,1)
U(0,0) = 1.0d0

! If max_nOrd > 0 then set up U(n,0).
if (max_nOrd > 0) then
  do kOsc=1,nOsc
    U(nInc(0,kOsc),0) = Work(ipd2+kOsc-1)
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

            U(nInc(iOrd,kOsc),0) = U(nInc(iOrd,kOsc),0)+Work(ipsqr+nMat(iOrd,lOsc))*Work(ipA2B2T+kOsc+nOsc*(lOsc-1)-1)* &
                                   U(nDec(iOrd,lOsc),0)
          end if
        end do

        U(nInc(iOrd,kOsc),0) = (U(nInc(iOrd,kOsc),0)+WOrk(ipd2+kOsc-1)*U(iOrd,0))/Work(ipsqr+nMat(nInc(iOrd,kOsc),kOsc))
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
          U(iOrd,jOrd) = U(iOrd,jOrd)+(Work(ipsqr+nMat(iOrd,kOsc))/Work(ipsqr+nMat(jOrd,lOsc)))*Work(ipA2+kOsc+nOsc*(lOsc-1)-1)* &
                         U(nDec(iOrd,kOsc),nDec(jOrd,lOsc))
        end if
      end do
    end do
  end do
end if

call GetMem('A1','Free','Real',ipA1,nOscSqr)
call GetMem('B1','Free','Real',ipB1,nOscSqr)
call GetMem('A2','Free','Real',ipA2,nOscSqr)
call GetMem('B2','Free','Real',ipB2,nOscSqr)
call GetMem('A1B1T','Free','Real',ipA1B1T,nOscSqr)
call GetMem('A2B2T','Free','Real',ipA2B2T,nOscSqr)
call GetMem('d1','Free','Real',ipd1,nOsc)
call GetMem('d2','Free','Real',ipd2,nOsc)

! Calculate Franck-Condon factors.
!FC = 0.0d0
call dcopy_((max_mord+1)*(max_nord+1),[0.0d0],0,FC,1)
do jOrd=0,max_nOrd
  do iOrd=0,max_mOrd
    sum = 0.0d0
    do kOrd=0,min(max_mOrd,max_nOrd)

      sum = sum+L(iOrd,kOrd)*U(jOrd,kOrd)
    end do
    FC(iOrd,jOrd) = FC00*sum
  end do
end do

call GetMem('sqr','Free','Real',ipsqr,n+1)

! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(r00)
  call Unused_integer(nnsiz)
end if

end subroutine FCval

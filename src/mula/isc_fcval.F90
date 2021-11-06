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
! Copyright (C) 2008,2009, Giovanni Ghigo                              *
!***********************************************************************

!  InterSystem Crossing rate evaluation: Multidimensional Franck-Condon
!  Modified copy of  FCval by Giovanni Ghigo.
!  Dip. Chimica Generale e Chimica Organica, Torino (ITALY)
!  28-Dec-08 - 06-Jan-09 ; June 2009

subroutine ISC_FCval(iPrint,iMaxYes,nTabDim,C1,W1,det1,r01,C2,W2,det2,r02,max_mOrd,max_nOrd,max_nOrd2,max_mInc,max_nInc,max_nInc2, &
                     mMat,nMat,mInc,nInc,mDec,nDec,C,W,det0,r00,L,U,FC00,Alpha1,Alpha2,Beta,nOsc,nnsiz,iMx_nOrd,nYes,VibWind2, &
                     FCWind2)

implicit real*8(a-h,o-z)
implicit integer(i-n)
#include "dims.fh"
real*8 C1(nOsc,nOsc), C2(nOsc,nOsc), W1(nOsc,nOsc)
real*8 W2(nOsc,nOsc), C(nOsc,nOsc), W(nOsc,nOsc)
real*8 L(0:max_mOrd,0:max_nInc2)
real*8 U(0:max_nOrd,0:max_nOrd2)
real*8 Alpha1(nOsc,nOsc), Alpha2(nOsc,nOsc), Beta(nOsc,nOsc)
real*8 r00(nOsc), r01(nOsc), r02(nOsc)
real*8 FCWind2(nYes)
integer mMat(0:mdim1,mdim2), mInc(0:mdim1,mdim2), mDec(0:mdim1,mdim2)
integer nMat(0:ndim1,ndim2), nInc(0:iMaxYes,nOsc), nDec(0:iMaxYes,nOsc)
integer VibWind2(nYes)
#include "WrkSpc.fh"
#include "inout.fh"

!write(6,*) 'CGGt[ISC_FCval] Enter '
!write(6,*) '                iMx_nOrd, nYes = ',iMx_nOrd, nYes
!write(6,*) '                VibWind2 :',(VibWind2(i),i=1,nYes)
!write(6,*) '            L matrix:',max_mOrd,max_nInc2
!write(6,*) '            U matrix:',max_nOrd,max_nOrd2
!do i=0,iMaxYes
!  write(6,*) (nInc(i,j),j=1,nOsc)
!end do
!call XFlush(6)

! Initialize.
nMaxMat = max(max_mOrd+1,max_nOrd+1)
!write(6,*) '            nMaxMat=',nMaxMat
nTabDim = max(nMaxMat,8)
!write(6,*) '            nTabDim=',nTabDim
nOscSqr = nOsc**2
!write(6,*) '            nOscSqr=',nOscSqr
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
!Write(6,*) 'CGGt[FCVal] Calculate alpha(s)'
!Call XFlush(6)
call GetMem('Alpha','Allo','Real',ipAlpha,nOscSqr)
call DGEMM_('T','N',nOsc,nOsc,nOsc,1.0d0,C1,nOsc,C1,nOsc,0.0d0,Alpha1,nOsc)
call dscal_(nOscSqr,0.5d0,Alpha1,1)
call DGEMM_('T','N',nOsc,nOsc,nOsc,1.0d0,C2,nOsc,C2,nOsc,0.0d0,Alpha2,nOsc)
call dscal_(nOscSqr,0.5d0,Alpha2,1)
!temp = alpha1+alpha2
call dcopy_(nOscSqr,Alpha1,1,Work(iptemp),1)
call Daxpy_(nOscSqr,1.0d0,Alpha2,1,Work(iptemp),1)
!alpha = 0.5d0*temp
call dcopy_(nOscSqr,[0.0d0],0,Work(ipAlpha),1)
call Daxpy_(nOscSqr,0.5d0,Work(iptemp),1,Work(ipAlpha),1)

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
!write(6,*) 'CGGt[FCVal] Calculate beta.'
!call XFlush(6)
do i=1,nOsc
  do j=1,nOsc
    Work(iptemp1+j+nOsc*(i-1)-1) = C1(i,j)
  end do
  !write(6,*) 'CGGt C1(',i,',j)=',(C1(i,jj),jj=1,nOsc)
end do
!call XFlush(6)
!temp1 = alpha1
call dcopy_(nOscSqr,Alpha1,1,Work(iptemp1),1)
!temp  = 2.0d0*alpha
call dcopy_(nOscSqr,[0.0d0],0,Work(iptemp),1)
call Daxpy_(nOscSqr,2.0d0,Work(ipAlpha),1,Work(iptemp),1)

call Dool_MULA(Work(iptemp),nOsc,nOsc,Work(iptemp1),nOsc,nOsc,det)
call DGEMM_('N','N',nOsc,nOsc,nOsc,1.0d0,Alpha2,nOsc,Work(iptemp1),nOsc,0.0d0,Beta,nOsc)

! Calculate FC00.
!r_temp1 = r02-r01
call dcopy_(nOsc,r01,1,Work(ipr_temp1),1)
call Daxpy_(nOsc,-1.0d0,r02,1,Work(ipr_temp1),1)

call DGEMM_('N','N',nOsc,1,nOsc,1.0d0,Beta,nOsc,Work(ipr_temp1),nOsc,0.0d0,Work(ipr_temp2),nOsc)
FC00_exp = Ddot_(nOsc,Work(ipr_temp1),1,Work(ipr_temp2),1)
FC00 = (sqrt(det1)*sqrt(det2)/det0)*exp(-FC00_exp)
!write(6,*) 'CGGt[FCVal] FC00_exp,FC00=',FC00_exp,FC00
!call XFlush(6)

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

! Initialize L matrix.
call dcopy_((max_mord+1)*(max_ninc2+1),[0.0d0],0,L,1)
!L = 0.0d0
L(0,0) = 1.0d0

! If max_mOrd > 0 then set up L(m,0).
if (max_mOrd > 0) then
  write(6,*) '*****************************************'
  write(6,*) ' Hot initial state not implemented yet !'
  write(6,*) '*****************************************'
  write(6,*)
  call Quit_OnUserError()
end if

! Initialize U matrix.
!write(6,*) 'CGGt[FCVal] Initialize U matrix.'
!call XFlush(6)
call dcopy_((max_nord+1)*(max_nOrd2+1),[0.0d0],0,U,1)
U(0,0) = 1.0d0
!GGt -------------------------------------------------------------------
!do kOsc=1,nOsc
!  write(6,*) 'CGGt[FCVal] Work(ipd2..)=',Work(ipd2+kOsc-1)
!  call XFlush(6)
!end do
!GGt -------------------------------------------------------------------

! If max_nOrd > 0 then set up U(n,0).
!write(6,*) 'CGGt[FCVal] max_nOrd > 0 then set up U(n,0).'
!call XFlush(6)
do kOsc=1,nOsc
  !write(6,*) 'CGGt[FCVal] nInc(0,kOsc)=',nInc(0,kOsc)
  !call XFlush(6)
  !write(6,*) 'CGGt[FCVal] Work(ipd2..)=',Work(ipd2+kOsc-1)
  !call XFlush(6)
  U(nInc(0,kOsc),0) = Work(ipd2+kOsc-1)
end do
!write(6,*) 'CGGt[FCVal] max_nInc=',max_nInc
!write(6,*) 'CGGt[FCVal]  iMaxYes=',iMaxYes
!call XFlush(6)
max_nInc = min(max_nInc,iMaxYes)
if (max_nInc > 0) then
  do iOrd=1,max_nInc
    !write(6,*) '              iOrd=',iOrd
    !call XFlush(6)
    kOsc_start = nOsc
    do while ((nMat(iOrd,kOsc_start) == 0) .and. (kOsc_start > 1))
      kOsc_start = kOsc_start-1
    end do
    do kOsc=kOsc_start,nOsc
      do lOsc=1,nOsc
        !write(6,*) 'iOrd,kOsc_start,kOsc,lOsc==',iOrd,kOsc_start,kOsc,lOsc
        !call XFlush(6)
        if (nMat(iOrd,lOsc) > 0) then

          U(nInc(iOrd,kOsc),0) = U(nInc(iOrd,kOsc),0)+Work(ipsqr+nMat(iOrd,lOsc))*Work(ipA2B2T+kOsc+nOsc*(lOsc-1)-1)* &
                                 U(nDec(iOrd,lOsc),0)
        end if
      end do

      U(nInc(iOrd,kOsc),0) = (U(nInc(iOrd,kOsc),0)+Work(ipd2+kOsc-1)*U(iOrd,0))/Work(ipsqr+nMat(nInc(iOrd,kOsc),kOsc))
    end do
  end do
end if

! Use recursion formula to obtain the rest of U.
!write(6,*) '            Use recursion ... rest of U.'
!write(6,*) '            max_nOrd2=',max_nOrd2
!call XFlush(6)
do jOrd=1,max_nOrd2
  lOsc = nOsc
  do while ((nMat(jOrd,lOsc) == 0) .and. (lOsc > 1))
    lOsc = lOsc-1
  end do
  do iOrd=0,max_nOrd
    do kOsc=1,nOsc
      if (nMat(iOrd,kOsc) > 0) then
        !write(6,*) '              ',iOrd,kOsc,nMat(iOrd,kOsc)
        U(iOrd,jOrd) = U(iOrd,jOrd)+(Work(ipsqr+nMat(iOrd,kOsc))/Work(ipsqr+nMat(jOrd,lOsc)))*Work(ipA2+kOsc+nOsc*(lOsc-1)-1)* &
                       U(nDec(iOrd,kOsc),nDec(jOrd,lOsc))
      end if
    end do
  end do
end do

call GetMem('A2B2T','Free','Real',ipA2B2T,nOscSqr)
call GetMem('A1B1T','Free','Real',ipA1B1T,nOscSqr)
call GetMem('d2','Free','Real',ipd2,nOsc)
call GetMem('d1','Free','Real',ipd1,nOsc)
call GetMem('B2','Free','Real',ipB2,nOscSqr)
call GetMem('A2','Free','Real',ipA2,nOscSqr)
call GetMem('B1','Free','Real',ipB1,nOscSqr)
call GetMem('A1','Free','Real',ipA1,nOscSqr)

! Calculate Franck-Condon factors.
if (iPrint >= 3) then
  write(6,*) ' Franck-Condon factors for States in the Window: '
  write(6,'(a,36a)') '  ',('=',i=1,36)
  write(6,*) '     #     jOrd   FC factor     jSum '
  write(6,'(a,36a)') '  ',('-',i=1,36)
end if
do ii=1,nYes
  jOrd = VibWind2(ii)
  dFC = FC00*L(0,0)*U(jOrd,0)
  FCWind2(ii) = dFC
  if (iPrint >= 3) then
    loc_n_max = 0
    do j=1,nOsc
      loc_n_max = loc_n_max+nMat(jOrd,j)
    end do
    write(6,'(a2,i5,i9,e15.6,a2,i4)') ' ',ii,jOrd,FCWind2(ii),' ',loc_n_max
  end if
end do
if (iPrint >= 3) then
  write(6,'(a,36a)') '  ',('-',i=1,36)
  write(6,*) ' FC_00 =',FC00
  write(6,*)
end if

if (iPrint >= 4) then
  write(6,*)
  write(6,*) ' Full Franck-Condon factors (FC_00=',FC00,'):'
  write(6,*) ' =================================================='
  write(6,*) '    jOrd   FC            level                     '
  write(6,*) ' --------------------------------------------------'
  do jOrd=0,max_nOrd
    loc_n_max = 0
    do j=1,nOsc
      loc_n_max = loc_n_max+nMat(jOrd,j)
    end do
    dFC = FC00*L(0,0)*U(jOrd,0)
    write(6,'(a,i8,e15.6,a2,i4,a2,24i3)') ' ',jOrd,dFC,' ',loc_n_max,' ',(nMat(jOrd,j),j=1,nOsc)
  end do
  write(6,*) ' --------------------------------------------------'
end if
!GGt -------------------------------------------------------------------

call GetMem('r_temp2','Free','Real',ipr_temp2,nOsc)
call GetMem('r_temp1','Free','Real',ipr_temp1,nOsc)
call GetMem('Alpha','Free','Real',ipAlpha,nOscSqr)
call GetMem('sqr','Free','Real',ipsqr,n+1)
call GetMem('temp2','Free','Real',iptemp2,nOscSqr)
call GetMem('temp1','Free','Real',iptemp1,nOscSqr)
call GetMem('temp','Free','Real',iptemp,nOscSqr)

!write(6,*) 'CGGt[FCVal] Exit'
!call XFlush(6)

! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(max_mInc)
  call Unused_integer_array(mMat)
  call Unused_integer_array(mInc)
  call Unused_integer_array(mDec)
  call Unused_real_array(r00)
  call Unused_integer(nnsiz)
  call Unused_integer(iMx_nOrd)
end if

end subroutine ISC_FCval

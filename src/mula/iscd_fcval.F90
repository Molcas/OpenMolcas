!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine ISCD_FCval(iPrint,iMaxYes,lnTabDim,nnTabDim,lNMAT0,lNMAT,lNINC,lNDEC,lBatch,nBatch,leftBatch,nIndex,C1,W1,det1,r01,C2, &
                      W2,det2,r02,max_mOrd,max_nOrd,max_nOrd2,max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,nInc,mDec,nDec,C,W,det0, &
                      r00,L,U,FC00,Alpha1,Alpha2,Beta,nOsc,nnsiz,iMx_nOrd,nYes,VibWind2,FCWind2,nMat0)

use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, u6

implicit real*8(a-h,o-z)
implicit integer(i-n)
#include "dims.fh"
#include "io_mula.fh"
integer nIndex(3,0:maxMax_n)
real*8 C1(nOsc,nOsc), C2(nOsc,nOsc), W1(nOsc,nOsc)
real*8 W2(nOsc,nOsc), C(nOsc,nOsc), W(nOsc,nOsc)
real*8 L(0:max_mOrd,0:max_nInc2)
real*8 U(0:max_nOrd,0:max_nOrd2)
real*8 Alpha1(nOsc,nOsc), Alpha2(nOsc,nOsc), Beta(nOsc,nOsc)
real*8 r00(nOsc), r01(nOsc), r02(nOsc)
real*8 FCWind2(nYes)
integer mMat(0:mdim1,mdim2), mInc(0:mdim1,mdim2), mDec(0:mdim1,mdim2)
integer nMat(nOsc,lBatch), nInc(nOsc,lBatch), nDec(nOsc,lBatch)
integer VibWind2(nYes), nMat0(nOsc), nnTabDim(0:lnTabDim)
#include "WrkSpc.fh"
#include "inout.fh"

!GGt -------------------------------------------------------------------
!write(u6,*) 'CGGt[ISCD_FCval] Enter '
!write(u6,*) '     iMx_nOrd, nYes = ',iMx_nOrd, nYes
!write(u6,*) '     VibWind2 :',(VibWind2(i),i=1,nYes)
!write(u6,*) '     L matrix:',max_mOrd,max_nInc2
!write(u6,*) '     U matrix:',max_nOrd,max_nOrd2
!do i=0,iMaxYes
!  write(u6,*) (nInc(i,j),j=1,nOsc)
!end do
!write(u6,*) '     lnTabDim+1=',lnTabDim+1,':'
!do i=0,lnTabDim
!  iIndex0 = nnTabDim(i)
!  !call iDaFile(lNMAT0,2,nMat0,nOsc,iIndex0)
!  write(u6,*) i,' read at',nnTabDim(i),'  M:',(nMat0(j),j=1,nOsc)
!end do
!write(u6,*) '-----------------------------------------------'
!call XFlush(u6)
!GGt -------------------------------------------------------------------

! Initialize.
rewind(lnMAT0)
rewind(lnMAT)
rewind(lnINC)
rewind(lnDEC)
nMaxMat = max(max_mOrd+1,max_nOrd+1)
!write(u6,*) '            nMaxMat=',nMaxMat
nTabDim = max(nMaxMat,8)
!write(u6,*) '            nTabDim=',nTabDim
nOscSqr = nOsc**2
!write(u6,*) '            nOscSqr=',nOscSqr
call GetMem('temp','Allo','Real',iptemp,nOscSqr)
call GetMem('temp1','Allo','Real',iptemp1,nOscSqr)
call GetMem('temp2','Allo','Real',iptemp2,nOscSqr)

! Setup sqr table.
n = nTabDim+1
call GetMem('sqr','Allo','Real',ipsqr,n+1)
do i=0,nTabDim+1
  Work(ipsqr+i) = sqrt(real(i,kind=wp))
end do

! Calculate alpha1, alpha2 and alpha.
!write(u6,*) 'CGGt[FCVal] Calculate alpha(s)'
!call XFlush(u6)
call GetMem('Alpha','Allo','Real',ipAlpha,nOscSqr)
call DGEMM_('T','N',nOsc,nOsc,nOsc,One,C1,nOsc,C1,nOsc,Zero,Alpha1,nOsc)
call dscal_(nOscSqr,Half,Alpha1,1)
call DGEMM_('T','N',nOsc,nOsc,nOsc,One,C2,nOsc,C2,nOsc,Zero,Alpha2,nOsc)
call dscal_(nOscSqr,Half,Alpha2,1)
!temp = alpha1+alpha2
call dcopy_(nOscSqr,Alpha1,1,Work(iptemp),1)
call Daxpy_(nOscSqr,One,Alpha2,1,Work(iptemp),1)
!alpha = Half*temp
call dcopy_(nOscSqr,[Zero],0,Work(ipAlpha),1)
call Daxpy_(nOscSqr,Half,Work(iptemp),1,Work(ipAlpha),1)

!call xxDgemul(C,nOsc,'T',C,nOsc,'N',alpha,nOsc,nOsc,nOsc,nOsc)
!call dscal_(nOscSqr,Half,alpha,1)

! Calculate C using a Cholesky factorization of 2*alpha.
!call Cholesky(temp,C)

! Calculate W.
!call dcopy_(nOscSqr,[Zero],0,W,1)
!call dcopy_(nOsc,[One],0,W,nOsc+1)
!temp = C
!call Dool_MULA(temp,W,det0)

! Calculate r00.
call GetMem('r_temp1','Allo','Real',ipr_temp1,nOsc)
call GetMem('r_temp2','Allo','Real',ipr_temp2,nOsc)

! Calculate beta.
!write(u6,*) 'CGGt[FCVal] Calculate beta.'
!call XFlush(u6)
do i=1,nOsc
  do j=1,nOsc
    Work(iptemp1+j+nOsc*(i-1)-1) = C1(i,j)
  end do
  !write(u6,*) 'CGGt C1(',i,',j)=',(C1(i,jj),jj=1,nOsc)
end do
!call XFlush(u6)
!temp1 = alpha1
call dcopy_(nOscSqr,Alpha1,1,Work(iptemp1),1)
!temp  = Two*alpha
call dcopy_(nOscSqr,[Zero],0,Work(iptemp),1)
call Daxpy_(nOscSqr,Two,Work(ipAlpha),1,Work(iptemp),1)

call Dool_MULA(Work(iptemp),nOsc,nOsc,Work(iptemp1),nOsc,nOsc,det)
call DGEMM_('N','N',nOsc,nOsc,nOsc,One,Alpha2,nOsc,Work(iptemp1),nOsc,Zero,Beta,nOsc)

! Calculate FC00.
!r_temp1 = r02-r01
call dcopy_(nOsc,r01,1,Work(ipr_temp1),1)
call Daxpy_(nOsc,-One,r02,1,Work(ipr_temp1),1)

call DGEMM_('N','N',nOsc,1,nOsc,One,Beta,nOsc,Work(ipr_temp1),nOsc,Zero,Work(ipr_temp2),nOsc)
FC00_exp = Ddot_(nOsc,Work(ipr_temp1),1,Work(ipr_temp2),1)
FC00 = (sqrt(det1)*sqrt(det2)/det0)*exp(-FC00_exp)
!write(u6,*) 'CGGt[FCVal] FC00_exp,FC00=',FC00_exp,FC00
!call XFlush(u6)

! Calculate A, B and d matrices.
call GetMem('A1','Allo','Real',ipA1,nOscSqr)
call GetMem('B1','Allo','Real',ipB1,nOscSqr)
call GetMem('A2','Allo','Real',ipA2,nOscSqr)
call GetMem('B2','Allo','Real',ipB2,nOscSqr)
call GetMem('d1','Allo','Real',ipd1,nOsc)
call GetMem('d2','Allo','Real',ipd2,nOsc)
call DGEMM_('N','N',nOsc,nOsc,nOsc,One,C1,nOsc,W,nOsc,Zero,Work(ipA1),nOsc)
call DGEMM_('T','T',nOsc,nOsc,nOsc,One,W1,nOsc,C,nOsc,Zero,Work(iptemp),nOsc)
!B1 = A1-temp
call dcopy_(nOscSqr,Work(ipA1),1,Work(ipB1),1)
call Daxpy_(nOscSqr,-One,Work(iptemp),1,Work(ipB1),1)

call DGEMM_('T','N',nOsc,1,nOsc,One,W1,nOsc,Work(ipr_temp2),nOsc,Zero,Work(ipd1),nOsc)
const = Work(ipsqr+8)
call dscal_(nOsc,const,Work(ipd1),1)
call DGEMM_('N','N',nOsc,nOsc,nOsc,One,C2,nOsc,W,nOsc,Zero,Work(ipA2),nOsc)
call DGEMM_('T','T',nOsc,nOsc,nOsc,One,W2,nOsc,C,nOsc,Zero,Work(iptemp),nOsc)
!B2 = A2-temp
call dcopy_(nOscSqr,Work(ipA2),1,Work(ipB2),1)
call Daxpy_(nOscSqr,-One,Work(iptemp),1,Work(ipB2),1)

call DGEMM_('T','N',nOsc,1,nOsc,One,W2,nOsc,Work(ipr_temp2),nOsc,Zero,Work(ipd2),nOsc)
const = -Work(ipsqr+8)
call dscal_(nOsc,const,Work(ipd2),1)

! Calculate A1B1T and A2B2T.
call GetMem('A1B1T','Allo','Real',ipA1B1T,nOscSqr)
call GetMem('A2B2T','Allo','Real',ipA2B2T,nOscSqr)

call DGEMM_('N','T',nOsc,nOsc,nOsc,One,Work(ipA1),nOsc,Work(ipB1),nOsc,Zero,Work(ipA1B1T),nOsc)
call DGEMM_('N','T',nOsc,nOsc,nOsc,One,Work(ipA2),nOsc,Work(ipB2),nOsc,Zero,Work(ipA2B2T),nOsc)

! Initialize L matrix.
call dcopy_((max_mord+1)*(max_ninc2+1),[Zero],0,L,1)
!L = Zero
L(0,0) = One

! If max_mOrd > 0 then set up L(m,0).
if (max_mOrd > 0) then
  write(u6,*) '*****************************************'
  write(u6,*) ' Hot initial state not implemented yet !'
  write(u6,*) '*****************************************'
  write(u6,*)
  call Quit_OnUserError()
end if

! Initialize U matrix.
!write(u6,*) 'CGGt[FCVal] Initialize U matrix.'
!call XFlush(u6)
call dcopy_((max_nord+1)*(max_nOrd2+1),[Zero],0,U,1)
U(0,0) = One
do kOsc=1,nOsc
!write(u6,*) 'CGGt[FCVal] Work(ipd2..)=',Work(ipd2+kOsc-1)
!call XFlush(u6)
end do

! Reading for first batch nMat, nInc, nDec

iIndex = nIndex(1,1)
jIndex = nIndex(2,1)
kIndex = nIndex(3,1)
!write(u6,*) 'CGGt[FCVal] Reading nMat at ',iIndex
!call XFlush(u6)
call iDaFile(lNMAT,2,nMat,nOsc*lBatch,iIndex)
!GGt -------------------------------------------------------------------
!write(u6,*) 'CGGt-nMat:'
!do i=0,lBatch-1
!  write(u6,*) i,'-M:',(nMat(j,i+1),j=1,nOsc)
!end do
!write(u6,*) '-----------------------------------------------'
!GGt -------------------------------------------------------------------

!write(u6,*) 'CGGt[FCVal] Reading nInc at ',jIndex
!call XFlush(u6)
call iDaFile(lNINC,2,nInc,nOsc*lBatch,jIndex)
!GGt -------------------------------------------------------------------
!write(u6,*) 'CGGt-nInc:'
!do i=0,lBatch-1
!  write(u6,*) i,'-I:',(nInc(j,i+1),j=1,nOsc)
!end do
!write(u6,*) '-----------------------------------------------'
!GGt -------------------------------------------------------------------

!write(u6,*) 'CGGt[FCVal] Reading nDec at ',kIndex
!call XFlush(u6)
call iDaFile(lNDEC,2,nDec,nOsc*lBatch,kIndex)
!GGt -------------------------------------------------------------------
!write(u6,*) 'CGGt-nDec:'
!do i=0,lBatch-1
!  write(u6,*) i,'-D:',(nDec(j,i+1),j=1,nOsc)
!end do
!write(u6,*) '-----------------------------------------------'
!GGt -------------------------------------------------------------------

! If max_nOrd > 0 then set up U(n,0).

!GGt -------------------------------------------------------------------
!write(u6,*) 'CGGt[FCVal] max_nOrd > 0 then set up U(n,0).'
!Call XFlush(u6)
!GGt -------------------------------------------------------------------
do kOsc=1,nOsc
  U(nInc(kOsc,1),0) = Work(ipd2+kOsc-1)
end do
max_nInc = min(max_nInc,iMaxYes)
!GGt -------------------------------------------------------------------
!do kOsc=1,nOsc
!  write(u6,*) 'CGGt[FCVal] U(',nInc(kOsc,1),',0) =',U(nInc(kOsc,1),0)
!end do
!write(u6,*) 'CGGt[FCVal] max_nInc=',max_nInc
!write(u6,*) 'CGGt[FCVal]  iMaxYes=',iMaxYes
!call XFlush(u6)
!GGt -------------------------------------------------------------------

iBatch = 1
if (max_nInc > 0) then
  do iOrd=1,max_nInc
    !write(u6,*) '              iOrd=',iOrd,'    iBatch=',iBatch
    !call XFlush(u6)
    iiOrd = iOrd-(iBatch-1)*lBatch+1
    if (iiOrd == (lBatch+1)) then
      iBatch = iBatch+1
      iIndex = nIndex(1,iBatch)
      jIndex = nIndex(2,iBatch)
      kIndex = nIndex(3,iBatch)
      !write(u6,*) 'CGGt[] iBatch=',iBatch,': Reading nMat at ',iIndex
      !call XFlush(u6)
      call iDaFile(lNMAT,2,nMat,nOsc*lBatch,iIndex)
      !write(u6,*) 'CGGt[] iBatch=',iBatch,': Reading nInc at ',iIndex
      !call XFlush(u6)
      call iDaFile(lNINC,2,nInc,nOsc*lBatch,jIndex)
      !write(u6,*) 'CGGt[] iBatch=',iBatch,': Reading nDec at ',iIndex
      !call XFlush(u6)
      call iDaFile(lNDEC,2,nDec,nOsc*lBatch,kIndex)
      iiOrd = iOrd-(iBatch-1)*lBatch+1
    end if
      !write(u6,*) '              iOrd=',iOrd,'  iiOrd=',iiOrd
      !write(u6,*) iOrd,'-MID:',(nMat(j,iiOrd),j=1,nOsc),' /',(nInc(j,iiOrd),j=1,nOsc),' /',(nDec(j,iiOrd),j=1,nOsc)
      !call XFlush(u6)
    kOsc_start = nOsc
    do while ((nMat(kOsc_start,iiOrd) == 0) .and. (kOsc_start > 1))
      kOsc_start = kOsc_start-1
    end do
    do kOsc=kOsc_start,nOsc
      do lOsc=1,nOsc
        !write(u6,*) 'iOrd,kOsc_start,kOsc,lOsc==',iOrd,kOsc_start,kOsc,lOsc
        !call XFlush(u6)
        if (nMat(lOsc,iiOrd) > 0) then
          U(nInc(kOsc,iiOrd),0) = U(nInc(kOsc,iiOrd),0)+Work(ipsqr+nMat(lOsc,iiOrd))*Work(ipA2B2T+kOsc+nOsc*(lOsc-1)-1)* &
                                  U(nDec(lOsc,iiOrd),0)
          !write(u6,*) iOrd,' U(',nInc(kOsc,iiOrd),')=',U(nInc(kOsc,iiOrd),0)
          !call XFlush(u6)
        end if
      end do
      !write(u6,*) iOrd,' jOrd=',nInc(kOsc,iiOrd)
      jOrd = nInc(kOsc,iiOrd)
      jjOrd = jOrd-(iBatch-1)*lBatch+1
      if ((jjOrd < 1) .or. (jjOrd > lBatch)) then
        iIndex0 = nnTabDim(jOrd)
        call iDaFile(lNMAT0,2,nMat0,nOsc,iIndex0)
        !write(u6,*) jOrd,' read at',nnTabDim(jOrd),'  M:',(nMat0(j),j=1,nOsc),'  jjOrd=',jjOrd
        kDelta = nMat0(kOsc)
      else
        !GGn kDelta = nMat(kOsc,nInc(kOsc,iiOrd)+1)
        kDelta = nMat(kOsc,jjOrd)
      end if
      !Write(u6,*) '         ',iOrd,' >',kDelta
      !Write(u6,*) iOrd,' jOrd=',jOrd  ,'  kDelta=',kDelta
      U(nInc(kOsc,iiOrd),0) = (U(nInc(kOsc,iiOrd),0)+Work(ipd2+kOsc-1)*U(iOrd,0))/Work(ipsqr+kDelta) ! nMat(kOsc,nInc(kOsc,iiOrd)))
    end do
  end do
end if

! Use recursion formula to obtain the rest of U.
!write(u6,*) '            Use recursion ... rest of U.'
!write(u6,*) '            max_nOrd2=',max_nOrd2
!call XFlush(u6)
!!do jOrd=1,max_nOrd2
!!  lOsc = nOsc
!!  do while ((nMat(jOrd,lOsc) == 0) .and. (lOsc > 1))
!!    lOsc = lOsc-1
!!  end do
!!  do iOrd=0,max_nOrd
!!    do kOsc=1,nOsc
!!      if (nMat(iOrd,kOsc) > 0) then
!!        !write(u6,*) '              ',iOrd,kOsc,nMat(iOrd,kOsc)
!!        U(iOrd,jOrd) = U(iOrd,jOrd)+(Work(ipsqr+nMat(iOrd,kOsc))/Work(ipsqr+nMat(jOrd,lOsc)))*Work(ipA2+kOsc+nOsc*(lOsc-1)-1)* &
!!                       U(nDec(iOrd,kOsc),nDec(jOrd,lOsc))
!!      end if
!!    end do
!!  end do
!!end do

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
  write(u6,*) ' Franck-Condon factors for States in the Window: '
  write(u6,'(a,36a)') '  ',('=',i=1,36)
  write(u6,*) '     #     jOrd   FC factor     jSum '
  write(u6,'(a,36a)') '  ',('-',i=1,36)
end if
do ii=1,nYes
  jOrd = VibWind2(ii)
  dFC = FC00*L(0,0)*U(jOrd,0)
  FCWind2(ii) = dFC
  if (iPrint >= 3) then
    loc_n_max = 0
    kIndex = nnTabDim(jOrd)
    call iDaFile(lNMAT0,2,nMat0,nOsc,kIndex)
    do j=1,nOsc
      loc_n_max = loc_n_max+nMat0(j)
    end do
    write(u6,'(a2,i5,i9,e15.6,a2,i4)') ' ',ii,jOrd,FCWind2(ii),' ',loc_n_max
  end if
end do
if (iPrint >= 3) then
  write(u6,'(a,36a)') '  ',('-',i=1,36)
  write(u6,*) ' FC_00 =',FC00
  write(u6,*)
end if

if (iPrint >= 4) then
  write(u6,*)
  write(u6,*) ' Full Franck-Condon factors (FC_00=',FC00,'):'
  write(u6,*) ' =================================================='
  write(u6,*) '    jOrd   FC            level                     '
  write(u6,*) ' --------------------------------------------------'
  do jOrd=0,max_nInc ! max_nOrd
    loc_n_max = 0
    kIndex = nnTabDim(jOrd)
    call iDaFile(lNMAT0,2,nMat0,nOsc,kIndex)
    do j=1,nOsc
      loc_n_max = loc_n_max+nMat0(j)
    end do
    dFC = FC00*L(0,0)*U(jOrd,0)
    write(u6,'(a,i8,e15.6,a2,i4,a2,24i3)') ' ',jOrd,dFC,' ',loc_n_max,' ',(nMat0(j),j=1,nOsc)
  end do
  write(u6,*) ' --------------------------------------------------'
end if

call GetMem('r_temp2','Free','Real',ipr_temp2,nOsc)
call GetMem('r_temp1','Free','Real',ipr_temp1,nOsc)
call GetMem('Alpha','Free','Real',ipAlpha,nOscSqr)
call GetMem('sqr','Free','Real',ipsqr,n+1)
call GetMem('temp2','Free','Real',iptemp2,nOscSqr)
call GetMem('temp1','Free','Real',iptemp1,nOscSqr)
call GetMem('temp','Free','Real',iptemp,nOscSqr)

!write(u6,*) 'CGGt[FCVal] Exit'
!call XFlush(u6)

! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(nBatch)
  call Unused_integer(leftBatch)
  call Unused_integer(max_mInc)
  call Unused_integer_array(mMat)
  call Unused_integer_array(mInc)
  call Unused_integer_array(mDec)
  call Unused_real_array(r00)
  call Unused_integer(nnsiz)
  call Unused_integer(iMx_nOrd)
end if

end subroutine ISCD_FCval

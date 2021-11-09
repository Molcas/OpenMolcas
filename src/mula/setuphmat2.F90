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

subroutine SetUpHmat2(energy1,energy2,C,W,det,r1,r2,max_mOrd,max_nOrd,max_nOrd2,max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,nInc, &
                      mDec,nDec,H,S,Hess,G0,Base,rOrigin,nnsiz,nDimTot,nOsc)
!  Purpose:
!
!  Input:
!
!  Output:
!
!  Uses:
!    Linalg
!    OptMod
!    FCMod
!    VibMod
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

use Constants, only: Zero

!use Linalg
!use OptMod
!use FCMod
!use TabMod
implicit real*8(a-h,o-z)
#include "Constants_mula.fh"
#include "dims.fh"
real*8 r1(nosc)
real*8 r2(nosc)
real*8 rOrigin(nosc)
real*8 C(nOsc,nOsc)
real*8 W(nOsc,nOsc)
integer mMat(0:mdim1,mdim2), mInc(0:mdim1,mdim2), mDec(0:mdim1,mdim2)
integer nMat(0:ndim1,ndim2), nInc(0:ndim1,ndim2), nDec(0:ndim1,ndim2)
real*8 H(nDimTot,nDimTot)
real*8 S(nDimTot,nDimTot)
real*8 Hess(nOsc,nOsc)
real*8 G0(nOsc,nOsc)
real*8 Base(nosc,nOsc)
#include "WrkSpc.fh"

! Initialize.
! arrays for setuphmat2
call GetMem('Hij','Allo','Real',ipHij,(max_mOrd+1)*(max_nOrd+1))
call GetMem('Sij','Allo','Real',ipSij,(max_mOrd+1)*(max_nOrd+1))
call GetMem('r0','Allo','Real',ipr0,nOsc)
call GetMem('r_diff','Allo','Real',ipr_diff,nOsc)
call GetMem('alpha1','Allo','Real',ipalpha1,nOsc*nOsc)
call GetMem('alpha2','Allo','Real',ipalpha2,nOsc*nOsc)
call GetMem('beta','Allo','Real',ipbeta,nOsc*nOsc)
call GetMem('L','Allo','Real',ipL,(max_mOrd+1)*(max_mOrd+1))
call GetMem('U','Allo','Real',ipU,(max_nOrd+1)*(max_nOrd+1))
call GetMem('C0','Allo','Real',ipC0,nOsc*nOsc)
call GetMem('W0','Allo','Real',ipW0,nOsc*nOsc)
call GetMem('grad','Allo','Real',ipgrad,nOsc)
call GetMem('D3','Allo','Real',ipD3,nOsc*nOsc*nOsc)
call GetMem('Gprime','Allo','Real',ipGprime,nOsc*nOsc*nOsc)
call GetMem('D4','Allo','Real',ipD4,nOsc*nOsc*nOsc*nOsc)
call GetMem('Gdble','Allo','Real',ipGdblePrime,nOsc*nOsc*nOsc*nOsc)

! - Call Franck-Condon routine.
! - Calculate matrix elements.

max_term = 2
!grad = Zero
call dcopy_(nOsc,[Zero],0,Work(ipgrad),1)
!D3 = Zero
call dcopy_(nOsc*nOsc*nOsc,[Zero],0,Work(ipD3),1)
call dcopy_(nOsc*nOsc*nOsc,[Zero],0,Work(ipGprime),1)
!D4 = Zero
call dcopy_(nOsc*nOsc*nOsc*nOsc,[Zero],0,Work(ipD4),1)
call dcopy_(nOsc*nOsc*nOsc*nOsc,[Zero],0,Work(ipGdblePrime),1)
!Gprime = Zero
!GdblePrime = Zero
!Base = Zero
!do i=1,nOsc
!  Base(i,i) = One
!end do
call Calc_r00(C,C,W,W,Work(ipC0),Work(ipW0),Work(ipalpha1),Work(ipalpha2),Work(ipr0),r1,r1,det0,det,det,FC00,nOsc)
call FCval(C,W,det0,Work(ipr0),C,W,det0,Work(ipr0),Work(ipSij),max_mOrd,max_nOrd,max_nOrd,max_mInc,max_nInc,max_nInc,mMat,nMat, &
           mInc,nInc,mDec,nDec,Work(ipC0),Work(ipW0),det0,Work(ipr0),Work(ipL),Work(ipU),FC00,Work(ipalpha1),Work(ipalpha2), &
           Work(ipbeta),nOsc,nnsiz)
!r_diff = Zero
call dcopy_(nOsc,[Zero],0,Work(ipr_diff),1)
call MatrixElements(Work(ipL),Work(ipU),FC00,Work(ipHij),Work(ipC0),Work(ipW0),Work(ipr_diff),mMat,nMat,nInc,nDec,max_nOrd, &
                    max_mOrd,nOsc,energy1,Work(ipgrad),Hess,Work(ipD3),Work(ipD4),G0,Work(ipGprime),Work(ipGdbleprime), &
                    Work(ipalpha1),Work(ipalpha2),Work(ipbeta),max_term,Base)

!H(1:max_mOrd+1,1:max_nOrd+1) = Hij
!S(1:max_mOrd+1,1:max_nOrd+1) = Sij
n_H = (max_mOrd+1)*(max_nOrd+1)
call dcopy_(n_H,Work(ipHij),1,H,1)
call dcopy_(n_H,Work(ipSij),1,S,1)

call GetMem('Hij','Free','Real',ipHij,(max_mOrd+1)*(max_nOrd+1))
call GetMem('Sij','Free','Real',ipSij,(max_mOrd+1)*(max_nOrd+1))
call GetMem('r0','Free','Real',ipr0,nOsc)
call GetMem('r_diff','Free','Real',ipr_diff,nOsc)
call GetMem('alpha1','Free','Real',ipalpha1,nOsc*nOsc)
call GetMem('alpha2','Free','Real',ipalpha2,nOsc*nOsc)
call GetMem('beta','Free','Real',ipbeta,nOsc*nOsc)
call GetMem('L','Free','Real',ipL,(max_mOrd+1)*(max_mOrd+1))
call GetMem('U','Free','Real',ipU,(max_nOrd+1)*(max_nOrd+1))
call GetMem('C0','Free','Real',ipC0,nOsc*nOsc)
call GetMem('W0','Free','Real',ipW0,nOsc*nOsc)
call GetMem('grad','Free','Real',ipgrad,nOsc)
call GetMem('D3','Free','Real',ipD3,nOsc*nOsc*nOsc)
call GetMem('Gprime','Free','Real',ipGprime,nOsc*nOsc*nOsc)
call GetMem('D4','Free','Real',ipD4,nOsc*nOsc*nOsc*nOsc)
call GetMem('Gdble','Free','Real',ipGdblePrime,nOsc*nOsc*nOsc*nOsc)

! Avoid unused argument warnings
if (.false.) then
  call Unused_real(energy2)
  call Unused_real_array(r2)
  call Unused_integer(max_nInc2)
  call Unused_integer(max_nOrd2)
  call Unused_real_array(rOrigin)
end if

end subroutine SetUpHmat2

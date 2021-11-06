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

subroutine Calc_r00(C1,C2,W1,W2,C,W,alpha1,alpha2,r00,r01,r02,det0,det1,det2,FC00,nOsc)
!  Purpose:
!    Calculate geometry of the intermediate oscillator.

implicit real*8(a-h,o-z)
real*8 C1(nOsc,nOsc), C2(nOsc,nOsc), C(nOsc,nOsc)
real*8 W1(nOsc,nOsc), W2(nOsc,nOsc), W(nOsc,nOsc)
real*8 alpha1(nOsc,nOsc), alpha2(nOsc,nOsc)
real*8 r00(nOsc), r01(nOsc), r02(nOsc)
#include "WrkSpc.fh"

! Initialize.
my1 = nOsc
my2 = nOsc
nOscSqr = nOsc**2

call GetMem('temp','Allo','Real',iptemp,nOscSqr)

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

! Calculate C using a Cholesky factorization of 2*alpha.
call Cholesky(Work(iptemp),C,nOsc)

! Calculate W.
call dcopy_(nOscSqr,[0.0d0],0,W,1)
call dcopy_(nOsc,[1.0d0],0,W,nOsc+1)
call dcopy_(nOscSqr,C,1,Work(iptemp),1)
!temp = C
call Dool_MULA(Work(iptemp),nOsc,nOsc,W,my1,my2,det0)
det0 = abs(det0)

! Calculate r00.
call GetMem('r_temp1','Allo','Real',ipr_temp1,nOsc)
call GetMem('r_temp2','Allo','Real',ipr_temp2,nOsc)
call GetMem('r_temp','Allo','Real',ipr_temp,nOsc)
call DGEMM_('N','N',nOsc,1,nOsc,1.0d0,alpha1,nOsc,r01,nOsc,0.0d0,Work(ipr_temp1),nOsc)
call DGEMM_('N','N',nOsc,1,nOsc,1.0d0,alpha2,nOsc,r02,nOsc,0.0d0,Work(ipr_temp2),nOsc)
!r_temp(:,1) = r_temp1+r_temp2
do i=1,nOsc
  Work(ipr_temp+i-1) = Work(ipr_temp1+i-1)+Work(ipr_temp2+i-1)
end do
!temp = 2.0*alpha
call dcopy_(nOscSqr,[0.0d0],0,Work(iptemp),1)
call Daxpy_(nOscSqr,2.0d0,Work(ipalpha),1,Work(iptemp),1)

call Dool_MULA(Work(iptemp),nOsc,nOsc,Work(ipr_temp),nOsc,1,det)
call GetMem('beta','Allo','Real',ipbeta,nOscSqr)
!r00 = r_temp(:,1)
do i=1,nOsc
  r00(i) = Work(ipr_temp+i-1)
end do

! Calculate beta.
call GetMem('temp1','Allo','Real',iptemp1,nOscSqr)
!temp1 = alpha1
call dcopy_(nOscSqr,alpha1,1,Work(iptemp1),1)

!temp = 2.0d0*alpha
call dcopy_(nOscSqr,[0.0d0],0,Work(iptemp),1)
call Daxpy_(nOscSqr,2.0d0,Work(ipalpha),1,Work(iptemp),1)

call Dool_MULA(Work(iptemp),nOsc,nOsc,Work(iptemp1),nOsc,nOsc,det)
call DGEMM_('N','N',nOsc,nOsc,nOsc,1.0d0,alpha2,nOsc,Work(iptemp1),nOsc,0.0d0,Work(ipbeta),nOsc)

! Calculate FC00.
!r_temp1 = r01-r02
call dcopy_(nOsc,r01,1,Work(ipr_temp1),1)
call Daxpy_(nOsc,-1.0d0,r02,1,Work(ipr_temp1),1)

call DGEMM_('N','N',nOsc,1,nOsc,1.0d0,Work(ipbeta),nOsc,Work(ipr_temp1),nOsc,0.0d0,Work(ipr_temp2),nOsc)
FC00_exp = Ddot_(nOsc,Work(ipr_temp1),1,Work(ipr_temp2),1)
FC00 = (sqrt(det1)*sqrt(det2)/det0)*exp(-FC00_exp)

call GetMem('r_temp1','Free','Real',ipr_temp1,nOsc)
call GetMem('r_temp2','Free','Real',ipr_temp2,nOsc)
call GetMem('r_temp','Free','Real',ipr_temp,nOsc)
call GetMem('alpha','Free','Real',ipalpha,nOscSqr)
call GetMem('beta','Free','Real',ipbeta,nOscSqr)
call GetMem('temp','Free','Real',iptemp,nOscSqr)
call GetMem('temp1','Free','Real',iptemp1,nOscSqr)

! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(W1)
  call Unused_real_array(W2)
end if

end subroutine Calc_r00

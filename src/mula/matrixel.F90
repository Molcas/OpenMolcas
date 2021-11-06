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

subroutine MatrixElements(L,U,FC00,Hmat,C,W,r_diff,mMat,nMat,iCre,iann,max_nOrd,max_mOrd,nOsc,energy,grad,Hess,D3,D4,G,Gprime, &
                          Gdbleprime,alpha1,alpha2,beta,max_term,Base)
!  Purpose:
!    Set up Hamilton matrix at a given center.
!
!  Input:
!
!  Output:
!
!  Uses:
!    Linalg
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

!use Linalg
!use Potkin
implicit real*8(a-h,o-z)
#include "dims.fh"
real*8 L(0:max_mOrd,0:max_mOrd)
real*8 U(0:max_nOrd,0:max_nOrd)
real*8 Hmat(0:max_mOrd,0:max_nOrd)
real*8 C(nosc,nosc)
real*8 W(nosc,nosc)
real*8 r_diff(nosc)
integer mMat(0:mdim1,mdim2)
integer nMat(0:ndim1,ndim2)
integer icre(0:ndim1,ndim2)
integer iann(0:ndim1,ndim2)
real*8 grad(nosc)
real*8 Hess(nosc,nosc)
real*8 D3(nosc,nosc,nosc)
real*8 D4(nosc,nosc,nosc,nosc)
real*8 G(nosc,nosc)
real*8 Gprime(nosc,nosc,nosc)
real*8 Gdbleprime(nosc,nosc,nosc,nosc)
real*8 alpha1(nosc,nosc), alpha2(nosc,nosc), beta(nosc,nosc)
real*8 Base(nosc,nosc)
#include "WrkSpc.fh"

! Initialize.
noscOld = nOsc
mPlus = max_mOrd+1
nPlus = max_nOrd+1
n_A = (max_mOrd+1)*(max_nOrd+1)
call GetMem('A','Allo','Real',ipA,n_A)
call dcopy_(n_A,[0.0d0],0,Work(ipA),1)
!A = 0.0d0

call GetMem('Wtemp','Allo','Real',ipWtemp,nOscold*nOsc)
call GetMem('Ctemp','Allo','Real',ipCtemp,nOsc*nOsc)
call GetMem('temp','Allo','Real',iptemp,nOsc*nOsc)

call DGEMM_('N','N',nOscold,nOsc,nOsc,1.0d0,Base,nOscOld,W,nOsc,0.0d0,Work(ipWtemp),nOscold)
call dcopy_(nOsc**2,[0.0d0],0,Work(ipCtemp),1)
call dcopy_(nOsc,[1.0d0],0,Work(ipCtemp),nOsc+1)
!temp = W
call dcopy_(nOsc*nOsc,W,1,Work(iptemp),1)
call Dool_MULA(Work(iptemp),nOsc,nOsc,Work(ipCtemp),nOsc,nOsc,det)
call GetMem('temp','Free','Real',iptemp,nOsc*nOsc)

call GetMem('rtemp1','Allo','Real',iprtemp1,nOsc)
call DGEMM_('N','N',nOsc,1,nOsc,1.0d0,Work(ipCtemp),nOsc,r_diff,nOsc,0.0d0,Work(iprtemp1),nOsc)

call PotEnergy(Work(ipA),nMat,iCre,iAnn,energy,grad,Hess,D3,D4,max_term,Work(ipWtemp),ndim1,ndim2,nOscOld)
!l_nMat_1 = ndim1
!l_nMat_2 = ndim2
call KinEnergy_drv(Work(ipA),nMat,iCre,iAnn,G,Gprime,Gdbleprime,max_term,C,W,alpha1,alpha2,beta,Work(iprtemp1),ndim1,ndim2,nOscOld)

! Calculate Hamilton matrix.
call GetMem('temp','Allo','Real',iptemp,mPlus*nPlus)
call DGEMM_('N','T',mPlus,nPlus,nPlus,1.0d0,Work(ipA),mPlus,U,nPlus,0.0d0,Work(ipTemp),mPlus)
call DGEMM_('N','N',mPlus,nPlus,mPlus,1.0d0,L,mPlus,Work(ipTemp),mPlus,0.0d0,Hmat,mPlus)
call dscal_((max_mOrd+1)*(max_nOrd+1),FC00,Hmat,1)
! Hmat = FC00*Hmat

call GetMem('A','Free','Real',ipA,n_A)
call GetMem('temp','Free','Real',iptemp,mPlus*nPlus)
call GetMem('Wtemp','Free','Real',ipWtemp,nOscold*nOsc)
call GetMem('Ctemp','Free','Real',ipCtemp,nOsc*nOsc)
call GetMem('rtemp1','Free','Real',iprtemp1,nOsc)

! Avoid unused argument warnings
if (.false.) call Unused_integer_array(mMat)

end subroutine MatrixElements

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

use mula_global, only: mdim1, mdim2, ndim1, ndim2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One

!use Linalg
!use Potkin
implicit real*8(a-h,o-z)
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
real*8, allocatable :: A(:,:), Ctemp(:,:), rtemp1(:), temp(:,:), Wtemp(:,:)

! Initialize.
noscOld = nOsc
mPlus = max_mOrd+1
nPlus = max_nOrd+1
call mma_allocate(A,[0,max_mOrd],[0,max_nOrd],label='A')
A(:,:) = Zero

call mma_allocate(Wtemp,nOscOld,nOsc,label='Wtemp')
call mma_allocate(Ctemp,nOsc,nOsc,label='Ctemp')
call mma_allocate(temp,nOsc,nOsc,label='temp')

call DGEMM_('N','N',nOscold,nOsc,nOsc,One,Base,nOscOld,W,nOsc,Zero,Wtemp,nOscold)
Ctemp(:,:) = Zero
do i=1,nOsc
  Ctemp(i,i) = One
end do
temp(:,:) = W
call Dool_MULA(temp,nOsc,nOsc,Ctemp,nOsc,nOsc,det)
call mma_deallocate(temp)

call mma_allocate(rtemp1,nOsc,label='rtemp1')
call DGEMM_('N','N',nOsc,1,nOsc,One,Ctemp,nOsc,r_diff,nOsc,Zero,rtemp1,nOsc)

call PotEnergy(A,nMat,iCre,iAnn,energy,grad,Hess,D3,D4,max_term,Wtemp,ndim1,ndim2,nOscOld)
!l_nMat_1 = ndim1
!l_nMat_2 = ndim2
call KinEnergy(A,nMat,iCre,iAnn,G,Gprime,Gdbleprime,max_term,C,W,alpha1,alpha2,beta,rtemp1,ndim1,ndim2,nOscOld)

! Calculate Hamilton matrix.
call mma_allocate(temp,[0,max_mOrd],[0,max_nOrd],label='temp')
call DGEMM_('N','T',mPlus,nPlus,nPlus,One,A,mPlus,U,nPlus,Zero,Temp,mPlus)
call DGEMM_('N','N',mPlus,nPlus,mPlus,FC00,L,mPlus,Temp,mPlus,Zero,Hmat,mPlus)

call mma_deallocate(A)
call mma_deallocate(Ctemp)
call mma_deallocate(Wtemp)
call mma_deallocate(temp)
call mma_deallocate(rtemp1)

! Avoid unused argument warnings
if (.false.) call Unused_integer_array(mMat)

end subroutine MatrixElements

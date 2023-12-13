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

subroutine SetUpHmat2(energy1,C,W,det,r1,max_mOrd,max_nOrd,max_mInc,max_nInc,mMat,nMat,mInc,nInc,mDec,nDec,H,S,Hess,G0,Base, &
                      nDimTot,nOsc)
!  Purpose:
!
!  Input:
!
!  Output:
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

use mula_global, only: mdim1, mdim2, ndim1, ndim2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: max_mOrd, max_nOrd, max_mInc, max_nInc, mMat(0:mdim1,mdim2), nMat(0:ndim1,ndim2), &
                                 mInc(0:mdim1,mdim2), nInc(0:ndim1,ndim2), mDec(0:mdim1,mdim2), nDec(0:ndim1,ndim2), nDimTot, nOsc
real(kind=wp), intent(in) :: energy1, C(nOsc,nOsc), W(nOsc,nOsc), det, r1(nOsc), Hess(nOsc,nOsc), G0(nOsc,nOsc), Base(nOsc,nOsc)
real(kind=wp), intent(out) :: H(nDimTot,nDimTot), S(nDimTot,nDimTot)
integer(kind=iwp) :: max_term
real(kind=wp) :: det0, FC00
real(kind=wp), allocatable :: alpha1(:,:), alpha2(:,:), beta(:,:), C0(:,:), D3(:,:,:), D4(:,:,:,:), Gdbleprime(:,:,:,:), &
                              Gprime(:,:,:), grad(:), Hij(:,:), L(:,:), r0(:), r_diff(:), Sij(:,:), U(:,:), W0(:,:)

! Initialize.
! arrays for setuphmat2
call mma_allocate(Hij,[0,max_mOrd],[0,max_nOrd],label='Hij')
call mma_allocate(Sij,[0,max_mOrd],[0,max_nOrd],label='Sij')
call mma_allocate(r0,nOsc,label='r0')
call mma_allocate(r_diff,nOsc,label='r_diff')
call mma_allocate(alpha1,nOsc,nOsc,label='alpha1')
call mma_allocate(alpha2,nOsc,nOsc,label='alpha2')
call mma_allocate(beta,nOsc,nOsc,label='beta')
call mma_allocate(L,[0,max_mOrd],[0,max_mOrd],label='L')
call mma_allocate(U,[0,max_nOrd],[0,max_nOrd],label='U')
call mma_allocate(C0,nOsc,nOsc,label='C0')
call mma_allocate(W0,nOsc,nOsc,label='W0')
call mma_allocate(grad,nOsc,label='grad')
call mma_allocate(Gprime,nOsc,nOsc,nOsc,label='Gprime')
call mma_allocate(D3,nOsc,nOsc,nOsc,label='D3')
call mma_allocate(Gdbleprime,nOsc,nOsc,nOsc,nOsc,label='Gdbleprime')
call mma_allocate(D4,nOsc,nOsc,nOsc,nOsc,label='D4')

! - Call Franck-Condon routine.
! - Calculate matrix elements.

max_term = 2
grad(:) = Zero
D3(:,:,:) = Zero
D4(:,:,:,:) = Zero
Gprime(:,:,:) = Zero
Gdbleprime(:,:,:,:) = Zero
!call unitmat(Base,nOsc)
call Calc_r00(C,C,C0,W0,alpha1,alpha2,r0,r1,r1,det0,det,det,FC00,nOsc)
call FCval(C,W,det0,r0,C,W,det0,r0,Sij,max_mOrd,max_nOrd,max_nOrd,max_mInc,max_nInc,max_nInc,mMat,nMat,mInc,nInc,mDec,nDec,C0,W0, &
           det0,L,U,FC00,alpha1,alpha2,beta,nOsc)
r_diff(:) = Zero
call MatrixElements(L,U,FC00,Hij,C0,W0,r_diff,nMat,nInc,nDec,max_nOrd,max_mOrd,nOsc,energy1,grad,Hess,D3,D4,G0,Gprime,Gdbleprime, &
                    alpha1,alpha2,beta,max_term,Base)

H(1:max_mOrd+1,1:max_nOrd+1) = Hij
S(1:max_mOrd+1,1:max_nOrd+1) = Sij

call mma_deallocate(Hij)
call mma_deallocate(Sij)
call mma_deallocate(r0)
call mma_deallocate(r_diff)
call mma_deallocate(alpha1)
call mma_deallocate(alpha2)
call mma_deallocate(beta)
call mma_deallocate(L)
call mma_deallocate(U)
call mma_deallocate(C0)
call mma_deallocate(W0)
call mma_deallocate(grad)
call mma_deallocate(D3)
call mma_deallocate(Gprime)
call mma_deallocate(D4)
call mma_deallocate(Gdbleprime)

end subroutine SetUpHmat2

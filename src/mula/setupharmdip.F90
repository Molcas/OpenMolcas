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

subroutine SetUpHarmDip(DipMat,max_term,m_max,n_max,mMat,mInc,mDec,nMat,nInc,nDec,C1,W1,det1,r01,C2,W2,det2,r02,C,W,det0,r00, &
                        TranDip,TranDipGrad,FC00,nnsiz,max_mOrd,max_nOrd,nOsc)
!  Purpose:
!    Calculate the matrix elements of the transition dipole moment
!    at the location of the intermediate oscillator.
!
!  Input:
!    max_term   : Integer - maximum order of the transition dipole terms.
!    W1,W2      : Real*8 two dimensional arrays - eigenvectors
!                 scaled by the square root of the eigenvalues.
!    C1,C2      : Real*8 two dimensional arrays - inverses
!                 of W1 and W2.
!    det1,det2  : Real*8 variables - determinants of C1 and C2.
!    r01,r02    : Real*8 arrays - coordinates of the two
!                 oscillators.
!    Forcefield : Logical variable - whether or not to use transition
!                 dipole from input.
!
!  Output:
!    DipMat     : Real*8 two dimensional array - contains the
!                 matrix elements of the transition dipole.
!
!  Uses:
!    TabMod
!    FCMod

use mula_global, only: mdim1, mdim2, ndim1, ndim2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One

!use TabMod
!use FCMod
implicit real*8(a-h,o-z)
real*8 DipMat(0:max_mOrd,0:max_nOrd,0:3)
integer mMat(0:mdim1,mdim2), mInc(0:mdim1,mdim2), mDec(0:mdim1,mdim2)
integer nMat(0:ndim1,ndim2), nInc(0:ndim1,ndim2), nDec(0:ndim1,ndim2)
real*8 C1(nosc,nosc), C2(nosc,nosc), W1(nosc,nosc), W2(nosc,nosc), C(nosc,nosc), W(nosc,nosc)
real*8 r01(nosc), r02(nosc), r00(nosc)
real*8 TranDipGrad(3,nosc)
real*8 TranDip(3)
integer nvTabDim
real*8, allocatable :: alpha1(:,:), alpha2(:,:), beta(:,:), F(:,:,:), L(:,:), Sij(:,:), Temp1(:,:), Temp2(:,:), U(:,:)

! Initialize.
call TabDim(m_max,nosc,nvTabDim)
m_max_ord = nvTabDim-1
call TabDim(min(n_max,m_max+1),nosc,nvTabDim)

mx_max_ord = nvTabDim-1
call TabDim(min(m_max,n_max+1),nosc,nvTabDim)
nx_max_ord = nvTabDim-1
call TabDim(m_max-1,nosc,nvTabDim)
max_mInc = nvTabDim-1
call TabDim(n_max,nosc,nvTabDim)
n_max_ord = nvTabDim-1
call TabDim(n_max-1,nosc,nvTabDim)
max_nInc = nvTabDim-1
call mma_allocate(L,[0,m_max_ord],[0,nx_max_ord],label='L')
call mma_allocate(U,[0,n_max_ord],[0,mx_max_ord],label='U')
call mma_allocate(Sij,[0,m_max_ord],[0,n_max_ord],label='Sij')

call mma_allocate(alpha1,nOsc,nOsc,label='alpha1')
call mma_allocate(alpha2,nOsc,nOsc,label='alpha2')
call mma_allocate(beta,nOsc,nOsc,label='beta')

! Calculate Franck-Condon factors.
call FCval(C1,W1,det1,r01,C2,W2,det2,r02,Sij,m_max_ord,n_max_ord,mx_max_ord,max_mInc,max_nInc,nx_max_ord,mMat,nMat,mInc,nInc,mDec, &
           nDec,C,W,det0,r00,L,U,FC00,alpha1,alpha2,beta,nOsc,nnsiz)

call mma_deallocate(alpha1)
call mma_deallocate(alpha2)
call mma_deallocate(beta)

! Get the zeroth order contribution from transition dipole.
Dipmat(:,:,:) = Zero
do iCar=1,3
  DipMat(:,:,iCar) = TranDip(iCar)*Sij
end do
call mma_deallocate(Sij)

!rt = sqrt(TranDip(1)**2+TranDip(2)**2+TranDip(3)**2)
!if (rt > Zero) DipMat(:,:,0) = rt*Sij
!end if

! Calculate LFU (just valid for tdm first derivatives)
if (max_term == 1) then
  call mma_allocate(Temp1,[0,m_max_ord],[0,mx_max_ord],label='Temp1')
  call mma_allocate(Temp2,[0,m_max_ord],[0,n_max_ord],label='Temp2')
  if (n_max > m_max) then
    call mma_allocate(F,[0,m_max_ord],[0,mx_max_ord],[1,3],label='F')
    call Fgenerator(nmat,F,nInc,nDec,TranDipGrad,m_max_ord,mx_max_ord,nosc)
    do iCar=1,3
      call DGEMM_('N','N',m_max_ord+1,mx_max_ord+1,m_max_ord+1,One,L,m_max_ord+1,F(:,:,iCar),m_max_ord+1,Zero,Temp1,m_max_ord+1)
      call DGEMM_('N','T',m_max_ord+1,n_max_ord+1,mx_max_ord+1,One,Temp1,m_max_ord+1,U,n_max_ord+1,Zero,Temp2,m_max_ord+1)

      DipMat(:,:,iCar) = DipMat(:,:,iCar)+Temp2*FC00
    end do
  else
    call mma_allocate(F,[0,n_max_ord],[0,nx_max_ord],[1,3],label='F')
    call Fgenerator(mmat,F,mInc,mDec,trandipgrad,n_max_ord,nx_max_ord,nosc)
    do iCar=1,3
      call DGEMM_('N','T',m_max_ord+1,n_max_ord+1,nx_max_ord+1,One,L,m_max_ord+1,F(:,:,iCar),n_max_ord+1,Zero,Temp1,m_max_ord+1)
      call DGEMM_('N','T',m_max_ord+1,n_max_ord+1,n_max_ord+1,One,Temp1,m_max_ord+1,U,n_max_ord+1,Zero,Temp2,m_max_ord+1)
      DipMat(:,:,iCar) = DipMat(:,:,iCar)+Temp2*FC00
    end do
  end if
  call mma_deallocate(Temp1)
  call mma_deallocate(Temp2)
  call mma_deallocate(F)

end if

call mma_deallocate(L)
call mma_deallocate(U)

end subroutine SetUpHarmDip

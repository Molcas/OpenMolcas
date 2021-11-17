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

subroutine SetUpDipMat2(DipMat,max_term,C1,W1,det1,r01,C2,W2,det2,r02,max_mOrd,max_nOrd,max_nOrd2,max_mInc,max_nInc,max_nInc2, &
                        mMat,nMat,mInc,nInc,mDec,nDec,det0,base,TranDip,TranDipGrad,nOsc,nDimTot)
!  Purpose:
!    Performs a least squares fit of the transition dipole at the two
!    centers and at the intermediate oscillator. Calculates the matrix
!    elements of the transition dipole at these centers.
!
!  Input:
!    ipow       : Two dimensional integer array - terms of the polynomial.
!    var        : Real two dimensional array - coordinates to be used in the fit.
!    dip        : Real array - values of dipole at the coordinates contained in var.
!    trfName    : Character array - transformation associated with each internal coordinate.
!    use_weight : Logical
!    max_term   : Integer - maximum order of the transition dipole terms.
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
!    DipMat     : Real two dimensional array - contains the matrix elements of the transition dipole.

use mula_global, only: mdim1, mdim2, ndim1, ndim2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: max_term, max_mOrd, max_nOrd, max_nOrd2, max_mInc, max_nInc, max_nInc2, mMat(0:mdim1,mdim2), &
                                 nMat(0:ndim1,ndim2), mInc(0:mdim1,mdim2), nInc(0:ndim1,ndim2), mDec(0:mdim1,mdim2), &
                                 nDec(0:ndim1,ndim2), nOsc, nDimTot
real(kind=wp), intent(out) :: DipMat(0:nDimTot,0:nDimTot), det0
real(kind=wp), intent(in) :: C1(nOsc,nOsc), W1(nOsc,nOsc), det1, r01(nOsc), C2(nOsc,nOsc), W2(nOsc,nOsc), det2, r02(nOsc), &
                             Base(nOsc,nOsc), TranDip(3), TranDipGrad(nOsc)
integer(kind=iwp) :: l_C1
real(kind=wp) :: D0(3), FC00
real(kind=wp), allocatable :: alpha1(:,:), alpha2(:,:), beta(:,:), C(:,:), D1(:), D2(:,:), D3(:,:,:), D4(:,:,:,:), Dij(:,:), &
                              L(:,:), r0vec(:), Sij(:,:), U(:,:), W(:,:)

! Initialize.
call mma_allocate(Dij,[0,max_mOrd],[0,max_mOrd],label='Dij')

call mma_allocate(C,nOsc,nOsc,label='C')
call mma_allocate(W,nOsc,nOsc,label='W')

call mma_allocate(L,[0,max_mOrd],[0,max_mOrd],label='L')
call mma_allocate(U,[0,max_nOrd2],[0,max_nOrd2],label='U')
call mma_allocate(Sij,[0,max_mOrd],[0,max_nOrd],label='Sij')
call mma_allocate(r0vec,nOsc,label='r0vec')
call mma_allocate(alpha1,nOsc,nOsc,label='alpha1')
call mma_allocate(alpha2,nOsc,nOsc,label='alpha2')
call mma_allocate(beta,nOsc,nOsc,label='beta')
call mma_allocate(D1,nOsc,label='D1')
call mma_allocate(D2,nOsc,nOsc,label='D2')
call mma_allocate(D3,nOsc,nOsc,nOsc,label='D3')
call mma_allocate(D4,nOsc,nOsc,nOsc,nOsc,label='D4')

! Calculate terms of type
!
!                    dM  |
!                    --  | < i | Q  | j >,
!                    dQ  |        k
!                      k  Q
!                          0
! where M is the (transition) dipole moment, |i> and |j> are harmonic
! oscillator states, Q_0 is the equilibrium geometry and Q_k is the
! k:th normal coordinate.

l_C1 = nOsc
call Calc_r00(C1,C2,C,W,alpha1,alpha2,r0vec,r01,r02,det0,det1,det2,FC00,l_C1)
call FCval(C1,W1,det1,r01,C2,W2,det2,r02,Sij,max_mOrd,max_nOrd,max_nOrd2,max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,nInc,mDec, &
           nDec,C,W,det0,L,U,FC00,alpha1,alpha2,beta,l_C1)
D0(1) = TranDip(1)
D0(2) = TranDip(2)
D0(3) = TranDip(3)
D1(:) = TranDipGrad
D2(:,:) = Zero
D3(:,:,:) = Zero
D4(:,:,:,:) = Zero

call DipMatEl(Dij,W,L,U,FC00,nMat,nInc,nDec,D0(1),D1,D2,D3,D4,max_term,Base,ndim1,ndim2,max_mOrd,max_nOrd2)
DipMat(0:max_mOrd,0:max_mOrd) = Dij

call mma_deallocate(Dij)
call mma_deallocate(C)
call mma_deallocate(W)
call mma_deallocate(L)
call mma_deallocate(U)
call mma_deallocate(r0vec)
call mma_deallocate(alpha1)
call mma_deallocate(alpha2)
call mma_deallocate(beta)
call mma_deallocate(D1)
call mma_deallocate(D2)
call mma_deallocate(D3)
call mma_deallocate(D4)

end subroutine SetUpDipMat2

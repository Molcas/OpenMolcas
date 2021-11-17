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
! Copyright (C) 1996,1999, Niclas Forsberg                             *
!               1996,1999, Anders Bernhardsson                         *
!***********************************************************************

!module MatElMod

!  Contains:
!    MatrixElements (L,U,FC00,Hmat,C,W,r_diff,mMat,nMat,
!                    max_nOrd,energy,grad,Hess,D3,D4,G,
!                    Gprime,Gdbleprime,alpha1,alpha2,beta,max_term)
!    LSPotFit       (r01,energy1,grad1,Hess1,D3_1,D4_1,
!                    r02,energy2,grad2,Hess2,D3_2,D4_2,
!                    r00,energy0,r_min,FitCoef,mMat,stand_dev,max_err,
!                    use_weight,max_term,pot)
!    SetUpHmat      (energy0,r_min,ipow,var,yin,coef,r00,trfName,max_term,
!                    C1,W1,det1,r01,C2,W2,det2,r02,max_mOrd,max_nOrd,
!                    max_nOrd2,max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,
!                    nInc,mDec,nDec,L,U,H,S,G1,G2,G0,Gprime1,Gprime2,
!                    Gprime0,Gdbleprime1,Gdbleprime2,Gdbleprime0,
!                    C0,W0,det0,Mass,rOrigin)
!
!  Written by:
!    Niclas Forsberg & Anders Bernhardsson,
!    Dept. of Theoretical Chemistry, Lund University, 1996.
!    Dept. of Theoretical Chemistry, Lund University, 1999.

!vv private

!contains

subroutine SetUpHmat(energy0,r_min,ipow,var,yin,trfName,max_term,C1,W1,det1,r01,C2,W2,det2,r02,max_mOrd,max_nOrd,max_nOrd2, &
                     max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,nInc,mDec,nDec,H,S,G1,G2,G0,Gprime1,Gprime2,Gprime0,Gdbleprime1, &
                     Gdbleprime2,Gdbleprime0,det0,Base,r0,r1,r2,nterm,nvar,ndata,nOsc,nDimTot)
!  Purpose:
!    Set up Hamilton matrix.
!
!  Input:
!    Energy
!    r_min
!    ipow
!    var
!    yin
!    coeff
!    trfname
!    Max_term
!    C1,W1,det1,r01
!    C2,W2,det2,r02
!    H S
!    G1 G2 G0 1' G2' g0' g0'' g1'' g2''
!    det0
!
!  The expansion point for lspotfit is r00!!!!!

use mula_global, only: mdim1, mdim2, ndim1, ndim2, ngdim, nPolyTerm
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nterm, nvar, ipow(nterm,nvar), max_term, max_mOrd, max_nOrd, max_nOrd2, max_mInc, max_nInc, &
                                 max_nInc2, mMat(0:mdim1,mdim2), nMat(0:ndim1,ndim2), mInc(0:mdim1,mdim2), nInc(0:ndim1,ndim2), &
                                 mDec(0:mdim1,mdim2), nDec(0:ndim1,ndim2), ndata, nOsc, nDimTot
real(kind=wp), intent(in) :: energy0, r_min(nOsc), var(ndata,nvar), yin(ndata), C1(nOsc,nOsc), W1(nOsc,nOsc), det1, r01(nOsc), &
                             C2(nOsc,nOsc), W2(nOsc,nOsc), det2, r02(nOsc), G1(nOsc,nOsc), G2(nOsc,nOsc), G0(nOsc,nOsc), &
                             Gprime1(ngdim,ngdim,ngdim), Gprime2(ngdim,ngdim,ngdim), Gprime0(ngdim,ngdim,ngdim), &
                             Gdbleprime1(ngdim,ngdim,ngdim,ngdim), Gdbleprime2(ngdim,ngdim,ngdim,ngdim), &
                             Gdbleprime0(ngdim,ngdim,ngdim,ngdim), Base(nOsc,nOsc), r0(nOsc)
real(kind=wp), intent(out) :: H(nDimTot,nDimTot), S(nDimTot,nDimTot), det0, r1(nOsc), r2(nOsc)
character(len=80), intent(in) :: trfName(nvar)
integer(kind=iwp) :: i, j, l_C1, l_C2, l_r2, nOscOld, numCoef
real(kind=wp) :: energy, energy1, energy2, FC00, max_err, stand_dev
logical(kind=iwp) :: find_minimum, pot, use_weight
integer(kind=iwp), allocatable :: jPow(:,:)
real(kind=wp), allocatable :: alpha1(:,:), alpha2(:,:), beta(:,:), C(:,:), Coef(:), D3(:,:,:), D3_1(:,:,:), D3_2(:,:,:), &
                              D4(:,:,:,:), D4_1(:,:,:,:), D4_2(:,:,:,:), FitCoef(:), grad(:), grad1(:), grad2(:), &
                              Gdbleprimetemp(:,:,:,:), Gprimetemp(:,:,:), Gtemp(:,:), Hess(:,:), Hess1(:,:), Hess2(:,:), Hij(:,:), &
                              HijTrans(:,:), L(:,:), r_diff(:), r0vec(:), Sij(:,:), SijTrans(:,:), U(:,:), W(:,:)

! Initialize.
l_r2 = nOsc
find_minimum = .false.
use_weight = .false.
nOscOld = nOsc

call mma_allocate(Coef,nPolyTerm,label='Coef')
call mma_allocate(grad1,nOscOld,label='grad1')
call mma_allocate(grad2,nOscOld,label='grad2')
call mma_allocate(Hess1,nOscOld,nOscOld,label='Hess1')
call mma_allocate(Hess2,nOscOld,nOscOld,label='Hess2')
call mma_allocate(D3_1,nOscOld,nOscOld,nOscOld,label='D3_1')
call mma_allocate(D3_1,nOscOld,nOscOld,nOscOld,label='D3_1')
call mma_allocate(D4_2,nOscOld,nOscOld,nOscOld,nOscOld,label='D4_2')
call mma_allocate(D4_2,nOscOld,nOscOld,nOscOld,nOscOld,label='D4_2')

! Fit polynomial and calculate energy, gradient, Hessian and third
! and fourth order force constants around r01 and r02.
call PotFit(nterm,nvar,ndata,ipow,var,yin,Coef,r1,nOscOld,energy1,grad1,Hess1,D3_1,D4_1,trfName,stand_dev,max_err,find_minimum, &
            max_term,use_weight,nOscOld,nOscOld,nOscOld)
call PotFit(nterm,nvar,ndata,ipow,var,yin,Coef,r2,l_r2,energy2,grad2,Hess2,D3_2,D4_2,trfName,stand_dev,max_err,find_minimum, &
            max_term,use_weight,nOscOld,nOscOld,nOscOld)
call mma_deallocate(Coef)

! Perform a least squares fit.
call TabDim(max_term,nOsc,numCoef)
call mma_allocate(FitCoef,numCoef,label='FitCoef')

call mma_allocate(jPow,numCoef,nOsc,label='jPow')
pot = .true.
call LSPotFit(r1,energy1,grad1,Hess1,D3_1,D4_1,r2,energy2,grad2,Hess2,D3_2,D4_2,r0,energy0,r_min,FitCoef,jPow,max_term,pot,nOsc, &
              numcoef)

call mma_deallocate(grad1)
call mma_deallocate(grad2)
call mma_deallocate(Hess1)
call mma_deallocate(Hess2)
call mma_deallocate(D3_1)
call mma_deallocate(D3_2)
call mma_deallocate(D4_1)
call mma_deallocate(D4_2)

! For each of the centers:
! - Call Franck-Condon routine.
! - Calculate matrix elements.

call mma_allocate(Hij,[0,max_mOrd],[0,max_nOrd],label='Hij')
call mma_allocate(Sij,[0,max_mOrd],[0,max_nOrd],label='Sij')
call mma_allocate(r0vec,nOscOld,label='r0vec')
call mma_allocate(r_diff,nOscOld,label='r_diff')
call mma_allocate(C,nOsc,nOsc,label='C')
call mma_allocate(W,nOsc,nOsc,label='W')
call mma_allocate(grad,nOscOld,label='grad')
call mma_allocate(D3,nOscOld,nOscOld,nOscOld,label='D3')
call mma_allocate(D4,nOscOld,nOscOld,nOscOld,nOscOld,label='D4')
call mma_allocate(Gtemp,nOsc,nOsc,label='Gtemp')
call mma_allocate(Gprimetemp,nOsc,nOsc,nOsc,label='Gprimetemp')
call mma_allocate(Gdbleprimetemp,nOsc,nOsc,nOsc,nOsc,label='Gdbleprimetemp')
call mma_allocate(alpha1,nOsc,nOsc,label='alpha1')
call mma_allocate(alpha2,nOsc,nOsc,label='alpha2')
call mma_allocate(beta,nOsc,nOsc,label='beta')
call mma_allocate(L,[0,max_mOrd],[0,max_mOrd],label='L')
call mma_allocate(U,[0,max_nOrd],[0,max_nOrd],label='U')

! Block 11.
l_C1 = nOsc
call Calc_r00(C1,C1,C,W,alpha1,alpha2,r0vec,r01,r01,det0,det1,det1,FC00,l_C1)
call FCval(C1,W1,det1,r01,C1,W1,det1,r01,Sij,max_mOrd,max_nOrd,max_nOrd2,max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,nInc,mDec, &
           nDec,C,W,det1,L,U,FC00,alpha1,alpha2,beta,l_C1)
r0vec(1:nOsc) = r1-r0
call funcval(r0vec,FitCoef,jPow,energy,nterm,nvar)
call gradient(r0vec,FitCoef,jPow,grad,nterm,nvar)
call Hessian(r0vec,FitCoef,jPow,Hess,nterm,nvar)
call thirdDer(r0vec,FitCoef,jPow,D3,nterm,nvar)
call fourthDer(r0vec,FitCoef,jPow,D4,nterm,nvar)
energy = energy-energy0
r0vec(:) = Zero
r_diff(:) = Zero
Gtemp(:,:) = G1
GprimeTemp(:,:,:) = Gprime1
GdbleprimeTemp(:,:,:,:) = Gdbleprime1
call MatrixElements(L,U,FC00,Hij,C,W,r_diff,nMat,ninc,ndec,max_nOrd,max_mOrd,nOsc,energy,grad,Hess,D3,D4,Gtemp,GprimeTemp, &
                    GdbleprimeTemp,alpha1,alpha2,beta,max_term,Base)

H(1:max_mOrd+1,1:max_nOrd+1) = Hij
S(1:max_mOrd+1,1:max_nOrd+1) = Sij

! Block 22.
l_C2 = nOsc
call Calc_r00(C2,C2,C,W,alpha1,alpha2,r0vec,r02,r02,det0,det2,det2,FC00,l_C2)
call FCval(C2,W2,det2,r02,C2,W2,det2,r02,Sij,max_mOrd,max_nOrd,max_nOrd2,max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,nInc,mDec, &
           nDec,C,W,det2,L,U,FC00,alpha1,alpha2,beta,l_C2)
r0vec(1:nOsc) = r2-r0
call funcval(r0vec,FitCoef,jPow,energy,nterm,nvar)
call gradient(r0vec,FitCoef,jPow,grad,nterm,nvar)
call Hessian(r0vec,FitCoef,jPow,Hess,nterm,nvar)
call thirdDer(r0vec,FitCoef,jPow,D3,nterm,nvar)
call fourthDer(r0vec,FitCoef,jPow,D4,nterm,nvar)
energy = energy-energy0
r0vec(:) = Zero
r_diff(:) = Zero
Gtemp(:,:) = G2
GprimeTemp(:,:,:) = Gprime2
GdbleprimeTemp(:,:,:,:) = Gdbleprime2
call MatrixElements(L,U,FC00,Hij,C,W,r_diff,nMat,ninc,ndec,max_nOrd,max_mOrd,nOsc,energy,grad,Hess,D3,D4,Gtemp,GprimeTemp, &
                    GdbleprimeTemp,alpha1,alpha2,beta,max_term,Base)
H(max_mOrd+2:2*max_mOrd+2,max_nOrd+2:2*max_nOrd+2) = Hij
S(max_mOrd+2:2*max_mOrd+2,max_nOrd+2:2*max_nOrd+2) = Sij

! Block 12 and 21.
l_C1 = nOsc
call Calc_r00(C1,C2,C,W,alpha1,alpha2,r0vec,r01,r02,det0,det1,det2,FC00,l_C1)
call FCval(C1,W1,det1,r01,C2,W2,det2,r02,Sij,max_mOrd,max_nOrd,max_nOrd2,max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,nInc,mDec, &
           nDec,C,W,det0,L,U,FC00,alpha1,alpha2,beta,l_C1)
r0vec(:) = Zero
call funcval(r0vec,FitCoef,jPow,energy,nterm,nvar)
call gradient(r0vec,FitCoef,jPow,grad,nterm,nvar)
call Hessian(r0vec,FitCoef,jPow,Hess,nterm,nvar)
call thirdDer(r0vec,FitCoef,jPow,D3,nterm,nvar)
call fourthDer(r0vec,FitCoef,jPow,D4,nterm,nvar)
energy = energy-energy0
r0vec(:) = Zero
r_diff(1:nOsc) = r01-r02
Gtemp(:,:) = G0
GprimeTemp(:,:,:) = Gprime0
GdbleprimeTemp(:,:,:,:) = Gdbleprime0
call MatrixElements(L,U,FC00,Hij,C,W,r_diff,nMat,nInc,nDec,max_nOrd,max_mOrd,nOsc,energy,grad,Hess,D3,D4,Gtemp,GprimeTemp, &
                    GdbleprimeTemp,alpha1,alpha2,beta,max_term,Base)
H(1:max_mOrd+1,max_nOrd+2:2*max_nOrd+2) = Hij
S(1:max_mOrd+1,max_nOrd+2:2*max_nOrd+2) = Sij

call mma_allocate(HijTrans,[0,max_nOrd],[0,max_mOrd],label='HijTrans')
call mma_allocate(SijTrans,[0,max_nOrd],[0,max_mOrd],label='SijTrans')
do i=0,max_mOrd
  do j=0,max_nOrd
    HijTrans(j,i) = Hij(i,j)
    SijTrans(j,i) = Sij(i,j)
  end do
end do
H(max_mOrd+2:2*max_mOrd+2,1:max_nOrd+1) = HijTrans
S(max_mOrd+2:2*max_mOrd+2,1:max_nOrd+1) = SijTrans

call mma_deallocate(jPow)
call mma_deallocate(FitCoef)
call mma_deallocate(Hij)
call mma_deallocate(Sij)
call mma_deallocate(HijTrans)
call mma_deallocate(SijTrans)
call mma_deallocate(r0vec)
call mma_deallocate(r_diff)
call mma_deallocate(C)
call mma_deallocate(W)
call mma_deallocate(grad)
call mma_deallocate(Hess)
call mma_deallocate(D3)
call mma_deallocate(D4)
call mma_deallocate(Gtemp)
call mma_deallocate(Gprimetemp)
call mma_deallocate(Gdbleprimetemp)
call mma_deallocate(alpha1)
call mma_deallocate(alpha2)
call mma_deallocate(beta)
call mma_deallocate(L)
call mma_deallocate(U)

end subroutine SetUpHmat

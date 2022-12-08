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
!    MatrixElements (L,U,FC00,Hmat,C,W,r_diff,mMat,nMat,max_nOrd,energy,grad,Hess,D3,D4,G,Gprime,Gdbleprime,alpha1,alpha2,beta,
!                    max_term)
!    LSPotFit       (r01,energy1,grad1,Hess1,D3_1,D4_1,r02,energy2,grad2,Hess2,D3_2,D4_2,r00,energy0,r_min,FitCoef,mMat,stand_dev,
!                    max_err,use_weight,max_term,pot)
!    SetUpHmat      (energy0,r_min,ipow,var,yin,coef,r00,trfName,max_term,C1,W1,det1,r01,C2,W2,det2,r02,max_mOrd,max_nOrd,max_nOrd2,
!                    max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,nInc,mDec,nDec,L,U,H,S,G1,G2,G0,Gprime1,Gprime2,Gprime0,Gdbleprime1,
!                    Gdbleprime2,Gdbleprime0,C0,W0,det0,Mass,rOrigin)
!
!  Written by:
!    Niclas Forsberg & Anders Bernhardsson,
!    Dept. of Theoretical Chemistry, Lund University, 1996.
!    Dept. of Theoretical Chemistry, Lund University, 1999.

!contains

subroutine MatrixElements(L,U,FC00,Hmat,C,W,r_diff,nMat,iCre,iann,max_nOrd,max_mOrd,nOsc,energy,grad,Hess,D3,D4,G,Gprime, &
                          Gdbleprime,alpha1,alpha2,beta,max_term,Base)
!  Purpose:
!    Set up Hamilton matrix at a given center.
!
!  Input:
!
!  Output:
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

use mula_global, only: ndim1, ndim2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nMat(0:ndim1,ndim2), icre(0:ndim1,ndim2), iann(0:ndim1,ndim2), max_nOrd, max_mOrd, nOsc, max_term
real(kind=wp), intent(out) :: Hmat(0:max_mOrd,0:max_nOrd)
real(kind=wp), intent(in) :: L(0:max_mOrd,0:max_mOrd), U(0:max_nOrd,0:max_nOrd), FC00, C(nOsc,nOsc), W(nOsc,nOsc), r_diff(nOsc), &
                             energy, grad(nOsc), Hess(nOsc,nOsc), D3(nOsc,nOsc,nOsc), D4(nOsc,nOsc,nOsc,nOsc), G(nOsc,nOsc), &
                             Gprime(nOsc,nOsc,nOsc), Gdbleprime(nOsc,nOsc,nOsc,nOsc), alpha1(nOsc,nOsc), alpha2(nOsc,nOsc), &
                             beta(nOsc,nOsc), Base(nOsc,nOsc)
integer(kind=iwp) :: mPlus, nOscOld, nPlus
real(kind=wp) :: det
real(kind=wp), allocatable :: A(:,:), Ctemp(:,:), rtemp1(:), temp(:,:), Wtemp(:,:)

! Initialize.
nOscOld = nOsc
mPlus = max_mOrd+1
nPlus = max_nOrd+1
call mma_allocate(A,[0,max_mOrd],[0,max_nOrd],label='A')
A(:,:) = Zero

call mma_allocate(Wtemp,nOscOld,nOsc,label='Wtemp')
call mma_allocate(Ctemp,nOsc,nOsc,label='Ctemp')
call mma_allocate(temp,nOsc,nOsc,label='temp')

call DGEMM_('N','N',nOscold,nOsc,nOsc,One,Base,nOscOld,W,nOsc,Zero,Wtemp,nOscold)
call unitmat(Ctemp,nOsc)
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

end subroutine MatrixElements

subroutine LSPotFit(r01,energy1,grad1,Hess1,D3_1,D4_1,r02,energy2,grad2,Hess2,D3_2,D4_2,r00,energy0,r_min,FitCoef,mMat,max_term, &
                    pot,nOsc,numcoef)
!  Purpose:
!    Perform a least squares fit of the potential at two different centers, r01 and r02.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Ten, Pi, auTocm
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nOsc, numcoef, max_term
integer(kind=iwp), intent(out) :: mMat(0:numcoef-1,nOsc)
real(kind=wp), intent(in) :: r01(nOsc), energy1, grad1(nOsc), Hess1(nOsc,nOsc), D3_1(nOsc,nOsc,nOsc), D4_1(nOsc,nOsc,nOsc,nOsc), &
                             r02(nOsc), energy2, grad2(nOsc), Hess2(nOsc,nOsc), D3_2(nOsc,nOsc,nOsc), D4_2(nOsc,nOsc,nOsc,nOsc), &
                             r00(nOsc), energy0, r_min(nOsc)
real(kind=wp), intent(out) :: FitCoef(numcoef,1)
logical(kind=iwp), intent(in) :: pot
integer(kind=iwp) :: i, ip, irow, iterm, ivar, j, jrow, jterm, jvar, k, krow, kvar, l, l_mMat, lrow, lvar, m, max_mInc, max_mOrd, &
                     mrow, mTabDim, mvar, n, nDim, nrow, nTabDim, nterm, nvar
real(kind=wp) :: det, pow, t
integer(kind=iwp), allocatable :: mDec(:,:), mInc(:,:)
real(kind=wp), allocatable :: Equmat(:,:), rhs(:), Temp(:,:), Tmat(:,:), vpow(:,:), weight(:,:), x(:)
integer(kind=iwp), parameter :: mxdeg = 6

! Initialize.
nvar = nOsc
mTabDim = numcoef-1
call TabDim(max_term,nOsc,nTabDim)
max_mOrd = nTabDim-1
call TabDim(max_term-1,nOsc,nTabDim)
max_mInc = nTabDim-1
call mma_allocate(mDec,[0,mTabDim],[1,nvar],label='mDec')
call mma_allocate(mInc,[0,mTabDim],[1,nvar],label='mInc')
l_mMat = nOsc
call MakeTab(max_term,max_mOrd,max_mInc,mMat,mInc,mDec,l_mMat)
nterm = max_mOrd+1

! Set up weight matrix and right hand side.
nDim = 2*nterm
!nDim = nVar*nterm !!!!
if (pot) nDim = nDim+1
call mma_allocate(rhs,nDim,label='rhs')

call mma_allocate(weight,nDim,nDim,label='weight')

weight(:,:) = Zero
n = 1
weight(n,n) = 1.0e4_wp
weight(n+nterm,n+nterm) = 1.0e4_wp
rhs(n) = energy1
rhs(n+nterm) = energy2
n = n+1
if (max_term > 0) then
  do i=1,nvar
    weight(n,n) = 1.0e3_wp
    weight(n+nterm,n+nterm) = 1.0e3_wp
    rhs(n) = grad1(i)
    rhs(n+nterm) = grad2(i)
    n = n+1
  end do
  if (max_term > 1) then
    do i=1,nvar
      do j=i,nvar
        weight(n,n) = 1.0e2_wp
        weight(n+nterm,n+nterm) = 1.0e2_wp
        rhs(n) = Hess1(i,j)
        rhs(n+nterm) = Hess2(i,j)
        n = n+1
      end do
    end do
    if (max_term > 2) then
      do i=1,nvar
        do j=i,nvar
          do k=j,nvar
            weight(n,n) = Ten
            weight(n+nterm,n+nterm) = Ten
            rhs(n) = D3_1(i,j,k)
            rhs(n+nterm) = D3_2(i,j,k)
            n = n+1
          end do
        end do
      end do
      if (max_term > 3) then
        do i=1,nvar
          do j=i,nvar
            do k=j,nvar
              do l=k,nvar
                weight(n,n) = One
                weight(n+nterm,n+nterm) = One
                rhs(n) = D4_1(i,j,k,l)
                rhs(n+nterm) = D4_2(i,j,k,l)
                n = n+1
              end do
            end do
          end do
        end do
      end if
    end if
  end if
end if
if (pot) then
  rhs(nDim) = energy0+2.0e4_wp/auTocm
  weight(nDim,nDim) = 1.0e4_wp
end if

call mma_allocate(vpow,[0,mxdeg],[1,nvar],label='vpow')
call mma_allocate(Tmat,nDim,nterm,label='Tmat')
Tmat(:,:) = Zero
call mma_allocate(x,nvar,label='x')

nrow = 1
do m=1,2
  mrow = 1
  if (m == 1) then
    x(:) = r01-r00
  else
    x(:) = r02-r00
  end if

  ! Calculate powers of individual variable values.
  do ivar=1,nvar
    pow = One
    vpow(0,ivar) = One
    do i=1,mxdeg
      pow = pow*x(ivar)
      vpow(i,ivar) = pow
    end do
  end do

  ! Calculate value of each polynomial term at this point.
  do iterm=1,nterm
    ip = mMat(iterm-1,1)
    t = vpow(ip,1)
    do ivar=2,nvar
      ip = mMat(iterm-1,ivar)
      t = t*vpow(ip,ivar)
    end do
    Tmat(nrow,iterm) = t
  end do
  mrow = mrow+1
  nrow = nrow+1

  ! First derivatives.
  if (max_term > 0) then
    do ivar=1,nvar
      do iterm=2,nterm
        irow = mDec(mrow-1,ivar)+1+((m-1)*(nDim/2))
        jvar = nvar
        do while ((mMat(iterm-1,jvar) == 0) .and. (jvar > 1))
          jvar = jvar-1
        end do
        jterm = mDec(iterm-1,jvar)+1
        Tmat(nrow,iterm) = x(jvar)*Tmat(nrow,jterm)
        if (ivar == jvar) then
          Tmat(nrow,iterm) = Tmat(nrow,iterm)+Tmat(irow,jterm)
        end if
      end do
      mrow = mrow+1
      nrow = nrow+1
    end do
  end if

  ! Second derivatives.
  if (max_term > 1) then
    do ivar=1,nvar
      do jvar=ivar,nvar
        do iterm=2,nterm
          irow = mDec(mrow-1,ivar)+1+((m-1)*(nDim/2))
          jrow = mDec(mrow-1,jvar)+1+((m-1)*(nDim/2))
          kvar = nvar
          do while ((mMat(iterm-1,kvar) == 0) .and. (kvar > 1))
            kvar = kvar-1
          end do
          jterm = mDec(iterm-1,kvar)+1
          Tmat(nrow,iterm) = x(kvar)*Tmat(nrow,jterm)
          if (ivar == kvar) then
            Tmat(nrow,iterm) = Tmat(nrow,iterm)+Tmat(irow,jterm)
          end if
          if (jvar == kvar) then
            Tmat(nrow,iterm) = Tmat(nrow,iterm)+Tmat(jrow,jterm)
          end if
        end do
        mrow = mrow+1
        nrow = nrow+1
      end do
    end do
  end if

  ! Third derivatives.
  if (max_term > 2) then
    do ivar=1,nvar
      do jvar=ivar,nvar
        do kvar=jvar,nvar
          do iterm=2,nterm
            irow = mDec(mrow-1,ivar)+1+((m-1)*(nDim/2))
            jrow = mDec(mrow-1,jvar)+1+((m-1)*(nDim/2))
            krow = mDec(mrow-1,kvar)+1+((m-1)*(nDim/2))
            lvar = nvar
            do while ((mMat(iterm-1,lvar) == 0) .and. (lvar > 1))
              lvar = lvar-1
            end do
            jterm = mDec(iterm-1,lvar)+1
            Tmat(nrow,iterm) = x(lvar)*Tmat(nrow,jterm)
            if (ivar == lvar) then
              Tmat(nrow,iterm) = Tmat(nrow,iterm)+Tmat(irow,jterm)
            end if
            if (jvar == lvar) then
              Tmat(nrow,iterm) = Tmat(nrow,iterm)+Tmat(jrow,jterm)
            end if
            if (kvar == lvar) then
              Tmat(nrow,iterm) = Tmat(nrow,iterm)+Tmat(krow,jterm)
            end if
          end do
          mrow = mrow+1
          nrow = nrow+1
        end do
      end do
    end do
  end if

  ! Fourth derivatives.
  if (max_term > 3) then
    do ivar=1,nvar
      do jvar=ivar,nvar
        do kvar=jvar,nvar
          do lvar=kvar,nvar
            do iterm=2,nterm
              irow = mDec(mrow-1,ivar)+1+((m-1)*(nDim/2))
              jrow = mDec(mrow-1,jvar)+1+((m-1)*(nDim/2))
              krow = mDec(mrow-1,kvar)+1+((m-1)*(nDim/2))
              lrow = mDec(mrow-1,lvar)+1+((m-1)*(nDim/2))
              mvar = nvar
              do while ((mMat(iterm-1,mvar) == 0) .and. (mvar > 1))
                mvar = mvar-1
              end do
              jterm = mDec(iterm-1,mvar)+1
              Tmat(nrow,iterm) = x(mvar)*Tmat(nrow,jterm)
              if (ivar == mvar) then
                Tmat(nrow,iterm) = Tmat(nrow,iterm)+Tmat(irow,jterm)
              end if
              if (jvar == mvar) then
                Tmat(nrow,iterm) = Tmat(nrow,iterm)+Tmat(jrow,jterm)
              end if
              if (kvar == mvar) then
                Tmat(nrow,iterm) = Tmat(nrow,iterm)+Tmat(krow,jterm)
              end if
              if (lvar == mvar) then
                Tmat(nrow,iterm) = Tmat(nrow,iterm)+Tmat(lrow,jterm)
              end if
            end do
            mrow = mrow+1
            nrow = nrow+1
          end do
        end do
      end do
    end do
  end if
end do

if (pot) then
  x(1) = r_min(1)-r00(1)
  x(2) = r_min(2)-r00(2)
  x(3) = Two*Pi-r_min(3)-r00(3)

  ! Calculate powers of individual variable values.
  do ivar=1,nvar
    pow = One
    vpow(0,ivar) = One
    do i=1,mxdeg
      pow = pow*x(ivar)
      vpow(i,ivar) = pow
    end do
  end do

  ! Calculate value of each polynomial term at this point.
  do iterm=1,nterm
    ip = mMat(iterm-1,1)
    t = vpow(ip,1)
    do ivar=2,nvar
      ip = mMat(iterm-1,ivar)
      t = t*vpow(ip,ivar)
    end do
    Tmat(nrow,iterm) = t
  end do
end if

! Calculate equation matrix, T(t)*weight*T, and T(t)*weight*rhs.
call mma_allocate(Temp,nterm,nDim,label='Temp')
call mma_allocate(Equmat,nterm,nterm,label='Equmat')

call DGEMM_('T','N',nterm,nDim,nDim,One,Tmat,nDim,weight,nDim,Zero,Temp,nterm)
call DGEMM_('N','N',nterm,nterm,nDim,One,Temp,nterm,Tmat,nDim,Zero,Equmat,nterm)

n = 2
if (max_term > 0) then
  do i=1,nvar
    n = n+1
  end do
  if (max_term > 1) then
    do i=1,nvar
      do j=i,nvar
        n = n+1
      end do
    end do
    if (max_term > 2) then
      do i=1,nvar
        do j=i,nvar
          do k=j,nvar
            Equmat(n,n) = Equmat(n,n)+0.1_wp
            n = n+1
          end do
        end do
      end do
      if (max_term > 3) then
        do i=1,nvar
          do j=i,nvar
            do k=j,nvar
              do l=k,nvar
                Equmat(n,n) = Equmat(n,n)+0.1_wp
                n = n+1
              end do
            end do
          end do
        end do
      end if
    end if
  end if
end if

call DGEMM_('N','N',nterm,1,nDim,One,Temp,nterm,rhs,nDim,Zero,FitCoef,nterm)

! Solve the resulting equation system.
call Dool_MULA(Equmat,nterm,nterm,FitCoef,nterm,nterm,det)

call mma_deallocate(weight)
call mma_deallocate(mDec)
call mma_deallocate(mInc)
call mma_deallocate(vpow)
call mma_deallocate(x)
call mma_deallocate(Tmat)
call mma_deallocate(rhs)
call mma_deallocate(Equmat)
call mma_deallocate(Temp)

end subroutine LSPotFit

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

!end module MatElMod

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

subroutine LSPotFit(r01,energy1,grad1,Hess1,D3_1,D4_1,r02,energy2,grad2,Hess2,D3_2,D4_2,r00,energy0,r_min,FitCoef,mMat,max_term, &
                    pot,nOsc,numcoef)
!  Purpose:
!    Perform a least squares fit of the potential at two different centers, r01 and r02.
!
!  Uses:
!    Constants
!    Linalg
!    FCMod
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

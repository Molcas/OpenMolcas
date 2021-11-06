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

subroutine LSPotFit(r01,energy1,grad1,Hess1,D3_1,D4_1,r02,energy2,grad2,Hess2,D3_2,D4_2,r00,energy0,r_min,FitCoef,mMat,stand_dev, &
                    max_err,use_weight,max_term,pot,nosc,numcoef)
!  Purpose:
!    Perform a least squares fit of the potentiag at two different
!    centra, r01 and r02.
!
!  Uses:
!    Constants
!    Linalg
!    FCMod
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

!use Linalg
!use FCMod
!use TabMod
implicit real*8(a-h,o-z)
#include "Constants_mula.fh"
parameter(mxdeg=6)
real*8 r01(nosc), r02(nosc), r00(nosc), r_min(nosc)
real*8 grad1(nosc), grad2(nosc)
real*8 Hess1(nosc,nosc), Hess2(nosc,nosc)
dimension D3_1(nosc,nosc,nosc), D3_2(nosc,nosc,nosc)
dimension D4_1(nosc,nosc,nosc,nosc), D4_2(nosc,nosc,nosc,nosc)
real*8 FitCoef(numcoef,1)
integer mMat(0:numcoef-1,nosc)
real*8 stand_dev, max_err
logical use_weight, pot
integer nTabDim
#include "WrkSpc.fh"

! Initialize.
nvar = nosc
mTabDim = numcoef-1
call TabDim_drv(max_term,nOsc,nTabDim)
max_mOrd = nTabDim-1
call TabDim_drv(max_term-1,nOsc,nTabDim)
max_mInc = nTabDim-1
n_mDec = mTabDim+1
call GetMem('mDec','Allo','Inte',ipmDec,n_mDec*nvar)
call GetMem('mInc','Allo','Inte',ipmInc,n_mDec*nvar)
l_mMat = nOsc
call MakeTab(max_term,max_mOrd,max_mInc,mMat,iWork(ipmInc),iWork(ipmDec),l_mMat)
nterm = max_mOrd+1

! Set up weight matrix and right hand side.
nDim = 2*nterm
!nDim = nVar*nterm !!!!
if (pot) nDim = nDim+1
call GetMem('rhs','Allo','Real',iprhs,nDim)

call GetMem('weight','Allo','Real',ipweight,nDim*nDim)

call dcopy_(ndim**2,[0.0d0],0,Work(ipweight),1)
n = 1
Work(ipweight+n+nDim*(n-1)-1) = 1.0d4
Work(ipweight+n+nterm+nDim*(n+nterm-1)-1) = 1.0d4
Work(iprhs+n-1) = energy1
Work(iprhs+n+nterm-1) = energy2
n = n+1
if (max_term > 0) then
  do i=1,nvar
    Work(ipweight+n+nDim*(n-1)-1) = 1.0d3
    Work(ipweight+n+nterm+nDim*(n+nterm-1)-1) = 1.0d3
    Work(iprhs+n-1) = grad1(i)
    Work(iprhs+n+nterm-1) = grad2(i)
    n = n+1
  end do
  if (max_term > 1) then
    do i=1,nvar
      do j=i,nvar
        Work(ipweight+n+nDim*(n-1)-1) = 1.0d2
        Work(ipweight+n+nterm+nDim*(n+nterm-1)-1) = 1.0d2
        Work(iprhs+n-1) = Hess1(i,j)
        Work(iprhs+n+nterm-1) = Hess2(i,j)
        n = n+1
      end do
    end do
    if (max_term > 2) then
      do i=1,nvar
        do j=i,nvar
          do k=j,nvar
            Work(ipweight+n+nDim*(n-1)-1) = 1.0d1
            Work(ipweight+n+nterm+nDim*(n+nterm-1)-1) = 1.0d1
            Work(iprhs+n-1) = D3_1(i,j,k)
            Work(iprhs+n+nterm-1) = D3_2(i,j,k)
            n = n+1
          end do
        end do
      end do
      if (max_term > 3) then
        do i=1,nvar
          do j=i,nvar
            do k=j,nvar
              do l=k,nvar
                Work(ipweight+n+nDim*(n-1)-1) = 1.0d0
                Work(ipweight+n+nterm+nDim*(n+nterm-1)-1) = 1.0d0
                Work(iprhs+n-1) = D4_1(i,j,k,l)
                Work(iprhs+n+nterm-1) = D4_2(i,j,k,l)
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
  Work(iprhs+nDim-1) = energy0+20000.0d0/HarToRcm
  Work(ipweight+nDim+nDim*(nDim-1)-1) = 1.0d4
end if

l_vpow = mxdeg+1
call GetMem('vpow','Allo','Real',ipvpow,l_vpow*nvar)
call GetMem('Tmat','Allo','Real',ipTmat,nDim*nterm)
call dcopy_(nDim*nterm,[0.0d0],0,work(ipTmat),1)
!Tmat = 0.0d0
call GetMem('x','Allo','Real',ipx,nvar)

nrow = 1
do m=1,2
  mrow = 1
  if (m == 1) then
    do iv=1,nvar
      Work(ipx+iv-1) = r01(iv)-r00(iv)
    end do
  else
    do iv=1,nvar
      Work(ipx+iv-1) = r02(iv)-r00(iv)
    end do
    !x = r02-r00
  end if

  ! Calculate powers of individual variable values.
  do ivar=1,nvar
    pow = 1.0d0
    work(ipvpow+1+l_vpow*(ivar-1)-1) = 1.0d0
    do i=1,mxdeg
      pow = pow*Work(ipx+ivar-1)
      Work(ipvpow+i+1+l_vpow*(ivar-1)-1) = pow
    end do
  end do

  ! Calculate value of each polynomial term at this point.
  do iterm=1,nterm
    ip = mMat(iterm-1,1)
    t = Work(ipvpow+ip)
    do ivar=2,nvar
      ip = mMat(iterm-1,ivar)
      t = t*Work(ipvpow+ip+1+l_vpow*(ivar-1)-1)
    end do
    Work(ipTmat+nrow+nDim*(iterm-1)-1) = t
  end do
  mrow = mrow+1
  nrow = nrow+1

  ! First derivatives.
  if (max_term > 0) then
    do ivar=1,nvar
      do iterm=2,nterm
        irow = iWork(ipmDec+mrow+n_mDec*(ivar-1)-1)+1+((m-1)*(nDim/2))
        jvar = nvar
        do while ((mMat(iterm-1,jvar) == 0) .and. (jvar > 1))
          jvar = jvar-1
        end do
        jterm = iWork(ipmDec+iterm+n_mDec*(jvar-1)-1)+1
        Work(ipTmat+nrow+nDim*(iterm-1)-1) = work(ipx+jvar-1)*Work(ipTmat+nrow+nDim*(jterm-1)-1)
        if (ivar == jvar) then
          Work(ipTmat+nrow+nDim*(iterm-1)-1) = Work(ipTmat+nrow+nDim*(iterm-1)-1)+Work(ipTmat+irow+nDim*(jterm-1)-1)
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
          irow = iWork(ipmDec+mrow+n_mDec*(ivar-1)-1)+1+((m-1)*(nDim/2))
          jrow = iWork(ipmDec+mrow+n_mDec*(jvar-1)-1)+1+((m-1)*(nDim/2))
          kvar = nvar
          do while ((mMat(iterm-1,kvar) == 0) .and. (kvar > 1))
            kvar = kvar-1
          end do
          jterm = iWork(ipmDec+iterm+n_mDec*(kvar-1)-1)+1
          Work(ipTmat+nrow+nDim*(iterm-1)-1) = Work(ipx+kvar-1)*Work(ipTmat+nrow+nDim*(jterm-1)-1)
          if (ivar == kvar) then
            Work(ipTmat+nrow+nDim*(iterm-1)-1) = Work(ipTmat+nrow+nDim*(iterm-1)-1)+Work(ipTmat+irow+nDim*(jterm-1)-1)
          end if
          if (jvar == kvar) then
            Work(ipTmat+nrow+nDim*(iterm-1)-1) = Work(ipTmat+nrow+nDim*(iterm-1)-1)+Work(ipTmat+jrow+nDim*(jterm-1)-1)
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
            irow = iWork(ipmDec+mrow+n_mDec*(ivar-1)-1)+1+((m-1)*(nDim/2))
            jrow = iWork(ipmDec+mrow+n_mDec*(jvar-1)-1)+1+((m-1)*(nDim/2))
            krow = iWork(ipmDec+mrow+n_mDec*(kvar-1)-1)+1+((m-1)*(nDim/2))
            lvar = nvar
            do while ((mMat(iterm-1,lvar) == 0) .and. (lvar > 1))
              lvar = lvar-1
            end do
            jterm = iWork(ipmDec+iterm+n_mDec*(lvar-1)-1)+1
            Work(ipTmat+nrow+nDim*(iterm-1)-1) = Work(ipx+lvar-1)*Work(ipTmat+nrow+nDim*(jterm-1)-1)
            if (ivar == lvar) then
              Work(ipTmat+nrow+nDim*(iterm-1)-1) = Work(ipTmat+nrow+nDim*(iterm-1)-1)+Work(ipTmat+irow+nDim*(jterm-1)-1)
            end if
            if (jvar == lvar) then
              Work(ipTmat+nrow+nDim*(iterm-1)-1) = Work(ipTmat+nrow+nDim*(iterm-1)-1)+Work(ipTmat+jrow+nDim*(jterm-1)-1)
            end if
            if (kvar == lvar) then
              Work(ipTmat+nrow+nDim*(iterm-1)-1) = Work(ipTmat+nrow+nDim*(iterm-1)-1)+Work(ipTmat+krow+nDim*(jterm-1)-1)
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
              irow = iWork(ipmDec+mrow+n_mDec*(ivar-1)-1)+1+((m-1)*(nDim/2))
              jrow = iWork(ipmDec+mrow+n_mDec*(jvar-1)-1)+1+((m-1)*(nDim/2))
              krow = iWork(ipmDec+mrow+n_mDec*(kvar-1)-1)+1+((m-1)*(nDim/2))
              lrow = iWork(ipmDec+mrow+n_mDec*(lvar-1)-1)+1+((m-1)*(nDim/2))
              mvar = nvar
              do while ((mMat(iterm-1,mvar) == 0) .and. (mvar > 1))
                mvar = mvar-1
              end do
              jterm = iWork(ipmDec+iterm+n_mDec*(mvar-1)-1)+1
              Work(ipTmat+nrow+nDim*(iterm-1)-1) = Work(ipx+mvar-1)*Work(ipTmat+nrow+nDim*(jterm-1)-1)
              if (ivar == mvar) then
                Work(ipTmat+nrow+nDim*(iterm-1)-1) = Work(ipTmat+nrow+nDim*(iterm-1)-1)+Work(ipTmat+irow+ndim*(jterm-1)-1)
              end if
              if (jvar == mvar) then
                work(ipTmat+nrow+nDim*(iterm-1)-1) = Work(ipTmat+nrow+nDim*(iterm-1)-1)+Work(ipTmat+jrow+nDim*(jterm-1)-1)
              end if
              if (kvar == mvar) then
                Work(ipTmat+nrow+nDim*(iterm-1)-1) = Work(ipTmat+nrow+nDim*(iterm-1)-1)+Work(ipTmat+krow+nDim*(jterm-1)-1)
              end if
              if (lvar == mvar) then
                Work(ipTmat+nrow+nDim*(iterm-1)-1) = Work(ipTmat+nrow+nDim*(iterm-1)-1)+Work(ipTmat+lrow+nDim*(jterm-1)-1)
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
  Work(ipx) = r_min(1)-r00(1)
  Work(ipx+1) = r_min(2)-r00(2)
  Work(ipx+2) = 2.0d0*rpi-r_min(3)-r00(3)

  ! Calculate powers of individual variable values.
  do ivar=1,nvar
    pow = 1.0d0
    Work(ipvpow+1+l_vpow*(ivar-1)-1) = 1.0d0
    do i=1,mxdeg
      pow = pow*Work(ipx+ivar-1)
      Work(ipvpow+i+1+l_vpow*(ivar-1)-1) = pow
    end do
  end do

  ! Calculate value of each polynomial term at this point.
  do iterm=1,nterm
    ip = mMat(iterm-1,1)
    t = Work(ipvpow+ip)
    do ivar=2,nvar
      ip = mMat(iterm-1,ivar)
      t = t*Work(ipvpow+ip+1+l_vpow*(ivar-1)-1)
    end do
    Work(ipTmat+nrow+nDim*(iterm-1)-1) = t
  end do
end if

! Calculate equation matrix, T(t)*weight*T, and T(t)*weight*rhs.
call GetMem('Temp','Allo','Real',ipTemp,nterm*nDim)
call GetMem('Equmat','Allo','Real',ipEqumat,nterm*nterm)

call DGEMM_('T','N',nterm,nDim,nDim,1.0d0,Work(ipTmat),nDim,Work(ipweight),nDim,0.0d0,Work(ipTemp),nterm)
call DGEMM_('N','N',nterm,nterm,nDim,1.0d0,Work(ipTemp),nterm,Work(ipTmat),nDim,0.0d0,Work(ipEqumat),nterm)

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
            Work(ipEqumat+n+nterm*(n-1)-1) = Work(ipEqumat+n+nterm*(n-1)-1)+1.0d-1
            n = n+1
          end do
        end do
      end do
      if (max_term > 3) then
        do i=1,nvar
          do j=i,nvar
            do k=j,nvar
              do l=k,nvar
                Work(ipEqumat+n+nterm*(n-1)-1) = Work(ipEqumat+n+nterm*(n-1)-1)+1.0d-1
                n = n+1
              end do
            end do
          end do
        end do
      end if
    end if
  end if
end if

call DGEMM_('N','N',nterm,1,nDim,1.0d0,Work(ipTemp),nterm,Work(iprhs),nDim,0.0d0,FitCoef,nterm)

! Solve the resulting equation system.
call Dool_MULA(Work(ipEqumat),nterm,nterm,FitCoef,nterm,nterm,det)

call GetMem('weight','Free','Real',ipweight,nDim*nDim)
call GetMem('mDec','Free','Inte',ipmDec,n_mDec*nvar)
call GetMem('mInc','Free','Inte',ipmInc,n_mDec*nvar)
call GetMem('vpow','Allo','Real',ipvpow,l_vpow*nvar)
call GetMem('x','Free','Real',ipx,nvar)
call GetMem('Tmat','Free','Real',ipTmat,nDim*nterm)
call GetMem('rhs','Free','Real',iprhs,nDim)
call GetMem('Equmat','Free','Real',ipEqumat,nterm*nterm)
call GetMem('Temp','Free','Real',ipTemp,nterm*nDim)

! Avoid unused argument warnings
if (.false.) then
  call Unused_real(stand_dev)
  call Unused_real(max_err)
  call Unused_logical(use_weight)
end if

end subroutine LSPotFit

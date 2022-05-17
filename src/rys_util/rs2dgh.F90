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
! Copyright (C) 1995, Roland Lindh                                     *
!               1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine Rs2Dgh(xyz2D0,nT,nRys,la,lb,lc,ld,xyz2D1,xyz2D2,IfHss,IndHss,IfGrad,IndGrd,IfG,Coora,Alpha,Beta,Gamma,Delta,nZeta,nEta, &
                  Scrtch,Scrtch2,Temp,Index1,Index2,Index3,Index4,ng,nh,ExpX,ExpY,mZeta,mEta,nIrrep,Tr)
!***********************************************************************
!                                                                      *
! Object:To compute the gradients and the Hessians of the 2D-integrals.*
!                                                                      *
!     Author: Anders Bernhardsson & Roland Lindh,                      *
!             Dept. of Theoretical Chemistry,                          *
!             University of Lund, SWEDEN                               *
!             Februar '95                                              *
!***********************************************************************

implicit real*8(A-H,O-Z)
external ExpX, ExpY
#include "real.fh"
real*8 xyz2D0(nRys*nT,0:la+2,0:lb+2,0:lc+2,0:ld+2,3), xyz2D1(nRys*nT,0:la,0:lb,0:lc,0:ld,3,3), &
       xyz2D2(nRys*nT,0:la,0:lb,0:lc,0:ld,3,6), Coora(3,4), Alpha(nZeta), Beta(nZeta), gamma(nEta), Delta(nEta), Scrtch2(nRys*nT), &
       Scrtch(nRys*nT), Temp(nT)
logical IfGrad(3,4), IfHss(4,3,4,3), IfG(4), EQ, Tr(4)
integer IndGrd(3,4,0:nIrrep-1), Ind1(3), Ind2(3), Ind3(3), Ind4(3), Index2(3,4,4), Index1(3,4), Index3(3,3), Index4(2,6,3), ng(3), &
        nh(3), IndHss(4,3,4,3,0:nIrrep-1)
#ifdef NAGFOR
save Ind1, Ind2, Ind3, Ind4
#endif

nx = 0
ny = 0
nz = 0
mx = 0
my = 0
mz = 0
call ICopy(12,[0],0,Index1,1)
call ICopy(48,[0],0,Index2,1)

! Differentiate with respect to the first center

if (IfG(1)) then
  call ExpX(Temp,mZeta,mEta,Alpha,sqrt(Two))
  call Exp_2(Scrtch,nRys,nT,Temp,sqrt(Two))
  nVec = 0
  if (IfGrad(1,1)) then
    nx = nx+1
    nVec = nVec+1
    Ind1(nVec) = nx
    Ind2(nVec) = 1
    Index1(1,1) = nx
    Index3(nx,1) = 1
  end if
  if (IfGrad(2,1)) then
    ny = ny+1
    nVec = nVec+1
    Ind1(nVec) = ny
    Ind2(nVec) = 2
    Index1(2,1) = ny
    Index3(ny,2) = 1
  end if
  if (IfGrad(3,1)) then
    nz = nz+1
    nVec = nVec+1
    Ind1(nVec) = nz
    Ind2(nVec) = 3
    Index1(3,1) = nz
    Index3(nz,3) = 1
  end if
  do i=nvec+1,3
    Ind1(i) = 0
  end do

  mvec = 0
  if (IfHss(1,1,1,1)) then
    mx = mx+1
    mvec = mvec+1
    Ind3(mVec) = mx
    Ind4(mVec) = 1
    Index2(1,1,1) = mx
    Index4(1,mx,1) = 1
    Index4(2,mx,1) = 1
  end if
  if (IfHss(1,2,1,2)) then
    my = my+1
    mvec = mvec+1
    Ind3(mVec) = my
    Ind4(mVec) = 2
    Index2(2,1,1) = my
    Index4(1,my,2) = 1
    Index4(2,my,2) = 1
  end if
  if (IfHss(1,3,1,3)) then
    mz = mz+1
    mvec = mvec+1
    Ind3(mVec) = mz
    Ind4(mVec) = 3
    Index2(3,1,1) = mz
    Index4(1,mz,3) = 1
    Index4(2,mz,3) = 1
  end if
  do i=mvec+1,3
    Ind3(i) = 0
  end do
  nvecx = max(nvec,mvec)
  if (nVecx /= 0) then

    ! Here we go with center 1

    do n=1,nvec
      do id=0,ld
        do ic=0,lc
          do ib=0,lb
            ra = -One
            do ia=0,la
              ra = ra+Two
              do iVec=1,nRys*nT
                xyz2D1(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) = Scrtch(iVec)*xyz2D0(iVec,ia+1,ib,ic,id,Ind2(n))
              end do
            end do
          end do
        end do
      end do
    end do
    do n=1,mvec
      do id=0,ld
        do ic=0,lc
          do ib=0,lb
            ra = -One
            do ia=0,la
              ra = ra+Two
              do iVec=1,nRys*nT
                xyz2D2(iVec,ia,ib,ic,id,Ind4(n),Ind3(n)) = Scrtch(iVec)**2*xyz2D0(iVec,ia+2,ib,ic,id,Ind4(n))- &
                                                           ra*Scrtch(iVec)*xyz2D0(iVec,ia,ib,ic,id,Ind4(n))
              end do
            end do
          end do
        end do
      end do
    end do
    if (la >= 1) then
      do n=1,nvec
        do id=0,ld
          do ic=0,lc
            do ib=0,lb
              ra = Zero
              do ia=1,la
                ra = ra+One
                call DaXpY_inline(nRys*nT,-ra,xyz2D0(1,ia-1,ib,ic,id,Ind2(n)),xyz2D1(1,ia,ib,ic,id,Ind2(n),Ind1(n)))
              end do
            end do
          end do
        end do
      end do
    end if
    if (la >= 2) then
      do n=1,mvec
        do id=0,ld
          do ic=0,lc
            do ib=0,lb
              ra = One
              do ia=2,la
                ra = ra+One
                Fact = ra*ra-ra
                call Daxpy_inline(nRys*nT,Fact,xyz2D0(1,ia-2,ib,ic,id,Ind4(n)),xyz2D2(1,ia,ib,ic,id,Ind4(n),Ind3(n)))
              end do
            end do
          end do
        end do
      end do
    end if
  end if
end if

! Cross term center 1 and center 2
if (IfG(2)) then
  call ExpX(Temp,mZeta,mEta,Beta,sqrt(Two))
  call Exp_2(Scrtch2,nRys,nT,Temp,sqrt(Two))
end if
if (IfG(2) .and. Ifg(1)) then
  nVec = 0
  if (ifHss(2,1,1,1)) then
    mx = mx+1
    nVec = nVec+1
    Ind1(nvec) = mx
    Ind2(nVec) = 1
    Index2(1,2,1) = mx
    Index4(1,mx,1) = 2
    Index4(2,mx,1) = 1
  end if
  if (ifHss(2,2,1,2)) then
    my = my+1
    nVec = nVec+1
    Ind1(nvec) = my
    Ind2(nVec) = 2
    Index2(2,2,1) = my
    Index4(1,my,2) = 2
    Index4(2,my,2) = 1
  end if
  if (ifHss(2,3,1,3)) then
    mz = mz+1
    nVec = nVec+1
    Ind1(nvec) = mz
    Ind2(nVec) = 3
    Index2(3,2,1) = mz
    Index4(1,mz,3) = 2
    Index4(2,mz,3) = 1
  end if
  do i=nVec+1,3
    Ind1(i) = 0
  end do
  if (nvec /= 0) then
    do n=1,nvec
      do id=0,ld
        do ic=0,lc
          do ib=0,lb
            do ia=0,la
              do iVec=1,nRys*nT
                xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) = Scrtch(iVec)*Scrtch2(iVec)*xyz2D0(iVec,ia+1,ib+1,ic,id,Ind2(n))
              end do
            end do
          end do
        end do
      end do
    end do
    if (la >= 1) then
      do n=1,nvec
        do id=0,ld
          do ic=0,lc
            do ib=0,lb
              ra = Zero
              do ia=1,la
                ra = ra+One
                do iVec=1,nRys*nT
                  xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n))- &
                                                             ra*Scrtch2(iVec)*xyz2D0(iVec,ia-1,ib+1,ic,id,Ind2(n))
                end do
              end do
            end do
          end do
        end do
      end do
    end if
    if (lb >= 1) then
      do n=1,nvec
        do id=0,ld
          do ic=0,lc
            rb = Zero
            do ib=1,lb
              rb = rb+One
              do ia=0,la
                do iVec=1,nRys*nT
                  xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n))- &
                                                             rb*Scrtch(iVec)*xyz2D0(iVec,ia+1,ib-1,ic,id,Ind2(n))
                end do
              end do
            end do
          end do
        end do
      end do
    end if
    if ((la >= 1) .and. (lb >= 1)) then
      do n=1,nvec
        do id=0,ld
          do ic=0,lc
            rb = Zero
            do ib=1,lb
              rb = rb+One
              ra = Zero
              do ia=1,la
                ra = ra+One
                Fact = ra*rb
                call DaXpY_inline(nRys*nT,Fact,xyz2D0(1,ia-1,ib-1,ic,id,Ind2(n)),xyz2D2(1,ia,ib,ic,id,Ind2(n),Ind1(n)))
              end do
            end do
          end do
        end do
      end do
    end if
  end if
end if

! Differentiate with respect to the second center

if (IfG(2)) then
  call ExpX(Temp,mZeta,mEta,Beta,sqrt(Two))
  call Exp_2(Scrtch2,nRys,nT,Temp,sqrt(Two))
  nVec = 0
  if (IfGrad(1,2)) then
    nx = nx+1
    nVec = nVec+1
    Ind1(nVec) = nx
    Ind2(nVec) = 1
    Index1(1,2) = nx
    Index3(nx,1) = 2
  end if
  if (IfGrad(2,2)) then
    ny = ny+1
    nVec = nVec+1
    Ind1(nVec) = ny
    Ind2(nVec) = 2
    Index3(ny,2) = 2
    Index1(2,2) = ny
  end if
  if (IfGrad(3,2)) then
    nz = nz+1
    nVec = nVec+1
    Ind1(nVec) = nz
    Ind2(nVec) = 3
    Index1(3,2) = nz
    Index3(nz,3) = 2
  end if
  do i=nvec+1,3
    Ind1(i) = 0
  end do

  mvec = 0
  if (IfHss(2,1,2,1)) then
    mx = mx+1
    mvec = mvec+1
    Ind3(mVec) = mx
    Ind4(mVec) = 1
    Index2(1,2,2) = mx
    Index4(1,mx,1) = 2
    Index4(2,mx,1) = 2
  end if
  if (IfHss(2,2,2,2)) then
    my = my+1
    mvec = mvec+1
    Ind3(mVec) = my
    Ind4(mvec) = 2
    Index2(2,2,2) = my
    Index4(1,my,2) = 2
    Index4(2,my,2) = 2
  end if
  if (IfHss(2,3,2,3)) then
    mz = mz+1
    mvec = mvec+1
    Ind3(mVec) = mz
    Ind4(mVec) = 3
    Index2(3,2,2) = mz
    Index4(1,mz,3) = 2
    Index4(2,mz,3) = 2
  end if
  do i=mvec+1,3
    Ind3(i) = 0
  end do
  nvecx = max(nvec,mvec)
  if (nVecx /= 0) then

    do n=1,nVec
      do id=0,ld
        do ic=0,lc
          rb = -One
          do ib=0,lb
            rb = rb+Two
            do ia=0,la
              do iVec=1,nRys*nT
                xyz2D1(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) = Scrtch2(iVec)*xyz2D0(iVec,ia,ib+1,ic,id,Ind2(n))
              end do
            end do
          end do
        end do
      end do
    end do
    do n=1,mVec
      do id=0,ld
        do ic=0,lc
          rb = -One
          do ib=0,lb
            rb = rb+Two
            do ia=0,la
              do iVec=1,nRys*nT
                xyz2D2(iVec,ia,ib,ic,id,Ind4(n),Ind3(n)) = Scrtch2(iVec)**2*xyz2D0(iVec,ia,ib+2,ic,id,Ind4(n))- &
                                                           rb*Scrtch2(iVec)*xyz2D0(iVec,ia,ib,ic,id,Ind4(n))
              end do
            end do
          end do
        end do
      end do
    end do

    if (lb >= 1) then
      do n=1,nvec
        do id=0,ld
          do ic=0,lc
            rb = Zero
            do ib=1,lb
              rb = rb+One
              do ia=0,la
                call DaXpy_inline(nRys*nT,-rb,xyz2D0(1,ia,ib-1,ic,id,Ind2(n)),xyz2D1(1,ia,ib,ic,id,Ind2(n),Ind1(n)))
              end do
            end do
          end do
        end do
      end do
    end if
    if (lb >= 2) then
      do n=1,mvec
        do id=0,ld
          do ic=0,lc
            rb = One
            do ib=2,lb
              rb = rb+One
              Fact = rb*rb-rb
              do ia=0,la
                call DaXpy_inline(nRys*nT,Fact,xyz2D0(1,ia,ib-2,ic,id,Ind4(n)),xyz2D2(1,ia,ib,ic,id,Ind4(n),Ind3(n)))
              end do
            end do
          end do
        end do
      end do
    end if
  end if
end if

! Cross Term center 2 and 3

if (IfG(2) .and. IfG(3)) then
  call ExpX(Temp,mZeta,mEta,Beta,sqrt(Two))
  call Exp_2(Scrtch2,nRys,nT,Temp,sqrt(Two))
  call ExpY(Temp,mZeta,mEta,Gamma,sqrt(Two))
  call Exp_2(Scrtch,nRys,nT,Temp,sqrt(Two))
  nVec = 0
  if (ifHss(3,1,2,1)) then
    mx = mx+1
    nVec = nVec+1
    Ind1(nvec) = mx
    Ind2(nVec) = 1
    Index2(1,3,2) = mx
    Index4(1,mx,1) = 3
    Index4(2,mx,1) = 2
  end if
  if (ifHss(3,2,2,2)) then
    my = my+1
    nVec = nVec+1
    Ind1(nvec) = my
    Ind2(nVec) = 2
    Index2(2,3,2) = my
    Index4(1,my,2) = 3
    Index4(2,my,2) = 2
  end if
  if (ifHss(3,3,2,3)) then
    mz = mz+1
    nVec = nVec+1
    Ind1(nvec) = mz
    Ind2(nVec) = 3
    Index2(3,3,2) = mz
    Index4(1,mz,3) = 3
    Index4(2,mz,3) = 2
  end if
  do i=nVec+1,3
    Ind1(i) = 0
  end do
  if (nvec /= 0) then

    do n=1,nvec
      do id=0,ld
        do ic=0,lc
          do ib=0,lb
            do ia=0,la
              do iVec=1,nRys*nT
                xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) = Scrtch(iVec)*Scrtch2(iVec)*xyz2D0(iVec,ia,ib+1,ic+1,id,Ind2(n))
              end do
            end do
          end do
        end do
      end do
    end do
    if (lb >= 1) then
      do n=1,nvec
        do id=0,ld
          do ic=0,lc
            rb = Zero
            do ib=1,lb
              rb = rb+One
              do ia=0,la
                do iVec=1,nRys*nT
                  xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n))- &
                                                             rb*Scrtch(iVec)*xyz2D0(iVec,ia,ib-1,ic+1,id,Ind2(n))
                end do
              end do
            end do
          end do
        end do
      end do
    end if
    if (lc >= 1) then
      do n=1,nvec
        do id=0,ld
          rc = Zero
          do ic=1,lc
            rc = rc+One
            do ib=0,lb
              do ia=0,la
                do iVec=1,nRys*nT
                  xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n))- &
                                                             rc*Scrtch2(iVec)*xyz2D0(iVec,ia,ib+1,ic-1,id,Ind2(n))
                end do
              end do
            end do
          end do
        end do
      end do
    end if
    if ((lb >= 1) .and. (lc >= 1)) then
      do n=1,nvec
        do id=0,ld
          rc = Zero
          do ic=1,lc
            rc = rc+One
            rb = Zero
            do ib=1,lb
              rb = rb+One
              Fact = rb*rc
              do ia=0,la
                call DaxPy_inline(nRys*nT,Fact,xyz2D0(1,ia,ib-1,ic-1,id,Ind2(n)),xyz2D2(1,ia,ib,ic,id,Ind2(n),Ind1(n)))
              end do
            end do
          end do
        end do
      end do
    end if
  end if
end if

! Differentiate with respect to the third center

if (IfG(3)) then
  call ExpY(Temp,mZeta,mEta,Gamma,sqrt(Two))
  call Exp_2(Scrtch,nRys,nT,Temp,sqrt(Two))
  nvec = 0
  if (IfGrad(1,3)) then
    nx = nx+1
    nVec = nVec+1
    Ind1(nVec) = nx
    Ind2(nVec) = 1
    Index1(1,3) = nx
    Index3(nx,1) = 3
  end if
  if (IfGrad(2,3)) then
    ny = ny+1
    nVec = nVec+1
    Ind1(nVec) = ny
    Ind2(nVec) = 2
    Index1(2,3) = ny
    Index3(ny,2) = 3
  end if
  if (IfGrad(3,3)) then
    nz = nz+1
    nVec = nVec+1
    Ind1(nVec) = nz
    Ind2(nVec) = 3
    Index1(3,3) = nz
    Index3(nz,3) = 3
  end if
  do i=nvec+1,3
    Ind1(i) = 0
  end do
!
  mvec = 0
  if (IfHss(3,1,3,1)) then
    mx = mx+1
    mvec = mvec+1
    Ind3(mVec) = mx
    Ind4(mVec) = 1
    Index2(1,3,3) = mx
    Index4(1,mx,1) = 3
    Index4(2,mx,1) = 3
  end if
  if (IfHss(3,2,3,2)) then
    my = my+1
    mvec = mvec+1
    Ind3(mVec) = my
    Ind4(mVec) = 2
    Index2(2,3,3) = my
    Index4(1,my,2) = 3
    Index4(2,my,2) = 3
  end if
  if (IfHss(3,3,3,3)) then
    mz = mz+1
    mvec = mvec+1
    Ind3(mVec) = mz
    Ind4(mVec) = 3
    Index2(3,3,3) = mz
    Index4(1,mz,3) = 3
    Index4(2,mz,3) = 3
  end if
  do i=mvec+1,3
    Ind3(i) = 0
  end do
  nvecx = max(nvec,mvec)
  if (nVecx /= 0) then
    do n=1,nvec
      do id=0,ld
        rc = -One
        do ic=0,lc
          rc = rc+Two
          do ib=0,lb
            do ia=0,la
              do iVec=1,nRys*nT
                xyz2D1(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) = Scrtch(iVec)*xyz2D0(iVec,ia,ib,ic+1,id,Ind2(n))
              end do
            end do
          end do
        end do
      end do
    end do
    do n=1,mvec
      do id=0,ld
        rc = -One
        do ic=0,lc
          rc = rc+Two
          do ib=0,lb
            do ia=0,la
              do iVec=1,nRys*nT
                xyz2D2(iVec,ia,ib,ic,id,Ind4(n),Ind3(n)) = Scrtch(iVec)**2*xyz2D0(iVec,ia,ib,ic+2,id,Ind4(n))- &
                                                           rc*Scrtch(iVec)*xyz2D0(iVec,ia,ib,ic,id,Ind4(n))
              end do
            end do
          end do
        end do
      end do
    end do

    if (lc >= 1) then
      do n=1,nvec
        do id=0,ld
          rc = Zero
          do ic=1,lc
            rc = rc+One
            do ib=0,lb
              do ia=0,la
                call DaXpY_inline(nRys*nT,-rc,xyz2D0(1,ia,ib,ic-1,id,Ind2(n)),xyz2D1(1,ia,ib,ic,id,Ind2(n),Ind1(n)))
              end do
            end do
          end do
        end do
      end do
    end if
    if (lc >= 2) then
      do n=1,mvec
        do id=0,ld
          rc = One
          do ic=2,lc
            rc = rc+One
            Fact = rc*rc-rc
            do ib=0,lb
              do ia=0,la
                call DaXpY_inline(nt*nrys,Fact,xyz2D0(1,ia,ib,ic-2,id,Ind4(n)),xyz2D2(1,ia,ib,ic,id,Ind4(n),Ind3(n)))
              end do
            end do
          end do
        end do
      end do
    end if
  end if
end if

! Cross term 1 3

if (IfG(1) .and. IfG(3)) then
  call ExpX(Temp,mZeta,mEta,Alpha,sqrt(Two))
  call Exp_2(Scrtch2,nRys,nT,Temp,sqrt(Two))
  call ExpY(Temp,mZeta,mEta,Gamma,sqrt(Two))
  call Exp_2(Scrtch,nRys,nT,Temp,sqrt(Two))
  nVec = 0
  if (ifHss(3,1,1,1)) then
    mx = mx+1
    nVec = nVec+1
    Ind1(nvec) = mx
    Ind2(nVec) = 1
    Index4(1,mx,1) = 3
    Index4(2,mx,1) = 1
    Index2(1,3,1) = mx
  end if
  if (ifHss(3,2,1,2)) then
    my = my+1
    nVec = nVec+1
    Ind1(nvec) = my
    Ind2(nVec) = 2
    Index4(1,my,2) = 3
    Index4(2,my,2) = 1
    Index2(2,3,1) = my
  end if
  if (ifHss(3,3,1,3)) then
    mz = mz+1
    nVec = nVec+1
    Ind1(nvec) = mz
    Ind2(nVec) = 3
    Index4(1,mz,3) = 3
    Index4(2,mz,3) = 1
    Index2(3,3,1) = mz
  end if
  do i=nVec+1,3
    Ind1(i) = 0
  end do
  if (nVec /= 0) then
    do n=1,nvec
      do id=0,ld
        do ic=0,lc
          do ib=0,lb
            do ia=0,la
              do iVec=1,nRys*nT
                xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) = Scrtch(iVec)*Scrtch2(iVec)*xyz2D0(iVec,ia+1,ib,ic+1,id,Ind2(n))
              end do
            end do
          end do
        end do
      end do
    end do
    if (la >= 1) then
      do n=1,nvec
        do id=0,ld
          do ic=0,lc
            do ib=0,lb
              ra = Zero
              do ia=1,la
                ra = ra+One
                do iVec=1,nRys*nT
                  xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n))- &
                                                             ra*Scrtch(iVec)*xyz2D0(iVec,ia-1,ib,ic+1,id,Ind2(n))
                end do
              end do
            end do
          end do
        end do
      end do
    end if
    if (lc >= 1) then
      do n=1,nVec
        do id=0,ld
          rc = Zero
          do ic=1,lc
            rc = rc+One
            do ib=0,lb
              do ia=0,la
                do iVec=1,nRys*nT
                  xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n))- &
                                                             rc*Scrtch2(iVec)*xyz2D0(iVec,ia+1,ib,ic-1,id,Ind2(n))
                end do
              end do
            end do
          end do
        end do
      end do
    end if
    if ((la >= 1) .and. (lc >= 1)) then
      do n=1,nvec
        do id=0,ld
          rc = Zero
          do ic=1,lc
            rc = rc+One
            do ib=0,lb
              ra = Zero
              do ia=1,la
                ra = ra+One
                Fact = rc*ra
                call DaXpy_inline(nt*nrys,Fact,xyz2D0(1,ia-1,ib,ic-1,id,Ind2(n)),xyz2D2(1,ia,ib,ic,id,Ind2(n),Ind1(n)))
              end do
            end do
          end do
        end do
      end do
    end if
  end if
end if

! 1 4

if (IfG(4)) then
  call ExpY(Temp,mZeta,mEta,Delta,sqrt(Two))
  call Exp_2(Scrtch,nRys,nT,Temp,sqrt(Two))
  if (IfG(1)) then
    call ExpX(Temp,mZeta,mEta,Alpha,sqrt(Two))
    call Exp_2(Scrtch2,nRys,nT,Temp,sqrt(Two))
    nVec = 0
    if (ifHss(4,1,1,1)) then
      mx = mx+1
      nVec = nVec+1
      Ind1(nvec) = mx
      Ind2(nVec) = 1
      Index2(1,4,1) = mx
      Index4(1,mx,1) = 4
      Index4(2,mx,1) = 1
    end if
    if (ifHss(4,2,1,2)) then
      my = my+1
      nVec = nVec+1
      Ind1(nvec) = my
      Ind2(nVec) = 2
      Index2(2,4,1) = my
      Index4(1,my,2) = 4
      Index4(2,my,2) = 1
    end if
    if (ifHss(4,3,1,3)) then
      mz = mz+1
      nVec = nVec+1
      Ind1(nvec) = mz
      Ind2(nVec) = 3
      Index2(3,4,1) = mz
      Index4(1,mz,3) = 4
      Index4(2,mz,3) = 1
    end if
    do i=nVec+1,3
      Ind1(i) = 0
    end do

    if (nVec /= 0) then

      do n=1,nvec
        do id=0,ld
          do ic=0,lc
            do ib=0,lb
              do ia=0,la
                do iVec=1,nRys*nT
                  xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) = Scrtch(iVec)*Scrtch2(iVec)*xyz2D0(iVec,ia+1,ib,ic,id+1,Ind2(n))
                end do
              end do
            end do
          end do
        end do
      end do

      if (la >= 1) then
        do n=1,nvec
          do id=0,ld
            do ic=0,lc
              do ib=0,lb
                ra = Zero
                do ia=1,la
                  ra = ra+One
                  do iVec=1,nRys*nT
                    xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n))- &
                                                               ra*Scrtch(iVec)*xyz2D0(iVec,ia-1,ib,ic,id+1,Ind2(n))
                  end do
                end do
              end do
            end do
          end do
        end do
      end if

      if (ld >= 1) then
        do n=1,nvec
          rd = Zero
          do id=1,ld
            rd = rd+One
            do ic=0,lc
              do ib=0,lb
                do ia=0,la
                  do iVec=1,nRys*nT
                    xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n))- &
                                                               rd*Scrtch2(iVec)*xyz2D0(iVec,ia+1,ib,ic,id-1,Ind2(n))
                  end do
                end do
              end do
            end do
          end do
        end do
      end if
      if ((la >= 1) .and. (ld >= 1)) then
        do n=1,nvec
          rd = Zero
          do id=1,ld
            rd = rd+One
            do ic=0,lc
              do ib=0,lb
                ra = Zero
                do ia=1,la
                  ra = ra+One
                  Fact = rd*ra
                  call DaxPy_inline(nRys*nT,Fact,xyz2D0(1,ia-1,ib,ic,id-1,Ind2(n)),xyz2D2(1,ia,ib,ic,id,Ind2(n),Ind1(n)))
                end do
              end do
            end do
          end do
        end do
      end if
    end if
  end if

  ! Cross terms between 2 4

  if (IfG(2) .and. Ifg(4)) then
    call ExpX(Temp,mZeta,mEta,Beta,sqrt(Two))
    call Exp_2(Scrtch2,nRys,nT,Temp,sqrt(Two))
    nVec = 0
    if (ifHss(4,1,2,1)) then
      mx = mx+1
      nVec = nVec+1
      Ind1(nvec) = mx
      Ind2(nVec) = 1
      Index2(1,4,2) = mx
      Index4(1,mx,1) = 4
      Index4(2,mx,1) = 2
    end if
    if (ifHss(4,2,2,2)) then
      my = my+1
      nVec = nVec+1
      Ind1(nvec) = my
      Ind2(nVec) = 2
      Index2(2,4,2) = my
      Index4(1,my,2) = 4
      Index4(2,my,2) = 2
    end if
    if (ifHss(4,3,2,3)) then
      mz = mz+1
      nVec = nVec+1
      Ind1(nvec) = mz
      Ind2(nVec) = 3
      Index2(3,4,2) = mz
      Index4(1,mz,3) = 4
      Index4(2,mz,3) = 2
    end if
    do i=nVec+1,3
      Ind1(i) = 0
    end do

    if (nvec /= 0) then

      do n=1,nvec
        do id=0,ld
          do ic=0,lc
            do ib=0,lb
              do ia=0,la
                do iVec=1,nRys*nT
                  xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) = Scrtch(iVec)*Scrtch2(iVec)*xyz2D0(iVec,ia,ib+1,ic,id+1,Ind2(n))
                end do
              end do
            end do
          end do
        end do
      end do

      if (lb >= 1) then
        do n=1,nvec
          do id=0,ld
            do ic=0,lc
              rb = Zero
              do ib=1,lb
                rb = rb+One
                do ia=0,la
                  do iVec=1,nRys*nT
                    xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n))- &
                                                               rb*Scrtch(iVec)*xyz2D0(iVec,ia,ib-1,ic,id+1,Ind2(n))
                  end do
                end do
              end do
            end do
          end do
        end do
      end if

      if (ld >= 1) then
        do n=1,nvec
          rd = Zero
          do id=1,ld
            rd = rd+One
            do ic=0,lc
              do ib=0,lb
                do ia=0,la
                  do iVec=1,nRys*nT
                    xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n))- &
                                                               rd*Scrtch2(iVec)*xyz2D0(iVec,ia,ib+1,ic,id-1,Ind2(n))
                  end do
                end do
              end do
            end do
          end do
        end do
      end if
      if ((lb >= 1) .and. (ld >= 1)) then
        do n=1,nvec
          rd = Zero
          do id=1,ld
            rd = rd+One
            do ic=0,lc
              rb = Zero
              do ib=1,lb
                rb = rb+One
                do ia=0,la
                  Fact = rb*rd
                  call DaxPy_inline(nt*nrys,Fact,xyz2D0(1,ia,ib-1,ic,id-1,Ind2(n)),xyz2D2(1,ia,ib,ic,id,Ind2(n),Ind1(n)))
                end do
              end do
            end do
          end do
        end do
      end if
    end if
  end if

  ! Cross Term 3 4

  if (IfG(3) .and. IfG(4)) then
    call ExpY(Temp,mZeta,mEta,Gamma,sqrt(Two))
    call Exp_2(Scrtch2,nRys,nT,Temp,sqrt(Two))
    nVec = 0
    if (ifHss(4,1,3,1)) then
      mx = mx+1
      nVec = nVec+1
      Ind1(nvec) = mx
      Ind2(nVec) = 1
      Index2(1,4,3) = mx
      Index4(1,mx,1) = 4
      Index4(2,mx,1) = 3
    end if
    if (ifHss(4,2,3,2)) then
      my = my+1
      nVec = nVec+1
      Ind1(nvec) = my
      Ind2(nVec) = 2
      Index2(2,4,3) = my
      Index4(1,my,2) = 4
      Index4(2,my,2) = 3

    end if
    if (ifHss(4,3,3,3)) then
      mz = mz+1
      nVec = nVec+1
      Ind1(nvec) = mz
      Ind2(nVec) = 3
      Index2(3,4,3) = mz
      Index4(1,mz,3) = 4
      Index4(2,mz,3) = 3
    end if
    do i=nVec+1,3
      Ind1(i) = 0
    end do

    if (nvec /= 0) then

      do n=1,nvec
        do id=0,ld
          do ic=0,lc
            do ib=0,lb
              do ia=0,la
                do iVec=1,nRys*nT
                  xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) = Scrtch(iVec)*Scrtch2(iVec)*xyz2D0(iVec,ia,ib,ic+1,id+1,Ind2(n))
                end do
              end do
            end do
          end do
        end do
      end do

      if (lc >= 1) then
        do n=1,nvec
          do id=0,ld
            rc = Zero
            do ic=1,lc
              rc = rc+One
              do ib=0,lb
                do ia=0,la
                  do iVec=1,nRys*nT
                    xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n))- &
                                                               rc*Scrtch(iVec)*xyz2D0(iVec,ia,ib,ic-1,id+1,Ind2(n))
                  end do
                end do
              end do
            end do
          end do
        end do
      end if
      if (ld >= 1) then
        do n=1,nvec
          rd = Zero
          do id=1,ld
            rd = rd+One
            do ic=0,lc
              do ib=0,lb
                do ia=0,la
                  do iVec=1,nRys*nT
                    xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n))- &
                                                               rd*Scrtch2(iVec)*xyz2D0(iVec,ia,ib,ic+1,id-1,Ind2(n))
                  end do
                end do
              end do
            end do
          end do
        end do
      end if
      if ((lc >= 1) .and. (ld >= 1)) then
        do n=1,nvec
          rd = Zero
          do id=1,ld
            rd = rd+One
            rc = Zero
            do ic=1,lc
              rc = rc+One
              do ib=0,lb
                do ia=0,la
                  Fact = rc*rd
                  call DaxPy_inline(nt*nRys,fact,xyz2D0(1,ia,ib,ic-1,id-1,Ind2(n)),xyz2D2(1,ia,ib,ic,id,Ind2(n),Ind1(n)))
                end do
              end do
            end do
          end do
        end do
      end if
    end if
  end if

  ! Differentiate with respect to the fourth center

  nvec = 0
  if (IfGrad(1,4)) then
    nx = nx+1
    nVec = nVec+1
    Ind1(nVec) = nx
    Ind2(nVec) = 1
    Index1(1,4) = nx
    Index3(nx,1) = 4
  end if
  if (IfGrad(2,4)) then
    ny = ny+1
    nVec = nVec+1
    Ind1(nVec) = ny
    Ind2(nVec) = 2
    Index1(2,4) = ny
    Index3(ny,2) = 4
  end if
  if (IfGrad(3,4)) then
    nz = nz+1
    nVec = nVec+1
    Ind1(nVec) = nz
    Ind2(nVec) = 3
    Index1(3,4) = nz
    Index3(nz,3) = 4
  end if
  do i=nvec+1,3
    Ind1(i) = 0
  end do

  mvec = 0
  if (IfHss(4,1,4,1)) then
    mx = mx+1
    mvec = mvec+1
    Ind3(mVec) = mx
    Ind4(mVec) = 1
    Index2(1,4,4) = mx
    Index4(1,mx,1) = 4
    Index4(2,mx,1) = 4
  end if
  if (IfHss(4,2,4,2)) then
    my = my+1
    mvec = mvec+1
    Ind3(mVec) = my
    Ind4(mVec) = 2
    Index2(2,4,4) = my
    Index4(1,my,2) = 4
    Index4(2,my,2) = 4
  end if
  if (IfHss(4,3,4,3)) then
    mz = mz+1
    mvec = mvec+1
    Ind3(mVec) = mz
    Ind4(mVec) = 3
    Index2(3,4,4) = mz
    Index4(1,mz,3) = 4
    Index4(2,mz,3) = 4
  end if
  do i=mvec+1,3
    Ind3(i) = 0
  end do
  nvecx = max(nvec,mvec)

  if (nVecx /= 0) then

    do n=1,nvec
      rd = -One
      do id=0,ld
        rd = rd+Two
        do ic=0,lc
          do ib=0,lb
            do ia=0,la
              do iVec=1,nRys*nT
                xyz2D1(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) = Scrtch(iVec)*xyz2D0(iVec,ia,ib,ic,id+1,Ind2(n))
              end do
            end do
          end do
        end do
      end do
    end do
    do n=1,mvec
      rd = -One
      do id=0,ld
        rd = rd+Two
        do ic=0,lc
          do ib=0,lb
            do ia=0,la
              do iVec=1,nRys*nT
                xyz2D2(iVec,ia,ib,ic,id,Ind4(n),Ind3(n)) = Scrtch(iVec)**2*xyz2D0(iVec,ia,ib,ic,id+2,Ind4(n))- &
                                                           rd*Scrtch(iVec)*xyz2D0(iVec,ia,ib,ic,id,Ind4(n))
              end do
            end do
          end do
        end do
      end do
    end do
    if (ld >= 1) then
      do n=1,nvec
        rd = Zero
        do id=1,ld
          rd = rd+One
          do ic=0,lc
            do ib=0,lb
              do ia=0,la
                call Daxpy_inline(nt*nrys,-rd,xyz2D0(1,ia,ib,ic,id-1,Ind2(n)),xyz2D1(1,ia,ib,ic,id,Ind2(n),Ind1(n)))
              end do
            end do
          end do
        end do
      end do
    end if
    if (ld >= 2) then
      do n=1,mvec
        rd = One
        do id=2,ld
          rd = rd+One
          Fact = rd*rd-rd
          do ic=0,lc
            do ib=0,lb
              do ia=0,la
                call Daxpy_inline(nt*nrys,Fact,xyz2D0(1,ia,ib,ic,id-2,Ind4(n)),xyz2D2(1,ia,ib,ic,id,Ind4(n),Ind3(n)))
              end do
            end do
          end do
        end do
      end do
    end if

  end if
end if

! Sum over common centers

do iCent=1,3
  if (IfG(iCent)) then
    do jCent=iCent+1,4
      if (EQ(Coora(1,iCent),Coora(1,jCent))) then
        if (IfG(jCent)) then
          do iCar=1,3
            i1 = Index2(iCar,iCent,iCent)
            i2 = Index2(iCar,jCent,jCent)
            i3 = Index2(iCar,jCent,iCent)
            j4 = Index1(iCar,jCent)
            j5 = Index1(iCar,iCent)
            if (IfHss(jCent,iCar,jCent,iCar) .and. IfHss(iCent,iCar,iCent,iCar)) then
              call DaXpY_inline(nRys*nT*(la+1)*(lb+1)*(lc+1)*(ld+1),One,xyz2D2(1,0,0,0,0,iCar,i2),xyz2D2(1,0,0,0,0,iCar,i1))
            end if
            if (IfHss(jCent,iCar,iCent,iCar) .and. IfHss(iCent,iCar,iCent,iCar)) then
              call DaXpY_inline(nRys*nT*(la+1)*(lb+1)*(lc+1)*(ld+1),Two,xyz2D2(1,0,0,0,0,iCar,i3),xyz2D2(1,0,0,0,0,iCar,i1))
            end if
            if ((j4 /= 0) .and. (j5 /= 0) .and. ifgrad(iCar,iCent) .and. ifgrad(iCar,jCent)) &
              call DaXpY_inline(nRys*nT*(la+1)*(lb+1)*(lc+1)*(ld+1),One,xyz2D1(1,0,0,0,0,iCar,j4),xyz2D1(1,0,0,0,0,iCar,j5))
            do kCent=1,4
              if (IfG(kCent)) then
                if ((kcent /= iCent) .and. (kcent /= jCent)) then
                  if (ifHss(kCent,iCar,jCent,iCar) .or. ifHss(jCent,iCar,kCent,iCar)) then
                    i4 = Index2(iCar,max(kCent,jCent),min(jCent,kCent))
                    i5 = Index2(iCar,max(kCent,iCent),min(iCent,kCent))
                    call DaXpY_inline(nRys*nT*(la+1)*(lb+1)*(lc+1)*(ld+1),One,xyz2D2(1,0,0,0,0,iCar,i4),xyz2D2(1,0,0,0,0,iCar,i5))
                  end if
                end if
              end if
            end do ! kCent
          end do ! iCar

          IfG(jCent) = .false.
          Tr(jCent) = .false.
          do jCar=1,3
            IfGrad(jcar,jCent) = .false.
            do iIrrep=0,nIrrep-1
              IndGrd(jCar,jcent,iIrrep) = 0
            end do
            do kCent=1,4
              do kCar=1,3
                IfHss(jCent,jCar,kCent,kCar) = .false.
                IfHss(kCent,kCar,jCent,jCar) = .false.
                do iIrrep=0,nIrrep-1
                  IndHss(jCent,jCar,kCent,kCar,iIrrep) = 0
                  IndHss(kCent,kCar,jCent,jCar,iIrrep) = 0
                end do
              end do
            end do
          end do

        end if
      end if ! end eq
    end do ! jCent
  end if
end do ! iCent

nh(1) = mx
nh(2) = my
nh(3) = mz
ng(1) = nx
ng(2) = ny
ng(3) = nz

return

end subroutine Rs2Dgh

subroutine Daxpy_inline(nt,r,A,B)

implicit real*8(A-H,O-Z)
real*8 A(*), B(*)

do i=1,nt
  B(i) = B(i)+r*A(i)
end do

return

end subroutine Daxpy_inline

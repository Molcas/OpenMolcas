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

subroutine Rs2Dgh(xyz2D0,nT,nRys,la,lb,lc,ld,xyz2D1,xyz2D2,IfHss,IndHss,IfGrad,IndGrd,IfG,Coora,Alpha,Beta,Gmma,Delta,nZeta,nEta, &
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

use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nT, nRys, la, lb, lc, ld, nZeta, nEta, mZeta, mEta, nIrrep
real(kind=wp), intent(in) :: xyz2D0(nRys*nT,0:la+2,0:lb+2,0:lc+2,0:ld+2,3), Coora(3,4), Alpha(nZeta), Beta(nZeta), Gmma(nEta), &
                             Delta(nEta)
real(kind=wp), intent(out) :: xyz2D1(nRys*nT,0:la,0:lb,0:lc,0:ld,3,3), xyz2D2(nRys*nT,0:la,0:lb,0:lc,0:ld,3,6), Scrtch(nRys*nT), &
                              Scrtch2(nRys*nT), Temp(nT)
logical(kind=iwp), intent(inout) :: IfHss(4,3,4,3), IfGrad(3,4), IfG(4), Tr(4)
integer(kind=iwp), intent(inout) :: IndHss(4,3,4,3,0:nIrrep-1), IndGrd(3,4,0:nIrrep-1)
integer(kind=iwp), intent(out) :: Index1(3,4), Index2(3,4,4), Index3(3,3), Index4(2,6,3), ng(3), nh(3)
external :: ExpX, ExpY
integer(kind=iwp) :: i1, i2, i3, i4, i5, ia, ib, ic, iCar, iCent, id, Ind1(3), Ind2(3), Ind3(3), Ind4(3), j4, j5, jCent, kCent, &
                     mVec, mx, my, mz, n, nVec, nvecx, nx, ny, nz
real(kind=wp) :: Fact, ra, rb, rc, rd
logical(kind=iwp), external :: EQ

nx = 0
ny = 0
nz = 0
mx = 0
my = 0
mz = 0
Index1(:,:) = 0
Index2(:,:,:) = 0

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
  Ind1(nVec+1:) = 0

  mVec = 0
  if (IfHss(1,1,1,1)) then
    mx = mx+1
    mVec = mVec+1
    Ind3(mVec) = mx
    Ind4(mVec) = 1
    Index2(1,1,1) = mx
    Index4(1,mx,1) = 1
    Index4(2,mx,1) = 1
  end if
  if (IfHss(1,2,1,2)) then
    my = my+1
    mVec = mVec+1
    Ind3(mVec) = my
    Ind4(mVec) = 2
    Index2(2,1,1) = my
    Index4(1,my,2) = 1
    Index4(2,my,2) = 1
  end if
  if (IfHss(1,3,1,3)) then
    mz = mz+1
    mVec = mVec+1
    Ind3(mVec) = mz
    Ind4(mVec) = 3
    Index2(3,1,1) = mz
    Index4(1,mz,3) = 1
    Index4(2,mz,3) = 1
  end if
  Ind3(mVec+1:) = 0
  nvecx = max(nVec,mVec)
  if (nVecx /= 0) then

    ! Here we go with center 1

    do n=1,nVec
      do id=0,ld
        do ic=0,lc
          do ib=0,lb
            ra = -One
            do ia=0,la
              ra = ra+Two
              xyz2D1(:,ia,ib,ic,id,Ind2(n),Ind1(n)) = Scrtch(:)*xyz2D0(:,ia+1,ib,ic,id,Ind2(n))
            end do
          end do
        end do
      end do
    end do
    do n=1,mVec
      do id=0,ld
        do ic=0,lc
          do ib=0,lb
            ra = -One
            do ia=0,la
              ra = ra+Two
              xyz2D2(:,ia,ib,ic,id,Ind4(n),Ind3(n)) = Scrtch(:)**2*xyz2D0(:,ia+2,ib,ic,id,Ind4(n))- &
                                                      ra*Scrtch(:)*xyz2D0(:,ia,ib,ic,id,Ind4(n))
            end do
          end do
        end do
      end do
    end do
    if (la >= 1) then
      do n=1,nVec
        do id=0,ld
          do ic=0,lc
            do ib=0,lb
              ra = Zero
              do ia=1,la
                ra = ra+One
                xyz2D1(:,ia,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D1(:,ia,ib,ic,id,Ind2(n),Ind1(n))-ra*xyz2D0(:,ia-1,ib,ic,id,Ind2(n))
              end do
            end do
          end do
        end do
      end do
    end if
    if (la >= 2) then
      do n=1,mVec
        do id=0,ld
          do ic=0,lc
            do ib=0,lb
              ra = One
              do ia=2,la
                ra = ra+One
                Fact = ra*ra-ra
                xyz2D2(:,ia,ib,ic,id,Ind4(n),Ind3(n)) = xyz2D2(:,ia,ib,ic,id,Ind4(n),Ind3(n))+Fact*xyz2D0(:,ia-2,ib,ic,id,Ind4(n))
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
if (IfG(2) .and. IfG(1)) then
  nVec = 0
  if (IfHss(2,1,1,1)) then
    mx = mx+1
    nVec = nVec+1
    Ind1(nVec) = mx
    Ind2(nVec) = 1
    Index2(1,2,1) = mx
    Index4(1,mx,1) = 2
    Index4(2,mx,1) = 1
  end if
  if (IfHss(2,2,1,2)) then
    my = my+1
    nVec = nVec+1
    Ind1(nVec) = my
    Ind2(nVec) = 2
    Index2(2,2,1) = my
    Index4(1,my,2) = 2
    Index4(2,my,2) = 1
  end if
  if (IfHss(2,3,1,3)) then
    mz = mz+1
    nVec = nVec+1
    Ind1(nVec) = mz
    Ind2(nVec) = 3
    Index2(3,2,1) = mz
    Index4(1,mz,3) = 2
    Index4(2,mz,3) = 1
  end if
  Ind1(nVec+1:) = 0
  if (nVec /= 0) then
    do n=1,nVec
      do id=0,ld
        do ic=0,lc
          do ib=0,lb
            do ia=0,la
              xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n)) = Scrtch(:)*Scrtch2(:)*xyz2D0(:,ia+1,ib+1,ic,id,Ind2(n))
            end do
          end do
        end do
      end do
    end do
    if (la >= 1) then
      do n=1,nVec
        do id=0,ld
          do ic=0,lc
            do ib=0,lb
              ra = Zero
              do ia=1,la
                ra = ra+One
                xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n))- &
                                                        ra*Scrtch2(:)*xyz2D0(:,ia-1,ib+1,ic,id,Ind2(n))
              end do
            end do
          end do
        end do
      end do
    end if
    if (lb >= 1) then
      do n=1,nVec
        do id=0,ld
          do ic=0,lc
            rb = Zero
            do ib=1,lb
              rb = rb+One
              do ia=0,la
                xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n))- &
                                                        rb*Scrtch(:)*xyz2D0(:,ia+1,ib-1,ic,id,Ind2(n))
              end do
            end do
          end do
        end do
      end do
    end if
    if ((la >= 1) .and. (lb >= 1)) then
      do n=1,nVec
        do id=0,ld
          do ic=0,lc
            rb = Zero
            do ib=1,lb
              rb = rb+One
              ra = Zero
              do ia=1,la
                ra = ra+One
                Fact = ra*rb
                xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n))+Fact*xyz2D0(:,ia-1,ib-1,ic,id,Ind2(n))
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
  Ind1(nVec+1:) = 0

  mVec = 0
  if (IfHss(2,1,2,1)) then
    mx = mx+1
    mVec = mVec+1
    Ind3(mVec) = mx
    Ind4(mVec) = 1
    Index2(1,2,2) = mx
    Index4(1,mx,1) = 2
    Index4(2,mx,1) = 2
  end if
  if (IfHss(2,2,2,2)) then
    my = my+1
    mVec = mVec+1
    Ind3(mVec) = my
    Ind4(mVec) = 2
    Index2(2,2,2) = my
    Index4(1,my,2) = 2
    Index4(2,my,2) = 2
  end if
  if (IfHss(2,3,2,3)) then
    mz = mz+1
    mVec = mVec+1
    Ind3(mVec) = mz
    Ind4(mVec) = 3
    Index2(3,2,2) = mz
    Index4(1,mz,3) = 2
    Index4(2,mz,3) = 2
  end if
  Ind3(mVec+1:) = 0
  nvecx = max(nVec,mVec)
  if (nVecx /= 0) then

    do n=1,nVec
      do id=0,ld
        do ic=0,lc
          rb = -One
          do ib=0,lb
            rb = rb+Two
            do ia=0,la
              xyz2D1(:,ia,ib,ic,id,Ind2(n),Ind1(n)) = Scrtch2(:)*xyz2D0(:,ia,ib+1,ic,id,Ind2(n))
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
              xyz2D2(:,ia,ib,ic,id,Ind4(n),Ind3(n)) = Scrtch2(:)**2*xyz2D0(:,ia,ib+2,ic,id,Ind4(n))- &
                                                      rb*Scrtch2(:)*xyz2D0(:,ia,ib,ic,id,Ind4(n))
            end do
          end do
        end do
      end do
    end do

    if (lb >= 1) then
      do n=1,nVec
        do id=0,ld
          do ic=0,lc
            rb = Zero
            do ib=1,lb
              rb = rb+One
              xyz2D1(:,:,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D1(:,:,ib,ic,id,Ind2(n),Ind1(n))-rb*xyz2D0(:,0:la,ib-1,ic,id,Ind2(n))
            end do
          end do
        end do
      end do
    end if
    if (lb >= 2) then
      do n=1,mVec
        do id=0,ld
          do ic=0,lc
            rb = One
            do ib=2,lb
              rb = rb+One
              Fact = rb*rb-rb
              xyz2D2(:,:,ib,ic,id,Ind4(n),Ind3(n)) = xyz2D2(:,:,ib,ic,id,Ind4(n),Ind3(n))+Fact*xyz2D0(:,0:la,ib-2,ic,id,Ind4(n))
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
  call ExpY(Temp,mZeta,mEta,Gmma,sqrt(Two))
  call Exp_2(Scrtch,nRys,nT,Temp,sqrt(Two))
  nVec = 0
  if (IfHss(3,1,2,1)) then
    mx = mx+1
    nVec = nVec+1
    Ind1(nVec) = mx
    Ind2(nVec) = 1
    Index2(1,3,2) = mx
    Index4(1,mx,1) = 3
    Index4(2,mx,1) = 2
  end if
  if (IfHss(3,2,2,2)) then
    my = my+1
    nVec = nVec+1
    Ind1(nVec) = my
    Ind2(nVec) = 2
    Index2(2,3,2) = my
    Index4(1,my,2) = 3
    Index4(2,my,2) = 2
  end if
  if (IfHss(3,3,2,3)) then
    mz = mz+1
    nVec = nVec+1
    Ind1(nVec) = mz
    Ind2(nVec) = 3
    Index2(3,3,2) = mz
    Index4(1,mz,3) = 3
    Index4(2,mz,3) = 2
  end if
  Ind1(nVec+1:) = 0
  if (nVec /= 0) then

    do n=1,nVec
      do id=0,ld
        do ic=0,lc
          do ib=0,lb
            do ia=0,la
              xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n)) = Scrtch(:)*Scrtch2(:)*xyz2D0(:,ia,ib+1,ic+1,id,Ind2(n))
            end do
          end do
        end do
      end do
    end do
    if (lb >= 1) then
      do n=1,nVec
        do id=0,ld
          do ic=0,lc
            rb = Zero
            do ib=1,lb
              rb = rb+One
              do ia=0,la
                xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n))- &
                                                        rb*Scrtch(:)*xyz2D0(:,ia,ib-1,ic+1,id,Ind2(n))
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
                xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n))- &
                                                        rc*Scrtch2(:)*xyz2D0(:,ia,ib+1,ic-1,id,Ind2(n))
              end do
            end do
          end do
        end do
      end do
    end if
    if ((lb >= 1) .and. (lc >= 1)) then
      do n=1,nVec
        do id=0,ld
          rc = Zero
          do ic=1,lc
            rc = rc+One
            rb = Zero
            do ib=1,lb
              rb = rb+One
              Fact = rb*rc
              xyz2D2(:,:,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D2(:,:,ib,ic,id,Ind2(n),Ind1(n))+Fact*xyz2D0(:,0:la,ib-1,ic-1,id,Ind2(n))
            end do
          end do
        end do
      end do
    end if
  end if
end if

! Differentiate with respect to the third center

if (IfG(3)) then
  call ExpY(Temp,mZeta,mEta,Gmma,sqrt(Two))
  call Exp_2(Scrtch,nRys,nT,Temp,sqrt(Two))
  nVec = 0
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
  Ind1(nVec+1:) = 0

  mVec = 0
  if (IfHss(3,1,3,1)) then
    mx = mx+1
    mVec = mVec+1
    Ind3(mVec) = mx
    Ind4(mVec) = 1
    Index2(1,3,3) = mx
    Index4(1,mx,1) = 3
    Index4(2,mx,1) = 3
  end if
  if (IfHss(3,2,3,2)) then
    my = my+1
    mVec = mVec+1
    Ind3(mVec) = my
    Ind4(mVec) = 2
    Index2(2,3,3) = my
    Index4(1,my,2) = 3
    Index4(2,my,2) = 3
  end if
  if (IfHss(3,3,3,3)) then
    mz = mz+1
    mVec = mVec+1
    Ind3(mVec) = mz
    Ind4(mVec) = 3
    Index2(3,3,3) = mz
    Index4(1,mz,3) = 3
    Index4(2,mz,3) = 3
  end if
  Ind3(mVec+1:) = 0
  nvecx = max(nVec,mVec)
  if (nVecx /= 0) then
    do n=1,nVec
      do id=0,ld
        rc = -One
        do ic=0,lc
          rc = rc+Two
          do ib=0,lb
            do ia=0,la
              xyz2D1(:,ia,ib,ic,id,Ind2(n),Ind1(n)) = Scrtch(:)*xyz2D0(:,ia,ib,ic+1,id,Ind2(n))
            end do
          end do
        end do
      end do
    end do
    do n=1,mVec
      do id=0,ld
        rc = -One
        do ic=0,lc
          rc = rc+Two
          do ib=0,lb
            do ia=0,la
              xyz2D2(:,ia,ib,ic,id,Ind4(n),Ind3(n)) = Scrtch(:)**2*xyz2D0(:,ia,ib,ic+2,id,Ind4(n))- &
                                                      rc*Scrtch(:)*xyz2D0(:,ia,ib,ic,id,Ind4(n))
            end do
          end do
        end do
      end do
    end do

    if (lc >= 1) then
      do n=1,nVec
        do id=0,ld
          rc = Zero
          do ic=1,lc
            rc = rc+One
            xyz2D1(:,:,:,ic,id,Ind2(n),Ind1(n)) = xyz2D1(:,:,:,ic,id,Ind2(n),Ind1(n))-rc*xyz2D0(:,0:la,0:lb,ic-1,id,Ind2(n))
          end do
        end do
      end do
    end if
    if (lc >= 2) then
      do n=1,mVec
        do id=0,ld
          rc = One
          do ic=2,lc
            rc = rc+One
            Fact = rc*rc-rc
            xyz2D2(:,:,:,ic,id,Ind4(n),Ind3(n)) = xyz2D2(:,:,:,ic,id,Ind4(n),Ind3(n))+Fact*xyz2D0(:,0:la,0:lb,ic-2,id,Ind4(n))
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
  call ExpY(Temp,mZeta,mEta,Gmma,sqrt(Two))
  call Exp_2(Scrtch,nRys,nT,Temp,sqrt(Two))
  nVec = 0
  if (IfHss(3,1,1,1)) then
    mx = mx+1
    nVec = nVec+1
    Ind1(nVec) = mx
    Ind2(nVec) = 1
    Index4(1,mx,1) = 3
    Index4(2,mx,1) = 1
    Index2(1,3,1) = mx
  end if
  if (IfHss(3,2,1,2)) then
    my = my+1
    nVec = nVec+1
    Ind1(nVec) = my
    Ind2(nVec) = 2
    Index4(1,my,2) = 3
    Index4(2,my,2) = 1
    Index2(2,3,1) = my
  end if
  if (IfHss(3,3,1,3)) then
    mz = mz+1
    nVec = nVec+1
    Ind1(nVec) = mz
    Ind2(nVec) = 3
    Index4(1,mz,3) = 3
    Index4(2,mz,3) = 1
    Index2(3,3,1) = mz
  end if
  Ind1(nVec+1:) = 0
  if (nVec /= 0) then
    do n=1,nVec
      do id=0,ld
        do ic=0,lc
          do ib=0,lb
            do ia=0,la
              xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n)) = Scrtch(:)*Scrtch2(:)*xyz2D0(:,ia+1,ib,ic+1,id,Ind2(n))
            end do
          end do
        end do
      end do
    end do
    if (la >= 1) then
      do n=1,nVec
        do id=0,ld
          do ic=0,lc
            do ib=0,lb
              ra = Zero
              do ia=1,la
                ra = ra+One
                xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n))- &
                                                        ra*Scrtch(:)*xyz2D0(:,ia-1,ib,ic+1,id,Ind2(n))
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
                xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n))- &
                                                        rc*Scrtch2(:)*xyz2D0(:,ia+1,ib,ic-1,id,Ind2(n))
              end do
            end do
          end do
        end do
      end do
    end if
    if ((la >= 1) .and. (lc >= 1)) then
      do n=1,nVec
        do id=0,ld
          rc = Zero
          do ic=1,lc
            rc = rc+One
            do ib=0,lb
              ra = Zero
              do ia=1,la
                ra = ra+One
                Fact = rc*ra
                xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n))+Fact*xyz2D0(:,ia-1,ib,ic-1,id,Ind2(n))
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
    if (IfHss(4,1,1,1)) then
      mx = mx+1
      nVec = nVec+1
      Ind1(nVec) = mx
      Ind2(nVec) = 1
      Index2(1,4,1) = mx
      Index4(1,mx,1) = 4
      Index4(2,mx,1) = 1
    end if
    if (IfHss(4,2,1,2)) then
      my = my+1
      nVec = nVec+1
      Ind1(nVec) = my
      Ind2(nVec) = 2
      Index2(2,4,1) = my
      Index4(1,my,2) = 4
      Index4(2,my,2) = 1
    end if
    if (IfHss(4,3,1,3)) then
      mz = mz+1
      nVec = nVec+1
      Ind1(nVec) = mz
      Ind2(nVec) = 3
      Index2(3,4,1) = mz
      Index4(1,mz,3) = 4
      Index4(2,mz,3) = 1
    end if
    Ind1(nVec+1:) = 0

    if (nVec /= 0) then

      do n=1,nVec
        do id=0,ld
          do ic=0,lc
            do ib=0,lb
              do ia=0,la
                xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n)) = Scrtch(:)*Scrtch2(:)*xyz2D0(:,ia+1,ib,ic,id+1,Ind2(n))
              end do
            end do
          end do
        end do
      end do

      if (la >= 1) then
        do n=1,nVec
          do id=0,ld
            do ic=0,lc
              do ib=0,lb
                ra = Zero
                do ia=1,la
                  ra = ra+One
                  xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n))- &
                                                          ra*Scrtch(:)*xyz2D0(:,ia-1,ib,ic,id+1,Ind2(n))
                end do
              end do
            end do
          end do
        end do
      end if

      if (ld >= 1) then
        do n=1,nVec
          rd = Zero
          do id=1,ld
            rd = rd+One
            do ic=0,lc
              do ib=0,lb
                do ia=0,la
                  xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n))- &
                                                          rd*Scrtch2(:)*xyz2D0(:,ia+1,ib,ic,id-1,Ind2(n))
                end do
              end do
            end do
          end do
        end do
      end if
      if ((la >= 1) .and. (ld >= 1)) then
        do n=1,nVec
          rd = Zero
          do id=1,ld
            rd = rd+One
            do ic=0,lc
              do ib=0,lb
                ra = Zero
                do ia=1,la
                  ra = ra+One
                  Fact = rd*ra
                  xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n))+ &
                                                          Fact*xyz2D0(:,ia-1,ib,ic,id-1,Ind2(n))
                end do
              end do
            end do
          end do
        end do
      end if
    end if
  end if

  ! Cross terms between 2 4

  if (IfG(2) .and. IfG(4)) then
    call ExpX(Temp,mZeta,mEta,Beta,sqrt(Two))
    call Exp_2(Scrtch2,nRys,nT,Temp,sqrt(Two))
    nVec = 0
    if (IfHss(4,1,2,1)) then
      mx = mx+1
      nVec = nVec+1
      Ind1(nVec) = mx
      Ind2(nVec) = 1
      Index2(1,4,2) = mx
      Index4(1,mx,1) = 4
      Index4(2,mx,1) = 2
    end if
    if (IfHss(4,2,2,2)) then
      my = my+1
      nVec = nVec+1
      Ind1(nVec) = my
      Ind2(nVec) = 2
      Index2(2,4,2) = my
      Index4(1,my,2) = 4
      Index4(2,my,2) = 2
    end if
    if (IfHss(4,3,2,3)) then
      mz = mz+1
      nVec = nVec+1
      Ind1(nVec) = mz
      Ind2(nVec) = 3
      Index2(3,4,2) = mz
      Index4(1,mz,3) = 4
      Index4(2,mz,3) = 2
    end if
    Ind1(nVec+1:) = 0

    if (nVec /= 0) then

      do n=1,nVec
        do id=0,ld
          do ic=0,lc
            do ib=0,lb
              do ia=0,la
                xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n)) = Scrtch(:)*Scrtch2(:)*xyz2D0(:,ia,ib+1,ic,id+1,Ind2(n))
              end do
            end do
          end do
        end do
      end do

      if (lb >= 1) then
        do n=1,nVec
          do id=0,ld
            do ic=0,lc
              rb = Zero
              do ib=1,lb
                rb = rb+One
                do ia=0,la
                  xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n))- &
                                                          rb*Scrtch(:)*xyz2D0(:,ia,ib-1,ic,id+1,Ind2(n))
                end do
              end do
            end do
          end do
        end do
      end if

      if (ld >= 1) then
        do n=1,nVec
          rd = Zero
          do id=1,ld
            rd = rd+One
            do ic=0,lc
              do ib=0,lb
                do ia=0,la
                  xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n))- &
                                                          rd*Scrtch2(:)*xyz2D0(:,ia,ib+1,ic,id-1,Ind2(n))
                end do
              end do
            end do
          end do
        end do
      end if
      if ((lb >= 1) .and. (ld >= 1)) then
        do n=1,nVec
          rd = Zero
          do id=1,ld
            rd = rd+One
            do ic=0,lc
              rb = Zero
              do ib=1,lb
                rb = rb+One
                Fact = rb*rd
                xyz2D2(:,:,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D2(:,:,ib,ic,id,Ind2(n),Ind1(n))+Fact*xyz2D0(:,0:la,ib-1,ic,id-1,Ind2(n))
              end do
            end do
          end do
        end do
      end if
    end if
  end if

  ! Cross Term 3 4

  if (IfG(3) .and. IfG(4)) then
    call ExpY(Temp,mZeta,mEta,Gmma,sqrt(Two))
    call Exp_2(Scrtch2,nRys,nT,Temp,sqrt(Two))
    nVec = 0
    if (IfHss(4,1,3,1)) then
      mx = mx+1
      nVec = nVec+1
      Ind1(nVec) = mx
      Ind2(nVec) = 1
      Index2(1,4,3) = mx
      Index4(1,mx,1) = 4
      Index4(2,mx,1) = 3
    end if
    if (IfHss(4,2,3,2)) then
      my = my+1
      nVec = nVec+1
      Ind1(nVec) = my
      Ind2(nVec) = 2
      Index2(2,4,3) = my
      Index4(1,my,2) = 4
      Index4(2,my,2) = 3

    end if
    if (IfHss(4,3,3,3)) then
      mz = mz+1
      nVec = nVec+1
      Ind1(nVec) = mz
      Ind2(nVec) = 3
      Index2(3,4,3) = mz
      Index4(1,mz,3) = 4
      Index4(2,mz,3) = 3
    end if
    Ind1(nVec+1:) = 0

    if (nVec /= 0) then

      do n=1,nVec
        do id=0,ld
          do ic=0,lc
            do ib=0,lb
              do ia=0,la
                xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n)) = Scrtch(:)*Scrtch2(:)*xyz2D0(:,ia,ib,ic+1,id+1,Ind2(n))
              end do
            end do
          end do
        end do
      end do

      if (lc >= 1) then
        do n=1,nVec
          do id=0,ld
            rc = Zero
            do ic=1,lc
              rc = rc+One
              do ib=0,lb
                do ia=0,la
                  xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n))- &
                                                          rc*Scrtch(:)*xyz2D0(:,ia,ib,ic-1,id+1,Ind2(n))
                end do
              end do
            end do
          end do
        end do
      end if
      if (ld >= 1) then
        do n=1,nVec
          rd = Zero
          do id=1,ld
            rd = rd+One
            do ic=0,lc
              do ib=0,lb
                do ia=0,la
                  xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n)) = xyz2D2(:,ia,ib,ic,id,Ind2(n),Ind1(n))- &
                                                          rd*Scrtch2(:)*xyz2D0(:,ia,ib,ic+1,id-1,Ind2(n))
                end do
              end do
            end do
          end do
        end do
      end if
      if ((lc >= 1) .and. (ld >= 1)) then
        do n=1,nVec
          rd = Zero
          do id=1,ld
            rd = rd+One
            rc = Zero
            do ic=1,lc
              rc = rc+One
              Fact = rc*rd
              xyz2D2(:,:,:,ic,id,Ind2(n),Ind1(n)) = xyz2D2(:,:,:,ic,id,Ind2(n),Ind1(n))+Fact*xyz2D0(:,0:la,0:lb,ic-1,id-1,Ind2(n))
            end do
          end do
        end do
      end if
    end if
  end if

  ! Differentiate with respect to the fourth center

  nVec = 0
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
  Ind1(nVec+1:) = 0

  mVec = 0
  if (IfHss(4,1,4,1)) then
    mx = mx+1
    mVec = mVec+1
    Ind3(mVec) = mx
    Ind4(mVec) = 1
    Index2(1,4,4) = mx
    Index4(1,mx,1) = 4
    Index4(2,mx,1) = 4
  end if
  if (IfHss(4,2,4,2)) then
    my = my+1
    mVec = mVec+1
    Ind3(mVec) = my
    Ind4(mVec) = 2
    Index2(2,4,4) = my
    Index4(1,my,2) = 4
    Index4(2,my,2) = 4
  end if
  if (IfHss(4,3,4,3)) then
    mz = mz+1
    mVec = mVec+1
    Ind3(mVec) = mz
    Ind4(mVec) = 3
    Index2(3,4,4) = mz
    Index4(1,mz,3) = 4
    Index4(2,mz,3) = 4
  end if
  Ind3(mVec+1:) = 0
  nvecx = max(nVec,mVec)

  if (nVecx /= 0) then

    do n=1,nVec
      rd = -One
      do id=0,ld
        rd = rd+Two
        do ic=0,lc
          do ib=0,lb
            do ia=0,la
              xyz2D1(:,ia,ib,ic,id,Ind2(n),Ind1(n)) = Scrtch(:)*xyz2D0(:,ia,ib,ic,id+1,Ind2(n))
            end do
          end do
        end do
      end do
    end do
    do n=1,mVec
      rd = -One
      do id=0,ld
        rd = rd+Two
        do ic=0,lc
          do ib=0,lb
            do ia=0,la
              xyz2D2(:,ia,ib,ic,id,Ind4(n),Ind3(n)) = Scrtch(:)**2*xyz2D0(:,ia,ib,ic,id+2,Ind4(n))- &
                                                      rd*Scrtch(:)*xyz2D0(:,ia,ib,ic,id,Ind4(n))
            end do
          end do
        end do
      end do
    end do
    if (ld >= 1) then
      do n=1,nVec
        rd = Zero
        do id=1,ld
          rd = rd+One
          xyz2D1(:,:,:,:,id,Ind2(n),Ind1(n)) = xyz2D1(:,:,:,:,id,Ind2(n),Ind1(n))-rd*xyz2D0(:,0:la,0:lb,0:lc,id-1,Ind2(n))
        end do
      end do
    end if
    if (ld >= 2) then
      do n=1,mVec
        rd = One
        do id=2,ld
          rd = rd+One
          Fact = rd*rd-rd
          xyz2D2(:,:,:,:,id,Ind4(n),Ind3(n)) = xyz2D2(:,:,:,:,id,Ind4(n),Ind3(n))+Fact*xyz2D0(:,0:la,0:lb,0:lc,id-2,Ind4(n))
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
            if (IfHss(jCent,iCar,jCent,iCar) .and. IfHss(iCent,iCar,iCent,iCar)) &
              xyz2D2(:,:,:,:,:,iCar,i1) = xyz2D2(:,:,:,:,:,iCar,i1)+xyz2D2(:,:,:,:,:,iCar,i2)
            if (IfHss(jCent,iCar,iCent,iCar) .and. IfHss(iCent,iCar,iCent,iCar)) &
              xyz2D2(:,:,:,:,:,iCar,i1) = xyz2D2(:,:,:,:,:,iCar,i1)+Two*xyz2D2(:,:,:,:,:,iCar,i3)
            if ((j4 /= 0) .and. (j5 /= 0) .and. IfGrad(iCar,iCent) .and. IfGrad(iCar,jCent)) &
              xyz2D1(:,:,:,:,:,iCar,j5) = xyz2D1(:,:,:,:,:,iCar,j5)+xyz2D1(:,:,:,:,:,iCar,j4)
            do kCent=1,4
              if (IfG(kCent)) then
                if ((kCent /= iCent) .and. (kCent /= jCent)) then
                  if (IfHss(kCent,iCar,jCent,iCar) .or. IfHss(jCent,iCar,kCent,iCar)) then
                    i4 = Index2(iCar,max(kCent,jCent),min(jCent,kCent))
                    i5 = Index2(iCar,max(kCent,iCent),min(iCent,kCent))
                    xyz2D2(:,:,:,:,:,iCar,i5) = xyz2D2(:,:,:,:,:,iCar,i5)+xyz2D2(:,:,:,:,:,iCar,i4)
                  end if
                end if
              end if
            end do ! kCent
          end do ! iCar

          IfG(jCent) = .false.
          Tr(jCent) = .false.
          IfGrad(:,jCent) = .false.
          IndGrd(:,jCent,:) = 0
          IfHss(jCent,:,:,:) = .false.
          IfHss(:,:,jCent,:) = .false.
          IndHss(jCent,:,:,:,:) = 0
          IndHss(:,:,jCent,:,:) = 0

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

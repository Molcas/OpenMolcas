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
! Copyright (C) 1991,1992, Roland Lindh                                *
!***********************************************************************

subroutine HRR2Da_mck(Arr1,nVec,nabMax,ncdMax,Arr2,A,B,la,lb,lc,ld,IfHss,IfGrd)
!***********************************************************************
!                                                                      *
! Object: to apply the transfer equation to the 2D-integrals.          *
!         The transformation is in place and the recursion             *
!         is replaced with the indentity when applicable.              *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             September '91                                            *
!             Modified to recurrence algorithm, February '92.          *
!             Improved algorithm, June '92.                            *
!***********************************************************************

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nVec, nabMax, ncdMax, la, lb, lc, ld
real(kind=wp), intent(in) :: Arr1(nVec,3,0:nabMax,0:ncdMax), A(3), B(3)
real(kind=wp), intent(out) :: Arr2(nVec,0:la+2,0:lb+2,0:ncdMax,3)
logical(kind=iwp), intent(in) :: IfHss(4,3,4,3), IfGrd(3,4)
integer(kind=iwp) :: ia, iab, ib, iCar, icd, ja, jb, ka, kb, lla, llb, llcd, ma, mb
real(kind=wp) :: AB

!iRout = 233
!iPrint = nPrint(iRout)

do iCar=1,3
  lla = 0
  llb = 0
  llcd = 0
  if (IfGrd(iCar,1)) lla = max(1,lla)
  if (IfHss(1,iCar,1,iCar)) lla = 2
  if (IfGrd(iCar,2)) llb = max(1,llb)
  if (IfHss(2,iCar,2,iCar)) llb = 2
  if (IfGrd(iCar,3) .or. IfGrd(iCar,4)) llcd = max(1,llcd)
  if (IfHss(3,iCar,3,iCar) .or. IfHss(4,iCar,4,iCar)) llcd = 2
  AB = A(iCar)-B(iCar)
  if (AB == Zero) then
    do icd=0,lc+ld+llcd
      ! Using the identity
      do ia=0,la+lla
        do ib=0,lb+llb
          iab = ia+ib
          if (iab > la+lb+max(lla,llb)) cycle
          Arr2(:,ia,ib,icd,iCar) = Arr1(:,iCar,iab,icd)
        end do
      end do
    end do
  else if (la >= lb) then
    do icd=0,lc+ld+llcd
      ! Move the first row I(ia,0)
      do ia=0,la+lb+max(lla,llb)
        ja = ia
        jb = 0
        if (ja > la+2) then
          ja = ja-(la+3)
          jb = 1
        end if
        Arr2(:,ja,jb,icd,iCar) = Arr1(:,iCar,ia,icd)
      end do
      ! Generate I(ia,ib) for fixed ib
      do ib=1,lb+llb
        do ia=la+lb+max(lla,llb)-ib,0,-1
          ja = ia
          jb = ib
          mb = ib-1
          if (ja > la+2) then
            ja = ja-(la+3)
            jb = jb+1
            mb = mb+1
          end if
          ma = ja
          ka = ia+1
          kb = ib-1
          if (ka > la+2) then
            ka = ka-(la+3)
            kb = kb+1
          end if
          Arr2(:,ja,jb,icd,iCar) = AB*Arr2(:,ma,mb,icd,iCar)+Arr2(:,ka,kb,icd,iCar)
        end do
      end do
    end do
  else
    AB = -AB
    do icd=0,lc+ld+llcd
      ! Move the first row I(0,ib)
      do ib=0,la+lb+max(lla,llb)
        jb = ib
        ja = 0
        if (jb > lb+2) then
          jb = jb-(lb+3)
          ja = 1
        end if
        Arr2(:,ja,jb,icd,iCar) = Arr1(:,iCar,ib,icd)
      end do
      ! Generate I(ia,ib) for fixed ia
      do ia=1,la+lla
        do ib=la+lb+max(lla,llb)-ia,0,-1
          jb = ib
          ja = ia
          ma = ia-1
          if (jb > lb+2) then
            jb = jb-(lb+3)
            ja = ja+1
            ma = ma+1
          end if
          mb = jb
          kb = ib+1
          ka = ia-1
          if (kb > lb+2) then
            kb = kb-(lb+3)
            ka = ka+1
          end if
          Arr2(:,ja,jb,icd,iCar) = AB*Arr2(:,ma,mb,icd,iCar)+Arr2(:,ka,kb,icd,iCar)
        end do
      end do
    end do
  end if
end do

return

end subroutine HRR2Da_mck

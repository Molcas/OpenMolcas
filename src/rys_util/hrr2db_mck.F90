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

subroutine Hrr2Db_mck(Arr1,nVec,ncdMax,Arr2,C,D,la,lb,lc,ld,IfHss,IfGrd,nt,nrys)
!***********************************************************************
!                                                                      *
! Object: to apply the transfer equation to the 2D-integrals.          *
!         The transformation is in place and the recursion             *
!         is replaced with the indentity when applicable.              *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             September '91                                            *
!             Modified to recurrence algorithm, February '92           *
!             Improved algorithm, June '92.                            *
!***********************************************************************

implicit real*8(A-H,O-Z)
!#include "print.fh"
#include "real.fh"
real*8 C(3), D(3), Arr1(nVec,0:la+2,0:lb+2,0:ncdMax,3), Arr2(nVec,0:la+2,0:lb+2,0:lc+2,0:ld+2,3)
logical IfHss(4,3,4,3), ifgrd(3,4)

!iRout = 233
!iPrint = nPrint(iRout)

do iCar=1,3
  llc = 0
  lld = 0
  lla = 0
  llb = 0
  if (IfGrd(iCar,3)) llc = max(llc,1)
  if (IfHss(3,iCar,3,iCar)) llc = 2
  if (IfGrd(iCar,4)) lld = max(1,lld)
  if (IfHss(4,iCar,4,iCar)) lld = 2
  if (IfGrd(iCar,1)) lla = max(1,lla)
  if (IfHss(1,iCar,1,iCar)) lla = 2
  if (IfGrd(iCar,2)) llb = max(llb,1)
  if (IfHss(2,iCar,2,iCar)) llb = 2
  CD = C(iCar)-D(iCar)
  if (CD == Zero) then
    do ia=0,la+lla
      do ib=0,lb+llb
        if (ia+ib > la+lb+max(lla,llb)) Go To 101
        ! Using the identity
        do ic=0,lc+llc
          do id=0,ld+lld
            icd = ic+id
            if (icd > lc+ld+max(llc,lld)) Go To 210
            do i=1,nVec
              Arr2(i,ia,ib,ic,id,iCar) = Arr1(i,ia,ib,icd,iCar)
            end do
210         continue
          end do
        end do
101     continue
      end do
    end do
  else
    if (lc >= ld) then
      do ia=0,la+lla
        do ib=0,lb+llb
          if (ia+ib > la+lb+max(lla,llb)) Go To 103
          ! Move the first row I(ic,0)
          do ic=0,lc+ld+max(llc,lld)
            jc = ic
            jd = 0
            if (jc > lc+2) then
              jc = jc-(lc+3)
              jd = 1
            end if
            do i=1,nVec
              Arr2(i,ia,ib,jc,jd,iCar) = Arr1(i,ia,ib,ic,iCar)
            end do
          end do
          ! Generate I(ic,id) for fixed id
          do id=1,ld+lld
            do ic=lc+ld+max(llc,lld)-id,0,-1
              jc = ic
              jd = id
              md = id-1
              if (jc > lc+2) then
                jc = jc-(lc+3)
                jd = jd+1
                md = md+1
              end if
              mc = jc
              kc = ic+1
              kd = id-1
              if (kc > lc+2) then
                kc = kc-(lc+3)
                kd = kd+1
              end if
              call DZaXpY(nVec,CD,Arr2(1,ia,ib,mc,md,iCar),1,Arr2(1,ia,ib,kc,kd,iCar),1,Arr2(1,ia,ib,jc,jd,iCar),1)
            end do
          end do
103       continue
        end do
      end do
    else
      CD = -CD
      do ia=0,la+lla
        do ib=0,lb+llb
          if (ia+ib > la+lb+max(lla,llb)) Go To 105
          ! Move the first row I(0,id)
          do id=0,lc+ld+max(llc,lld)
            jd = id
            jc = 0
            if (jd > ld+2) then
              jd = jd-(ld+3)
              jc = 1
            end if
            do i=1,nVec
              Arr2(i,ia,ib,jc,jd,iCar) = Arr1(i,ia,ib,id,iCar)
            end do
          end do
          ! Generate I(ic,id) for fixed ic
          do ic=1,lc+llc
            do id=lc+ld+max(llc,lld)-ic,0,-1
              jd = id
              jc = ic
              mc = ic-1
              if (jd > ld+2) then
                jd = jd-(ld+3)
                jc = jc+1
                mc = mc+1
              end if
              md = jd
              kd = id+1
              kc = ic-1
              if (kd > ld+2) then
                kd = kd-(ld+3)
                kc = kc+1
              end if
              call DZaXpY(nVec,CD,Arr2(1,ia,ib,mc,md,iCar),1,Arr2(1,ia,ib,kc,kd,iCar),1,Arr2(1,ia,ib,jc,jd,iCar),1)
            end do
          end do
105       continue
        end do
      end do
    end if
  end if
end do

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(nt)
  call Unused_integer(nrys)
end if

end subroutine Hrr2Db_mck

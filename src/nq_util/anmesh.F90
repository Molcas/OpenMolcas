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
! Copyright (C) 2001, Roland Lindh                                     *
!               2001, Laura Gagliardi                                  *
!***********************************************************************

subroutine AnMesh(nscheme,pa,rPt,wPt)

use Constants, only: Zero, One, Two, Three
use Definitions, only: wp, iwp
!define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nscheme(8)
real(kind=wp), intent(in) :: pa(*)
real(kind=wp), intent(_OUT_) :: rPt(3,*), wPt(*)
integer(kind=iwp) :: i, ii, ip, ix, iy, iz, j, j1, jj, n1
real(kind=wp) :: c, pp, qq, rr, ss, tt, uu, vv

!                                                                      *
!***********************************************************************
!                                                                      *

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' ******** The Angular Lebedev Grid ********'
write(u6,*)
write(u6,*)
#endif

i = 0
ip = 0

! nscheme(2) -> 6 points

!lg write(u6,*) 'nscheme',(nscheme(i),i=1,8)
if (nscheme(2) > 0) then
  !lg write(u6,*) 'nscheme(2)',nscheme(2)
  ip = ip+1
  do ix=1,3
    do iy=1,-1,-2
      i = i+1
      wPt(i) = pa(ip)
      do j=1,3
        rPt(j,i) = Zero
      end do
      rPt(ix,i) = real(iy,kind=wp)
      !lg write(u6,*) rPt(ix,i),wPt(i)
    end do
  end do
end if

! nscheme(3) -> 8 points

if (nscheme(3) > 0) then   !
  c = One/sqrt(Three)
  ip = ip+1
  do ix=1,-1,-2
    do iy=1,-1,-2
      do iz=1,-1,-2
        i = i+1
        wPt(i) = pa(ip)
        rPt(1,i) = real(ix,kind=wp)*c
        rPt(2,i) = real(iy,kind=wp)*c
        rPt(3,i) = real(iz,kind=wp)*c
      end do
    end do
  end do
end if

! nscheme(4) -> 12 points

if (nscheme(4) > 0) then
  c = One/sqrt(Two)
  ip = ip+1
  do ix=1,-1,-2
    do iy=1,-1,-2
      do iz=1,3
        i = i+1
        wPt(i) = pa(ip)
        rPt(iz,i) = real(ix,kind=wp)*c
        j = mod(iz,3)+1
        rPt(j,i) = real(iy,kind=wp)*c
        j = 6-iz-j
        rPt(j,i) = Zero
      end do
    end do
  end do
end if

! 24a points

n1 = nscheme(5)
do jj=1,n1
  ip = ip+1
  uu = pa(ip)
  vv = sqrt(One-Two*uu*uu)
  ip = ip+1
  do ix=1,-1,-2
    do iy=1,-1,-2
      do iz=1,-1,-2
        do j=1,3
          i = i+1
          wPt(i) = pa(ip)
          do j1=1,3
            rPt(j1,i) = uu
          end do
          rPt(j,i) = vv
          rPt(1,i) = rPt(1,i)*real(ix,kind=wp)
          rPt(2,i) = rPt(2,i)*real(iy,kind=wp)
          rPt(3,i) = rPt(3,i)*real(iz,kind=wp)
        end do
      end do
    end do
  end do
end do

! 24b points

n1 = nscheme(6)
do jj=1,n1
  ip = ip+1
  pp = pa(ip)
  qq = sqrt(One-pp*pp)
  ip = ip+1
  do ix=1,-1,-2
    do iy=1,-1,-2
      do ii=0,1
        do j=1,3
          i = i+1
          wPt(i) = pa(ip)
          j1 = mod(j+ii,3)+1
          rPt(j1,i) = pp*real(ix,kind=wp)
          j1 = mod(j+1-ii,3)+1
          rPt(j1,i) = qq*real(iy,kind=wp)
          rPt(j,i) = Zero
        end do
      end do
    end do
  end do
end do

! 48 points

n1 = nscheme(7)
!lg write(u6,*) 'i, n1 =',i,n1
do jj=1,n1
  ip = ip+1
  rr = pa(ip)
  ip = ip+1
  ss = pa(ip)
  tt = sqrt(One-rr*rr-ss*ss)
  ip = ip+1
  do ix=1,-1,-2
    do iy=1,-1,-2
      do iz=1,-1,-2
        do j=1,3
          do ii=0,1
            i = i+1
            wPt(i) = pa(ip)
            rPt(j,i) = rr*real(ix,kind=wp)
            j1 = mod(j+ii,3)+1
            rPt(j1,i) = ss*real(iy,kind=wp)
            j1 = mod(j+1-ii,3)+1
            rPt(j1,i) = tt*real(iz,kind=wp)
            !lg write (u6,*) rPt(j1,i),wPt(i),j1,i
          end do
          !lg write (u6,*) 'Enddo1',i
        end do
        !lg write (u6,*) 'Enddo2',i
      end do
      !lg write (u6,*) 'Enddo3',i
    end do
    !lg write (u6,*) 'Enddo4',i
  end do
  !lg write (u6,*) 'Enddo5',i
end do
!lg write (u6,*) 'enddo',n1
!                                                                      *
!***********************************************************************
!                                                                      *

!lg write(u6,*) 'End of AnMesh'

return

end subroutine AnMesh

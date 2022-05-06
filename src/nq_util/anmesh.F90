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

use nq_Info
implicit real*8(a-h,o-z)
implicit integer(i-n)
#include "real.fh"
#include "debug.fh"
#include "WrkSpc.fh"
integer nscheme(8)
real*8 pa(*), rPt(3,*), wPt(*)

!                                                                      *
!***********************************************************************
!                                                                      *

if (Debug) then
  write(6,*)
  write(6,*) ' ******** The Angular Lebedev Grid ********'
  write(6,*)
  write(6,*)
end if

i = 0
ip = 0

! nscheme(2) -> 6 points

!lg write(6,*) 'nscheme',(nscheme(i),i=1,8)
if (nscheme(2) > 0) then
  !lg write(6,*) 'nscheme(2)',nscheme(2)
  ip = ip+1
  do ix=1,3
    do iy=1,-1,-2
      i = i+1
      wPt(i) = pa(ip)
      do j=1,3
        rPt(j,i) = Zero
      end do
      rPt(ix,i) = dble(iy)
      !lg write(6,*) rPt(ix,i),wPt(i)
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
        rPt(1,i) = dble(ix)*c
        rPt(2,i) = dble(iy)*c
        rPt(3,i) = dble(iz)*c
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
        rPt(iz,i) = dble(ix)*c
        j = mod(iz,3)+1
        rPt(j,i) = dble(iy)*c
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
          rPt(1,i) = rPt(1,i)*dble(ix)
          rPt(2,i) = rPt(2,i)*dble(iy)
          rPt(3,i) = rPt(3,i)*dble(iz)
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
          rPt(j1,i) = pp*dble(ix)
          j1 = mod(j+1-ii,3)+1
          rPt(j1,i) = qq*dble(iy)
          rPt(j,i) = zero
        end do
      end do
    end do
  end do
end do

! 48 points

n1 = nscheme(7)
!lg write(6,*) 'i, n1 =',i,n1
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
            rPt(j,i) = rr*dble(ix)
            j1 = mod(j+ii,3)+1
            rPt(j1,i) = ss*dble(iy)
            j1 = mod(j+1-ii,3)+1
            rPt(j1,i) = tt*dble(iz)
            !lg write (6,*) rPt(j1,i),wPt(i),j1,i
          end do
          !lg write (6,*) 'Enddo1',i
        end do
        !lg write (6,*) 'Enddo2',i
      end do
      !lg write (6,*) 'Enddo3',i
    end do
    !lg write (6,*) 'Enddo4',i
  end do
  !lg write (6,*) 'Enddo5',i
end do
!lg write (6,*) 'enddo',n1
!                                                                      *
!***********************************************************************
!                                                                      *

!lg write(6,*) 'End of AnMesh'

return

end subroutine AnMesh

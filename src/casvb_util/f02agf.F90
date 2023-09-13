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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine f02agf(a,ia,n,rr,ri,vr,ivr,vi,ivi,intger,ifail)

implicit real*8(a-h,o-z)
logical pair
dimension a(ia,*), rr(*), ri(*), vr(ivr,*), vi(ivi,*), intger(*)
save thresh
data thresh/1.d-8/

if (ifail /= 0) call SysHalt('ifail f02agf')

if ((ia /= ivr) .or. (ia /= ivi)) call SysHalt('f02agf dim')
call rg(ia,n,a,rr,ri,1,vr,intger,vi,info)
if (info /= 0) call SysHalt('info f02agf')
call fzero(vi,n*ivi)
pair = .false.
do k=1,n-1
  if ((ri(k) /= 0d0) .and. (.not. pair)) then
    if (rr(k) /= rr(k+1)) call SysHalt('rr trouble')
    if (abs(ri(k)+ri(k+1)) > 1d-12) call SysHalt('ri trouble')
    pair = .true.
    ! If eig value almost real: return real value & vectors:
    if (abs(ri(k)) <= thresh) then
      ri(k) = 0.d0
      ri(k+1) = 0.d0
    else
      do l=1,n
        vi(l,k) = vr(l,k+1)
        vi(l,k+1) = -vr(l,k+1)
      end do
      do l=1,n
        vr(l,k+1) = vr(l,k)
      end do
    end if
  else
    pair = .false.
  end if
end do

end subroutine f02agf

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

use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ia, n, ivr, ivi, ifail
real(kind=wp), intent(inout) :: a(ia,*)
real(kind=wp), intent(_OUT_) :: rr(*), ri(*), vr(ivr,*), vi(ivi,*)
integer(kind=iwp), intent(_OUT_) :: intger(*)
integer(kind=iwp) :: info, k
logical(kind=iwp) :: pair
real(kind=wp), parameter :: thresh = 1.0e-8_wp

if (ifail /= 0) call SysHalt('ifail f02agf')

if ((ia /= ivr) .or. (ia /= ivi)) call SysHalt('f02agf dim')
call rg(ia,n,a,rr,ri,1,vr,intger,vi,info)
if (info /= 0) call SysHalt('info f02agf')
vi(:,1:n) = Zero
pair = .false.
do k=1,n-1
  if ((ri(k) /= Zero) .and. (.not. pair)) then
    if (rr(k) /= rr(k+1)) call SysHalt('rr trouble')
    if (abs(ri(k)+ri(k+1)) > 1.0e-12_wp) call SysHalt('ri trouble')
    pair = .true.
    ! If eig value almost real: return real value & vectors:
    if (abs(ri(k)) <= thresh) then
      ri(k) = Zero
      ri(k+1) = Zero
    else
      vi(1:n,k) = vr(1:n,k+1)
      vi(1:n,k+1) = -vr(1:n,k+1)
      vr(1:n,k+1) = vr(1:n,k)
    end if
  else
    pair = .false.
  end if
end do

end subroutine f02agf

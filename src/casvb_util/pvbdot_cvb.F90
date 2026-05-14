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

subroutine pvbdot_cvb(cfrom,cto,ret)

use casvb_global, only: iapr, iform_ci, ixapr, ndet
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: cfrom(0:ndet)
real(kind=wp), intent(inout) :: cto(0:ndet)
real(kind=wp), intent(out) :: ret
integer(kind=iwp) :: icfrom, icto

icfrom = nint(cfrom(0))
icto = nint(cto(0))
if ((iform_ci(icfrom) /= 0) .or. (iform_ci(icto) /= 0)) then
  write(u6,*) ' Unsupported format in PVBDOT'
  call abend_cvb()
end if
call pvbcopy2_cvb(cfrom(1:),cto(1:),iapr,ixapr,ret,1)

return

end subroutine pvbdot_cvb

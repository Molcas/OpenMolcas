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

subroutine pvbcopy_cvb(cfrom,cto)

use casvb_global, only: iapr, icnt_ci, iform_ci, ixapr, ndet
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: cfrom(0:ndet)
real(kind=wp), intent(inout) :: cto(0:ndet)
integer(kind=iwp) :: icfrom, icto
real(kind=wp) :: dum

icfrom = nint(cfrom(0))
icto = nint(cto(0))
if ((iform_ci(icfrom) /= 0) .or. (iform_ci(icto) /= 0)) then
  write(u6,*) ' Unsupported format in PVBCOPY'
  call abend_cvb()
end if
call pvbcopy2_cvb(cfrom(1:),cto(1:),iapr,ixapr,dum,0)
icnt_ci(icto) = 0

return

end subroutine pvbcopy_cvb

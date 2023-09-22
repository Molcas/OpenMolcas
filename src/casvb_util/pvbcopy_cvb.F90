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

use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: cfrom(*), cto(*)
#include "main_cvb.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: icfrom, icto
real(kind=wp) :: dum

icfrom = nint(cfrom(1))
icto = nint(cto(1))
if ((iform_ci(icfrom) /= 0) .or. (iform_ci(icto) /= 0)) then
  write(u6,*) ' Unsupported format in PVBCOPY'
  call abend_cvb()
end if
call pvbcopy2_cvb(work(iaddr_ci(icfrom)),work(iaddr_ci(icto)),iwork(ll(11)),iwork(ll(12)),dum,0)
call setcnt2_cvb(icto,0)

return

end subroutine pvbcopy_cvb

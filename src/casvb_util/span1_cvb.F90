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

subroutine span1_cvb(c,nvec,s,n,metr)

use casvb_global, only: iaddr, nvecmx, nvtot
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nvec, n, metr
real(kind=wp) :: c(n,nvec), s(*)
#include "WrkSpc.fh"
integer(kind=iwp) :: nvmove, nvremain

nvremain = nvec
do
  nvmove = min(nvremain,nvecmx-nvtot)
  if ((nvmove == 0) .and. (nvremain > 0)) then
    write(u6,*) ' Fatal error in SPAN_CVB!',nvmove,nvremain
    call abend_cvb()
  end if
  call fmove_cvb(c(1,1+nvec-nvremain),work(nvtot*n+iaddr),n*nvmove)
  nvtot = nvtot+nvmove
  if (nvtot == nvecmx) call span_cvb(work(iaddr),nvtot,nvtot,s,n,metr)
  nvremain = nvremain-nvmove
  if (nvremain <= 0) exit
end do

return

end subroutine span1_cvb

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

use casvb_global, only: nvecmx, nvtot, span
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nvec, n, metr
real(kind=wp), intent(in) :: c(n,nvec), s(*)
integer(kind=iwp) :: nvmove, nvremain

nvremain = nvec
do
  nvmove = min(nvremain,nvecmx-nvtot)
  if ((nvmove == 0) .and. (nvremain > 0)) then
    write(u6,*) ' Fatal error in SPAN_CVB!',nvmove,nvremain
    call abend_cvb()
  end if
  span(:,nvtot+1:nvtot+nvmove) = c(:,nvec-nvremain+1:nvec-nvremain+nvmove)
  nvtot = nvtot+nvmove
  if (nvtot == nvecmx) call span_cvb(span,nvtot,nvtot,s,n,metr)
  nvremain = nvremain-nvmove
  if (nvremain <= 0) exit
end do

return

end subroutine span1_cvb

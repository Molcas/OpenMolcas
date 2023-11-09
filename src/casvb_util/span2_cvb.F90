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

subroutine span2_cvb(c,nvec,s,n,metr)

use casvb_global, only: nvtot, span
use stdalloc, only: mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(inout) :: nvec
integer(kind=iwp), intent(in) :: n, metr
real(kind=wp), intent(inout) :: c(n,nvec)
real(kind=wp), intent(in) :: s(*)
integer(kind=iwp) :: nvtot_

if (nvtot /= 0) then
  call span_cvb(span,nvtot,nvtot_,s,n,metr)
  nvtot = nvtot_
  c(:,1:nvtot) = span(:,1:nvtot)
end if
nvec = nvtot
call mma_deallocate(span)

return

end subroutine span2_cvb

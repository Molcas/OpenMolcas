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

subroutine span0_cvb(nvecmx1,n)

use casvb_global, only: nvecmx, nvtot, span
use stdalloc, only: mma_allocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nvecmx1, n
integer(kind=iwp) :: mavailr
integer(kind=iwp), parameter :: nmult = 5

call mma_maxDBLE(mavailr)
nvecmx = min(nmult*nvecmx1,mavailr/n)
if (nvecmx <= 0) then
  write(u6,*) ' Not enough vectors in SPAN0_CVB!',nvecmx
  write(u6,*) ' Remaining memory :',mavailr
  write(u6,*) ' Max number of vectors :',nvecmx1
  call abend_cvb()
end if
call mma_allocate(span,n,nvecmx,label='span')
nvtot = 0

return

end subroutine span0_cvb

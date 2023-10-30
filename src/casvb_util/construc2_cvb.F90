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

subroutine construc2_cvb(tconstr)

use casvb_global, only: iconstruc, nconstr, nvb
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: tconstr(nvb,nvb)
integer(kind=iwp) :: iconstruc_kp, ivb
real(kind=wp) :: dum(1)
real(kind=wp), allocatable :: repm(:)

iconstruc_kp = iconstruc
iconstruc = 1
call mma_allocate(repm,nvb,label='nvb')

call span0_cvb(nvb,nvb)
do ivb=1,nvb
  repm(:) = Zero
  repm(ivb) = -One
  call symtrizcvb_cvb(repm)
  repm(ivb) = repm(ivb)+One
  call span1_cvb(repm,1,dum,nvb,0)
end do
call span2_cvb(tconstr,nconstr,dum,nvb,0)

call mma_deallocate(repm)
iconstruc = iconstruc_kp

return

end subroutine construc2_cvb

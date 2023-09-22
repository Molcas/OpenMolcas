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

use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
real(kind=wp) :: tconstr(nvb,nvb)
#include "WrkSpc.fh"
integer(kind=iwp) :: iconstruc_kp, irepm, ivb
real(kind=wp) :: dum(1)
integer(kind=iwp), external :: mstackr_cvb

iconstruc_kp = iconstruc
iconstruc = 1
irepm = mstackr_cvb(nvb)

call span0_cvb(nvb,nvb)
do ivb=1,nvb
  call fzero(work(irepm),nvb)
  work(ivb+irepm-1) = -One
  call symtrizcvb_cvb(work(irepm))
  work(ivb+irepm-1) = work(ivb+irepm-1)+One
  call span1_cvb(work(irepm),1,dum,nvb,0)
end do
call span2_cvb(tconstr,nconstr,dum,nvb,0)

call mfreer_cvb(irepm)
iconstruc = iconstruc_kp

return

end subroutine construc2_cvb

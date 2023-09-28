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

subroutine o12ea_cvb( &
#                    define _CALLING_
#                    include "opta_interface.fh"
                    )

use casvb_global, only: have_solved_it
use Definitions, only: iwp

implicit none
#include "opta_interface.fh"
#include "main_cvb.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: i1, i2, i3, i4, iv, ivuse2, ivuse_h, ivuse_s
integer(kind=iwp), external :: mstackr_cvb
logical(kind=iwp), external :: tstcnt_cvb ! ... Content of CI vectors ...

call ddnewopt_cvb()
have_solved_it = .false.

! Find CIVBS & CIVBH:
ivuse_s = 0
ivuse_h = 0
do iv=1,nv
  if (tstcnt_cvb(work(lc(iv)),4)) ivuse_s = iv
  if (tstcnt_cvb(work(lc(iv)),5)) ivuse_h = iv
end do
ivuse2 = 3
if ((ivuse_h == 3) .or. (ivuse_s == 3)) ivuse2 = 2
if ((ivuse_h == 2) .or. (ivuse_s == 2)) ivuse2 = 4
if (ivuse2 > nv) ivuse2 = 1
if ((ivuse_s /= 0) .and. (ivuse_h /= 0)) then
  i1 = mstackr_cvb(nparam)
  i2 = mstackr_cvb(nparam)
  i3 = mstackr_cvb(nparam)
  i4 = mstackr_cvb(nvb+nprorb)
  call o12ea2_cvb(work(i1),work(i2),work(i3),nparam,work(lc(ivuse2)),work(lc(ivuse_s)),work(lc(ivuse_h)),work(lw(9)),work(lv(2)), &
                  work(i4))
  call mfreer_cvb(i1)
else
  if (strucopt) then
    call ddguess_cvb(work(lv(2)),nvb,nprorb)
  else
    call ddguess_cvb([one],1,0)
  end if
end if
call str2vbc_cvb(work(lv(2)),work(lw(9)))
call vb2cic_cvb(work(lw(9)),work(lc(3)))

return

end subroutine o12ea_cvb

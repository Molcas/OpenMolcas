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

logical function firsttime_cvb()

use casvb_global, only: icode, inputmode, iopt2step, ioptim, ioptstep, istackrep, joptstep, loopstep

implicit real*8(a-h,o-z)
logical begbracket, second_time_round
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
external istkprobe_cvb
logical istkprobe_cvb

if (inputmode == 2) then
  if (loopstep == 0) then
    begbracket = .false.
  else
    begbracket = ((icode(loopstep) == 1) .and. (icode(loopstep+1) == 2)) .or. &
                 ((icode(loopstep) == 3) .and. (icode(loopstep+1) == 4))
  end if

  firsttime_cvb = ((joptstep == ioptstep-1) .or. ((ioptstep == 0) .and. (joptstep == 0))) .or. &
                  ((joptstep == ioptstep) .and. begbracket)

  if ((ioptim > 1) .and. (iopt2step(ioptim) == iopt2step(ioptim-1))) firsttime_cvb = .false.

  if (istkprobe_cvb(istackrep)) then
    call istkpop_cvb(istackrep,nc_zeroed_l)
    call istkpop_cvb(istackrep,nconvinone_l)
    call istkpop_cvb(istackrep,italter_l)
    second_time_round = (italter_l > 1)
    call istkpush_cvb(istackrep,italter_l)
    call istkpush_cvb(istackrep,nconvinone_l)
    call istkpush_cvb(istackrep,nc_zeroed_l)
  else
    second_time_round = .false.
  end if

  if (second_time_round) firsttime_cvb = .false.

  if (nmcscf >= 2) firsttime_cvb = .false.
else
  firsttime_cvb = .false.
end if

return

end function firsttime_cvb

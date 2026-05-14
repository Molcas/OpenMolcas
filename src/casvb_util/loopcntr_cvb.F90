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

subroutine loopcntr_cvb(icode1)

use casvb_global, only: icnt, icode, inputmode, ioptstep, ipos, joptstep, loopstep, loopstepmx, mxstep, ncnt, noptstep
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: icode1
logical(kind=iwp) :: begbracket

loopstep = loopstep+1
if ((inputmode == 2) .and. ((icode1 == 5) .or. (icode1 == 6))) return
if (loopstep > mxstep) then
  write(u6,*) ' Loop structure too complicated -- mxstep :',mxstep
  call abend_cvb()
end if
if (inputmode == 1) then
  icode(loopstep) = icode1
  ipos(loopstep) = ncnt
end if
if ((icode(loopstep) == 1) .or. (icode(loopstep) == 3)) joptstep = joptstep+1
if (inputmode == 2) then
  if (joptstep == ioptstep) call setifinish_cvb(icode(loopstep))
  begbracket = ((icode(loopstep) == 1) .and. (icode(loopstep+1) == 2)) .or. ((icode(loopstep) == 3) .and. (icode(loopstep+1) == 4))
  if ((joptstep >= ioptstep+1) .or. ((joptstep == ioptstep) .and. (.not. begbracket))) then
    ! Go to end of input:
    icnt = ncnt
    loopstep = loopstepmx
    joptstep = noptstep
  else if ((joptstep < ioptstep) .and. begbracket) then
    ! Go to next closing bracket:
    icnt = ipos(loopstep+1)
    loopstep = loopstep+1
  end if
end if

return

end subroutine loopcntr_cvb

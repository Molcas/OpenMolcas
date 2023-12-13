!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine MandatoryInp(YesNo)

use Definitions, only: iwp, u6

implicit none
logical(kind=iwp), intent(in) :: YesNo(*)
#include "warnings.h"

if ((.not. YesNo(8)) .and. (.not. YesNo(7))) then
  write(u6,*)
  write(u6,*) ' You have not specified what type of calculation this is.'
  write(u6,*) ' Use either the RUN keyword or the SINGle-point keyword.'
  call Quit(_RC_INPUT_ERROR_)
end if

if (YesNo(3) .and. YesNo(4)) then
  write(u6,*)
  write(u6,*) ' You have specified both a SCFSection and a RASSisection.'
  write(u6,*) ' They are mutually exclusive. Remove one.'
  call Quit(_RC_INPUT_ERROR_)
end if

if (YesNo(7) .and. (.not. YesNo(6))) then
  write(u6,*)
  write(u6,*) ' You have requested a single-point calculation, but no input coordinates were given.'
  write(u6,*) ' Provide these in the SOLVent section.'
  call Quit(_RC_INPUT_ERROR_)
end if

if (YesNo(5) .and. (.not. YesNo(6))) then
  write(u6,*)
  write(u6,*) ' You have specified that initial coordinates are to be given in input, but no coordinates are found.'
  write(u6,*) ' Provide these in the SOLVent section.'
  call Quit(_RC_INPUT_ERROR_)
end if

if ((.not. YesNo(2)) .and. (.not. YesNo(7))) then
  write(u6,*)
  write(u6,*) ' You fail to specify where from initial configuration should be collected.'
  write(u6,*) ' Do this with the CONFiguration keyword.'
  call Quit(_RC_INPUT_ERROR_)
end if

if (YesNo(9) .and. (.not. YesNo(10))) then
  write(u6,*)
  write(u6,*) ' Your file specification implies that an extraction file is to be generated.'
  write(u6,*) ' However, you have no EXTRact section.'
  call Quit(_RC_INPUT_ERROR_)
end if

if ((.not. YesNo(9)) .and. YesNo(10)) then
  write(u6,*)
  write(u6,*) ' You have a EXTRact section, but the file to read from is not a sampfile.'
  write(u6,*) ' Change this after the FILE keyword.'
  call Quit(_RC_INPUT_ERROR_)
end if

return

end subroutine MandatoryInp

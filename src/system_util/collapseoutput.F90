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
! Copyright (C) Valera Veryazov                                        *
!***********************************************************************
!  CollapseOutput
!
!> @brief
!>   Place markers for collapsed parts in output
!> @author V. Veryazov
!>
!> @details
!> If \p iSw = ``1``, print Message preceded by a marker to
!> start collapsed part in the output.
!>
!> If \p iSw = ``0``, place a marker to terminate collapsed part.
!> \p STR is not used, but recommended.
!>
!> @param[in] iSw Begin/End switch
!> @param[in] STR Message
!***********************************************************************

subroutine CollapseOutput(iSw,STR)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iSw
character(len=*), intent(in) :: STR
#include "print.fh"

if (icolorize == 1) then
  if (iSw == 1) then
    write(u6,'(A,A)') '++ ',trim(STR)
  else
    write(u6,'(A)') '--'
  end if
else
  if (iSw == 1) write(u6,'(A)') trim(STR)
end if

return

end subroutine CollapseOutput

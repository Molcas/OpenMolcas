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

subroutine Lu2Lu(Filename,LuInput)

use Definitions, only: iwp, u6

implicit none
character(len=*), intent(in) :: FileName
integer(kind=iwp), intent(in) :: LuInput
#include "warnings.h"
integer(kind=iwp) :: istatus, LuSpool2
character(len=180) :: Line
logical(kind=iwp) :: Exists
integer(kind=iwp), external :: IsFreeUnit

call f_inquire(Filename,Exists)
if (.not. Exists) then
  write(u6,*) 'SuperMac: Missing ',Filename
  call Finish(_RC_INTERNAL_ERROR_)
end if

LuSpool2 = 77
LuSpool2 = IsFreeUnit(LuSpool2)
call Molcas_Open(LuSpool2,Filename)

do
  read(LuSpool2,'(A)',iostat=istatus) Line
  if (istatus < 0) exit
  write(LuInput,'(A)') Line
end do

close(LuSpool2)

return

end subroutine Lu2Lu

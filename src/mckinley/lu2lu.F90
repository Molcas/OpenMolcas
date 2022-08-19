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

character FileName*(*), Line*180
logical Exist
#include "warnings.h"

call f_inquire(Filename,Exist)
if (.not. Exist) then
  write(6,*) 'SuperMac: Missing ',Filename
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

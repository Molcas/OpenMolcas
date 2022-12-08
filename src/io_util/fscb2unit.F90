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
! Copyright (C) Victor P. Vysotskiy                                    *
!***********************************************************************
!  FSCB2UNIT
!
!> @brief
!>   Translate system (C)file descriptor into internal Molcas's one
!> @author V. Vysotskiy
!>
!> @details
!> Translate system (C)file descriptor into internal Molcas's one
!>
!> @param[in]  cunit System (C)file descriptor
!> @param[out] LuP
!***********************************************************************

subroutine FSCB2UNIT(cunit,LuP)

use Fast_IO, only: FSCB, LuName, LuNameProf, MxFile, NProfFiles
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: cunit
integer(kind=iwp), intent(out) :: LuP
integer(kind=iwp) :: i, Lu

Lu = -1
do i=1,MxFile
  !write(u6,*) i,FSCB(i),cunit
  if (FSCB(i) == cunit) then
    Lu = i
  end if
end do
#ifndef _I8_
Lu = MPUnit(0,Lu)
#endif
LuP = -1
if (Lu == -1) call Abend()
do i=1,NProfFiles
  !write(u6,*) i,LuNameProf(i),LuName(Lu)
  if (LuNameProf(i) == LuName(Lu)) then
    LuP = i
  end if
end do

if (LuP == -1) call Abend()

return

end subroutine FSCB2UNIT

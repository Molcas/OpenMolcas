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
! Copyright (C) 2001, Valera Veryazov                                  *
!***********************************************************************

subroutine SysValueWarnMsg(Text1,N1)

use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: Text1
integer(kind=iwp), intent(in) :: N1
character(len=20) :: s

write(s,'(a,i16)') ' = ',N1
call SysPuts('Value: ',Text1,s)

return

end subroutine SysValueWarnMsg

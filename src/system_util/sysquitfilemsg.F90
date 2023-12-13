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

subroutine SysQuitFileMsg(rc,Location,TheFile,Text1,Text2)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: rc
character(len=*), intent(in) :: Location, Text1, Text2, TheFile

call SysWarnFileMsg(Location,TheFile,Text1,Text2)
call Quit(rc)

return

end subroutine SysQuitFileMsg

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
!***********************************************************************
!  SysAbendMsg
!
!> @brief
!>   Stop calculation
!> @author V. Veryazov
!>
!> @details
!> Print formatted message and stop the calculation.
!>
!> @param[in] Location routine name
!> @param[in] Text1    message text
!> @param[in] Text2    message text
!***********************************************************************

subroutine SysAbendMsg(Location,Text1,Text2)

implicit none
character(len=*), intent(in) :: Location, Text1, Text2

call SysWarnMsg(Location,Text1,Text2)
call Abend()

return

end subroutine SysAbendMsg

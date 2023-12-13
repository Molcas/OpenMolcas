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
!  StatusLine
!
!> @brief
!>   Update status file
!> @author V. Veryazov
!>
!> @details
!> Print parameters into status file. There is no
!> difference between params.
!>
!> \code
!> call StatusLine(ModuleName,': just started')
!> \endcode
!>
!> @param[in] STR  Status
!> @param[in] STR1 Status
!***********************************************************************

subroutine StatusLine(STR,STR1)

#ifdef _MOLCAS_MPP_
use Para_Info, only: King
#endif
use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: STR, STR1
integer(kind=iwp) :: Lu

#ifdef _MOLCAS_MPP_
if (.not. King()) return
#endif
Lu = 2
call molcas_open(Lu,'status')
write(Lu,'(A,A)') STR,STR1
close(Lu)

return

end subroutine StatusLine

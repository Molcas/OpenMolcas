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

subroutine CHO_QUIT(STRING,IERR)
!
! Purpose: echo message STRING and abort execution.

use Cholesky, only: LuPri
use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: STRING
integer(kind=iwp), intent(in) :: IERR
integer(kind=iwp) :: MOLCASCODE

write(LUPRI,'(//,1X,A)') '***'
if ((IERR == 0) .or. (IERR == 100)) then
  write(LUPRI,'(1X,A)') '*** Execution stopped by Cholesky Decomposition Utility'
  write(LUPRI,'(1X,A,A)') '*** Message: ',STRING
else
  write(LUPRI,'(1X,A)') '*** Error in Cholesky Core Routine'
  write(LUPRI,'(1X,A,A)') '*** Message: ',STRING
  write(LUPRI,'(1X,A,I5)') '*** Code   : ',IERR
end if
write(LUPRI,'(1X,A,//)') '***'
call CHO_TRANSLATEERRORCODE(IERR,MOLCASCODE)
call QUIT(MOLCASCODE)

end subroutine CHO_QUIT

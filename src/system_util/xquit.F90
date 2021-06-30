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

subroutine xquit(rc)
!SVC: routine that terminates Molcas properly

use warnings, only: rc_msg
#ifdef _MOLCAS_MPP_
use Para_Info, only: King
#endif
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: rc
integer(kind=iwp) :: lb, ub
character(len=128) :: msg
logical(kind=iwp), external :: bomb_on_error
#include "warnings.h"

call xflush(u6)

!SVC: write rc to stderr if not 0 (all is well)
if (rc /= _RC_ALL_IS_WELL_) then
  !IFG: do not write message if rc is "out of bounds"
  !     (this avoids it, e.g., when a slave process quits with -2,
  !     also avoids writing garbage from rc_msg)
  lb = lbound(rc_msg,dim=1)
  ub = ubound(rc_msg,dim=1)
  if ((rc >= lb) .and. (rc <= ub)) then
    write(msg,'(a,i6,2a)') 'xquit (rc = ',rc,'): ',rc_msg(rc)
    call write_stderr(msg)
  end if
end if

!SVC: write return code to file
call write_rc(rc)

!SVC: critical errors result in backtrace + immediate abort, while
!     regular errors only do this if MOLCAS_BOMB has been set too.
if ((rc >= _RC_GROUP_CRITICAL_) .or. &
#   if _MOLCAS_MPP_
!IFG: in a parallel (real or fake) run, we have to be more strict
!     with errors, or a deadlock may occur when some slave process
!     quits (e.g., a glitch or bug causes a missing file). Of course,
!     it could be argued that the error raised for those cases should
!     be critical... but it's often an "INPUT ERROR"
    (rc >= _RC_GROUP_USER_ERROR_ .and. (.not. King())) .or. &
#   endif
    (rc >= _RC_GROUP_ERROR_ .and. bomb_on_error())) then
  call xabort(rc)
end if

!SVC: terminate GA/MPI gracefully
call GATerminate()

!SVC: eventual exit normally. We should not exit with rc here because
!     some MPI implementations (like MPICH) will treat a non-zero return
!     code similar to an abort.
stop

end subroutine xquit

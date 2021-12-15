************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      subroutine xquit(rc)
CSVC: routine that terminates Molcas properly
      implicit none
      integer rc, lb, ub
#include "warnings.fh"
      character(len=128) :: msg
      logical, external :: bomb_on_error
#ifdef _MOLCAS_MPP_
      logical, external :: King
#endif

      call xflush(6)

CSVC: write rc to stderr if not 0 (all is well)
      if (rc .ne. _RC_ALL_IS_WELL_) then
CIFG: do not write message if rc is "out of bounds"
C     (this avoids it, e.g., when a slave process quits with -2,
C     also avoids writing garbage from rc_msg)
        lb = lbound(rc_msg,dim=1)
        ub = ubound(rc_msg,dim=1)
        if (rc .ge. lb .and. rc .le. ub) then
          write(msg,'(a,i6,2a)') 'xquit (rc = ', rc, '): ', rc_msg(rc)
          call write_stderr(msg)
        end if
      end if

CSVC: write return code to file
      call write_rc(rc)

CSVC: critical errors result in backtrace + immediate abort, while
C     regular errors only do this if MOLCAS_BOMB has been set too.
      if (        (rc .ge. _RC_GROUP_CRITICAL_)
#if _MOLCAS_MPP_
CIFG: in a parallel (real or fake) run, we have to be more strict
C     with errors, or a deadlock may occur when some slave process
C     quits (e.g., a glitch or bug causes a missing file). Of course,
C     it could be argued that the error raised for those cases should
C     be critical... but it's often an "INPUT ERROR"
     &        .or.(rc .ge. _RC_GROUP_USER_ERROR_ .and. (.not.King()))
#endif
     &        .or.(rc .ge. _RC_GROUP_ERROR_ .and. bomb_on_error())) then
        call xabort(rc)
      end if

CSVC: terminate GA/MPI gracefully
      call GATerminate

CSVC: eventual exit normally. We should not exit with rc here because
C     some MPI implementations (like MPICH) will treat a non-zero return
C     code similar to an abort.
      stop
      end

      function bomb_on_error() result(rc)
      implicit none
      logical :: rc
      character(len=16) :: bomb, env
      bomb=' '
      env='MOLCAS_BOMB'
      call getenvf(env,bomb)
      if (bomb(1:1).eq.'Y'.or.bomb(1:1).eq.'y'.or.bomb(1:1).eq.'1') then
        rc = .true.
      else
        rc = .false.
      end if
      end
c-----------------------------------------------------------------------
c     convenient wrapper routines for xquit
c-----------------------------------------------------------------------
      subroutine quit(rc)
      implicit none
      integer rc
      call xquit(rc)
      End
      Subroutine Quit_OnUserError
      implicit none
#include "warnings.fh"
      Call xQuit(_RC_INPUT_ERROR_)
      End
      Subroutine Quit_OnConvError
      implicit none
#include "warnings.fh"
      Call xQuit(_RC_NOT_CONVERGED_)
      End
      Subroutine Quit_OnInstError
      implicit none
#include "warnings.fh"
      Call xQuit(_RC_INSTALL_ERROR_)
      End
      Subroutine Quit_OnTimeOut
      implicit none
#include "warnings.fh"
      Call xQuit(_RC_TIMEOUT_)
      End

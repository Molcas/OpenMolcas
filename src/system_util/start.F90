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

subroutine Start(ModuleName)
! Initialization procedure for a program module.

use UnixInfo, only: init_UnixInfo
use warnings, only: rc_msg_init
use Para_Info, only: MyRank, King
#ifdef _MOLCAS_MPP_
use Para_Info, only: Set_Do_Parallel
use Definitions, only: iwp
#endif
#ifdef _DEBUGPRINT_
use Para_Info, only: nProcs
#endif
use Definitions, only: u5, u6

implicit none
character(len=*), intent(in) :: ModuleName
character(len=8) :: Prin
#ifdef _MOLCAS_MPP_
logical(kind=iwp) :: parallelized
#endif
#include "standard_iounits.fh"

#ifdef _DEBUGPRINT_
write(u6,*) ' Start.'
#endif

! fill the returncode database with human-readable messages
call rc_msg_init()

! Initialize delayed BLAS+LAPACK

call Init_LinAlg()
#ifdef _DEBUGPRINT_
write(u6,*) ' BLAS+LAPACK initialized '
#endif

! Initialize Timing

call SetTim()
#ifdef _DEBUGPRINT_
write(u6,*) ' TIMING set'
#endif

! Initialize the parallel environment

call GAInit()
#ifdef _DEBUGPRINT_
write(u6,*) ' GA initialized: MyRank, nProcs = ',MyRank,nProcs
#endif

! Now that each process is in its own directory, we can read/write
! files that are local to the' work directory of each process.

call write_rc(-1)
call getenvinit()

! Enable proper functions to handle various signals (needs environment!)

call set_sighandlers(myrank)

! Unix-related information must be set or checked here or later, PID and
! master/slave status may not have been set before now.

call write_pid()

! Have processes print their PID and then sleep for an amount of seconds
! specified by the environment variable MOLCAS_NAP. This is to allow for
! attaching a debugger to multiple processes.
call nap_time()

! Check that binary was executed via molcas script (check $MOLCAS_STAMP)
! (DO NOT MOVE FROM HERE)

#ifdef _HAVE_EXTRA_
call CheckRun()
#ifdef _DEBUGPRINT_
write(u6,*) ' CHECKRUN passed'
#endif
#endif

!  Initialize Memory manager (MOLCAS_MEM)
!  (DO NOT MOVE FROM HERE)
!
call IniMem()
#ifdef _DEBUGPRINT_
write(u6,*) ' MEMORY MANAGER initialized'
#endif
!
!  Get various unix-related information :
!  (DO NOT MOVE FROM HERE)
!
call init_UnixInfo(ModuleName,ModuleName)
#ifdef _DEBUGPRINT_
write(u6,*) ' UNIXINFO passed '
#endif

call prgminit(ModuleName)
#ifdef _DEBUGPRINT_
write(u6,*) ' PRGM initialized'
#endif

! Now that prgm has been initialized, we know if a module is supposed to
! support parallel execution. If not, then we can serialize the
! processes, after which the slaves jump to the finish and wait there.

#ifdef _MOLCAS_MPP_
parallelized = .true.
if (.not. parallelized) then
  call Set_Do_Parallel(.false.)
  call GATerminate()
  if (.not. King()) then
    call Finish(-2)
  end if
end if
#endif

! Redirect standard input to a file 'stdin'.

LuRd = u5
close(LuRd)
call molcas_open(LuRd,'stdin')
#ifdef _DEBUGPRINT_
write(u6,*) ' STDIN opened '
#endif

! Redirect standard output for slaves to a file 'stdout'.

LuWr = u6
if (.not. King()) then
  close(LuWr)
  call molcas_open(LuWr,'stdout')
  call Append_file(LuWr)
end if
#ifdef _DEBUGPRINT_
write(u6,*) ' STDOUT opened '
#endif

! Initialize XML

call ColorizeInit()
call xml_Open('module',' ',' ',0,ModuleName)
#ifdef _DEBUGPRINT_
write(u6,*) ' XML initialized '
#endif

! Initiate spool mode to read from standard input

Spool = .true.

! Initiate I/O

call FIOInit()
#ifdef _DEBUGPRINT_
write(u6,*) ' I/O initialized '
#endif

! Initialize timings and statistics

#ifdef _DEBUGPRINT_
write(u6,*) ' TIMINGS/STATISTICS initialized'
#endif

! Set standard name for the runfile

call NameRun('RUNFILE')
#ifdef _DEBUGPRINT_
write(u6,*) ' RUNFILE initialized '
#endif
call Init_Run_Use()

! Initialize the peek/poke utility

call Init_ppu(.true.)
call poke_iScalar('xml opened',0)

! Set Grid Set

call NQ_Init()
#ifdef _DEBUGPRINT_
write(u6,*) ' GRID initialized '
#endif

!***********************************************************************
! At this point, everything should have been properly initialized!     *
!***********************************************************************
! Finally, print useful runtime information about the module
call getenvf('MOLCAS_PRINT',Prin)
if ((Prin(1:1) /= '0') .and. (Prin(1:1) /= 'S')) then
  call Print_Module_Header(ModuleName)
  ! Force output to be written to stdout, such that in case of problems
  ! later on, at least all the start info has been printed.
  call xFlush(u6)
end if
! Write to the status file that the moduel has started
call StatusLine(ModuleName,' properly started!')

return

end subroutine Start

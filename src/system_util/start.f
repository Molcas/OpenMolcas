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
      Subroutine Start(ModuleName)
* Initialization procedure for a program module.
      Implicit None
      Character*(*) ModuleName
      Character*8 Prin
#include "para_info.fh"
#ifdef _MOLCAS_MPP_
      Logical parallelized
#endif
#include "standard_iounits.fh"

#ifdef _DEBUG_
      Write (6,*) ' Start.'
#endif

* fill the returncode database with human-readable messages
      Call rc_msg_init
* initialize global counter for the maximum warning level
      Call WarningInit

*
* Initialize delayed BLAS+LAPACK
*
      Call Init_LinAlg
#  ifdef _DEBUG_
      Write (6,*) ' BLAS+LAPACK initialized '
#  endif
*
*  Initialize Timing
*
      Call SetTim
#ifdef _DEBUG_
      Write (6,*) ' TIMING set'
#endif
*
* Initialize the parallel environment
*
      Call GAInit
#ifdef _DEBUG_
      Write(6,*) ' GA initialized: MyRank, nProcs = ', MyRank, nProcs
#endif
*
* Now that each process is in its own directory, we can read/write
* files that are local to the' work directory of each process.
*
      call write_rc(-1)
      Call getenvinit()
*
* Enable proper functions to handle various signals (needs environment!)
*
      call set_sighandlers(myrank)
*
*  Unix-related information must be set or checked here or later, PID and
*  master/slave status may not have been set before now.
*
      call write_pid
*
* Have processes print their PID and then sleep for an amount of seconds
* specified by the environment variable MOLCAS_NAP. This is to allow for
* attaching a debugger to multiple processes.
      call nap_time
*
*  Check that binary was executed via molcas script (check $MOLCAS_STAMP)
*  (DO NOT MOVE FROM HERE)
*
      Call CheckRun()
#  ifdef _DEBUG_
      Write (6,*) ' CHECKRUN passed'
#  endif

*  Initialize Memory manager (MOLCAS_MEM)
*  (DO NOT MOVE FROM HERE)
*
      Call IniMem()
#  ifdef _DEBUG_
      Write (6,*) ' MEMORY MANAGER initialized'
#  endif
*
*  Get various unix-related information :
*  (DO NOT MOVE FROM HERE)
*
      Call UnixInfo(ModuleName,ModuleName)
#  ifdef _DEBUG_
      Write (6,*) ' UNIXINFO passed '
#  endif

      call prgminit(ModuleName)
#  ifdef _DEBUG_
      Write (6,*) ' PRGM initialized'
#  endif

*  Now that prgm has been initialized, we know if a module is supposed to
*  support parallel execution. If not, then we can serialize the
*  processes, after which the slaves jump to the finish and wait there.

#ifdef _MOLCAS_MPP_
      parallelized = .True.
      If (.not.parallelized) Then
        Call Set_Do_Parallel(.False.)
        Call GATerminate
        If (.not.King()) Then
          Call Finish(-2)
        End If
      End If
#endif
*
* Redirect standard input to a file 'stdin'.
*
      LuRd=5
      close(LuRd)
      call molcas_open(LuRd,'stdin')
#  ifdef _DEBUG_
      Write (6,*) ' STDIN opened '
#  endif
*
*  Redirect standard output for slaves to a file 'stdout'.
*
      LuWr=6
      If(.Not.King())Then
        Close(LuWr)
        call molcas_open(LuWr,'stdout')
        call Append_file(LuWr)
      Endif
#  ifdef _DEBUG_
      Write (6,*) ' STDOUT opened '
#  endif
*
* Initialize XML
*
      Call ColorizeInit()
      call xml_Open('module',' ',' ',0,ModuleName)
#  ifdef _DEBUG_
      Write (6,*) ' XML initialized '
#  endif
*
*  Initiate spool mode to read from standard input
*
      Spool=.True.
*
*  Initiate I/O
*
      Call FIOInit()
#  ifdef _DEBUG_
      Write (6,*) ' I/O initialized '
#  endif
*
*  Initialize timings and statistics
*
      Call IniTim()
      Call IniStat()
#  ifdef _DEBUG_
      Write (6,*) ' TIMINGS/STATISTICS initialized'
#  endif
*
* Set standard name for the runfile
*
      Call NameRun('RUNFILE')
#  ifdef _DEBUG_
      Write (6,*) ' RUNFILE initialized '
#  endif
      Call Init_Run_Use()
*
* Initialize the peek/poke utility
*
      Call Init_ppu(.true.)
      Call poke_iScalar("xml opened",0)
*
* Set Grid Set
*
      Call NQ_Init()
#  ifdef _DEBUG_
      Write (6,*) ' GRID initialized '
#  endif

************************************************************************
* At this point, everything should have been properly initialized!     *
************************************************************************
* Finally, print useful runtime information about the module
      call getenvf("MOLCAS_PRINT", Prin)
      if(Prin(1:1).ne.'0'.and.Prin(1:1).ne.'S') then
      Call Print_Module_Header(ModuleName)
* Force output to be written to stdout, such that in case of problems
* later on, at least all the start info has been printed.
      Call xFlush(6)
      endif
* Write to the status file that the moduel has started
      Call StatusLine(ModuleName,' properly started!')
      Return
      End

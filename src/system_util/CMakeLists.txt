#***********************************************************************
# This file is part of OpenMolcas.                                     *
#                                                                      *
# OpenMolcas is free software; you can redistribute it and/or modify   *
# it under the terms of the GNU Lesser General Public License, v. 2.1. *
# OpenMolcas is distributed in the hope that it will be useful, but it *
# is provided "as is" and without any express or implied warranties.   *
# For more details see the full text of the license in the file        *
# LICENSE or in <http://www.gnu.org/licenses/>.                        *
#***********************************************************************

set (sources
     abendc.c
     abend.F90
     bomb_on_error.F90
     center_text.F90
     checkrun.F90
     collapseoutput.F90
     colorizeinit.F90
     constants.F90
     cwtime.F90
     data_structures.F90
     datimx.F90
     datimxc.c
     definitions.F90
     filesystem.F90
     finish.F90
     fndlnk.c
     fortran_strings.F90
     getenvc.c
     getenvf.F90
     get_genome.c
     init_linalg.F90
     iprintlevel.F90
     isstructure.F90
     locase.F90
     molcascontrol.F90
     molcascontrolinit.F90
     normal.F90
     para_info.F90
     print_module_header.F90
     printresult.F90
     putenvc.c
     quit.F90
     quit_onconverror.F90
     quit_onusererror.F90
     seconds.F90
     set_sighandlers.c
     settim.F90
     sighandler.c
     start.F90
     startlight.F90
     statusline.F90
     stdin_name.F90
     sysabendmsg.F90
     syscondmsg.F90
     sysdumpstr.F90
     sysexpand.F90
     sysfilemsg.F90
     syshalt.F90
     sysputsend.F90
     sysputs.F90
     sysputsstart.F90
     sysquitfilemsg.F90
     sysquitmsg.F90
     systemc.c
     systemf.F90
     sysvaluemsg.F90
     sysvaluewarnmsg.F90
     syswarnfilemsg.F90
     syswarnmsg.F90
     timingc.c
     timing.F90
     unixinfo.F90
     unixinfoc.c
     untested.F90
     upcase.F90
     warningmessage.F90
     warnings.F90
     whichmolcas.F90
     write_cc.c
     write_pid.c
     write_rc.c
     write_stderr.F90
     xabort.F90
     xbacktrace.F90
     xflush.F90
     xquit.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list
     constants.F90
     data_structures.F90
     definitions.F90
     filesystem.F90
     fortran_strings.F90
     para_info.F90
     unixinfo.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)

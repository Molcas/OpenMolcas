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

module warnings

use Definitions, only: iwp

implicit none
private

integer(kind=iwp) :: MaxWarnMess = -1
character(len=22), parameter :: rc_msg(0:255) = '_UKNOWN_ERROR_CODE_'
! Is there a better way to set the values?
#include "warnings.h"
data rc_msg(_RC_ALL_IS_WELL_)          / '_ALL_IS_WELL_' /
data rc_msg(_RC_JOB_KILLED_)           / '_JOB_KILLED_' /
data rc_msg(_RC_NOT_AVAILABLE_)        / '_NOT_AVAILABLE_' /
data rc_msg(_RC_CONTINUE_LOOP_)        / '_CONTINUE_LOOP_' /
data rc_msg(_RC_INVOKED_OTHER_MODULE_) / '_INVOKED_OTHER_MODULE_' /
data rc_msg(_RC_CONTINUE_UNIX_LOOP_)   / '_CONTINUE_UNIX_LOOP_' /
data rc_msg(_RC_CHO_DUM_)              / '_CHO_DUM_' /
data rc_msg(_RC_EXIT_)                 / '_EXIT_' /
data rc_msg(_RC_EXIT_EXPECTED_)        / '_EXIT_EXPECTED_' /
data rc_msg(_RC_DO_TASKS_)             / '_DO_TASKS_' /
data rc_msg(_RC_GENERAL_WARNING_)      / '_GENERAL_WARNING_' /
data rc_msg(_RC_NOT_CONVERGED_)        / '_NOT_CONVERGED_' /
data rc_msg(_RC_TIMEOUT_)              / '_TIMEOUT_' /
data rc_msg(_RC_INPUT_ERROR_)          / '_INPUT_ERROR_' /
data rc_msg(_RC_INPUT_EMIL_ERROR_)     / '_INPUT_EMIL_ERROR_' /
data rc_msg(_RC_LICENSE_)              / '_LICENSE_' /
data rc_msg(_RC_CHO_INP_)              / '_CHO_INP_' /
data rc_msg(_RC_CHECK_ERROR_)          / '_CHECK_ERROR_' /
data rc_msg(_RC_INSTALL_ERROR_)        / '_INSTALL_ERROR_' /
data rc_msg(_RC_INTERNAL_ERROR_)       / '_INTERNAL_ERROR_' /
data rc_msg(_RC_EXTERNAL_TERMINATION_) / '_EXTERNAL_TERMINATION_' /
data rc_msg(_RC_GENERAL_ERROR_)        / '_GENERAL_ERROR_' /
data rc_msg(_RC_FLOATING_EXCEPTION_)   / '_FLOATING_EXCEPTION_' /
data rc_msg(_RC_EXTERNAL_TERM_)        / '_EXTERNAL_TERM_' /
data rc_msg(_RC_MEMORY_ERROR_)         / '_MEMORY_ERROR_' /
data rc_msg(_RC_IO_ERROR_WRITE_)       / '_IO_ERROR_WRITE_' /
data rc_msg(_RC_IO_ERROR_READ_)        / '_IO_ERROR_READ_' /
data rc_msg(_RC_CHO_MEM_)              / '_CHO_MEM_' /
data rc_msg(_RC_CHO_INI_)              / '_CHO_INI_' /
data rc_msg(_RC_CHO_LOG_)              / '_CHO_LOG_' /
data rc_msg(_RC_CHO_RUN_)              / '_CHO_RUN_' /

public :: MaxWarnMess, rc_msg

end module warnings

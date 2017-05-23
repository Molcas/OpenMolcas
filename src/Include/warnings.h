/***********************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
***********************************************************************/
#ifndef _WARNINGS_H_
#define _WARNINGS_H_
#ifdef _WARNINGS_H_DOC_

This file contains return codes issued by the modules in MOLCAS.

_RC_GROUP_AND_            Mask (hex F0) to get group of return code
_RC_GROUP_INTERNAL_       Group internal rc (hex 40)
_RC_GROUP_WARNING_        Group of warnings (hex 50)
_RC_GROUP_ERROR_          Group of errors (hex 60)
_RC_GROUP_USER_ERROR_     Group of user errors (hex 70)
_RC_GROUP_CRITICAL_       Group of critical errors (hex 80)

Do an and of the RC with _RC_GROUP_AND_ to obtain the group of error or
warnings.

_RC_ALL_IS_WELL_          No return code
_RC_DO_TASKS_             Run a bunch of MOLCAS tasks
_RC_CONTINUE_LOOP_        Do another iteration in geometry optimization
_RC_INVOKED_OTHER_MODULE_ The requested module invoked another module
_RC_EXIT_                 Terminate the execution (e.g. for non
                          implemented) but with return
                          code _RC_ALL_IS_WELL_
_RC_EXIT_EXPECTED_        Same as _RC_EXIT_ but without warnings
_RC_CONTINUE_UNIX_LOOP_   Do another iteration in a unix loop.
_RC_GENERAL_WARNING_      Warning of general nature was issued
_RC_NOT_CONVERGED_        The calculation did not converge
_RC_INTERNAL_ERROR_       Some internal inconsistency of the code was
                          detected
_RC_EXTERNAL_TERMINATION_ The module was terminated by an external
                          error, such as disk full
_RC_INPUT_ERROR_          There is an error in the user input
_RC_INPUT_EMIL_ERROR_     Error in EMIL input
_RC_LICENSE_              Error in license
_RC_GENERAL_ERROR_        Error of a general or unclassified nature.
_RC_IO_ERROR_WRITE_       I/O write error
_RC_IO_ERROR_READ_        I/O read error
_RC_CHO_DUM_              Cholesky code: not really an error; e.g. stop
                          during debug.
_RC_CHO_MEM_              Cholesky code: memory allocation problems.
_RC_CHO_INI_              Cholesky code: initialization problems.
_RC_CHO_LOG_              Cholesky code: logical error (usually a bug).
_RC_CHO_RUN_              Cholesky code: run-time error (or a bug).
_RC_CHO_INP_              Cholesky code: input error.
_RC_MEMORY_ERROR_         Memory error.
_RC_JOB_KILLED_           Job terminated by user
_RC_INSTALL_ERROR_        Installation error
_RC_TIMEOUT_              Job has reached user-specified time limit.
#endif

#define _RC_GROUP_AND_            240
#define _RC_GROUP_INTERNAL_        64
#define _RC_GROUP_WARNING_         80
#define _RC_GROUP_ERROR_           96
#define _RC_GROUP_USER_ERROR_     112
#define _RC_GROUP_CRITICAL_       128

#define _RC_ALL_IS_WELL_            0
#define _RC_JOB_KILLED_             1
#define _RC_NOT_AVAILABLE_         36

#define _RC_CONTINUE_LOOP_         64
#define _RC_INVOKED_OTHER_MODULE_  65
#define _RC_CONTINUE_UNIX_LOOP_    66
#define _RC_CHO_DUM_               67
#define _RC_EXIT_                  68
#define _RC_EXIT_EXPECTED_         69
#define _RC_DO_TASKS_              70
#define _RC_GENERAL_WARNING_       80

#define _RC_NOT_CONVERGED_         96
#define _RC_TIMEOUT_              100

#define _RC_INPUT_ERROR_          112
#define _RC_INPUT_EMIL_ERROR_     113
#define _RC_LICENSE_              114
#define _RC_CHO_INP_              115
#define _RC_CHECK_ERROR_          116
#define _RC_INSTALL_ERROR_        117

#define _RC_INTERNAL_ERROR_       128
#define _RC_EXTERNAL_TERMINATION_ 129
#define _RC_GENERAL_ERROR_        130
#define _RC_FLOATING_EXCEPTION_   134
#define _RC_EXTERNAL_TERM_        137
#define _RC_MEMORY_ERROR_         139
#define _RC_IO_ERROR_WRITE_       161
#define _RC_IO_ERROR_READ_        162
#define _RC_CHO_MEM_              163
#define _RC_CHO_INI_              164
#define _RC_CHO_LOG_              165
#define _RC_CHO_RUN_              166
#endif

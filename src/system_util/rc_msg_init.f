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
      subroutine rc_msg_init
      implicit none
#include "warnings.fh"
      rc_msg                           ='_UNKNOWN_ERROR_CODE_          '

      rc_msg(_RC_ALL_IS_WELL_         )='_ALL_IS_WELL_                 '
      rc_msg(_RC_JOB_KILLED_          )='_JOB_KILLED_                  '
      rc_msg(_RC_CONTINUE_LOOP_       )='_CONTINUE_LOOP_               '
      rc_msg(_RC_INVOKED_OTHER_MODULE_)='_INVOKED_OTHER_MODULE_        '
      rc_msg(_RC_CONTINUE_UNIX_LOOP_  )='_CONTINUE_UNIX_LOOP_          '
      rc_msg(_RC_CHO_DUM_             )='_CHO_DUM_                     '
      rc_msg(_RC_EXIT_                )='_EXIT_                        '
      rc_msg(_RC_EXIT_EXPECTED_       )='_EXIT_EXPECTED_               '
      rc_msg(_RC_DO_TASKS_            )='_DO_TASKS_                    '
      rc_msg(_RC_GENERAL_WARNING_     )='_GENERAL_WARNING_             '
      rc_msg(_RC_NOT_CONVERGED_       )='_NOT_CONVERGED_               '
      rc_msg(_RC_TIMEOUT_             )='_TIMEOUT_                     '
      rc_msg(_RC_INPUT_ERROR_         )='_INPUT_ERROR_                 '
      rc_msg(_RC_INPUT_EMIL_ERROR_    )='_INPUT_EMIL_ERROR_            '
      rc_msg(_RC_LICENSE_             )='_LICENSE_                     '
      rc_msg(_RC_CHO_INP_             )='_CHO_INP_                     '
      rc_msg(_RC_CHECK_ERROR_         )='_CHECK_ERROR_                 '
      rc_msg(_RC_INSTALL_ERROR_       )='_INSTALL_ERROR_               '
      rc_msg(_RC_INTERNAL_ERROR_      )='_INTERNAL_ERROR_              '
      rc_msg(_RC_EXTERNAL_TERMINATION_)='_EXTERNAL_TERMINATION_        '
      rc_msg(_RC_GENERAL_ERROR_       )='_GENERAL_ERROR_               '
      rc_msg(_RC_FLOATING_EXCEPTION_  )='_FLOATING_EXCEPTION_          '
      rc_msg(_RC_EXTERNAL_TERM_       )='_EXTERNAL_TERM_               '
      rc_msg(_RC_MEMORY_ERROR_        )='_MEMORY_ERROR_                '
      rc_msg(_RC_IO_ERROR_WRITE_      )='_IO_ERROR_WRITE_              '
      rc_msg(_RC_IO_ERROR_READ_       )='_IO_ERROR_READ_               '
      rc_msg(_RC_CHO_MEM_             )='_CHO_MEM_                     '
      rc_msg(_RC_CHO_INI_             )='_CHO_INI_                     '
      rc_msg(_RC_CHO_LOG_             )='_CHO_LOG_                     '
      rc_msg(_RC_CHO_RUN_             )='_CHO_RUN_                     '
      end

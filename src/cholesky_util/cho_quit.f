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
      SUBROUTINE CHO_QUIT(STRING,IERR)
C
C     Purpose: echo message STRING and abort execution.
C
#include "implicit.fh"
      CHARACTER*(*) STRING
#include "cholesky.fh"
      WRITE(LUPRI,'(//,1X,A)') '***'
      IF (IERR.EQ.0 .OR. IERR.EQ.100) THEN
         WRITE(LUPRI,'(1X,A)')
     &   '*** Execution stopped by Cholesky Decomposition Utility'
         WRITE(LUPRI,'(1X,A,A)')
     &   '*** Message: ',STRING
      ELSE
         WRITE(LUPRI,'(1X,A)')    '*** Error in Cholesky Core Routine'
         WRITE(LUPRI,'(1X,A,A)')  '*** Message: ',STRING
         WRITE(LUPRI,'(1X,A,I5)') '*** Code   : ',IERR
      END IF
      WRITE(LUPRI,'(1X,A,//)') '***'
      CALL CHO_TRANSLATEERRORCODE(IERR,MOLCASCODE)
      CALL QUIT(MOLCASCODE)

      END

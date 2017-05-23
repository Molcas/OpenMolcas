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
      SUBROUTINE CHO_GETRSTC()
C
C     Purpose: read and check decomposition restart info.
C
#include "implicit.fh"
#include "cholesky.fh"

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CHO_GETRSTC')

C     Read restart file, populating restart common block.
C     ---------------------------------------------------

      IFAIL = 0
      CALL CHO_RDRSTC(IFAIL)
      IF (IFAIL .NE. 0) THEN
         WRITE(LUPRI,'(A,A)')
     &   SECNAM,': error reading decomposition restart file.'
         WRITE(LUPRI,'(A,A,I10)')
     &   SECNAM,': return code from reading routine:',IFAIL
         CALL CHO_QUIT('Error reading decomposition restart file',104)
      END IF

C     Check system info .
C     -------------------

      IFAIL = 0
      CALL CHO_RSTMOL(IFAIL)
      IF (IFAIL .NE. 0) THEN
         WRITE(LUPRI,'(A,A)')
     &   SECNAM,': decomposition restart failure.'
         CALL CHO_QUIT('Decomposition restart failure in '//SECNAM,105)
      END IF

C     Check decomposition configuration.
C     ----------------------------------

      IFAIL = 0
      CALL CHO_RSTCNF(IFAIL)
      IF (IFAIL .NE. 0) THEN
         WRITE(LUPRI,'(A,A,I6,A)')
     &   SECNAM,':',IFAIL,' configuration discrepancies detected.'
         IF (MODRST .EQ. -1) THEN
            WRITE(LUPRI,'(A)')
     &      'Recovery: using configuration from restart file.'
            CALL CHO_RESETCNF
         ELSE IF (MODRST .EQ. 0) THEN
            WRITE(LUPRI,'(A)')
     &      'Recovery: none, program stops.'
            CALL CHO_QUIT('Restart configuration error in '//SECNAM,105)
         ELSE IF (MODRST .EQ. 1) THEN
            WRITE(LUPRI,'(A)')
     &      'Recovery: using input configuration.'
         ELSE
            WRITE(LUPRI,'(A,A,I6,A)')
     &      SECNAM,': restart model,',MODRST,', not recognized.'
            CALL CHO_QUIT('Error in '//SECNAM,103)
         END IF
      END IF

      END

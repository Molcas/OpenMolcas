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
      SUBROUTINE CHO_MCA_INT_1_DBG(DIAG,LEVEL)
C
C     Purpose: debug seward interface routine CHO_MCA_INT_1.
C
C     LEVEL =  1: test diagonal, reduced set 1 (i.e. initial).
C              2: test diagonal, reduced set 2 (i.e. current).
C              3: test symmetry of integral matrix (shell quadruple-based)
C
#include "implicit.fh"
      DIMENSION DIAG(*)
#include "cholesky.fh"

      CHARACTER*17 SECNAM
      PARAMETER (SECNAM = 'CHO_MCA_INT_1_DBG')

      LOGICAL LOCDIAG, LOCSYM

      CALL CHO_HEAD('Debugging CHO_MCA_INT_1','=',80,LUPRI)
      WRITE(LUPRI,'(A,I2)') 'Debug level',LEVEL

      IF (LEVEL .EQ. 1) THEN
         LOCDIAG = .TRUE.
         LOCSYM  = .FALSE.
         IRED    = 1
      ELSE IF (LEVEL .EQ. 2) THEN
         LOCDIAG = .TRUE.
         LOCSYM  = .FALSE.
         IRED    = 2
      ELSE IF (LEVEL .EQ. 3) THEN
         LOCDIAG = .FALSE.
         LOCSYM  = .TRUE.
      ELSE
         LOCDIAG = .FALSE.
         LOCSYM  = .FALSE.
         WRITE(LUPRI,'(A)') 'Debug level not recognized ---',
     &                      ' debug cancelled!'
      END IF

      IF (LOCDIAG) THEN
         CALL CHO_MCA_INT_1_DBG1(DIAG,IRED)
      END IF

      IF (LOCSYM) THEN
         CALL CHO_MCA_INT_1_DBG2()
      END IF

      END

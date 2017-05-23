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
      SUBROUTINE CHO_RSTMOL(NERR)
C
C     Purpose: check restart molecular info.
C
      IMPLICIT NONE
      INTEGER NERR
#include "cholesky.fh"
#include "choorb.fh"

      INTEGER ISYM

      NERR = 0

      IF (XNSYM .NE. NSYM) THEN
         WRITE(LUPRI,'(A,I3,A,I3)')
     &   'RESTART ERROR: #irreps from restart file:',XNSYM,
     &   ' Expected:',NSYM
         NERR = NERR + 1
      ELSE
         DO ISYM = 1,NSYM
            IF (XNBAS(ISYM) .NE. NBAS(ISYM)) THEN
               WRITE(LUPRI,'(A,I2,A,I9,A,I9)')
     &         'RESTART ERROR: #basis functions (sym.',ISYM,
     &         ') from restart file:',XNBAS(ISYM),
     &         ' Expected:',NBAS(ISYM)
               NERR = NERR + 1
            END IF
         END DO
      END IF

      IF (XNSHELL .NE. NSHELL) THEN
         WRITE(LUPRI,'(A,I9,A,I9)')
     &   'RESTART ERROR: #shells from restart file:',XNSHELL,
     &   ' Expected:',NSHELL
         NERR = NERR + 1
      END IF

      IF (XNNSHL .NE. NNSHL) THEN
         WRITE(LUPRI,'(A,I9,A,I9)')
     &   'RESTART ERROR: #shell pairs from restart file:',XNNSHL,
     &   ' Expected:',NNSHL
         NERR = NERR + 1
      END IF

      END

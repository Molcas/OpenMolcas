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
      SUBROUTINE CHO_PRTTIM(SECTION,TCPU2,TCPU1,TWALL2,TWALL1,IOPT)
C
C     Purpose: print timing for a section.
C
      IMPLICIT NONE
      CHARACTER*(*) SECTION
      REAL*8        TCPU1, TCPU2, TWALL1, TWALL2
      INTEGER       IOPT
#include "cholesky.fh"

      INTEGER IHRC, IMNC
      INTEGER IHRW, IMNW
      INTEGER LENSEC
      REAL*8  TCPUT, TWALLT, SECC, SECW
      CHARACTER*80 STRING

      TCPUT  = TCPU2  - TCPU1
      TWALLT = TWALL2 - TWALL1
      CALL CHO_CNVTIM(TCPUT,IHRC,IMNC,SECC)
      CALL CHO_CNVTIM(TWALLT,IHRW,IMNW,SECW)

      IF (IOPT .EQ. 0) THEN
         LENSEC = LEN(SECTION)
         WRITE(LUPRI,'(/,A,A,A)')
     &   '***** ',SECTION(1:LENSEC),' completed *****'
         WRITE(LUPRI,'(A,I8,A,I2,A,F6.2,A)')
     &   'Total CPU  time:',IHRC,' hours ',IMNC,' minutes ',SECC,
     &   ' seconds'
         WRITE(LUPRI,'(A,I8,A,I2,A,F6.2,A,/)')
     &   'Total wall time:',IHRW,' hours ',IMNW,' minutes ',SECW,
     &   ' seconds'
      ELSE IF (IOPT .EQ. 1) THEN
         LENSEC = LEN(SECTION)
         WRITE(LUPRI,'(///,A,A,A)')
     &   '***** ',SECTION(1:LENSEC),' completed *****'
         WRITE(LUPRI,'(A,I8,A,I2,A,F6.2,A)')
     &   'Total CPU  time:',IHRC,' hours ',IMNC,' minutes ',SECC,
     &   ' seconds'
         WRITE(LUPRI,'(A,I8,A,I2,A,F6.2,A,//)')
     &   'Total wall time:',IHRW,' hours ',IMNW,' minutes ',SECW,
     &   ' seconds'
      ELSE IF (IOPT .EQ. 2) THEN
         LENSEC = MIN(LEN(SECTION),70)
         WRITE(STRING,'(A10,A)') 'Timing of ',SECTION(1:LENSEC)
         LENSEC = LENSEC + 10
         CALL CHO_HEAD(STRING(1:LENSEC),'=',80,LUPRI)
         WRITE(LUPRI,'(/,A,I8,A,I2,A,F6.2,A)')
     &   'Total CPU  time:',IHRC,' hours ',IMNC,' minutes ',SECC,
     &   ' seconds'
         WRITE(LUPRI,'(A,I8,A,I2,A,F6.2,A)')
     &   'Total wall time:',IHRW,' hours ',IMNW,' minutes ',SECW,
     &   ' seconds'
      ELSE
         WRITE(LUPRI,'(/,A,I8,A,I2,A,F6.2,A)')
     &   'Total CPU  time:',IHRC,' hours ',IMNC,' minutes ',SECC,
     &   ' seconds'
         WRITE(LUPRI,'(A,I8,A,I2,A,F6.2,A)')
     &   'Total wall time:',IHRW,' hours ',IMNW,' minutes ',SECW,
     &   ' seconds'
      END IF

      CALL CHO_FLUSH(LUPRI)

      END

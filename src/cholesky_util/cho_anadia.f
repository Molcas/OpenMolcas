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
      SUBROUTINE CHO_ANADIA(DIAG,BIN1,STEP,NUMBIN,FULL)
C
C     Purpose: analyze diagonal (histogram).
C
#include "implicit.fh"
      DIMENSION DIAG(*)
      LOGICAL   FULL
#include "cholesky.fh"

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_ANADIA')

      PARAMETER (MAXBIN = 50, NUMSTA = 7)
      DIMENSION BIN(MAXBIN), STAT(NUMSTA)

      LOGICAL FOUND

#if defined (_DEBUGPRINT_)
      CALL QENTER('_ANADIA')
#endif

C     Print header.
C     -------------

      CALL CHO_HEAD('Histogram of Diagonal Elements','=',80,LUPRI)

C     Set up size bins for analysis of diagonal.
C     ------------------------------------------

      IF (NUMBIN .LT. 1) THEN
         MBIN = MIN(10,MAXBIN)
         BINLOC = 1.0D2
         STPLOC = 1.0D-2
      ELSE
         MBIN = MIN(NUMBIN,MAXBIN)
         BINLOC = BIN1
         STPLOC = STEP
      END IF

      BIN(1) = BINLOC
      DO IBIN = 2,MBIN
         BIN(IBIN) = BIN(IBIN-1)*STPLOC
      END DO

C     Set smallest BIN according to full.
C     -----------------------------------

      IF (FULL) THEN
         NBIN = MBIN
      ELSE
         NBIN  = MBIN
         IBIN  = MBIN
         FOUND = .FALSE.
         DO WHILE ((IBIN.GT.1) .AND. (.NOT.FOUND))
            IBIN = IBIN - 1
            IF (THRCOM .GE. BIN(IBIN)) THEN
               NBIN = IBIN + 1
            ELSE
               FOUND = .TRUE.
            END IF
         END DO
      END IF

C     Histogram.
C     ----------

      CALL CHO_ANASIZE(DIAG,NNBSTRT(1),BIN,NBIN,LUPRI)

C     Count converged.
C     ----------------

      NCONV = 0
      DO IAB = 1,NNBSTRT(1)
         IF (DIAG(IAB) .LE. THRCOM) NCONV = NCONV + 1
      END DO
      WRITE(LUPRI,'(/,1X,A,I10,/,1X,A,I10)')
     & 'Converged  : ',NCONV,'Unconverged: ',NNBSTRT(1)-NCONV

C     Print total number of negative zeroed diagonal as well as the most
C     negative one.
C     ------------------------------------------------------------------

      WRITE(LUPRI,'(/,1X,A,5X,I10)')
     & 'Total number of zeroed negative diagonals: ',NNZTOT
      IF (NNZTOT .GT. 0) THEN
         IF (IABMNZ .LT. 1) THEN
            WRITE(LUPRI,'(1X,A)')
     &     'WARNING: most negative zeroed diagonal has not been stored!'
         ELSE
            WRITE(LUPRI,'(1X,A,1P,D15.6)')
     &      '- most negative zeroed diagonal          : ',DIAMNZ
         END IF
      END IF

C     Print statistics.
C     -----------------

      CALL STATISTICS(DIAG,NNBSTRT(1),STAT,1,2,3,4,5,6,7)
      WRITE(LUPRI,'(/,1X,A,1P,D15.6)')
     & 'Minimum diagonal: ',STAT(3)
      WRITE(LUPRI,'(1X,A,1P,D15.6)')
     & 'Maximum diagonal: ',STAT(4)
      WRITE(LUPRI,'(1X,A,1P,D15.6)')
     & 'Mean value      : ',STAT(1)
      WRITE(LUPRI,'(1X,A,1P,D15.6)')
     & 'Mean abs. value : ',STAT(2)
      WRITE(LUPRI,'(1X,A,1P,D15.6)')
     & 'Biased variance : ',STAT(6)
      WRITE(LUPRI,'(1X,A,1P,D15.6,A)')
     & 'Standard dev.   : ',STAT(7),' (unbiased variance)'

#if defined (_DEBUGPRINT_)
      CALL QEXIT('_ANADIA')
#endif

      END

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
      SUBROUTINE CHO_CHKDIA(DIAG,ISYM,XM,YM,ZM,NNEGT,NNEG,NCONV)
C
C     Purpose: 1) find min. in updated diagonal, XM (sym. ISYM only)
C              2) find max. in updated diagonal, YM (sym. ISYM only)
C              3) find abs. max. in updated diagonal, ZM (sym. ISYM only)
C              4) count #diagonals < 0.0D0, NNEGT
C              5) count #diagonals < THRNEG, NNEGT
C              6) count #screenable diagonals, NCONV
C
C     -- also: a) zero negative diagonals < THRNEG (from cholesky.fh)
C              b) screen diagonal if requested (flag SCDIAG from cholesky.fh)
C              c) Keep track of most negative zeroed diagonal.
C
      use ChoSwp, only: IndRed
#include "implicit.fh"
      DIMENSION DIAG(*)
#include "cholesky.fh"

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_CHKDIA')

      PARAMETER (ZERO = 0.0D0)

C     Initialization.
C     ---------------

      NNEGT = 0
      NNEG  = 0
      NCONV = 0

      IF (NNBSTR(ISYM,2) .GT. 0) THEN
         JAB1 = IIBSTR(ISYM,2) + 1
         JAB2 = JAB1 + NNBSTR(ISYM,2) - 1
         XM   = DIAG(INDRED(JAB1,2))
         YM   = DIAG(INDRED(JAB1,2))
         ZM   = ABS(YM)
      ELSE
         XM = ZERO
         YM = ZERO
         ZM = ZERO
         RETURN
      END IF

C     Find min. and max. diagonal and zero too negative diagonals.
C     ------------------------------------------------------------

      DO JAB = JAB1,JAB2
         IAB = INDRED(JAB,2)  ! get address in first red. set
         XM  = MIN(XM,DIAG(IAB))
         YM  = MAX(YM,DIAG(IAB))
         IF (DIAG(IAB) .LT. ZERO) THEN
            NNEGT = NNEGT + 1
            IF (DIAG(IAB) .LT. THRNEG) THEN
               NNEG   = NNEG + 1
               IF (DIAG(IAB) .LT. TOONEG) THEN
                  WRITE(LUPRI,'(A,A,I12,1X,1P,D16.8)')
     &            SECNAM,': diagonal too negative: ',IAB,DIAG(IAB)
                  WRITE(LUPRI,'(A,A)')
     &            SECNAM,': shutting down Cholesky decomposition!'
                  CALL CHO_QUIT('Diagonal too negative in '//SECNAM,104)
               END IF
               IF (DIAG(IAB) .LT. WARNEG) THEN
                  WRITE(LUPRI,'(A,A,I12,1X,1P,D16.8,A)')
     &            SECNAM,': Negative diagonal: ',IAB,DIAG(IAB),
     &            ' (zeroed)'
               END IF
               IF (DIAG(IAB) .LT. DIAMNZ) THEN
                  DIAMNZ = DIAG(IAB)
                  IABMNZ = IAB
               END IF
               DIAG(IAB) = ZERO
            END IF
         END IF
      END DO
      ZM = MAX(ABS(XM),ABS(YM))

C     Screen diagonal (if requested) and count the screenables as
C     converged.
C     NOTE: the screening is actually identical to setting up
C           reduced sets. However, doing the screening here will
C           force entries in later vectors of this integral pass
C           to have zero entries at screened diagonals.
C     -----------------------------------------------------------

      IF (SCDIAG) THEN
         DO JAB = JAB1,JAB2
            IAB = INDRED(JAB,2)  ! get address in first red. set
            TST = SQRT(ABS(DIAG(IAB))*ZM)*DAMP(2)
            IF (TST .LE. THRCOM) THEN
               NCONV     = NCONV + 1
               DIAG(IAB) = ZERO
            END IF
         END DO
      ELSE
         DO JAB = JAB1,JAB2
            IAB = INDRED(JAB,2)  ! get address in first red. set
            TST = SQRT(ABS(DIAG(IAB))*ZM)*DAMP(2)
            IF (TST .LE. THRCOM) THEN
               NCONV = NCONV + 1
            END IF
         END DO
      END IF

      END

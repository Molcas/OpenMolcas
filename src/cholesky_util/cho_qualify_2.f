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
      SUBROUTINE CHO_QUALIFY_2(DIAG,ISYM,ISHLAB,MEM,MEM0,LEFT)
C
C     Purpose: qualify diagonals ("qualify until full, then largest").
C
      use ChoSwp, only: iQuAB, nnBstRSh, iiBstRSh, IndRed
#include "implicit.fh"
      DIMENSION DIAG(*)
#include "cholesky.fh"

      CHARACTER*13 SECNAM
      PARAMETER (SECNAM = 'CHO_QUALIFY_2')

      NDIM = NNBSTRSH(ISYM,ISHLAB,2)
      IF (NDIM .GT. 0) THEN
         MAXQ = MIN(MAXQUAL-NQUAL(ISYM),LEFT/NNBSTR(ISYM,2))
         NUMQ = 0
         IF (MAXQ .GT. 0) THEN
            I1 = IIBSTR(ISYM,2) + IIBSTRSH(ISYM,ISHLAB,2) + 1
            I2 = I1 + NDIM - 1
            IF (MAXQ .EQ. 1) THEN ! qualify the largest > DIAMIN
               XMAX = DIAMIN(ISYM)
               IMAX = -1
               DO I = I1,I2
                  J = INDRED(I,2)
                  IF (DIAG(J) .GE. XMAX) THEN
                     XMAX = DIAG(J)
                     IMAX = I
                  END IF
               END DO
               IF (IMAX .GT. 0) THEN
                  NUMQ = NUMQ + 1
                  iQuAB(IOFFQ(ISYM)+NUMQ,ISYM) = IMAX
               END IF
            ELSE ! full search
               DO I = I1,I2
                  J = INDRED(I,2)
                  IF (DIAG(J) .GE. DIAMIN(ISYM)) THEN
                     IF (NUMQ .LT. MAXQ) THEN
                        NUMQ = NUMQ + 1
                        iQuAB(IOFFQ(ISYM)+NUMQ,ISYM) = I
                     ELSE IF (NUMQ .EQ. MAXQ) THEN
                        K1  = IOFFQ(ISYM) + 1
                        K2  = K1 + NUMQ - 1
                        II1 = IQUAB(K1,ISYM)
                        JJ1 = INDRED(II1,2)
                        XMIN = DIAG(JJ1)
                        KKMN = K1
                        DO K = K1+1,K2 ! find min. among qualified
                           II = IQUAB(K,ISYM)
                           JJ = INDRED(II,2)
                           IF (DIAG(JJ) .LT. XMIN) THEN
                              XMIN = DIAG(JJ)
                              KKMN = K
                           END IF
                        END DO
                        IF (DIAG(J) .GT. XMIN) THEN ! replace
                           iQuAB(KKMN,ISYM) = I
                        END IF
                     ELSE
                        CALL CHO_QUIT('Logical error in '//SECNAM,
     &                                104)
                     END IF
                  END IF
               END DO
            END IF
         END IF
         NQUAL(ISYM) = NQUAL(ISYM) + NUMQ
         MEM0 = MEM0 + NUMQ*NNBSTR(ISYM,2)
         LEFT = MEM  - MEM0
      END IF

      END

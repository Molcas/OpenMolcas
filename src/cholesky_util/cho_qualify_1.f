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
      SUBROUTINE CHO_QUALIFY_1(DIAG,ISYM,ISHLAB,MEM,MEM0,LEFT)
C
C     Purpose: qualify diagonals ("qualify until full").
C
      use ChoSwp, only: iQuAB, nnBstRSh, iiBstRSh, IndRed
#include "implicit.fh"
      DIMENSION DIAG(*)
#include "cholesky.fh"

      IF (NNBSTRSH(ISYM,ISHLAB,2) .GT. 0) THEN
         I  = IIBSTR(ISYM,2) + IIBSTRSH(ISYM,ISHLAB,2)
         I2 = I + NNBSTRSH(ISYM,ISHLAB,2)
         MAXQ = MIN(MAXQUAL-NQUAL(ISYM),LEFT/NNBSTR(ISYM,2))
         NUMQ = 0
         DO WHILE ((I.LT.I2) .AND. (NUMQ.LT.MAXQ))
            I = I + 1
            J = INDRED(I,2)
            IF (DIAG(J) .GE. DIAMIN(ISYM)) THEN
               NUMQ = NUMQ + 1
               iQuAB(IOFFQ(ISYM)+NUMQ,ISYM)=I
            END IF
         END DO
         NQUAL(ISYM) = NQUAL(ISYM) + NUMQ
         MEM0 = MEM0 + NUMQ*NNBSTR(ISYM,2)
         LEFT = MEM  - MEM0
      END IF

      END

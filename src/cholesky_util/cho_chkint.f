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
      SUBROUTINE CHO_CHKINT(XINT,DIAG,ISYM,NERR,TOL,REPORT)
C
C     Purpose: check diagonals in qualified integral columns.
C
      use ChoSwp, only: iQuAB, IndRed
#include "implicit.fh"
      DIMENSION XINT(*), DIAG(*)
      LOGICAL   REPORT
#include "cholesky.fh"
#include "choptr.fh"

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_CHKINT')

      NERR = 0
      DO I = 1,NQUAL(ISYM)
         II = IQUAB(I,ISYM)
         JJ = INDRED(II,2)
         IK = II - IIBSTR(ISYM,2)
         KK = NNBSTR(ISYM,2)*(I - 1) + IK
         DF = DIAG(JJ) - XINT(KK)
         IF (ABS(DF) .GT. TOL) THEN
            NERR = NERR + 1
            IF (REPORT) THEN
               WRITE(LUPRI,*) SECNAM,': diag error: ',DIAG(JJ),XINT(KK)
               WRITE(LUPRI,*) '            diagonal elm    : ',JJ,
     &                        ' (rs1) ',II,' (rs2)'
               WRITE(LUPRI,*) '            integral row,col: ',IK,I
            END IF
         END IF
      END DO

      END

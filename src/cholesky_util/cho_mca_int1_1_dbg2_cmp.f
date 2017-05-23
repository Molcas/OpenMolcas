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
      SUBROUTINE CHO_MCA_INT1_1_DBG2_CMP(XINT1,XINT2,NI,NJ,
     &                                   ERRMIN,IMN,JMN,
     &                                   ERRMAX,IMX,JMX,
     &                                   ITST,IERR,THR,PRTERR,LUPRI)
#include "implicit.fh"
      DIMENSION XINT1(NI,NJ), XINT2(NJ,NI)
      LOGICAL PRTERR

      IF ((NI.LT.1) .OR. (NJ.LT.1)) THEN
         ERRMAX = 0.0D0
         ERRMIN = 0.0D0
         IMN = 0
         JMN = 0
         IMX = 0
         JMX = 0
         RETURN
      END IF

      ERRMAX = XINT1(1,1) - XINT2(1,1)
      ERRMIN = XINT1(1,1) - XINT2(1,1)
      IMN = 1
      JMN = 1
      IMX = 1
      JMX = 1

      JERR = 0
      DO J = 1,NJ
         DO I = 1,NI
            ITST = ITST + 1
            DIFF = XINT1(I,J) - XINT2(J,I)
            IF (ABS(DIFF) .GT. THR) THEN
               JERR = JERR + 1
               IF (PRTERR) THEN
                  WRITE(LUPRI,*) '      Error: ',I,J,DIFF
               END IF
            END IF
            IF (DIFF .LT. ERRMIN) THEN
               ERRMIN = DIFF
               IMN = I
               JMN = J
            END IF
            IF (DIFF .GT. ERRMAX) THEN
               ERRMAX = DIFF
               IMX = I
               JMX = J
            END IF
         END DO
      END DO
      IERR = IERR + JERR

      IF ((JERR.NE.0) .AND. (NI.EQ.NJ)) THEN
         JERR = 0
         WRITE(LUPRI,*) '         Checking for identity...'
         DO J = 1,NJ
            DO I = 1,NI
               DIFF = XINT1(I,J) - XINT2(I,J)
               IF (ABS(DIFF) .GT. 1.0D-14) JERR = JERR + 1
            END DO
         END DO
         IF (JERR .NE. 0) THEN
            WRITE(LUPRI,*) '      ...not identical!!'
         ELSE
            WRITE(LUPRI,*) '      ...identical!!'
         END IF
      END IF

      END

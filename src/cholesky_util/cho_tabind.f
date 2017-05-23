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
      INTEGER FUNCTION CHO_TABIND(TABLE,LKEY,NTABLE,EOINP,LEOINP,NEOINP,
     &                            WORD)
C
C     Purpose: table lookup.
C              First, try to find WORD in TABLE. If success, return ID,
C              else, check if WORD is a special string
C              in EOINP (if any supplied). If success, return NTABLE+1,
C              else, return -1.
C
      IMPLICIT NONE
      INTEGER LKEY, NTABLE, LEOINP, NEOINP
      CHARACTER*(*) TABLE(NTABLE)  ! <-- character*(lkey)
      CHARACTER*(*) EOINP(NEOINP)  ! <-- character*(leoinp)
      CHARACTER*(*) WORD           ! <-- character*(lkey)

      INTEGER IJUMP, LCMP

C     Find entry.
C     -----------

      IF (LKEY.GT.0 .AND. NTABLE.GT.0) THEN
         IJUMP = 1
         DO WHILE (IJUMP.LE.NTABLE .AND. TABLE(IJUMP).NE.WORD)
            IJUMP = IJUMP + 1
         END DO
         IF (IJUMP .GT. NTABLE) THEN
            IF (LEOINP.GT.0 .AND. NEOINP.GT.0) THEN
               LCMP  = MIN(LEOINP,LKEY)
               IJUMP = 1
               DO WHILE (IJUMP.LE.NEOINP .AND.
     &                   EOINP(IJUMP)(1:LCMP).NE.WORD(1:LCMP))
                  IJUMP = IJUMP + 1
               END DO
               IF (IJUMP .GT. NEOINP) THEN
                  CHO_TABIND = -1
               ELSE
                  CHO_TABIND = NTABLE + 1
               END IF
            ELSE
               CHO_TABIND = -1
            END IF
         ELSE
            CHO_TABIND = IJUMP
         END IF
      ELSE
         CHO_TABIND = -1
      END IF

      END

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
      SUBROUTINE CHO_SETSH2(ISHLSO,ISOSHL,NBSTSH,NBAST,NSHELL)
      IMPLICIT NONE
      INTEGER NBAST, NSHELL
      INTEGER ISHLSO(NBAST), ISOSHL(NBAST), NBSTSH(NSHELL)

      INTEGER I, ICOUNT, ISHL

      DO ISHL = 1,NSHELL
         ICOUNT = 0
         I      = 0
         DO WHILE ((I.LT.NBAST) .AND. (ICOUNT.LT.NBSTSH(ISHL)))
            I = I + 1
            IF (ISOSHL(I) .EQ. ISHL) THEN
               ICOUNT = ICOUNT + 1
               ISHLSO(I) = ICOUNT
            END IF
         END DO
      END DO

      END

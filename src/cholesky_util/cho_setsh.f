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
      SUBROUTINE CHO_SETSH(IBASSH,NBASSH,NBSTSH,IBAS,NBAS,ISOSHL,
     &                     NSYM,NSHELL,NBAST)
      IMPLICIT NONE
      INTEGER NSYM, NSHELL, NBAST
      INTEGER IBASSH(NSYM,NSHELL), NBASSH(NSYM,NSHELL)
      INTEGER NBSTSH(NSHELL)
      INTEGER ISOSHL(NBAST)
      INTEGER IBAS(NSYM), NBAS(NSYM)

      INTEGER ISYM, IA, ISHL

      CALL CHO_IZERO(NBASSH,NSYM*NSHELL)
      DO ISYM = 1,NSYM
         DO IA = 1,NBAS(ISYM)
            ISHL = ISOSHL(IBAS(ISYM)+IA)
            NBASSH(ISYM,ISHL) = NBASSH(ISYM,ISHL) + 1
         END DO
      END DO

      DO ISHL = 1,NSHELL
         IBASSH(1,ISHL) = 0
         NBSTSH(ISHL) = NBASSH(1,ISHL)
         DO ISYM = 2,NSYM
            IBASSH(ISYM,ISHL) = NBSTSH(ISHL)
            NBSTSH(ISHL) = NBSTSH(ISHL) + NBASSH(ISYM,ISHL)
         END DO
      END DO

      END

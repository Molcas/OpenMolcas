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
      SUBROUTINE CHO_SETREDIND(IIBSTRSH,NNBSTRSH,MSYM,MMSHL,IRED)
C
C     Purpose: set index arrays for reduced set IRED. The counter
C              array NNBSTRSH must be set on entry.
C
#include "implicit.fh"
      INTEGER IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM,MMSHL,3)
#include "cholesky.fh"

      CHARACTER*13 SECNAM
      PARAMETER (SECNAM = 'CHO_SETREDIND')

      J = IRED

#if defined (_DEBUGPRINT_)
      IF ((NNSHL.NE.MMSHL) .OR. (NSYM.NE.MSYM))
     & CALL CHO_QUIT('[1] Dimension error in '//SECNAM,104)
      IF ((J.LT.1) .OR. (J.GT.3))
     & CALL CHO_QUIT('[2] Dimension error in '//SECNAM,104)
#endif

      IF (NNSHL .LT. 1) THEN ! may occur in parallel runs
         NNBSTRT(J) = 0
         CALL CHO_IZERO(IIBSTR(1,J),NSYM)
         CALL CHO_IZERO(NNBSTR(1,J),NSYM)
         RETURN
      END IF

      NNBSTRT(J) = 0
      DO ISYM = 1,NSYM
         IIBSTRSH(ISYM,1,J) = 0
         NNBSTR(ISYM,J) = NNBSTRSH(ISYM,1,J)
         DO ISHLAB = 2,NNSHL
            IIBSTRSH(ISYM,ISHLAB,J) = NNBSTR(ISYM,J)
            NNBSTR(ISYM,J) = NNBSTR(ISYM,J) + NNBSTRSH(ISYM,ISHLAB,J)
         END DO
         IIBSTR(ISYM,J) = NNBSTRT(J)
         NNBSTRT(J) = NNBSTRT(J) + NNBSTR(ISYM,J)
      END DO

      END

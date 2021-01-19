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
      SUBROUTINE CHO_GETSTOR(VCSTOR)
C
C     Purpose: get total vector storage (in words).
C
#include "implicit.fh"
      DIMENSION VCSTOR(*)
#include "cholesky.fh"

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CHO_GETSTOR')

      DO ISYM = 1,NSYM
         IF (NUMCHO(ISYM) .GT. MAXVEC) THEN
            WRITE(LUPRI,*) SECNAM,': too many Cholesky vectors ',
     &                     'in symmetry ',ISYM,': ',NUMCHO(ISYM)
            CALL CHO_QUIT('Error in '//SECNAM,103)
            VCSTOR(ISYM) = 0.0D0
         ELSE IF (NUMCHO(ISYM) .LT. 0) THEN
            WRITE(LUPRI,*) SECNAM,': negative #Cholesky vectors ',
     &                     'in symmetry ',ISYM,': ',NUMCHO(ISYM)
            CALL CHO_QUIT('Error in '//SECNAM,103)
            VCSTOR(ISYM) = 0.0D0
         ELSE
            CALL CHO_GETSTOR_S(VCSTOR(ISYM),ISYM)
         END IF
      END DO

      END
      SUBROUTINE CHO_GETSTOR_S(VCSTOR,ISYM)
C
C     Purpose: get total vector storage (in words), symmetry ISYM.
C
      use ChoArr, only: nDimRS
      use ChoSwp, only: nnBstRSh, iiBstRSh, InfVec
#include "implicit.fh"
#include "cholesky.fh"

      IF (NUMCHO(ISYM) .LT. 1) THEN
         VCSTOR = 0.0D0
      ELSE
         IF (.NOT.Allocated(nDimRS)) Then
            IRED = INFVEC(NUMCHO(ISYM),2,ISYM)
            JRED = 3
            CALL CHO_GETRED(IRED,JRED,.FALSE.)
            CALL CHO_SETREDIND(IIBSTRSH,NNBSTRSH,NSYM,NNSHL,JRED)
            VCSTOR = DBLE(INFVEC(NUMCHO(ISYM),4,ISYM))
     &             + DBLE(NNBSTR(ISYM,JRED))
         ELSE
            IRED = INFVEC(NUMCHO(ISYM),2,ISYM)
            VCSTOR = DBLE(INFVEC(NUMCHO(ISYM),4,ISYM))
     &             + DBLE(NDIMRS(ISYM,IRED))
         END IF
      END IF

      END

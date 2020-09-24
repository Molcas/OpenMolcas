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
      SUBROUTINE CHO_SETPASS(DIAG,DIASH,ISYSH,IRED,CONV,NPOTSH)
C
C     Purpose: Check convergence and, if not converged, set up
C              integral pass.
C
#include "implicit.fh"
      DIMENSION DIAG(*), DIASH(*)
      INTEGER   ISYSH(*)
      LOGICAL   CONV
#include "cholesky.fh"

#if defined (_DEBUGPRINT_)
      CALL QENTER('_SETPASS')
#endif

C     Initialize the potential number of shell pairs that can
C     contribute.
C     -------------------------------------------------------

      NPOTSH = 0

C     Find max. abs. diagonal in each symmetry and the global max.
C     ------------------------------------------------------------

      DGMAX = -1.0D15
      CALL CHO_MAXABSDIAG(DIAG,IRED,DGMAX)

C     If not converged, set next integral pass.
C     -----------------------------------------

      CONV = DGMAX .LT. THRCOM
      IF (.NOT. CONV) THEN
         CALL CHO_SETMAXSHL(DIAG,DIASH,ISYSH,IRED)
         DO ISYM = 1,NSYM
            DIAMIN(ISYM) = MAX(DIAMAX(ISYM)*SPAN,THRCOM)
         END DO
         DO ISHLAB = 1,NNSHL
            IF (DIASH(ISHLAB) .GT. THRCOM) THEN
               NPOTSH = NPOTSH + 1
            ELSE
               DIASH(ISHLAB) = 0.0D0
            END IF
         END DO
      END IF

#if defined (_DEBUGPRINT_)
      CALL QEXIT('_SETPASS')
#endif

      END

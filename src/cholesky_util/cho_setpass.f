!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      SUBROUTINE CHO_SETPASS(DIAG,DIASH,ISYSH,IRED,CONV,NPOTSH)
!
!     Purpose: Check convergence and, if not converged, set up
!              integral pass.
!
      Implicit Real*8 (a-h,o-z)
      DIMENSION DIAG(*), DIASH(*)
      INTEGER   ISYSH(*)
      LOGICAL   CONV
#include "cholesky.fh"

!     Initialize the potential number of shell pairs that can
!     contribute.
!     -------------------------------------------------------

      NPOTSH = 0

!     Find max. abs. diagonal in each symmetry and the global max.
!     ------------------------------------------------------------

      DGMAX = -1.0D15
      CALL CHO_MAXABSDIAG(DIAG,IRED,DGMAX)

!     If not converged, set next integral pass.
!     -----------------------------------------

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

      END

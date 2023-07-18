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
      SUBROUTINE CHO_INIT1()
!
!     Purpose: initialize counter arrays.
!
      use ChoSwp, only: InfRed, InfVec
      Implicit Real*8 (a-h,o-z)
#include "cholesky.fh"
#include "choglob.fh"
#include "cho_para_info.fh"

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CHO_INIT1')

      INTEGER  CHO_ISUMELM
      EXTERNAL CHO_ISUMELM

      IF (RSTCHO) THEN

!        Read restart information.
!        -------------------------

         CALL CHO_GETRSTC()
         NUMCHT = CHO_ISUMELM(NUMCHO,NSYM)

      ELSE

!        Initialize vector info and counters.
!        ------------------------------------

         CALL IZERO(INFVEC,SIZE(INFVEC))
         CALL IZERO(NUMCHO,NSYM)
         NUMCHT = 0

!        Initialize reduced set info.
!        ----------------------------

         CALL IZERO(INFRED,SIZE(INFRED))

!        Initialize global integral pass counter.
!        ----------------------------------------

         XNPASS = 0

      END IF

!     Parallel init.
!     --------------

      IF (Cho_Real_Par) CALL IZERO(MYNUMCHO,NSYM)

      END

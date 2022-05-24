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
      SUBROUTINE CHO_RSCOPY(IRS1,IRS2)
C
C     Purpose: copy red. set info from location IRS1 to IRS2.
C              Special action is taken with INDRED if IRS1=1 so that it
C              will point as expected for the "current" reduced set.
C
      use ChoSwp, only: IndRed, iiBstRSh, nnBstRSh
      IMPLICIT NONE
      INTEGER IRS1, IRS2
      INTEGER IAB
      INTEGER MSYM
#include "cholesky.fh"

      MSYM=SIZE(iiBstRSh,1)
      nnBstRSh(:,:,IRS2) = nnBstRSh(:,:,IRS1)
      iiBstRSh(:,:,IRS2) = iiBstRSh(:,:,IRS1)
      iiBstR    (1:MSYM,IRS2) = iiBstR    (1:MSYM,IRS1)
      nnBstR    (1:MSYM,IRS2) = nnBstR    (1:MSYM,IRS1)
      IF (IRS1 .EQ. 1) THEN
         DO IAB = 1,SIZE(INDRED,1)
            INDRED(IAB,IRS2) = IAB
         END DO
      ELSE
         IndRed(:,iRS2) = IndRed(:,iRS1)
      END IF
      NNBSTRT(IRS2) = NNBSTRT(IRS1)

      END

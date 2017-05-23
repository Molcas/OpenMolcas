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
      SUBROUTINE CHO_RSCOPY(IIBSTRSH,NNBSTRSH,INDRED,IRS1,IRS2,
     &                      MSYM,MMSHL,LMMBSTRT,NRS)
C
C     Purpose: copy red. set info from location IRS1 to IRS2.
C              Special action is taken with INDRED if IRS1=1 so that it
C              will point as expected for the "current" reduced set.
C
      IMPLICIT NONE
      INTEGER IRS1, IRS2, MSYM, MMSHL, LMMBSTRT, NRS
      INTEGER IIBSTRSH(MSYM,MMSHL,NRS), NNBSTRSH(MSYM,MMSHL,NRS)
      INTEGER INDRED(LMMBSTRT,NRS)
      INTEGER IAB
#include "cholesky.fh"

      CALL ICOPY(MSYM*MMSHL,IIBSTRSH(1,1,IRS1),1,IIBSTRSH(1,1,IRS2),1)
      CALL ICOPY(MSYM*MMSHL,NNBSTRSH(1,1,IRS1),1,NNBSTRSH(1,1,IRS2),1)
      CALL ICOPY(MSYM,IIBSTR(1,IRS1),1,IIBSTR(1,IRS2),1)
      CALL ICOPY(MSYM,NNBSTR(1,IRS1),1,NNBSTR(1,IRS2),1)
      IF (IRS1 .EQ. 1) THEN
         DO IAB = 1,MMBSTRT
            INDRED(IAB,IRS2) = IAB
         END DO
      ELSE
         CALL ICOPY(MMBSTRT,INDRED(1,IRS1),1,INDRED(1,IRS2),1)
      END IF
      NNBSTRT(IRS2) = NNBSTRT(IRS1)

      END

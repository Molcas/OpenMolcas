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
      INTEGER FUNCTION CHO_ISAO(IAO)
C
C     Purpose: return symmetry of AO number IAO (in global list).
C
#include "implicit.fh"
#include "cholesky.fh"
#include "choorb.fh"

      INTEGER  CHO_IRANGE
      EXTERNAL CHO_IRANGE

#if defined (_DEBUGPRINT_)
      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'CHO_ISAO')

      IF ((IAO.GT.NBAST) .OR. (IAO.LT.1)) THEN
         WRITE(LUPRI,'(//,1X,A,A,I10)')
     &   SECNAM,': AO index out of bounds: ',IAO
         WRITE(LUPRI,'(A,I10,A,/)')
     &   'Maximum possible: NBAST = ',NBAST,'(from common block)'
         IF (NBAST .LT. 1) THEN
            CALL CHO_QUIT('Initialization error detected in '//SECNAM,
     &                    102)
         ELSE
            CALL CHO_QUIT('Internal error detected in '//SECNAM,
     &                    103)
         END IF
      END IF
#endif

      CHO_ISAO = CHO_IRANGE(IAO,IBAS,NSYM,.FALSE.)

      END

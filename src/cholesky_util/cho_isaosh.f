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
      INTEGER FUNCTION CHO_ISAOSH(IAO,ISHL)
C
C     Purpose: return symmetry of AO number IAO in shell ISHL.
C
#if defined (_DEBUGPRINT_)
      use ChoArr, only: nBstSh
#endif
#include "implicit.fh"
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      INTEGER  CHO_IRANGE
      EXTERNAL CHO_IRANGE

#if defined (_DEBUGPRINT_)
      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_ISAOSH')

      IF ((ISHL.GT.NSHELL) .OR. (ISHL.LT.1)) THEN
         WRITE(LUPRI,'(//,1X,A,A,I10)')
     &   SECNAM,': shell index out of bounds: ',ISHL
         WRITE(LUPRI,'(A,I10,A,/)')
     &   'Maximum possible: NSHELL = ',NSHELL,'(from common block)'
         IF (NSHELL .LT. 1) THEN
            CALL CHO_QUIT('Initialization error detected in '//SECNAM,
     &                    102)
         ELSE
            CALL CHO_QUIT('Internal error detected in '//SECNAM,
     &                    103)
         END IF
      ELSE IF ((IAO.GT.NBSTSH(ISHL)) .OR. (IAO.LT.1)) THEN
         WRITE(LUPRI,'(//,1X,A,A,I10)')
     &   SECNAM,': AO index out of bounds: ',IAO,' shell: ',ISHL
         WRITE(LUPRI,'(A,I10,A,/)')
     &   'Maximum possible: NBSTSH(ISHL) = ',NBSTSH(ISHL),
     &   '(from common block)'
         IF (NBSTSH(ISHL) .LT. 1) THEN
            CALL CHO_QUIT('Initialization error detected in '//SECNAM,
     &                    102)
         ELSE
            CALL CHO_QUIT('Internal error detected in '//SECNAM,
     &                    103)
         END IF
      END IF
#endif

      KIBASSH = ip_IBASSH + NSYM*(ISHL-1)
      CHO_ISAOSH = CHO_IRANGE(IAO,IWORK(KIBASSH),NSYM,.FALSE.)

      END

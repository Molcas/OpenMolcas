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
      SUBROUTINE CHO_INTCHK_ID_OF(LABEL,ID,IOPT)
C
C     Purpose: for minimal integral checking,
C              IOPT=-1 : return label corresponding to id ID.
C              else    : return index of shell quadruple corresponding
C                        to check label LABEL.
C
      IMPLICIT NONE
      CHARACTER*8 LABEL
      INTEGER ID, IOPT
#include "cholesky.fh"

      INTEGER  CHO_TABIND
      EXTERNAL CHO_TABIND

      INTEGER     NTABLE
      PARAMETER   (NTABLE = 12)
      CHARACTER*8 TABLE(NTABLE)
      DATA TABLE  /'EXCL RS1','MAX|XRS1','MIN|XRS1',
     &             'NEG DIAG','MAX|NEG ','MIN|NEG ',
     &             'NEG->ZER','MAX|NEGZ','MIN|NEGZ',
     &             'MAX DIAG','MIN DIAG','MAX|MIN '/

      IF (IOPT .EQ. -1) THEN
         IF (ID.LT.1 .OR. ID.GT.NTABLE) THEN
            LABEL = 'UNKNOWN '
         ELSE
            LABEL = TABLE(ID)
         END IF
      ELSE
         ID = CHO_TABIND(TABLE,8,NTABLE,' ',0,0,LABEL)
      END IF

      END

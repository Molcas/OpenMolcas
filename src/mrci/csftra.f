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
      SUBROUTINE CSFTRA(KEY,CI,AREF)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*4 KEY
      DIMENSION CI(NCONF),AREF(NREF,NREF)

#include "SysDef.fh"

#include "mrci.fh"
      DIMENSION TMP(MXREF)
      IF(NREF.EQ.1) RETURN
      IF(KEY.EQ.' CSF') THEN
        DO 20 I=1,NREF
          SUM=0.0D00
          DO 10 J=1,NREF
            SUM=SUM+AREF(I,J)*CI(IREFX(J))
10        CONTINUE
          TMP(I)=SUM
20      CONTINUE
      ELSE
        DO 50 I=1,NREF
          SUM=0.0D00
          DO 40 J=1,NREF
            SUM=SUM+AREF(J,I)*CI(IREFX(J))
40        CONTINUE
          TMP(I)=SUM
50      CONTINUE
      END IF
      DO 60 I=1,NREF
        CI(IREFX(I))=TMP(I)
60    CONTINUE
      RETURN
      END

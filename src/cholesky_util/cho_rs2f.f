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
      INTEGER FUNCTION CHO_RS2F(LAB,ISHLAB,ISYMAB,IRED)
C
C     Purpose: return index in reduced set IRED (1,2,3) of
C              element LAB in shell pair ISHLAB (sym. ISYMAB).
C              If not included in this reduced set, 0 is returned.
C
#include "implicit.fh"
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'CHO_RS2F')

      INTEGER K, K2

      IIBSTRSH(I,J,K)=IWORK(ip_IIBSTRSH-1+NSYM*NNSHL*(K-1)+NSYM*(J-1)+I)
      NNBSTRSH(I,J,K)=IWORK(ip_NNBSTRSH-1+NSYM*NNSHL*(K-1)+NSYM*(J-1)+I)
      INDRED(I,J)=IWORK(ip_INDRED-1+MMBSTRT*(J-1)+I)

      CHO_RS2F = 0

      K  = IIBSTR(ISYMAB,IRED) + IIBSTRSH(ISYMAB,ISHLAB,IRED)
      K2 = K + NNBSTRSH(ISYMAB,ISHLAB,IRED)
      IF (IRED .EQ. 1) THEN
         DO WHILE (K.LT.K2 .AND. CHO_RS2F.EQ.0)
            K = K + 1
            IF (INDRED(K,1) .EQ. LAB) THEN
               CHO_RS2F = K
            END IF
         END DO
      ELSE IF (IRED.EQ.2 .OR. IRED.EQ.3) THEN
         DO WHILE (K.LT.K2 .AND. CHO_RS2F.EQ.0)
            K = K + 1
            IF (INDRED(INDRED(K,IRED),1) .EQ. LAB) THEN
               CHO_RS2F = K
            END IF
         END DO
      ELSE
         CALL CHO_QUIT('IRED error in '//SECNAM,104)
      END IF

      END

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
      INTEGER FUNCTION CHO_LREAD(ISYM,LWRK)
C
C     Purpose: return a reasonable scratch space dimension for reading
C              previous vectors using cho_getvec.
C
#include "implicit.fh"
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      INTEGER MNVECRS1
      PARAMETER (MNVECRS1 = 5)

      PARAMETER (N2 = INFVEC_N2)

      INFVEC(I,J,K)=IWORK(ip_INFVEC-1+MAXVEC*N2*(K-1)+MAXVEC*(J-1)+I)

      IF (CHO_IOVEC .EQ. 1) THEN
         IF (NVECRS1(ISYM).LT.1 .AND. NUMCHO(ISYM).GT.0) THEN
            NVECRS1(ISYM) = 1
            JVEC = 1
            IRED = INFVEC(JVEC,2,ISYM)
            DO WHILE (JVEC .LT. NUMCHO(ISYM))
               JVEC = JVEC + 1
               JRED = INFVEC(JVEC,2,ISYM)
               IF (JRED .EQ. IRED) THEN
                  NVECRS1(ISYM) = NVECRS1(ISYM) + 1
               ELSE
                  JVEC = NUMCHO(ISYM)
               END IF
            END DO
         END IF
         LEN1 = LWRK/3 - 1
         LEN2 = MAX(NVECRS1(ISYM),MNVECRS1)*NNBSTR(ISYM,1)
         LEN3 = MIN(LEN1,LEN2)
         LMIN = 2*NNBSTR(ISYM,1)
         CHO_LREAD = MAX(LEN3,LMIN) + 1
      ELSE IF (CHO_IOVEC.EQ.2 .OR. CHO_IOVEC.EQ.3 .OR. CHO_IOVEC.EQ.4)
     & THEN
         LEN1 = LWRK/3 - 1
         LMIN = 2*NNBSTR(ISYM,1)
         CHO_LREAD = MAX(LEN1,LMIN) + 1
      ELSE
         CHO_LREAD = 2*NNBSTR(ISYM,1)
      END IF

      END

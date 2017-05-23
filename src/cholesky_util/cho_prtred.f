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
      SUBROUTINE CHO_PRTRED(IOPT)
C
C     Purpose: print information about reduced set.
C
#include "implicit.fh"
#include "choorb.fh"
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      DIMENSION XBAS(8), XXBAS(8)

      INTEGER NSHP(2)
      LOGICAL CONTRIB(2)

      MULD2H(I,J)=IEOR(I-1,J-1)+1
      NNBSTRSH(I,J,K)=IWORK(ip_NNBSTRSH-1+NSYM*NNSHL*(K-1)+NSYM*(J-1)+I)

      DO ISYM = 1,NSYM
         XBAS(ISYM) = DBLE(NBAS(ISYM))
      END DO

      XXBAST = 0.0D0
      DO ISYM = 1,NSYM
         XXBAS(ISYM) = 0.0D0
         DO ISYMB = 1,NSYM
            ISYMA = MULD2H(ISYMB,ISYM)
            IF (ISYMA .EQ. ISYMB) THEN
               XXBAS(ISYM) = XXBAS(ISYM)
     &                     + XBAS(ISYMA)*(XBAS(ISYMA) + 1.0D0)/2.0D0
            ELSE IF (ISYMA .GT. ISYMB) THEN
               XXBAS(ISYM) = XXBAS(ISYM) + XBAS(ISYMA)*XBAS(ISYMB)
            END IF
         END DO
         XXBAST = XXBAST + XXBAS(ISYM)
      END DO

      IF (IOPT .EQ. 1) THEN
         NRED = 1
      ELSE
         NRED = 2
      END IF

      DO IRED = 1,NRED
         NSHP(IRED) = 0
         DO ISHLAB = 1,NNSHL
            CONTRIB(IRED) = .FALSE.
            ISYM = 0
            DO WHILE ((.NOT.CONTRIB(IRED)) .AND. ISYM.LT.NSYM)
               ISYM = ISYM + 1
               IF (NNBSTRSH(ISYM,ISHLAB,IRED) .GT. 0) THEN
                  CONTRIB(IRED) = .TRUE.
               END IF
            END DO
            IF (CONTRIB(IRED)) THEN
               NSHP(IRED) = NSHP(IRED) + 1
            END IF
         END DO
      END DO

      CALL CHO_HEAD('Reduced Set Information','=',80,LUPRI)

      IF (NNSHL_TOT .EQ. 0) THEN
         PCT1 = 9.9D9
      ELSE
         PCT1 = 1.0D2*DBLE(NSHP(1))/DBLE(NNSHL_TOT)
      END IF
      IF (IOPT .EQ. 1) THEN
         WRITE(LUPRI,'(/,A,/,A)')
     &   'Sym.          Full   First Red. Set',
     &   '-----------------------------------'
         DO ISYM = 1,NSYM
            WRITE(LUPRI,'(I3,3X,F12.1,7X,I10)')
     &      ISYM,XXBAS(ISYM),NNBSTR(ISYM,1)
         END DO
         WRITE(LUPRI,'(A)')
     &   '-----------------------------------'
         WRITE(LUPRI,'(A,F12.1,7X,I10)')
     &   'Total:',XXBAST,NNBSTRT(1)
         WRITE(LUPRI,'(A)')
     &   '-----------------------------------'
         WRITE(LUPRI,'(/,A,I10,A,I10,A,F7.2,A)')
     &   'First Reduced Set:',NSHP(1),' of',NNSHL_TOT,
     &   ' shell pairs contribute (',PCT1,'%)'
      ELSE
         IF (NNSHL_TOT .EQ. 0) THEN
            PCT2 = 9.9D9
         ELSE
            PCT2 = 1.0D2*DBLE(NSHP(2))/DBLE(NNSHL_TOT)
         END IF
         WRITE(LUPRI,'(/,A,/,A,/,A)')
     &   '                          Reduced Set',
     &   'Sym.          Full      First    Current',
     &   '----------------------------------------'
         DO ISYM = 1,NSYM
            WRITE(LUPRI,'(I3,3X,F12.1,1X,I10,1X,I10)')
     &      ISYM,XXBAS(ISYM),NNBSTR(ISYM,1),NNBSTR(ISYM,2)
         END DO
         WRITE(LUPRI,'(A)')
     &   '----------------------------------------'
         WRITE(LUPRI,'(A,F12.1,1X,I10,1X,I10)')
     &   'Total:',XXBAST,NNBSTRT(1),NNBSTRT(2)
         WRITE(LUPRI,'(A)')
     &   '----------------------------------------'
         WRITE(LUPRI,'(/,A,I10,A,I10,A,F7.2,A)')
     &   'First Reduced Set:',NSHP(1),' of',NNSHL_TOT,
     &   ' shell pairs contribute (',PCT1,'%)'
         WRITE(LUPRI,'(A,I10,A,I10,A,F7.2,A)')
     &   'Curr. Reduced Set:',NSHP(2),' of',NNSHL_TOT,
     &   ' shell pairs contribute (',PCT2,'%)'
      END IF

      END

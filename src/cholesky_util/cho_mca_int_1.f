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
      SUBROUTINE CHO_MCA_INT_1(IJ,KL,XINT,LINT,LOCPRT)
C
C     Purpose: calculate shell-quadruple (IJ|KL) and return
C              them in XINT.
C
C     Notes:
C        LOCPRT: flag for printing the shell quadruple to output;
C                output format differs depending on IFCSEW.
C
      use ChoArr, only: nBstSh, iSP2F
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL Integral_WrOut_Cho
      REAL*8   XINT(LINT)
      LOGICAL  LOCPRT
#include "itmax.fh"
#include "cholesky.fh"

      CHARACTER*13 SECNAM
      PARAMETER (SECNAM = 'CHO_MCA_INT_1')

      ITRI(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J

C     Initializations.
C     ----------------

      CALL CHO_INVPCK(ISP2F(IJ),I,J,.TRUE.)
      CALL CHO_INVPCK(ISP2F(KL),K,L,.TRUE.)

      SHCD = IJ
      SHAB = KL
      SHC = I
      SHD = J
      SHA = K
      SHB = L

C     Calculate integrals.
C     --------------------

      CALL EVAL_IJKL(I,J,K,L,XINT,LINT,Integral_WrOut_Cho)

C     Print integrals.
C     ----------------

      IF (LOCPRT) THEN

         IF (IFCSEW.EQ.1) THEN

            WRITE(LUPRI,'(//,5X,A,A,4I5,A)')
     &      SECNAM,': shell quadruple ',I,J,K,L,':'

            NUMI = NBSTSH(I)
            NUMJ = NBSTSH(J)
            NUMK = NBSTSH(K)
            NUML = NBSTSH(L)
            IF (I .EQ. J) THEN
               NUMIJ = NUMI*(NUMJ + 1)/2
            ELSE
               NUMIJ = NUMI*NUMJ
            END IF

            IF (K .EQ. L) THEN
               DO LL = 1,NUML
                  DO KK = 1,LL
                     KKLL = ITRI(KK,LL)
                     IF (I .EQ. J) THEN
                        DO JJ = 1,NUMJ
                           DO II = 1,JJ
                              IIJJ = ITRI(II,JJ)
                              KOFF = NUMIJ*(KKLL - 1) + IIJJ
                              WRITE(LUPRI,*)
     &                        '(',I,J,K,L,') [',II,JJ,KK,LL,'] = ',
     &                        XINT(KOFF)
                           END DO
                        END DO
                     ELSE
                        DO JJ = 1,NUMJ
                           DO II = 1,NUMI
                              IIJJ = NUMI*(JJ - 1) + II
                              KOFF = NUMIJ*(KKLL - 1) + IIJJ
                              WRITE(LUPRI,*)
     &                        '(',I,J,K,L,') [',II,JJ,KK,LL,'] = ',
     &                        XINT(KOFF)
                           END DO
                        END DO
                     END IF
                  END DO
               END DO
            ELSE
               DO LL = 1,NUML
                  DO KK = 1,NUMK
                     KKLL = NUMK*(LL - 1) + KK
                     IF (I .EQ. J) THEN
                        DO JJ = 1,NUMJ
                           DO II = 1,JJ
                              IIJJ = ITRI(II,JJ)
                              KOFF = NUMIJ*(KKLL - 1) + IIJJ
                              WRITE(LUPRI,*)
     &                        '(',I,J,K,L,') [',II,JJ,KK,LL,'] = ',
     &                        XINT(KOFF)
                           END DO
                        END DO
                     ELSE
                        DO JJ = 1,NUMJ
                           DO II = 1,NUMI
                              IIJJ = NUMI*(JJ - 1) + II
                              KOFF = NUMIJ*(KKLL - 1) + IIJJ
                              WRITE(LUPRI,*)
     &                        '(',I,J,K,L,') [',II,JJ,KK,LL,'] = ',
     &                        XINT(KOFF)
                           END DO
                        END DO
                     END IF
                  END DO
               END DO
            END IF

         ELSE IF (IFCSEW.EQ.2 .OR. IFCSEW.EQ.3) THEN

            CALL CHO_PRTINT(IJ,KL,XINT,LINT)

         ELSE

            WRITE(LUPRI,*) SECNAM,': IFCSEW=',IFCSEW
            CALL CHO_QUIT(SECNAM//': IFCSEW out of bounds!',103)

         END IF

      END IF

      END

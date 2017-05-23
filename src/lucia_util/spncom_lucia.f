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
      SUBROUTINE SPNCOM_LUCIA(  NOPEN,    MS2,   NDET, IABDET, IABUPP,
     &                          IFLAG, PSSIGN, IPRCSF)
*
* Combinations of nopen unpaired electrons.Required
* spin projection MS2/2.
*
#include "implicit.fh"
*. MXPDIM is included to have access to MXPORB
#include "mxpdim.fh"
      INTEGER ADD
      DIMENSION IABDET(NOPEN,*),IABUPP(NOPEN,*)
*. Should have length of max number of open orbitals
      DIMENSION IWORK(MXPORB)
*
* LENGTH OF IWORK MUST BE AT LEAST NOPEN
*
      NTEST = 0
      NTEST = MAX(NTEST,IPRCSF)
      NDET=0
      NUPPER = 0
*
* Determinants are considered as binary numbers,1=alpha,0=beta
*
      MX=2 ** NOPEN
      CALL ISETVC(IWORK,0,NOPEN)
      IFIRST = 1
* Loop over all possible binary numbers
      DO 200 I=1, MX
        IF(IFIRST.EQ.1) THEN
*. Initial number
          IZERO = 0
          CALL ISETVC(IWORK,IZERO,NOPEN)
          IFIRST = 0
        ELSE
*. Next number
          ADD=1
          J=0
  190     CONTINUE
          J=J+1
          IF(IWORK(J).EQ.1) THEN
            IWORK(J)=0
          ELSE
            IWORK(J)=1
            ADD=0
          END IF
          IF( ADD .EQ. 1 ) GOTO 190
        END IF
C
C.. 2 :  CORRECT SPIN PROJECTION ?
        NALPHA=0
        DO J=1,NOPEN
          NALPHA=NALPHA+IWORK(J)
        END DO
C
        IF(2*NALPHA-NOPEN.EQ.MS2.AND.
     &    .NOT.(PSSIGN.NE.0.0D0 .AND. IWORK(1).EQ.0)) THEN
          IF (IFLAG .LT. 3 ) THEN
            NDET=NDET+1
            CALL ICOPVE(IWORK,IABDET(1,NDET),NOPEN)
          END IF
*
          IF (IFLAG .GT. 1 ) THEN
C UPPER DET ?
            MS2L = 0
            LUPPER = 1
C
            DO IEL = 1,NOPEN
              IF (IWORK(IEL).EQ.1) THEN
                 MS2L = MS2L + 1
              ELSE
                 MS2L = MS2L - 1
              END IF
              IF( MS2L .LT. 0 ) LUPPER = 0
            END DO
*
            IF( LUPPER .EQ. 1 ) THEN
              NUPPER = NUPPER + 1
              CALL ICOPVE(IWORK,IABUPP(1,NUPPER),NOPEN)
            END IF
          END IF
        END  IF
C
  200 CONTINUE
C
      XMSD2=0.5D0*DBLE(MS2)
C
      IF(NTEST.GE.5.AND.IFLAG .NE.3) THEN
         WRITE(6,1010) NOPEN,NDET,XMSD2
 1010    FORMAT(1H0,2X,I3,' Unpaired electrons give ',I5,/,
     +'           combinations with spin projection ',F12.7)
         WRITE(6,*)
         WRITE(6,'(A)') '  Combinations : '
         WRITE(6,'(A)') '  ============== '
         WRITE(6,*)
         DO 20 J=1,NDET
           WRITE(6,1020) J,(IABDET(K,J),K=1,NOPEN)
  20     CONTINUE
 1020    FORMAT(1H0,I5,2X,30I2,/,(1H ,7X,30I2))
      END IF
C
      IF( IFLAG.GT.1.AND.NTEST.GE.5) THEN
         WRITE(6,*)
         WRITE(6,'(A)') ' Upper determinants '
         WRITE(6,'(A)') ' ================== '
         WRITE(6,*)
         DO 22 J=1,NUPPER
           WRITE(6,1020) J,(IABUPP(K,J),K=1,NOPEN)
  22     CONTINUE
      END IF
C
      RETURN
      END

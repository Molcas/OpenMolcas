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
      SUBROUTINE CSDTVC_MCLR(CSFVEC,DETVEC,IWAY,DTOCMT,ICTSDT,
     &                  IREFSM,ICOPY,IPRNT)
*
* IWAY = 1 : CSF to DETERMINANT TRANSFORMATION
* IWAY = 2 : DETERMINANT TO CSF TRANSFORMATION
*
* ICOPY .NE. 0 : Copy output into input
*                so input becomes output while
*                output remains output
*
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CSFVEC(*),DETVEC(*)
      DIMENSION DTOCMT(*),ICTSDT(*)
#include "detdim.fh"
#include "spinfo_mclr.fh"
*

      ZERO = 0.0D0
      ONE = 1.0D0
      IOFFCS = 0 ! dummy initialize
      IOFFDT = 0 ! dummy initialize
      IOFFCD = 0 ! dummy initialize
*

      NDET = NDTASM(IREFSM)
      NCSF = NCSASM(IREFSM)
*
      IF(IWAY .EQ. 1 ) THEN
*
* ===========================
*.. CSF to DET transformation
* ===========================
*
        CALL SETVEC(DETVEC,ZERO ,NDET)
*. Multiply with  expansion matrix
        DO 100 ITYP = 1,NTYP
          call xflush(6)
          IDET = NDPCNT(ITYP)
          ICSF = NCPCNT(ITYP)
          ICNF = NCNATS(ITYP,IREFSM)
          IF(ITYP .EQ. 1 ) THEN
            IOFFCS = 1
            IOFFDT = 1
            IOFFCD = 1
          ELSE
            IOFFCS = IOFFCS+NCNATS(ITYP-1,IREFSM)*NCPCNT(ITYP-1)
            IOFFDT = IOFFDT+NCNATS(ITYP-1,IREFSM)*NDPCNT(ITYP-1)
            IOFFCD = IOFFCD + NDPCNT(ITYP-1)*NCPCNT(ITYP-1)
          END IF
          IF( IDET*ICNF*ICSF .GT. 0 ) THEN
          CALL  DGEMM_('N','N',IDET,ICNF,ICSF,ONE ,
     &                DTOCMT(IOFFCD),IDET,
     &                CSFVEC(IOFFCS),ICSF,ZERO,
     &                DETVEC(IOFFDT),IDET)
         END IF
  100   CONTINUE
*. Sign changes
        CALL COPVEC(DETVEC,CSFVEC,NDET)
*. Change to string ordering
        CALL SCAVCS(DETVEC,CSFVEC,ICTSDT,NDET)
        IF(ICOPY.NE.0) CALL COPVEC(DETVEC,CSFVEC,NDET)
      ELSE
*
* ====================================
*  Determinant to csf transformation
* ====================================
*
C. To CSF ordering

        CALL GATVCS(CSFVEC,DETVEC,ICTSDT,NDET)
        CALL COPVEC(CSFVEC,DETVEC,NDET)
C. Multiply with CIND expansion matrix
        DO 200 ITYP = 1,NTYP
          IDET = NDPCNT(ITYP)
          ICSF = NCPCNT(ITYP)
          ICNF = NCNATS(ITYP,IREFSM)
          IF(ITYP .EQ. 1 ) THEN
            IOFFCS = 1
            IOFFDT = 1
            IOFFCD = 1
          ELSE
            IOFFCS = IOFFCS+NCNATS(ITYP-1,IREFSM)*NCPCNT(ITYP-1)
            IOFFDT = IOFFDT+NCNATS(ITYP-1,IREFSM)*NDPCNT(ITYP-1)
            IOFFCD = IOFFCD + NDPCNT(ITYP-1)*NCPCNT(ITYP-1)
          END IF
          IF( IDET*ICNF*ICSF .GT. 0 ) THEN
          CALL  DGEMM_('T','N',ICSF,ICNF,IDET,ONE ,
     &                DTOCMT(IOFFCD),IDET,
     &                DETVEC(IOFFDT),IDET,ZERO,
     &                CSFVEC(IOFFCS),ICSF)
          END IF
  200   CONTINUE
        IF( ICOPY .NE. 0 ) CALL COPVEC(CSFVEC,DETVEC,NCSF)
       END IF

      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(IPRNT)
      END

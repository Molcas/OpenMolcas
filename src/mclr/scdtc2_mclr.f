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
      SUBROUTINE SCDTC2_MCLR(RASVEC,ISMOST,ICBLTP,NSMST,NOCTPA,NOCTPB,
     &                  NSASO,NSBSO,IOCOC,IDC,IWAY,IMMLST,IPRNT)
*
* Scale elements of a RAS vector to transfer between
* combinations and packed determinants
* IWAY = 1 : dets to combs
* IWAY = 2 : combs to dets
* Combination storage mode is defined BY IDC
*
* General symmetry version , Feb 1991
*
      IMPLICIT real*8(A-H,O-Z)
      DIMENSION RASVEC(*),NSASO(NOCTPA,*),NSBSO(NOCTPB,*)
      DIMENSION IOCOC(NOCTPA,NOCTPB)
      DIMENSION ISMOST(*),ICBLTP(*),IMMLST(*)
*

      NTEST = 0000
      NTEST = MAX(NTEST,IPRNT)
      IF( NTEST .GT. 10 ) THEN
        WRITE(6,*) ' Information from SCDTC2 '
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' Input vector '
        CALL WRTRS2_MCLR(RASVEC,ISMOST,ICBLTP,IOCOC,
     &              NOCTPA,NOCTPB,NSASO,NSBSO,NSMST)
      END IF
*
      SQ2 = SQRT(2.0D0)
      SQ2I = 1.0D0/SQ2
*
      IBASE = 1
      DO 200 IASM = 1, NSMST
*
        IBSM = ISMOST(IASM)
        IF(IBSM.EQ.0.OR.ICBLTP(IASM).EQ.0) GOTO 200
        DO  100 IATP = 1, NOCTPA
          IF(ICBLTP(IASM).EQ.2) THEN
            IBTPMX = IATP
          ELSE
            IBTPMX = NOCTPB
          END IF
          NIA   = NSASO(IATP,IASM)
          DO 50 IBTP = 1,IBTPMX
            IF(IOCOC(IATP,IBTP).EQ.0) GOTO   50
*. Number of elements in this block
          call xflush(6)
            NIB = NSBSO(IBTP,IBSM)
            IF(ICBLTP(IASM).EQ.2.AND.IATP.EQ.IBTP) THEN
                NELMNT =  NIA*(NIA+1)/2
            ELSE
                NELMNT =  NIA*NIB
            END IF

          IF(IDC.EQ.2) THEN
            IF(IWAY.EQ.1) THEN
              FACTOR = SQ2
            ELSE
              FACTOR = SQ2I
            END IF
            CALL DSCAL_(NELMNT,FACTOR,RASVEC(IBASE),1)
            IF(IASM.EQ.IBSM.AND.IATP.EQ.IBTP) THEN
              FACTOR = 1.0D0/FACTOR
              CALL SCLDIA(RASVEC(IBASE),FACTOR,NIA,1)
            END IF
          ELSE IF(IDC.EQ.3.AND.IMMLST(IASM).NE.IASM) THEN
            IF(IWAY.EQ.1) THEN
              FACTOR = SQ2
            ELSE
              FACTOR = SQ2I
            END IF
            CALL DSCAL_(NELMNT,FACTOR,RASVEC(IBASE),1)
*Ml Ms combinations
          call xflush(6)
          ELSE IF(IDC.EQ.4) THEN
            IF(IWAY.EQ.1) THEN
              IF(IASM.EQ.IBSM) THEN
                FACTOR = SQ2
              ELSE
                FACTOR = 2.0D0
              END IF
            ELSE IF(IWAY.EQ.2) THEN
              IF(IASM.EQ.IBSM) THEN
                FACTOR = SQ2I
              ELSE
                FACTOR = 0.5D0
              END IF
            END IF
            CALL DSCAL_(NELMNT,FACTOR,RASVEC(IBASE),1)
            IF(IATP.EQ.IBTP) THEN
              IF(IWAY.EQ.1) THEN
                FACTOR = SQ2I
              ELSE IF(IWAY.EQ.2) THEN
                FACTOR = SQ2
              END IF
              CALL SCLDIA(RASVEC(IBASE),FACTOR,NIA,1)
            END IF
          END IF
*
          IBASE = IBASE + NELMNT
  50      CONTINUE
 100    CONTINUE
 200  CONTINUE

*
      IF( NTEST .GT. 10 ) THEN
        WRITE(6,*) ' Scaled vector '
        call xflush(6)
        CALL WRTRS2_MCLR(RASVEC,ISMOST,ICBLTP,IOCOC,
     &              NOCTPA,NOCTPB,NSASO,NSBSO,NSMST)
      END IF
*
      RETURN
      END


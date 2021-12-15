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
      SUBROUTINE NGASDT(  IOCCMN,  IOCCMX,    NGAS,  ITOTSM,   NSMST,
     &                    NOCTPA,  NOCTPB,   NSSOA,   NSSOB,   IAOCC,
     &                     IBOCC, MXPNGAS,     NCOMB,  XNCOMB,
     &                      MXSB,  MXSOOB,   IBLTP,  NTTSBL,    LCOL,
     &                     IOCOC,MXSOOB_AS)
*
*
* Number of combimations with symmetry ITOTSM and
* occupation between IOCCMN and IOCCMX
*
* In view of the limited range of I*4, the number of dets
* is returned as integer and  real*8
*
* MXSB is largest UNPACKED symmetry block
* MXSOOB is largest UNPACKED symmetry-type-type block
* NTTSBL is number of TTS blocks in vector
* LCOL is the sum of the number of columns in each block
*
*
* Winter 94/95
* May 1999 : Loops restructrured to sym,type,type (leftmost innerst)
*            MXSB not calculated
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Allowed combinations of alpha and beta types
      INTEGER IOCOC(NOCTPA,NOCTPB)
*. Occupation constraints
      DIMENSION IOCCMN(NGAS),IOCCMX(NGAS)
*. Occupation of alpha and beta strings
      DIMENSION IAOCC(MXPNGAS,*),IBOCC(MXPNGAS,*)
*. Number of strings per supergroup and symmetry
      DIMENSION NSSOA(NSMST,*),NSSOB(NSMST,*)
*. block types
      DIMENSION IBLTP(*)
*
      NTEST = 0
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' NGASDT speaking'
        WRITE(6,*) ' ==============='
        WRITE(6,*) ' NGAS NOCTPA,NOCTPB ',NGAS,NOCTPA,NOCTPB
        WRITE(6,*) ' ITOTSM ', ITOTSM
        WRITE(6,*) ' Upper and lower occupation constraints'
        CALL IWRTMA(IOCCMN,1,NGAS,1,NGAS)
        CALL IWRTMA(IOCCMX,1,NGAS,1,NGAS)
        WRITE(6,*) ' IOCOC matrix '
        CALL IWRTMA(IOCOC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
        WRITE(6,*) ' Number of alpha and beta strings '
        CALL IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
        CALL IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
      END IF
*
      MXSB = 0
      MXSOOB = 0
      MXSOOB_AS = 0
      NCOMB = 0
      XNCOMB = 0.0D0
      NTTSBL = 0
      LCOL = 0
*
      DO 200 IATP = 1, NOCTPA
        DO 100 IBTP = 1, NOCTPB
*
          IF(NTEST.GE.10) THEN
            WRITE(6,*) ' Alpha super group and beta super group'
            CALL IWRTMA(IAOCC(1,IATP),1,NGAS,1,NGAS)
            CALL IWRTMA(IBOCC(1,IBTP),1,NGAS,1,NGAS)
          END IF
*
          IF(IOCOC(IATP,IBTP).EQ.1) THEN
*
            LTTS_AS = 0
            DO 300 IASM = 1, NSMST
              IF(IBLTP(IASM).EQ.0) GOTO 300
              CALL SYMCOM(2,1,IASM,IBSM,ITOTSM)
              IF(IBSM.NE.0) THEN
                IF(IBLTP(IASM).EQ.2) THEN
                  ISYM = 1
                ELSE
                  ISYM = 0
                END IF
                IF(ISYM.EQ.1.AND.IBTP.GT.IATP) GOTO 300
                LASTR = NSSOA(IASM,IATP)
                LBSTR = NSSOB(IBSM,IBTP)
*. Size of unpacked block
                LTTSUP =  LASTR*LBSTR
*. Size of packed block
                IF(ISYM.EQ.0.OR.IATP.NE.IBTP) THEN
                  LTTSBL = LASTR*LBSTR
                  XNCOMB = XNCOMB + dble(LASTR)*dble(LBSTR)
                ELSE
                  LTTSBL = LASTR*(LASTR+1)/2
                  XNCOMB = XNCOMB + 0.5D0*dble(LASTR+1)*dble(LASTR)
                END IF
                LTTS_AS = LTTS_AS + LTTSUP
                NCOMB = NCOMB + LTTSBL
                MXSOOB = MAX(MXSOOB,LTTSUP)
                NTTSBL = NTTSBL + 1
                LCOL = LCOL + NSSOB(IBSM,IBTP)
              END IF
  300       CONTINUE
            MXSOOB_AS = MAX(MXSOOB_AS,LTTS_AS)
          END IF
  100   CONTINUE
  200 CONTINUE
*
      IF(NTEST.GE.1) THEN
        WRITE(6,*) ' NGASDT : NCOMB XNCOMB ,NTTSBL',
     &               NCOMB,XNCOMB,NTTSBL
      END IF
*
*
      RETURN
      END

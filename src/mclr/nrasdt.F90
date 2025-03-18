!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1984,1989-1993, Jeppe Olsen                            *
!***********************************************************************
      SUBROUTINE NRASDT(MNRS1,MXRS1,MNRS3,MXRS3,ITOTSM,                 &
     &                  NSMST,NOCTPA,NOCTPB,IEL1A,IEL1B,                &
     &                  NSSOA,NSSOB,                                    &
     &                  IEL3A,IEL3B,NCOMB,XNCOMB,MXSB,MXSOOB,           &
     &                  IBLTP)
!
! Number of combimations with symmetry ITOTSM and
!       MNRS1 - MXRS1 elecs in RAS1
!       MNRS3 - MXRS3 elecs in RAS3
!
! In view of the limited range of I*4, the number of dets
! is returned as integer and  real*8
!
! MXSB is largest UNPACKED symmetry block
! MXSOOB is largest UNPACKED symmetry-type-type block
!
! Updated with IBLTP, Summer of 93
!
      use Symmetry_Info, only: Mul
      IMPLICIT None
      Integer MNRS1,MXRS1,MNRS3,MXRS3,ITOTSM,                           &
     &                  NSMST,NOCTPA,NOCTPB
      Integer IEL1A(*),IEL1B(*)
      Integer NSSOA(NOCTPA,*),NSSOB(NOCTPB,*)
      Integer IEL3A(*),IEL3B(*)
      Integer NCOMB
      Real*8 XNCOMB
      Integer MXSB,MXSOOB
      Integer IBLTP(*)

! local variables
      Integer IASM,LSB,IBSM,ISYM,IATP,MXBTP,IBTP,IEL1,IEL3,LTTSBL,      &
     &        LTTSUP,NTEST
!
      MXSB = 0
      MXSOOB = 0
      NCOMB = 0
      XNCOMB = 0.0D0
      DO 300 IASM = 1, NSMST
        IF(IBLTP(IASM).EQ.0) GOTO 300
        IBSM = Mul(IASM,ITOTSM)
        LSB = 0
        IF(IBSM.NE.0) THEN
          IF(IBLTP(IASM).EQ.2) THEN
            ISYM = 1
          ELSE
            ISYM = 0
          END IF
          DO 200 IATP = 1, NOCTPA
           IF(ISYM.EQ.1) THEN
             MXBTP = IATP
           ELSE
             MXBTP = NOCTPB
           END IF
           DO 100 IBTP = 1, MXBTP
             IEL1 = IEL1A(IATP)+IEL1B(IBTP)
             IEL3 = IEL3A(IATP)+IEL3B(IBTP)
             IF(IEL1.GE.MNRS1.AND.IEL1.LE.MXRS1.AND.                    &
     &       IEL3.GE.MNRS3.AND.IEL3.LE.MXRS3 ) THEN
!. Size of unpacked block
               LTTSUP =  NSSOA(IATP,IASM)*NSSOB(IBTP,IBSM)
!. Size of packed block
               IF(ISYM.EQ.0.OR.IATP.NE.IBTP) THEN
                 LTTSBL = NSSOA(IATP,IASM)*NSSOB(IBTP,IBSM)
               ELSE
                 LTTSBL = NSSOA(IATP,IASM)*(NSSOA(IATP,IASM)+1)/2
               END IF
               NCOMB = NCOMB + LTTSBL
               LSB = LSB + LTTSUP
               MXSOOB = MAX(MXSOOB,LTTSUP)
               IF(ISYM.EQ.0.OR.IATP.NE.IBTP) THEN
                 XNCOMB = XNCOMB +                                      &
     &         DBLE(NSSOA(IATP,IASM))*DBLE(NSSOB(IBTP,IBSM))
               ELSE
                 XNCOMB = XNCOMB +                                      &
     &           DBLE(NSSOA(IATP,IASM))*                                &
     &           (DBLE(NSSOB(IBTP,IBSM))+1.0D0)/2.0D0
               END IF
             END IF
  100      CONTINUE
  200     CONTINUE
          MXSB = MAX(MXSB,LSB)
        END IF
  300 CONTINUE
!
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' NCOMB and XNCOMB ', NCOMB,XNCOMB
      END IF
!
      END SUBROUTINE NRASDT

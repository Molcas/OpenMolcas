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
      SUBROUTINE SS_HEAD_DBL_TAIL_ACT_G(LRA)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      COMMON/ONEPL/LINE,JPH,JPEL,JPER,LRG,LRS,JWL,JWR,W0,W1
      w0ss1=0.d0
      w1ss1=0.d0
      w0ss2=0.d0
      w1ss2=0.d0
      w0ss3=0.d0
      w1ss3=0.d0
      w0ss4=0.d0
      w1ss4=0.d0
      w0ss5=0.d0
      w1ss5=0.d0
      w0ss6=0.d0
      w1ss6=0.d0
      w0ss7=0.d0
      w1ss7=0.d0
      w0ss8=0.d0
      w1ss8=0.d0
      w0ss9=0.d0
      w1ss9=0.d0
      w0ss10=0.d0
      w1ss10=0.d0
      w0ss11=0.d0
      w1ss11=0.d0
      w0ss12=0.d0
      w1ss12=0.d0
      w0ss13=0.d0
      w1ss13=0.d0
      w0ss14=0.d0
      w1ss14=0.d0
      w0ss15=0.d0
      w1ss15=0.d0
      w0ss16=0.d0
      w1ss16=0.d0
      w0ss18=0.d0
      w1ss18=0.d0

!SS(1-1)  Ar(01)-Bl(32)-
!SS(1-2)  Ar(02)-Bl(31)-
!SS(1-3)  Ar(13)-Bl(20)-
!SS(1-4)  Ar(23)-Bl(10)-
!SS(1-5)  (22)-Ar(13)-Bl(31)-
!SS(1-6)  (11)-Ar(23)-Bl(32)-
!SS(1-7)  Ar(13)-C'(21)-Bl(32)-
!SS(1-8)  Ar(13)-C'(22)-Bl(31)-
!SS(1-9)  Ar(23)-C'(11)-Bl(32)-
!SS(1-10) Ar(23)-C'(12)-Bl(31)-
!SS(1-11) Ar(13)-Bl(31)-C"(22)-
!SS(1-12) Ar(13)-Bl(32)-C"(21)-
!SS(1-13) Ar(23)-Bl(31)-C"(12)-
!SS(1-14) Ar(23)-Bl(32)-C"(11)-
!SS(1-15) (22)-Drl(11)-
!SS(1-16) (11)-Drl(22)-
!SS(1-17) Drl(22)-C"(11)-
!SS(1-18) Drl(11)-C"(22)-
!SS(1-19) Drl(12)-C"(21)-
!SS(1-20) Drl(33)-C"(00)-
!SS(1-20) Drl(33)-C"(11)-C"(22)-
!SS(1-20) (11)Drl(33)-C"(22)-
!SS(1-20) (11)(22)Drl(33)-
!SS(1-20) Drl(33)-C"(22)-C"(11)-
!SS(1-20) (22)Drl(33)-C"(11)-
!SS(1-20) (22)(11)Drl(33)-
!ss(1-20) (21)(12)Drl(33)-


      DO LRI=NORB_FRZ+1,NORB_DZ-1
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          JMLR=MUL_TAB(JML,JMR)
          IF(LMIJ.NE.JMLR) CYCLE
          W0SS2=W0_SS(2)
          W1SS2=W1_SS(2)
          W0SS4=W0_SS(4)
          W1SS4=W1_SS(4)
          W0SS5=W0_SS(5)
          W1SS5=W1_SS(5)
          W0SS10=-W0_SS(10)
          W1SS10=-W1_SS(10)
          W0SS14=W0_SS(14)
          W1SS14=W1_SS(14)
          IF(JB_SYS.GT.0) THEN
            W0SS1=W0_SS(1)
            W1SS1=W1_SS(1)
            W0SS3=W0_SS(3)
            W1SS3=W1_SS(3)
            W0SS6=W0_SS(6)
            W1SS6=W1_SS(6)
            W0SS7=W0_SS(7)
            W1SS7=W1_SS(7)
            W0SS8=W0_SS(8)
            W1SS8=W1_SS(8)
            W0SS9=W0_SS(9)
            W1SS9=W1_SS(9)
            W0SS11=W0_SS(11)
            W1SS11=W1_SS(11)
            W0SS12=W0_SS(12)
            W1SS12=W1_SS(12)
            W0SS13=W0_SS(13)
            W1SS13=W1_SS(13)
          ENDIF
          NI=MOD(LRJ-LRI,2)
          IF(NI.EQ.0) THEN
            W0SS2=-W0SS2
            W1SS2=-W1SS2
            W0SS4=-W0SS4
            W1SS4=-W1SS4
            W0SS5=-W0SS5
            W1SS5=-W1SS5
            W0SS10=-W0SS10
            W1SS10=-W1SS10
            W0SS14=-W0SS14
            W1SS14=-W1SS14
            IF(JB_SYS.GT.0) THEN
              W0SS1=-W0SS1
              W1SS1=-W1SS1
              W0SS3=-W0SS3
              W1SS3=-W1SS3
              W0SS6=-W0SS6
              W1SS6=-W1SS6
              W0SS7=-W0SS7
              W1SS7=-W1SS7
              W0SS8=-W0SS8
              W1SS8=-W1SS8
              W0SS9=-W0SS9
              W1SS9=-W1SS9
              W0SS11=-W0SS11
              W1SS11=-W1SS11
              W0SS12=-W0SS12
              W1SS12=-W1SS12
              W0SS13=-W0SS13
              W1SS13=-W1SS13
            ENDIF
          ENDIF

          IF(JML.EQ.1.AND.LMIJ.EQ.JMR) THEN
            IWDL=JUST(LRI,LRI)
!SS(1-1)   Ar(01)-Bl(32)-
            IF(JB_SYS.GT.0) THEN
              IWDR=JUST(LRJ,LRI)
              VLOP0=W0*W0SS1
              VLOP1=W1*W1SS1
      IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
        CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
      ENDIF
      IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
        CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
      ENDIF
      IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
        CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
      ENDIF

            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :      CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

            ENDIF
!SS(1-2)  Ar(02)-Bl(31)-
            IWDR=JUST(LRI,LRJ)
            VLOP0=W0*W0SS2
            VLOP1=W1*W1SS2
      IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
        CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
      ENDIF
      IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
        CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
      ENDIF
      IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
        CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
      ENDIF
        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :  CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

         ENDIF

          IF(JMR.EQ.1.AND.LMIJ.EQ.JML) THEN
            IWDR=JUST(LRJ,LRJ)
!SS(1-3)  Ar(13)-Bl(20)-
            IF(JB_SYS.GT.0) THEN
              IWDL=JUST(LRJ,LRI)
              VLOP0=W0*W0SS3
              VLOP1=W1*W1SS3
      IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
        CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
      ENDIF
      IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
        CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
      ENDIF
      IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
        CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
      ENDIF
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :         CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

            ENDIF
!SS(1-4)  Ar(23)-Bl(10)-        ACT -C"-                         ! IPRAD
            IWDL=JUST(LRI,LRJ)
            VLOP0=W0*W0SS4
            VLOP1=W1*W1SS4
      IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
        CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
      ENDIF
      IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
        CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
      ENDIF
      IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
        CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
      ENDIF
        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :  CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)
         ENDIF

          VLOP0=W0*W0SS5
          VLOP1=W1*W1SS5
          IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
            CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
          IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
            CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
          IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
            CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
        IF(JB_SYS.GT.0) THEN
          VLOP0=W0*W0SS6
          VLOP1=W1*W1SS6
          IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
           CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
          ENDIF
          IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
           CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
          ENDIF
          IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
           CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
          ENDIF
        ENDIF
          DO LRK=NORB_FRZ+1,LRI-1
            LMK=LSM_INN(LRK)
            LMKI=MUL_TAB(LMK,LMI)
            LMKJ=MUL_TAB(LMK,LMJ)
            IF(LMKI.EQ.JML.AND.LMKJ.EQ.JMR) THEN
!SS(1-5)  (22)-Ar(13)-Bl(31)-
              IWDL=JUST(LRK,LRI)
              IWDR=JUST(LRK,LRJ)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :         CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

              IF(JB_SYS.GT.0) THEN
!SS(1-6)  (11)-Ar(23)-Bl(32)-
                IWDL=JUST(LRI,LRK)
                IWDR=JUST(LRJ,LRK)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0_1,JPER,LIST0_1)
      IF(LIST1_1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1_1,JPER,LIST1_1)

              ENDIF
            ENDIF
          ENDDO

       IF(JB_SYS.GT.0) THEN
          VLOP0=W0*W0SS7
          VLOP1=W1*W1SS7
          IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
           CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
          ENDIF
          IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
           CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
          ENDIF
          IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
           CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
          ENDIF
          VLOP0=W0*W0SS8
          VLOP1=W1*W1SS8
          IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
           CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_2,LIST0_2,
     :                     WL1_2,LIST1_2)
          ENDIF
          IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
           CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_2,LIST0_2,
     :                     WL1_2,LIST1_2)
          ENDIF
          IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
           CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_2,LIST0_2,
     :                     WL1_2,LIST1_2)
          ENDIF
          VLOP0=W0*W0SS9
          VLOP1=W1*W1SS9
          IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
           CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_3,LIST0_3,
     :                     WL1_3,LIST1_3)
          ENDIF
          IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
           CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_3,LIST0_3,
     :                     WL1_3,LIST1_3)
          ENDIF
          IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
           CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_3,LIST0_3,
     :                     WL1_3,LIST1_3)
          ENDIF
        ENDIF
          VLOP0=W0*W0SS10
          VLOP1=W1*W1SS10
          IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
           CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_4,LIST0_4,
     :                     WL1_4,LIST1_4)
          ENDIF
          IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
           CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_4,LIST0_4,
     :                     WL1_4,LIST1_4)
          ENDIF
          IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
           CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_4,LIST0_4,
     :                     WL1_4,LIST1_4)
          ENDIF

          DO LRK=LRI+1,LRJ-1
            LMK=LSM_INN(LRK)
            LMKI=MUL_TAB(LMK,LMI)
            LMKJ=MUL_TAB(LMK,LMJ)
            IF(LMKI.EQ.JML.AND.LMKJ.EQ.JMR) THEN
              IF(JB_SYS.GT.0) THEN
!SS(1-7)  Ar(13)-C'(21)-Bl(32)-
                IWDL=JUST(LRK,LRI)
                IWDR=JUST(LRJ,LRK)
C                CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,-WL1,JPER)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL0_1,JPER,LIST0_1)
      IF(LIST1_1.NE.0)
     :      CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL1_1,JPER,LIST1_1)

!SS(1-8)  Ar(13)-C'(22)-Bl(31)-
                IWDL=JUST(LRK,LRI)
                IWDR=JUST(LRK,LRJ)
C                CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,-WL2,JPER)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL0_2,JPER,LIST0_2)
      IF(LIST1_2.NE.0)
     :      CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL1_2,JPER,LIST1_2)

!SS(1-9)  Ar(23)-C'(11)-Bl(32)-
                IWDL=JUST(LRI,LRK)
                IWDR=JUST(LRJ,LRK)
!               WL=(VLOP0-VLOP1)*vint_ci(LIST)-2*VLOP0*vint_ci(LIST+1)
C                CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,-WL3,JPER)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL0_3,JPER,LIST0_3)
      IF(LIST1_3.NE.0)
     :      CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL1_3,JPER,LIST1_3)

              ENDIF
!SS(1-10) Ar(23)-C'(12)-Bl(31)-
              IWDL=JUST(LRI,LRK)
              IWDR=JUST(LRK,LRJ)
C              CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL4,JPER)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0_4,JPER,LIST0_4)
      IF(LIST1_4.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1_4,JPER,LIST1_4)

            ENDIF
          ENDDO

        IF(JB_SYS.GT.0) THEN
          VLOP0=W0*W0SS11
          VLOP1=W1*W1SS11
          IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
           CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
          ENDIF
          IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
           CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
          ENDIF
          IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
           CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
          ENDIF
          VLOP0=W0*W0SS12
          VLOP1=W1*W1SS12
          IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
           CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_2,LIST0_2,
     :                     WL1_2,LIST1_2)
          ENDIF
          IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
           CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_2,LIST0_2,
     :                     WL1_2,LIST1_2)
          ENDIF
          IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
           CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_2,LIST0_2,
     :                     WL1_2,LIST1_2)
          ENDIF
          VLOP0=W0*W0SS13
          VLOP1=W1*W1SS13
          IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
           CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_3,LIST0_3,
     :                     WL1_3,LIST1_3)
          ENDIF
          IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
           CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_3,LIST0_3,
     :                     WL1_3,LIST1_3)
          ENDIF
          IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
           CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_3,LIST0_3,
     :                     WL1_3,LIST1_3)
          ENDIF
        ENDIF
          VLOP0=W0*W0SS14
          VLOP1=W1*W1SS14
          IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
           CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_4,LIST0_4,
     :                     WL1_4,LIST1_4)
          ENDIF
          IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
           CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_4,LIST0_4,
     :                     WL1_4,LIST1_4)
          ENDIF
          IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
           CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_4,LIST0_4,
     :                     WL1_4,LIST1_4)
          ENDIF

          DO LRK=LRJ+1,NORB_DZ
            LMK=LSM_INN(LRK)
            LMKI=MUL_TAB(LMK,LMI)
            LMKJ=MUL_TAB(LMK,LMJ)
            IF(LMKI.EQ.JML.AND.LMKJ.EQ.JMR) THEN
              IF(JB_SYS.GT.0) THEN
!SS(1-11) Ar(13)-Bl(31)-C"(22)-
                IWDL=JUST(LRK,LRI)
                IWDR=JUST(LRK,LRJ)
C                CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0_1,JPER,LIST0_1)
      IF(LIST1_1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1_1,JPER,LIST1_1)

!SS(1-12) Ar(13)-Bl(32)-C"(21)-
                IWDL=JUST(LRK,LRI)
                IWDR=JUST(LRJ,LRK)
C                CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0_2,JPER,LIST0_2)
      IF(LIST1_2.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1_2,JPER,LIST1_2)

!SS(1-13) Ar(23)-Bl(31)-C"(12)-
                IWDL=JUST(LRI,LRK)
                IWDR=JUST(LRK,LRJ)
C                CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL3,JPER)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0_3,JPER,LIST0_3)
      IF(LIST1_3.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1_3,JPER,LIST1_3)

              ENDIF
!SS(1-14) Ar(23)-Bl(32)-C"(11)- ACT -C"-
              IWDL=JUST(LRI,LRK)
              IWDR=JUST(LRJ,LRK)
C              CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL4,JPER)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0_4,JPER,LIST0_4)
      IF(LIST1_4.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1_4,JPER,LIST1_4)

            ENDIF
          ENDDO
        ENDDO
      ENDDO

      IF(JPAD.NE.JPADL) RETURN

!      IF(JB_SYS.GT.0.OR.JWL.GE.JWR) THEN
      IF(JB_SYS.GT.0) THEN
      IF(JPAD.GT.17.AND.JPAD.LT.25) THEN
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JML.OR.JML.NE.JMR) CYCLE
            IWDL=JUST(LRJ,LRI)
            IWDR=JUST(LRI,LRJ)
!SS(1-19) Drl(12)-C"(21)-
            VLOP0=W0*W0_SS(19)
            VLOP1=W1*W1_SS(19)
            IF(LINE.EQ.26) THEN    !LRI,LRA
              CALL COMP_LOOP_G(9,LRI,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
            IF(LINE.EQ.28) THEN    !LRI,LRS,LRA
              CALL COMP_LOOP_G(12,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
            IF(LINE.EQ.29) THEN    !LRI,LRS,LRA
              CALL COMP_LOOP_G(11,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
!         WL=(VLOP0-VLOP1)*VOINT(LRI,LRB)
C           CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

        ENDDO
      ENDDO
      ENDIF
      ENDIF
      IF(JWL.GE.JWR) RETURN

      W0SS15=W0_SS(15)
      W1SS15=W1_SS(15)
      W0SS17=W0_SS(17)
      W1SS17=W1_SS(17)
      W0SS20=W0_SS(20)
      IF(JB_SYS.GT.0) THEN
        W0SS16=W0_SS(16)
        W1SS16=W1_SS(16)
        W0SS18=W0_SS(18)
        W1SS18=W1_SS(18)
      ENDIF

      IF(JML.EQ.1.AND.JMR.EQ.1) THEN
!SS(1-20) Drl(33)-C"(00)-                            ! IPL(R)AD=1 or =NS
        DO LR0=NORB_FRZ+1,NORB_DZ
          IWDL=JUST(LR0,LR0)
          IWDR=IWDL
          VLOP0=W0*W0_SS(20)
          VLOP1=0.D0
C         WL=0.D0
          DO LRK=1,NORB_DZ
            IF(LRK.EQ.LR0) CYCLE
            IF(LINE.EQ.26) THEN    !LRK,LRA
             CALL COMP_LOOP_G(9,LRK,0,0,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
            ENDIF
            IF(LINE.EQ.28) THEN    !LRK,LRS,LRA    LRS,LRA,LRK
        CALL COMP_LOOP_G(12,LRK,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
            ENDIF
            IF(LINE.EQ.29) THEN    !LRK,LRS,LRA
        CALL COMP_LOOP_G(11,LRK,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
            ENDIF
C            WL=WL+WL1
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0_1,JPER,LIST0_1)
      IF(LIST1_1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1_1,JPER,LIST1_1)

          ENDDO
C          CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
        ENDDO
      ENDIF
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JML.OR.JML.NE.JMR) CYCLE
          IF(JWL.GE.JWR.AND.JB_SYS.EQ.0) CYCLE
          WL=0.D0
          IWDL=JUST(LRI,LRJ)
          IWDR=IWDL
!SS(1-15) (22)-Drl(11)-
          VLOP0=W0*W0SS15
          VLOP1=W1*W1SS15
C          WL=0.D0
          IF(LINE.EQ.26) THEN    !LRJ,LRA
            CALL COMP_LOOP_G(9,LRJ,0,0,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
          ENDIF
          IF(LINE.EQ.28) THEN    !LRJ,LRS,LRA
            CALL COMP_LOOP_G(12,LRJ,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
          ENDIF
          IF(LINE.EQ.29) THEN    !LRJ,LRS,LRA
            CALL COMP_LOOP_G(11,LRJ,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
          ENDIF
C          WL=WL+WL1
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0_1,JPER,LIST0_1)
      IF(LIST1_1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1_1,JPER,LIST1_1)

!         WL=WL+(VLOP0-VLOP1)*VOINT(LRJ,LRB)
!SS(1-17) Drl(22)-C"(11)-
          VLOP0=W0*W0SS17
         VLOP1=W1*W1SS17
          IF(LINE.EQ.26) THEN    !LRI,LRA
            CALL COMP_LOOP_G(9,LRI,0,0,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
          ENDIF
          IF(LINE.EQ.28) THEN    !LRI,LRS,LRA
            CALL COMP_LOOP_G(12,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
          ENDIF
          IF(LINE.EQ.29) THEN    !LRI,LRS,LRA
            CALL COMP_LOOP_G(11,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
          ENDIF
C          WL=WL+WL1
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0_1,JPER,LIST0_1)
      IF(LIST1_1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1_1,JPER,LIST1_1)

!         WL=WL+(VLOP0-VLOP1)*VOINT(LRI,LRB)
!SS(1-20) (22)(11)Drl(33)-
!SS(1-20) (22)Drl(33)-C"(11)-
!SS(1-20) Drl(33)-C"(22)-C"(11)-
          VLOP0=W0*W0SS20
          VLOP1=0.D0
         DO LRK=1,NORB_DZ
            IF(LRK.EQ.LRI) CYCLE
            IF(LRK.EQ.LRJ) CYCLE
            IF(LINE.EQ.26) THEN    !LRK,LRA
              CALL COMP_LOOP_G(9,LRK,0,0,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
            ENDIF
            IF(LINE.EQ.28) THEN    !LRK,LRS,LRA
            CALL COMP_LOOP_G(12,LRK,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
            ENDIF
            IF(LINE.EQ.29) THEN    !LRK,LRS,LRA
            CALL COMP_LOOP_G(11,LRK,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
            ENDIF
C            WL=WL+WL1
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0_1,JPER,LIST0_1)
      IF(LIST1_1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1_1,JPER,LIST1_1)

!           WL=WL+VLOP0*VOINT(LRK,LRB)
          ENDDO
C          CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
          IF(JB_SYS.GT.0) THEN
            IWDL=JUST(LRJ,LRI)
            IWDR=IWDL
!SS(1-16) (11)-Drl(22)-
            VLOP0=W0*W0SS16
           VLOP1=W1*W1SS16
C          WL=0.D0
            IF(LINE.EQ.26) THEN    !LRJ,LRA
            CALL COMP_LOOP_G(9,LRJ,0,0,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
            ENDIF
            IF(LINE.EQ.28) THEN    !LRJ,LRS,LRA
          CALL COMP_LOOP_G(12,LRJ,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
            ENDIF
            IF(LINE.EQ.29) THEN    !LRJ,LRS,LRA
          CALL COMP_LOOP_G(11,LRJ,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
            ENDIF
c            WL=WL+WL1
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0_1,JPER,LIST0_1)
      IF(LIST1_1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1_1,JPER,LIST1_1)

!            WL=(VLOP0-VLOP1)*VOINT(LRJ,LRB)
!SS(1-18) Drl(11)-C"(22)-
            VLOP0=W0*W0SS18
            VLOP1=W1*W1SS18
            IF(LINE.EQ.26) THEN    !LRI,LRA
            CALL COMP_LOOP_G(9,LRI,0,0,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
            ENDIF
            IF(LINE.EQ.28) THEN    !LRI,LRS,LRA
            CALL COMP_LOOP_G(12,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
            ENDIF
            IF(LINE.EQ.29) THEN    !LRI,LRS,LRA
            CALL COMP_LOOP_G(11,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
            ENDIF
C            WL=WL+WL1
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0_1,JPER,LIST0_1)
      IF(LIST1_1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1_1,JPER,LIST1_1)

!            WL=WL+(VLOP0-VLOP1)*VOINT(LRI,LRB)
!SS(1-20) Drl(33)-C"(11)-C"(22)-
!SS(1-20) (11)Drl(33)-C"(22)-
!SS(1-20) (11)(22)Drl(33)-
            VLOP0=W0*W0SS20
            VLOP1=0.D0
            DO LRK=1,NORB_DZ
              IF(LRK.EQ.LRI) CYCLE
              IF(LRK.EQ.LRJ) CYCLE
              IF(LINE.EQ.26) THEN    !LRK,LRA
           CALL COMP_LOOP_G(9,LRK,0,0,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
              ENDIF
              IF(LINE.EQ.28) THEN    !LRK,LRS,LRA
           CALL COMP_LOOP_G(12,LRK,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
              ENDIF
              IF(LINE.EQ.29) THEN    !LRK,LRS,LRA
           CALL COMP_LOOP_G(11,LRK,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
              ENDIF
C              WL=WL+WL1
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0_1,JPER,LIST0_1)
      IF(LIST1_1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1_1,JPER,LIST1_1)

!              WL=WL+VLOP0*VOINT(LRK,LRB)
            ENDDO
C            CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
           ENDIF
        ENDDO
      ENDDO
      RETURN

      END

      SUBROUTINE ST_HEAD_DBL_TAIL_ACT_G(LRA)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      COMMON/ONEPL/LINE,JPH,JPEL,JPER,LRG,LRS,JWL,JWR,W0,W1
!ST(2-1) Ar(02)-Bl(32)-
!ST(2-2) (22)Ar(13)-Bl(32)-
!ST(2-4) Ar(23)-C'(12)-Bl(32)-
!ST(2-4) Ar(23)-Bl(32)-C'(12)-
!ST(2-5) (22)Drl(12)-
!ST(2-6) Drl(22)-C"(12)-
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          LMIJ=MUL_TAB(LMIJ,1)
          IF(JML.EQ.JMR.AND.LMIJ.EQ.JML) THEN
            IWDS=JUST(LRI,LRJ)
            IWDT=IWDS              !
!ST(2-5) (22)Drl(12)-
            VLOP1=W1*W1_ST(5)             !D2-5
            VLOP0=0.D0
            IF(LINE.EQ.26) THEN    !LRJ,LRA
              CALL COMP_LOOP_G(9,LRJ,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
            IF(LINE.EQ.28) THEN    !LRJ,LRS,LRA
              CALL COMP_LOOP_G(12,LRJ,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
            IF(LINE.EQ.29) THEN    !LRJ,LRS,LRA
              CALL COMP_LOOP_G(11,LRJ,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
C           CALL PRODAB_1(3,JPEL,IWDS,IWDT,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,WL1,JPER,LIST1)

!ST(2-6) Drl(22)-C"(12)-
            VLOP1=W1*W1_ST(6)             !D2-6
            VLOP0=0.D0
            IF(LINE.EQ.26) THEN    !LRI,LRA
          CALL COMP_LOOP_G(9,LRI,0,0,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
            ENDIF
            IF(LINE.EQ.28) THEN    !LRI,LRS,LRA
          CALL COMP_LOOP_G(12,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
            ENDIF
            IF(LINE.EQ.29) THEN    !LRI,LRS,LRA
          CALL COMP_LOOP_G(11,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
            ENDIF

C          CALL PRODAB_1(3,JPEL,IWDS,IWDT,JWL,JWR,WL2_1,JPER,LIST2_1,
C    :                                                       LIST3_1)
              CALL PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,WL0_1,JPER,LIST0_1)
      IF(LIST1_1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,WL1_1,JPER,LIST1_1)

C           WL=WL+WL1
!            LIST=LIST3(LRA,LRB,LRI)
!            WL=WL-VLOP1*vint_ci(LIST)    !4.3 Vlop0=0        !!!!!
C            CALL PRODAB(3,JPEL,IWDS,IWDT,JWL,JWR,WL,JPER)
!ST(2-7) Drl(12)-C"(22)-
            IF(JB_SYS.GT.0) THEN
              IWDS=JUST(LRJ,LRI)
              IWDT=JUST(LRI,LRJ)
            VLOP1=W1*W1_ST(7)
            VLOP0=0.D0             !D2-6
!              LIST=LIST3(LRA,LRB,LRI)
!              WL=WL-VLOP1*vint_ci(LIST)    !4.3 Vlop0=0        !!!!!
              IF(LINE.EQ.26) THEN    !LRI,LRA
                CALL COMP_LOOP_G(9,LRI,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
              ENDIF
              IF(LINE.EQ.28) THEN    !LRI,LRS,LRA
                CALL COMP_LOOP_G(12,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
              ENDIF
              IF(LINE.EQ.29) THEN    !LRI,LRS,LRA
                CALL COMP_LOOP_G(11,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
              ENDIF
C           CALL PRODAB_1(3,JPEL,IWDS,IWDT,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,WL1,JPER,LIST1)

C              CALL PRODAB(3,JPEL,IWDS,IWDT,JWL,JWR,WL,JPER)
            ENDIF
          ENDIF

          JMLR=MUL_TAB(JML,JMR)
         IF(LMIJ.NE.JMLR) CYCLE
          W1ST1=W1_ST(1)
          W1ST2=W1_ST(2)
          W1ST3=W1_ST(3)
          W1ST4=W1_ST(4)
          NI=MOD(LRJ-LRI,2)
          IF(NI.EQ.0) THEN
            W1ST1=-W1ST1
            W1ST2=-W1ST2
            W1ST3=-W1ST3
            W1ST4=-W1ST4
          ENDIF
          IF(JML.EQ.1)THEN
!ST(2-1) Ar(02)-Bl(32)-
            IWDS=JUST(LRI,LRI)
            IWDT=JUST(LRI,LRJ)      !
            VLOP1=W1*W1ST1
           VLOP0=0.D0
!          LIST=LIST4(LRI,LRJ,LRA,LRB)
            IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
              CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
            IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
               CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
            IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
              CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
!          WL=-VLOP1*vint_ci(LIST)    !1.1 VLOP0=0        !!!!!
C            CALL PRODAB(3,JPEL,IWDS,IWDT,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDS,IWDT,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,WL1,JPER,LIST1)

          ENDIF

!ST(2-2) (22)Ar(13)-Bl(32)-
          VLOP1=W1*W1ST2
!       WL=-VLOP1*vint_ci(LIST)    !1.1
          VLOP0=0.D0
          IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
            CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
          IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
            CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
          IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
            CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
         DO LRK=NORB_FRZ+1,LRI-1
            LMK=LSM_INN(LRK)
            IF(MUL_TAB(LMK,LMI).NE.JML) CYCLE
            IF(MUL_TAB(LMK,LMJ).NE.JMR) CYCLE
            IWDS=JUST(LRK,LRI)
            IWDT=JUST(LRK,LRJ)      !
C            CALL PRODAB(3,JPEL,IWDS,IWDT,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDS,IWDT,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,WL1,JPER,LIST1)
          ENDDO
!ST(2-3) Ar(13)-Bl(32)-C'(22)-
!ST(2-4) Ar(23)-Bl(32)-C'(12)-
          VLOP1=W1*W1ST4
         VLOP0=0.D0
          IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
            CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
          IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
            CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
          IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
            CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
!        WL=-VLOP1*vint_ci(LIST)    !1.1 VLOP0=0
        IF(JB_SYS.GT.0) THEN
          VLOP1=W1*W1ST3
          VLOP0=0.D0
          IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
        CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
          ENDIF
          IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
        CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
          ENDIF
          IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
        CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,
     :                     WL1_1,LIST1_1)
          ENDIF
        ENDIF
!         WL1=-VLOP1*vint_ci(LIST)
         DO LRK=LRJ+1,NORB_DZ
            LMK=LSM_INN(LRK)
            IF(MUL_TAB(LMK,LMI).NE.JML) CYCLE
            IF(MUL_TAB(LMK,LMJ).NE.JMR) CYCLE
            IWDS=JUST(LRI,LRK)
            IWDT=JUST(LRJ,LRK)     !
C            CALL PRODAB(3,JPEL,IWDS,IWDT,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDS,IWDT,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,WL1,JPER,LIST1)

            IF(JB_SYS.GT.0) THEN
              IWDS=JUST(LRK,LRI)
C              CALL PRODAB(3,JPEL,IWDS,IWDT,JWL,JWR,WL1,JPER)
C          CALL PRODAB_1(3,JPEL,IWDS,IWDT,JWL,JWR,WL2_1,JPER,LIST2_1,
C    :                                                       LIST3_1)
              CALL PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,WL0_1,JPER,LIST0_1)
      IF(LIST1_1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,WL1_1,JPER,LIST1_1)

            ENDIF
          ENDDO
!ST(2-4) Ar(23)-C'(12)-Bl(32)-
         DO LRK=LRI+1,LRJ-1
            LMK=LSM_INN(LRK)
            IF(MUL_TAB(LMK,LMI).NE.JML) CYCLE
            IF(MUL_TAB(LMK,LMJ).NE.JMR) CYCLE
            IWDS=JUST(LRI,LRK)
            IWDT=JUST(LRK,LRJ)    !
C            CALL PRODAB(3,JPEL,IWDS,IWDT,JWL,JWR,-WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDS,IWDT,JWL,JWR,-WL2,JPER,LIST2,LIST3
          CALL PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,-WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :       CALL PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,-WL1,JPER,LIST1)

!ST(2-3) Ar(13)-C'(22)-Bl(32)-
            IF(JB_SYS.GT.0) THEN
              IWDS=JUST(LRK,LRI)
C              CALL PRODAB(3,JPEL,IWDS,IWDT,JWL,JWR,-WL1,JPER)
C          CALL PRODAB_1(3,JPEL,IWDS,IWDT,JWL,JWR,-WL2_1,JPER,LIST2_1,
C    :                                                        LIST3_1)
           CALL PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,-WL0_1,JPER,LIST0_1)
      IF(LIST1_1.NE.0)
     :     CALL PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,-WL1_1,JPER,LIST1_1)

           ENDIF
          ENDDO
        ENDDO
      ENDDO
      RETURN

      END

      SUBROUTINE TS_HEAD_DBL_TAIL_ACT_G(LRA)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      COMMON/ONEPL/LINE,JPH,JPEL,JPER,LRG,LRS,JWL,JWR,W0,W1
!TS(3-1) Ar(23)-Bl(20)-
!TS(3-2) (22)Ar(23)-Bl(31)-
!TS(3-2) Ar(23)-C'(22)-Bl(31)-
!TS(3-3) Ar(23)-Bl(31)-C"(22)-
!TS(3-4) Ar(23)-Bl(32)-C"(21)-
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
      DO LRJ=LRI+1,NORB_DZ
        LMJ=LSM_INN(LRJ)
        LMIJ=MUL_TAB(LMI,LMJ)
        W1TS1=W1_TS(1)
        W1TS2=W1_TS(2)
        W1TS3=W1_TS(3)
        W1TS4=W1_TS(4)
        NI=MOD(LRJ-LRI,2)
        IF(NI.eq.0) THEN
          W1TS1=-W1TS1
          W1TS2=-W1TS2
          W1TS3=-W1TS3
          W1TS4=-W1TS4
        ENDIF
!        LIST=LIST3(LRI,LRJ,LRB)
!-------------------------------------------------------------------
!TS(3-1) Ar(23)-Bl(20)-
          IF(JMR.EQ.1.and.LMIJ.eq.JML) THEN
            IWDT=JUST(LRI,LRJ)    !
            IWDS=JUST(LRJ,LRJ)
            VLOP1=W1*W1TS1             !D3-1
            VLOP0=0.D0
            IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
              CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
            IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
              CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
            IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
              CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
!           WL=-VLOP1*(VINT_CI(LIST))    !2.2 vlop0=0
C            CALL PRODAB(3,JPEL,IWDT,IWDS,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDT,IWDS,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDT,IWDS,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDT,IWDS,JWL,JWR,WL1,JPER,LIST1)
          ENDIF
!-------------------------------------------------------------------
          VLOP1=W1*W1TS2             !D3-2
          VLOP0=0.D0
          IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
            CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
          IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
            CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
          IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
            CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
!          WL=-VLOP1*(VINT_CI(LIST))    !2.2 vlop0=0
!TS(3-2) (22)Ar(23)-Bl(31)-
          DO LRK=NORB_FRZ+1,LRI-1
            LMK=LSM_INN(LRK)
            IF(MUL_TAB(LMK,LMI).NE.JML.or.MUL_TAB(LMK,LMJ).NE.JMR) CYCLE
           IWDT=JUST(LRK,LRI)   !
            IWDS=JUST(LRK,LRJ)
C            CALL PRODAB(3,JPEL,IWDT,IWDS,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDT,IWDS,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDT,IWDS,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDT,IWDS,JWL,JWR,WL1,JPER,LIST1)
          ENDDO

!TS(3-2) Ar(23)-C'(22)-Bl(31)-
         DO LRK=LRI+1,LRJ-1
            LMK=LSM_INN(LRK)
            IF(MUL_TAB(LMK,LMI).NE.JML.or.MUL_TAB(LMK,LMJ).NE.JMR) CYCLE
           IWDT=JUST(LRI,LRK)   !
            IWDS=JUST(LRK,LRJ)
C            CALL PRODAB(3,JPEL,IWDT,IWDS,JWL,JWR,-WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDT,IWDS,JWL,JWR,-WL2,JPER,LIST2,LIST3
              CALL PRODAB_2(3,JPEL,IWDT,IWDS,JWL,JWR,-WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDT,IWDS,JWL,JWR,-WL1,JPER,LIST1)
          ENDDO

!-------------------------------------------------------------------
!TS(3-4) Ar(23)-Bl(32)-C"(21)-
          VLOP1=W1*W1TS4             !D3-4
          VLOP0=0.D0
          IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
            CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
          IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
            CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
          IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
            CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
          DO LRK=LRJ+1,NORB_DZ
            LMK=LSM_INN(LRK)
            IF(MUL_TAB(LMK,LMI).NE.JML.or.MUL_TAB(LMK,LMJ).NE.JMR) CYCLE
            IWDT=JUST(LRI,LRK)   !
            IWDS=JUST(LRJ,LRK)
C            CALL PRODAB(3,JPEL,IWDT,IWDS,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDT,IWDS,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDT,IWDS,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDT,IWDS,JWL,JWR,WL1,JPER,LIST1)

          ENDDO

!-------------------------------------------------------------------
!TS(3-3) Ar(23)-Bl(31)-C"(22)-
          IF(JB_SYS.GT.0) THEN
            VLOP1=W1*W1TS3             !D3-4
            VLOP0=0.D0
            IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
              CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
            IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
              CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
            IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
              CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
!            WL=-VLOP1*VINT_CI(LIST)    !2.2 vlop0=0
            DO LRK=LRJ+1,NORB_DZ
              LMK=LSM_INN(LRK)
       IF(MUL_TAB(LMK,LMI).NE.JML.or.MUL_TAB(LMK,LMJ).NE.JMR) CYCLE
              IWDT=JUST(LRI,LRK)   !
              IWDS=JUST(LRK,LRJ)
C              CALL PRODAB(3,JPEL,IWDT,IWDS,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDT,IWDS,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDT,IWDS,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDT,IWDS,JWL,JWR,WL1,JPER,LIST1)

            ENDDO
          ENDIF
        ENDDO
      ENDDO
      RETURN

      END

      SUBROUTINE STT_HEAD_DBL_TAIL_ACT_G(LRA)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      COMMON/ONEPL/LINE,JPH,JPEL,JPER,LRG,LRS,JWL,JWR,W0,W1
!ST1(4-1) Ar(01)-Bl(31)-
!ST1(4-2) (11)Ar(23)-Bl(31)-
!ST1(4-3) Ar(13)-C'(21)-Bl(31)-
!ST1(4-3) Ar(13)-Bl(31)-C"(21)-
!ST1(4-4) Ar(23)-C'(11)-Bl(31)-
!ST1(4-4) Ar(23)-Bl(31)-C"(11)-
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          W1ST1=W1_ST1(1)
          W1ST2=W1_ST1(2)
          W1ST3=W1_ST1(3)
          W1ST4=W1_ST1(4)
          NI=MOD(LRJ-LRI,2)
         IF(NI.EQ.0) THEN
            W1ST1=-W1ST1
            W1ST2=-W1ST2
            W1ST3=-W1ST3
            W1ST4=-W1ST4
          ENDIF
!          LIST=LIST3(LRI,LRJ,LRB)
         IF(JML.EQ.1.AND.LMIJ.EQ.JMR) THEN
!ST1(4-1) Ar(01)-Bl(31)-
             IWDL=JUST(LRI,LRI)
             IWDR=JUST(LRI,LRJ)
            VLOP1=W1*W1ST1
             VLOP0=0.D0
             IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
               CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
             ENDIF
             IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
               CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
             ENDIF
             IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
               CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
             ENDIF
!             WL=-VLOP1*(VINT_CI(LIST))    !2.2 vlop0=0
C             CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)
           ENDIF
!ST1(4-2) (11)Ar(23)-Bl(31)-
           VLOP1=W1*W1ST2
           VLOP0=0.D0
           IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
             CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
           ENDIF
           IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
             CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
           ENDIF
           IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
             CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
           ENDIF
           DO LRK=NORB_FRZ+1,LRI-1
             LMK=LSM_INN(LRK)
             IF(MUL_TAB(LMK,LMI).NE.JML.OR.
     :            MUL_TAB(LMK,LMJ).NE.JMR) CYCLE
!            WL=-VLOP1*VINT_CI(LIST)
             IWDL=JUST(LRI,LRK)
             IWDR=JUST(LRK,LRJ)
C            CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)
          ENDDO
           VLOP1=W1*W1ST4
           VLOP0=0.D0
           IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
              CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
           ENDIF
           IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
             CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
           ENDIF
           IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
             CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
           ENDIF
!ST1(4-4) Ar(23)-C'(11)-Bl(31)-
           DO LRK=LRI+1,LRJ-1
             LMK=LSM_INN(LRK)
             IF(MUL_TAB(LMI,LMK).NE.JML.OR.
     :             MUL_TAB(LMK,LMJ).NE.JMR) CYCLE
             IWDL=JUST(LRI,LRK)
             IWDR=JUST(LRK,LRJ)
C            CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,-WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,-WL2,JPER,LIST2,LIST3
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL1,JPER,LIST1)
          ENDDO
!ST1(4-4) Ar(23)-Bl(31)-C"(11)-
           DO LRK=LRJ+1,NORB_DZ
             LMK=LSM_INN(LRK)
            IF(MUL_TAB(LMI,LMK).NE.JML.OR.
     :          MUL_TAB(LMJ,LMK).NE.JMR) CYCLE
             IWDL=JUST(LRI,LRK)
             IWDR=JUST(LRJ,LRK)
C            CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)
          ENDDO

           VLOP1=W1*W1ST3
           VLOP0=0.D0
           IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
             CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
           ENDIF
           IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
             CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
           ENDIF
           IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
             CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
           ENDIF
!ST1(4-3) Ar(13)-C'(21)-Bl(31)-
           DO LRK=LRI+1,LRJ-1
             LMK=LSM_INN(LRK)
             IF(MUL_TAB(LMI,LMK).NE.JML.OR.
     :         MUL_TAB(LMK,LMJ).NE.JMR) CYCLE
             IWDL=JUST(LRK,LRI)
             IWDR=JUST(LRK,LRJ)
C            CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,-WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,-WL2,JPER,LIST2,LIST3
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL1,JPER,LIST1)
          ENDDO
!ST1(4-3) Ar(13)-Bl(31)-C"(21)-
           DO LRK=LRJ+1,NORB_DZ
             LMK=LSM_INN(LRK)
             IF(MUL_TAB(LMI,LMK).NE.JML.OR.
     :            MUL_TAB(LMJ,LMK).NE.JMR) CYCLE
             IWDL=JUST(LRK,LRI)
             IWDR=JUST(LRJ,LRK)
C            CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)
          ENDDO
        ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE TTS_HEAD_DBL_TAIL_ACT_G(LRA)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      COMMON/ONEPL/LINE,JPH,JPEL,JPER,LRG,LRS,JWL,JWR,W0,W1
!T1S(5-1)   Ar(13)-Bl(10)-
!T1S(5-2)   Ar(13)-Bl(32)-
!T1S(5-2)   Ar(13)-C'(11)-Bl(32)-
!T1S(5-3)   Ar(13)-Bl(31)-C"(12)-
!T1S(5-4)   Ar(13)-Bl(32)-C"(11)-
!T1S(5-5)   Drl(12)-
!T1S(5-6)   Drl(12)-C"(12)-
!T1S(5-7)   Drl(12)-C"(11)-
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
      DO LRJ=LRI+1,NORB_DZ
        LMJ=LSM_INN(LRJ)
        LMIJ=MUL_TAB(LMI,LMJ)
        W1T1S1=W1_T1S(1)
        W1T1S2=W1_T1S(2)
        W1T1S3=W1_T1S(3)
        W1T1S4=W1_T1S(4)
        W1T1S5=W1_T1S(5)
        W1T1S6=W1_T1S(6)
        W1T1S7=W1_T1S(7)
        NI=MOD(LRJ-LRI,2)
        IF(NI.EQ.0) THEN
          W1T1S1=-W1T1S1
          W1T1S2=-W1T1S2
          W1T1S3=-W1T1S3
          W1T1S4=-W1T1S4
!         W1T1S5=-W1T1S5
!         W1T1S6=-W1T1S6
!         W1T1S7=-W1T1S7
        ENDIF
!       LIST=LIST3(LRI,LRJ,LRB)
        IF(JML.EQ.LMIJ.AND.JMR.EQ.1) THEN
!T1S(5-1)   Ar(13)-Bl(10)-
          VLOP1=W1*W1T1S1
          VLOP0=0.D0
          IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
            CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
          IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
            CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
          IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
            CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
!         WL=-VLOP1*VINT_CI(LIST)
         IWDL=JUST(LRI,LRJ)
          IWDR=JUST(LRJ,LRJ)
C        CALL  PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)
        ENDIF
        VLOP1=W1*W1T1S2
        VLOP0=0.D0
        IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
          CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
        ENDIF
        IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
          CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
        ENDIF
        IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
          CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
        ENDIF
!         WL=-VLOP1*VINT_CI(LIST)
        DO LRK=NORB_FRZ+1,LRI-1
          LMK=LSM_INN(LRK)
!T1S(5-2)   (11)Ar(13)-Bl(32)-
          IF(JML.NE.MUL_TAB(LMK,LMI).OR.JMR.NE.
     :        MUL_TAB(LMK,LMJ)) CYCLE
!          VLOP1=W1*W1T1S2
!          WL=-VLOP1*VINT_CI(LIST)
           IWDL=JUST(LRK,LRI)
           IWDR=JUST(LRJ,LRK)
C         CALL  PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

        ENDDO
        DO LRK=LRI+1,LRJ-1
          LMK=LSM_INN(LRK)
!T1S(5-2)   Ar(13)-C'(11)-Bl(32)-
          IF(JML.NE.MUL_TAB(LMI,LMK).OR.JMR.NE.
     :        MUL_TAB(LMK,LMJ)) CYCLE
!          VLOP1=W1*W1T1S2
!          WL=-VLOP1*VINT_CI(LIST)
           IWDL=JUST(LRI,LRK)
           IWDR=JUST(LRJ,LRK)
C         CALL  PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,-WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,-WL2,JPER,LIST2,LIST3
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL1,JPER,LIST1)

        ENDDO
        VLOP1=W1*W1T1S3
        VLOP0=0.D0
        IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
          CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
        ENDIF
        IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
          CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
        ENDIF
        IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
          CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
        ENDIF
        DO LRK=LRJ+1,NORB_DZ
          LMK=LSM_INN(LRK)
!T1S(5-3)   Ar(13)-Bl(31)-C"(12)-
          IF(JML.NE.MUL_TAB(LMI,LMK).OR.JMR.NE.
     :        MUL_TAB(LMJ,LMK)) CYCLE
!          VLOP1=W1*W1T1S3
!          WL=-VLOP1*VINT_CI(LIST)
           IWDL=JUST(LRI,LRK)
           IWDR=JUST(LRK,LRJ)
C         CALL  PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

        ENDDO
        VLOP1=W1*W1T1S4
        VLOP0=0.D0
        IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
          CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
        ENDIF
        IF(LINE.EQ.28) THEN     !LRI,LRK,LRS,LRA
          CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
        ENDIF
        IF(LINE.EQ.29) THEN     !LRI,LRK,LRS,LRA
          CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
        ENDIF
!       WL=-VLOP1*VINT_CI(LIST)
        DO LRK=LRJ+1,NORB_DZ
          LMK=LSM_INN(LRK)
!T1S(5-4)   Ar(13)-Bl(32)-C"(11)-
          IF(JML.NE.MUL_TAB(LMI,LMK).OR.JMR.NE.
     :        MUL_TAB(LMJ,LMK)) CYCLE
           IWDL=JUST(LRI,LRK)
           IWDR=JUST(LRJ,LRK)
C         CALL  PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)
        ENDDO
        IF(LMIJ.EQ.JML.AND.LMIJ.EQ.JMR) THEN
!T1S(5-5)   (11)Drl(12)-
          VLOP1=W1*W1T1S5
          VLOP0=0.D0
          IF(LINE.EQ.26) THEN    !LRJ,LRA
            CALL COMP_LOOP_G(9,LRJ,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
          IF(LINE.EQ.28) THEN    !LRJ,LRS,LRA
            CALL COMP_LOOP_G(12,LRJ,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
          IF(LINE.EQ.29) THEN    !LRJ,LRS,LRA
            CALL COMP_LOOP_G(11,LRJ,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
!         WL=-VLOP1*VOINT(LRB,LRJ)
          IWDL=JUST(LRI,LRJ)
          IWDR=JUST(LRJ,LRI)
C         CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

!T1S(5-6)   Drl(11)-C"(12)-
          VLOP1=W1*W1T1S6
          VLOP0=0.D0
          IF(LINE.EQ.26) THEN    !LRI,LRA
            CALL COMP_LOOP_G(9,LRI,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
          IF(LINE.EQ.28) THEN    !LRI,LRS,LRA
            CALL COMP_LOOP_G(12,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
          IF(LINE.EQ.29) THEN    !LRI,LRS,LRA
            CALL COMP_LOOP_G(11,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
!         WL=-VLOP1*VOINT(LRB,LRI)
          IWDL=JUST(LRI,LRJ)
          IWDR=JUST(LRJ,LRI)
C         CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

!T1S(5-7)   Drl(12)-C"(11)-
          VLOP1=W1*W1T1S7
          VLOP0=0.D0
          IF(LINE.EQ.26) THEN    !LRI,LRA
            CALL COMP_LOOP_G(9,LRI,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
          IF(LINE.EQ.28) THEN    !LRI,LRS,LRA
            CALL COMP_LOOP_G(12,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
          IF(LINE.EQ.29) THEN    !LRI,LRS,LRA
            CALL COMP_LOOP_G(11,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
!         WL=-VLOP1*VOINT(LRB,LRI)
          IWDL=JUST(LRI,LRJ)
          IWDR=IWDL
C         CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

        ENDIF
      ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE TT_HEAD_DBL_TAIL_ACT_G(LRA)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      COMMON/ONEPL/LINE,JPH,JPEL,JPER,LRG,LRS,JWL,JWR,W0,W1
!TT(11-1) (22)Ar(23)-Bl(32)-
!TT(11-1) Ar(23)-C'(22)-Bl(32)-
!TT(11-1) Ar(23)-Bl(32)-C"(22)-
!TT(11-2) (22)Drl(22)-
!TT(11-2) Drl(22)-C"(22)-
!TT(11-3) (22)Drl(33)-
!TT(11-3) Drl(33)-C"(22)-
!TT(11-3) Drl(33)-C"(22)-C"(22)-
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          JMLR=MUL_TAB(JML,JMR)
          IF(LMIJ.NE.JMLR) CYCLE
          W0TT1=W0_TT(1)
          W1TT1=W1_TT(1)
          NI=MOD(LRJ-LRI,2)
          IF(NI.EQ.0) THEN
            W0TT1=-W0TT1
            W1TT1=-W1TT1
          ENDIF
!          LIST=LIST3(LRI,LRJ,LRB)
          VLOP0=W0*W0TT1
          VLOP1=W1*W1TT1
          IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
            CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
          IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
            CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
          IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
            CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
!          WL=(VLOP0-VLOP1)*vint_ci(LIST)-2*VLOP0*vint_ci(LIST+1)
!TT(11-1) (22)Ar(23)-Bl(32)-
          DO LRK=NORB_FRZ+1,LRI-1
            LMK=LSM_INN(LRK)
            LMKI=MUL_TAB(LMK,LMI)
            LMKJ=MUL_TAB(LMK,LMJ)
            IF(LMKI.EQ.JML.AND.LMKJ.EQ.JMR) THEN
              IWDL=JUST(LRK,LRI)     !
              IWDR=JUST(LRK,LRJ)     !
C              CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

             ENDIF
          ENDDO
!TT(11-1) Ar(23)-Bl(32)-C"(22)-    ACT -C"-
          DO LRK=LRJ+1,NORB_DZ
            LMK=LSM_INN(LRK)
            LMKI=MUL_TAB(LMK,LMI)
            LMKJ=MUL_TAB(LMK,LMJ)
            IF(LMKI.EQ.JML.AND.LMKJ.EQ.JMR) THEN
              IWDL=JUST(LRI,LRK)     !
              IWDR=JUST(LRJ,LRK)     !
C              CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

            ENDIF
          ENDDO
!TT(11-1) Ar(23)-C'(22)-Bl(32)-    ACT -C"-
          DO LRK=LRI+1,LRJ-1
           LMK=LSM_INN(LRK)
            LMKI=MUL_TAB(LMK,LMI)
            LMKJ=MUL_TAB(LMK,LMJ)
            IF(LMKI.EQ.JML.AND.LMKJ.EQ.JMR) THEN
              IWDL=JUST(LRI,LRK)      !
              IWDR=JUST(LRK,LRJ)      !
C              CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,-WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,-WL2,JPER,LIST2,LIST3
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL1,JPER,LIST1)

            ENDIF
          ENDDO
        ENDDO
      ENDDO

      IF(JPAD.NE.JPADL) RETURN
      IF(JWL .GE. JWR ) RETURN

      W0TT2=W0_TT(2)
      W1TT2=W1_TT(2)
      W0TT3=W0_TT(3)

      DO LRI=NORB_FRZ+1,NORB_DZ-1
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JML.OR.LMIJ.NE.JMR) CYCLE
          IF(JWL.GE.JWR) CYCLE
!TT(11-2) (22)Drl(22)-
!TT(11-2) Drl(22)-C"(22)-
          IWDL=JUST(LRI,LRJ)     !
          IWDR=IWDL
          VLOP0=W0*W0TT2
          VLOP1=W1*W1TT2
          IF(LINE.EQ.26) THEN    !LRJ,LRA
            CALL COMP_LOOP_G(9,LRJ,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
          IF(LINE.EQ.28) THEN    !LRJ,LRS,LRA
            CALL COMP_LOOP_G(12,LRJ,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
          IF(LINE.EQ.29) THEN    !LRJ,LRS,LRA
            CALL COMP_LOOP_G(11,LRJ,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF

C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

          IF(LINE.EQ.26) THEN    !LRI,LRA
            CALL COMP_LOOP_G(9,LRI,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
          IF(LINE.EQ.28) THEN    !LRI,LRS,LRA
            CALL COMP_LOOP_G(12,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
          IF(LINE.EQ.29) THEN    !LRI,LRS,LRA
            CALL COMP_LOOP_G(11,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)


C         WL=WL+WLTMP
!          WL=(VLOP0-VLOP1)*(VOINT(LRB,LRI)+VOINT(LRB,LRJ))
          VLOP0=W0*W0TT3
          VLOP1=0.D0
          DO LRK=1,NORB_DZ
            IF(LRK.EQ.LRI) CYCLE
            IF(LRK.EQ.LRJ) CYCLE
!TT(11-3) Drl(33)-C"(22)-C"(22)-
!TT(11-3) (22)Drl(33)-C"(22)-
!TT(11-3) (22)(22)Drl(33)-
            IF(LINE.EQ.26) THEN    !LRK,LRA
              CALL COMP_LOOP_G(9,LRK,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
            IF(LINE.EQ.28) THEN    !LRK,LRS,LRA
              CALL COMP_LOOP_G(12,LRK,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
            IF(LINE.EQ.29) THEN    !LRK,LRS,LRA
              CALL COMP_LOOP_G(11,LRK,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
C            WL=WL+WLTMP
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

          ENDDO
C          CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
        ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE TTTT_HEAD_DBL_TAIL_ACT_G(LRA)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      COMMON/ONEPL/LINE,JPH,JPEL,JPER,LRG,LRS,JWL,JWR,W0,W1
!T1T1(12-1)  Ar(13)-Bl(31)-
!T1T1(12-1)  Ar(13)-C'(11)-Bl(31)-
!T1T1(12-1)  Ar(13)-Bl(31)-C"(11)-
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        LMI=LSM_INN(LRI)
      DO LRJ=LRI+1,NORB_DZ
        LMJ=LSM_INN(LRJ)
        W0TT1=W0_T1T1(1)
        W1TT1=W1_T1T1(1)
        NI=MOD(LRJ-LRI,2)
        IF(NI.EQ.0) THEN
          W0TT1=-W0TT1
          W1TT1=-W1TT1
        ENDIF
!        LIST=LIST4(LRI,LRJ,LRA,LRB)
        VLOP0=W0*W0TT1
        VLOP1=W1*W1TT1
        IF(LINE.EQ.26) THEN    !LRI,LRJ,LRA
          CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
        ENDIF
        IF(LINE.EQ.28) THEN     !LRI,LRJ,LRS,LRA
          CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
        ENDIF
        IF(LINE.EQ.29) THEN     !LRI,LRJ,LRS,LRA
          CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
        ENDIF
!       WL=(VLOP0-VLOP1)*VINT_CI(LIST)-2*VLOP0*VINT_CI(LIST+1)
!T1T1(12-1)  (11)Ar(13)-Bl(31)-
        DO LRM=NORB_FRZ+1,LRI-1
          LMM=LSM_INN(LRM)
          IF(JML.NE.MUL_TAB(LMI,LMM).OR.JMR.NE.MUL_TAB(LMM,LMJ)) CYCLE
          IWDL=JUST(LRM,LRI)
          IWDR=JUST(LRM,LRJ)
C         CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

        ENDDO
!T1T1(12-1)  Ar(13)-Bl(31)-C"(11)-
        DO LRM=LRJ+1,NORB_DZ
          LMM=LSM_INN(LRM)
          IF(JML.NE.MUL_TAB(LMI,LMM).OR.JMR.NE.MUL_TAB(LMM,LMJ)) CYCLE
          IWDL=JUST(LRI,LRM)
          IWDR=JUST(LRJ,LRM)
C         CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

        ENDDO
!T1T1(12-1)  Ar(13)-C'(11)-Bl(31)-
        DO LRM=LRI+1,LRJ-1
          LMM=LSM_INN(LRM)
          IF(JML.NE.MUL_TAB(LMI,LMM).OR.JMR.NE.MUL_TAB(LMM,LMJ)) CYCLE
          IWDL=JUST(LRI,LRM)
          IWDR=JUST(LRM,LRJ)
C         CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,-WL,JPER)
C          CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,-WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL1,JPER,LIST1)

        ENDDO
      ENDDO
      ENDDO

      IF(JPAD.NE.JPADL) RETURN
      IF(JWL.GE.JWR) RETURN

      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
      DO LRJ=LRI+1,NORB_DZ               !BBS_TMP
        LMJ=LSM_INN(LRJ)
        LMIJ=MUL_TAB(LMI,LMJ)
        IF(JML.NE.LMIJ) CYCLE
        VLOP0=W0*W0_T1T1(2)
        VLOP1=W1*W1_T1T1(2)
C----------------------------------------
C lyb
        IWDL=JUST(LRI,LRJ)
        IWDR=IWDL

!T1T1(12-2)  (11)Drl(11)-
        IF(LINE.EQ.26) THEN    !LRI,LRA
           CALL COMP_LOOP_G(9,LRI,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
        ENDIF
        IF(LINE.EQ.28) THEN    !LRK,LRS,LRA
          CALL COMP_LOOP_G(12,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
        ENDIF
        IF(LINE.EQ.29) THEN    !LRK,LRS,LRA
          CALL COMP_LOOP_G(11,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
        ENDIF
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

!T1T1(12-2)  Drl(11)-C"(11)-
        IF(LINE.EQ.26) THEN    !LRJ,LRA
           CALL COMP_LOOP_G(9,LRJ,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
        ENDIF
        IF(LINE.EQ.28) THEN    !LRJ,LRS,LRA
          CALL COMP_LOOP_G(12,LRJ,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
        ENDIF
        IF(LINE.EQ.29) THEN    !LRJ,LRS,LRA
          CALL COMP_LOOP_G(11,LRJ,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
        ENDIF
C       WL=WL+WLTMP
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

!T1T1(12-3)  (11)(11)Drl(33)-
!T1T1(12-3)  (11)Drl(33)-C"(11)-
!T1T1(12-3)  Drl(33)-C"(11)-C"(11)-
        DO LRK=1,NORB_DZ
          IF(LRK.EQ.LRI) CYCLE
          IF(LRK.EQ.LRJ) CYCLE
          VLOP0=W0*W0_T1T1(3)
          VLOP1=0.D0
          IF(LINE.EQ.26) THEN    !LRK,LRA
            CALL COMP_LOOP_G(9,LRK,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
          IF(LINE.EQ.28) THEN    !LRK,LRS,LRA
            CALL COMP_LOOP_G(12,LRK,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
          IF(LINE.EQ.29) THEN    !LRK,LRS,LRA
            CALL COMP_LOOP_G(11,LRK,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
C          WL=WL+WLTMP
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

        ENDDO
C        IWDL=JUST(LRI,LRJ)
C        IWDR=IWDL
C       CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
      ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE DD_HEAD_DBL_TAIL_ACT_G(LRA)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      COMMON/ONEPL/LINE,JPH,JPEL,JPER,LRG,LRS,JWL,JWR,W0,W1
!DD(19-1) Ar(23)-Bl(32)-
!DD(19-2) Drl(22)-
!DD(19-3) Drl(33)-
!DD(19-3) Drl(33)-C"(22)-
      DO LRIL=NORB_FRZ+1,NORB_DZ
        IMIL=LSM_INN(LRIL)
        IF(IMIL.NE.JML) CYCLE
        IWDL=JUD(LRIL)
        DO LRIR=LRIL,NORB_DZ
          IMIR=LSM_INN(LRIR)
          IF(IMIR.NE.JMR) CYCLE
          IWDR=JUD(LRIR)

         W0DD1=W0_DD(1)
          W1DD1=W1_DD(1)
          NI=MOD(LRIR-LRIL,2)
          IF(NI.EQ.0) THEN
            W0DD1=-W0DD1
            W1DD1=-W1DD1
         ENDIF

         IF(LRIL.EQ.LRIR.and.JWL.LT.JWR) THEN
!DD(19-2) Drl(22)-
            VLOP0=W0*W0_DD(2)
            VLOP1=W1*W1_DD(2)
            IF(LINE.EQ.26) THEN    !LRIL,LRA
              CALL COMP_LOOP_G(9,LRIL,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
            IF(LINE.EQ.28) THEN    !LRIL,LRS,LRA
              CALL COMP_LOOP_G(12,LRIL,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
            IF(LINE.EQ.29) THEN    !LRIL,LRS,LRA
              CALL COMP_LOOP_G(11,LRIL,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

!            WL=(VLOP0-VLOP1)*VOINT(LRA,LRIL)
!            CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
!DD(19-3) Drl(33)-
            VLOP0=W0*W0_DD(3)
           VLOP1=0.D0
            DO LRK=1,NORB_DZ
             IF(LRK.EQ.LRIL) CYCLE
!             LIST=LIST3(LRS,LRA,LRK)
             IF(LINE.EQ.26) THEN    !LRK,LRA
               CALL COMP_LOOP_G(9,LRK,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)     !WYB
             ENDIF
             IF(LINE.EQ.28) THEN    !LRK,LRS,LRA
               CALL COMP_LOOP_G(12,LRK,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)   !WYB
             ENDIF
             IF(LINE.EQ.29) THEN    !LRK,LRS,LRA
               CALL COMP_LOOP_G(11,LRK,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)   !WYB
             ENDIF
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)
!              WL=WL+VLOP0*VOINT(LRK,LRA)
C             WL=WLTMP+WL
            ENDDO
C             CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
          ENDIF
          IF(LRIL.NE.LRIR)THEN
!DD(19-1) Ar(23)-Bl(32)-
            VLOP0=W0*W0DD1
            VLOP1=W1*W1DD1
            IF(LINE.EQ.26) THEN   !LRIL,LRIR,LRA
             CALL COMP_LOOP_G(5,LRIL,LRIR,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
            IF(LINE.EQ.28) THEN   !LRIL,LRIR,LRS,LRA
             CALL COMP_LOOP_G(7,LRIL,LRIR,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
            IF(LINE.EQ.29) THEN   !LRIL,LRIR,LRS,LRA
             CALL COMP_LOOP_G(6,LRIL,LRIR,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
!            LIST=LIST4(LRIL,LRIR,LRS,LRA)
!            WL=VLOP0*(vint_ci(LIST)-2*vint_ci(LIST+1)) !1.1       !!!!!
!     :       -VLOP1*vint_ci(LIST)
C            CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

          ENDIF
        ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE DDDD_HEAD_DBL_TAIL_ACT_G(LRA)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      COMMON/ONEPL/LINE,JPH,JPEL,JPER,LRG,LRS,JWL,JWR,W0,W1
!D1D1(20-1) Ar(13)-BL(31)-
!D1D1(20-1) Drl(11)-
!D1D1(20-1) Drl(33)-
!D1D1(20-1) Drl(33)-C"(11)-
      DO LRIL=NORB_FRZ+1,NORB_DZ
        IMIL=LSM_INN(LRIL)
        IF(IMIL.NE.JML) CYCLE
        IWDL=JUD(LRIL)
        DO LRIR=LRIL,NORB_DZ
          IMIR=LSM_INN(LRIR)
          IF(IMIR.NE.JMR) CYCLE
          W0DD1=W0_D1D1(1)
         W1DD1=W1_D1D1(1)
         NI=MOD(LRIR-LRIL,2)
         IF(NI.EQ.0) THEN
            W0DD1=-W0DD1
            W1DD1=-W1DD1
         ENDIF
         IWDR=JUD(LRIR)
          IF(LRIL.EQ.LRIR.and.JWL.LT.JWR)THEN
!D1D1(20-1) Drl(11)-
           VLOP0=W0*W0_D1D1(2)
            VLOP1=W1*W1_D1D1(2)
            IF(LINE.EQ.26) THEN    !LRIL,LRA
              CALL COMP_LOOP_G(9,LRIL,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
            IF(LINE.EQ.28) THEN    !LRIL,LRS,LRA
              CALL COMP_LOOP_G(12,LRIL,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
            IF(LINE.EQ.29) THEN    !LRIL,LRS,LRA
              CALL COMP_LOOP_G(11,LRIL,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

!D1D1(20-1) Drl(33)-
!D1D1(20-1) Drl(33)-C"(11)-
            VLOP0=W0*W0_D1D1(3)
           VLOP1=0.D0
           DO LRK=1,NORB_DZ
              IF(LRK.EQ.LRIL) CYCLE
              IF(LINE.EQ.26) THEN    !LRK,LRA
                CALL COMP_LOOP_G(9,LRK,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
              ENDIF
              IF(LINE.EQ.28) THEN    !LRK,LRS,LRA
                CALL COMP_LOOP_G(12,LRK,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
              ENDIF
              IF(LINE.EQ.29) THEN    !LRK,LRS,LRA
                CALL COMP_LOOP_G(11,LRK,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
              ENDIF
C             WL=WL+WLTMP
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

            ENDDO
C          CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
          ENDIF
          IF(LRIL.NE.LRIR)THEN
!D1D1(20-1) Ar(13)-BL(31)-
            VLOP0=W0*W0DD1
            VLOP1=W1*W1DD1
            IF(LINE.EQ.26) THEN   !LRIL,LRIR,LRA
             CALL COMP_LOOP_G(5,LRIL,LRIR,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
            IF(LINE.EQ.28) THEN   !LRIL,LRIR,LRS,LRA
             CALL COMP_LOOP_G(7,LRIL,LRIR,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
            IF(LINE.EQ.29) THEN   !LRIL,LRIR,LRS,LRA
             CALL COMP_LOOP_G(6,LRIL,LRIR,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
!            LIST=LIST3(LRIL,LRIR,LRA)
!            WL=(VLOP0-VLOP1)*vint_ci(LIST)-2*VLOP0*vint_ci(LIST+1)   !2
C            CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

         ENDIF
        ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE DD1_HEAD_DBL_TAIL_ACT_G(LRA)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      COMMON/ONEPL/LINE,JPH,JPEL,JPER,LRG,LRS,JWL,JWR,W0,W1
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        IWDL=JUD(LRI)
      DO LRJ=LRI+1,NORB_DZ
!DD1(21) Ar(23)-Bl(31)-
        LMJ=LSM_INN(LRJ)
        IF(JML.NE.LMI.OR.JMR.NE.LMJ) CYCLE
        VLOP1=W1*W1_DD1
        VLOP0=0.D0
        NI=MOD(LRJ-LRI,2)
        IF(NI.EQ.0) THEN
          VLOP1=-VLOP1
        ENDIF

        IF(LINE.EQ.26) THEN   !LRI,LRJ,LRA
          CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
        ENDIF
        IF(LINE.EQ.28) THEN   !LRIL,LRJ,LRS,LRA
          CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
        ENDIF
        IF(LINE.EQ.29) THEN   !LRI,LRJ,LRS,LRA
          CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
        ENDIF
!       LIST=LIST3(LRI,LRJ,LRA)
!       WL=-VLOP1*VINT_CI(LIST)
        IWDR=JUD(LRJ)
C       CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

      ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE D1D_HEAD_DBL_TAIL_ACT_G(LRA)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      COMMON/ONEPL/LINE,JPH,JPEL,JPER,LRG,LRS,JWL,JWR,W0,W1
!D1D(22-1)   Ar(13)-Bl(32)-
!D1D(22-2)   Drl(12)-
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        LMI=LSM_INN(LRI)
        IWDL=JUD(LRI)
      DO LRJ=LRI+1,NORB_DZ
        LMJ=LSM_INN(LRJ)
        IF(JML.NE.LMI.OR.JMR.NE.LMJ) CYCLE
        VLOP1=W1*W1_D1D(1)
        VLOP0=0.D0
        IF(MOD(LRJ-LRI,2).EQ.0) THEN
          VLOP1=-VLOP1
        ENDIF
        IF(LINE.EQ.26) THEN   !LRI,LRJ,LRA
          CALL COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
        ENDIF
        IF(LINE.EQ.28) THEN   !LRIL,LRJ,LRS,LRA
          CALL COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
        ENDIF
        IF(LINE.EQ.29) THEN   !LRI,LRJ,LRS,LRA
          CALL COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
        ENDIF
!       LIST=LIST3(LRI,LRJ,LRA)
!       WL=-VLOP1*VINT_CI(LIST)
        IWDR=JUD(LRJ)
C       CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

      ENDDO
      ENDDO
      VLOP1=W1*W1_D1D(2)
      VLOP0=0.D0
      DO LRI=NORB_FRZ+1,NORB_DZ
!D1D(22-2)   Drl(12)-
        LMI=LSM_INN(LRI)
        IF(JML.NE.LMI.OR.JMR.NE.LMI) CYCLE
        IF(LINE.EQ.26) THEN    !LRI,LRA
          CALL COMP_LOOP_G(9,LRI,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
        ENDIF
        IF(LINE.EQ.28) THEN    !LRI,LRS,LRA
          CALL COMP_LOOP_G(12,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
        ENDIF
        IF(LINE.EQ.29) THEN    !LRI,LRS,LRA
          CALL COMP_LOOP_G(11,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
        ENDIF
!       VLOP1=W1*W1_D1D(2)
!       WL=-VLOP1*VOINT(LRI,LRA)
        IWDL=JUD(LRI)
        IWDR=IWDL
C       CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

      ENDDO

      RETURN
      END

      SUBROUTINE SV_HEAD_DBL_TAIL_ACT_G(LRA)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      COMMON/ONEPL/LINE,JPH,JPEL,JPER,LRG,LRS,JWL,JWR,W0,W1
!SV(10-1) Ar(13)-Br(23)-
!SV(10-2) Ar(23)-Br(13)-
!SV(10-3) Drr(03)-
      IWDR=0
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
      DO LRJ=LRI,NORB_DZ
        LMJ=LSM_INN(LRJ)
        LMIJ=MUL_TAB(LMI,LMJ)
        IF(LMIJ.NE.JML) CYCLE
        W0SV1=W0_SV(1)
        W1SV1=W1_SV(1)
        W0SV2=W0_SV(2)
        W1SV2=W1_SV(2)
        NI=MOD(LRJ-LRI,2)
        IF(NI.EQ.0) THEN
          W0SV1=-W0SV1
          W1SV1=-W1SV1
          W0SV2=-W0SV2
          W1SV2=-W1SV2
        ENDIF
        IWDL=JUST(LRI,LRJ)
        IF(LRI.EQ.LRJ) THEN
          VLOP0=W0*W0_SV(3)            !D10-3
          VLOP1=0.D0
!         WL=VLOP0*VOINT(LRA,LRI)/2
          IF(LINE.EQ.25) THEN    !LRI,LRA
Cwyb         CALL COMP_LOOP_G(8,LRI,LRG,LRS,LRA,VLOP0,VLOP1,WL)
             CALL COMP_LOOP_G(8,LRI,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
          IF(LINE.EQ.27) THEN    !LRI,LRS,LRA
Cwyb         CALL COMP_LOOP_G(10,LRI,LRS,LRA,LRA,VLOP0,VLOP1,WL)
            CALL COMP_LOOP_G(10,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
C          CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

        ELSE
          VLOP0=W0*W0SV2             !D10-2
          VLOP1=W1*W1SV2
!         LIST=LIST3(LRI,LRJ,LRA)
!          WL=(VLOP0+VLOP1)*vint_ci(LIST)        !2.1          !!!!!
          IF(LINE.EQ.25) THEN    !LRI,LRJ,LRA
            CALL COMP_LOOP_G(3,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
          ENDIF
          IF(LINE.EQ.27) THEN      !LRI,LRJ,LRS,LRA
            CALL COMP_LOOP_G(4,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
         ENDIF
C          CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

          IF(JB_SYS.GT.0) THEN
            IWDL=JUST(LRJ,LRI)
            VLOP0=W0*W0SV1
            VLOP1=W1*W1SV1
!           LIST=LIST3(LRI,LRJ,LRA)
!            WL=(VLOP0+VLOP1)*vint_ci(LIST)        !2.1          !!!!!
            IF(LINE.EQ.25) THEN    !LRI,LRJ,LRA
              CALL COMP_LOOP_G(3,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
            IF(LINE.EQ.27) THEN      !LRI,LRJ,LRS,LRA
              CALL COMP_LOOP_G(4,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,
     :                     WL1,LIST1)
            ENDIF
C            CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
C           CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
              CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      IF(LIST1.NE.0)
     :        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

          ENDIF
        ENDIF
      ENDDO
      ENDDO
      RETURN

      END

      SUBROUTINE SD_HEAD_DBL_TAIL_ACT_G(LRA,LPCOE)
!**********************************************
!     LRA, ....partial loop.......
!**********************************************
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      DIMENSION LPCOE(NORB_DZ+1:NORB_INN)
      COMMON/ONEPL/LINE,JPH,JPEL,JPER,LRG,LRS,JWL,JWR,W0,W1
      JMLR=MUL_TAB(JML,JMR)
!SD(6-1) A&r(02)-
!SD(6-2) C(22)A&r(13)-
!SD(6-4) A&r(23)C'(12)-
!SD(6-5) A&r(23)B&r(13)B^r(32)
!SD(6-8) A&r(23)B&l(32)B^l(13)
!SD(6-9) D&r&r(03)B^r(32)
!SD(6-11) D&r&l(22)B^l(13)
!SD(6-12) D&r&l(33)B^l(02)
!SD(6-13) (22)D&r&l(33)B^l(13)
!SD(6-14) D&r&l(33)C"(22)B^l(13)
!SD(6-16) D&r&l(33)B^l(23)C'(12)

!SD(6-3) A&r(13)C'(22)-
!SD(6-6) A&r(13)B&r(23)B^r(32)
!SD(6-7) A&r(13)B&l(32)B^l(23)
!SD(6-10) D&r&l(12)B^l(23)
!SD(6-15) D&r&l(33)B^l(13)C'(22)

      IF(JML.NE.1) GOTO 207

!SD(6-1) A&r(02)-
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        IWDL=JUST(LRI,LRI)
        W0SD1 =W0_SD(1)
        W0SD12=W0_SD(12)
        NI=MOD(NORB_DZ-LRI,2)
        IF(NI.EQ.1) THEN
          W0SD1 =-W0SD1
          W0SD12=-W0SD12
        ENDIF
        IF(LMI.EQ.JMLR) THEN
          IWDR=JUD(LRI)
          VLOP0=W0*W0SD1
C          WL=VLOP0*VOINT(LRI,LRA)
          WL=VLOP0
          CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,LRI,LRA)

          DO LR=LRI+1,NORB_DZ
C           LIST =LIST3(LRI,LRA,LR)
C            WL=WL+VLOP0*(2*VINT_CI(LIST+1)-VINT_CI(LIST)) !  310:NEOC=2
          WL=VLOP0*2
          CALL TRANS_IJKL_INTPOS(LRA,LRI,LR,LR,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
          WL=-VLOP0
          CALL TRANS_IJKL_INTPOS(LRA,LR,LRI,LR,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
          ENDDO
          DO LRK=norb_dz+1,LRA
C            LIST=LIST3(LRI,LRA,LRK)
            KCOE=LPCOE(LRK)
            CALL NEOC(KCOE,NOCC,TCOE)
C            WL=WL+VLOP0*NOCC*(VINT_CI(LIST+1)+TCOE*VINT_CI(LIST))
            WL=VLOP0*NOCC
          CALL TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=VLOP0*NOCC*TCOE
          CALL TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

          ENDDO
C          WL=WL*VLOP0
!SD(6-12) D&rl(33)B^l(02)
          VLOP0=W0*W0SD12
         DO LRK=1,LRI-1
C            LIST=LIST3(LRI,LRA,LRK)
C            WL=WL+VLOP0*(vint_ci(LIST)-2*VINT_CI(LIST+1))
            WL=VLOP0
          CALL TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
           WL=-VLOP0*2
          CALL TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

          ENDDO
C          CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
        ENDIF

!SD(6-9) D&r&r(03)B^r(32)
        DO LRD=LRI+1,NORB_DZ
         LMD=LSM_INN(LRD)
          IF(LMD.NE.JMR) CYCLE
          W0SD9=W0_SD(9)
          NI=MOD(NORB_DZ-LRD,2)
          IF(NI.EQ.1)   W0SD9=-W0SD9
          IWDR=JUD(LRD)
          VLOP0=W0*W0SD9
C          LIST=LIST3(LRD,LRA,LRI)
C          WL=VINT_CI(LIST)*VLOP0
           WL=VLOP0
          CALL TRANS_IJKL_INTPOS(LRA,LRI,LRD,LRI,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
        ENDDO
      ENDDO

207   DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JML) CYCLE
          W0SD2 =W0_SD(2)
          W1SD2 =W1_SD(2)
          W0SD11=W0_SD(11)
          W1SD11=W1_SD(11)
          W0SD14=W0_SD(14)
          W1SD14=0.D0
          NI=MOD(NORB_DZ-LRJ,2)
          IF(NI.EQ.1) THEN
            W0SD2 =-W0SD2
            W1SD2 =-W1SD2
            W0SD11=-W0SD11
            W1SD11=-W1SD11
            W0SD14=-W0SD14
          ENDIF
          IWDL=JUST(LRI,LRJ)
          IWDL1=JUST(LRJ,LRI)
C**********************************************************
!SD(6-2) C(22)-A&r(13)-
          IF(LMI.EQ.JMR) THEN
            IWDR=JUD(LRI)
            VLOP0=W0*W0SD2
C            LIST=LIST3(LRJ,LRA,LRJ)
C            WL=VLOP0*(VOINT(LRJ,LRA)+VINT_CI(LIST))          !310,act_c
            WL=VLOP0
          CALL TRANS_IJKL_INTPOS(LRA,LRJ,LRJ,LRJ,NXO)
          CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,LRJ,LRA)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

            DO LR=LRJ+1,NORB_DZ
C              LIST =LIST3(LRJ,LRA,LR)
C              WL=WL+VLOP0*(2*VINT_CI(LIST+1)-VINT_CI(LIST)) !  310:NEOC
            WL=VLOP0*2
          CALL TRANS_IJKL_INTPOS(LRA,LRJ,LR,LR,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=-VLOP0
          CALL TRANS_IJKL_INTPOS(LRA,LR,LRJ,LR,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

            ENDDO
            DO LRK=norb_dz+1,LRA
C              LIST=LIST3(LRJ,LRA,LRK)
              KCOE=LPCOE(LRK)
              CALL NEOC(KCOE,NOCC,TCOE)
C              WL=WL+VLOP0*NOCC*(VINT_CI(LIST+1)+TCOE*VINT_CI(LIST))
            WL=VLOP0*NOCC
          CALL TRANS_IJKL_INTPOS(LRA,LRJ,LRK,LRK,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=VLOP0*NOCC*TCOE
          CALL TRANS_IJKL_INTPOS(LRA,LRK,LRJ,LRK,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

            ENDDO
C            WL=WL*VLOP0
!SD(6-11) D&r&l(22)B^l(13)
            VLOP0=W0*W0SD11
            VLOP1=W1*W1SD11
C            LIST=LIST3(LRJ,LRA,LRI)
C          WL=WL+(VLOP0-VLOP1)*vint_ci(LIST)-2*VLOP0*VINT_CI(LIST+1)
            WL=VLOP0-VLOP1
          CALL TRANS_IJKL_INTPOS(LRA,LRI,LRJ,LRI,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
           WL=-VLOP0*2
          CALL TRANS_IJKL_INTPOS(LRA,LRJ,LRI,LRI,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

!SD(6-13) (22)D&r&l(33)B^l(13)
!SD(6-14) D&r&l(33)C"(22)B^l(13)
            VLOP0=W0*W0SD14
            DO LRK=1,LRJ-1
              IF(LRK.EQ.LRI) CYCLE
C              LIST=LIST3(LRJ,LRA,LRK)
C              WL=WL+VLOP0*(vint_ci(LIST)-2*VINT_CI(LIST+1))
            WL=VLOP0
          CALL TRANS_IJKL_INTPOS(LRA,LRK,LRJ,LRK,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
           WL=-VLOP0*2
          CALL TRANS_IJKL_INTPOS(LRA,LRJ,LRK,LRK,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

            ENDDO
C            CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
          ENDIF
        IF(JB_SYS.GT.0) THEN
          IF(LMJ.EQ.JMR) THEN
            IWDR=JUD(LRJ)
            W0SD3=W0_SD(3)
            W0SD15=W0_SD(15)
            W1SD10=W1_SD(10)
           NI=MOD(NORB_DZ-LRI,2)
            IF(NI.EQ.0) THEN
             W0SD3=-W0SD3
             W0SD15=-W0SD15
              W1SD10=-W1SD10
            ENDIF
!SD(6-3) A&r(13)-C'(22)-
           VLOP0=W0*W0SD3
C            LIST=LIST3(LRI,LRA,LRI)
C            WL=VOINT(LRI,LRA)+VINT_CI(LIST)
          WL=VLOP0
          CALL PRODAB_1(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,LRI,LRA)
          CALL TRANS_IJKL_INTPOS(LRA,LRI,LRI,LRI,NXO)
          CALL PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)
C            LIST=LIST3(LRI,LRA,LRJ)
C            WL=WL+VINT_CI(LIST+1)-VINT_CI(LIST)   !310 C'(22)NEOC=1,COE
            WL=VLOP0
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRJ,LRJ,NXO)
            CALL PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=-VLOP0
            CALL TRANS_IJKL_INTPOS(LRA,LRJ,LRI,LRJ,NXO)
            CALL PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)

            DO LR=LRI+1,NORB_DZ
              IF(LR.EQ.LRJ) CYCLE
C           LIST =LIST3(LRI,LRA,LR)
C              WL=WL+2*VINT_CI(LIST+1)-VINT_CI(LIST)       !310:NEOC=2,C
            WL=VLOP0*2
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LR,LR,NXO)
            CALL PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=-VLOP0
            CALL TRANS_IJKL_INTPOS(LRA,LR,LRI,LR,NXO)
            CALL PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)

            ENDDO
            DO LRK=norb_dz+1,LRA
C              LIST=LIST3(LRI,LRA,LRK)
              KCOE=LPCOE(LRK)
              CALL NEOC(KCOE,NOCC,TCOE)
C              WL=WL+NOCC*(VINT_CI(LIST+1)+TCOE*VINT_CI(LIST))
            WL=VLOP0*NOCC
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
            CALL PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=VLOP0*NOCC*TCOE
            CALL TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
            CALL PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)

            ENDDO
C            WL=WL*VLOP0
!SD(6-15) D&r&l(33)B^l(13)C'(22)
            VLOP0=W0*W0SD15
            DO LRK=1,LRI-1
C              LIST=LIST3(LRI,LRA,LRK)
C              WL=WL-VLOP0*(2*VINT_CI(LIST+1)-vint_ci(LIST))
            WL=-VLOP0*2
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
            CALL PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=VLOP0
            CALL TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
            CALL PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)

            ENDDO
C            CALL PRODAB(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER)
          ENDIF
         IF(LMIJ.EQ.JML.AND.LMI.EQ.JMR) THEN
            IWDR=JUD(LRI)
            IWDL1=JUST(LRJ,LRI)
            W1SD10=W1_SD(10)
            IF(MOD(NORB_DZ-LRJ,2).EQ.1) THEN
              W1SD10=-W1SD10
            ENDIF
!SD(6-10) D&r&l(12)B^l(23)
            VLOP1=W1*W1SD10
C            LIST=LIST3(LRJ,LRA,LRI)
C          WL=-VLOP1*vint_ci(LIST)      !4.3
            WL=-VLOP1
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRJ,LRI,NXO)
            CALL PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)
          ENDIF
        ENDIF
!SD(6-4) A&r(23)-C'(12)-
          IF(LMJ.EQ.JMR) THEN
            IWDR=JUD(LRJ)
            W0SD4=W0_SD(4)
            W0SD16=W0_SD(16)
           NI=MOD(NORB_DZ-LRI,2)
            IF(NI.EQ.0) W0SD4=-W0SD4
            IF(NI.EQ.0) W0SD16=-W0SD16
            VLOP0=W0*W0SD4
C            LIST=LIST3(LRI,LRA,LRI)
C            WL=VOINT(LRI,LRA)+VINT_CI(LIST)             !310,act_coe,61
            WL=VLOP0
          CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,LRI,LRA)
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRI,LRI,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

C            LIST=LIST3(LRI,LRA,LRJ)
C            WL=WL+VINT_CI(LIST+1)-(JB_SYS+2)*1.d0*VINT_CI(LIST)
            WL=VLOP0
              CALL TRANS_IJKL_INTPOS(LRA,LRI,LRJ,LRJ,NXO)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=-VLOP0*(JB_SYS+2)
              CALL TRANS_IJKL_INTPOS(LRA,LRJ,LRI,LRJ,NXO)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

            DO LR=LRI+1,NORB_DZ
              IF(LR.EQ.LRJ) CYCLE
C           LIST =LIST3(LRI,LRA,LR)
C              WL=WL+2*VINT_CI(LIST+1)-VINT_CI(LIST)       !310:NEOC=2,C
            WL=VLOP0*2
              CALL TRANS_IJKL_INTPOS(LRA,LRI,LR,LR,NXO)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=-VLOP0
              CALL TRANS_IJKL_INTPOS(LRA,LR,LRI,LR,NXO)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

            ENDDO
            DO LRK=norb_dz+1,LRA
C              LIST=LIST3(LRI,LRA,LRK)
              KCOE=LPCOE(LRK)
              CALL NEOC(KCOE,NOCC,TCOE)
C              WL=WL+NOCC*(VINT_CI(LIST+1)+TCOE*VINT_CI(LIST))
            WL=VLOP0*NOCC
              CALL TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=VLOP0*NOCC*TCOE
              CALL TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

            ENDDO
C            WL=WL*VLOP0
!SD(6-16) D&r&l(33)B^l(23)C'(12)
            VLOP0=W0*W0SD16
            DO LRK=1,LRI-1
C              LIST=LIST3(LRI,LRA,LRK)
C              WL=WL-VLOP0*(2*VINT_CI(LIST+1)-vint_ci(LIST))
            WL=-VLOP0*2
              CALL TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=VLOP0
              CALL TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

            ENDDO
C            CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
          ENDIF

!SD(6-5) A&r(23)B&r(13)B^r(32)
          DO LRD=LRJ+1,NORB_DZ
           LMD=LSM_INN(LRD)
            IF(LMD.NE.JMR) CYCLE
           IWDR=JUD(LRD)
            W0SD5=W0_SD(5)
            W1SD5=W1_SD(5)
            NI=MOD(LRJ-LRI+NORB_DZ-LRD,2)
            IF(NI.EQ.0)   W0SD5=-W0SD5
            IF(NI.EQ.0)   W1SD5=-W1SD5
            VLOP0=W0*W0SD5
            VLOP1=W1*W1SD5
C            LIST=LIST4(LRI,LRJ,LRD,LRA)
C            WL=(VLOP0-VLOP1)*vint_ci(LIST+2)+
C     :         (VLOP0+VLOP1)*vint_ci(LIST)             !1.3
            WL=VLOP0-VLOP1
              CALL TRANS_IJKL_INTPOS(LRA,LRI,LRJ,LRD,NXO)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=VLOP0+VLOP1
              CALL TRANS_IJKL_INTPOS(LRA,LRJ,LRD,LRI,NXO)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

C         CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
         ENDDO
        IF(JB_SYS.GT.0) THEN
!SD(6-6) A&r(13)B&r(23)B^r(32)
          DO LRD=LRJ+1,NORB_DZ
           LMD=LSM_INN(LRD)
            IF(LMD.NE.JMR) CYCLE
           IWDR=JUD(LRD)
            W0SD6=W0_SD(6)
            W1SD6=W1_SD(6)
            NI=MOD(LRJ-LRI+NORB_DZ-LRD,2)
            IF(NI.EQ.0)   W0SD6=-W0SD6
            IF(NI.EQ.0)   W1SD6=-W1SD6
            VLOP0=W0*W0SD6
            VLOP1=W1*W1SD6
C            LIST=LIST4(LRI,LRJ,LRD,LRA)
C            WL=(VLOP0-VLOP1)*vint_ci(LIST+2)+
C     :         (VLOP0+VLOP1)*vint_ci(LIST)             !1.3
            WL=VLOP0-VLOP1
              CALL TRANS_IJKL_INTPOS(LRA,LRI,LRJ,LRD,NXO)
            CALL PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=VLOP0+VLOP1
              CALL TRANS_IJKL_INTPOS(LRA,LRJ,LRD,LRI,NXO)
            CALL PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)

C         CALL PRODAB(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER)
         ENDDO
!SD(6-7) A&r(13)B&l(32)B^l(23)
         DO LRD=LRI+1,LRJ-1
           LMD=LSM_INN(LRD)
           IF(LMD.NE.JMR) CYCLE
           IWDR=JUD(LRD)
           W0SD7=W0_SD(7)
           W1SD7=W1_SD(7)
           NI=MOD(LRD-LRI+NORB_DZ-LRJ,2)
           IF(NI.EQ.0)   W0SD7=-W0SD7
           IF(NI.EQ.0)   W1SD7=-W1SD7
           VLOP0=W0*W0SD7
           VLOP1=W1*W1SD7
C           LIST=LIST4(LRI,LRD,LRJ,LRA)
C          WL=(VLOP0-VLOP1)*vint_ci(LIST+2)-
C     :             2*VLOP0*vint_ci(LIST+1)      !1.2
            WL=VLOP0-VLOP1
              CALL TRANS_IJKL_INTPOS(LRA,LRI,LRD,LRJ,NXO)
            CALL PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=-2*VLOP0
              CALL TRANS_IJKL_INTPOS(LRA,LRJ,LRD,LRI,NXO)
            CALL PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)
C          CALL PRODAB(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER)
         ENDDO
        ENDIF
!SD(6-8) A&r(23)B&l(32)B^l(13)
         DO LRD=LRI+1,LRJ-1
           LMD=LSM_INN(LRD)
           IF(LMD.NE.JMR) CYCLE
           IWDR=JUD(LRD)
           W0SD8=W0_SD(8)
           W1SD8=W1_SD(8)
           NI=MOD(LRD-LRI+NORB_DZ-LRJ,2)
           IF(NI.EQ.0)   W0SD8=-W0SD8
           IF(NI.EQ.0)   W1SD8=-W1SD8
           VLOP0=W0*W0SD8
           VLOP1=W1*W1SD8
C           LIST=LIST4(LRI,LRD,LRJ,LRA)
C          WL=(VLOP0-VLOP1)*vint_ci(LIST+2)-
C     :             2*VLOP0*vint_ci(LIST+1)      !1.2
            WL=VLOP0-VLOP1
              CALL TRANS_IJKL_INTPOS(LRA,LRI,LRD,LRJ,NXO)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=-2*VLOP0
              CALL TRANS_IJKL_INTPOS(LRA,LRJ,LRD,LRI,NXO)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

C        CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
        ENDDO
      ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE SDD_HEAD_DBL_TAIL_ACT_G(LRA,LPCOE)
!**********************************************
!     LRA, ....partial loop.......
!**********************************************
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      DIMENSION LPCOE(NORB_DZ+1:NORB_INN)
      COMMON/ONEPL/LINE,JPH,JPEL,JPER,LRG,LRS,JWL,JWR,W0,W1
      JMLR=MUL_TAB(JML,JMR)
!SD1(8-1) A&r(01)-
!SD1(8-2) C(11)A&r(23)-
!SD1(8-3) A&r(13)C'(21)-
!SD1(8-4) A&r(23)C'(11)-
!SD1(8-5) A&r(13)B&r(23)B^r(31)
!SD1(8-6) A&r(23)B&r(13)B^r(31)
!SD1(8-7) A&r(13)B&l(31)B^l(23)
!SD1(8-8) A&r(23)B&l(31)B^l(13)
!SD1(8-9) D&r&r(03)B^r(31)
!SD1(8-10) D&r&l(11)B^l(23)
!SD1(8-11) D&r&l(33)B^l(01)
!SD1(8-12) D&r&l(33)B^l(23)
!SD1(8-13) D&r&l(33)C"(13)B^l(23)
!SD1(8-14) D&r&l(33)B^l(11)C'(23)
!SD1(8-15) D&r&l(33)B^l(23)C'(11)
      IF(JML.NE.1) GOTO 209
      DO LRI=NORB_FRZ+1,NORB_DZ
!SD1(8-1) A&r(01)-
        LMI=LSM_INN(LRI)
        IWDL=JUST(LRI,LRI)
        W0SD1=W0_SD1(1)
        W0SD11=W0_SD1(9)
        NI=MOD(NORB_DZ-LRI,2)
        IF(NI.EQ.1) THEN
          W0SD1 =-W0SD1
          W0SD11=-W0SD11
        ENDIF
        IF(LMI.EQ.JMLR) THEN
          IWDR=JUD(LRI)
          VLOP0=W0*W0SD1
C          WL=VLOP0*VOINT(LRI,LRA)
          WL=VLOP0
          CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,LRI,LRA)

          DO LR=LRI+1,NORB_DZ
C            LIST =LIST3(LRI,LRA,LR)
C            WL=WL+VLOP0*(2*VINT_CI(LIST+1)-VINT_CI(LIST)) !  310:NEOC=2
          WL=VLOP0*2
          CALL TRANS_IJKL_INTPOS(LRA,LRI,LR,LR,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
          WL=-VLOP0
          CALL TRANS_IJKL_INTPOS(LRA,LR,LRI,LR,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

          ENDDO
          DO LRK=norb_dz+1,LRA
C            LIST=LIST3(LRI,LRA,LRK)
            KCOE=LPCOE(LRK)
            CALL NEOC(KCOE,NOCC,TCOE)
C            WL=WL+VLOP0*NOCC*(VINT_CI(LIST+1)+TCOE*VINT_CI(LIST))
            WL=VLOP0*NOCC
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=VLOP0*NOCC*TCOE
            CALL TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

          ENDDO
C          WL=WL*VLOP0
!SD1(8-11) D&rl(33)B^l(01)
          VLOP0=W0*W0SD11
         DO LRK=1,LRI-1
C            LIST=LIST3(LRI,LRA,LRK)
C            WL=WL+VLOP0*(vint_ci(LIST)-2*VINT_CI(LIST+1))
            WL=VLOP0
            CALL TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=-VLOP0*2
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

          ENDDO
C          CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
        ENDIF

!SD1(8-9) D&r&r(03)B^r(31)
        DO LRD=LRI+1,NORB_DZ
         LMD=LSM_INN(LRD)
          IF(LMD.NE.JMR) CYCLE
          W0SD9=W0_SD1(9)
          NI=MOD(NORB_DZ-LRD,2)
          IF(NI.EQ.1)   W0SD9=-W0SD9
          IWDR=JUD(LRD)
          VLOP0=W0*W0SD9
C          LIST=LIST3(LRD,LRA,LRI)
C          WL=VINT_CI(LIST)*VLOP0
            WL=VLOP0
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRD,LRI,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

C          CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
        ENDDO
      ENDDO

209   DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          IF(LRI.EQ.LRJ) CYCLE
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JML) CYCLE
          W0SD2 =W0_SD1(2)
          W1SD2 =W1_SD1(2)
          W0SD3 =W0_SD1(3)
          W1SD3 =W1_SD1(3)
          W0SD4 =W0_SD1(4)
          W1SD4 =W1_SD1(4)
          W0SD10=W0_SD1(10)
          W1SD10=W1_SD1(10)
          W0SD11=W0_SD1(11)
          W1SD11=0.D0
          W0SD12=W0_SD1(12)
          W0SD13=W0_SD1(13)
          NI=MOD(NORB_DZ-LRJ,2)
          IF(NI.EQ.1) THEN
            W0SD2 =-W0SD2
            W1SD2 =-W1SD2
            W0SD10=-W0SD10
            W1SD10=-W1SD10
            W0SD11=-W0SD11
          ENDIF
          IF(MOD(NORB_DZ-LRI,2).EQ.1) THEN
            W0SD3 =-W0SD3
            W1SD3 =-W1SD3
            W0SD4 =-W0SD4
            W1SD4 =-W1SD4
            W0SD12=-W0SD12
            W0SD13=-W0SD13
          ENDIF
          IWDL=JUST(LRJ,LRI)
          IWDL1=JUST(LRI,LRJ)
C**********************************************************
!SD1(8-2) C(11)-A&r(23)-
          IF(LMI.EQ.JMR) THEN
            IWDL=JUST(LRJ,LRI)
            IWDR=JUD(LRI)
            VLOP0=W0*W0SD2
C            LIST=LIST3(LRJ,LRA,LRJ)
C            WL=VLOP0*(VOINT(LRJ,LRA)+VINT_CI(LIST))
              WL=VLOP0
        CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,LRJ,LRA)
            CALL TRANS_IJKL_INTPOS(LRA,LRJ,LRJ,LRJ,NXO)
        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

            DO LR=LRJ+1,NORB_DZ
C              LIST =LIST3(LRJ,LRA,LR)
C              WL=WL+VLOP0*(2*VINT_CI(LIST+1)-VINT_CI(LIST)) !  310:NEOC
            WL=VLOP0*2
            CALL TRANS_IJKL_INTPOS(LRA,LRJ,LR,LR,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=-VLOP0
            CALL TRANS_IJKL_INTPOS(LRA,LR,LRJ,LR,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

            ENDDO
            DO LRK=norb_dz+1,LRA
C              LIST=LIST3(LRJ,LRA,LRK)
              KCOE=LPCOE(LRK)
              CALL NEOC(KCOE,NOCC,TCOE)
C              WL=WL+VLOP0*NOCC*(VINT_CI(LIST+1)+TCOE*VINT_CI(LIST))
            WL=VLOP0*NOCC
            CALL TRANS_IJKL_INTPOS(LRA,LRJ,LRK,LRK,NXO)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=VLOP0*NOCC*TCOE
            CALL TRANS_IJKL_INTPOS(LRA,LRK,LRJ,LRK,NXO)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

            ENDDO
C            WL=WL*VLOP0
!SD1(8-10) D&r&l(11)B^l(23)
            VLOP0=W0*W0SD10
            VLOP1=W1*W1SD10
C            LIST=LIST3(LRJ,LRA,LRI)
C          WL=WL+(VLOP0-VLOP1)*vint_ci(LIST)-2*VLOP0*VINT_CI(LIST+1)

            WL=VLOP0-VLOP1
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRJ,LRI,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=-VLOP0*2
            CALL TRANS_IJKL_INTPOS(LRA,LRJ,LRI,LRI,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)


!SD1(8-11) (11)D&r&l(33)B^l(23)
!SD1(8-11) D&r&l(33)C"(11)B^l(23)
            VLOP0=W0*W0SD11
            DO LRK=1,LRJ-1
              IF(LRK.EQ.LRI) CYCLE
C              LIST=LIST3(LRJ,LRA,LRK)
C              WL=WL+VLOP0*(vint_ci(LIST)-2*VINT_CI(LIST+1))
            WL=VLOP0
            CALL TRANS_IJKL_INTPOS(LRA,LRK,LRJ,LRK,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=-VLOP0*2
            CALL TRANS_IJKL_INTPOS(LRA,LRJ,LRK,LRK,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

            ENDDO
C            CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
          ENDIF
!SD1(8-3) A&r(13)-C'(21)-
          IF(LMJ.EQ.JMR) THEN
            VLOP0=-W0*W0SD3
C--------------------------------------
C lyb
           IWDL=JUST(LRJ,LRI)
           IWDR=JUD(LRJ)

C            LIST=LIST3(LRI,LRA,LRI)
C            WL=VOINT(LRI,LRA)+VINT_CI(LIST)
              WL=VLOP0
C            write(nf2,'(a3,5I8)') 'sdd',JPEL,IWDL,IWDR,JWL,JWR
        CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,LRI,LRA)
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRI,LRI,NXO)
        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

C            LIST=LIST3(LRI,LRA,LRJ)
C            WL=WL+VINT_CI(LIST+1)+JB_SYS*1.d0*VINT_CI(LIST)
            WL=VLOP0
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRJ,LRJ,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=VLOP0*JB_SYS
            CALL TRANS_IJKL_INTPOS(LRA,LRJ,LRI,LRJ,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

            DO LR=LRI+1,NORB_DZ
              IF(LR.EQ.LRJ) CYCLE
C           LIST =LIST3(LRI,LRA,LR)
C              WL=WL+2*VINT_CI(LIST+1)-VINT_CI(LIST)       !310:NEOC=2,C
            WL=VLOP0*2
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LR,LR,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=-VLOP0
            CALL TRANS_IJKL_INTPOS(LRA,LR,LRI,LR,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

            ENDDO
            DO LRK=norb_dz+1,LRA
C              LIST=LIST3(LRI,LRA,LRK)
              KCOE=LPCOE(LRK)
              CALL NEOC(KCOE,NOCC,TCOE)
C              WL=WL+NOCC*(VINT_CI(LIST+1)+TCOE*VINT_CI(LIST))
            WL=VLOP0*NOCC
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=VLOP0*NOCC*TCOE
            CALL TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            ENDDO
C            WL=WL*VLOP0
!SD1(8-12) Drl(33)-BL(13)-C'(21)-
C            IWDL=JUST(LRJ,LRI)
C           IWDR=JUD(LRJ)
            VLOP0=-W0*W0SD12
            DO LRK=1,LRI-1
C              LIST=LIST3(LRI,LRA,LRK)
C              WL=WL-VLOP0*(2*VINT_CI(LIST+1)-vint_ci(LIST))
            WL=-VLOP0*2
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=VLOP0
            CALL TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

            ENDDO
C            CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
          ENDIF
!SD1(8-4) A&r(23)-C'(11)-
          IF(LMJ.EQ.JMR) THEN
C-----------------------------
C lyb
           IWDL=JUST(LRI,LRJ)
           IWDR=JUD(LRJ)
             VLOP0=-W0*W0SD4
C            LIST=LIST3(LRI,LRA,LRI)
C            WL=VOINT(LRI,LRA)+VINT_CI(LIST)
              WL=VLOP0
        CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,LRI,LRA)
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRI,LRI,NXO)
        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

C            LIST=LIST3(LRI,LRA,LRJ)
C            WL=WL+VINT_CI(LIST+1)-VINT_CI(LIST)
            WL=VLOP0
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRJ,LRJ,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=-VLOP0
            CALL TRANS_IJKL_INTPOS(LRA,LRJ,LRI,LRJ,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

            DO LR=LRI+1,NORB_DZ
              IF(LR.EQ.LRJ) CYCLE
C           LIST =LIST3(LRI,LRA,LR)
C              WL=WL+2*VINT_CI(LIST+1)-VINT_CI(LIST)       !310:NEOC=2,C
            WL=VLOP0*2
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LR,LR,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=-VLOP0
            CALL TRANS_IJKL_INTPOS(LRA,LR,LRI,LR,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

            ENDDO
            DO LRK=norb_dz+1,LRA
C              LIST=LIST3(LRI,LRA,LRK)
              KCOE=LPCOE(LRK)
              CALL NEOC(KCOE,NOCC,TCOE)
C              WL=WL+NOCC*(VINT_CI(LIST+1)+TCOE*VINT_CI(LIST))
            WL=VLOP0*NOCC
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=VLOP0*NOCC*TCOE
            CALL TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

            ENDDO
C            WL=WL*VLOP0
!SD1(8-13) Drl(33)-BL(23)-C'(11)-
C            IWDL=JUST(LRI,LRJ)
C           IWDR=JUD(LRJ)
            VLOP0=-W0*W0SD13
            DO LRK=1,LRI-1
C              LIST=LIST3(LRI,LRA,LRK)
C              WL=WL-VLOP0*(2*VINT_CI(LIST+1)-vint_ci(LIST))
            WL=-VLOP0*2
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=VLOP0
            CALL TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

            ENDDO
C            CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
          ENDIF
!SD1(8-5) A&r(13)B&r(23)B^r(31)
          DO LRD=LRJ+1,NORB_DZ
           LMD=LSM_INN(LRD)
            IF(LMD.NE.JMR) CYCLE
            IWDL=JUST(LRJ,LRI)
           IWDR=JUD(LRD)
            W0SD5=W0_SD1(5)
            W1SD5=W1_SD1(5)
            NI=MOD(LRJ-LRI+NORB_DZ-LRD,2)
            IF(NI.EQ.0)   W0SD5=-W0SD5
            IF(NI.EQ.0)   W1SD5=-W1SD5
            VLOP0=W0*W0SD5
            VLOP1=W1*W1SD5
C            LIST=LIST4(LRI,LRJ,LRD,LRA)
C            WL=(VLOP0-VLOP1)*vint_ci(LIST+2)+
C     :         (VLOP0+VLOP1)*vint_ci(LIST)             !1.3
            WL=VLOP0-VLOP1
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRJ,LRD,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=VLOP0+VLOP1
            CALL TRANS_IJKL_INTPOS(LRA,LRJ,LRD,LRI,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
C          CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
!SD1(8-6) A&r(23)B&r(13)B^r(31)
            IWDL=JUST(LRI,LRJ)
            W0SD6=W0_SD1(6)
            W1SD6=W1_SD1(6)
            NI=MOD(LRJ-LRI+NORB_DZ-LRD,2)
            IF(NI.EQ.0)   W0SD6=-W0SD6
            IF(NI.EQ.0)   W1SD6=-W1SD6
            VLOP0=W0*W0SD6
            VLOP1=W1*W1SD6
C            LIST=LIST4(LRI,LRJ,LRD,LRA)
C            WL=(VLOP0-VLOP1)*vint_ci(LIST+2)+
C     :         (VLOP0+VLOP1)*vint_ci(LIST)             !1.3
            WL=VLOP0-VLOP1
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRJ,LRD,NXO)
          CALL PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=VLOP0+VLOP1
            CALL TRANS_IJKL_INTPOS(LRA,LRJ,LRD,LRI,NXO)
          CALL PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)

C         CALL PRODAB(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER)
         ENDDO
!SD1(8-7) A&r(13)B&l(31)B^l(23)
         DO LRD=LRI+1,LRJ-1
           LMD=LSM_INN(LRD)
           IF(LMD.NE.JMR) CYCLE
           IWDR=JUD(LRD)
           IWDL=JUST(LRJ,LRI)
           W0SD7=W0_SD1(7)
           W1SD7=W1_SD1(7)
           NI=MOD(LRD-LRI+NORB_DZ-LRJ,2)
           IF(NI.EQ.0)   W0SD7=-W0SD7
           IF(NI.EQ.0)   W1SD7=-W1SD7
           VLOP0=W0*W0SD7
           VLOP1=W1*W1SD7
C           LIST=LIST4(LRI,LRD,LRJ,LRA)
C          WL=(VLOP0-VLOP1)*vint_ci(LIST+2)-
C     :             2*VLOP0*vint_ci(LIST+1)      !1.2
            WL=VLOP0-VLOP1
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRD,LRJ,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=-2*VLOP0
            CALL TRANS_IJKL_INTPOS(LRA,LRJ,LRD,LRI,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

C         CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
!SD1(8-8) A&r(23)B&l(31)B^l(13)
           IWDL=JUST(LRI,LRJ)
           W0SD8=W0_SD1(8)
           W1SD8=W1_SD1(8)
           NI=MOD(LRD-LRI+NORB_DZ-LRJ,2)
           IF(NI.EQ.0)   W0SD8=-W0SD8
           IF(NI.EQ.0)   W1SD8=-W1SD8
           VLOP0=W0*W0SD8
           VLOP1=W1*W1SD8
C           LIST=LIST4(LRI,LRD,LRJ,LRA)
C          WL=(VLOP0-VLOP1)*vint_ci(LIST+2)-
C     :             2*VLOP0*vint_ci(LIST+1)      !1.2
            WL=VLOP0-VLOP1
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRD,LRJ,NXO)
          CALL PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=-2*VLOP0
            CALL TRANS_IJKL_INTPOS(LRA,LRJ,LRD,LRI,NXO)
          CALL PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)

C         CALL PRODAB(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER)
         ENDDO
       ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE TD_HEAD_DBL_TAIL_ACT_G(LRA,LPCOE)
!**********************************************
!     LRA, ....partial loop.......
!**********************************************
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      DIMENSION LPCOE(NORB_DZ+1:NORB_INN)
      COMMON/ONEPL/LINE,JPH,JPEL,JPER,LRG,LRS,JWL,JWR,W0,W1
      JMLR=MUL_TAB(JML,JMR)
!TD(13-1) (22)A&(23)
!TD(13-1) A&(23)C'(22)
!TD(13-5) (22)D&&l(33)B^l(23)
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        ISMA=LMI
        IF(LMI.NE.JMLR) CYCLE
        W0TD1=W0_TD(1)
        W0TD4=W0_TD(4)
        W0TD5=W0_TD(5)
        NI=MOD(NORB_DZ-LRI,2)
        IF(NI.EQ.1)W0TD1=-W0TD1
        IF(NI.EQ.1)W0TD4=-W0TD4
        IF(NI.EQ.1)W0TD5=-W0TD5

!TD(13-1) A&(23)C'(22)
        DO LRD=LRI+1,NORB_DZ
          LMD=LSM_INN(LRD)
          IF(LMD.NE.JMR) CYCLE
         IWDL=JUST(LRI,LRD)      !
          IWDR=JUD(LRD)
          VLOP0=-W0*W0TD1
C          LIST=LIST3(LRI,LRA,LRI)
C          WL=VOINT(LRI,LRA)+VINT_CI(LIST)             !310,act_coe,610,
              WL=VLOP0
        CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,LRI,LRA)
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRI,LRI,NXO)
        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

C          LIST=LIST3(LRI,LRA,LRD)
C         WL=WL+VINT_CI(LIST+1)                          !310 C'(22) COE
              WL=VLOP0
            CALL TRANS_IJKL_INTPOS(LRA,LRD,LRI,LRD,NXO)
        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

         DO LR=LRI+1,NORB_DZ
            IF(LR.EQ.LRD) CYCLE
C          LIST =LIST3(LRI,LRA,LR)
C            WL=WL+2*VINT_CI(LIST+1)-VINT_CI(LIST)       !310:NEOC=2,COE
            WL=2*VLOP0
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LR,LR,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=-VLOP0
            CALL TRANS_IJKL_INTPOS(LRA,LR,LRI,LR,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

          ENDDO
          DO LRK=norb_dz+1,LRA
C            LIST=LIST3(LRI,LRA,LRK)
            KCOE=LPCOE(LRK)
            CALL NEOC(KCOE,NOCC,TCOE)
C            WL=WL+NOCC*(VINT_CI(LIST+1)+TCOE*VINT_CI(LIST))
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
            WL=VLOP0*NOCC
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=VLOP0*NOCC*TCOE
            CALL TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

          ENDDO
C          WL=WL*VLOP0
!TD(13-5) D&rl(33)B^l(23)C'(22)          !CC (22)D&rl(33)B^l(23)???
          VLOP0=-W0*W0TD5
          DO LRK=1,LRI-1
C            LIST=LIST3(LRI,LRA,LRK)
C            WL=WL-VLOP0*(2*VINT_CI(LIST+1)-vint_ci(LIST))
            WL=-2*VLOP0
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=VLOP0
            CALL TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

          ENDDO
C          CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
        ENDDO
!-------------------------------------------------------------------
        DO LRD=NORB_FRZ+1,LRI-1
          LMD=LSM_INN(LRD)
          IF(LMD.NE.JMR) CYCLE
         IWDL=JUST(LRD,LRI)      !
          IWDR=JUD(LRD)
!TD(13-1) (22)A&(23)
          VLOP0=W0*W0TD1
C          LIST=LIST3(LRI,LRA,LRI)
C          WL=VLOP0*(VOINT(LRI,LRA)+VINT_CI(LIST))             !310,act_
              WL=VLOP0
        CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,LRI,LRA)
              CALL TRANS_IJKL_INTPOS(LRA,LRI,LRI,LRI,NXO)
       CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

          DO LR=LRI+1,NORB_DZ
C            LIST =LIST3(LRI,LRA,LR)
C            WL=WL+VLOP0*(2*VINT_CI(LIST+1)-VINT_CI(LIST)) !  310:NEOC=2
            WL=2*VLOP0
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LR,LR,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=-VLOP0
            CALL TRANS_IJKL_INTPOS(LRA,LR,LRI,LR,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

          ENDDO
          DO LRK=norb_dz+1,LRA
C            LIST=LIST3(LRI,LRA,LRK)
            KCOE=LPCOE(LRK)
            CALL NEOC(KCOE,NOCC,TCOE)
C            WL=WL+VLOP0*NOCC*(VINT_CI(LIST+1)+TCOE*VINT_CI(LIST))
            WL=VLOP0*NOCC
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=VLOP0*NOCC*TCOE
            CALL TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
          ENDDO
C            WL=WL*VLOP0
!TD(13-4) D&r&l(22)B^l(23)
          VLOP0=W0*W0TD4
          VLOP1=W1*W0TD4
C          LIST=LIST3(LRI,LRA,LRD)
C          WL=WL+(VLOP0-VLOP1)*vint_ci(LIST)-2*VLOP0*VINT_CI(LIST+1)

            WL=VLOP0-VLOP1
            CALL TRANS_IJKL_INTPOS(LRA,LRD,LRI,LRD,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=-2*VLOP0
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRD,LRD,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
!TD(13-5) D&rl(33)C"(22)B^l(23)
          VLOP0=W0*W0TD5
          DO LRK=1,LRI-1
            IF(LRK.EQ.LRD) CYCLE
C          LIST=LIST3(LRI,LRA,LRK)
C            WL=WL+VLOP0*(vint_ci(LIST)-2*VINT_CI(LIST+1))      !4.3
            WL=VLOP0
            CALL TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=-2*VLOP0
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
          ENDDO
C          CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
        ENDDO
      ENDDO
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        LMI=LSM_INN(LRI)
      DO LRJ=LRI+1,NORB_DZ
        LMJ=LSM_INN(LRJ)
        LMIJ=MUL_TAB(LMI,LMJ)
        IF(LMIJ.NE.JML) CYCLE
        IWDL=JUST(LRI,LRJ)     !

!TD(13-2) A&(23)B&r(23)B^r(32)
        DO LRD=LRJ+1,NORB_DZ
          LMD=LSM_INN(LRD)
          IF(LMD.NE.JMR) CYCLE
          W0TD2=W0_TD(2)
          W1TD2=W1_TD(2)
          NI=MOD(LRJ-LRI+NORB_DZ-LRD,2)
        IF(NI.EQ.0) W0TD2=-W0TD2
        IF(NI.EQ.0) W1TD2=-W1TD2

          IWDR=JUD(LRD)
          VLOP0=W0*W0TD2
          VLOP1=W1*W1TD2
C          LIST=LIST4(LRI,LRJ,LRD,LRA)
C         WL=VLOP0*(vint_ci(LIST+2)+vint_ci(LIST))  !1.3
C    :       -VLOP1*(vint_ci(LIST+2)-vint_ci(LIST))
            WL=VLOP0-VLOP1
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRJ,LRD,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=VLOP0+VLOP1
              CALL TRANS_IJKL_INTPOS(LRA,LRJ,LRD,LRI,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

C       CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
        ENDDO
!TD(13-3) A&(23)B&l(32)B^l(23)
        DO LRD=LRI+1,LRJ-1
          LMD=LSM_INN(LRD)
          IF(LMD.NE.JMR) CYCLE
          IWDR=JUD(LRD)
          W0TD3=W0_TD(3)
          W1TD3=W1_TD(3)
          NI=MOD(LRD-LRI+NORB_DZ-LRJ,2)
          IF(NI.EQ.0)   W0TD3=-W0TD3
          IF(NI.EQ.0)   W1TD3=-W1TD3
          VLOP0=W0*W0TD3                !D6-8
          VLOP1=W1*W1TD3
C          LIST=LIST4(LRI,LRD,LRJ,LRA)
C       WL=VLOP0*(vint_ci(LIST+2)-2*vint_ci(LIST+1))      !1.2
C     :       -VLOP1*vint_ci(LIST+2)
            WL=VLOP0-VLOP1
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRD,LRJ,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=-2*VLOP0
              CALL TRANS_IJKL_INTPOS(LRA,LRJ,LRD,LRI,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

C       CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
        ENDDO
      ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE TTDD_HEAD_DBL_TAIL_ACT_G(LRA,LPCOE)
!**********************************************
!     LRA, ....partial loop.......
!**********************************************
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      DIMENSION LPCOE(NORB_DZ+1:NORB_INN)
      COMMON/ONEPL/LINE,JPH,JPEL,JPER,LRG,LRS,JWL,JWR,W0,W1
      JMLR=MUL_TAB(JML,JMR)
!T1D1(15-1) (11)A&(13)
!T1D1(15-1) A&(13)C'(11)
!T1D1(15-5) (11)D&&l(33)B^l(13)
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        ISMA=LMI
        IF(LMI.NE.JMLR) CYCLE
        W0TD1=W0_T1D1(1)
        W0TD4=W0_T1D1(4)
        W0TD5=W0_T1D1(5)
        NI=MOD(NORB_DZ-LRI,2)
        IF(NI.EQ.1) THEN
          W0TD1=-W0TD1
          W0TD4=-W0TD4
          W0TD5=-W0TD5
        ENDIF
!T1D1(15-2) A&(13)C'(11)
        DO LRD=LRI+1,NORB_DZ
          LMD=LSM_INN(LRD)
          IF(LMD.NE.JMR) CYCLE
         IWDL=JUST(LRI,LRD)      !
          IWDR=JUD(LRD)
          VLOP0=-W0*W0TD1
C          LIST=LIST3(LRI,LRA,LRI)
C          WL=VOINT(LRI,LRA)+VINT_CI(LIST)             !310,act_coe,610,
              WL=VLOP0
        CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,LRI,LRA)
              CALL TRANS_IJKL_INTPOS(LRA,LRI,LRI,LRI,NXO)
        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

C          LIST=LIST3(LRI,LRA,LRD)
C          WL=WL+VINT_CI(LIST+1)                          !310 C'(22) CO
              WL=VLOP0
              CALL TRANS_IJKL_INTPOS(LRA,LRI,LRD,LRD,NXO)
        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

         DO LR=LRI+1,NORB_DZ
            IF(LR.EQ.LRD) CYCLE
C          LIST =LIST3(LRI,LRA,LR)
C            WL=WL+2*VINT_CI(LIST+1)-VINT_CI(LIST)       !310:NEOC=2,COE
            WL=2*VLOP0
              CALL TRANS_IJKL_INTPOS(LRA,LRI,LR,LR,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=-VLOP0
              CALL TRANS_IJKL_INTPOS(LRA,LR,LRI,LR,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

          ENDDO
          DO LRK=norb_dz+1,LRA
C            LIST=LIST3(LRI,LRA,LRK)
            KCOE=LPCOE(LRK)
            CALL NEOC(KCOE,NOCC,TCOE)
C            WL=WL+NOCC*(VINT_CI(LIST+1)+TCOE*VINT_CI(LIST))
            WL=VLOP0*NOCC
              CALL TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=VLOP0*NOCC*TCOE
              CALL TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

          ENDDO
C          WL=WL*VLOP0
!T1D1(15-5) D&rl(33)B^l(13)C'(11)
          VLOP0=-W0*W0TD5
          DO LRK=1,LRI-1
C            LIST=LIST3(LRI,LRA,LRK)
C            WL=WL-VLOP0*(2*VINT_CI(LIST+1)-vint_ci(LIST))
            WL=-2*VLOP0
              CALL TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=VLOP0
              CALL TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

          ENDDO
C          CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
        ENDDO
!-------------------------------------------------------------------
        DO LRD=NORB_FRZ+1,LRI-1
          LMD=LSM_INN(LRD)
          IF(LMD.NE.JMR) CYCLE
         IWDL=JUST(LRD,LRI)      !
          IWDR=JUD(LRD)
!T1D1(15-1) (11)A&(13)
          VLOP0=W0*W0TD1
C          LIST=LIST3(LRI,LRA,LRI)
C          WL=VLOP0*(VOINT(LRI,LRA)+VINT_CI(LIST))             !310,act_
              WL=VLOP0
        CALL PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,LRI,LRA)
              CALL TRANS_IJKL_INTPOS(LRA,LRI,LRI,LRI,NXO)
        CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)


          DO LR=LRI+1,NORB_DZ
C            LIST =LIST3(LRI,LRA,LR)
C            WL=WL+VLOP0*(2*VINT_CI(LIST+1)-VINT_CI(LIST)) !  310:NEOC=2
            WL=2*VLOP0
              CALL TRANS_IJKL_INTPOS(LRA,LRI,LR,LR,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=-VLOP0
              CALL TRANS_IJKL_INTPOS(LRA,LR,LRI,LR,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

          ENDDO
          DO LRK=norb_dz+1,LRA
C            LIST=LIST3(LRI,LRA,LRK)
            KCOE=LPCOE(LRK)
            CALL NEOC(KCOE,NOCC,TCOE)
C            WL=WL+VLOP0*NOCC*(VINT_CI(LIST+1)+TCOE*VINT_CI(LIST))
            WL=VLOP0*NOCC
              CALL TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=VLOP0*NOCC*TCOE
              CALL TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
            CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

          ENDDO
C            WL=WL*VLOP0
!T1D1(15-4) D&r&l(11)B^l(13)
          VLOP0=W0*W0TD4
          VLOP1=W1*W0TD4
C          LIST=LIST3(LRI,LRA,LRD)
C          WL=WL+(VLOP0-VLOP1)*vint_ci(LIST)-2*VLOP0*VINT_CI(LIST+1)
            WL=VLOP0-VLOP1
              CALL TRANS_IJKL_INTPOS(LRA,LRD,LRI,LRD,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=-2*VLOP0
              CALL TRANS_IJKL_INTPOS(LRA,LRI,LRD,LRD,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)


!T1D1(15-5) D&rl(33)C"(11)B^l(13)
          VLOP0=W0*W0TD5
          DO LRK=1,LRI-1
            IF(LRK.EQ.LRD) CYCLE
C          LIST=LIST3(LRI,LRA,LRK)
C           WL=WL+VLOP0*(vint_ci(LIST)-2*VINT_CI(LIST+1))      !4.3
            WL=-2*VLOP0
              CALL TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=VLOP0
              CALL TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

          ENDDO
C          CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
        ENDDO
      ENDDO
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        LMI=LSM_INN(LRI)
      DO LRJ=LRI+1,NORB_DZ
        LMJ=LSM_INN(LRJ)
        LMIJ=MUL_TAB(LMI,LMJ)
        IF(LMIJ.NE.JML) CYCLE
        IWDL=JUST(LRI,LRJ)     !

!T1D1(15-2) A&(13)B&r(13)B^r(31)
        DO LRD=LRJ+1,NORB_DZ
          LMD=LSM_INN(LRD)
          IF(LMD.NE.JMR) CYCLE
          W0TD2=W0_T1D1(2)
          W1TD2=W1_T1D1(2)
          NI=MOD(LRJ-LRI+NORB_DZ-LRD,2)
        IF(NI.EQ.0) W0TD2=-W0TD2
        IF(NI.EQ.0) W1TD2=-W1TD2

          IWDR=JUD(LRD)
          VLOP0=W0*W0TD2
          VLOP1=W1*W1TD2
C          LIST=LIST4(LRI,LRJ,LRD,LRA)
C         WL=VLOP0*(vint_ci(LIST+2)+vint_ci(LIST))  !1.3
C     :       -VLOP1*(vint_ci(LIST+2)-vint_ci(LIST))
            WL=VLOP0-VLOP1
              CALL TRANS_IJKL_INTPOS(LRA,LRI,LRJ,LRD,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=VLOP0+VLOP1
              CALL TRANS_IJKL_INTPOS(LRA,LRJ,LRD,LRI,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

C       CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
        ENDDO
!T1D1(15-3) A&(13)B&l(31)B^l(13)
        DO LRD=LRI+1,LRJ-1
          LMD=LSM_INN(LRD)
          IF(LMD.NE.JMR) CYCLE
          IWDR=JUD(LRD)
          W0TD3=W0_T1D1(3)
          W1TD3=W1_T1D1(3)
          NI=MOD(LRD-LRI+NORB_DZ-LRJ,2)
          IF(NI.EQ.0)   W0TD3=-W0TD3
          IF(NI.EQ.0)   W1TD3=-W1TD3
          VLOP0=W0*W0TD3                !D6-8
          VLOP1=W1*W1TD3
C          LIST=LIST4(LRI,LRD,LRJ,LRA)
C       WL=VLOP0*(vint_ci(LIST+2)-2*vint_ci(LIST+1))      !1.2
C     :       -VLOP1*vint_ci(LIST+2)
            WL=VLOP0-VLOP1
              CALL TRANS_IJKL_INTPOS(LRA,LRI,LRD,LRJ,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
            WL=-2*VLOP0
              CALL TRANS_IJKL_INTPOS(LRA,LRJ,LRD,LRI,NXO)
          CALL PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

C       CALL PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
        ENDDO
      ENDDO
      ENDDO

      RETURN
      END

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
      SUBROUTINE SS_DRL_ACT_BL_SGT0(LIN,LRA)
!=======================================================================
!SS(1)    ACT -BL-
!SS(1-20) Drl(33)-C"(11)-C"(22)-
!SS(1-20) (11)Drl(33)-C"(22)-
!SS(1-20) (11)(22)Drl(33)-
!=======================================================================
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      JMLR=1
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
!SS(1-16) (11)-Drl(22)-
          IF(JML.NE.MUL_TAB(LMI,LMJ).OR.JMR.NE.MUL_TAB(LMI,LMJ)) CYCLE
          IWDL=JUST(LRJ,LRI)
          IWDR=IWDL
         DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0_SS(16)
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_SS(16)
          ENDDO
          CALL Drl_BL_EXT_AR_NEW(LIN,LRJ,LRA)
!SS(1-18) Drl(11)-C"(22)-
         DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0_SS(18)
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_SS(18)
          ENDDO
          CALL Drl_BL_EXT_AR_NEW(LIN,LRI,LRA)
!SS(1-20) Drl(33)-C"(11)-C"(22)-
!SS(1-20) (11)Drl(33)-C"(22)-
!SS(1-20) (11)(22)Drl(33)-
         DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0_SS(20)
            VPLP_W1(MPL)=0.D0
          ENDDO
          IF(LRA.GT.NORB_DZ) THEN
            call Drl_BL_SUM_AR_new(LIN,LRI,LRJ,LRA)
          ELSE
            DO LRK=1,NORB_DZ
              IF(LRK.EQ.LRI) CYCLE
              IF(LRK.EQ.LRJ) CYCLE
              CALL Drl_BL_EXT_AR_NEW(LIN,LRK,LRA)
            ENDDO
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE SS_ARBL_ACT_BL_SGT0(LIN,LRA)
!=======================================================================
!SS(1)    ACT -BL-
!SS(1-1)  Ar(01)-Bl(32)-
!SS(1-3)  Ar(13)-Bl(20)-
!SS(1-6)  (11)-Ar(23)-Bl(32)-
!SS(1-7)  Ar(13)-C'(21)-Bl(32)-
!SS(1-8)  Ar(13)-C'(22)-Bl(31)-
!SS(1-9)  Ar(23)-C'(11)-Bl(32)-
!SS(1-11) Ar(13)-Bl(31)-C"(22)-
!SS(1-12) Ar(13)-Bl(32)-C"(21)-
!SS(1-13) Ar(23)-Bl(31)-C"(12)-
!=======================================================================
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      JMLR=MUL_TAB(JML,JMR)
      ISMA=MUL_TAB(IML,IMR)
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JMLR) CYCLE
          IJK=LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ)
          intpos=INTIND_IJKA(IJK)
!-------------------------------------------------------------------
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
          NI=MOD(LRJ-LRI,2)
          IF(NI.EQ.0) THEN
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
!SS(1-1)  Ar(01)-Bl(32)-
          IF(JML.EQ.1) THEN
            IWDL=JUST(LRI,LRI)
            IWDR=JUST(LRJ,LRI)
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS1
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SS1
            ENDDO
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
          ENDIF
!SS(1-3)  Ar(13)-Bl(20)-
          IF(JMR.EQ.1) THEN
            IWDL=JUST(LRJ,LRI)
            IWDR=JUST(LRJ,LRJ)
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS3
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SS3
            ENDDO
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
          CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
          ENDIF
!-------------------------------------------------------------------
!SS(1-6)  (11)-Ar(23)-Bl(32)-
          DO LRK=NORB_FRZ+1,LRI-1
            LMK=LSM_INN(LRK)
            IF(MUL_TAB(LMK,LMI).NE.JML) CYCLE
            IWDL=JUST(LRI,LRK)
            IWDR=JUST(LRJ,LRK)
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS6
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SS6
            ENDDO
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
          ENDDO
!-------------------------------------------------------------------
!SS(1-7)  Ar(13)-C'(21)-Bl(32)-
          DO LRK=LRI+1,LRJ-1
            LMK=LSM_INN(LRK)
            IF(MUL_TAB(LMK,LMI).NE.JML) CYCLE
            IWDL=JUST(LRK,LRI)
            IWDR=JUST(LRJ,LRK)
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=-VPLPNEW_W0(MPL)*W0SS7
              VPLP_W1(MPL)=-VPLPNEW_W1(MPL)*W1SS7
            ENDDO
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
!SS(1-8)  Ar(13)-C'(22)-Bl(31)-
            IWDL=JUST(LRK,LRI)
            IWDR=JUST(LRK,LRJ)
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=-VPLPNEW_W0(MPL)*W0SS8
              VPLP_W1(MPL)=-VPLPNEW_W1(MPL)*W1SS8
            ENDDO
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
!SS(1-9)  Ar(23)-C'(11)-Bl(32)-
            IWDL=JUST(LRI,LRK)
            IWDR=JUST(LRJ,LRK)
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=-VPLPNEW_W0(MPL)*W0SS9
              VPLP_W1(MPL)=-VPLPNEW_W1(MPL)*W1SS9
            ENDDO
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
          ENDDO
!-------------------------------------------------------------------
!SS(1-11) Ar(13)-Bl(31)-C"(22)-
          DO LRK=LRJ+1,NORB_DZ
            LMK=LSM_INN(LRK)
            IF(MUL_TAB(LMK,LMI).NE.JML) CYCLE
            IWDL=JUST(LRK,LRI)
            IWDR=JUST(LRK,LRJ)
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS11
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SS11
            ENDDO
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
!SS(1-12) Ar(13)-Bl(32)-C"(21)-
            IWDL=JUST(LRK,LRI)
            IWDR=JUST(LRJ,LRK)
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS12
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SS12
            ENDDO
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
!SS(1-13) Ar(23)-Bl(31)-C"(12)-
            IWDL=JUST(LRI,LRK)
            IWDR=JUST(LRK,LRJ)
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS13
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SS13
            ENDDO
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
          ENDDO
!-------------------------------------------------------------------
        ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE ST_ARBL_ACT_BL_SGT0(LIN,LRA)
!*********************************
!ST(2-3) Ar(13)-C'(22)-Bl(32)-
!ST(2-3) Ar(13)-Bl(32)-C'(22)-
!*********************************
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      JMLR=MUL_TAB(JML,JMR)
      ISMA=MUL_TAB(IML,IMR)
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
      DO LRJ=LRI+1,NORB_DZ
        LMJ=LSM_INN(LRJ)
        LMIJ=MUL_TAB(LMI,LMJ)
        IF(LMIJ.NE.JMLR) CYCLE

        W1ST3=W1_ST(3)
        NI=MOD(LRJ-LRI,2)
        IF(NI.EQ.0) W1ST3=-W1ST3

        IJK=LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ)
        intpos=INTIND_IJKA(IJK)
!-------------------------------------------------------------------
!ST(2-3) Ar(13)-C'(22)-Bl(32)-
          DO LRK=LRI+1,LRJ-1
            LMK=LSM_INN(LRK)
            IF(MUL_TAB(LMK,LMI).NE.JML) CYCLE
            IWDL=JUST(LRK,LRI)
            IWDR=JUST(LRK,LRJ)      !
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=0.D0
              VPLP_W1(MPL)=-VPLPNEW_W1(MPL)*W1ST3
            ENDDO
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
          CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
          ENDDO
!ST(2-3) Ar(13)-Bl(32)-C'(22)-
          DO LRK=LRJ+1,NORB_DZ
            LMK=LSM_INN(LRK)
            IF(MUL_TAB(LMK,LMI).NE.JML) CYCLE
            IWDL=JUST(LRK,LRI)
            IWDR=JUST(LRJ,LRK)      !
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=0.D0
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1ST3
            ENDDO
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
          CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE ST_DRL_ACT_BL_SGT0(LIN,LRA)
!ST(2-7) Drl(12)-C"(22)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      JMLR=1
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
      DO LRJ=LRI+1,NORB_DZ
        LMJ=LSM_INN(LRJ)
        LMIJ=MUL_TAB(LMI,LMJ)
        IF(LMIJ.NE.JML) CYCLE
!-------------------------------------------------------------------
!ST(2-7) Drl(12)-C"(22)-
          IWDL=JUST(LRJ,LRI)
          IWDR=JUST(LRI,LRJ)             !
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=0.D0
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_ST(7)
          ENDDO
          CALL Drl_BL_EXT_AR_NEW(LIN,LRI,LRA)
        ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE TS_ARBL_ACT_BL_SGT0(LIN,LRA)
!=======================================================================
!TS(3) A&R-B^L-  ACT -B&L ............................................
!TS(3-3) Ar(23)-Bl(31)-C"(22)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      JMLR=MUL_TAB(JML,JMR)
      ISMA=MUL_TAB(IML,IMR)
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JMLR) CYCLE
          IJK=LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ)
          intpos=INTIND_IJKA(IJK)
!-------------------------------------------------------------------
          W1TS3=W1_TS(3)
          NI=MOD(LRJ-LRI,2)
          IF(NI.EQ.0) THEN
            W1TS3=-W1TS3
          ENDIF
!TS(3-3) Ar(23)-Bl(31)-C"(22)-
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=0.D0
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1TS3
          ENDDO
          DO LRK=LRJ+1,NORB_DZ
            LMK=LSM_INN(LRK)
            IF(MUL_TAB(LMK,LMI).NE.JML) CYCLE
            IWDL=JUST(LRI,LRK)             !
            IWDR=JUST(LRK,LRJ)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
          ENDDO
        ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE STT_ARBL_ACT_BL_SGT1(LIN,LRA)
!=======================================================================
!STT(4) A&R-B^L-  ACT -B&L ............................................
!ST1(4-1) Ar(01)-Bl(31)-
!ST1(4-2) Ar(23)-Bl(31)-
!ST1(4-3) Ar(13)-C'(21)-Bl(31)-
!ST1(4-3) Ar(13)-Bl(31)-C"(21)-
!ST1(4-4) Ar(23)-C'(11)-Bl(31)-
!ST1(4-4) Ar(23)-Bl(31)-C"(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      JMLR=MUL_TAB(JML,JMR)
      ISMA=MUL_TAB(IML,IMR)
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
      DO LRJ=LRI+1,NORB_DZ
        LMJ=LSM_INN(LRJ)
        LMIJ=MUL_TAB(LMI,LMJ)
        W1ST1=W1_ST1(1)
        W1ST2=W1_ST1(2)
        W1ST3=W1_ST1(3)
        W1ST4=W1_ST1(4)
        IF(MOD(LRJ-LRI,2).EQ.0) THEN
          W1ST1=-W1ST1
          W1ST2=-W1ST2
          W1ST3=-W1ST3
          W1ST4=-W1ST4
        ENDIF
        IJK=LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ) !???
        intpos=INTIND_IJKA(IJK)                         !???
!ST1(4-1) Ar(01)-Bl(31)-
!ST1(4-2) Ar(23)-Bl(31)-
        IF(JML.EQ.1.AND.JMR.EQ.LMIJ) THEN
          IWDL=JUST(LRI,LRI)
          IWDR=JUST(LRI,LRJ)
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=0.D0
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1ST1
          ENDDO
          CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
        ENDIF
!ST1(4-2) (11)Ar(23)-Bl(31)-
        DO LRK=NORB_FRZ+1,LRI-1
          LMK=LSM_INN(LRK)
         IF(JML.EQ.MUL_TAB(LMK,LMI).AND.JMR.EQ.MUL_TAB(LMK,LMJ)) THEN
            IWDL=JUST(LRI,LRK)
            IWDR=JUST(LRK,LRJ)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=0.D0
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1ST2
            ENDDO
            CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
          ENDIF
        ENDDO
!ST1(4-3) Ar(13)-C'(21)-Bl(31)-
        DO LRK=LRI+1,LRJ-1
          LMK=LSM_INN(LRK)
         IF(JML.EQ.MUL_TAB(LMI,LMK).AND.JMR.EQ.MUL_TAB(LMK,LMJ)) THEN
            IWDL=JUST(LRK,LRI)
            IWDR=JUST(LRK,LRJ)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=0.D0
              VPLP_W1(MPL)=-VPLPNEW_W1(MPL)*W1ST3
            ENDDO
            CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
!ST1(4-4) Ar(23)-C'(11)-Bl(31)-
            IWDL=JUST(LRI,LRK)
            IWDR=JUST(LRK,LRJ)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=0.D0
              VPLP_W1(MPL)=-VPLPNEW_W1(MPL)*W1ST4
            ENDDO
            CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
          ENDIF
        ENDDO
!ST1(4-3) Ar(13)-Bl(31)-C"(21)-
        DO LRK=LRJ+1,NORB_DZ
          LMK=LSM_INN(LRK)
         IF(JML.EQ.MUL_TAB(LMI,LMK).AND.JMR.EQ.MUL_TAB(LMJ,LMK)) THEN
            IWDL=JUST(LRK,LRI)
            IWDR=JUST(LRJ,LRK)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=0.D0
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1ST3
            ENDDO
            CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
!ST1(4-4) Ar(23)-Bl(31)-C"(11)-
            IWDL=JUST(LRI,LRK)
            IWDR=JUST(LRJ,LRK)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=0.D0
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1ST4
            ENDDO
            CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
          ENDIF
        ENDDO
      ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE TTS_ARBL_ACT_BL_SGT1(LIN,LRA)
!=======================================================================
!TTS(5) A&R-B^L-  ACT -B&L ............................................
!T1S(5-1)   Ar(13)-Bl(10)-
!T1S(5-2)   (11)Ar(13)-Bl(32)-
!T1S(5-2)   Ar(13)-C'(11)-Bl(32)-
!T1S(5-3)   Ar(13)-Bl(31)-C"(12)-
!T1S(5-4)   Ar(13)-Bl(32)-C"(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      JMLR=MUL_TAB(JML,JMR)
      ISMA=MUL_TAB(IML,IMR)
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
      DO LRJ=LRI+1,NORB_DZ
        LMJ=LSM_INN(LRJ)
        LMIJ=MUL_TAB(LMI,LMJ)
        W1TS1=W1_T1S(1)
        W1TS2=W1_T1S(2)
        W1TS3=W1_T1S(3)
        W1TS4=W1_T1S(4)
        IF(MOD(LRJ-LRI,2).EQ.0) THEN
          W1TS1=-W1TS1
          W1TS2=-W1TS2
          W1TS3=-W1TS3
          W1TS4=-W1TS4
        ENDIF
        IJK=LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ) !???
        intpos=INTIND_IJKA(IJK)                         !???
!T1S(5-1)   Ar(13)-Bl(10)-
        IF(JMR.EQ.1.AND.JML.EQ.LMIJ) THEN
          IWDL=JUST(LRI,LRJ)
          IWDR=JUST(LRJ,LRJ)
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=0.D0
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1TS1
          ENDDO
          CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
        ENDIF
!T1S(5-2)   (11)Ar(13)-Bl(32)-
        DO LRK=NORB_FRZ+1,LRI-1
          LMK=LSM_INN(LRK)
         IF(JML.EQ.MUL_TAB(LMK,LMI).AND.JMR.EQ.MUL_TAB(LMK,LMJ)) THEN
            IWDL=JUST(LRK,LRI)
            IWDR=JUST(LRJ,LRK)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=0.D0
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1TS2
            ENDDO
            CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
          ENDIF
        ENDDO
!T1S(5-2)   Ar(13)-C'(11)-Bl(32)-
        DO LRK=LRI+1,LRJ-1
          LMK=LSM_INN(LRK)
          IF(JML.EQ.MUL_TAB(LMI,LMK).AND.JMR.EQ.MUL_TAB(LMK,LMJ)) THEN
            IWDL=JUST(LRI,LRK)
            IWDR=JUST(LRJ,LRK)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=0.D0
              VPLP_W1(MPL)=-VPLPNEW_W1(MPL)*W1TS2
            ENDDO
            CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
          ENDIF
        ENDDO
!T1S(5-3)   Ar(13)-Bl(31)-C"(12)-
        DO LRK=LRJ+1,NORB_DZ
          LMK=LSM_INN(LRK)
          IF(JML.EQ.MUL_TAB(LMI,LMK).AND.JMR.EQ.MUL_TAB(LMJ,LMK)) THEN
            IWDL=JUST(LRI,LRK)
            IWDR=JUST(LRK,LRJ)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=0.D0
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1TS3
            ENDDO
            CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
!T1S(5-4)   Ar(13)-Bl(32)-C"(11)-
            IWDL=JUST(LRI,LRK)
            IWDR=JUST(LRJ,LRK)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=0.D0
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1TS4
            ENDDO
            CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
          ENDIF
        ENDDO
      ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE TTS_Drl_ACT_BL_SGT1(LIN,LRA)
!T1S(5-5)   (11)Drl(12)-
!T1S(5-6)   Drl(12)-C"(12)-
!T1S(5-7)   Drl(12)-C"(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          IF(JML.NE.MUL_TAB(LMI,LMJ).OR.JMR.NE.MUL_TAB(LMI,LMJ)) CYCLE
!T1S(5-5)   (11)Drl(12)-
          IWDL=JUST(LRI,LRJ)
          IWDR=JUST(LRJ,LRI)
         DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=0.D0
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_T1S(5)
          ENDDO
          CALL Drl_BL_EXT_AR_NEW(LIN,LRJ,LRA)
!T1S(5-6)   Drl(11)-C"(12)-
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=0.D0
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_T1S(6)
          ENDDO
          CALL Drl_BL_EXT_AR_NEW(LIN,LRI,LRA)
!T1S(5-7)   Drl(12)-C"(11)-
          IWDL=JUST(LRI,LRJ)
          IWDR=IWDL
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=0.D0
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_T1S(7)
          ENDDO
          CALL Drl_BL_EXT_AR_NEW(LIN,LRI,LRA)
        ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE SDD_ABB_ACT_C_SGT0(LIN)
!SD1(8-5)    Ar(13)-Br(23)-BR(31)-
!SD1(8-6)    Ar(23)-Br(13)-BR(31)-
!SD1(8-7)    Ar(13)-Bl(31)-BL(23)-
!SD1(8-8)    Ar(23)-Bl(31)-BL(13)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      ISMA=MUL_TAB(IML,IMR)
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
      DO LRJ=LRI+1,NORB_DZ
        LMJ=LSM_INN(LRJ)
        LMIJ=MUL_TAB(LMI,LMJ)
        DO LRK=LRJ+1,NORB_DZ
          LMK=LSM_INN(LRK)
          W0SD5=W0_SD1(5)
          W1SD5=W1_SD1(5)
          W0SD6=W0_SD1(6)
          W1SD6=W1_SD1(6)
          W0SD7=W0_SD1(7)
          W1SD7=W1_SD1(7)
          W0SD8=W0_SD1(8)
          W1SD8=W1_SD1(8)
          NI=MOD(LRJ-LRI+NORB_DZ-LRK,2)
          IF(NI.EQ.0) THEN
            W0SD5=-W0SD5
            W1SD5=-W1SD5
            W0SD6=-W0SD6
            W1SD6=-W1SD6
            W0SD7=-W0SD7
            W1SD7=-W1SD7
            W0SD8=-W0SD8
            W1SD8=-W1SD8
         ENDIF
          IJK=LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRK-NORB_FRZ)
          INTPOS=INTIND_IJKA(IJK)

          IF(JML.EQ.LMIJ.AND.JMR.EQ.LMK) THEN
!SD1(8-5)    Ar(13)-Br(23)-BR(31)-
            IWDL=JUST(LRJ,LRI)
            IWDR=JUD(LRK)
           DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SD5
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SD5
            ENDDO
            CALL Ar_Br_Br_EXT_AR_NEW(LIN,INTPOS,ISMA)
!SD1(8-6)    Ar(23)-Br(13)-BR(31)-
            IWDL=JUST(LRI,LRJ)
            IWDR=JUD(LRK)
           DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SD6
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SD6
            ENDDO
            CALL Ar_Br_Br_EXT_AR_NEW(LIN,INTPOS,ISMA)
          ENDIF
!SD1(8-7)    Ar(13)-Bl(31)-BL(23)-
          IF(JML.EQ.MUL_TAB(LMI,LMK).AND.JMR.EQ.LMJ) THEN
            IWDL=JUST(LRK,LRI)
            IWDR=JUD(LRJ)
           DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SD7
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SD7
            ENDDO
            CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
!SD1(8-8)    Ar(23)-Bl(31)-BL(13)-
            IWDL=JUST(LRI,LRK)
            IWDR=JUD(LRJ)
           DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SD8
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SD8
            ENDDO
            CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
          ENDIF
        ENDDO
      ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE SDD_DrlBl_ACT_C_SGT0(LIN)
!SD1(8-9)    Drl(33)-BL(01)-
!SD1(8-10) Drl(11)-BL(23)-
!SD1(8-11) (11)Drl(33)-BL(23)-
!SD1(8-11) Drl(33)-C"(11)-BL(23)-
!SD1(8-12) Drl(33)-BL(13)-C'(21)-
!SD1(8-13) Drl(33)-BL(23)-C'(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        IF(JML.NE.1.OR.JMR.NE.LMI) CYCLE
        W0SD9=W0_SD1(9)
        W1SD9=W1_SD1(9)
        IF(MOD(NORB_DZ-LRI,2).EQ.1) THEN
          W0SD9=-W0SD9
          W1SD9=-W1SD9
        ENDIF
        DO LRK=1,LRI-1
!SD1(8-9)    Drl(33)-BL(01)-
          IWDL=JUST(LRI,LRI)
          IWDR=JUD(LRI)
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SD9
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SD9
          ENDDO
          CALL Drl_BL_EXT_AR_NEW(LIN,LRK,LRI)
        ENDDO
      ENDDO
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
      DO LRJ=LRI+1,NORB_DZ
        LMJ=LSM_INN(LRJ)
          W0SD10=W0_SD1(10)
          W1SD10=W1_SD1(10)
          W0SD11=W0_SD1(11)
          W1SD11=W1_SD1(11)
          W0SD12=W0_SD1(12)
          W1SD12=W1_SD1(12)
          W0SD13=W0_SD1(13)
          W1SD13=W1_SD1(13)
          IF(MOD(NORB_DZ-LRJ,2).EQ.1) THEN
            W0SD10=-W0SD10
            W1SD10=-W1SD10
            W0SD11=-W0SD11
            W1SD11=-W1SD11
          ENDIF
          IF(MOD(NORB_DZ-LRI,2).EQ.1) THEN
            W0SD12=-W0SD12
            W1SD12=-W1SD12
            W0SD13=-W0SD13
            W1SD13=-W1SD13
          ENDIF
!SD1(8-10) Drl(11)-BL(23)-
        IF(JML.EQ.MUL_TAB(LMI,LMJ).AND.JMR.EQ.LMI) THEN
          IWDL=JUST(LRJ,LRI)
         IWDR=JUD(LRI)
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SD10
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SD10
          ENDDO
          CALL Drl_BL_EXT_AR_NEW(LIN,LRI,LRJ)
        ENDIF
        DO LRK=1,LRI-1
          LMK=LSM_INN(LRK)
!SD1(8-12) Drl(33)-BL(13)-C'(21)-
          IF(JML.EQ.MUL_TAB(LMI,LMJ).AND.JMR.EQ.LMJ) THEN
             IWDL=JUST(LRJ,LRI)
            IWDR=JUD(LRJ)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=-VPLPNEW_W0(MPL)*W0SD12
              VPLP_W1(MPL)=-VPLPNEW_W1(MPL)*W1SD12
            ENDDO
            CALL Drl_BL_EXT_AR_NEW(LIN,LRK,LRI)
!SD1(8-13) Drl(33)-BL(23)-C'(11)-
            IWDL=JUST(LRI,LRJ)
            IWDR=JUD(LRJ)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=-VPLPNEW_W0(MPL)*W0SD13
              VPLP_W1(MPL)=-VPLPNEW_W1(MPL)*W1SD13
            ENDDO
            CALL Drl_BL_EXT_AR_NEW(LIN,LRK,LRI)
          ENDIF
         IF(JML.EQ.MUL_TAB(LMI,LMJ).AND.JMR.EQ.LMI) THEN
!SD1(8-11) Drl(33)-C"(11)-BL(23)-
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SD11
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SD11   !?????????
            ENDDO
            IWDL=JUST(LRJ,LRI)
            IWDR=JUD(LRI)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
            ENDDO
            CALL Drl_BL_EXT_AR_NEW(LIN,LRK,LRJ)
          ENDIF
        ENDDO
        DO LRK=NORB_FRZ+1,LRI-1
          LMK=LSM_INN(LRK)
!SD1(8-11) (11)Drl(33)-BL(23)-
          IF(JML.EQ.MUL_TAB(LMK,LMJ).AND.JMR.EQ.LMK) THEN
            IWDL=JUST(LRJ,LRK)
            IWDR=JUD(LRK)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SD11
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SD11
            ENDDO
            CALL Drl_BL_EXT_AR_NEW(LIN,LRI,LRJ)
          ENDIF
        ENDDO
      ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE SDD_DrrBR_ACT_C_SGT0(LIN)
!SD1(8-9)    Drr(03)-BR(31)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      DO LRI=NORB_FRZ+1,NORB_DZ
      DO LRJ=LRI+1,NORB_DZ
        LMJ=LSM_INN(LRJ)
        IF(JML.NE.1.OR.JMR.NE.LMJ) CYCLE
        W0SD=W0_SD1(9)
        IF(MOD(NORB_DZ-LRJ,2).EQ.1) W0SD=-W0SD
        IWDL=JUST(LRI,LRI)
        IWDR=JUD(LRJ)
        DO MPL=1,MHLP
           IWAL=LPNEW_LWEI(MPL)
           IWAR=LPNEW_RWEI(MPL)
           LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
           LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        ENDDO
        DO MPL=1,MTYPE
          VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SD
          VPLP_W1(MPL)=0.D0
        ENDDO
        CALL Drr_BR_EXT_AR(LIN,LRI,LRJ)
      ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE SDD_Ar_ACT_BrBR_SGT0(LIN,LRA)
!SD1(8-1)    Ar(01)-
!SD1(8-2)    Ar(23)-
!SD1(8-3)    Ar(13)-C'(21)-
!SD1(8-4)    Ar(23)-C'(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      ISMA=MUL_TAB(IML,IMR)
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        W0SD1=W0_SD1(1)
        W1SD1=W1_SD1(1)
        W0SD2=W0_SD1(2)
        W1SD2=W1_SD1(2)
        W0SD3=W0_SD1(3)
        W1SD3=W1_SD1(3)
        W0SD4=W0_SD1(4)
        W1SD4=W1_SD1(4)
        NI=MOD(NORB_DZ-LRI,2)
        IF(NI.EQ.1) THEN
          W0SD1=-W0SD1
          W1SD1=-W1SD1
          W0SD2=-W0SD2
          W1SD2=-W1SD2
          W0SD3=-W0SD3
          W1SD3=-W1SD3
          W0SD4=-W0SD4
          W1SD4=-W1SD4
        ENDIF
!SD1(8-1)    Ar(01)-
        IF(JML.EQ.1.AND.JMR.EQ.LMI) THEN
          IWDL=JUST(LRI,LRI)
          IWDR=JUD(LRI)
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SD1
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SD1
          ENDDO
          IJK=LRI-NORB_FRZ+LRA
          INTPOS=INTIND_IJKA(IJK)
          CALL Ar_Br_BR_EXT_AR_NEW(LIN,INTPOS,ISMA)     !!!!!!
        ENDIF
        DO LRJ=NORB_FRZ+1,LRI-1
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(JML.NE.LMIJ) CYCLE
          IF(JMR.EQ.LMJ) THEN
!SD1(8-2)    (11)Ar(23)-
            IWDL=JUST(LRI,LRJ)
            IWDR=JUD(LRJ)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SD2
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SD2
            ENDDO
            IJK=LRI-NORB_FRZ+LRA
            INTPOS=INTIND_IJKA(IJK)
            CALL Ar_Br_BR_EXT_AR_NEW(LIN,INTPOS,ISMA)     !!!!!!
          ENDIF
        ENDDO
!SD1(8-3)    Ar(13)-C'(21)-
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(JML.NE.LMIJ) CYCLE
          IF(JMR.EQ.LMJ) THEN
            IWDL=JUST(LRJ,LRI)
            IWDR=JUD(LRJ)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=-VPLPNEW_W0(MPL)*W0SD3
              VPLP_W1(MPL)=-VPLPNEW_W1(MPL)*W1SD3
            ENDDO
            IJK=LRI-NORB_FRZ+LRA
            INTPOS=INTIND_IJKA(IJK)
            CALL Ar_Br_BR_EXT_AR_NEW(LIN,INTPOS,ISMA)     !!!!!!
!SD1(8-4)    Ar(23)-C'(11)-
            IWDL=JUST(LRI,LRJ)
            IWDR=JUD(LRJ)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=-VPLPNEW_W0(MPL)*W0SD4
              VPLP_W1(MPL)=-VPLPNEW_W1(MPL)*W1SD4
            ENDDO
            IJK=LRI-NORB_FRZ+LRA
            INTPOS=INTIND_IJKA(IJK)
            CALL Ar_Br_BR_EXT_AR_NEW(LIN,INTPOS,ISMA)     !!!!!!
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END

!-----------------------------------
      SUBROUTINE SDD_Ar_ACT_BlBL_SGT0(LIN,LRA)
!SD1(8-1)    Ar(01)-
!SD1(8-2)    Ar(23)-
!SD1(8-3)    Ar(13)-C'(21)-
!SD1(8-4)    Ar(23)-C'(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      ISMA=MUL_TAB(IML,IMR)
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        W0SD1=W0_SD1(1)
        W1SD1=W1_SD1(1)
        W0SD2=W0_SD1(2)
        W1SD2=W1_SD1(2)
        W0SD3=W0_SD1(3)
        W1SD3=W1_SD1(3)
        W0SD4=W0_SD1(4)
        W1SD4=W1_SD1(4)
        NI=MOD(NORB_DZ-LRI,2)
        IF(NI.EQ.1) THEN
          W0SD1=-W0SD1
          W1SD1=-W1SD1
          W0SD2=-W0SD2
          W1SD2=-W1SD2
          W0SD3=-W0SD3
          W1SD3=-W1SD3
          W0SD4=-W0SD4
          W1SD4=-W1SD4
        ENDIF
!SD1(8-1)    Ar(01)-
        IF(JML.EQ.1.AND.JMR.EQ.LMI) THEN
          IWDL=JUST(LRI,LRI)
          IWDR=JUD(LRI)
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SD1
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SD1
          ENDDO
          IJK=LRI-NORB_FRZ+LRA
          INTPOS=INTIND_IJKA(IJK)
          CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
        ENDIF
        DO LRJ=NORB_FRZ+1,LRI-1
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(JML.NE.LMIJ) CYCLE
          IF(JMR.EQ.LMJ) THEN
!SD1(8-2)    (11)Ar(23)-
            IWDL=JUST(LRI,LRJ)
            IWDR=JUD(LRJ)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SD2
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SD2
            ENDDO
            IJK=LRI-NORB_FRZ+LRA
            INTPOS=INTIND_IJKA(IJK)
            CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
          ENDIF
        ENDDO
!SD1(8-3)    Ar(13)-C'(21)-
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(JML.NE.LMIJ) CYCLE
          IF(JMR.EQ.LMJ) THEN
            IWDL=JUST(LRJ,LRI)
            IWDR=JUD(LRJ)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=-VPLPNEW_W0(MPL)*W0SD3
              VPLP_W1(MPL)=-VPLPNEW_W1(MPL)*W1SD3
            ENDDO
            IJK=LRI-NORB_FRZ+LRA
            INTPOS=INTIND_IJKA(IJK)
            CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
!SD1(8-4)    Ar(23)-C'(11)-
            IWDL=JUST(LRI,LRJ)
            IWDR=JUD(LRJ)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=-VPLPNEW_W0(MPL)*W0SD4
              VPLP_W1(MPL)=-VPLPNEW_W1(MPL)*W1SD4
            ENDDO
            IJK=LRI-NORB_FRZ+LRA
            INTPOS=INTIND_IJKA(IJK)
            CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE DDS_ABB_ACT_C_SGT0(LIN)
!D1S(9-2)    Ar(13)-Bl(31)-BR(32)-
!D1S(9-3)    Ar(13)-Bl(32)-BR(31)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      ISMA=MUL_TAB(IML,IMR)
      DO LRI=NORB_FRZ+1,NORB_DZ
         LMI=LSM_INN(LRI)
      DO LRJ=LRI+1,NORB_DZ
         LMJ=LSM_INN(LRJ)
         DO LRK=LRJ+1,NORB_DZ
           LMK=LSM_INN(LRK)
          W1DS2=W1_D1S(2)
           W1DS3=W1_D1S(3)
           NI=MOD(LRJ-LRI+NORB_DZ-LRK,2)
           IF(NI.EQ.1) THEN
            W1DS2=-W1DS2
             W1DS3=-W1DS3
           ENDIF
           IJK=LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRK-NORB_FRZ)
           INTPOS=INTIND_IJKA(IJK)
!D1S(9-2)    Ar(13)-Bl(31)-BR(32)-
           IF(JML.EQ.LMI.AND.JMR.EQ.MUL_TAB(LMJ,LMK)) THEN
             IWDL=JUD(LRI)
             IWDR=JUST(LRK,LRJ)
            DO MPL=1,MHLP
               IWAL=LPNEW_LWEI(MPL)
               IWAR=LPNEW_RWEI(MPL)
               LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
               LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
             ENDDO
             DO MPL=1,MTYPE
               VPLP_W0(MPL)=0.D0
               VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1DS2
             ENDDO
             CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
!D1S(9-3)    Ar(13)-Bl(32)-BR(31)-
             IWDL=JUD(LRI)
             IWDR=JUST(LRJ,LRK)
            DO MPL=1,MHLP
               IWAL=LPNEW_LWEI(MPL)
               IWAR=LPNEW_RWEI(MPL)
               LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
               LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
             ENDDO
             DO MPL=1,MTYPE
               VPLP_W0(MPL)=0.D0
               VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1DS3
             ENDDO
             CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
          ENDIF
         ENDDO
      ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE TTTT_Drl_ACT_BL_SGT1(LIN,LRA)
!T1T1(12-2)  (11)Drl(11)-
!T1T1(12-2)  Drl(11)-C"(11)-
!T1T1(12-3)  (11)Drl(33)-
!T1T1(12-3)  (11)Drl(33)-C"(11)-
!T1T1(12-3)  Drl(33)-C"(11)-C"(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      IF(JML.NE.JMR) RETURN
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
      DO LRJ=LRI+1,NORB_DZ
        LMJ=LSM_INN(LRJ)
        LMIJ=MUL_TAB(LMI,LMJ)
        IF(JML.NE.LMIJ) CYCLE
!T1T1(12-2)  (11)Drl(11)-
!T1T1(12-2)  Drl(11)-C"(11)-
        DO MPL=1,MTYPE
          VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0_T1T1(2)
          VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_T1T1(2)
        ENDDO
        IWDL=JUST(LRI,LRJ)
        IWDR=IWDL
        DO MPL=1,MHLP
          IWAL=LPNEW_LWEI(MPL)
          IWAR=LPNEW_RWEI(MPL)
          LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
          LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
        ENDDO
        CALL Drl_BL_EXT_AR_NEW(LIN,LRI,LRA)
        CALL Drl_BL_EXT_AR_NEW(LIN,LRJ,LRA)
!T1T1(12-3)  (11)Drl(33)-
!T1T1(12-3)  (11)Drl(33)-C"(11)-
!T1T1(12-3)  Drl(33)-C"(11)-C"(11)-
        DO MPL=1,MTYPE
          VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0_T1T1(3)
          VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_T1T1(3)
        ENDDO
        IF(LRA.GT.NORB_DZ) THEN
          CALL Drl_BL_SUM_AR_new(LIN,LRI,LRJ,LRA)
        ELSE
          DO LRK=1,NORB_DZ
            IF(LRK.EQ.LRI) CYCLE
            IF(LRK.EQ.LRJ) CYCLE
            CALL Drl_BL_EXT_AR_NEW(LIN,LRK,LRA)
          ENDDO
        ENDIF
      ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE TTTT_ArBl_ACT_BL_SGT1(LIN,LRA)
!T1T1(12-1)  (11)Ar(13)-Bl(31)-
!T1T1(12-1)  Ar(13)-C'(11)-Bl(31)-
!T1T1(12-1)  Ar(13)-Bl(31)-C"(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      ISMA=MUL_TAB(IML,IMR)
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
      DO LRJ=LRI+1,NORB_DZ
        LMJ=LSM_INN(LRJ)
        W0TT1=W0_T1T1(1)
        W1TT1=W1_T1T1(1)
        IF(MOD(LRJ-LRI,2).EQ.0) THEN
          W0TT1=-W0TT1
          W1TT1=-W1TT1
        ENDIF
        DO LRK=NORB_FRZ+1,LRI-1
          LMK=LSM_INN(LRK)
!T1T1(12-1)  (11)Ar(13)-Bl(31)-
          IF(JML.EQ.MUL_TAB(LMK,LMI).AND.JMR.EQ.MUL_TAB(LMK,LMJ)) THEN
            IWDL=JUST(LRK,LRI)
            IWDR=JUST(LRK,LRJ)
            IJK=LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ)
            INTPOS=INTIND_IJKA(IJK)
            DO MPL=1,MTYPE
             VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0TT1
             VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1TT1
            ENDDO
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
            ENDDO
            CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
          ENDIF
        ENDDO
        DO LRK=LRJ+1,NORB_DZ
          LMK=LSM_INN(LRK)
          IF(JML.EQ.MUL_TAB(LMI,LMk).AND.JMR.EQ.MUL_TAB(LMJ,LMK)) THEN
!T1T1(12-1)  Ar(13)-Bl(31)-C"(11)-
            IWDL=JUST(LRI,LRK)
            IWDR=JUST(LRJ,LRK)
            IJK=LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ)
            INTPOS=INTIND_IJKA(IJK)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
             VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0TT1
             VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1TT1
            ENDDO
            CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
         ENDIF
        ENDDO
        DO LRK=LRI+1,LRJ-1
          LMK=LSM_INN(LRK)
          IF(JML.EQ.MUL_TAB(LMI,LMK).AND.JMR.EQ.MUL_TAB(LMK,LMJ)) THEN
!T1T1(12-1)  Ar(13)-C'(11)-Bl(31)-
            IJK=LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ)
            INTPOS=INTIND_IJKA(IJK)
            IWDL=JUST(LRI,LRK)
            IWDR=JUST(LRK,LRJ)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
             VPLP_W0(MPL)=-VPLPNEW_W0(MPL)*W0TT1
             VPLP_W1(MPL)=-VPLPNEW_W1(MPL)*W1TT1
            ENDDO
            CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
         ENDIF
        ENDDO
      ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE TTDD_DrlBL_ACT_C_SGT1(LIN)
!T1D1(15-4)  Drl(11)-BL(13)-
!T1D1(15-5)  (11)Drl(33)-BL(13)-
!T1D1(15-5)  Drl(33)-C"(11)-BL(13)-
!T1D1(15-5)  Drl(33)-BL(13)-C'(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      DO LRI=NORB_FRZ+1,NORB_DZ
         LMI=LSM_INN(LRI)
      DO LRJ=LRI+1,NORB_DZ
         LMJ=LSM_INN(LRJ)
         IF(JML.EQ.MUL_TAB(LMI,LMJ).AND.JMR.EQ.LMI) THEN
           IWDL=JUST(LRI,LRJ)
           IWDR=JUD(LRI)
           W0TD4=W0_T1D1(4)
           W1TD4=W1_T1D1(4)
           W0TD5=W0_T1D1(5)
           W1TD5=W1_T1D1(5)
           IF(MOD(NORB_DZ-LRJ,2).EQ.1) THEN
             W0TD4=-W0TD4
             W1TD4=-W1TD4
             W0TD5=-W0TD5
             W1TD5=-W1TD5
           ENDIF
!T1D1(15-4)  Drl(11)-BL(13)-
!T1D1(15-5)  (11)Drl(33)-BL(13)-
           DO MPL=1,MTYPE
             VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0TD4
             VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1TD4
           ENDDO
           DO MPL=1,MHLP
             IWAL=LPNEW_LWEI(MPL)
             IWAR=LPNEW_RWEI(MPL)
             LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
             LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
           ENDDO
           CALL Drl_BL_EXT_AR_NEW(LIN,LRI,LRJ)

           DO MPL=1,MTYPE
             VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0TD5
             VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1TD5
           ENDDO
           DO LRK=1,LRI-1
!T1D1(15-5)  Drl(33)-C"(11)-BL(13)-
             DO MPL=1,MHLP
               IWAL=LPNEW_LWEI(MPL)
               IWAR=LPNEW_RWEI(MPL)
               LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
               LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
             ENDDO
             CALL Drl_BL_EXT_AR_NEW(LIN,LRK,LRJ)
           ENDDO
!T1D1(15-5)  (11)Drl(33)-BL(13)-
           DO LRK=LRI+1,LRJ-1
             DO MPL=1,MHLP
               IWAL=LPNEW_LWEI(MPL)
               IWAR=LPNEW_RWEI(MPL)
               LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
               LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
             ENDDO
             CALL Drl_BL_EXT_AR_NEW(LIN,LRK,LRJ)
           ENDDO
         ENDIF
         W0TD5=W0_T1D1(5)
         W1TD5=W1_T1D1(5)
         IF(MOD(NORB_DZ-LRI,2).EQ.0) THEN
           W0TD5=-W0TD5
           W1TD5=-W1TD5
         ENDIF
         IF(JML.EQ.MUL_TAB(LMI,LMJ).AND.JMR.EQ.LMJ) THEN
           DO MPL=1,MTYPE
             VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0TD5
             VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1TD5
           ENDDO
!T1D1(15-5)  Drl(33)-BL(13)-C'(11)-
           DO LRK=1,LRI-1
             IWDL=JUST(LRI,LRJ)
             IWDR=JUD(LRJ)
             DO MPL=1,MHLP
               IWAL=LPNEW_LWEI(MPL)
               IWAR=LPNEW_RWEI(MPL)
               LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
               LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
             ENDDO
             CALL Drl_BL_EXT_AR_NEW(LIN,LRK,LRI)
           ENDDO
         ENDIF
      ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE TTDD_ABB_ACT_C_SGT1(LIN)
!T1D1(15-2)  Ar(13)-Br(13)-BR(31)-
!T1D1(15-3)  Ar(13)-Bl(31)-Bl(13)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      ISMA=MUL_TAB(IML,IMR)
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
      DO LRJ=LRI+1,NORB_DZ
        LMJ=LSM_INN(LRJ)
        DO LRK=LRJ+1,NORB_DZ
           LMK=LSM_INN(LRK)
           W1TD2=W1_T1D1(2)
           W0TD3=W0_T1D1(3)
           W1TD3=W1_T1D1(3)
           NI=MOD(LRJ-LRI+NORB_DZ-LRK,2)
           IF(NI.EQ.0) THEN
             W1TD2=-W1TD2
             W0TD3=-W0TD3
             W1TD3=-W1TD3
           ENDIF
           IJK=LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRK-NORB_FRZ)
           INTPOS=INTIND_IJKA(IJK)
           IF(JML.EQ.MUL_TAB(LMI,LMJ).AND.JMR.EQ.LMK) THEN
!T1D1(15-2)  Ar(13)-Br(13)-BR(31)-
               IWDL=JUST(LRI,LRJ)
               IWDR=JUD(LRK)
               DO MPL=1,MHLP
                IWAL=LPNEW_LWEI(MPL)
                IWAR=LPNEW_RWEI(MPL)
                LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
                LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
               ENDDO
               DO MPL=1,MTYPE
                 VPLP_W0(MPL)=0.D0
                 VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1TD2
               ENDDO
               CALL Ar_Br_BR_EXT_AR_NEW(LIN,INTPOS,ISMA)
           ENDIF
           IF(JML.EQ.MUL_TAB(LMI,LMK).AND.JMR.EQ.LMJ) THEN
!T1D1(15-3)  Ar(13)-Bl(31)-Bl(13)-
               IWDL=JUST(LRI,LRK)
               IWDR=JUD(LRJ)
               DO MPL=1,MHLP
                IWAL=LPNEW_LWEI(MPL)
                IWAR=LPNEW_RWEI(MPL)
                LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
                LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
               ENDDO
               DO MPL=1,MTYPE
                 VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0TD3
                 VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1TD3
               ENDDO
               CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
          ENDIF
        ENDDO
      ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE TTDD_Ar_ACT_BrBR_SGT0(LIN,LRA)
!T1D1(15-1)  (11)Ar(13)-
!T1D1(15-1)  Ar(13)-C'(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      ISMA=MUL_TAB(IML,IMR)
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        W0TD1=W0_T1D1(1)
        W1TD1=W1_T1D1(1)
        IF(MOD(NORB_DZ-LRI,2).EQ.1) THEN
          W0TD1=-W0TD1
          W1TD1=-W1TD1
        ENDIF
      DO LRJ=NORB_FRZ+1,LRI-1
!T1D1(15-1)  (11)Ar(13)-
        LMJ=LSM_INN(LRJ)
        LMIJ=MUL_TAB(LMI,LMJ)
        IF(JML.EQ.LMIJ.AND.JMR.EQ.LMJ) THEN
          IWDL=JUST(LRJ,LRI)
          IWDR=JUD(LRJ)
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0TD1
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1TD1
          ENDDO
          IJK=LRI-NORB_FRZ+LRA
          INTPOS=INTIND_IJKA(IJK)
          CALL Ar_Br_Br_EXT_AR_NEW(LIN,INTPOS,ISMA)     !!!!!!
        ENDIF
      ENDDO
      DO LRJ=LRI+1,NORB_DZ
        LMJ=LSM_INN(LRJ)
        LMIJ=MUL_TAB(LMI,LMJ)
!T1D1(15-1)  Ar(13)-C'(11)-
        IF(JML.EQ.LMIJ.AND.JMR.EQ.LMJ) THEN
          IWDL=JUST(LRI,LRJ)
          IWDR=JUD(LRJ)
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=-VPLPNEW_W0(MPL)*W0TD1
            VPLP_W1(MPL)=-VPLPNEW_W1(MPL)*W1TD1
          ENDDO
          IJK=LRI-NORB_FRZ+LRA
          INTPOS=INTIND_IJKA(IJK)
          CALL Ar_Br_Br_EXT_AR_NEW(LIN,INTPOS,ISMA)     !!!!!!
        ENDIF
      ENDDO
      ENDDO

      RETURN
      END

!------------------------------------------
      SUBROUTINE TTDD_Ar_ACT_BlBL_SGT0(LIN,LRA)
!T1D1(15-1)  (11)Ar(13)-
!T1D1(15-1)  Ar(13)-C'(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      ISMA=MUL_TAB(IML,IMR)
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        W0TD1=W0_T1D1(1)
        W1TD1=W1_T1D1(1)
        IF(MOD(NORB_DZ-LRI,2).EQ.1) THEN
          W0TD1=-W0TD1
          W1TD1=-W1TD1
        ENDIF
      DO LRJ=NORB_FRZ+1,LRI-1
!T1D1(15-1)  (11)Ar(13)-
        LMJ=LSM_INN(LRJ)
        LMIJ=MUL_TAB(LMI,LMJ)
        IF(JML.EQ.LMIJ.AND.JMR.EQ.LMJ) THEN
          IWDL=JUST(LRJ,LRI)
          IWDR=JUD(LRJ)
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0TD1
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1TD1
          ENDDO
          IJK=LRI-NORB_FRZ+LRA
          INTPOS=INTIND_IJKA(IJK)
          CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
        ENDIF
      ENDDO
      DO LRJ=LRI+1,NORB_DZ
!T1D1(15-1)  Ar(13)-C'(11)-
        LMJ=LSM_INN(LRJ)
        LMIJ=MUL_TAB(LMI,LMJ)
        IF(JML.EQ.LMIJ.AND.JMR.EQ.LMJ) THEN
          IWDL=JUST(LRI,LRJ)
          IWDR=JUD(LRJ)
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=-VPLPNEW_W0(MPL)*W0TD1
            VPLP_W1(MPL)=-VPLPNEW_W1(MPL)*W1TD1
          ENDDO
          IJK=LRI-NORB_FRZ+LRA
          INTPOS=INTIND_IJKA(IJK)
          CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
        ENDIF
      ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE TTV_ArBr_ACT_C_SGT1(LIN,LRA)
!T1V(18) Ar(13)-Br(13)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      isma=mul_tab(iml,imr)
      DO LRI=NORB_FRZ+1,NORB_DZ
         LMI=LSM_INN(LRI)
      DO LRJ=LRI+1,NORB_DZ
         LMJ=LSM_INN(LRJ)
         LMIJ=MUL_TAB(LMI,LMJ)
         IF(JML.NE.LMIJ.OR.JMR.NE.1) CYCLE
         W1TV1=W1_T1V
         IF(MOD(LRJ-LRI,2).EQ.0) W1TV1=-W1TV1
         IJK=LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ)
         INTPOS=INTIND_IJKA(IJK)
         IWDL=JUST(LRI,LRJ)
         IWDR=0
         DO MPL=1,MHLP
           IWAL=LPNEW_LWEI(MPL)
           IWAR=LPNEW_RWEI(MPL)
           LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
           LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
         ENDDO
         DO MPL=1,MTYPE
           VPLP_W0(MPL)=0.D0
           VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1TV1
         ENDDO
         CALL Ar_Br_BR_EXT_AR_NEW(LIN,INTPOS,ISMA)
      ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE DDDD_ArBl_ACT_BL_SGT0(LIN,LRA)
!D1D1(20-1) Ar(13)-BL(31)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      ISMA=MUL_TAB(IML,IMR)
      DO LRI=NORB_FRZ+1,NORB_DZ
         LMI=LSM_INN(LRI)
      DO LRJ=LRI+1,NORB_DZ
         LMJ=LSM_INN(LRJ)
         IF(JML.NE.LMI.OR.JMR.NE.LMJ) CYCLE
         W0DD1=W0_D1D1(1)
         W1DD1=W1_D1D1(1)
         IF(MOD(LRJ-LRI,2).EQ.0) THEN
           W0DD1=-W0DD1
           W1DD1=-W1DD1
         ENDIF
         IJK=LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ)
         INTPOS=INTIND_IJKA(IJK)
         IWDL=JUD(LRI)
         IWDR=JUD(LRJ)
         DO MPL=1,MHLP
           IWAL=LPNEW_LWEI(MPL)
           IWAR=LPNEW_RWEI(MPL)
           LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
           LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
         ENDDO
         DO MPL=1,MTYPE
           VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0DD1
           VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1DD1
         ENDDO
         CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
      ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE DDDD_Drl_ACT_BL_SGT0(LIN,LRA)
!D1D1(20-2) Drl(11)-
!D1D1(20-3) (11)Drl(33)-
!D1D1(20-3) Drl(33)-C"(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      IF(JML.NE.JMR) RETURN
      DO LRI=NORB_FRZ+1,NORB_DZ
         LMI=LSM_INN(LRI)
         IF(JML.NE.LMI) CYCLE
!D1D1(20-2) Drl(11)-
         IWDL=JUD(LRI)
         IWDR=IWDL
         DO MPL=1,MHLP
           IWAL=LPNEW_LWEI(MPL)
           IWAR=LPNEW_RWEI(MPL)
           LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
           LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
         ENDDO
         DO MPL=1,MTYPE
           VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0_D1D1(2)
           VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_D1D1(2)
         ENDDO
         CALL Drl_BL_EXT_AR_NEW(LIN,LRI,LRA)
         DO MPL=1,MTYPE
           VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0_D1D1(3)
           VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_D1D1(3)
         ENDDO
         DO LRJ=1,NORB_DZ
           IF(LRJ.EQ.LRI) CYCLE
!D1D1(20-3) Drl(33)-C"(11)-
!D1D1(20-3) (11)Drl(33)-
           CALL Drl_BL_EXT_AR_NEW(LIN,LRJ,LRA)
         ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE DD1_ArBl_ACT_BL_SGT0(LIN,LRA)
!DD1(21) Ar(23)-Bl(31)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      ISMA=MUL_TAB(IML,IMR)
      DO LRI=NORB_FRZ+1,NORB_DZ
         LMI=LSM_INN(LRI)
      DO LRJ=LRI+1,NORB_DZ
         LMJ=LSM_INN(LRJ)
         IF(JML.NE.LMI.OR.JMR.NE.LMJ) CYCLE
         W1DD1=W1_DD1
         IF(MOD(LRJ-LRI,2).EQ.0) THEN
           W1DD1=-W1DD1
         ENDIF
         IJK=LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ)
         INTPOS=INTIND_IJKA(IJK)
         IWDL=JUD(LRI)
         IWDR=JUD(LRJ)
         DO MPL=1,MHLP
           IWAL=LPNEW_LWEI(MPL)
           IWAR=LPNEW_RWEI(MPL)
           LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
           LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
         ENDDO
         DO MPL=1,MTYPE
           VPLP_W0(MPL)=0.D0
           VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1DD1
         ENDDO
         CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
      ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE D1D_ArBl_ACT_BL_SGT0(LIN,LRA)
!DD1(21) Ar(13)-Bl(32)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      ISMA=MUL_TAB(IML,IMR)
      DO LRI=NORB_FRZ+1,NORB_DZ
         LMI=LSM_INN(LRI)
      DO LRJ=LRI+1,NORB_DZ
         LMJ=LSM_INN(LRJ)
         IF(JML.NE.LMI.OR.JMR.NE.LMJ) CYCLE
         W1DD1=W1_D1D(1)
         IF(MOD(LRJ-LRI,2).EQ.0) THEN
           W1DD1=-W1DD1
         ENDIF
         IJK=LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ)
         INTPOS=INTIND_IJKA(IJK)
         IWDL=JUD(LRI)
         IWDR=JUD(LRJ)
         DO MPL=1,MHLP
           IWAL=LPNEW_LWEI(MPL)
           IWAR=LPNEW_RWEI(MPL)
           LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
           LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
         ENDDO
         DO MPL=1,MTYPE
           VPLP_W0(MPL)=0.D0
           VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1DD1
         ENDDO
         CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
      ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE D1D_Drl_ACT_BL_SGT0(LIN,LRA)
!D1D(22-2) Drl(12)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      IF(JML.NE.JMR) RETURN
      DO LRI=NORB_FRZ+1,NORB_DZ
         LMI=LSM_INN(LRI)
         IF(JML.NE.LMI) CYCLE
!D1D(22-2) Drl(12)-
         IWDL=JUD(LRI)
         IWDR=IWDL
         DO MPL=1,MHLP
           IWAL=LPNEW_LWEI(MPL)
           IWAR=LPNEW_RWEI(MPL)
           LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
           LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
         ENDDO
         DO MPL=1,MTYPE
           VPLP_W0(MPL)=0.D0
           VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_D1D(2)
         ENDDO
         CALL Drl_BL_EXT_AR_NEW(LIN,LRI,LRA)
      ENDDO

      RETURN
      END

      SUBROUTINE  D1V_Drl_BL_ACT_C_SGT0(LIN)
!D1V(24-1)  Ar(13)-
!D1V(24-2)  Drl(33)-BL(13)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      ISMA=MUL_TAB(IML,IMR)
      IF(JMR.NE.1) RETURN
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        IF(JML.NE.LMI) CYCLE
        NI=MOD(NORB_DZ-LRI,2)
!D1V(24-1)  Ar(13)-
        W0=W0_D1V(1)
        IF(NI.EQ.1)W0=-W0_D1V(1)
        IWDL=JUD(LRI)
        IWDR=0
        DO MPL=1,MTYPE
          VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0
          VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W0
        ENDDO
        DO MPL=1,MHLP
          IWAL=LPNEW_LWEI(MPL)
          IWAR=LPNEW_RWEI(MPL)
          LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        ENDDO
        CALL Ar_DV_EXT_AR(51,isma,LRI,0)   !ar_dv
!D1V(24-2)  Drl(33)-BL(13)-
        W0=W0_D1V(2)
        IF(NI.EQ.1)W0=-W0_D1V(2)
        DO MPL=1,MTYPE
          VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0
          VPLP_W1(MPL)=0.D0
        ENDDO
        DO LRJ=1,LRI-1
          CALL Drl_BL_EXT_AR_NEW(LIN,LRJ,LRI)
        ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE D1V_Ar_ACT_BrBR_SGT0(LIN,LRA)
!D1V(24-1)  Ar(13)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      ISMA=MUL_TAB(IML,IMR)
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        IF(JML.NE.LMI.OR.JMR.NE.1) CYCLE
        W0DV1=W0_D1V(1)
        W1DV1=W1_D1V(1)
        IF(MOD(NORB_DZ-LRI,2).EQ.1) THEN
          W0DV1=-W0_D1V(1)
          W1DV1=-W1_D1V(1)
        ENDIF
        IJK=LRI-NORB_FRZ+LRA
        INTPOS=INTIND_IJKA(IJK)
        IWDL=JUD(LRI)
        IWDR=0
        DO MPL=1,MHLP
          IWAL=LPNEW_LWEI(MPL)
          IWAR=LPNEW_RWEI(MPL)
          LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
          LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
        ENDDO
        DO MPL=1,MTYPE
          VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0DV1
          VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1DV1
        ENDDO
        CALL Ar_Br_Br_EXT_AR_NEW(LIN,INTPOS,ISMA)
      ENDDO

      RETURN
      END

      SUBROUTINE D1V_Ar_ACT_BlBL_SGT0(LIN,LRA)
!D1V(24-1)  Ar(13)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      ISMA=MUL_TAB(IML,IMR)
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        IF(JML.NE.LMI.OR.JMR.NE.1) CYCLE
        W0DV1=W0_D1V(1)
        W1DV1=W1_D1V(1)
        IF(MOD(NORB_DZ-LRI,2).EQ.1) THEN
          W0DV1=-W0_D1V(1)
          W1DV1=-W1_D1V(1)
        ENDIF
        IJK=LRI-NORB_FRZ+LRA
        INTPOS=INTIND_IJKA(IJK)
        IWDL=JUD(LRI)
        IWDR=0
        DO MPL=1,MHLP
          IWAL=LPNEW_LWEI(MPL)
          IWAR=LPNEW_RWEI(MPL)
          LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
          LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
        ENDDO
        DO MPL=1,MTYPE
          VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0DV1
          VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1DV1
        ENDDO
        CALL Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
      ENDDO

      RETURN
      END


!=======================================================
! END OF DV, NEXT VD
!=======================================================
      SUBROUTINE SS_DRL_ACT_BR_SGT0(LIN,LRA)
!SS(1-16) (11)-Drl(22)-
!SS(1-19) Drl(12)-C"(21)-
!SS(1-20) Drl(33)-C"(11)-C"(22)-
!SS(1-20) (11)Drl(33)-C"(22)-
!SS(1-20) (11)(22)Drl(33)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      IF(JML.NE.JMR) RETURN
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JML) CYCLE
!SS(1-16) (11)-Drl(22)-
          IWDL=JUST(LRJ,LRI)
          IWDR=IWDL
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0_SS(16)
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_SS(16)
          ENDDO
          CALL Drl_BR_EXT_AL_NEW(LIN,LRJ,LRA)
!SS(1-18) Drl(11)-C"(22)-
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0_SS(18)
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_SS(18)
          ENDDO
          CALL Drl_BR_EXT_AL_NEW(LIN,LRI,LRA)
!SS(1-20) Drl(33)-C"(11)-C"(22)-
!SS(1-20) (11)Drl(33)-C"(22)-
!SS(1-20) (11)(22)Drl(33)-
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0_SS(20)
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_SS(20)
          ENDDO
          IF(LRA.GT.NORB_DZ) THEN
            call Drl_BR_SUM_AL_new(LIN,LRI,LRJ,LRA)
          ELSE
            DO LRK=1,norb_dz
              if(LRK.EQ.LRI) CYCLE
              if(LRK.EQ.LRJ) CYCLE
              CALL Drl_BR_EXT_AL_NEW(LIN,LRK,LRA)
            ENDDO
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END

!============================================
      SUBROUTINE SS_S_DRL_ACT_BR_SGT0(LIN,LRA)
!SS(1-19) Drl(12)-C"(21)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      IF(JML.NE.JMR) RETURN
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JML) CYCLE
!SS(1-19) Drl(12)-C"(21)-
          IWDL=JUST(LRJ,LRI)
          IWDR=JUST(LRI,LRJ)
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0_SS(19)
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_SS(19)
          ENDDO
          CALL Drl_BR_EXT_AL_NEW(LIN,LRI,LRA)
        ENDDO
      ENDDO

      RETURN
      END

!==============================================
      SUBROUTINE SS_S_DRL_ACT_BL_SGT0(LIN,LRA)
!SS(1-19) Drl(12)-C"(21)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      IF(JML.NE.JMR) RETURN
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JML) CYCLE
!SS(1-19) Drl(12)-C"(21)-
          IWDL=JUST(LRJ,LRI)
          IWDR=JUST(LRI,LRJ)
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0_SS(19)
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_SS(19)
          ENDDO
          CALL Drl_BL_EXT_AR_NEW(LIN,LRI,LRA)
        ENDDO
      ENDDO

      RETURN
      END


      SUBROUTINE SS_ARBL_ACT_BR_SGT0(LIN,LRA)
!=======================================================================
!SS(1-1)  Ar(01)-Bl(32)-
!SS(1-3)  Ar(13)-Bl(20)-
!SS(1-6)  (11)-Ar(23)-Bl(32)-
!SS(1-7)  Ar(13)-C'(21)-Bl(32)-
!SS(1-8)  Ar(13)-C'(22)-Bl(31)-
!SS(1-9)  Ar(23)-C'(11)-Bl(32)-
!SS(1-11) Ar(13)-Bl(31)-C"(22)-
!SS(1-12) Ar(13)-Bl(32)-C"(21)-
!SS(1-13) Ar(23)-Bl(31)-C"(12)-
!=======================================================================
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      JMLR=MUL_TAB(JML,JMR)
      ISMA=MUL_TAB(IML,IMR)
        DO LRI=NORB_FRZ+1,NORB_DZ
          LMI=LSM_INN(LRI)
          DO LRJ=LRI+1,NORB_DZ
            LMJ=LSM_INN(LRJ)
            LMIJ=MUL_TAB(LMI,LMJ)
            IF(LMIJ.NE.JMLR) CYCLE
           IJK=LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ)
            intpos=INTIND_IJKA(IJK)
!-------------------------------------------------------------------
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
          NI=MOD(LRJ-LRI,2)
          IF(NI.EQ.0) THEN
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
!-------------------------------------------------------------------
!SS(1-1)  Ar(01)-Bl(32)-
          IF(JML.EQ.1) THEN
            IWDL=JUST(LRI,LRI)
            IWDR=JUST(LRJ,LRI)
            DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS1
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SS1
            ENDDO
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            CALL Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
          ENDIF
!SS(1-3)  Ar(13)-Bl(20)-
          IF(JMR.EQ.1) THEN
            IWDL=JUST(LRJ,LRI)
            IWDR=JUST(LRJ,LRJ)
            DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS3
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SS3
            ENDDO
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
          CALL Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
          ENDIF
!-------------------------------------------------------------------
!SS(1-6)  (11)-Ar(23)-Bl(32)-
          DO LRK=NORB_FRZ+1,LRI-1
            LMK=LSM_INN(LRK)
            IF(MUL_TAB(LMK,LMI).NE.JML) CYCLE
            IWDL=JUST(LRI,LRK)
            IWDR=JUST(LRJ,LRK)
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS6
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SS6
            ENDDO
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            CALL Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
          ENDDO
!-------------------------------------------------------------------
          DO LRK=LRI+1,LRJ-1
            LMK=LSM_INN(LRK)
            IF(MUL_TAB(LMK,LMI).NE.JML) CYCLE
!SS(1-7)  Ar(13)-C'(21)-Bl(32)-
            IWDL=JUST(LRK,LRI)
            IWDR=JUST(LRJ,LRK)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=-VPLPNEW_W0(MPL)*W0SS7
              VPLP_W1(MPL)=-VPLPNEW_W1(MPL)*W1SS7
            ENDDO
            CALL Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
!SS(1-8)  Ar(13)-C'(22)-Bl(31)-
            IWDL=JUST(LRK,LRI)
            IWDR=JUST(LRK,LRJ)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=-VPLPNEW_W0(MPL)*W0SS8
              VPLP_W1(MPL)=-VPLPNEW_W1(MPL)*W1SS8
            ENDDO
            CALL Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
!SS(1-9)  Ar(23)-C'(11)-Bl(32)-
            IWDL=JUST(LRI,LRK)
            IWDR=JUST(LRJ,LRK)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=-VPLPNEW_W0(MPL)*W0SS9
              VPLP_W1(MPL)=-VPLPNEW_W1(MPL)*W1SS9
            ENDDO
            CALL Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
          ENDDO
!-------------------------------------------------------------------
          DO LRK=LRJ+1,NORB_DZ
            LMK=LSM_INN(LRK)
            IF(MUL_TAB(LMK,LMI).NE.JML) CYCLE
!SS(1-11) Ar(13)-Bl(31)-C"(22)-
            IWDL=JUST(LRK,LRI)
            IWDR=JUST(LRK,LRJ)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS11
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SS11
            ENDDO
            CALL Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
!SS(1-12) Ar(13)-Bl(32)-C"(21)-
            IWDL=JUST(LRK,LRI)
            IWDR=JUST(LRJ,LRK)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS12
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SS12
            ENDDO
            CALL Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
!SS(1-13) Ar(23)-Bl(31)-C"(12)-
            IWDL=JUST(LRI,LRK)
            IWDR=JUST(LRK,LRJ)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS13
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SS13
            ENDDO
            CALL Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
          ENDDO
!-------------------------------------------------------------------
        ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE ST_ARBL_ACT_BR_SGT0(LIN,LRA)
!ST(2-3) Ar(13)-C'(22)-Bl(32)-
!ST(2-3) Ar(13)-Bl(32)-C'(22)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      JMLR=MUL_TAB(JML,JMR)
      ISMA=MUL_TAB(IML,IMR)
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
      DO LRJ=LRI+1,NORB_DZ
        LMJ=LSM_INN(LRJ)
        LMIJ=MUL_TAB(LMI,LMJ)
        IF(LMIJ.NE.JMLR) CYCLE

        W1ST3=W1_ST(3)
        NI=MOD(LRJ-LRI,2)
        IF(NI.EQ.0) THEN
          W1ST3=-W1ST3
        ENDIF
        IJK=LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ)
        intpos=INTIND_IJKA(IJK)
!-------------------------------------------------------------------
!ST(2-3) Ar(13)-C'(22)-Bl(32)-
          DO LRK=LRI+1,LRJ-1
            LMK=LSM_INN(LRK)
            IF(MUL_TAB(LMK,LMI).NE.JML) CYCLE
            IWDL=JUST(LRK,LRI)
            IWDR=JUST(LRK,LRJ)       !
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=0.D0
              VPLP_W1(MPL)=-VPLPNEW_W1(MPL)*W1ST3
            ENDDO
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
          CALL Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
          ENDDO
!-------------------------------------------------------------------
!ST(2-3) Ar(13)-Bl(32)-C'(22)-
          DO LRK=LRJ+1,NORB_DZ
            LMK=LSM_INN(LRK)
            IF(MUL_TAB(LMK,LMI).NE.JML) CYCLE
            IWDL=JUST(LRK,LRI)
            IWDR=JUST(LRJ,LRK)   !
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=0.D0
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1ST3
            ENDDO
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
          CALL Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE ST_DRL_ACT_BR_SGT0(LIN,LRA)
!ST(2-7) Drl(12)-C"(22)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      JMLR=1
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
      DO LRJ=LRI+1,NORB_DZ
        LMJ=LSM_INN(LRJ)
        LMIJ=MUL_TAB(LMI,LMJ)
        IF(LMIJ.NE.JML) CYCLE
!-------------------------------------------------------------------
!ST(2-7) Drl(12)-C"(22)-
          IWDL=JUST(LRJ,LRI)
          IWDR=JUST(LRI,LRJ)       !
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=0.D0
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_ST(7)
          ENDDO
          CALL Drl_BR_EXT_AL_NEW(LIN,LRI,LRA)
        ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE TS_ARBL_ACT_BR_SGT0(LIN,LRA)
!=======================================================================
!TS(3) A&R-B^L-  ACT -B&R ............................................
!TS(3-3) Ar(23)-Bl(31)-C"(22)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      JMLR=MUL_TAB(JML,JMR)
      ISMA=MUL_TAB(IML,IMR)
        DO LRI=NORB_FRZ+1,NORB_DZ
          LMI=LSM_INN(LRI)
          DO LRJ=LRI+1,NORB_DZ
            LMJ=LSM_INN(LRJ)
            IJK=LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ)
            intpos=INTIND_IJKA(IJK)
!-------------------------------------------------------------------
            W1TS3=W1_TS(3)
            NI=MOD(LRJ-LRI,2)
            IF(NI.EQ.0) THEN
              W1TS3=-W1TS3
            ENDIF
!-------------------------------------------------------------------
!TS(3-3) Ar(23)-Bl(31)-C"(22)-
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=0.D0
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1TS3
            ENDDO
            DO LRK=LRJ+1,NORB_DZ
              LMK=LSM_INN(LRK)
              LMIK=MUL_TAB(LMI,LMK)
              LMJK=MUL_TAB(LMJ,LMK)
            IF(LMIK.NE.JML.OR.LMJK.NE.JMR) CYCLE
              IWDL=JUST(LRI,LRK)    !
              IWDR=JUST(LRK,LRJ)
              DO MPL=1,MHLP
                IWAL=LPNEW_LWEI(MPL)
                IWAR=LPNEW_RWEI(MPL)
                LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
                LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
              ENDDO
              CALL Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
            ENDDO
!-------------------------------------------------------------------
          ENDDO
        ENDDO
      RETURN
      END


      SUBROUTINE STT_ARBL_ACT_BR_SGT1(LIN,LRA)
!=======================================================================
!STT(4) A&R-B^L-  ACT -BR ............................................
!ST1(4-1) Ar(01)-Bl(31)-
!ST1(4-2) Ar(23)-Bl(31)-
!ST1(4-3) Ar(13)-C'(21)-Bl(31)-
!ST1(4-3) Ar(13)-Bl(31)-C"(21)-
!ST1(4-4) Ar(23)-C'(11)-Bl(31)-
!ST1(4-4) Ar(23)-Bl(31)-C"(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      JMLR=MUL_TAB(JML,JMR)
      ISMA=MUL_TAB(IML,IMR)
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
      DO LRJ=LRI+1,NORB_DZ
        LMJ=LSM_INN(LRJ)
        LMIJ=MUL_TAB(LMI,LMJ)
        W1ST1=W1_ST1(1)
        W1ST2=W1_ST1(2)
        W1ST3=W1_ST1(3)
        W1ST4=W1_ST1(4)
        IF(MOD(LRJ-LRI,2).EQ.0) THEN
          W1ST1=-W1ST1
          W1ST2=-W1ST2
          W1ST3=-W1ST3
          W1ST4=-W1ST4
        ENDIF
        IJK=LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ) !???
        intpos=INTIND_IJKA(IJK)                         !???
!ST1(4-1) Ar(01)-Bl(31)-
!ST1(4-2) Ar(23)-Bl(31)-
        IF(JML.EQ.1.AND.JMR.EQ.LMIJ) THEN
          IWDL=JUST(LRI,LRI)
          IWDR=JUST(LRI,LRJ)
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=0.D0
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1ST1
          ENDDO
          CALL Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
        ENDIF
!ST1(4-2) (11)Ar(23)-Bl(31)-
        DO LRK=NORB_FRZ+1,LRI-1
          LMK=LSM_INN(LRK)
         IF(JML.EQ.MUL_TAB(LMK,LMI).AND.JMR.EQ.MUL_TAB(LMK,LMJ)) THEN
            IWDL=JUST(LRI,LRK)
            IWDR=JUST(LRK,LRJ)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=0.D0
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1ST2
            ENDDO
            CALL Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
          ENDIF
        ENDDO
!ST1(4-3) Ar(13)-C'(21)-Bl(31)-
        DO LRK=LRI+1,LRJ-1
          LMK=LSM_INN(LRK)
         IF(JML.EQ.MUL_TAB(LMI,LMK).AND.JMR.EQ.MUL_TAB(LMK,LMJ)) THEN
            IWDL=JUST(LRK,LRI)
            IWDR=JUST(LRK,LRJ)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=0.D0
              VPLP_W1(MPL)=-VPLPNEW_W1(MPL)*W1ST3
            ENDDO
            CALL Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
!ST1(4-4) Ar(23)-C'(11)-Bl(31)-
            IWDL=JUST(LRI,LRK)
            IWDR=JUST(LRK,LRJ)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=0.D0
              VPLP_W1(MPL)=-VPLPNEW_W1(MPL)*W1ST4
            ENDDO
            CALL Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
          ENDIF
        ENDDO
!ST1(4-3) Ar(13)-Bl(31)-C"(21)-
        DO LRK=LRJ+1,NORB_DZ
          LMK=LSM_INN(LRK)
         IF(JML.EQ.MUL_TAB(LMI,LMK).AND.JMR.EQ.MUL_TAB(LMJ,LMK)) THEN
            IWDL=JUST(LRK,LRI)
            IWDR=JUST(LRJ,LRK)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=0.D0
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1ST3
            ENDDO
            CALL Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
!ST1(4-4) Ar(23)-Bl(31)-C"(11)-
            IWDL=JUST(LRI,LRK)
            IWDR=JUST(LRJ,LRK)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=0.D0
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1ST4
            ENDDO
            CALL Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
          ENDIF
        ENDDO
      ENDDO
      ENDDO

      RETURN
      END


      SUBROUTINE TTS_Drl_ACT_BR_SGT1(LIN,LRA)
!T1S(5-5)   (11)Drl(12)-
!T1S(5-6)   Drl(12)-C"(12)-
!T1S(5-7)   Drl(12)-C"(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          IF(JML.NE.MUL_TAB(LMI,LMJ).OR.JMR.NE.MUL_TAB(LMI,LMJ)) CYCLE
!T1S(5-5)   (11)Drl(12)-
          IWDL=JUST(LRI,LRJ)
          IWDR=JUST(LRJ,LRI)
         DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=0.D0
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_T1S(5)
          ENDDO
          CALL Drl_BR_EXT_AL_NEW(LIN,LRJ,LRA)
!T1S(5-6)   Drl(11)-C"(12)-
          IWDL=JUST(LRI,LRJ)
          IWDR=JUST(LRJ,LRI)
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=0.D0
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_T1S(6)
          ENDDO
          CALL Drl_BR_EXT_AL_NEW(LIN,LRI,LRA)
!T1S(5-7)   Drl(12)-C"(11)-
          IWDL=JUST(LRI,LRJ)
          IWDR=IWDL
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,IPAE,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=0.D0
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_T1S(7)
          ENDDO
          CALL Drl_BR_EXT_AL_NEW(LIN,LRI,LRA)
        ENDDO
      ENDDO

      RETURN
      END

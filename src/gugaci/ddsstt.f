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
C calculate dd, ss and tt space
      subroutine dd_drt_ci_new()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      common/lpdisk/idisk_lp,idisk_array(13)
      DATA VSQ2/0.7071067811865950D0/
      w0gdd=vsq2
!      w1gdd=-sq3vsq2
      w1gdd=-sqrt(3.d0)/sqrt(2.d0)
      idisk_lp=idisk_array(3)

      DO lpblock=1,lpblock_dd
        call read_lp()
        IpaeL=iml+1
        Ipae =imr+1
        int_dd_drl=int_dd_offset(iml,imr)
        call logicg_dd(iml,imr)
        call get_jpty(jpadlr,jptyl,jptyr)
        call get_jp(jptyl,jml,jpadl,1)
        call get_jp(jptyr,jmr,jpad,1)
C        JMLR=MUL_TAB(JML,JMR)
        if(linelp.le.12)   then
          call dd_ext_head_in_act()
        else
          call dd_ext_head_in_dbl()
        endif
      enddo
      return
      end

      subroutine dd_ext_head_in_dbl()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      LOGIC_DH=.TRUE.
      JMLR=MUL_TAB(JML,JMR)
      LPOK=JPADLR
      GOTO(101,102,103,104,105,106,10,108,10,10,111,112,
     :     113,10,115,10,10,10,119,120,121,122,123,124,125,10),LPOK
!====================================================================
!SD(6)    ACT: -B&L-
106   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      IF(JB_SYS.GT.0) THEN
        CALL SD_Ar_ACT_Bl_DD_EXT_SGT0(LRA)
      ENDIF
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        IF(LMI.NE.JMLR) CYCLE
        W0SD1=W0_SD(1)
        W0SD2=W0_SD(2)
        W0SD3=-W0_SD(3)
        W0SD4=-W0_SD(4)
        NI=MOD(NORB_DZ-LRI,2)
        IF(NI.EQ.1)W0SD1=-W0SD1
        IF(NI.EQ.1)W0SD2=-W0SD2
        IF(NI.EQ.1)W0SD3=-W0SD3
        IF(NI.EQ.1)W0SD4=-W0SD4
        IF(JML.EQ.1.AND.LMI.EQ.JMR) THEN
!SD(6-1) A&r(02)-
          IWDL=JUST(LRI,LRI)
          IWDR=JUD(LRI)
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SD1
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W0SD1
          ENDDO
          CALL Ar_BL_DD_EXT(LRI,LRA,1)
        ENDIF
!SD(6-2) C(22)-A&r(13)-
        DO LRK=NORB_FRZ+1,LRI-1
          LMK=LSM_INN(LRK)
          LMKI=MUL_TAB(LMK,LMI)
          IF(LMKI.EQ.JML.AND.LMK.EQ.JMR) THEN
            IWDL=JUST(LRK,LRI)
            IWDR=JUD(LRK)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
             ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SD2
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W0SD2
            ENDDO
            CALL Ar_BL_DD_EXT(LRI,LRA,1)
          ENDIF
        ENDDO
!SD(6-4) A&r(23)-C'(12)-
        DO LRK=LRI+1,NORB_DZ
          LMK=LSM_INN(LRK)
          LMKI=MUL_TAB(LMK,LMI)
          IF(LMKI.NE.JML.OR.LMK.NE.JMR) CYCLE
C..........................03_01.......................
!          IF(jroute_sys.GT.1) THEN
!SD(6-3) A&r(13)-C'(22)-
!            IWDL=JUST(LRK,LRI)
!            IWDR=JUD(LRK)
!            DO MPL=1,MHLP
!              IWAL=LPNEW_LWEI(MPL)
!              IWAR=LPNEW_RWEI(MPL)
!              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
!              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
!            ENDDO
!            DO MPL=1,MTYPE
!              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SD3
!              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W0SD3
!            ENDDO
!            CALL Ar_BL_DD_EXT(LRI,LRA,1)
!          ENDIF
C...........................03_01......................
          IWDL=JUST(LRI,LRK)
          IWDR=JUD(LRK)
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SD4
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W0SD4
          ENDDO
          CALL Ar_BL_DD_EXT(LRI,LRA,1)
        ENDDO
      ENDDO
      GOTO 10
!====================================================================
!TD(13) ACT -B&L-
113   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        IF(LMI.NE.JMLR) CYCLE
        W0TD1=W0_TD(1)
        NI=MOD(NORB_DZ-LRI,2)
        IF(NI.EQ.1)W0TD1=-W0TD1

        DO LRK=NORB_FRZ+1,LRI-1
          LMK=LSM_INN(LRK)
          IF(LMK.EQ.JMR) THEN
!TD(13-1) C(22)-A&r(23)-
            IWDL=JUST(LRK,LRI)    !
            IWDR=JUD(LRK)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0TD1
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W0TD1
            ENDDO
          CALL Ar_BL_DD_EXT(LRI,LRA,1)
          ENDIF
        ENDDO
        DO LRK=LRI+1,NORB_DZ
          LMK=LSM_INN(LRK)
          IF(LMK.EQ.JMR) THEN
!TD(13-1) A&r(23)-C'(22)-
            IWDL=JUST(LRI,LRK)    !
            IWDR=JUD(LRK)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=-VPLPNEW_W0(MPL)*W0TD1
              VPLP_W1(MPL)=-VPLPNEW_W1(MPL)*W0TD1
            ENDDO
            CALL Ar_BL_DD_EXT(LRI,LRA,1)
          ENDIF
        ENDDO
      ENDDO
      GOTO 10
!=======================================================================
!DV(23) ACT -C'-..................................................
123   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        IF(LMI.NE.JMLR) CYCLE
        W0DV1=W0_DV(1)
        NI=MOD(NORB_DZ-LRI,2)
        IF(NI.EQ.1)W0DV1=-W0DV1
!DV(23-1) A&r(23)-
        IWDL=JUD(LRI)
        IWDR=0
        DO MPL=1,MHLP
          IWAL=LPNEW_LWEI(MPL)
          IWAR=LPNEW_RWEI(MPL)
          LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        ENDDO
        DO MPL=1,MTYPE
          VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0DV1
          VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W0DV1
        ENDDO
        CALL Ar_BL_DD_EXT(LRI,LRA,1)
      ENDDO
      GOTO 10
!=======================================================================
!SS(1)   ACT -C"-
!-------------------------------------------------------------------
!SS(1-9)  Ar(23)-C'(11)-Bl(32)- ACT -C"-
!SS(1-11) Ar(13)-Bl(31)-C"(22)- ACT -C"-
!SS(1-12) Ar(13)-Bl(32)-C"(21)- ACT -C"-
!SS(1-13) Ar(23)-Bl(31)-C"(12)- ACT -C"-
!SS(1-16) (11)-Drl(22)-         ACT -C"-
!SS(1-18) Drl(11)-C"(22)-       ACT -C"-
!SS(1-19) Drl(12)-C"(21)-       ACT -C"-
!SS(1-20) (11)-Drl(33)-C"(22)-  ACT -C"-
!SS(1-20) Drl(33)-C"(11)-C"(22)-ACT -C"-
!-------------------------------------------------------------------
!SS(1)   ACT 14: -C"-
101   IF(LINELP.NE.14) RETURN
      IF(NLG2.EQ.1) GOTO 1011
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        DO LRJ=LRI+1,NORB_DZ
          CALL SS2_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_DD_EXT(LRI,LRJ,1)
          CALL SS4_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_DD_EXT(LRI,LRJ,1)
          CALL SS5_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_DD_EXT(LRI,LRJ,NK)
          CALL SS10_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_DD_EXT(LRI,LRJ,NK)
          CALL SS14_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_DD_EXT(LRI,LRJ,NK)
!         IF(JROUTE_SYS.GT.1) THEN
!           CALL SS1_EXT(LRI,LRJ)
!           CALL Ar_BL_DD_EXT(LRI,LRJ,1)
!           CALL SS3_EXT(LRI,LRJ)
!           CALL Ar_BL_DD_EXT(LRI,LRJ,1)
!         ENDIF
        ENDDO
      ENDDO
      IF(JB_SYS.GT.0) THEN
        CALL SS_ArBr_ACT_C_DD_EXT_SGT0()
        CALL SS_S_Drl_ACT_C_DD_EXT_SGT0()
      ENDIF
      RETURN

1011  IF(JB_SYS.GT.0) THEN
      LRA=NLG1
      CALL SS_Drl_ACT_C_DD_EXT_SGT0()
      ENDIF
      W0SS15=W0_SS(15)
      W1SS15=W1_SS(15)
      W0SS17=W0_SS(17)
      W1SS17=W1_SS(17)
      W0SS20=W0_SS(20)
      IF(JML.EQ.1.AND.JMR.EQ.1) THEN
!SS(1-20) Drl(33)-C"(00)-       ACT -C"-
        DO LR0=NORB_FRZ+1,NORB_DZ
          IWDL=JUST(LR0,LR0)
          IWDR=IWDL
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS20
            VPLP_W1(MPL)=0.D0
          ENDDO
          DO LRK=1,NORB_DZ
            IF(LRK.EQ.LR0) CYCLE
            CALL Drl_DD_EXT(LRK)
          ENDDO
        ENDDO
      ENDIF
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JML.OR.LMIJ.NE.JMR) CYCLE
          IWDL=JUST(LRI,LRJ)
          IWDR=IWDL
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          ENDDO
!SS(1-15) (22)-Drl(11)-         ACT -C"-
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS15
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SS15
          ENDDO
          CALL Drl_DD_EXT(LRJ)
!SS(1-17) Drl(22)-C"(11)-       ACT -C"-
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS17
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SS17
          ENDDO
          CALL Drl_DD_EXT(LRI)
!SS(1-20) (22)(11)Drl(33)-      ACT -C"-
!SS(1-20) (22)Drl(33)-C"(11)-   ACT -C"-
!SS(1-20) Drl(33)-C"(22)-C"(11)-ACT -C"-
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS20
            VPLP_W1(MPL)=0.D0
          ENDDO
          DO LRK=1,NORB_DZ
            IF(LRK.EQ.LRI) CYCLE
            IF(LRK.EQ.LRJ) CYCLE
            CALL Drl_DD_EXT(LRK)
          ENDDO
        ENDDO
      ENDDO
      RETURN
!=======================================================================
!ST(2)   ACT -C"-
      IF(JB_SYS.GT.0) CALL ST_Drl_ACT_C_DD_EXT()
      IF(JB_SYS.GT.0) CALL ST_ArBl_ACT_C_DD_EXT_SGT0()

102   IF(LINELP.NE.14.OR.NLG2.EQ.1) RETURN
      IF(JB_SYS.GT.0) CALL ST_Drl_ACT_C_DD_EXT_SGT0()
      IF(JB_SYS.GT.0) CALL ST_ArBl_ACT_C_DD_EXT_SGT0()
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        DO LRJ=LRI+1,NORB_DZ
          CALL ST1_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_DD_EXT(LRI,LRJ,1)
          CALL ST2_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_DD_EXT(LRI,LRJ,NK)
          CALL ST4_EXT(LRI,LRJ,NK,1)
          IF(NK.NE.0) CALL Ar_BL_DD_EXT(LRI,LRJ,NK)
          CALL ST4_EXT(LRI,LRJ,NK,-1)
          IF(NK.NE.0) CALL Ar_BL_DD_EXT(LRI,LRJ,NK)
        ENDDO
      ENDDO

      DO LRI=NORB_FRZ+1,NORB_DZ-1
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JML.OR.LMIJ.NE.JMR) CYCLE
!ST(2-5) (22)Drl(12)-          ACT -C"-
            IWDL=JUST(LRI,LRJ)
            IWDR=JUST(LRI,LRJ)      !
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=0.D0
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_ST(5)
            ENDDO
            CALL Drl_DD_EXT(LRJ)
!ST(2-6) Drl(22)-C"(12)-       ACT -C"-
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=0.D0
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_ST(6)
            ENDDO
            CALL Drl_DD_EXT(LRI)
!----------------------------------------------------------
!ST(2-3) Ar(13)-C'(22)-Bl(32)-   ACT -C"-
!ST(2-3) Ar(13)-Bl(32)-C'(22)-   ACT -C"-
!ST(2-7) Drl(12)-C"(22)-         ACT -C"-
!----------------------------------------------------------
          ENDDO
        ENDDO
      GOTO 10
!=======================================================================
!TS(3) D&R^L-  ACT -C"-
103   IF(LINELP.NE.14.OR.NLG2.EQ.1) RETURN
      IF(JB_SYS.GT.0) THEN
        CALL TS_ArBl_ACT_C_DD_EXE_SGT0()
      ENDIF
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        DO LRJ=LRI+1,NORB_DZ
          CALL TS1_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_DD_EXT(LRI,LRJ,1)
          CALL TS2_EXT(LRI,LRJ,NK,1)
          IF(NK.NE.0) CALL Ar_BL_DD_EXT(LRI,LRJ,NK)
          CALL TS2_EXT(LRI,LRJ,NK,-1)
          IF(NK.NE.0) CALL Ar_BL_DD_EXT(LRI,LRJ,NK)
          CALL TS4_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_DD_EXT(LRI,LRJ,NK)
        ENDDO
      ENDDO
      GOTO 10
!=======================================================================
!STT(4) ArBl-  ACT -C"-
104   IF(LINELP.NE.14.OR.NLG2.EQ.1) RETURN
      LRA=NLG1
      CALL STT_ArBl_ACT_C_DD_EXT_SGT1()
      RETURN
!=======================================================================
!TTS(3) ArBl-  ACT -C"-
105   IF(LINELP.NE.14.OR.NLG2.EQ.1) RETURN
      CALL TTS_ArBl_ACT_C_DD_EXT_SGT1()
      CALL TTS_Drl_ACT_C_DD_EXT_SGT1()
      RETURN
!=======================================================================
!SDD(8) Ar- ACT -Bl-
108   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      CALL SDD_Ar_ACT_Bl_DD_EXT_SGT0(LRA)
      RETURN
!=======================================================================
!TT(11-1) (22)Ar(23)-Bl(32)-      IGF=1   ACT -C"-
!TT(11-1) Ar(23)-Bl(32)-C"(22)-  IGF=1    ACT -C"-
!TT(11-1) Ar(23)-C'(22)-Bl(32)-  IGF=-1    ACT -C"-
111   IF(LINELP.NE.14) RETURN
      IF(NLG2.EQ.1) GOTO 1111
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        DO LRJ=LRI+1,NORB_DZ
          call TT1_EXT(LRI,LRJ,NK,1)
          CALL Ar_BL_DD_EXT(LRI,LRJ,NK)
          call TT1_EXT(LRI,LRJ,NK,-1)
          CALL Ar_BL_DD_EXT(LRI,LRJ,NK)
        ENDDO
      ENDDO
      RETURN

1111  W0TT2=W0_TT(2)
      W1TT2=W1_TT(2)
      W0TT3=W0_TT(3)
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JML.OR.LMIJ.NE.JMR) CYCLE
!TT(11-2) (22)Drl(22)-
!TT(11-2) Drl(22)-C"(22)-
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0TT2
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1TT2
            ENDDO
            IWDL=JUST(LRI,LRJ)   !
            IWDR=IWDL
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            CALL Drl_DD_EXT(LRI)
            CALL Drl_DD_EXT(LRJ)
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0TT3
              VPLP_W1(MPL)=0.D0
            ENDDO
            DO LRK=1,NORB_DZ
              IF(LRK.EQ.LRI) CYCLE
              IF(LRK.EQ.LRJ) CYCLE
              LMK=LSM_INN(LRK)
              LMKI=MUL_TAB(LMK,LMI)
!TT(11-3) Drl(33)-C"(22)-C"(22)-
!TT(11-3) (22)Drl(33)-C"(22)-
!TT(11-3) (22)(22)Drl(33)-
              CALL Drl_DD_EXT(LRK)
            ENDDO
          ENDDO
        ENDDO
      RETURN
!=======================================================================
!TTTT(12) Drl- ArBl-  ACT -C"-
112   IF(LINELP.NE.14) RETURN
      IF(NLG2.EQ.1) THEN
        CALL TTTT_Drl_ACT_C_DD_EXT_SGT1()
      ELSE
        CALL TTTT_ArBl_ACT_C_DD_EXT_SGT1()
      ENDIF
!=======================================================================
!T1D1(15) Ar- ACT -Bl-
115   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      CALL TTDD_Ar_ACT_Bl_DD_EXT_SGT1(LRA)
      RETURN
!=======================================================================
!DD(19) ACT -C"- ....................................................
119   IF(LINELP.NE.14) RETURN
      IF(NLG2.EQ.1) GOTO 1191
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JMLR) CYCLE
          W0DD1=W0_DD(1)
          W1DD1=W1_DD(1)
          NI=MOD(LRJ-LRI,2)
          IF(NI.EQ.0) THEN
            W0DD1=-W0DD1
            W1DD1=-W1DD1
          ENDIF
          IF(LMI.EQ.JML.AND.LMJ.EQ.JMR) THEN
!DD(19-1) Ar(23)-Bl(32)-      ACT -C"-
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0DD1
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1DD1
            ENDDO
            IWDL=JUD(LRI)
            IWDR=JUD(LRJ)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            CALL Ar_BL_DD_EXT(LRI,LRJ,1)
          ENDIF
        ENDDO
      ENDDO
      RETURN

1191  W0DD2=W0_DD(2)
      W1DD2=W1_DD(2)
      W0DD3=W0_DD(3)
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        IF(LMI.NE.JML) CYCLE
!DD(19-2) Drl(22)-
        DO MPL=1,MTYPE
          VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0DD2
          VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1DD2
        ENDDO
        IWDL=JUD(LRI)
        IWDR=IWDL
        DO MPL=1,MHLP
          IWAL=LPNEW_LWEI(MPL)
          IWAR=LPNEW_RWEI(MPL)
          LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        ENDDO
        CALL Drl_DD_EXT(LRI)
        DO MPL=1,MTYPE
          VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0DD3
          VPLP_W1(MPL)=0.D0
        ENDDO
!DD(19-3) (22)Drl(33)-
!DD(19-3) Drl(33)-C"(22)-
        DO LRK=1,NORB_DZ
          IF(LRK.EQ.LRI) CYCLE
          LMK=LSM_INN(LRK)
          LMKI=MUL_TAB(LMK,LMI)
          CALL Drl_DD_EXT(LRK)
        ENDDO
      ENDDO
      GOTO 10
!=======================================================================
!DDDD(19) ACT -C"- ....................................................
120   IF(LINELP.NE.14) RETURN
      IF(NLG2.EQ.1) THEN
        CALL DDDD_Drl_ACT_C_DD_EXT_SGT0()
      ELSE
        CALL DDDD_ArBl_ACT_C_DD_EXT_SGT0()
      ENDIF
      RETURN
!=======================================================================
!DD1(19) ACT -C"- ....................................................
121   IF(LINELP.NE.14.OR.NLG2.NE.2) RETURN
      CALL DD1_ArBl_ACT_C_DD_EXT_SGT0()
      RETURN
!=======================================================================
!D1D(19) ACT -C"- ....................................................
122   IF(LINELP.NE.14) RETURN
!      IF(NLG2.EQ.1) THEN
        CALL D1D_Drl_ACT_C_DD_EXT_SGT0()
!     ELSE
        CALL D1D_ArBl_ACT_C_DD_EXT_SGT0()
!     ENDIF
      RETURN
!=======================================================================
!D1V(24) Ar- ACT -Bl-..................................................
124   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      CALL D1V_Ar_ACT_Bl_DD_EXT_SGT0(LRA)
      RETURN
!=======================================================================
!VV(25) ACT -BL- ....................................................
125   IF(LINELP.NE.14.OR.NLG2.NE.1) RETURN
!VV(25) Drl(33)-
      IWDL=0
      IWDR=0
      DO MPL=1,MTYPE
        VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0_VV
        VPLP_W1(MPL)=0.D0
      ENDDO
      DO MPL=1,MHLP
        IWAL=LPNEW_LWEI(MPL)
        IWAR=LPNEW_RWEI(MPL)
        LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      ENDDO
      DO LRK=1,NORB_DZ
        CALL Drl_DD_EXT(LRK)
      ENDDO
      GOTO 10
C=======================================================================
10    return
      end

      subroutine DD_ext_head_in_act()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      LOGIC_DH=.FALSE.
      LRAI=NLG1
      LRAJ=NLG2
!line=5 A&r-B&l<-->EXT
!5    CONTINUE
      IF(LINELP.EQ.5) THEN
        CALL Ar_BL_DD_EXT(lrai,lraj,1)
      ENDIF
!line=9 D&rl<-->EXT
!9    CONTINUE
      IF(LINELP.EQ.9) THEN
        call Drl_DD_EXT(lrai)
      endif
      return
      end

      subroutine ss_drt_ci_new()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      common/lpdisk/idisk_lp,idisk_array(13)

      call external_space_plpmode_value_ss()
      idisk_lp=idisk_array(13)

      DO lpblock=1,lpblock_ss
        call read_lp()
        IpaeL=iml+17
        Ipae =imr+17
        imlr=mul_tab(iml,imr)
        nvalue_space_ss=iseg_downwei(9+imlr)
        idownwei_g131415=iseg_downwei(17+iml)
C        if(imlr.eq.1)imspace=iml
        call logicg_st(iml,imr,4,4)      ! irtype=4(S),lptype=5:ArBl-
        call get_jpty(jpadlr,jptyl,jptyr)
        call get_jp(jptyl,jml,jpadl,1)
        call get_jp(jptyr,jmr,jpad,1)
C        JMLR=MUL_TAB(JML,JMR)
        if(linelp.le.12)   then
          call ss_ext_head_in_act()
        else
          call ss_ext_head_in_dbl()
        endif
      enddo
      return
      end

      subroutine ss_ext_head_in_dbl()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      LOGIC_DH=.TRUE.
      JMLR=MUL_TAB(JML,JMR)
      LPOK=JPADLR
      GOTO(101,102,10,104,105,106,10,108,10,10,111,112,
     :     113,10,115,10,10,10,119,120,121,122,123,124,125,10),LPOK
!=======================================================================
!SD(6-1) ACT -B&L-
106   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      CALL SD_AR_ACT_BL(1,LRA)
      IF(JB_SYS.GT.0) CALL SD_AR_ACT_BL_SGT0(1,LRA)
      GOTO 10
!=======================================================================
!TD(13) ACT -BL-
113   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      CALL TD_AR_ACT_BL(1,LRA)
      RETURN
!=======================================================================
!DV(23) ACT -C'-..................................................
123   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        IF(LMI.NE.JMLR) CYCLE
        W0DV1=W0_DV(1)
        NI=MOD(NORB_DZ-LRI,2)
        IF(NI.EQ.1)W0DV1=-W0DV1
!DV(23-1) A&r(23)-
        IWDL=JUD(LRI)
        IWDR=0
        DO MPL=1,MHLP
          IWAL=LPNEW_LWEI(MPL)
          IWAR=LPNEW_RWEI(MPL)
          LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        ENDDO
        DO MPL=1,MTYPE
          VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0DV1
        ENDDO
        CALL Ar_BL_EXT_SS(LRI,LRA,1)
      ENDDO
      GOTO 10
!=======================================================================
!SS(1)   ACT -C"-
!-------------------------------------------------------------------
!SS(1-1)  Ar(01)-Bl(32)-        ACT -C"-
!SS(1-3)  Ar(13)-Bl(20)-        ACT -C"-
!SS(1-6)  (11)-Ar(23)-Bl(32)-   ACT -C"-
!SS(1-7)  Ar(13)-C'(21)-Bl(32)- ACT -C"-
!SS(1-8)  Ar(13)-C'(22)-Bl(31)- ACT -C"-
!SS(1-9)  Ar(23)-C'(11)-Bl(32)- ACT -C"-
!SS(1-11) Ar(13)-Bl(31)-C"(22)- ACT -C"-
!SS(1-12) Ar(13)-Bl(32)-C"(21)- ACT -C"-
!SS(1-13) Ar(23)-Bl(31)-C"(12)- ACT -C"-
!SS(1-16) (11)-Drl(22)-         ACT -C"-
!SS(1-18) Drl(11)-C"(22)-       ACT -C"-
!SS(1-19) Drl(12)-C"(21)-       ACT -C"-
!SS(1-20) (11)-Drl(33)-C"(22)-  ACT -C"-
!SS(1-20) Drl(33)-C"(11)-C"(22)-ACT -C"-
!-------------------------------------------------------------------
!SS(1)   ACT 14: -C"-
101   IF(LINELP.NE.14) RETURN
      IF(NLG2.EQ.1) GOTO 1011
      IF(JB_SYS.GT.0) THEN
       CALL SS_ArBl_ACT_C_EXT_AB_SGT0(1)
       CALL SS_S_Drl_ACT_C_EXT_AB_SGT0(1)
      ENDIF
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        DO LRJ=LRI+1,NORB_DZ
          CALL SS2_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_EXT_SS(LRI,LRJ,1)
          CALL SS4_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_EXT_SS(LRI,LRJ,1)
          CALL SS5_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_EXT_SS(LRI,LRJ,NK)
          CALL SS10_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_EXT_SS(LRI,LRJ,NK)
          CALL SS14_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_EXT_SS(LRI,LRJ,NK)
        ENDDO
      ENDDO
      RETURN

1011  IF(JB_SYS.GT.0) CALL SS_Drl_ACT_C_EXT_AB_SGT0(1)
      W0SS15=W0_SS(15)
      W1SS15=W1_SS(15)
      W0SS17=W0_SS(17)
      W1SS17=W1_SS(17)
      W0SS20=W0_SS(20)
      IF(JML.EQ.1.AND.JMR.EQ.1) THEN
!SS(1-20) Drl(33)-C"(00)-       ACT -C"-                     ! IPL(R)AD=
        DO LR0=NORB_FRZ+1,NORB_DZ
          IWDL=JUST(LR0,LR0)
          IWDR=IWDL
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS20
            VPLP_W1(MPL)=0.D0
          ENDDO
          CALL Drl_SS_SUM(LR0,0)
        ENDDO
      ENDIF
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JML.OR.LMIJ.NE.JMR) CYCLE
          IWDL=JUST(LRI,LRJ)
          IWDR=IWDL
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          ENDDO
!SS(1-15) (22)-Drl(11)-         ACT -C"-
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS15
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SS15
          ENDDO
          CALL Drl_SS_EXT(LRJ)
!SS(1-17) Drl(22)-C"(11)-       ACT -C"-
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS17
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SS17
          ENDDO
          CALL Drl_SS_EXT(LRI)
!SS(1-20) (22)(11)Drl(33)-      ACT -C"-
!SS(1-20) (22)Drl(33)-C"(11)-   ACT -C"-
!SS(1-20) Drl(33)-C"(22)-C"(11)-ACT -C"-
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS20
            VPLP_W1(MPL)=0.D0
          ENDDO
          CALL Drl_SS_SUM(LRI,LRJ)
        ENDDO
      ENDDO
      RETURN
!ST(2)   ACT -C"-
102   IF(LINELP.NE.14.OR.NLG2.EQ.1) RETURN
      IF(JB_SYS.GT.0) CALL ST_Drl_ACT_C_EXT_AB_SGT0(1)
      IF(JB_SYS.GT.0) CALL ST_ArBl_ACT_C_EXT_AB_SGT0(1)
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JML.OR.LMIJ.NE.JMR) CYCLE
!ST(2-5) (22)Drl(12)-          ACT -C"-
            IWDL=JUST(LRI,LRJ)
            IWDR=IWDL             !
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=0.D0
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_ST(5)
            ENDDO
            CALL Drl_SS_EXT(LRJ)
!ST(2-6) Drl(22)-C"(12)-       ACT -C"-
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=0.D0
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_ST(6)
            ENDDO
            CALL Drl_SS_EXT(LRI)
!----------------------------------------------------------
!ST(2-3) Ar(13)-C'(22)-Bl(32)-   ACT -C"-
!ST(2-3) Ar(13)-Bl(32)-C'(22)-   ACT -C"-
!ST(2-7) Drl(12)-C"(22)-         ACT -C"-
!----------------------------------------------------------
          ENDDO
        ENDDO
      GOTO 10
!=======================================================================
!TS(3)   ACT -C"-  no used
!======================================================
!ST1(4) Ar-Bl- Drl- ACT -C"-
104   IF(LINELP.NE.14.or.nlg2.ne.2) RETURN
        CALL STT_ArBl_ACT_C_EXT_AB_SGT1(1)
      RETURN
!=======================================================
!T1S(5) Ar-Bl- Drl ACT -C"-
105   IF(LINELP.NE.14.OR.NLG2.NE.2) RETURN
        CALL TTS_Drl_ACT_C_EXT_AB_SGT1(1)
        CALL TTS_ArBl_ACT_C_EXT_AB_SGT1(1)
      RETURN
!===========================================================
!SD1(8) Ar ACT -Bl-
108   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      IF(JB_SYS.GT.0) CALL SDD_AR_ACT_BL_SGT0(1,LRA)
      RETURN
!===========================================================
!TTTT(12) Ar-Bl- Drl- ACT -C"-
112   IF(LINELP.NE.14) RETURN
      IF(NLG2.EQ.1) THEN
        CALL TTTT_Drl_ACT_C_EXT_AB_SGT0(1)
      ELSE
        CALL TTTT_ArBl_ACT_C_EXT_AB_SGT0(1)
      ENDIF
      RETURN
!=============================================================
!T1D1(15) Ar- ACT -Bl-
115   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      CALL TTDD_AR_ACT_BL_SGT1(1,LRA)
      RETURN
!============================================================
!D1D1(20) Drl- Ar-Bl- ACT -C"-
120   IF(LINELP.NE.14) RETURN
      IF(NLG2.EQ.1) THEN
        CALL D1D1_Drl_ACT_C_EXT_AB_SGT0(1)
      ELSE
        CALL D1D1_ArBl_ACT_C_EXT_AB_SGT0(1)
      ENDIF
      RETURN
!=================================================================
!DD1(21) Ar-Bl- ACT -C"-
121   IF(LINELP.NE.14.or.nlg2.ne.2) RETURN
      CALL DD1_ArBl_ACT_C_EXT_AB_SGT0(1)
      RETURN
!=================================================================
!D1D(22) Ar-Bl- Drl- ACT -C"-
122   IF(LINELP.NE.14.or.nlg2.ne.2) RETURN
      CALL D1D_ArBl_ACT_C_EXT_AB_SGT0(1)
      CALL D1D_Drl_ACT_C_EXT_AB_SGT0(1)
      RETURN
!=================================================================
!D1V(24) Ar- ACT -Bl-
124   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      CALL D1V_Ar_ACT_Bl_EXT_AB_SGT0(1,LRA)
      RETURN
!=======================================================================
!TT(11) AR-BL   ACT -C"-
111   IF(LINELP.NE.14) RETURN
      IF(NLG2.EQ.1) GOTO 1111
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        DO LRJ=LRI+1,NORB_DZ
          call TT1_EXT(LRI,LRJ,NK,1)
          IF(NK.NE.0) CALL AR_BL_EXT_SS(LRI,LRJ,NK)
          call TT1_EXT(LRI,LRJ,NK,-1)
          IF(NK.NE.0) CALL AR_BL_EXT_SS(LRI,LRJ,NK)
        ENDDO
      ENDDO
      RETURN

1111  W0TT2=W0_TT(2)
      W1TT2=W1_TT(2)
      W0TT3=W0_TT(3)
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JML.OR.LMIJ.NE.JMR) CYCLE
!TT(11-2) (22)Drl(22)-
!TT(11-2) Drl(22)-C"(22)-
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0TT2
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1TT2
            ENDDO
            IWDL=JUST(LRI,LRJ)   !
            IWDR=IWDL
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            CALL Drl_SS_EXT(LRI)
            CALL Drl_SS_EXT(LRJ)
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0TT3
              VPLP_W1(MPL)=0.D0
            ENDDO
!TT(11-3) Drl(33)-C"(22)-C"(22)-
!TT(11-3) (22)Drl(33)-C"(22)-
!TT(11-3) (22)(22)Drl(33)-
            CALL Drl_SS_SUM(LRI,LRJ)
          ENDDO
        ENDDO
      RETURN
!=======================================================================
!DD(19) ACT -C"- ....................................................
119   IF(LINELP.NE.14) RETURN
      IF(NLG2.EQ.1) GOTO 1191
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JMLR) CYCLE
          W0DD1=W0_DD(1)
          W1DD1=W1_DD(1)
          NI=MOD(LRJ-LRI,2)
          IF(NI.EQ.0) THEN
            W0DD1=-W0DD1
            W1DD1=-W1DD1
          ENDIF
          IF(LMI.EQ.JML.AND.LMJ.EQ.JMR) THEN
!DD(19-1) Ar(23)-Bl(32)-      ACT -C"-
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0DD1
            ENDDO
            IWDL=JUD(LRI)
            IWDR=JUD(LRJ)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            CALL Ar_BL_EXT_SS(LRI,LRJ,1)
          ENDIF
        ENDDO
      ENDDO
      RETURN

1191  W0DD2=W0_DD(2)
      W1DD2=W1_DD(2)
      W0DD3=W0_DD(3)
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        IF(LMI.NE.JML) CYCLE
!DD(19-2) Drl(22)-
        DO MPL=1,MTYPE
          VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0DD2
          VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1DD2
        ENDDO
        IWDL=JUD(LRI)
        IWDR=IWDL
        DO MPL=1,MHLP
          IWAL=LPNEW_LWEI(MPL)
          IWAR=LPNEW_RWEI(MPL)
          LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        ENDDO
        CALL Drl_SS_EXT(LRI)
        DO MPL=1,MTYPE
          VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0DD3
          VPLP_W1(MPL)=0.D0
        ENDDO
!DD(19-3) (22)Drl(33)-
!DD(19-3) Drl(33)-C"(22)-
        CALL Drl_SS_SUM(LRI,0)
      ENDDO
      GOTO 10
!=======================================================================
!VV(25) ACT -C"- .false.
125   IF(LINELP.NE.14.OR.NLG2.NE.1) RETURN
!VV(25) Drl(33)-
      IWDL=0
      IWDR=0
      DO MPL=1,MTYPE
        VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0_VV
        VPLP_W1(MPL)=0.D0
      ENDDO
      DO MPL=1,MHLP
        IWAL=LPNEW_LWEI(MPL)
        IWAR=LPNEW_RWEI(MPL)
        LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      ENDDO
        CALL Drl_SS_SUM(0,0)
C      do lrk=1,norb_dz
C        CALL Drl_ss_ext(lrk)
C     enddo
      GOTO 10
C=======================================================================
10    return
      end

      subroutine ss_ext_head_in_act()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      LOGIC_DH=.FALSE.
      LRAI=NLG1
      LRAJ=NLG2
!line=5 A&r-B&l<-->EXT
!5    CONTINUE
      IF(LINELP.EQ.5) THEN
        CALL Ar_BL_EXT_SS(lrai,lraj,1)
      ENDIF
!line=9 D&rl<-->EXT
!9    CONTINUE
      IF(LINELP.EQ.9) THEN
        call Drl_SS_EXT(lrai)
      endif
      return
      end

      subroutine st_drt_ci_new()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      common/lpdisk/idisk_lp,idisk_array(13)

      call external_space_plpmode_value_ST()
      idisk_lp=idisk_array(12)

      DO lpblock=1,lpblock_st
        call read_lp()
        IpaeL=iml+17
        Ipae =imr+9
        imlr=mul_tab(iml,imr)
        idownwei_g131415=iseg_downwei(9+iml)    !(17+iml)???
        nvalue_space_ss=iseg_downwei(9+imlr)
        call logicg_st(iml,imr,4,3)             ! irtype=4(S),3(T)
        call get_jpty(jpadlr,jptyl,jptyr)
        call get_jp(jptyl,jml,jpadl,1)
        call get_jp(jptyr,jmr,jpad,1)
C        JMLR=MUL_TAB(JML,JMR)
        if(linelp.le.12)   then
          call st_ext_head_in_act()
        else
          call st_ext_head_in_dbl()
        endif
      enddo
      return
      end

      subroutine st_ext_head_in_dbl()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      LOGIC_DH=.TRUE.
      JMLR=MUL_TAB(JML,JMR)
      LPOK=JPADLR
!      GOTO(101,102,103,10,10,106,10,10,10,10,111,10,
!     :     113,10,10,10,10,10,119,10,10,10,123,10,125,10),LPOK
      GOTO(101,102,103,104,105,106,10,108,10,10,111,112,
     :     113,10,115,10,10,10,119,120,121,122,123,124,125,10),LPOK
!=======================================================================
!SD(6-1) ACT -B&L-
106   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      CALL SD_AR_ACT_BL(2,LRA)
      IF(JB_SYS.GT.0) CALL SD_AR_ACT_BL_SGT0(2,LRA)
      GOTO 10
!=======================================================================
!TD(13) ACT -BL-
113   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      CALL TD_AR_ACT_BL(2,LRA)
      RETURN
!=======================================================================
!DV(23) ACT -C'-..................................................
123   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        IF(LMI.NE.JMLR) CYCLE
        W0DV1=W0_DV(1)
        NI=MOD(NORB_DZ-LRI,2)
        IF(NI.EQ.1)W0DV1=-W0DV1
!DV(23-1) A&r(23)-
        IWDL=JUD(LRI)
        IWDR=0
        DO MPL=1,MHLP
          IWAL=LPNEW_LWEI(MPL)
          IWAR=LPNEW_RWEI(MPL)
          LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        ENDDO
        DO MPL=1,MTYPE
          VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W0DV1
        ENDDO
        CALL Ar_BL_EXT_ST(LRI,LRA,1)
      ENDDO
      GOTO 10
!=======================================================================
!SS(1)   ACT -C"-
!-------------------------------------------------------------------
!SS(1-1)  Ar(01)-Bl(32)-        ACT -C"-
!SS(1-3)  Ar(13)-Bl(20)-        ACT -C"-
!SS(1-6)  (11)-Ar(23)-Bl(32)-   ACT -C"-
!SS(1-7)  Ar(13)-C'(21)-Bl(32)- ACT -C"-
!SS(1-8)  Ar(13)-C'(22)-Bl(31)- ACT -C"-
!SS(1-9)  Ar(23)-C'(11)-Bl(32)- ACT -C"-
!SS(1-11) Ar(13)-Bl(31)-C"(22)- ACT -C"-
!SS(1-12) Ar(13)-Bl(32)-C"(21)- ACT -C"-
!SS(1-13) Ar(23)-Bl(31)-C"(12)- ACT -C"-
!SS(1-16) (11)-Drl(22)-         ACT -C"-
!SS(1-18) Drl(11)-C"(22)-       ACT -C"-
!SS(1-19) Drl(12)-C"(21)-       ACT -C"-
!SS(1-20) (11)-Drl(33)-C"(22)-  ACT -C"-
!SS(1-20) Drl(33)-C"(11)-C"(22)-ACT -C"-
!-------------------------------------------------------------------
!SS(1)   ACT 14: -C"-
101   IF(LINELP.NE.14) RETURN
      IF(NLG2.EQ.1) GOTO 1011
      LRA=NLG1
      IF(JB_SYS.GT.0) THEN
       CALL SS_ArBl_ACT_C_EXT_AB_SGT0(2)
       CALL SS_S_Drl_ACT_C_EXT_AB_SGT0(2)
      ENDIF
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        DO LRJ=LRI+1,NORB_DZ
          CALL SS2_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_EXT_ST(LRI,LRJ,1)
          CALL SS4_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_EXT_ST(LRI,LRJ,1)
          CALL SS5_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_EXT_ST(LRI,LRJ,NK)
          CALL SS10_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_EXT_ST(LRI,LRJ,NK)
          CALL SS14_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_EXT_ST(LRI,LRJ,NK)
        ENDDO
      ENDDO
      RETURN

1011  IF(JB_SYS.GT.0) CALL SS_Drl_ACT_C_EXT_AB_SGT0(2)
      W0SS15=W0_SS(15)
      W1SS15=W1_SS(15)
      W0SS17=W0_SS(17)
      W1SS17=W1_SS(17)
      W0SS20=W0_SS(20)
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JML.OR.LMIJ.NE.JMR) CYCLE
          IWDL=JUST(LRI,LRJ)
          IWDR=IWDL
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          ENDDO
!SS(1-15) (22)-Drl(11)-         ACT -C"-
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS15
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SS15
          ENDDO
          CALL Drl_ST_EXT(LRJ)
!SS(1-17) Drl(22)-C"(11)-       ACT -C"-
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS17
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SS17
          ENDDO
          CALL Drl_ST_EXT(LRI)
!SS(1-20) (22)(11)Drl(33)-      ACT -C"-
!SS(1-20) (22)Drl(33)-C"(11)-   ACT -C"-
!SS(1-20) Drl(33)-C"(22)-C"(11)-ACT -C"-
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS20
            VPLP_W1(MPL)=0.D0
          ENDDO
         DO LRK=1,NORB_DZ
            IF(LRK.EQ.LRI) CYCLE
            IF(LRK.EQ.LRJ) CYCLE
            CALL Drl_ST_EXT(LRK)
          ENDDO
        ENDDO
      ENDDO
      RETURN
!ST(2)   ACT -C"-
102   IF(LINELP.NE.14.OR.NLG2.EQ.1) RETURN
      IF(JB_SYS.GT.0) CALL ST_Drl_ACT_C_EXT_AB_SGT0(2)
      IF(JB_SYS.GT.0) CALL ST_ArBl_ACT_C_EXT_AB_SGT0(2)
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        DO LRJ=LRI+1,NORB_DZ
          CALL ST1_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_EXT_ST(LRI,LRJ,1)
          CALL ST2_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_EXT_ST(LRI,LRJ,NK)
          CALL ST4_EXT(LRI,LRJ,NK,1)
          IF(NK.NE.0) CALL Ar_BL_EXT_ST(LRI,LRJ,NK)
          CALL ST4_EXT(LRI,LRJ,NK,-1)
          IF(NK.NE.0) CALL Ar_BL_EXT_ST(LRI,LRJ,NK)
        ENDDO
      ENDDO

      DO LRI=NORB_FRZ+1,NORB_DZ-1
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JML.OR.LMIJ.NE.JMR) CYCLE
!ST(2-5) (22)Drl(12)-          ACT -C"-
            IWDL=JUST(LRI,LRJ)
            IWDR=IWDL                !
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=0.D0
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_ST(5)
            ENDDO
            CALL Drl_ST_EXT(LRJ)
!ST(2-6) Drl(22)-C"(12)-       ACT -C"-
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=0.D0
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_ST(6)
            ENDDO
            CALL Drl_ST_EXT(LRI)
!----------------------------------------------------------
!ST(2-3) Ar(13)-C'(22)-Bl(32)-   ACT -C"-
!ST(2-3) Ar(13)-Bl(32)-C'(22)-   ACT -C"-
!ST(2-7) Drl(12)-C"(22)-         ACT -C"-
!----------------------------------------------------------
          ENDDO
        ENDDO
      GOTO 10
!=======================================================================
!TS(3) A&R-B^L-  ACT -C"-
103   IF(LINELP.NE.14.OR.NLG2.EQ.1) RETURN
      IF(JB_SYS.GT.0) THEN
        CALL TS_ArBl_ACT_C_EXT_AB_SGT0(2)
      ENDIF
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        DO LRJ=LRI+1,NORB_DZ
          CALL TS1_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_EXT_ST(LRI,LRJ,1)
          CALL TS2_EXT(LRI,LRJ,NK,1)
          IF(NK.NE.0) CALL Ar_BL_EXT_ST(LRI,LRJ,NK)
          CALL TS2_EXT(LRI,LRJ,NK,-1)
          IF(NK.NE.0) CALL Ar_BL_EXT_ST(LRI,LRJ,NK)
          CALL TS4_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_EXT_ST(LRI,LRJ,NK)
        ENDDO
      ENDDO
      GOTO 10
!=======================================================================
!TT(11) Drl-  ACT -C"-
111   IF(LINELP.NE.14) RETURN
      IF(NLG2.EQ.1) GOTO 1111
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        DO LRJ=LRI+1,NORB_DZ
          call TT1_EXT(LRI,LRJ,NK,1)
          CALL AR_BL_EXT_TS(LRI,LRJ,NK)
          call TT1_EXT(LRI,LRJ,NK,-1)
          CALL AR_BL_EXT_ST(LRI,LRJ,NK)
        ENDDO
      ENDDO
      RETURN

1111  W0TT2=W0_TT(2)
      W1TT2=W1_TT(2)
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JML.OR.LMIJ.NE.JMR) CYCLE
!TT(11-2) (22)Drl(22)-
!TT(11-2) Drl(22)-C"(22)-
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0TT2
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1TT2
            ENDDO
            IWDL=JUST(LRI,LRJ)     !
            IWDR=IWDL
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            CALL Drl_ST_EXT(LRI)
            CALL Drl_ST_EXT(LRJ)
          ENDDO
        ENDDO
      RETURN
!=======================================================================
!DD(19) ACT -C"- ....................................................
119   IF(LINELP.NE.14) RETURN
      IF(NLG2.EQ.1) GOTO 1191
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JMLR) CYCLE
          W0DD1=W0_DD(1)
          W1DD1=W1_DD(1)
          NI=MOD(LRJ-LRI,2)
          IF(NI.EQ.0) THEN
            W0DD1=-W0DD1
            W1DD1=-W1DD1
          ENDIF
          IF(LMI.EQ.JML.AND.LMJ.EQ.JMR) THEN
!DD(19-1) Ar(23)-Bl(32)-      ACT -C"-
            DO MPL=1,MTYPE
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1DD1
            ENDDO
            IWDL=JUD(LRI)
            IWDR=JUD(LRJ)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            CALL Ar_BL_EXT_ST(LRI,LRJ,1)
          ENDIF
        ENDDO
      ENDDO
      RETURN

1191  W0DD2=W0_DD(2)
      W1DD2=W1_DD(2)
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        IF(LMI.NE.JML) CYCLE
!DD(19-2) Drl(22)-
        DO MPL=1,MTYPE
          VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0DD2
          VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1DD2
        ENDDO
        IWDL=JUD(LRI)
        IWDR=IWDL
        DO MPL=1,MHLP
          IWAL=LPNEW_LWEI(MPL)
          IWAR=LPNEW_RWEI(MPL)
          LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        ENDDO
        CALL Drl_ST_EXT(LRI)
      ENDDO
      GOTO 10
!======================================================
!ST1(4) Ar-Bl- Drl- ACT -C"-
104   IF(LINELP.NE.14.or.nlg2.ne.2) RETURN
        CALL STT_ArBl_ACT_C_EXT_AB_SGT1(2)
      RETURN
!=======================================================
!T1S(5) Ar-Bl- Drl ACT -C"-
105   IF(LINELP.NE.14.OR.NLG2.NE.2) RETURN
        CALL TTS_Drl_ACT_C_EXT_AB_SGT1(2)
        CALL TTS_ArBl_ACT_C_EXT_AB_SGT1(2)
      RETURN
!===========================================================
!SD1(8) Ar ACT -Bl-
108   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      IF(JB_SYS.GT.0) CALL SDD_AR_ACT_BL_SGT0(2,LRA)
      RETURN
!===========================================================
!TTTT(12) Ar-Bl- Drl- ACT -C"-
112   IF(LINELP.NE.14) RETURN
      IF(NLG2.EQ.1) THEN
        CALL TTTT_Drl_ACT_C_EXT_AB_SGT0(2)
      ELSE
        CALL TTTT_ArBl_ACT_C_EXT_AB_SGT0(2)
      ENDIF
      RETURN
!=============================================================
!T1D1(15) Ar- ACT -Bl-
115   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      CALL TTDD_AR_ACT_BL_SGT1(2,LRA)
      RETURN
!============================================================
!D1D1(20) Drl- Ar-Bl- ACT -C"-
120   IF(LINELP.NE.14) RETURN
      IF(NLG2.EQ.1) THEN
        CALL D1D1_Drl_ACT_C_EXT_AB_SGT0(2)
      ELSE
        CALL D1D1_ArBl_ACT_C_EXT_AB_SGT0(2)
      ENDIF
      RETURN
!=================================================================
!DD1(21) Ar-Bl- ACT -C"-
121   IF(LINELP.NE.14.or.nlg2.ne.2) RETURN
      CALL DD1_ArBl_ACT_C_EXT_AB_SGT0(2)
      RETURN
!=================================================================
!D1D(22) Ar-Bl- Drl- ACT -C"-
122   IF(LINELP.NE.14.or.nlg2.ne.2) RETURN
      CALL D1D_ArBl_ACT_C_EXT_AB_SGT0(2)
      CALL D1D_Drl_ACT_C_EXT_AB_SGT0(2)
      RETURN
!=================================================================
!D1V(24) Ar- ACT -Bl-
124   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      CALL D1V_Ar_ACT_Bl_EXT_AB_SGT0(2,LRA)
      RETURN
!=======================================================================
!VV(25) ACT -BL- ....................................................
125   RETURN
C=======================================================================
10    return
      end

      subroutine ST_ext_head_in_act()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      LOGIC_DH=.FALSE.
      LRAI=NLG1
      LRAJ=NLG2
!line=5 A&r-B&l<-->EXT
!5    CONTINUE
      IF(LINELP.EQ.5) THEN
        CALL Ar_BL_EXT_ST(lrai,lraj,1)
      ENDIF
!line=9 D&rl<-->EXT
!9    CONTINUE
      IF(LINELP.EQ.9) THEN
        call Drl_ST_EXT(lrai)
      endif
      return
      end

      subroutine ts_drt_ci_new()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      common/lpdisk/idisk_lp,idisk_array(13)

      call external_space_plpmode_value_TS()
      idisk_lp=idisk_array(9)
      DO lpblock=1,lpblock_ts
        call read_lp()
        IpaeL=iml+9
        Ipae =imr+17
        imlr=mul_tab(iml,imr)
        idownwei_g131415=iseg_downwei(9+iml)    !(17+iml)???
        nvalue_space_ss=iseg_downwei(9+imlr)
        call logicg_st(iml,imr,3,4)             ! irtype=4(S),3(T)
        call get_jpty(jpadlr,jptyl,jptyr)
        call get_jp(jptyl,jml,jpadl,1)
        call get_jp(jptyr,jmr,jpad,1)
C        JMLR=MUL_TAB(JML,JMR)
        if(linelp.le.12)   then
          call ts_ext_head_in_act()
        else
          call ts_ext_head_in_dbl()
        endif
      enddo
      return
      end

      subroutine ts_ext_head_in_dbl()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      LOGIC_DH=.TRUE.
      JMLR=MUL_TAB(JML,JMR)
      LPOK=JPADLR
      GOTO(101,102,103,104,105,106,10,108,10,10,111,112,
     :     113,10,115,10,10,10,119,120,121,122,123,124,125,10),LPOK
!=======================================================================
!SD(6-1) ACT -B&L-
106   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      CALL SD_AR_ACT_BL(3,LRA)
      IF(JB_SYS.GT.0) CALL SD_AR_ACT_BL_SGT0(3,LRA)
      GOTO 10
!=======================================================================
!TD(13) ACT -BL-
113   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      CALL TD_AR_ACT_BL(3,LRA)
      RETURN
!=======================================================================
!DV(23) ACT -C'-..................................................
123   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        IF(LMI.NE.JMLR) CYCLE
        W0DV1=W0_DV(1)
        NI=MOD(NORB_DZ-LRI,2)
        IF(NI.EQ.1)W0DV1=-W0DV1
!DV(23-1) A&r(23)-
        IWDL=JUD(LRI)
        IWDR=0
        DO MPL=1,MHLP
          IWAL=LPNEW_LWEI(MPL)
          IWAR=LPNEW_RWEI(MPL)
          LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        ENDDO
        DO MPL=1,MTYPE
          VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W0DV1
        ENDDO
        CALL Ar_BL_EXT_TS(LRI,LRA,1)
      ENDDO
      GOTO 10
!=======================================================================
!SS(1)   ACT -C"-
!-------------------------------------------------------------------
!SS(1-1)  Ar(01)-Bl(32)-        ACT -C"-
!SS(1-3)  Ar(13)-Bl(20)-        ACT -C"-
!SS(1-6)  (11)-Ar(23)-Bl(32)-   ACT -C"-
!SS(1-7)  Ar(13)-C'(21)-Bl(32)- ACT -C"-
!SS(1-8)  Ar(13)-C'(22)-Bl(31)- ACT -C"-
!SS(1-9)  Ar(23)-C'(11)-Bl(32)- ACT -C"-
!SS(1-11) Ar(13)-Bl(31)-C"(22)- ACT -C"-
!SS(1-12) Ar(13)-Bl(32)-C"(21)- ACT -C"-
!SS(1-13) Ar(23)-Bl(31)-C"(12)- ACT -C"-
!SS(1-16) (11)-Drl(22)-         ACT -C"-
!SS(1-18) Drl(11)-C"(22)-       ACT -C"-
!SS(1-19) Drl(12)-C"(21)-       ACT -C"-
!SS(1-20) (11)-Drl(33)-C"(22)-  ACT -C"-
!SS(1-20) Drl(33)-C"(11)-C"(22)-ACT -C"-
!-------------------------------------------------------------------
!SS(1)   ACT 14: -C"-
101   IF(LINELP.NE.14) RETURN
      IF(NLG2.EQ.1) GOTO 1011
      IF(JB_SYS.GT.0) THEN
       CALL SS_ArBl_ACT_C_EXT_AB_SGT0(3)
       CALL SS_S_Drl_ACT_C_EXT_AB_SGT0(3)
      ENDIF
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        DO LRJ=LRI+1,NORB_DZ
          CALL SS2_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL AR_BL_EXT_TS(LRI,LRJ,1)
          CALL SS4_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL AR_BL_EXT_TS(LRI,LRJ,1)
          CALL SS5_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL AR_BL_EXT_TS(LRI,LRJ,NK)
          CALL SS10_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL AR_BL_EXT_TS(LRI,LRJ,NK)
          CALL SS14_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL AR_BL_EXT_TS(LRI,LRJ,NK)
        ENDDO
      ENDDO
      RETURN

1011  IF(JB_SYS.GT.0) CALL SS_Drl_ACT_C_EXT_AB_SGT0(3)
      W0SS15=W0_SS(15)
      W1SS15=W1_SS(15)
      W0SS17=W0_SS(17)
      W1SS17=W1_SS(17)
      W0SS20=W0_SS(20)
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JML.OR.LMIJ.NE.JMR) CYCLE
          IWDL=JUST(LRI,LRJ)
          IWDR=IWDL
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          ENDDO
!SS(1-15) (22)-Drl(11)-         ACT -C"-
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS15
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SS15
          ENDDO
          CALL Drl_TS_EXT(LRJ)
!SS(1-17) Drl(22)-C"(11)-       ACT -C"-
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS17
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SS17
          ENDDO
          CALL Drl_TS_EXT(LRI)
!SS(1-20) (22)(11)Drl(33)-      ACT -C"-
!SS(1-20) (22)Drl(33)-C"(11)-   ACT -C"-
!SS(1-20) Drl(33)-C"(22)-C"(11)-ACT -C"-
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS20
            VPLP_W1(MPL)=0.D0
          ENDDO
         DO LRK=1,NORB_DZ
            IF(LRK.EQ.LRI) CYCLE
            IF(LRK.EQ.LRJ) CYCLE
            CALL Drl_TS_EXT(LRK)
          ENDDO
        ENDDO
      ENDDO
      RETURN
!ST(2)   ACT -C"-
102   IF(LINELP.NE.14.OR.NLG2.EQ.1) RETURN
      IF(JB_SYS.GT.0) CALL ST_Drl_ACT_C_EXT_AB_SGT0(3)
      IF(JB_SYS.GT.0) CALL ST_ArBl_ACT_C_EXT_AB_SGT0(3)
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        DO LRJ=LRI+1,NORB_DZ
          CALL ST1_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_EXT_TS(LRI,LRJ,1)
          CALL ST2_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_EXT_TS(LRI,LRJ,NK)
          CALL ST4_EXT(LRI,LRJ,NK,1)
          IF(NK.NE.0) CALL Ar_BL_EXT_TS(LRI,LRJ,NK)
          CALL ST4_EXT(LRI,LRJ,NK,-1)
          IF(NK.NE.0) CALL Ar_BL_EXT_TS(LRI,LRJ,NK)
        ENDDO
      ENDDO

      DO LRI=NORB_FRZ+1,NORB_DZ-1
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JML.OR.LMIJ.NE.JMR) CYCLE
!ST(2-5) (22)Drl(12)-          ACT -C"-
            IWDL=JUST(LRI,LRJ)
            IWDR=IWDL               !
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=0.D0
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_ST(5)
            ENDDO
            CALL Drl_TS_EXT(LRJ)
!ST(2-6) Drl(22)-C"(12)-       ACT -C"-
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=0.D0
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_ST(6)
            ENDDO
            CALL Drl_TS_EXT(LRI)
!----------------------------------------------------------
!ST(2-3) Ar(13)-C'(22)-Bl(32)-   ACT -C"-
!ST(2-3) Ar(13)-Bl(32)-C'(22)-   ACT -C"-
!ST(2-7) Drl(12)-C"(22)-         ACT -C"-
!----------------------------------------------------------
          ENDDO
        ENDDO
      GOTO 10
!=======================================================================
!TS(3) D&R^L-  ACT -C"-
103   IF(LINELP.NE.14.OR.NLG2.EQ.1) RETURN
      IF(JB_SYS.GT.0) THEN
        CALL TS_ArBl_ACT_C_EXT_AB_SGT0(3)
      ENDIF
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        DO LRJ=LRI+1,NORB_DZ
          CALL TS1_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_EXT_TS(LRI,LRJ,1)
          CALL TS2_EXT(LRI,LRJ,NK,1)
          IF(NK.NE.0) CALL Ar_BL_EXT_TS(LRI,LRJ,NK)
          CALL TS2_EXT(LRI,LRJ,NK,-1)
          IF(NK.NE.0) CALL Ar_BL_EXT_TS(LRI,LRJ,NK)
          CALL TS4_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_EXT_TS(LRI,LRJ,NK)
        ENDDO
      ENDDO
      GOTO 10
!======================================================
!ST1(4) Ar-Bl- Drl- ACT -C"-
104   IF(LINELP.NE.14.or.nlg2.ne.2) RETURN
        CALL STT_ArBl_ACT_C_EXT_AB_SGT1(3)
      RETURN
!=======================================================
!T1S(5) Ar-Bl- Drl ACT -C"-
105   IF(LINELP.NE.14.OR.NLG2.NE.2) RETURN
        CALL TTS_Drl_ACT_C_EXT_AB_SGT1(3)
        CALL TTS_ArBl_ACT_C_EXT_AB_SGT1(3)
      RETURN
!===========================================================
!SD1(8) Ar ACT -Bl-
108   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      IF(JB_SYS.GT.0) CALL SDD_AR_ACT_BL_SGT0(3,LRA)
      RETURN
!===========================================================
!TTTT(12) Ar-Bl- Drl- ACT -C"-
112   IF(LINELP.NE.14) RETURN
      IF(NLG2.EQ.1) THEN
        CALL TTTT_Drl_ACT_C_EXT_AB_SGT0(3)
      ELSE
        CALL TTTT_ArBl_ACT_C_EXT_AB_SGT0(3)
      ENDIF
      RETURN
!=============================================================
!T1D1(15) Ar- ACT -Bl-
115   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      CALL TTDD_AR_ACT_BL_SGT1(3,LRA)
      RETURN
!============================================================
!D1D1(20) Drl- Ar-Bl- ACT -C"-
120   IF(LINELP.NE.14) RETURN
      IF(NLG2.EQ.1) THEN
        CALL D1D1_Drl_ACT_C_EXT_AB_SGT0(3)
      ELSE
        CALL D1D1_ArBl_ACT_C_EXT_AB_SGT0(3)
      ENDIF
      RETURN
!=================================================================
!DD1(21) Ar-Bl- ACT -C"-
121   IF(LINELP.NE.14.or.nlg2.ne.2) RETURN
      CALL DD1_ArBl_ACT_C_EXT_AB_SGT0(3)
      RETURN
!=================================================================
!D1D(22) Ar-Bl- Drl- ACT -C"-
122   IF(LINELP.NE.14.or.nlg2.ne.2) RETURN
      CALL D1D_ArBl_ACT_C_EXT_AB_SGT0(3)
      CALL D1D_Drl_ACT_C_EXT_AB_SGT0(3)
      RETURN
!=================================================================
!D1V(24) Ar- ACT -Bl-
124   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      CALL D1V_Ar_ACT_Bl_EXT_AB_SGT0(3,LRA)
      RETURN
!=======================================================================
!TT(11) Drl-  ACT -C"-
111   IF(LINELP.NE.14) RETURN
      IF(NLG2.EQ.1) GOTO 1111
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        DO LRJ=LRI+1,NORB_DZ
          call TT1_EXT(LRI,LRJ,NK,1)
          CALL AR_BL_EXT_TS(LRI,LRJ,NK)
          call TT1_EXT(LRI,LRJ,NK,-1)
          CALL AR_BL_EXT_TS(LRI,LRJ,NK)
        ENDDO
      ENDDO
      RETURN

1111  W0TT2=W0_TT(2)
      W1TT2=W1_TT(2)
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JML.OR.LMIJ.NE.JMR) CYCLE
!TT(11-2) (22)Drl(22)-
!TT(11-2) Drl(22)-C"(22)-
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0TT2
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1TT2
            ENDDO
            IWDL=JUST(LRI,LRJ)
            IWDR=IWDL
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            CALL Drl_TS_EXT(LRI)
            CALL Drl_TS_EXT(LRJ)
          ENDDO
        ENDDO
      RETURN
!=======================================================================
!DD(19) ACT -C"- ....................................................
119   IF(LINELP.NE.14) RETURN
      IF(NLG2.EQ.1) GOTO 1191
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JMLR) CYCLE
          W0DD1=W0_DD(1)
          W1DD1=W1_DD(1)
          NI=MOD(LRJ-LRI,2)
          IF(NI.EQ.0) THEN
            W0DD1=-W0DD1
            W1DD1=-W1DD1
          ENDIF
          IF(LMI.EQ.JML.AND.LMJ.EQ.JMR) THEN
!DD(19-1) Ar(23)-Bl(32)-      ACT -C"-
            DO MPL=1,MTYPE
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1DD1
            ENDDO
            IWDL=JUD(LRI)
            IWDR=JUD(LRJ)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            CALL AR_BL_EXT_TS(LRI,LRJ,1)
          ENDIF
        ENDDO
      ENDDO
      RETURN

1191  W0DD2=W0_DD(2)
      W1DD2=W1_DD(2)
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        IF(LMI.NE.JML) CYCLE
!DD(19-2) Drl(22)-
        DO MPL=1,MTYPE
          VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0DD2
          VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1DD2
        ENDDO
        IWDL=JUD(LRI)
        IWDR=IWDL
        DO MPL=1,MHLP
          IWAL=LPNEW_LWEI(MPL)
          IWAR=LPNEW_RWEI(MPL)
          LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        ENDDO
        CALL Drl_TS_EXT(LRI)
      ENDDO
      GOTO 10
!=======================================================================
!VV(25) ACT -BL- ....................................................
125   GOTO 10
C=======================================================================
10    return
      end

      subroutine TS_ext_head_in_act()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      LOGIC_DH=.FALSE.
      LRAI=NLG1
      LRAJ=NLG2
!line=5 A&r-B&l<-->EXT
!5    CONTINUE
      IF(LINELP.EQ.5) THEN
        CALL AR_BL_EXT_TS(lrai,lraj,1)
      ENDIF
!line=9 D&rl<-->EXT
!9    CONTINUE
      IF(LINELP.EQ.9) THEN
        call Drl_TS_EXT(lrai)
      endif
      return
      end


      subroutine tt_drt_ci_new()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      common/lpdisk/idisk_lp,idisk_array(13)

      call external_space_plpmode_value_tt()
      idisk_lp=idisk_array(8)

      DO lpblock=1,lpblock_tt
        call read_lp()
        IpaeL=iml+9
        Ipae =imr+9
        call logicg_st(iml,imr,3,3)      ! irtype=3(T),lptype=5:ArBl-
        imlr=mul_tab(iml,imr)
C        if(imlr.eq.1)imspace=iml
        nvalue_space_ss=iseg_downwei(9+imlr)
        idownwei_g131415=iseg_downwei(9+iml)
        call get_jpty(jpadlr,jptyl,jptyr)
        call get_jp(jptyl,jml,jpadl,1)
        call get_jp(jptyr,jmr,jpad,1)
C        JMLR=MUL_TAB(JML,JMR)
        if(linelp.le.12)   then
          call tt_ext_head_in_act()
        else
          call tt_ext_head_in_dbl()
        endif
      enddo
      return
      end

      subroutine tt_ext_head_in_dbl()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      LOGIC_DH=.TRUE.
      JMLR=MUL_TAB(JML,JMR)
      LPOK=JPADLR
      GOTO(101,102,103,104,105,106,10,108,10,10,111,112,
     :     113,10,115,10,10,10,119,120,121,122,123,124,125,10),LPOK
!=======================================================================
!SD(6-1) ACT -B&L-
106   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      CALL SD_AR_ACT_BL(11,LRA)
      IF(JB_SYS.GT.0) CALL SD_AR_ACT_BL_SGT0(11,LRA)
      GOTO 10
!====================================================================
!TD(13) ACT -B&L-
113   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      CALL TD_AR_ACT_BL(11,LRA)
      GOTO 10
!=======================================================================
!DV(23) ACT -C'-..................................................
123   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        IF(LMI.NE.JMLR) CYCLE
        W0DV1=W0_DV(1)
        NI=MOD(NORB_DZ-LRI,2)
        IF(NI.EQ.1)W0DV1=-W0DV1
!DV(23-1) A&r(23)-
        IWDL=JUD(LRI)
        IWDR=0
        DO MPL=1,MHLP
          IWAL=LPNEW_LWEI(MPL)
          IWAR=LPNEW_RWEI(MPL)
          LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        ENDDO
        DO MPL=1,MTYPE
          VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0DV1
          VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W0DV1
        ENDDO
        CALL Ar_BL_EXT_TT(LRI,LRA,1)
      ENDDO
      GOTO 10
!=======================================================================
!SS(1)   ACT -C"-
!-------------------------------------------------------------------
!SS(1-1)  Ar(01)-Bl(32)-        ACT -C"-
!SS(1-3)  Ar(13)-Bl(20)-        ACT -C"-
!SS(1-6)  (11)-Ar(23)-Bl(32)-   ACT -C"-
!SS(1-7)  Ar(13)-C'(21)-Bl(32)- ACT -C"-
!SS(1-8)  Ar(13)-C'(22)-Bl(31)- ACT -C"-
!SS(1-9)  Ar(23)-C'(11)-Bl(32)- ACT -C"-
!SS(1-11) Ar(13)-Bl(31)-C"(22)- ACT -C"-
!SS(1-12) Ar(13)-Bl(32)-C"(21)- ACT -C"-
!SS(1-13) Ar(23)-Bl(31)-C"(12)- ACT -C"-
!SS(1-16) (11)-Drl(22)-         ACT -C"-
!SS(1-18) Drl(11)-C"(22)-       ACT -C"-
!SS(1-19) Drl(12)-C"(21)-       ACT -C"-
!SS(1-20) (11)-Drl(33)-C"(22)-  ACT -C"-
!SS(1-20) Drl(33)-C"(11)-C"(22)-ACT -C"-
!-------------------------------------------------------------------
!SS(1)   ACT 14: -C"-
101   IF(LINELP.NE.14) RETURN
      IF(NLG2.EQ.1) GOTO 1011
      IF(JB_SYS.GT.0) THEN
       CALL SS_ArBl_ACT_C_EXT_AB_SGT0(11)
       CALL SS_S_Drl_ACT_C_EXT_AB_SGT0(11)
      ENDIF
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        DO LRJ=LRI+1,NORB_DZ
          CALL SS2_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_EXT_TT(LRI,LRJ,1)
          CALL SS4_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_EXT_TT(LRI,LRJ,1)
          CALL SS5_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_EXT_TT(LRI,LRJ,NK)
          CALL SS10_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_EXT_TT(LRI,LRJ,NK)
          CALL SS14_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_EXT_TT(LRI,LRJ,NK)
        ENDDO
      ENDDO
      RETURN

1011  IF(JB_SYS.GT.0) CALL SS_Drl_ACT_C_EXT_AB_SGT0(11)
      W0SS15=W0_SS(15)
      W1SS15=W1_SS(15)
      W0SS17=W0_SS(17)
      W1SS17=W1_SS(17)
      W0SS20=W0_SS(20)
      IF(JML.EQ.1.AND.JMR.EQ.1) THEN
!SS(1-20) Drl(33)-C"(00)-       ACT -C"-                     ! IPL(R)AD=
        DO LR0=NORB_FRZ+1,NORB_DZ
          IWDL=JUST(LR0,LR0)
          IWDR=IWDL
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          ENDDO
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS20
            VPLP_W1(MPL)=0.D0
          ENDDO
C      do lrk=1,norb_dz
C        if(lrk.eq.lr0) cycle
C        CALL Drl_tt_ext(lrk)
C     enddo
          CALL Drl_TT_SUM(LR0,0)
        ENDDO
      ENDIF
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JML.OR.LMIJ.NE.JMR) CYCLE
          IWDL=JUST(LRI,LRJ)
          IWDR=IWDL
          DO MPL=1,MHLP
            IWAL=LPNEW_LWEI(MPL)
            IWAR=LPNEW_RWEI(MPL)
            LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          ENDDO
!SS(1-15) (22)-Drl(11)-         ACT -C"-
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS15
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SS15
          ENDDO
          CALL Drl_TT_EXT(LRJ)
!SS(1-17) Drl(22)-C"(11)-       ACT -C"-
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS17
            VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SS17
          ENDDO
          CALL Drl_TT_EXT(LRI)
!SS(1-20) (22)(11)Drl(33)-      ACT -C"-
!SS(1-20) (22)Drl(33)-C"(11)-   ACT -C"-
!SS(1-20) Drl(33)-C"(22)-C"(11)-ACT -C"-
          DO MPL=1,MTYPE
            VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS20
            VPLP_W1(MPL)=0.D0
          ENDDO
          CALL Drl_TT_SUM(LRI,LRJ)
        ENDDO
      ENDDO
      RETURN
!ST(2)   ACT -C"-
102   IF(LINELP.NE.14.OR.NLG2.EQ.1) RETURN
      IF(JB_SYS.GT.0) CALL ST_Drl_ACT_C_EXT_AB_SGT0(11)
      IF(JB_SYS.GT.0) CALL ST_ArBl_ACT_C_EXT_AB_SGT0(11)
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        DO LRJ=LRI+1,NORB_DZ
          CALL ST1_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_EXT_TT(LRI,LRJ,1)
          CALL ST2_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_EXT_TT(LRI,LRJ,NK)
          CALL ST4_EXT(LRI,LRJ,NK,1)
          IF(NK.NE.0) CALL Ar_BL_EXT_TT(LRI,LRJ,NK)
          CALL ST4_EXT(LRI,LRJ,NK,-1)
          IF(NK.NE.0) CALL Ar_BL_EXT_TT(LRI,LRJ,NK)
        ENDDO
      ENDDO

      DO LRI=NORB_FRZ+1,NORB_DZ-1
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JML.OR.LMIJ.NE.JMR) CYCLE
!ST(2-5) (22)Drl(12)-          ACT -C"-
            IWDL=JUST(LRI,LRJ)
            IWDR=IWDL           !
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=0.D0
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_ST(5)
            ENDDO
            CALL Drl_TT_EXT(LRJ)
!ST(2-6) Drl(22)-C"(12)-       ACT -C"-
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=0.D0
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1_ST(6)
            ENDDO
            CALL Drl_TT_EXT(LRI)
!----------------------------------------------------------
!ST(2-3) Ar(13)-C'(22)-Bl(32)-   ACT -C"-
!ST(2-3) Ar(13)-Bl(32)-C'(22)-   ACT -C"-
!ST(2-7) Drl(12)-C"(22)-         ACT -C"-
!----------------------------------------------------------
          ENDDO
        ENDDO
      GOTO 10
!=======================================================================
!TS(3) D&R^L-  ACT -C"-
103   IF(LINELP.NE.14.OR.NLG2.EQ.1) RETURN
      IF(JB_SYS.GT.0) THEN
        CALL TS_ArBl_ACT_C_EXT_AB_SGT0(11)
      ENDIF
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        DO LRJ=LRI+1,NORB_DZ
          CALL TS1_EXT(LRI,LRJ,NK)
         IF(NK.NE.0) CALL Ar_BL_EXT_TT(LRI,LRJ,1)
          CALL TS2_EXT(LRI,LRJ,NK,1)
          IF(NK.NE.0) CALL Ar_BL_EXT_TT(LRI,LRJ,NK)
          CALL TS2_EXT(LRI,LRJ,NK,-1)
          IF(NK.NE.0) CALL Ar_BL_EXT_TT(LRI,LRJ,NK)
          CALL TS4_EXT(LRI,LRJ,NK)
          IF(NK.NE.0) CALL Ar_BL_EXT_TT(LRI,LRJ,NK)
        ENDDO
      ENDDO
      GOTO 10
!======================================================
!ST1(4) Ar-Bl- Drl- ACT -C"-
104   IF(LINELP.NE.14.or.nlg2.ne.2) RETURN
        CALL STT_ArBl_ACT_C_EXT_AB_SGT1(11)
      RETURN
!=======================================================
!T1S(5) Ar-Bl- Drl ACT -C"-
105   IF(LINELP.NE.14.OR.NLG2.NE.2) RETURN
        CALL TTS_Drl_ACT_C_EXT_AB_SGT1(11)
        CALL TTS_ArBl_ACT_C_EXT_AB_SGT1(11)
      RETURN
!===========================================================
!SD1(8) Ar ACT -Bl-
108   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      IF(JB_SYS.GT.0) CALL SDD_AR_ACT_BL_SGT0(11,LRA)
      RETURN
!===========================================================
!TTTT(12) Ar-Bl- Drl- ACT -C"-
112   IF(LINELP.NE.14) RETURN
      IF(NLG2.EQ.1) THEN
        CALL TTTT_Drl_ACT_C_EXT_AB_SGT0(11)
      ELSE
        CALL TTTT_ArBl_ACT_C_EXT_AB_SGT0(11)
      ENDIF
      RETURN
!=============================================================
!T1D1(15) Ar- ACT -Bl-
115   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      CALL TTDD_AR_ACT_BL_SGT1(11,LRA)
      RETURN
!============================================================
!D1D1(20) Drl- Ar-Bl- ACT -C"-
120   IF(LINELP.NE.14) RETURN
      IF(NLG2.EQ.1) THEN
        CALL D1D1_Drl_ACT_C_EXT_AB_SGT0(11)
      ELSE
        CALL D1D1_ArBl_ACT_C_EXT_AB_SGT0(11)
      ENDIF
      RETURN
!=================================================================
!DD1(21) Ar-Bl- ACT -C"-
121   IF(LINELP.NE.14.or.nlg2.ne.2) RETURN
      CALL DD1_ArBl_ACT_C_EXT_AB_SGT0(11)
      RETURN
!=================================================================
!D1D(22) Ar-Bl- Drl- ACT -C"-
122   IF(LINELP.NE.14.or.nlg2.ne.2) RETURN
      CALL D1D_ArBl_ACT_C_EXT_AB_SGT0(11)
      CALL D1D_Drl_ACT_C_EXT_AB_SGT0(11)
      RETURN
!=================================================================
!D1V(24) Ar- ACT -Bl-
124   IF(LINELP.NE.17) RETURN
      LRA=NLG1
      CALL D1V_Ar_ACT_Bl_EXT_AB_SGT0(11,LRA)
      RETURN
!=======================================================================
!TT(11) Drl-  ACT -C"-
111   IF(LINELP.NE.14) RETURN
      IF(NLG2.EQ.1) GOTO 1111
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        DO LRJ=LRI+1,NORB_DZ
          call TT1_EXT(LRI,LRJ,NK,1)
          CALL Ar_BL_EXT_TT(LRI,LRJ,NK)
          call TT1_EXT(LRI,LRJ,NK,-1)
          CALL Ar_BL_EXT_TT(LRI,LRJ,NK)
        ENDDO
      ENDDO
      RETURN

1111  W0TT2=W0_TT(2)
      W1TT2=W1_TT(2)
      W0TT3=W0_TT(3)
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JML.OR.LMIJ.NE.JMR) CYCLE
!TT(11-2) (22)Drl(22)-
!TT(11-2) Drl(22)-C"(22)-
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0TT2
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1TT2
            ENDDO
            IWDL=JUST(LRI,LRJ)
            IWDR=IWDL
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            CALL Drl_TT_EXT(LRI)
            CALL Drl_TT_EXT(LRJ)
!TT(11-3) Drl(33)-C"(22)-C"(22)-
!TT(11-3) (22)Drl(33)-C"(22)-
!TT(11-3) (22)(22)Drl(33)-
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0TT3
              VPLP_W1(MPL)=0.D0
            ENDDO
C      do lrk=1,norb_dz
C        if(lrk.eq.lri) cycle
C        if(lrk.eq.lrj) cycle
C        CALL Drl_tt_ext(lrk)
C     enddo
            CALL Drl_TT_SUM(LRI,LRJ)
          ENDDO
        ENDDO
      RETURN
!=======================================================================
!DD(19) ACT -C"- ....................................................
119   IF(LINELP.NE.14) RETURN
      IF(NLG2.EQ.1) GOTO 1191
      DO LRI=NORB_FRZ+1,NORB_DZ-1
        LMI=LSM_INN(LRI)
        DO LRJ=LRI+1,NORB_DZ
          LMJ=LSM_INN(LRJ)
          LMIJ=MUL_TAB(LMI,LMJ)
          IF(LMIJ.NE.JMLR) CYCLE
          W0DD1=W0_DD(1)
          W1DD1=W1_DD(1)
          NI=MOD(LRJ-LRI,2)
          IF(NI.EQ.0) THEN
            W0DD1=-W0DD1
            W1DD1=-W1DD1
          ENDIF
          IF(LMI.EQ.JML.AND.LMJ.EQ.JMR) THEN
!DD(19-1) Ar(23)-Bl(32)-      ACT -C"-
            DO MPL=1,MTYPE
              VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0DD1
              VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1DD1
            ENDDO
            IWDL=JUD(LRI)
            IWDR=JUD(LRJ)
            DO MPL=1,MHLP
              IWAL=LPNEW_LWEI(MPL)
              IWAR=LPNEW_RWEI(MPL)
              LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            ENDDO
            CALL Ar_BL_EXT_TT(LRI,LRJ,1)
          ENDIF
        ENDDO
      ENDDO
      RETURN

1191  W0DD2=W0_DD(2)
      W1DD2=W1_DD(2)
      W0DD3=W0_DD(3)
      DO LRI=NORB_FRZ+1,NORB_DZ
        LMI=LSM_INN(LRI)
        IF(LMI.NE.JML) CYCLE
!DD(19-2) Drl(22)-
        DO MPL=1,MTYPE
          VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0DD2
          VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1DD2
        ENDDO
        IWDL=JUD(LRI)
        IWDR=IWDL
        DO MPL=1,MHLP
          IWAL=LPNEW_LWEI(MPL)
          IWAR=LPNEW_RWEI(MPL)
          LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        ENDDO
        CALL Drl_TT_EXT(LRI)
!DD(19-3) (22)Drl(33)-
!DD(19-3) Drl(33)-C"(22)-
        DO MPL=1,MTYPE
          VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0DD3
          VPLP_W1(MPL)=0.D0
        ENDDO
        CALL Drl_TT_SUM(LRI,0)
C      do lrk=1,norb_dz
C        if(lrk.eq.lri) cycle
C        CALL Drl_tt_ext(lrk)
C     enddo
      ENDDO
      GOTO 10
!=======================================================================
!VV(25) ACT -BL- ....................................................
125   IF(LINELP.NE.14.OR.NLG2.NE.1) RETURN
!VV(25) Drl(33)-
      IWDL=0
      IWDR=0
      DO MPL=1,MTYPE
        VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0_VV
        VPLP_W1(MPL)=0.D0
      ENDDO
      DO MPL=1,MHLP
        IWAL=LPNEW_LWEI(MPL)
        IWAR=LPNEW_RWEI(MPL)
        LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      ENDDO
!      CALL Drl_TT_SUM(0,0)
      do lrk=1,norb_dz
        CALL Drl_tt_ext(lrk)
      enddo
      GOTO 10
C=======================================================================
10    return
      end

      subroutine TT_ext_head_in_act()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      LOGIC_DH=.FALSE.
      LRAI=NLG1
      LRAJ=NLG2
!line=5 A&r-B&l<-->EXT
!5    CONTINUE
      IF(LINELP.EQ.5) THEN
        CALL Ar_BL_EXT_TT(lrai,lraj,1)
      ENDIF
!line=9 D&rl<-->EXT
!9    CONTINUE
      IF(LINELP.EQ.9) THEN
        call Drl_TT_EXT(lrai)
      endif
      return
      end


      subroutine logicg_st(ilnodesm,irnodesm,iltype,irtype)
#include "drt_h.fh"
#include "lpextmode_h.fh"
      ilrsm=mul_tab(ilnodesm,irnodesm)
      iii=1      !index to determine lwei rwei iposint and nlinkorb
!G2G4a G2G4b G1415 G13
      logic_g36a=.false.
      logic_g36b=.false.
      logic_g35a=.false.
      logic_g35b=.false.
      logic_g34a=.false.
      logic_g34b=.false.

      lpsta36a=iii
      call do_g36mode(ilrsm,ilnodesm,iii)
      lpend36a=iii-4
      if ( lpend36a .ge. lpsta36a ) logic_g36a=.true.
      lpsta35a=iii
      call do_g35mode(ilrsm,ilnodesm,iii)
      lpend35a=iii-4
      if ( lpend35a .ge. lpsta35a ) logic_g35a=.true.
      lpsta34a=iii
      call do_g34mode(ilrsm,ilnodesm,iii)
      lpend34a=iii-4
      if ( lpend34a .ge. lpsta34a ) logic_g34a=.true.
      if ( ilrsm .ne. 1 ) then
            lpsta36b=iii
            call do_g36mode(ilrsm,irnodesm,iii)
            lpend36b=iii-4
            if ( lpend36b .ge. lpsta36b ) logic_g36b=.true.
            lpsta35b=iii
            call do_g35mode(ilrsm,irnodesm,iii)
            lpend35b=iii-4
            if ( lpend35b .ge. lpsta35b ) logic_g35b=.true.
            lpsta34b=iii
            call do_g34mode(ilrsm,irnodesm,iii)
            lpend34b=iii-4
            if ( lpend34b .ge. lpsta34b ) logic_g34b=.true.
      else
            logic_g36b=logic_g36a
            lpsta36b=lpsta36a
            lpend36b=lpend36a
            logic_g35b=logic_g35a
            lpsta35b=lpsta35a
            lpend35b=lpend35a
            logic_g34b=logic_g34a
            lpsta34b=lpsta34a
            lpend34b=lpend34a
      endif
      logic_g2g4a=.false.
      logic_g2g4b=.false.
      logic_g1415=.false.
      logic_g13=.false.

      if ( irnodesm .eq. 1 .and. irtype .eq. 4 ) then        ! ir:S_1
            logic_g2g4a=.true.
            ism_g2g4=ilnodesm
      endif
      if ( ilnodesm .eq. 1 .and. iltype .eq. 4 ) then        ! il:S_1
            logic_g2g4b=.true.
            ism_g2g4=irnodesm
      endif
      if ( ilnodesm .eq. irnodesm) then        ! ir:S_1
        logic_g1415=.true.
        ism_g1415=ilnodesm
      endif
      end

!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

! calculate dd, ss and tt space
subroutine dd_drt_ci_new()

use gugaci_global, only: idisk_array, idisk_lp, iml, imr, int_dd_drl, int_dd_offset, ipae, ipael, jml, jmr, jpad, jpadl, jpadlr, &
                         linelp, lpblock_dd, v_onevsqtwo, v_sqthreevsqtwo, w0gdd, w1gdd
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: jptyl, jptyr, lpb

w0gdd = v_onevsqtwo
w1gdd = -v_sqthreevsqtwo
idisk_lp = idisk_array(3)

do lpb=1,lpblock_dd
  call read_lp()
  IpaeL = iml+1
  Ipae = imr+1
  int_dd_drl = int_dd_offset(iml,imr)
  call logicg_dd(iml,imr)
  call get_jpty(jpadlr,jptyl,jptyr)
  call get_jp(jptyl,jml,jpadl,1)
  call get_jp(jptyr,jmr,jpad,1)
  !JMLR = Mul(JML,JMR)
  if (linelp <= 12) then
    call dd_ext_head_in_act()
  else
    call dd_ext_head_in_dbl()
  end if
end do

return

end subroutine dd_drt_ci_new

subroutine dd_ext_head_in_dbl()

use gugaci_global, only: ipae, ipael, jb_sys, jml, jmr, jpad, jpadl, jpadlr, jud, just, linelp, logic_dh, lp_lwei, lp_rwei, &
                         lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, nlg1, nlg2, norb_dz, norb_frz, vplp_w0, vplp_w1, &
                         vplpnew_w0, vplpnew_w1, w0_dd, w0_dv, w0_sd, w0_ss, w0_td, w0_tt, w0_vv, w1_dd, w1_ss, w1_st, w1_tt
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, jmlr, lmi, lmij, lmj, lmk, lmki, lpok, lr0, lra, lri, lrj, lrk, mpl, ni, nk
real(kind=wp) :: w0dd1, w0dd2, w0dd3, w0dv1, w0sd1, w0sd2, w0sd3, w0sd4, w0ss15, w0ss17, w0ss20, w0td1, w0tt2, w0tt3, w1dd1, &
                 w1dd2, w1ss15, w1ss17, w1tt2
integer(kind=iwp), external :: iwalk_ad

LOGIC_DH = .true.
JMLR = Mul(JML,JMR)
LPOK = JPADLR
select case (LPOK)
  case (1)
    !===================================================================
    ! SS(1)   ACT -C"-
    !-------------------------------------------------------------------
    ! SS(1-9)  Ar(23)-C'(11)-Bl(32)- ACT -C"-
    ! SS(1-11) Ar(13)-Bl(31)-C"(22)- ACT -C"-
    ! SS(1-12) Ar(13)-Bl(32)-C"(21)- ACT -C"-
    ! SS(1-13) Ar(23)-Bl(31)-C"(12)- ACT -C"-
    ! SS(1-16) (11)-Drl(22)-         ACT -C"-
    ! SS(1-18) Drl(11)-C"(22)-       ACT -C"-
    ! SS(1-19) Drl(12)-C"(21)-       ACT -C"-
    ! SS(1-20) (11)-Drl(33)-C"(22)-  ACT -C"-
    ! SS(1-20) Drl(33)-C"(11)-C"(22)-ACT -C"-
    !-------------------------------------------------------------------
    ! SS(1)   ACT 14: -C"-
    if (LINELP /= 14) return
    if (NLG2 == 1) then
      if (JB_SYS > 0) then
        LRA = NLG1
        call SS_Drl_ACT_C_DD_EXT_SGT0()
      end if
      W0SS15 = W0_SS(15)
      W1SS15 = W1_SS(15)
      W0SS17 = W0_SS(17)
      W1SS17 = W1_SS(17)
      W0SS20 = W0_SS(20)
      if ((JML == 1) .and. (JMR == 1)) then
        ! SS(1-20) Drl(33)-C"(00)-       ACT -C"-
        do LR0=NORB_FRZ+1,NORB_DZ
          IWDL = JUST(LR0,LR0)
          IWDR = IWDL
          do MPL=1,MHLP
            IWAL = LPNEW_LWEI(MPL)
            IWAR = LPNEW_RWEI(MPL)
            LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          end do
          do MPL=1,MTYPE
            VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS20
            VPLP_W1(MPL) = Zero
          end do
          do LRK=1,NORB_DZ
            if (LRK == LR0) cycle
            call Drl_DD_EXT(LRK)
          end do
        end do
      end if
      do LRI=NORB_FRZ+1,NORB_DZ-1
        LMI = LSM_INN(LRI)
        do LRJ=LRI+1,NORB_DZ
          LMJ = LSM_INN(LRJ)
          LMIJ = Mul(LMI,LMJ)
          if ((LMIJ /= JML) .or. (LMIJ /= JMR)) cycle
          IWDL = JUST(LRI,LRJ)
          IWDR = IWDL
          do MPL=1,MHLP
            IWAL = LPNEW_LWEI(MPL)
            IWAR = LPNEW_RWEI(MPL)
            LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          end do
          ! SS(1-15) (22)-Drl(11)-         ACT -C"-
          do MPL=1,MTYPE
            VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS15
            VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SS15
          end do
          call Drl_DD_EXT(LRJ)
          ! SS(1-17) Drl(22)-C"(11)-       ACT -C"-
          do MPL=1,MTYPE
            VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS17
            VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SS17
          end do
          call Drl_DD_EXT(LRI)
          ! SS(1-20) (22)(11)Drl(33)-      ACT -C"-
          ! SS(1-20) (22)Drl(33)-C"(11)-   ACT -C"-
          ! SS(1-20) Drl(33)-C"(22)-C"(11)-ACT -C"-
          do MPL=1,MTYPE
            VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS20
            VPLP_W1(MPL) = Zero
          end do
          do LRK=1,NORB_DZ
            if (LRK == LRI) cycle
            if (LRK == LRJ) cycle
            call Drl_DD_EXT(LRK)
          end do
        end do
      end do
    else
      do LRI=NORB_FRZ+1,NORB_DZ-1
        do LRJ=LRI+1,NORB_DZ
          call SS2_EXT(LRI,LRJ,NK)
          if (NK /= 0) call Ar_BL_DD_EXT(LRI,LRJ,1)
          call SS4_EXT(LRI,LRJ,NK)
          if (NK /= 0) call Ar_BL_DD_EXT(LRI,LRJ,1)
          call SS5_EXT(LRI,LRJ,NK)
          if (NK /= 0) call Ar_BL_DD_EXT(LRI,LRJ,NK)
          call SS10_EXT(LRI,LRJ,NK)
          if (NK /= 0) call Ar_BL_DD_EXT(LRI,LRJ,NK)
          call SS14_EXT(LRI,LRJ,NK)
          if (NK /= 0) call Ar_BL_DD_EXT(LRI,LRJ,NK)
          !if (JROUTE_SYS > 1) then
          !  call SS1_EXT(LRI,LRJ)
          !  call Ar_BL_DD_EXT(LRI,LRJ,1)
          !  call SS3_EXT(LRI,LRJ)
          !  call Ar_BL_DD_EXT(LRI,LRJ,1)
          !end if
        end do
      end do
      if (JB_SYS > 0) then
        call SS_ArBr_ACT_C_DD_EXT_SGT0()
        call SS_S_Drl_ACT_C_DD_EXT_SGT0()
      end if
    end if

  case (2)
    !===================================================================
    ! ST(2)   ACT -C"-
    !if (JB_SYS > 0) call ST_Drl_ACT_C_DD_EXT()
    !if (JB_SYS > 0) call ST_ArBl_ACT_C_DD_EXT_SGT0()

    if ((LINELP /= 14) .or. (NLG2 == 1)) return
    if (JB_SYS > 0) call ST_Drl_ACT_C_DD_EXT_SGT0()
    if (JB_SYS > 0) call ST_ArBl_ACT_C_DD_EXT_SGT0()
    do LRI=NORB_FRZ+1,NORB_DZ-1
      do LRJ=LRI+1,NORB_DZ
        call ST1_EXT(LRI,LRJ,NK)
        if (NK /= 0) call Ar_BL_DD_EXT(LRI,LRJ,1)
        call ST2_EXT(LRI,LRJ,NK)
        if (NK /= 0) call Ar_BL_DD_EXT(LRI,LRJ,NK)
        call ST4_EXT(LRI,LRJ,NK,1)
        if (NK /= 0) call Ar_BL_DD_EXT(LRI,LRJ,NK)
        call ST4_EXT(LRI,LRJ,NK,-1)
        if (NK /= 0) call Ar_BL_DD_EXT(LRI,LRJ,NK)
      end do
    end do

    do LRI=NORB_FRZ+1,NORB_DZ-1
      LMI = LSM_INN(LRI)
      do LRJ=LRI+1,NORB_DZ
        LMJ = LSM_INN(LRJ)
        LMIJ = Mul(LMI,LMJ)
        if ((LMIJ /= JML) .or. (LMIJ /= JMR)) cycle
        ! ST(2-5) (22)Drl(12)-          ACT -C"-
        IWDL = JUST(LRI,LRJ)
        IWDR = JUST(LRI,LRJ)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = Zero
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_ST(5)
        end do
        call Drl_DD_EXT(LRJ)
        ! ST(2-6) Drl(22)-C"(12)-       ACT -C"-
        do MPL=1,MTYPE
          VPLP_W0(MPL) = Zero
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_ST(6)
        end do
        call Drl_DD_EXT(LRI)
        !---------------------------------------------------------------
        ! ST(2-3) Ar(13)-C'(22)-Bl(32)-   ACT -C"-
        ! ST(2-3) Ar(13)-Bl(32)-C'(22)-   ACT -C"-
        ! ST(2-7) Drl(12)-C"(22)-         ACT -C"-
        !---------------------------------------------------------------
      end do
    end do

  case (3)
    !===================================================================
    ! TS(3) D&R^L-  ACT -C"-
    if ((LINELP /= 14) .or. (NLG2 == 1)) return
    if (JB_SYS > 0) then
      call TS_ArBl_ACT_C_DD_EXE_SGT0()
    end if
    do LRI=NORB_FRZ+1,NORB_DZ-1
      do LRJ=LRI+1,NORB_DZ
        call TS1_EXT(LRI,LRJ,NK)
        if (NK /= 0) call Ar_BL_DD_EXT(LRI,LRJ,1)
        call TS2_EXT(LRI,LRJ,NK,1)
        if (NK /= 0) call Ar_BL_DD_EXT(LRI,LRJ,NK)
        call TS2_EXT(LRI,LRJ,NK,-1)
        if (NK /= 0) call Ar_BL_DD_EXT(LRI,LRJ,NK)
        call TS4_EXT(LRI,LRJ,NK)
        if (NK /= 0) call Ar_BL_DD_EXT(LRI,LRJ,NK)
      end do
    end do

  case (4)
    !===================================================================
    ! STT(4) ArBl-  ACT -C"-
    if ((LINELP /= 14) .or. (NLG2 == 1)) return
    LRA = NLG1
    call STT_ArBl_ACT_C_DD_EXT_SGT1()

  case (5)
    !===================================================================
    ! TTS(3) ArBl-  ACT -C"-
    if ((LINELP /= 14) .or. (NLG2 == 1)) return
    call TTS_ArBl_ACT_C_DD_EXT_SGT1()
    call TTS_Drl_ACT_C_DD_EXT_SGT1()

  case default ! (6)
    !===================================================================
    ! SD(6)    ACT: -B&L-
    if (LINELP /= 17) return
    LRA = NLG1
    if (JB_SYS > 0) then
      call SD_Ar_ACT_Bl_DD_EXT_SGT0(LRA)
    end if
    do LRI=NORB_FRZ+1,NORB_DZ
      LMI = LSM_INN(LRI)
      if (LMI /= JMLR) cycle
      W0SD1 = W0_SD(1)
      W0SD2 = W0_SD(2)
      W0SD3 = -W0_SD(3)
      W0SD4 = -W0_SD(4)
      NI = mod(NORB_DZ-LRI,2)
      if (NI == 1) W0SD1 = -W0SD1
      if (NI == 1) W0SD2 = -W0SD2
      if (NI == 1) W0SD3 = -W0SD3
      if (NI == 1) W0SD4 = -W0SD4
      if ((JML == 1) .and. (LMI == JMR)) then
        ! SD(6-1) A&r(02)-
        IWDL = JUST(LRI,LRI)
        IWDR = JUD(LRI)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SD1
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W0SD1
        end do
        call Ar_BL_DD_EXT(LRI,LRA,1)
      end if
      ! SD(6-2) C(22)-A&r(13)-
      do LRK=NORB_FRZ+1,LRI-1
        LMK = LSM_INN(LRK)
        LMKI = Mul(LMK,LMI)
        if ((LMKI == JML) .and. (LMK == JMR)) then
          IWDL = JUST(LRK,LRI)
          IWDR = JUD(LRK)
          do MPL=1,MHLP
            IWAL = LPNEW_LWEI(MPL)
            IWAR = LPNEW_RWEI(MPL)
            LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          end do
          do MPL=1,MTYPE
            VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SD2
            VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W0SD2
          end do
          call Ar_BL_DD_EXT(LRI,LRA,1)
        end if
      end do
      ! SD(6-4) A&r(23)-C'(12)-
      do LRK=LRI+1,NORB_DZ
        LMK = LSM_INN(LRK)
        LMKI = Mul(LMK,LMI)
        if ((LMKI /= JML) .or. (LMK /= JMR)) cycle
        !......................03_01....................................
        !if (jroute_sys > 1) then
        !  ! SD(6-3) A&r(13)-C'(22)-
        !  IWDL = JUST(LRK,LRI)
        !  IWDR = JUD(LRK)
        !  do MPL=1,MHLP
        !    IWAL = LPNEW_LWEI(MPL)
        !    IWAR = LPNEW_RWEI(MPL)
        !    LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        !    LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        !  end do
        !  do MPL=1,MTYPE
        !    VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SD3
        !    VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W0SD3
        !  end do
        !  call Ar_BL_DD_EXT(LRI,LRA,1)
        !end if
        !.......................03_01...................................
        IWDL = JUST(LRI,LRK)
        IWDR = JUD(LRK)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SD4
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W0SD4
        end do
        call Ar_BL_DD_EXT(LRI,LRA,1)
      end do
    end do

  case (8)
    !===================================================================
    ! SDD(8) Ar- ACT -Bl-
    if (LINELP /= 17) return
    LRA = NLG1
    call SDD_Ar_ACT_Bl_DD_EXT_SGT0(LRA)

  case (11)
    !===================================================================
    ! TT(11-1) (22)Ar(23)-Bl(32)-      IGF=1   ACT -C"-
    ! TT(11-1) Ar(23)-Bl(32)-C"(22)-  IGF=1    ACT -C"-
    ! TT(11-1) Ar(23)-C'(22)-Bl(32)-  IGF=-1    ACT -C"-
    if (LINELP /= 14) return
    if (NLG2 == 1) then
      W0TT2 = W0_TT(2)
      W1TT2 = W1_TT(2)
      W0TT3 = W0_TT(3)
      do LRI=NORB_FRZ+1,NORB_DZ-1
        LMI = LSM_INN(LRI)
        do LRJ=LRI+1,NORB_DZ
          LMJ = LSM_INN(LRJ)
          LMIJ = Mul(LMI,LMJ)
          if ((LMIJ /= JML) .or. (LMIJ /= JMR)) cycle
          ! TT(11-2) (22)Drl(22)-
          ! TT(11-2) Drl(22)-C"(22)-
          do MPL=1,MTYPE
            VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0TT2
            VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1TT2
          end do
          IWDL = JUST(LRI,LRJ)
          IWDR = IWDL
          do MPL=1,MHLP
            IWAL = LPNEW_LWEI(MPL)
            IWAR = LPNEW_RWEI(MPL)
            LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          end do
          call Drl_DD_EXT(LRI)
          call Drl_DD_EXT(LRJ)
          do MPL=1,MTYPE
            VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0TT3
            VPLP_W1(MPL) = Zero
          end do
          do LRK=1,NORB_DZ
            if (LRK == LRI) cycle
            if (LRK == LRJ) cycle
            LMK = LSM_INN(LRK)
            LMKI = Mul(LMK,LMI)
            ! TT(11-3) Drl(33)-C"(22)-C"(22)-
            ! TT(11-3) (22)Drl(33)-C"(22)-
            ! TT(11-3) (22)(22)Drl(33)-
            call Drl_DD_EXT(LRK)
          end do
        end do
      end do
    else
      do LRI=NORB_FRZ+1,NORB_DZ-1
        do LRJ=LRI+1,NORB_DZ
          call TT1_EXT(LRI,LRJ,NK,1)
          call Ar_BL_DD_EXT(LRI,LRJ,NK)
          call TT1_EXT(LRI,LRJ,NK,-1)
          call Ar_BL_DD_EXT(LRI,LRJ,NK)
        end do
      end do
    end if

  case (12)
    !===================================================================
    ! TTTT(12) Drl- ArBl-  ACT -C"-
    if (LINELP /= 14) return
    if (NLG2 == 1) then
      call TTTT_Drl_ACT_C_DD_EXT_SGT1()
    else
      call TTTT_ArBl_ACT_C_DD_EXT_SGT1()
    end if

  case (13)
    !===================================================================
    ! TD(13) ACT -B&L-
    if (LINELP /= 17) return
    LRA = NLG1
    do LRI=NORB_FRZ+1,NORB_DZ
      LMI = LSM_INN(LRI)
      if (LMI /= JMLR) cycle
      W0TD1 = W0_TD(1)
      NI = mod(NORB_DZ-LRI,2)
      if (NI == 1) W0TD1 = -W0TD1

      do LRK=NORB_FRZ+1,LRI-1
        LMK = LSM_INN(LRK)
        if (LMK == JMR) then
          ! TD(13-1) C(22)-A&r(23)-
          IWDL = JUST(LRK,LRI)
          IWDR = JUD(LRK)
          do MPL=1,MHLP
            IWAL = LPNEW_LWEI(MPL)
            IWAR = LPNEW_RWEI(MPL)
            LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          end do
          do MPL=1,MTYPE
            VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0TD1
            VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W0TD1
          end do
          call Ar_BL_DD_EXT(LRI,LRA,1)
        end if
      end do
      do LRK=LRI+1,NORB_DZ
        LMK = LSM_INN(LRK)
        if (LMK == JMR) then
          ! TD(13-1) A&r(23)-C'(22)-
          IWDL = JUST(LRI,LRK)
          IWDR = JUD(LRK)
          do MPL=1,MHLP
            IWAL = LPNEW_LWEI(MPL)
            IWAR = LPNEW_RWEI(MPL)
            LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          end do
          do MPL=1,MTYPE
            VPLP_W0(MPL) = -VPLPNEW_W0(MPL)*W0TD1
            VPLP_W1(MPL) = -VPLPNEW_W1(MPL)*W0TD1
          end do
          call Ar_BL_DD_EXT(LRI,LRA,1)
        end if
      end do
    end do

  case (15)
    !===================================================================
    ! T1D1(15) Ar- ACT -Bl-
    if (LINELP /= 17) return
    LRA = NLG1
    call TTDD_Ar_ACT_Bl_DD_EXT_SGT1(LRA)

  case (19)
    !===================================================================
    ! DD(19) ACT -C"- ..................................................
    if (LINELP /= 14) return
    if (NLG2 == 1) then
      W0DD2 = W0_DD(2)
      W1DD2 = W1_DD(2)
      W0DD3 = W0_DD(3)
      do LRI=NORB_FRZ+1,NORB_DZ
        LMI = LSM_INN(LRI)
        if (LMI /= JML) cycle
        ! DD(19-2) Drl(22)-
        do MPL=1,MTYPE
          VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0DD2
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1DD2
        end do
        IWDL = JUD(LRI)
        IWDR = IWDL
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        end do
        call Drl_DD_EXT(LRI)
        do MPL=1,MTYPE
          VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0DD3
          VPLP_W1(MPL) = Zero
        end do
        ! DD(19-3) (22)Drl(33)-
        ! DD(19-3) Drl(33)-C"(22)-
        do LRK=1,NORB_DZ
          if (LRK == LRI) cycle
          LMK = LSM_INN(LRK)
          LMKI = Mul(LMK,LMI)
          call Drl_DD_EXT(LRK)
        end do
      end do
    else
      do LRI=NORB_FRZ+1,NORB_DZ-1
        LMI = LSM_INN(LRI)
        do LRJ=LRI+1,NORB_DZ
          LMJ = LSM_INN(LRJ)
          LMIJ = Mul(LMI,LMJ)
          if (LMIJ /= JMLR) cycle
          W0DD1 = W0_DD(1)
          W1DD1 = W1_DD(1)
          NI = mod(LRJ-LRI,2)
          if (NI == 0) then
            W0DD1 = -W0DD1
            W1DD1 = -W1DD1
          end if
          if ((LMI == JML) .and. (LMJ == JMR)) then
            ! DD(19-1) Ar(23)-Bl(32)-      ACT -C"-
            do MPL=1,MTYPE
              VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0DD1
              VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1DD1
            end do
            IWDL = JUD(LRI)
            IWDR = JUD(LRJ)
            do MPL=1,MHLP
              IWAL = LPNEW_LWEI(MPL)
              IWAR = LPNEW_RWEI(MPL)
              LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            end do
            call Ar_BL_DD_EXT(LRI,LRJ,1)
          end if
        end do
      end do
    end if

  case (20)
    !===================================================================
    ! DDDD(19) ACT -C"- ................................................
    if (LINELP /= 14) return
    if (NLG2 == 1) then
      call DDDD_Drl_ACT_C_DD_EXT_SGT0()
    else
      call DDDD_ArBl_ACT_C_DD_EXT_SGT0()
    end if

  case (21)
    !===================================================================
    ! DD1(19) ACT -C"- .................................................
    if ((LINELP /= 14) .or. (NLG2 /= 2)) return
    call DD1_ArBl_ACT_C_DD_EXT_SGT0()

  case (22)
    !===================================================================
    ! D1D(19) ACT -C"- .................................................
    if (LINELP /= 14) return
    !if (NLG2 == 1) then
    call D1D_Drl_ACT_C_DD_EXT_SGT0()
    !else
    call D1D_ArBl_ACT_C_DD_EXT_SGT0()
    !end if

  case (23)
    !===================================================================
    ! DV(23) ACT -C'-..................................................
    if (LINELP /= 17) return
    LRA = NLG1
    do LRI=NORB_FRZ+1,NORB_DZ
      LMI = LSM_INN(LRI)
      if (LMI /= JMLR) cycle
      W0DV1 = W0_DV(1)
      NI = mod(NORB_DZ-LRI,2)
      if (NI == 1) W0DV1 = -W0DV1
      ! DV(23-1) A&r(23)-
      IWDL = JUD(LRI)
      IWDR = 0
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      do MPL=1,MTYPE
        VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0DV1
        VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W0DV1
      end do
      call Ar_BL_DD_EXT(LRI,LRA,1)
    end do

  case (24)
    !===================================================================
    ! D1V(24) Ar- ACT -Bl-..............................................
    if (LINELP /= 17) return
    LRA = NLG1
    call D1V_Ar_ACT_Bl_DD_EXT_SGT0(LRA)

  case (25)
    !===================================================================
    ! VV(25) ACT -BL- ..................................................
    if ((LINELP /= 14) .or. (NLG2 /= 1)) return
    ! VV(25) Drl(33)-
    IWDL = 0
    IWDR = 0
    do MPL=1,MTYPE
      VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0_VV
      VPLP_W1(MPL) = Zero
    end do
    do MPL=1,MHLP
      IWAL = LPNEW_LWEI(MPL)
      IWAR = LPNEW_RWEI(MPL)
      LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
      LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
    end do
    do LRK=1,NORB_DZ
      call Drl_DD_EXT(LRK)
    end do

  case (7,9:10,14,16:18,26)
end select

return

end subroutine dd_ext_head_in_dbl

subroutine dd_ext_head_in_act()

use gugaci_global, only: linelp, logic_dh, nlg1, nlg2
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: lrai, lraj

LOGIC_DH = .false.
LRAI = NLG1
LRAJ = NLG2
! line=5 A&r-B&l<-->EXT
if (LINELP == 5) then
  call Ar_BL_DD_EXT(lrai,lraj,1)
end if
! line=9 D&rl<-->EXT
if (LINELP == 9) then
  call Drl_DD_EXT(lrai)
end if

return

end subroutine dd_ext_head_in_act

subroutine ss_drt_ci_new()

use gugaci_global, only: idisk_array, idisk_lp, idownwei_g131415, iml, imr, ipae, ipael, iseg_downwei, jml, jmr, jpad, jpadl, &
                         jpadlr, linelp, lpblock_ss, nvalue_space_ss
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: imlr, jptyl, jptyr, lpb

call external_space_plpmode_value_ss()
idisk_lp = idisk_array(13)

do lpb=1,lpblock_ss
  call read_lp()
  IpaeL = iml+17
  Ipae = imr+17
  imlr = Mul(iml,imr)
  nvalue_space_ss = iseg_downwei(9+imlr)
  idownwei_g131415 = iseg_downwei(17+iml)
  !if (imlr == 1) imspace= i ml
  call logicg_st(iml,imr,4,4)      ! irtype=4(S),lptype=5:ArBl-
  call get_jpty(jpadlr,jptyl,jptyr)
  call get_jp(jptyl,jml,jpadl,1)
  call get_jp(jptyr,jmr,jpad,1)
  !JMLR = Mul(JML,JMR)
  if (linelp <= 12) then
    call ss_ext_head_in_act()
  else
    call ss_ext_head_in_dbl()
  end if
end do

return

end subroutine ss_drt_ci_new

subroutine ss_ext_head_in_dbl()

use gugaci_global, only: ipae, ipael, jb_sys, jml, jmr, jpad, jpadl, jpadlr, jud, just, linelp, logic_dh, lp_lwei, lp_rwei, &
                         lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, nlg1, nlg2, norb_dz, norb_frz, vplp_w0, vplp_w1, &
                         vplpnew_w0, vplpnew_w1, w0_dd, w0_dv, w0_ss, w0_tt, w0_vv, w1_dd, w1_ss, w1_st, w1_tt
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, jmlr, lmi, lmij, lmj, lpok, lr0, lra, lri, lrj, mpl, ni, nk
real(kind=wp) :: w0dd1, w0dd2, w0dd3, w0dv1, w0ss15, w0ss17, w0ss20, w0tt2, w0tt3, w1dd1, w1dd2, w1ss15, w1ss17, w1tt2
integer(kind=iwp), external :: iwalk_ad

LOGIC_DH = .true.
JMLR = Mul(JML,JMR)
LPOK = JPADLR
select case (LPOK)
  case (1)
    !===================================================================
    ! SS(1)   ACT -C"-
    !-------------------------------------------------------------------
    ! SS(1-1)  Ar(01)-Bl(32)-        ACT -C"-
    ! SS(1-3)  Ar(13)-Bl(20)-        ACT -C"-
    ! SS(1-6)  (11)-Ar(23)-Bl(32)-   ACT -C"-
    ! SS(1-7)  Ar(13)-C'(21)-Bl(32)- ACT -C"-
    ! SS(1-8)  Ar(13)-C'(22)-Bl(31)- ACT -C"-
    ! SS(1-9)  Ar(23)-C'(11)-Bl(32)- ACT -C"-
    ! SS(1-11) Ar(13)-Bl(31)-C"(22)- ACT -C"-
    ! SS(1-12) Ar(13)-Bl(32)-C"(21)- ACT -C"-
    ! SS(1-13) Ar(23)-Bl(31)-C"(12)- ACT -C"-
    ! SS(1-16) (11)-Drl(22)-         ACT -C"-
    ! SS(1-18) Drl(11)-C"(22)-       ACT -C"-
    ! SS(1-19) Drl(12)-C"(21)-       ACT -C"-
    ! SS(1-20) (11)-Drl(33)-C"(22)-  ACT -C"-
    ! SS(1-20) Drl(33)-C"(11)-C"(22)-ACT -C"-
    !-------------------------------------------------------------------
    ! SS(1)   ACT 14: -C"-
    if (LINELP /= 14) return
    if (NLG2 == 1) then
      if (JB_SYS > 0) call SS_Drl_ACT_C_EXT_AB_SGT0(1)
      W0SS15 = W0_SS(15)
      W1SS15 = W1_SS(15)
      W0SS17 = W0_SS(17)
      W1SS17 = W1_SS(17)
      W0SS20 = W0_SS(20)
      if ((JML == 1) .and. (JMR == 1)) then
        ! SS(1-20) Drl(33)-C"(00)-       ACT -C"-                  ! IPL(R)AD=
        do LR0=NORB_FRZ+1,NORB_DZ
          IWDL = JUST(LR0,LR0)
          IWDR = IWDL
          do MPL=1,MHLP
            IWAL = LPNEW_LWEI(MPL)
            IWAR = LPNEW_RWEI(MPL)
            LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          end do
          do MPL=1,MTYPE
            VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS20
            VPLP_W1(MPL) = Zero
          end do
          call Drl_SS_SUM(LR0,0)
        end do
      end if
      do LRI=NORB_FRZ+1,NORB_DZ-1
        LMI = LSM_INN(LRI)
        do LRJ=LRI+1,NORB_DZ
          LMJ = LSM_INN(LRJ)
          LMIJ = Mul(LMI,LMJ)
          if ((LMIJ /= JML) .or. (LMIJ /= JMR)) cycle
          IWDL = JUST(LRI,LRJ)
          IWDR = IWDL
          do MPL=1,MHLP
            IWAL = LPNEW_LWEI(MPL)
            IWAR = LPNEW_RWEI(MPL)
            LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          end do
          ! SS(1-15) (22)-Drl(11)-         ACT -C"-
          do MPL=1,MTYPE
            VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS15
            VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SS15
          end do
          call Drl_SS_EXT(LRJ)
          ! SS(1-17) Drl(22)-C"(11)-       ACT -C"-
          do MPL=1,MTYPE
            VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS17
            VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SS17
          end do
          call Drl_SS_EXT(LRI)
          ! SS(1-20) (22)(11)Drl(33)-      ACT -C"-
          ! SS(1-20) (22)Drl(33)-C"(11)-   ACT -C"-
          ! SS(1-20) Drl(33)-C"(22)-C"(11)-ACT -C"-
          do MPL=1,MTYPE
            VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS20
            VPLP_W1(MPL) = Zero
          end do
          call Drl_SS_SUM(LRI,LRJ)
        end do
      end do
    else
      if (JB_SYS > 0) then
        call SS_ArBl_ACT_C_EXT_AB_SGT0(1)
        call SS_S_Drl_ACT_C_EXT_AB_SGT0(1)
      end if
      do LRI=NORB_FRZ+1,NORB_DZ-1
        do LRJ=LRI+1,NORB_DZ
          call SS2_EXT(LRI,LRJ,NK)
          if (NK /= 0) call Ar_BL_EXT_SS(LRI,LRJ,1)
          call SS4_EXT(LRI,LRJ,NK)
          if (NK /= 0) call Ar_BL_EXT_SS(LRI,LRJ,1)
          call SS5_EXT(LRI,LRJ,NK)
          if (NK /= 0) call Ar_BL_EXT_SS(LRI,LRJ,NK)
          call SS10_EXT(LRI,LRJ,NK)
          if (NK /= 0) call Ar_BL_EXT_SS(LRI,LRJ,NK)
          call SS14_EXT(LRI,LRJ,NK)
          if (NK /= 0) call Ar_BL_EXT_SS(LRI,LRJ,NK)
        end do
      end do
    end if

  case (2)
    ! ST(2)   ACT -C"-
    if ((LINELP /= 14) .or. (NLG2 == 1)) return
    if (JB_SYS > 0) call ST_Drl_ACT_C_EXT_AB_SGT0(1)
    if (JB_SYS > 0) call ST_ArBl_ACT_C_EXT_AB_SGT0(1)
    do LRI=NORB_FRZ+1,NORB_DZ-1
      LMI = LSM_INN(LRI)
      do LRJ=LRI+1,NORB_DZ
        LMJ = LSM_INN(LRJ)
        LMIJ = Mul(LMI,LMJ)
        if ((LMIJ /= JML) .or. (LMIJ /= JMR)) cycle
        ! ST(2-5) (22)Drl(12)-          ACT -C"-
        IWDL = JUST(LRI,LRJ)
        IWDR = IWDL
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = Zero
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_ST(5)
        end do
        call Drl_SS_EXT(LRJ)
        ! ST(2-6) Drl(22)-C"(12)-       ACT -C"-
        do MPL=1,MTYPE
          VPLP_W0(MPL) = Zero
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_ST(6)
        end do
        call Drl_SS_EXT(LRI)
        !---------------------------------------------------------------
        ! ST(2-3) Ar(13)-C'(22)-Bl(32)-   ACT -C"-
        ! ST(2-3) Ar(13)-Bl(32)-C'(22)-   ACT -C"-
        ! ST(2-7) Drl(12)-C"(22)-         ACT -C"-
        !---------------------------------------------------------------
      end do
    end do

  case (4)
    !===================================================================
    ! TS(3)   ACT -C"-  no used
    !===================================================================
    ! ST1(4) Ar-Bl- Drl- ACT -C"-
    if ((LINELP /= 14) .or. (nlg2 /= 2)) return
    call STT_ArBl_ACT_C_EXT_AB_SGT1(1)

  case (5)
    !===================================================================
    ! T1S(5) Ar-Bl- Drl ACT -C"-
    if ((LINELP /= 14) .or. (NLG2 /= 2)) return
    call TTS_Drl_ACT_C_EXT_AB_SGT1(1)
    call TTS_ArBl_ACT_C_EXT_AB_SGT1(1)

  case default ! (6)
    !===================================================================
    ! SD(6-1) ACT -B&L-
    if (LINELP /= 17) return
    LRA = NLG1
    call SD_AR_ACT_BL(1,LRA)
    if (JB_SYS > 0) call SD_AR_ACT_BL_SGT0(1,LRA)

  case (8)
    !===================================================================
    ! SD1(8) Ar ACT -Bl-
    if (LINELP /= 17) return
    LRA = NLG1
    if (JB_SYS > 0) call SDD_AR_ACT_BL_SGT0(1,LRA)

  case (11)
    !===================================================================
    ! TT(11) AR-BL   ACT -C"-
    if (LINELP /= 14) return
    if (NLG2 == 1) then
      W0TT2 = W0_TT(2)
      W1TT2 = W1_TT(2)
      W0TT3 = W0_TT(3)
      do LRI=NORB_FRZ+1,NORB_DZ-1
        LMI = LSM_INN(LRI)
        do LRJ=LRI+1,NORB_DZ
          LMJ = LSM_INN(LRJ)
          LMIJ = Mul(LMI,LMJ)
          if ((LMIJ /= JML) .or. (LMIJ /= JMR)) cycle
          ! TT(11-2) (22)Drl(22)-
          ! TT(11-2) Drl(22)-C"(22)-
          do MPL=1,MTYPE
            VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0TT2
            VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1TT2
          end do
          IWDL = JUST(LRI,LRJ)
          IWDR = IWDL
          do MPL=1,MHLP
            IWAL = LPNEW_LWEI(MPL)
            IWAR = LPNEW_RWEI(MPL)
            LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          end do
          call Drl_SS_EXT(LRI)
          call Drl_SS_EXT(LRJ)
          do MPL=1,MTYPE
            VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0TT3
            VPLP_W1(MPL) = Zero
          end do
          ! TT(11-3) Drl(33)-C"(22)-C"(22)-
          ! TT(11-3) (22)Drl(33)-C"(22)-
          ! TT(11-3) (22)(22)Drl(33)-
          call Drl_SS_SUM(LRI,LRJ)
        end do
      end do
    else
      do LRI=NORB_FRZ+1,NORB_DZ-1
        do LRJ=LRI+1,NORB_DZ
          call TT1_EXT(LRI,LRJ,NK,1)
          if (NK /= 0) call AR_BL_EXT_SS(LRI,LRJ,NK)
          call TT1_EXT(LRI,LRJ,NK,-1)
          if (NK /= 0) call AR_BL_EXT_SS(LRI,LRJ,NK)
        end do
      end do
    end if

  case (12)
    !===================================================================
    ! TTTT(12) Ar-Bl- Drl- ACT -C"-
    if (LINELP /= 14) return
    if (NLG2 == 1) then
      call TTTT_Drl_ACT_C_EXT_AB_SGT0(1)
    else
      call TTTT_ArBl_ACT_C_EXT_AB_SGT0(1)
    end if

  case (13)
    !===================================================================
    ! TD(13) ACT -BL-
    if (LINELP /= 17) return
    LRA = NLG1
    call TD_AR_ACT_BL(1,LRA)

  case (15)
    !===================================================================
    ! T1D1(15) Ar- ACT -Bl-
    if (LINELP /= 17) return
    LRA = NLG1
    call TTDD_AR_ACT_BL_SGT1(1,LRA)

  case (19)
    !===================================================================
    ! DD(19) ACT -C"- ..................................................
    if (LINELP /= 14) return
    if (NLG2 == 1) then
      W0DD2 = W0_DD(2)
      W1DD2 = W1_DD(2)
      W0DD3 = W0_DD(3)
      do LRI=NORB_FRZ+1,NORB_DZ
        LMI = LSM_INN(LRI)
        if (LMI /= JML) cycle
        ! DD(19-2) Drl(22)-
        do MPL=1,MTYPE
          VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0DD2
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1DD2
        end do
        IWDL = JUD(LRI)
        IWDR = IWDL
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        end do
        call Drl_SS_EXT(LRI)
        do MPL=1,MTYPE
          VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0DD3
          VPLP_W1(MPL) = Zero
        end do
        ! DD(19-3) (22)Drl(33)-
        ! DD(19-3) Drl(33)-C"(22)-
        call Drl_SS_SUM(LRI,0)
      end do
    else
      do LRI=NORB_FRZ+1,NORB_DZ-1
        LMI = LSM_INN(LRI)
        do LRJ=LRI+1,NORB_DZ
          LMJ = LSM_INN(LRJ)
          LMIJ = Mul(LMI,LMJ)
          if (LMIJ /= JMLR) cycle
          W0DD1 = W0_DD(1)
          W1DD1 = W1_DD(1)
          NI = mod(LRJ-LRI,2)
          if (NI == 0) then
            W0DD1 = -W0DD1
            W1DD1 = -W1DD1
          end if
          if ((LMI == JML) .and. (LMJ == JMR)) then
            ! DD(19-1) Ar(23)-Bl(32)-      ACT -C"-
            do MPL=1,MTYPE
              VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0DD1
            end do
            IWDL = JUD(LRI)
            IWDR = JUD(LRJ)
            do MPL=1,MHLP
              IWAL = LPNEW_LWEI(MPL)
              IWAR = LPNEW_RWEI(MPL)
              LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            end do
            call Ar_BL_EXT_SS(LRI,LRJ,1)
          end if
        end do
      end do
    end if

  case (20)
    !===================================================================
    ! D1D1(20) Drl- Ar-Bl- ACT -C"-
    if (LINELP /= 14) return
    if (NLG2 == 1) then
      call D1D1_Drl_ACT_C_EXT_AB_SGT0(1)
    else
      call D1D1_ArBl_ACT_C_EXT_AB_SGT0(1)
    end if

  case (21)
    !===================================================================
    ! DD1(21) Ar-Bl- ACT -C"-
    if ((LINELP /= 14) .or. (nlg2 /= 2)) return
    call DD1_ArBl_ACT_C_EXT_AB_SGT0(1)

  case (22)
    !===================================================================
    ! D1D(22) Ar-Bl- Drl- ACT -C"-
    if ((LINELP /= 14) .or. (nlg2 /= 2)) return
    call D1D_ArBl_ACT_C_EXT_AB_SGT0(1)
    call D1D_Drl_ACT_C_EXT_AB_SGT0(1)

  case (23)
    !===================================================================
    ! DV(23) ACT -C'-..................................................
    if (LINELP /= 17) return
    LRA = NLG1
    do LRI=NORB_FRZ+1,NORB_DZ
      LMI = LSM_INN(LRI)
      if (LMI /= JMLR) cycle
      W0DV1 = W0_DV(1)
      NI = mod(NORB_DZ-LRI,2)
      if (NI == 1) W0DV1 = -W0DV1
      ! DV(23-1) A&r(23)-
      IWDL = JUD(LRI)
      IWDR = 0
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      do MPL=1,MTYPE
        VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0DV1
      end do
      call Ar_BL_EXT_SS(LRI,LRA,1)
    end do

  case (24)
    !===================================================================
    ! D1V(24) Ar- ACT -Bl-
    if (LINELP /= 17) return
    LRA = NLG1
    call D1V_Ar_ACT_Bl_EXT_AB_SGT0(1,LRA)

  case (25)
    !===================================================================
    ! VV(25) ACT -C"- .false.
    if ((LINELP /= 14) .or. (NLG2 /= 1)) return
    !VV(25) Drl(33)-
    IWDL = 0
    IWDR = 0
    do MPL=1,MTYPE
      VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0_VV
      VPLP_W1(MPL) = Zero
    end do
    do MPL=1,MHLP
      IWAL = LPNEW_LWEI(MPL)
      IWAR = LPNEW_RWEI(MPL)
      LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
      LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
    end do
    call Drl_SS_SUM(0,0)
    !do lrk=1,norb_dz
    !  call Drl_ss_ext(lrk)
    !end do

  case (3,7,9:10,14,16:18,26)
end select

return

end subroutine ss_ext_head_in_dbl

subroutine ss_ext_head_in_act()

use gugaci_global, only: linelp, logic_dh, nlg1, nlg2
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: lrai, lraj

LOGIC_DH = .false.
LRAI = NLG1
LRAJ = NLG2
! line=5 A&r-B&l<-->EXT
if (LINELP == 5) then
  call Ar_BL_EXT_SS(lrai,lraj,1)
end if
! line=9 D&rl<-->EXT
if (LINELP == 9) then
  call Drl_SS_EXT(lrai)
end if

return

end subroutine ss_ext_head_in_act

subroutine st_drt_ci_new()

use gugaci_global, only: idisk_array, idisk_lp, idownwei_g131415, iml, imr, ipae, ipael, iseg_downwei, jml, jmr, jpad, jpadl, &
                         jpadlr, linelp, lpblock_st, nvalue_space_ss
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: imlr, jptyl, jptyr, lpb

call external_space_plpmode_value_ST()
idisk_lp = idisk_array(12)

do lpb=1,lpblock_st
  call read_lp()
  IpaeL = iml+17
  Ipae = imr+9
  imlr = Mul(iml,imr)
  idownwei_g131415 = iseg_downwei(9+iml)  !(17+iml)???
  nvalue_space_ss = iseg_downwei(9+imlr)
  call logicg_st(iml,imr,4,3)             ! irtype=4(S),3(T)
  call get_jpty(jpadlr,jptyl,jptyr)
  call get_jp(jptyl,jml,jpadl,1)
  call get_jp(jptyr,jmr,jpad,1)
  !JMLR = Mul(JML,JMR)
  if (linelp <= 12) then
    call st_ext_head_in_act()
  else
    call st_ext_head_in_dbl()
  end if
end do

return

end subroutine st_drt_ci_new

subroutine st_ext_head_in_dbl()

use gugaci_global, only: ipae, ipael, jb_sys, jml, jmr, jpad, jpadl, jpadlr, jud, just, linelp, logic_dh, lp_lwei, lp_rwei, &
                         lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, nlg1, nlg2, norb_dz, norb_frz, vplp_w0, vplp_w1, &
                         vplpnew_w0, vplpnew_w1, w0_dd, w0_dv, w0_ss, w0_tt, w1_dd, w1_ss, w1_st, w1_tt
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, jmlr, lmi, lmij, lmj, lpok, lra, lri, lrj, lrk, mpl, ni, nk
real(kind=wp) :: w0dd1, w0dd2, w0dv1, w0ss15, w0ss17, w0ss20, w0tt2, w1dd1, w1dd2, w1ss15, w1ss17, w1tt2
integer(kind=iwp), external :: iwalk_ad

LOGIC_DH = .true.
JMLR = Mul(JML,JMR)
LPOK = JPADLR
select case (LPOK)
  case (1)
    !===================================================================
    ! SS(1)   ACT -C"-
    !-------------------------------------------------------------------
    ! SS(1-1)  Ar(01)-Bl(32)-        ACT -C"-
    ! SS(1-3)  Ar(13)-Bl(20)-        ACT -C"-
    ! SS(1-6)  (11)-Ar(23)-Bl(32)-   ACT -C"-
    ! SS(1-7)  Ar(13)-C'(21)-Bl(32)- ACT -C"-
    ! SS(1-8)  Ar(13)-C'(22)-Bl(31)- ACT -C"-
    ! SS(1-9)  Ar(23)-C'(11)-Bl(32)- ACT -C"-
    ! SS(1-11) Ar(13)-Bl(31)-C"(22)- ACT -C"-
    ! SS(1-12) Ar(13)-Bl(32)-C"(21)- ACT -C"-
    ! SS(1-13) Ar(23)-Bl(31)-C"(12)- ACT -C"-
    ! SS(1-16) (11)-Drl(22)-         ACT -C"-
    ! SS(1-18) Drl(11)-C"(22)-       ACT -C"-
    ! SS(1-19) Drl(12)-C"(21)-       ACT -C"-
    ! SS(1-20) (11)-Drl(33)-C"(22)-  ACT -C"-
    ! SS(1-20) Drl(33)-C"(11)-C"(22)-ACT -C"-
    !-------------------------------------------------------------------
    ! SS(1)   ACT 14: -C"-
    if (LINELP /= 14) return
    if (NLG2 == 1) then
      if (JB_SYS > 0) call SS_Drl_ACT_C_EXT_AB_SGT0(2)
      W0SS15 = W0_SS(15)
      W1SS15 = W1_SS(15)
      W0SS17 = W0_SS(17)
      W1SS17 = W1_SS(17)
      W0SS20 = W0_SS(20)
      do LRI=NORB_FRZ+1,NORB_DZ-1
        LMI = LSM_INN(LRI)
        do LRJ=LRI+1,NORB_DZ
          LMJ = LSM_INN(LRJ)
          LMIJ = Mul(LMI,LMJ)
          if ((LMIJ /= JML) .or. (LMIJ /= JMR)) cycle
          IWDL = JUST(LRI,LRJ)
          IWDR = IWDL
          do MPL=1,MHLP
            IWAL = LPNEW_LWEI(MPL)
            IWAR = LPNEW_RWEI(MPL)
            LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          end do
          ! SS(1-15) (22)-Drl(11)-         ACT -C"-
          do MPL=1,MTYPE
            VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS15
            VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SS15
          end do
          call Drl_ST_EXT(LRJ)
          ! SS(1-17) Drl(22)-C"(11)-       ACT -C"-
          do MPL=1,MTYPE
            VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS17
            VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SS17
          end do
          call Drl_ST_EXT(LRI)
          ! SS(1-20) (22)(11)Drl(33)-      ACT -C"-
          ! SS(1-20) (22)Drl(33)-C"(11)-   ACT -C"-
          ! SS(1-20) Drl(33)-C"(22)-C"(11)-ACT -C"-
          do MPL=1,MTYPE
            VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS20
            VPLP_W1(MPL) = Zero
          end do
          do LRK=1,NORB_DZ
            if (LRK == LRI) cycle
            if (LRK == LRJ) cycle
            call Drl_ST_EXT(LRK)
          end do
        end do
      end do
    else
      LRA = NLG1
      if (JB_SYS > 0) then
        call SS_ArBl_ACT_C_EXT_AB_SGT0(2)
        call SS_S_Drl_ACT_C_EXT_AB_SGT0(2)
      end if
      do LRI=NORB_FRZ+1,NORB_DZ-1
        do LRJ=LRI+1,NORB_DZ
          call SS2_EXT(LRI,LRJ,NK)
          if (NK /= 0) call Ar_BL_EXT_ST(LRI,LRJ,1)
          call SS4_EXT(LRI,LRJ,NK)
          if (NK /= 0) call Ar_BL_EXT_ST(LRI,LRJ,1)
          call SS5_EXT(LRI,LRJ,NK)
          if (NK /= 0) call Ar_BL_EXT_ST(LRI,LRJ,NK)
          call SS10_EXT(LRI,LRJ,NK)
          if (NK /= 0) call Ar_BL_EXT_ST(LRI,LRJ,NK)
          call SS14_EXT(LRI,LRJ,NK)
          if (NK /= 0) call Ar_BL_EXT_ST(LRI,LRJ,NK)
        end do
      end do
    end if

  case (2)
    ! ST(2)   ACT -C"-
    if ((LINELP /= 14) .or. (NLG2 == 1)) return
    if (JB_SYS > 0) call ST_Drl_ACT_C_EXT_AB_SGT0(2)
    if (JB_SYS > 0) call ST_ArBl_ACT_C_EXT_AB_SGT0(2)
    do LRI=NORB_FRZ+1,NORB_DZ-1
      do LRJ=LRI+1,NORB_DZ
        call ST1_EXT(LRI,LRJ,NK)
        if (NK /= 0) call Ar_BL_EXT_ST(LRI,LRJ,1)
        call ST2_EXT(LRI,LRJ,NK)
        if (NK /= 0) call Ar_BL_EXT_ST(LRI,LRJ,NK)
        call ST4_EXT(LRI,LRJ,NK,1)
        if (NK /= 0) call Ar_BL_EXT_ST(LRI,LRJ,NK)
        call ST4_EXT(LRI,LRJ,NK,-1)
        if (NK /= 0) call Ar_BL_EXT_ST(LRI,LRJ,NK)
      end do
    end do

    do LRI=NORB_FRZ+1,NORB_DZ-1
      LMI = LSM_INN(LRI)
      do LRJ=LRI+1,NORB_DZ
        LMJ = LSM_INN(LRJ)
        LMIJ = Mul(LMI,LMJ)
        if ((LMIJ /= JML) .or. (LMIJ /= JMR)) cycle
        ! ST(2-5) (22)Drl(12)-          ACT -C"-
        IWDL = JUST(LRI,LRJ)
        IWDR = IWDL
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = Zero
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_ST(5)
        end do
        call Drl_ST_EXT(LRJ)
        ! ST(2-6) Drl(22)-C"(12)-       ACT -C"-
        do MPL=1,MTYPE
          VPLP_W0(MPL) = Zero
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_ST(6)
        end do
        call Drl_ST_EXT(LRI)
        !---------------------------------------------------------------
        ! ST(2-3) Ar(13)-C'(22)-Bl(32)-   ACT -C"-
        ! ST(2-3) Ar(13)-Bl(32)-C'(22)-   ACT -C"-
        ! ST(2-7) Drl(12)-C"(22)-         ACT -C"-
        !---------------------------------------------------------------
      end do
    end do

  case (3)
    !===================================================================
    ! TS(3) A&R-B^L-  ACT -C"-
    if ((LINELP /= 14) .or. (NLG2 == 1)) return
    if (JB_SYS > 0) then
      call TS_ArBl_ACT_C_EXT_AB_SGT0(2)
    end if
    do LRI=NORB_FRZ+1,NORB_DZ-1
      do LRJ=LRI+1,NORB_DZ
        call TS1_EXT(LRI,LRJ,NK)
        if (NK /= 0) call Ar_BL_EXT_ST(LRI,LRJ,1)
        call TS2_EXT(LRI,LRJ,NK,1)
        if (NK /= 0) call Ar_BL_EXT_ST(LRI,LRJ,NK)
        call TS2_EXT(LRI,LRJ,NK,-1)
        if (NK /= 0) call Ar_BL_EXT_ST(LRI,LRJ,NK)
        call TS4_EXT(LRI,LRJ,NK)
        if (NK /= 0) call Ar_BL_EXT_ST(LRI,LRJ,NK)
      end do
    end do

  case (4)
    !===================================================================
    ! ST1(4) Ar-Bl- Drl- ACT -C"-
    if ((LINELP /= 14) .or. (nlg2 /= 2)) return
    call STT_ArBl_ACT_C_EXT_AB_SGT1(2)

  case (5)
    !===================================================================
    ! T1S(5) Ar-Bl- Drl ACT -C"-
    if ((LINELP /= 14) .or. (NLG2 /= 2)) return
    call TTS_Drl_ACT_C_EXT_AB_SGT1(2)
    call TTS_ArBl_ACT_C_EXT_AB_SGT1(2)

  case default ! (6)
    !===================================================================
    ! SD(6-1) ACT -B&L-
    if (LINELP /= 17) return
    LRA = NLG1
    call SD_AR_ACT_BL(2,LRA)
    if (JB_SYS > 0) call SD_AR_ACT_BL_SGT0(2,LRA)

  case (8)
    !===================================================================
    ! SD1(8) Ar ACT -Bl-
    if (LINELP /= 17) return
    LRA = NLG1
    if (JB_SYS > 0) call SDD_AR_ACT_BL_SGT0(2,LRA)

  case (11)
    !===================================================================
    ! TT(11) Drl-  ACT -C"-
    if (LINELP /= 14) return
    if (NLG2 == 1) then
      W0TT2 = W0_TT(2)
      W1TT2 = W1_TT(2)
      do LRI=NORB_FRZ+1,NORB_DZ-1
        LMI = LSM_INN(LRI)
        do LRJ=LRI+1,NORB_DZ
          LMJ = LSM_INN(LRJ)
          LMIJ = Mul(LMI,LMJ)
          if ((LMIJ /= JML) .or. (LMIJ /= JMR)) cycle
          ! TT(11-2) (22)Drl(22)-
          ! TT(11-2) Drl(22)-C"(22)-
          do MPL=1,MTYPE
            VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0TT2
            VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1TT2
          end do
          IWDL = JUST(LRI,LRJ)
          IWDR = IWDL
          do MPL=1,MHLP
            IWAL = LPNEW_LWEI(MPL)
            IWAR = LPNEW_RWEI(MPL)
            LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          end do
          call Drl_ST_EXT(LRI)
          call Drl_ST_EXT(LRJ)
        end do
      end do
    else
      do LRI=NORB_FRZ+1,NORB_DZ-1
        do LRJ=LRI+1,NORB_DZ
          call TT1_EXT(LRI,LRJ,NK,1)
          call AR_BL_EXT_TS(LRI,LRJ,NK)
          call TT1_EXT(LRI,LRJ,NK,-1)
          call AR_BL_EXT_ST(LRI,LRJ,NK)
        end do
      end do
    end if

  case (12)
    !===================================================================
    ! TTTT(12) Ar-Bl- Drl- ACT -C"-
    if (LINELP /= 14) return
    if (NLG2 == 1) then
      call TTTT_Drl_ACT_C_EXT_AB_SGT0(2)
    else
      call TTTT_ArBl_ACT_C_EXT_AB_SGT0(2)
    end if

  case (13)
    !===================================================================
    ! TD(13) ACT -BL-
    if (LINELP /= 17) return
    LRA = NLG1
    call TD_AR_ACT_BL(2,LRA)

  case (15)
    !===================================================================
    ! T1D1(15) Ar- ACT -Bl-
    if (LINELP /= 17) return
    LRA = NLG1
    call TTDD_AR_ACT_BL_SGT1(2,LRA)

  case (19)
    !===================================================================
    ! DD(19) ACT -C"- ..................................................
    if (LINELP /= 14) return
    if (NLG2 == 1) then
      W0DD2 = W0_DD(2)
      W1DD2 = W1_DD(2)
      do LRI=NORB_FRZ+1,NORB_DZ
        LMI = LSM_INN(LRI)
        if (LMI /= JML) cycle
        ! DD(19-2) Drl(22)-
        do MPL=1,MTYPE
          VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0DD2
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1DD2
        end do
        IWDL = JUD(LRI)
        IWDR = IWDL
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        end do
        call Drl_ST_EXT(LRI)
      end do
    else
      do LRI=NORB_FRZ+1,NORB_DZ-1
        LMI = LSM_INN(LRI)
        do LRJ=LRI+1,NORB_DZ
          LMJ = LSM_INN(LRJ)
          LMIJ = Mul(LMI,LMJ)
          if (LMIJ /= JMLR) cycle
          W0DD1 = W0_DD(1)
          W1DD1 = W1_DD(1)
          NI = mod(LRJ-LRI,2)
          if (NI == 0) then
            W0DD1 = -W0DD1
            W1DD1 = -W1DD1
          end if
          if ((LMI == JML) .and. (LMJ == JMR)) then
            ! DD(19-1) Ar(23)-Bl(32)-      ACT -C"-
            do MPL=1,MTYPE
              VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1DD1
            end do
            IWDL = JUD(LRI)
            IWDR = JUD(LRJ)
            do MPL=1,MHLP
              IWAL = LPNEW_LWEI(MPL)
              IWAR = LPNEW_RWEI(MPL)
              LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            end do
            call Ar_BL_EXT_ST(LRI,LRJ,1)
          end if
        end do
      end do
    end if

  case (20)
    !===================================================================
    ! D1D1(20) Drl- Ar-Bl- ACT -C"-
    if (LINELP /= 14) return
    if (NLG2 == 1) then
      call D1D1_Drl_ACT_C_EXT_AB_SGT0(2)
    else
      call D1D1_ArBl_ACT_C_EXT_AB_SGT0(2)
    end if

  case (21)
    !===================================================================
    ! DD1(21) Ar-Bl- ACT -C"-
    if ((LINELP /= 14) .or. (nlg2 /= 2)) return
    call DD1_ArBl_ACT_C_EXT_AB_SGT0(2)

  case (22)
    !===================================================================
    ! D1D(22) Ar-Bl- Drl- ACT -C"-
    if ((LINELP /= 14) .or. (nlg2 /= 2)) return
    call D1D_ArBl_ACT_C_EXT_AB_SGT0(2)
    call D1D_Drl_ACT_C_EXT_AB_SGT0(2)

  case (23)
    !===================================================================
    ! DV(23) ACT -C'-..................................................
    if (LINELP /= 17) return
    LRA = NLG1
    do LRI=NORB_FRZ+1,NORB_DZ
      LMI = LSM_INN(LRI)
      if (LMI /= JMLR) cycle
      W0DV1 = W0_DV(1)
      NI = mod(NORB_DZ-LRI,2)
      if (NI == 1) W0DV1 = -W0DV1
      ! DV(23-1) A&r(23)-
      IWDL = JUD(LRI)
      IWDR = 0
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      do MPL=1,MTYPE
        VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W0DV1
      end do
      call Ar_BL_EXT_ST(LRI,LRA,1)
    end do

  case (24)
    !===================================================================
    ! D1V(24) Ar- ACT -Bl-
    if (LINELP /= 17) return
    LRA = NLG1
    call D1V_Ar_ACT_Bl_EXT_AB_SGT0(2,LRA)

  case (25)
    !===================================================================
    ! VV(25) ACT -BL- ..................................................
    return

  case (7,9:10,14,16:18,26)
end select
return

end subroutine st_ext_head_in_dbl

subroutine st_ext_head_in_act()

use gugaci_global, only: linelp, logic_dh, nlg1, nlg2
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: lrai, lraj

LOGIC_DH = .false.
LRAI = NLG1
LRAJ = NLG2
!line=5 A&r-B&l<-->EXT
if (LINELP == 5) then
  call Ar_BL_EXT_ST(lrai,lraj,1)
end if
!line=9 D&rl<-->EXT
if (LINELP == 9) then
  call Drl_ST_EXT(lrai)
end if

return

end subroutine st_ext_head_in_act

subroutine ts_drt_ci_new()

use gugaci_global, only: idisk_array, idisk_lp, idownwei_g131415, iml, imr, ipae, ipael, iseg_downwei, jml, jmr, jpad, jpadl, &
                         jpadlr, linelp, lpblock_ts, nvalue_space_ss
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: imlr, jptyl, jptyr, lpb

call external_space_plpmode_value_TS()
idisk_lp = idisk_array(9)
do lpb=1,lpblock_ts
  call read_lp()
  IpaeL = iml+9
  Ipae = imr+17
  imlr = Mul(iml,imr)
  idownwei_g131415 = iseg_downwei(9+iml)  !(17+iml)???
  nvalue_space_ss = iseg_downwei(9+imlr)
  call logicg_st(iml,imr,3,4)             ! irtype=4(S),3(T)
  call get_jpty(jpadlr,jptyl,jptyr)
  call get_jp(jptyl,jml,jpadl,1)
  call get_jp(jptyr,jmr,jpad,1)
  !JMLR = Mul(JML,JMR)
  if (linelp <= 12) then
    call ts_ext_head_in_act()
  else
    call ts_ext_head_in_dbl()
  end if
end do

return

end subroutine ts_drt_ci_new

subroutine ts_ext_head_in_dbl()

use gugaci_global, only: ipae, ipael, jb_sys, jml, jmr, jpad, jpadl, jpadlr, jud, just, linelp, logic_dh, lp_lwei, lp_rwei, &
                         lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, nlg1, nlg2, norb_dz, norb_frz, vplp_w0, vplp_w1, &
                         vplpnew_w0, vplpnew_w1, w0_dd, w0_dv, w0_ss, w0_tt, w1_dd, w1_ss, w1_st, w1_tt
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, jmlr, lmi, lmij, lmj, lpok, lra, lri, lrj, lrk, mpl, ni, nk
real(kind=wp) :: w0dd1, w0dd2, w0dv1, w0ss15, w0ss17, w0ss20, w0tt2, w1dd1, w1dd2, w1ss15, w1ss17, w1tt2
integer(kind=iwp), external :: iwalk_ad

LOGIC_DH = .true.
JMLR = Mul(JML,JMR)
LPOK = JPADLR
select case (LPOK)
  case (1)
    !===================================================================
    ! SS(1)   ACT -C"-
    !-------------------------------------------------------------------
    ! SS(1-1)  Ar(01)-Bl(32)-        ACT -C"-
    ! SS(1-3)  Ar(13)-Bl(20)-        ACT -C"-
    ! SS(1-6)  (11)-Ar(23)-Bl(32)-   ACT -C"-
    ! SS(1-7)  Ar(13)-C'(21)-Bl(32)- ACT -C"-
    ! SS(1-8)  Ar(13)-C'(22)-Bl(31)- ACT -C"-
    ! SS(1-9)  Ar(23)-C'(11)-Bl(32)- ACT -C"-
    ! SS(1-11) Ar(13)-Bl(31)-C"(22)- ACT -C"-
    ! SS(1-12) Ar(13)-Bl(32)-C"(21)- ACT -C"-
    ! SS(1-13) Ar(23)-Bl(31)-C"(12)- ACT -C"-
    ! SS(1-16) (11)-Drl(22)-         ACT -C"-
    ! SS(1-18) Drl(11)-C"(22)-       ACT -C"-
    ! SS(1-19) Drl(12)-C"(21)-       ACT -C"-
    ! SS(1-20) (11)-Drl(33)-C"(22)-  ACT -C"-
    ! SS(1-20) Drl(33)-C"(11)-C"(22)-ACT -C"-
    !-------------------------------------------------------------------
    ! SS(1)   ACT 14: -C"-
    if (LINELP /= 14) return
    if (NLG2 == 1) then
      if (JB_SYS > 0) call SS_Drl_ACT_C_EXT_AB_SGT0(3)
      W0SS15 = W0_SS(15)
      W1SS15 = W1_SS(15)
      W0SS17 = W0_SS(17)
      W1SS17 = W1_SS(17)
      W0SS20 = W0_SS(20)
      do LRI=NORB_FRZ+1,NORB_DZ-1
        LMI = LSM_INN(LRI)
        do LRJ=LRI+1,NORB_DZ
          LMJ = LSM_INN(LRJ)
          LMIJ = Mul(LMI,LMJ)
          if ((LMIJ /= JML) .or. (LMIJ /= JMR)) cycle
          IWDL = JUST(LRI,LRJ)
          IWDR = IWDL
          do MPL=1,MHLP
            IWAL = LPNEW_LWEI(MPL)
            IWAR = LPNEW_RWEI(MPL)
            LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          end do
          ! SS(1-15) (22)-Drl(11)-         ACT -C"-
          do MPL=1,MTYPE
            VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS15
            VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SS15
          end do
          call Drl_TS_EXT(LRJ)
          ! SS(1-17) Drl(22)-C"(11)-       ACT -C"-
          do MPL=1,MTYPE
            VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS17
            VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SS17
          end do
          call Drl_TS_EXT(LRI)
          ! SS(1-20) (22)(11)Drl(33)-      ACT -C"-
          ! SS(1-20) (22)Drl(33)-C"(11)-   ACT -C"-
          ! SS(1-20) Drl(33)-C"(22)-C"(11)-ACT -C"-
          do MPL=1,MTYPE
            VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS20
            VPLP_W1(MPL) = Zero
          end do
          do LRK=1,NORB_DZ
            if (LRK == LRI) cycle
            if (LRK == LRJ) cycle
            call Drl_TS_EXT(LRK)
          end do
        end do
      end do
    else
      if (JB_SYS > 0) then
        call SS_ArBl_ACT_C_EXT_AB_SGT0(3)
        call SS_S_Drl_ACT_C_EXT_AB_SGT0(3)
      end if
      do LRI=NORB_FRZ+1,NORB_DZ-1
        do LRJ=LRI+1,NORB_DZ
          call SS2_EXT(LRI,LRJ,NK)
          if (NK /= 0) call AR_BL_EXT_TS(LRI,LRJ,1)
          call SS4_EXT(LRI,LRJ,NK)
          if (NK /= 0) call AR_BL_EXT_TS(LRI,LRJ,1)
          call SS5_EXT(LRI,LRJ,NK)
          if (NK /= 0) call AR_BL_EXT_TS(LRI,LRJ,NK)
          call SS10_EXT(LRI,LRJ,NK)
          if (NK /= 0) call AR_BL_EXT_TS(LRI,LRJ,NK)
          call SS14_EXT(LRI,LRJ,NK)
          if (NK /= 0) call AR_BL_EXT_TS(LRI,LRJ,NK)
        end do
      end do
    end if

  case (2)
    ! ST(2)   ACT -C"-
    if ((LINELP /= 14) .or. (NLG2 == 1)) return
    if (JB_SYS > 0) call ST_Drl_ACT_C_EXT_AB_SGT0(3)
    if (JB_SYS > 0) call ST_ArBl_ACT_C_EXT_AB_SGT0(3)
    do LRI=NORB_FRZ+1,NORB_DZ-1
      do LRJ=LRI+1,NORB_DZ
        call ST1_EXT(LRI,LRJ,NK)
        if (NK /= 0) call Ar_BL_EXT_TS(LRI,LRJ,1)
        call ST2_EXT(LRI,LRJ,NK)
        if (NK /= 0) call Ar_BL_EXT_TS(LRI,LRJ,NK)
        call ST4_EXT(LRI,LRJ,NK,1)
        if (NK /= 0) call Ar_BL_EXT_TS(LRI,LRJ,NK)
        call ST4_EXT(LRI,LRJ,NK,-1)
        if (NK /= 0) call Ar_BL_EXT_TS(LRI,LRJ,NK)
      end do
    end do

    do LRI=NORB_FRZ+1,NORB_DZ-1
      LMI = LSM_INN(LRI)
      do LRJ=LRI+1,NORB_DZ
        LMJ = LSM_INN(LRJ)
        LMIJ = Mul(LMI,LMJ)
        if ((LMIJ /= JML) .or. (LMIJ /= JMR)) cycle
        ! ST(2-5) (22)Drl(12)-          ACT -C"-
        IWDL = JUST(LRI,LRJ)
        IWDR = IWDL
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = Zero
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_ST(5)
        end do
        call Drl_TS_EXT(LRJ)
        ! ST(2-6) Drl(22)-C"(12)-       ACT -C"-
        do MPL=1,MTYPE
          VPLP_W0(MPL) = Zero
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_ST(6)
        end do
        call Drl_TS_EXT(LRI)
        !---------------------------------------------------------------
        ! ST(2-3) Ar(13)-C'(22)-Bl(32)-   ACT -C"-
        ! ST(2-3) Ar(13)-Bl(32)-C'(22)-   ACT -C"-
        ! ST(2-7) Drl(12)-C"(22)-         ACT -C"-
        !---------------------------------------------------------------
      end do
    end do

  case (3)
    !===================================================================
    ! TS(3) D&R^L-  ACT -C"-
    if ((LINELP /= 14) .or. (NLG2 == 1)) return
    if (JB_SYS > 0) then
      call TS_ArBl_ACT_C_EXT_AB_SGT0(3)
    end if
    do LRI=NORB_FRZ+1,NORB_DZ-1
      do LRJ=LRI+1,NORB_DZ
        call TS1_EXT(LRI,LRJ,NK)
        if (NK /= 0) call Ar_BL_EXT_TS(LRI,LRJ,1)
        call TS2_EXT(LRI,LRJ,NK,1)
        if (NK /= 0) call Ar_BL_EXT_TS(LRI,LRJ,NK)
        call TS2_EXT(LRI,LRJ,NK,-1)
        if (NK /= 0) call Ar_BL_EXT_TS(LRI,LRJ,NK)
        call TS4_EXT(LRI,LRJ,NK)
        if (NK /= 0) call Ar_BL_EXT_TS(LRI,LRJ,NK)
      end do
    end do

  case (4)
    !===================================================================
    ! ST1(4) Ar-Bl- Drl- ACT -C"-
    if ((LINELP /= 14) .or. (nlg2 /= 2)) return
    call STT_ArBl_ACT_C_EXT_AB_SGT1(3)

  case (5)
    !===================================================================
    ! T1S(5) Ar-Bl- Drl ACT -C"-
    if ((LINELP /= 14) .or. (NLG2 /= 2)) return
    call TTS_Drl_ACT_C_EXT_AB_SGT1(3)
    call TTS_ArBl_ACT_C_EXT_AB_SGT1(3)

  case default ! (6)
    !===================================================================
    ! SD(6-1) ACT -B&L-
    if (LINELP /= 17) return
    LRA = NLG1
    call SD_AR_ACT_BL(3,LRA)
    if (JB_SYS > 0) call SD_AR_ACT_BL_SGT0(3,LRA)

  case (8)
    !===================================================================
    ! SD1(8) Ar ACT -Bl-
    if (LINELP /= 17) return
    LRA = NLG1
    if (JB_SYS > 0) call SDD_AR_ACT_BL_SGT0(3,LRA)

  case (11)
    !===================================================================
    ! TT(11) Drl-  ACT -C"-
    if (LINELP /= 14) return
    if (NLG2 == 1) then
      W0TT2 = W0_TT(2)
      W1TT2 = W1_TT(2)
      do LRI=NORB_FRZ+1,NORB_DZ-1
        LMI = LSM_INN(LRI)
        do LRJ=LRI+1,NORB_DZ
          LMJ = LSM_INN(LRJ)
          LMIJ = Mul(LMI,LMJ)
          if ((LMIJ /= JML) .or. (LMIJ /= JMR)) cycle
          ! TT(11-2) (22)Drl(22)-
          ! TT(11-2) Drl(22)-C"(22)-
          do MPL=1,MTYPE
            VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0TT2
            VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1TT2
          end do
          IWDL = JUST(LRI,LRJ)
          IWDR = IWDL
          do MPL=1,MHLP
            IWAL = LPNEW_LWEI(MPL)
            IWAR = LPNEW_RWEI(MPL)
            LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          end do
          call Drl_TS_EXT(LRI)
          call Drl_TS_EXT(LRJ)
        end do
      end do
    else
      do LRI=NORB_FRZ+1,NORB_DZ-1
        do LRJ=LRI+1,NORB_DZ
          call TT1_EXT(LRI,LRJ,NK,1)
          call AR_BL_EXT_TS(LRI,LRJ,NK)
          call TT1_EXT(LRI,LRJ,NK,-1)
          call AR_BL_EXT_TS(LRI,LRJ,NK)
        end do
      end do
    end if

  case (12)
    !===================================================================
    ! TTTT(12) Ar-Bl- Drl- ACT -C"-
    if (LINELP /= 14) return
    if (NLG2 == 1) then
      call TTTT_Drl_ACT_C_EXT_AB_SGT0(3)
    else
      call TTTT_ArBl_ACT_C_EXT_AB_SGT0(3)
    end if

  case (13)
    !===================================================================
    ! TD(13) ACT -BL-
    if (LINELP /= 17) return
    LRA = NLG1
    call TD_AR_ACT_BL(3,LRA)

  case (15)
    !===================================================================
    ! T1D1(15) Ar- ACT -Bl-
    if (LINELP /= 17) return
    LRA = NLG1
    call TTDD_AR_ACT_BL_SGT1(3,LRA)

  case (19)
    !===================================================================
    ! DD(19) ACT -C"- ..................................................
    if (LINELP /= 14) return
    if (NLG2 == 1) then
      W0DD2 = W0_DD(2)
      W1DD2 = W1_DD(2)
      do LRI=NORB_FRZ+1,NORB_DZ
        LMI = LSM_INN(LRI)
        if (LMI /= JML) cycle
        ! DD(19-2) Drl(22)-
        do MPL=1,MTYPE
          VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0DD2
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1DD2
        end do
        IWDL = JUD(LRI)
        IWDR = IWDL
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        end do
        call Drl_TS_EXT(LRI)
      end do
    else
      do LRI=NORB_FRZ+1,NORB_DZ-1
        LMI = LSM_INN(LRI)
        do LRJ=LRI+1,NORB_DZ
          LMJ = LSM_INN(LRJ)
          LMIJ = Mul(LMI,LMJ)
          if (LMIJ /= JMLR) cycle
          W0DD1 = W0_DD(1)
          W1DD1 = W1_DD(1)
          NI = mod(LRJ-LRI,2)
          if (NI == 0) then
            W0DD1 = -W0DD1
            W1DD1 = -W1DD1
          end if
          if ((LMI == JML) .and. (LMJ == JMR)) then
            ! DD(19-1) Ar(23)-Bl(32)-      ACT -C"-
            do MPL=1,MTYPE
              VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1DD1
            end do
            IWDL = JUD(LRI)
            IWDR = JUD(LRJ)
            do MPL=1,MHLP
              IWAL = LPNEW_LWEI(MPL)
              IWAR = LPNEW_RWEI(MPL)
              LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            end do
            call AR_BL_EXT_TS(LRI,LRJ,1)
          end if
        end do
      end do
    end if

  case (20)
    !===================================================================
    ! D1D1(20) Drl- Ar-Bl- ACT -C"-
    if (LINELP /= 14) return
    if (NLG2 == 1) then
      call D1D1_Drl_ACT_C_EXT_AB_SGT0(3)
    else
      call D1D1_ArBl_ACT_C_EXT_AB_SGT0(3)
    end if

  case (21)
    !===================================================================
    ! DD1(21) Ar-Bl- ACT -C"-
    if ((LINELP /= 14) .or. (nlg2 /= 2)) return
    call DD1_ArBl_ACT_C_EXT_AB_SGT0(3)

  case (22)
    !===================================================================
    ! D1D(22) Ar-Bl- Drl- ACT -C"-
    if ((LINELP /= 14) .or. (nlg2 /= 2)) return
    call D1D_ArBl_ACT_C_EXT_AB_SGT0(3)
    call D1D_Drl_ACT_C_EXT_AB_SGT0(3)

  case (23)
    !===================================================================
    ! DV(23) ACT -C'-..................................................
    if (LINELP /= 17) return
    LRA = NLG1
    do LRI=NORB_FRZ+1,NORB_DZ
      LMI = LSM_INN(LRI)
      if (LMI /= JMLR) cycle
      W0DV1 = W0_DV(1)
      NI = mod(NORB_DZ-LRI,2)
      if (NI == 1) W0DV1 = -W0DV1
      ! DV(23-1) A&r(23)-
      IWDL = JUD(LRI)
      IWDR = 0
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      do MPL=1,MTYPE
        VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W0DV1
      end do
      call Ar_BL_EXT_TS(LRI,LRA,1)
    end do

  case (24)
    !===================================================================
    ! D1V(24) Ar- ACT -Bl-
    if (LINELP /= 17) return
    LRA = NLG1
    call D1V_Ar_ACT_Bl_EXT_AB_SGT0(3,LRA)

  case (25)
    !===================================================================
    ! VV(25) ACT -BL- ..................................................

  case (7,9:10,14,16:18,26)
end select

return

end subroutine ts_ext_head_in_dbl

subroutine ts_ext_head_in_act()

use gugaci_global, only: linelp, logic_dh, nlg1, nlg2
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: lrai, lraj

LOGIC_DH = .false.
LRAI = NLG1
LRAJ = NLG2
!line=5 A&r-B&l<-->EXT
if (LINELP == 5) then
  call AR_BL_EXT_TS(lrai,lraj,1)
end if
!line=9 D&rl<-->EXT
if (LINELP == 9) then
  call Drl_TS_EXT(lrai)
end if

return

end subroutine ts_ext_head_in_act

subroutine tt_drt_ci_new()

use gugaci_global, only: idisk_array, idisk_lp, idownwei_g131415, iml, imr, ipae, ipael, iseg_downwei, jml, jmr, jpad, jpadl, &
                         jpadlr, linelp, lpblock_tt, nvalue_space_ss
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: imlr, jptyl, jptyr, lpb

call external_space_plpmode_value_tt()
idisk_lp = idisk_array(8)

do lpb=1,lpblock_tt
  call read_lp()
  IpaeL = iml+9
  Ipae = imr+9
  call logicg_st(iml,imr,3,3)      ! irtype=3(T),lptype=5:ArBl-
  imlr = Mul(iml,imr)
  !if (imlr == 1) imspace = iml
  nvalue_space_ss = iseg_downwei(9+imlr)
  idownwei_g131415 = iseg_downwei(9+iml)
  call get_jpty(jpadlr,jptyl,jptyr)
  call get_jp(jptyl,jml,jpadl,1)
  call get_jp(jptyr,jmr,jpad,1)
  !JMLR = Mul(JML,JMR)
  if (linelp <= 12) then
    call tt_ext_head_in_act()
  else
    call tt_ext_head_in_dbl()
  end if
end do

return

end subroutine tt_drt_ci_new

subroutine tt_ext_head_in_dbl()

use gugaci_global, only: ipae, ipael, jb_sys, jml, jmr, jpad, jpadl, jpadlr, jud, just, linelp, logic_dh, lp_lwei, lp_rwei, &
                         lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, nlg1, nlg2, norb_dz, norb_frz, vplp_w0, vplp_w1, &
                         vplpnew_w0, vplpnew_w1, w0_dd, w0_dv, w0_ss, w0_tt, w0_vv, w1_dd, w1_ss, w1_st, w1_tt
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, jmlr, lmi, lmij, lmj, lpok, lr0, lra, lri, lrj, lrk, mpl, ni, nk
real(kind=wp) :: w0dd1, w0dd2, w0dd3, w0dv1, w0ss15, w0ss17, w0ss20, w0tt2, w0tt3, w1dd1, w1dd2, w1ss15, w1ss17, w1tt2
integer(kind=iwp), external :: iwalk_ad

LOGIC_DH = .true.
JMLR = Mul(JML,JMR)
LPOK = JPADLR
select case (LPOK)
  case (1)
    !===================================================================
    ! SS(1)   ACT -C"-
    !-------------------------------------------------------------------
    ! SS(1-1)  Ar(01)-Bl(32)-        ACT -C"-
    ! SS(1-3)  Ar(13)-Bl(20)-        ACT -C"-
    ! SS(1-6)  (11)-Ar(23)-Bl(32)-   ACT -C"-
    ! SS(1-7)  Ar(13)-C'(21)-Bl(32)- ACT -C"-
    ! SS(1-8)  Ar(13)-C'(22)-Bl(31)- ACT -C"-
    ! SS(1-9)  Ar(23)-C'(11)-Bl(32)- ACT -C"-
    ! SS(1-11) Ar(13)-Bl(31)-C"(22)- ACT -C"-
    ! SS(1-12) Ar(13)-Bl(32)-C"(21)- ACT -C"-
    ! SS(1-13) Ar(23)-Bl(31)-C"(12)- ACT -C"-
    ! SS(1-16) (11)-Drl(22)-         ACT -C"-
    ! SS(1-18) Drl(11)-C"(22)-       ACT -C"-
    ! SS(1-19) Drl(12)-C"(21)-       ACT -C"-
    ! SS(1-20) (11)-Drl(33)-C"(22)-  ACT -C"-
    ! SS(1-20) Drl(33)-C"(11)-C"(22)-ACT -C"-
    !-------------------------------------------------------------------
    ! SS(1)   ACT 14: -C"-
    if (LINELP /= 14) return
    if (NLG2 == 1) then
      if (JB_SYS > 0) call SS_Drl_ACT_C_EXT_AB_SGT0(11)
      W0SS15 = W0_SS(15)
      W1SS15 = W1_SS(15)
      W0SS17 = W0_SS(17)
      W1SS17 = W1_SS(17)
      W0SS20 = W0_SS(20)
      if ((JML == 1) .and. (JMR == 1)) then
        ! SS(1-20) Drl(33)-C"(00)-       ACT -C"-                  ! IPL(R)AD=
        do LR0=NORB_FRZ+1,NORB_DZ
          IWDL = JUST(LR0,LR0)
          IWDR = IWDL
          do MPL=1,MHLP
            IWAL = LPNEW_LWEI(MPL)
            IWAR = LPNEW_RWEI(MPL)
            LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          end do
          do MPL=1,MTYPE
            VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS20
            VPLP_W1(MPL) = Zero
          end do
          !do lrk=1,norb_dz
          !  if (lrk == lr0) cycle
          !  call Drl_tt_ext(lrk)
          !end do
          call Drl_TT_SUM(LR0,0)
        end do
      end if
      do LRI=NORB_FRZ+1,NORB_DZ-1
        LMI = LSM_INN(LRI)
        do LRJ=LRI+1,NORB_DZ
          LMJ = LSM_INN(LRJ)
          LMIJ = Mul(LMI,LMJ)
          if ((LMIJ /= JML) .or. (LMIJ /= JMR)) cycle
          IWDL = JUST(LRI,LRJ)
          IWDR = IWDL
          do MPL=1,MHLP
            IWAL = LPNEW_LWEI(MPL)
            IWAR = LPNEW_RWEI(MPL)
            LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          end do
          ! SS(1-15) (22)-Drl(11)-         ACT -C"-
          do MPL=1,MTYPE
            VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS15
            VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SS15
          end do
          call Drl_TT_EXT(LRJ)
          ! SS(1-17) Drl(22)-C"(11)-       ACT -C"-
          do MPL=1,MTYPE
            VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS17
            VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SS17
          end do
          call Drl_TT_EXT(LRI)
          ! SS(1-20) (22)(11)Drl(33)-      ACT -C"-
          ! SS(1-20) (22)Drl(33)-C"(11)-   ACT -C"-
          ! SS(1-20) Drl(33)-C"(22)-C"(11)-ACT -C"-
          do MPL=1,MTYPE
            VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS20
            VPLP_W1(MPL) = Zero
          end do
          call Drl_TT_SUM(LRI,LRJ)
        end do
      end do
    else
      if (JB_SYS > 0) then
        call SS_ArBl_ACT_C_EXT_AB_SGT0(11)
        call SS_S_Drl_ACT_C_EXT_AB_SGT0(11)
      end if
      do LRI=NORB_FRZ+1,NORB_DZ-1
        do LRJ=LRI+1,NORB_DZ
          call SS2_EXT(LRI,LRJ,NK)
          if (NK /= 0) call Ar_BL_EXT_TT(LRI,LRJ,1)
          call SS4_EXT(LRI,LRJ,NK)
          if (NK /= 0) call Ar_BL_EXT_TT(LRI,LRJ,1)
          call SS5_EXT(LRI,LRJ,NK)
          if (NK /= 0) call Ar_BL_EXT_TT(LRI,LRJ,NK)
          call SS10_EXT(LRI,LRJ,NK)
          if (NK /= 0) call Ar_BL_EXT_TT(LRI,LRJ,NK)
          call SS14_EXT(LRI,LRJ,NK)
          if (NK /= 0) call Ar_BL_EXT_TT(LRI,LRJ,NK)
        end do
      end do
    end if

  case (2)
    !===================================================================
    ! ST(2)   ACT -C"-
    if ((LINELP /= 14) .or. (NLG2 == 1)) return
    if (JB_SYS > 0) call ST_Drl_ACT_C_EXT_AB_SGT0(11)
    if (JB_SYS > 0) call ST_ArBl_ACT_C_EXT_AB_SGT0(11)
    do LRI=NORB_FRZ+1,NORB_DZ-1
      do LRJ=LRI+1,NORB_DZ
        call ST1_EXT(LRI,LRJ,NK)
        if (NK /= 0) call Ar_BL_EXT_TT(LRI,LRJ,1)
        call ST2_EXT(LRI,LRJ,NK)
        if (NK /= 0) call Ar_BL_EXT_TT(LRI,LRJ,NK)
        call ST4_EXT(LRI,LRJ,NK,1)
        if (NK /= 0) call Ar_BL_EXT_TT(LRI,LRJ,NK)
        call ST4_EXT(LRI,LRJ,NK,-1)
        if (NK /= 0) call Ar_BL_EXT_TT(LRI,LRJ,NK)
      end do
    end do

    do LRI=NORB_FRZ+1,NORB_DZ-1
      LMI = LSM_INN(LRI)
      do LRJ=LRI+1,NORB_DZ
        LMJ = LSM_INN(LRJ)
        LMIJ = Mul(LMI,LMJ)
        if ((LMIJ /= JML) .or. (LMIJ /= JMR)) cycle
        ! ST(2-5) (22)Drl(12)-          ACT -C"-
        IWDL = JUST(LRI,LRJ)
        IWDR = IWDL
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = Zero
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_ST(5)
        end do
        call Drl_TT_EXT(LRJ)
        ! ST(2-6) Drl(22)-C"(12)-       ACT -C"-
        do MPL=1,MTYPE
          VPLP_W0(MPL) = Zero
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_ST(6)
        end do
        call Drl_TT_EXT(LRI)
        !---------------------------------------------------------------
        ! ST(2-3) Ar(13)-C'(22)-Bl(32)-   ACT -C"-
        ! ST(2-3) Ar(13)-Bl(32)-C'(22)-   ACT -C"-
        ! ST(2-7) Drl(12)-C"(22)-         ACT -C"-
        !---------------------------------------------------------------
      end do
    end do

  case (3)
    !===================================================================
    ! TS(3) D&R^L-  ACT -C"-
    if ((LINELP /= 14) .or. (NLG2 == 1)) return
    if (JB_SYS > 0) then
      call TS_ArBl_ACT_C_EXT_AB_SGT0(11)
    end if
    do LRI=NORB_FRZ+1,NORB_DZ-1
      do LRJ=LRI+1,NORB_DZ
        call TS1_EXT(LRI,LRJ,NK)
        if (NK /= 0) call Ar_BL_EXT_TT(LRI,LRJ,1)
        call TS2_EXT(LRI,LRJ,NK,1)
        if (NK /= 0) call Ar_BL_EXT_TT(LRI,LRJ,NK)
        call TS2_EXT(LRI,LRJ,NK,-1)
        if (NK /= 0) call Ar_BL_EXT_TT(LRI,LRJ,NK)
        call TS4_EXT(LRI,LRJ,NK)
        if (NK /= 0) call Ar_BL_EXT_TT(LRI,LRJ,NK)
      end do
    end do

  case (4)
    !===================================================================
    ! ST1(4) Ar-Bl- Drl- ACT -C"-
    if ((LINELP /= 14) .or. (nlg2 /= 2)) return
    call STT_ArBl_ACT_C_EXT_AB_SGT1(11)

  case (5)
    !===================================================================
    ! T1S(5) Ar-Bl- Drl ACT -C"-
    if ((LINELP /= 14) .or. (NLG2 /= 2)) return
    call TTS_Drl_ACT_C_EXT_AB_SGT1(11)
    call TTS_ArBl_ACT_C_EXT_AB_SGT1(11)

  case default ! (6)
    !===================================================================
    ! SD(6-1) ACT -B&L-
    if (LINELP /= 17) return
    LRA = NLG1
    call SD_AR_ACT_BL(11,LRA)
    if (JB_SYS > 0) call SD_AR_ACT_BL_SGT0(11,LRA)

  case (8)
    !===================================================================
    ! SD1(8) Ar ACT -Bl-
    if (LINELP /= 17) return
    LRA = NLG1
    if (JB_SYS > 0) call SDD_AR_ACT_BL_SGT0(11,LRA)

  case (11)
    !===================================================================
    ! TT(11) Drl-  ACT -C"-
    if (LINELP /= 14) return
    if (NLG2 == 1) then
      W0TT2 = W0_TT(2)
      W1TT2 = W1_TT(2)
      W0TT3 = W0_TT(3)
      do LRI=NORB_FRZ+1,NORB_DZ-1
        LMI = LSM_INN(LRI)
        do LRJ=LRI+1,NORB_DZ
          LMJ = LSM_INN(LRJ)
          LMIJ = Mul(LMI,LMJ)
          if ((LMIJ /= JML) .or. (LMIJ /= JMR)) cycle
          ! TT(11-2) (22)Drl(22)-
          ! TT(11-2) Drl(22)-C"(22)-
          do MPL=1,MTYPE
            VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0TT2
            VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1TT2
          end do
          IWDL = JUST(LRI,LRJ)
          IWDR = IWDL
          do MPL=1,MHLP
            IWAL = LPNEW_LWEI(MPL)
            IWAR = LPNEW_RWEI(MPL)
            LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
            LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
          end do
          call Drl_TT_EXT(LRI)
          call Drl_TT_EXT(LRJ)
          ! TT(11-3) Drl(33)-C"(22)-C"(22)-
          ! TT(11-3) (22)Drl(33)-C"(22)-
          ! TT(11-3) (22)(22)Drl(33)-
          do MPL=1,MTYPE
            VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0TT3
            VPLP_W1(MPL) = Zero
          end do
          !do lrk=1,norb_dz
          !  if (lrk == lri) cycle
          !  if (lrk == lrj) cycle
          !  call Drl_tt_ext(lrk)
          !end do
          call Drl_TT_SUM(LRI,LRJ)
        end do
      end do
    else
      do LRI=NORB_FRZ+1,NORB_DZ-1
        do LRJ=LRI+1,NORB_DZ
          call TT1_EXT(LRI,LRJ,NK,1)
          call Ar_BL_EXT_TT(LRI,LRJ,NK)
          call TT1_EXT(LRI,LRJ,NK,-1)
          call Ar_BL_EXT_TT(LRI,LRJ,NK)
        end do
      end do
    end if

  case (12)
    !===================================================================
    ! TTTT(12) Ar-Bl- Drl- ACT -C"-
    if (LINELP /= 14) return
    if (NLG2 == 1) then
      call TTTT_Drl_ACT_C_EXT_AB_SGT0(11)
    else
      call TTTT_ArBl_ACT_C_EXT_AB_SGT0(11)
    end if

  case (13)
    !===================================================================
    ! TD(13) ACT -B&L-
    if (LINELP /= 17) return
    LRA = NLG1
    call TD_AR_ACT_BL(11,LRA)

  case (15)
    !===================================================================
    ! T1D1(15) Ar- ACT -Bl-
    if (LINELP /= 17) return
    LRA = NLG1
    call TTDD_AR_ACT_BL_SGT1(11,LRA)

  case (19)
    !===================================================================
    ! DD(19) ACT -C"- ..................................................
    if (LINELP /= 14) return
    if (NLG2 == 1) then
      W0DD2 = W0_DD(2)
      W1DD2 = W1_DD(2)
      W0DD3 = W0_DD(3)
      do LRI=NORB_FRZ+1,NORB_DZ
        LMI = LSM_INN(LRI)
        if (LMI /= JML) cycle
        ! DD(19-2) Drl(22)-
        do MPL=1,MTYPE
          VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0DD2
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1DD2
        end do
        IWDL = JUD(LRI)
        IWDR = IWDL
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        end do
        call Drl_TT_EXT(LRI)
        ! DD(19-3) (22)Drl(33)-
        ! DD(19-3) Drl(33)-C"(22)-
        do MPL=1,MTYPE
          VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0DD3
          VPLP_W1(MPL) = Zero
        end do
        call Drl_TT_SUM(LRI,0)
        !do lrk=1,norb_dz
        !  if (lrk == lri) cycle
        !  call Drl_tt_ext(lrk)
        !end do
      end do
    else
      do LRI=NORB_FRZ+1,NORB_DZ-1
        LMI = LSM_INN(LRI)
        do LRJ=LRI+1,NORB_DZ
          LMJ = LSM_INN(LRJ)
          LMIJ = Mul(LMI,LMJ)
          if (LMIJ /= JMLR) cycle
          W0DD1 = W0_DD(1)
          W1DD1 = W1_DD(1)
          NI = mod(LRJ-LRI,2)
          if (NI == 0) then
            W0DD1 = -W0DD1
            W1DD1 = -W1DD1
          end if
          if ((LMI == JML) .and. (LMJ == JMR)) then
            ! DD(19-1) Ar(23)-Bl(32)-      ACT -C"-
            do MPL=1,MTYPE
              VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0DD1
              VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1DD1
            end do
            IWDL = JUD(LRI)
            IWDR = JUD(LRJ)
            do MPL=1,MHLP
              IWAL = LPNEW_LWEI(MPL)
              IWAR = LPNEW_RWEI(MPL)
              LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
              LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
            end do
            call Ar_BL_EXT_TT(LRI,LRJ,1)
          end if
        end do
      end do
    end if

  case (20)
    !===================================================================
    ! D1D1(20) Drl- Ar-Bl- ACT -C"-
    if (LINELP /= 14) return
    if (NLG2 == 1) then
      call D1D1_Drl_ACT_C_EXT_AB_SGT0(11)
    else
      call D1D1_ArBl_ACT_C_EXT_AB_SGT0(11)
    end if

  case (21)
    !===================================================================
    ! DD1(21) Ar-Bl- ACT -C"-
    if ((LINELP /= 14) .or. (nlg2 /= 2)) return
    call DD1_ArBl_ACT_C_EXT_AB_SGT0(11)

  case (22)
    !===================================================================
    ! D1D(22) Ar-Bl- Drl- ACT -C"-
    if ((LINELP /= 14) .or. (nlg2 /= 2)) return
    call D1D_ArBl_ACT_C_EXT_AB_SGT0(11)
    call D1D_Drl_ACT_C_EXT_AB_SGT0(11)

  case (23)
    !===================================================================
    ! DV(23) ACT -C'-..................................................
    if (LINELP /= 17) return
    LRA = NLG1
    do LRI=NORB_FRZ+1,NORB_DZ
      LMI = LSM_INN(LRI)
      if (LMI /= JMLR) cycle
      W0DV1 = W0_DV(1)
      NI = mod(NORB_DZ-LRI,2)
      if (NI == 1) W0DV1 = -W0DV1
      ! DV(23-1) A&r(23)-
      IWDL = JUD(LRI)
      IWDR = 0
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      do MPL=1,MTYPE
        VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0DV1
        VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W0DV1
      end do
      call Ar_BL_EXT_TT(LRI,LRA,1)
    end do

  case (24)
    !===================================================================
    ! D1V(24) Ar- ACT -Bl-
    if (LINELP /= 17) return
    LRA = NLG1
    call D1V_Ar_ACT_Bl_EXT_AB_SGT0(11,LRA)

  case (25)
    !===================================================================
    ! VV(25) ACT -BL- ..................................................
    if ((LINELP /= 14) .or. (NLG2 /= 1)) return
    ! VV(25) Drl(33)-
    IWDL = 0
    IWDR = 0
    do MPL=1,MTYPE
      VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0_VV
      VPLP_W1(MPL) = Zero
    end do
    do MPL=1,MHLP
      IWAL = LPNEW_LWEI(MPL)
      IWAR = LPNEW_RWEI(MPL)
      LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
      LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
    end do
    !call Drl_TT_SUM(0,0)
    do lrk=1,norb_dz
      call Drl_tt_ext(lrk)
    end do

  case (7,9:10,14,16:18,26)
end select

return

end subroutine tt_ext_head_in_dbl

subroutine tt_ext_head_in_act()

use gugaci_global, only: linelp, logic_dh, nlg1, nlg2
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: lrai, lraj

LOGIC_DH = .false.
LRAI = NLG1
LRAJ = NLG2
!line=5 A&r-B&l<-->EXT
if (LINELP == 5) then
  call Ar_BL_EXT_TT(lrai,lraj,1)
end if
!line=9 D&rl<-->EXT
if (LINELP == 9) then
  call Drl_TT_EXT(lrai)
end if

return

end subroutine tt_ext_head_in_act

subroutine logicg_st(ilnodesm,irnodesm,iltype,irtype)

use gugaci_global, only: ism_g1415, ism_g2g4, logic_g13, logic_g1415, logic_g2g4a, logic_g2g4b, logic_g34a, logic_g34b, &
                         logic_g35a, logic_g35b, logic_g36a, logic_g36b, lpend34a, lpend34b, lpend35a, lpend35b, lpend36a, &
                         lpend36b, lpsta34a, lpsta34b, lpsta35a, lpsta35b, lpsta36a, lpsta36b
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: ilnodesm, irnodesm, iltype, irtype
integer(kind=iwp) :: iii, ilrsm

ilrsm = Mul(ilnodesm,irnodesm)
iii = 1      !index to determine lwei rwei iposint and nlinkorb
! G2G4a G2G4b G1415 G13
logic_g36a = .false.
logic_g36b = .false.
logic_g35a = .false.
logic_g35b = .false.
logic_g34a = .false.
logic_g34b = .false.

lpsta36a = iii
call do_g36mode(ilrsm,ilnodesm,iii)
lpend36a = iii-4
if (lpend36a >= lpsta36a) logic_g36a = .true.
lpsta35a = iii
call do_g35mode(ilrsm,ilnodesm,iii)
lpend35a = iii-4
if (lpend35a >= lpsta35a) logic_g35a = .true.
lpsta34a = iii
call do_g34mode(ilrsm,ilnodesm,iii)
lpend34a = iii-4
if (lpend34a >= lpsta34a) logic_g34a = .true.
if (ilrsm /= 1) then
  lpsta36b = iii
  call do_g36mode(ilrsm,irnodesm,iii)
  lpend36b = iii-4
  if (lpend36b >= lpsta36b) logic_g36b = .true.
  lpsta35b = iii
  call do_g35mode(ilrsm,irnodesm,iii)
  lpend35b = iii-4
  if (lpend35b >= lpsta35b) logic_g35b = .true.
  lpsta34b = iii
  call do_g34mode(ilrsm,irnodesm,iii)
  lpend34b = iii-4
  if (lpend34b >= lpsta34b) logic_g34b = .true.
else
  logic_g36b = logic_g36a
  lpsta36b = lpsta36a
  lpend36b = lpend36a
  logic_g35b = logic_g35a
  lpsta35b = lpsta35a
  lpend35b = lpend35a
  logic_g34b = logic_g34a
  lpsta34b = lpsta34a
  lpend34b = lpend34a
end if
logic_g2g4a = .false.
logic_g2g4b = .false.
logic_g1415 = .false.
logic_g13 = .false.

if ((irnodesm == 1) .and. (irtype == 4)) then        ! ir:S_1
  logic_g2g4a = .true.
  ism_g2g4 = ilnodesm
end if
if ((ilnodesm == 1) .and. (iltype == 4)) then        ! il:S_1
  logic_g2g4b = .true.
  ism_g2g4 = irnodesm
end if
if (ilnodesm == irnodesm) then                       ! ir:S_1
  logic_g1415 = .true.
  ism_g1415 = ilnodesm
end if

end subroutine logicg_st

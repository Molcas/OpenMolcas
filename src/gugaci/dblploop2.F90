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

subroutine SS_DRL_ACT_BL_SGT0(LIN,LRA)
!=======================================================================
! SS(1)    ACT -BL-
! SS(1-20) Drl(33)-C"(11)-C"(22)-
! SS(1-20) (11)Drl(33)-C"(22)-
! SS(1-20) (11)(22)Drl(33)-
!=======================================================================

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_ss, w1_ss
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lmj, lri, lrj, lrk, mpl
integer(kind=iwp), external :: iwalk_ad

do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    ! SS(1-16) (11)-Drl(22)-
    if ((JML /= Mul(LMI,LMJ)) .or. (JMR /= Mul(LMI,LMJ))) cycle
    IWDL = JUST(LRJ,LRI)
    IWDR = IWDL
    do MPL=1,MHLP
      IWAL = LPNEW_LWEI(MPL)
      IWAR = LPNEW_RWEI(MPL)
      LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
      LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
    end do
    do MPL=1,MTYPE
      VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0_SS(16)
      VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_SS(16)
    end do
    call Drl_BL_EXT_AR_NEW(LIN,LRJ,LRA)
    ! SS(1-18) Drl(11)-C"(22)-
    do MPL=1,MHLP
      IWAL = LPNEW_LWEI(MPL)
      IWAR = LPNEW_RWEI(MPL)
      LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
      LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
    end do
    do MPL=1,MTYPE
      VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0_SS(18)
      VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_SS(18)
    end do
    call Drl_BL_EXT_AR_NEW(LIN,LRI,LRA)
    ! SS(1-20) Drl(33)-C"(11)-C"(22)-
    ! SS(1-20) (11)Drl(33)-C"(22)-
    ! SS(1-20) (11)(22)Drl(33)-
    do MPL=1,MHLP
      IWAL = LPNEW_LWEI(MPL)
      IWAR = LPNEW_RWEI(MPL)
      LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
      LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
    end do
    do MPL=1,MTYPE
      VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0_SS(20)
      VPLP_W1(MPL) = Zero
    end do
    if (LRA > NORB_DZ) then
      call Drl_BL_SUM_AR_new(LIN,LRI,LRJ,LRA)
    else
      do LRK=1,NORB_DZ
        if (LRK == LRI) cycle
        if (LRK == LRJ) cycle
        call Drl_BL_EXT_AR_NEW(LIN,LRK,LRA)
      end do
    end if
  end do
end do

return

end subroutine SS_DRL_ACT_BL_SGT0

subroutine SS_ARBL_ACT_BL_SGT0(LIN,LRA)
!=======================================================================
! SS(1)    ACT -BL-
! SS(1-1)  Ar(01)-Bl(32)-
! SS(1-3)  Ar(13)-Bl(20)-
! SS(1-6)  (11)-Ar(23)-Bl(32)-
! SS(1-7)  Ar(13)-C'(21)-Bl(32)-
! SS(1-8)  Ar(13)-C'(22)-Bl(31)-
! SS(1-9)  Ar(23)-C'(11)-Bl(32)-
! SS(1-11) Ar(13)-Bl(31)-C"(22)-
! SS(1-12) Ar(13)-Bl(32)-C"(21)-
! SS(1-13) Ar(23)-Bl(31)-C"(12)-
!=======================================================================

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, lsm_inn, mhlp, mtype, ngw2, ngw3, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, &
                         vplpnew_w1, w0_ss, w1_ss
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, jmlr, lmi, lmij, lmj, lmk, lri, lrj, lrk, mpl, ni
real(kind=wp) :: w0ss1, w0ss11, w0ss12, w0ss13, w0ss3, w0ss6, w0ss7, w0ss8, w0ss9, w1ss1, w1ss11, w1ss12, w1ss13, w1ss3, w1ss6, &
                 w1ss7, w1ss8, w1ss9
integer(kind=iwp), external :: iwalk_ad

JMLR = Mul(JML,JMR)
ISMA = Mul(IML,IMR)
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    if (LMIJ /= JMLR) cycle
    IJK = LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ)
    intpos = INTIND_IJKA(IJK)
    !-------------------------------------------------------------------
    W0SS1 = W0_SS(1)
    W1SS1 = W1_SS(1)
    W0SS3 = W0_SS(3)
    W1SS3 = W1_SS(3)
    W0SS6 = W0_SS(6)
    W1SS6 = W1_SS(6)
    W0SS7 = W0_SS(7)
    W1SS7 = W1_SS(7)
    W0SS8 = W0_SS(8)
    W1SS8 = W1_SS(8)
    W0SS9 = W0_SS(9)
    W1SS9 = W1_SS(9)
    W0SS11 = W0_SS(11)
    W1SS11 = W1_SS(11)
    W0SS12 = W0_SS(12)
    W1SS12 = W1_SS(12)
    W0SS13 = W0_SS(13)
    W1SS13 = W1_SS(13)
    NI = mod(LRJ-LRI,2)
    if (NI == 0) then
      W0SS1 = -W0SS1
      W1SS1 = -W1SS1
      W0SS3 = -W0SS3
      W1SS3 = -W1SS3
      W0SS6 = -W0SS6
      W1SS6 = -W1SS6
      W0SS7 = -W0SS7
      W1SS7 = -W1SS7
      W0SS8 = -W0SS8
      W1SS8 = -W1SS8
      W0SS9 = -W0SS9
      W1SS9 = -W1SS9
      W0SS11 = -W0SS11
      W1SS11 = -W1SS11
      W0SS12 = -W0SS12
      W1SS12 = -W1SS12
      W0SS13 = -W0SS13
      W1SS13 = -W1SS13
    end if
    ! SS(1-1)  Ar(01)-Bl(32)-
    if (JML == 1) then
      IWDL = JUST(LRI,LRI)
      IWDR = JUST(LRJ,LRI)
      do MPL=1,MTYPE
        VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS1
        VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SS1
      end do
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
    end if
    ! SS(1-3)  Ar(13)-Bl(20)-
    if (JMR == 1) then
      IWDL = JUST(LRJ,LRI)
      IWDR = JUST(LRJ,LRJ)
      do MPL=1,MTYPE
        VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS3
        VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SS3
      end do
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
    end if
    !-------------------------------------------------------------------
    ! SS(1-6)  (11)-Ar(23)-Bl(32)-
    do LRK=NORB_FRZ+1,LRI-1
      LMK = LSM_INN(LRK)
      if (Mul(LMK,LMI) /= JML) cycle
      IWDL = JUST(LRI,LRK)
      IWDR = JUST(LRJ,LRK)
      do MPL=1,MTYPE
        VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS6
        VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SS6
      end do
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
    end do
    !-------------------------------------------------------------------
    ! SS(1-7)  Ar(13)-C'(21)-Bl(32)-
    do LRK=LRI+1,LRJ-1
      LMK = LSM_INN(LRK)
      if (Mul(LMK,LMI) /= JML) cycle
      IWDL = JUST(LRK,LRI)
      IWDR = JUST(LRJ,LRK)
      do MPL=1,MTYPE
        VPLP_W0(MPL) = -VPLPNEW_W0(MPL)*W0SS7
        VPLP_W1(MPL) = -VPLPNEW_W1(MPL)*W1SS7
      end do
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
      ! SS(1-8)  Ar(13)-C'(22)-Bl(31)-
      IWDL = JUST(LRK,LRI)
      IWDR = JUST(LRK,LRJ)
      do MPL=1,MTYPE
        VPLP_W0(MPL) = -VPLPNEW_W0(MPL)*W0SS8
        VPLP_W1(MPL) = -VPLPNEW_W1(MPL)*W1SS8
      end do
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
      ! SS(1-9)  Ar(23)-C'(11)-Bl(32)-
      IWDL = JUST(LRI,LRK)
      IWDR = JUST(LRJ,LRK)
      do MPL=1,MTYPE
        VPLP_W0(MPL) = -VPLPNEW_W0(MPL)*W0SS9
        VPLP_W1(MPL) = -VPLPNEW_W1(MPL)*W1SS9
      end do
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
    end do
    !-------------------------------------------------------------------
    ! SS(1-11) Ar(13)-Bl(31)-C"(22)-
    do LRK=LRJ+1,NORB_DZ
      LMK = LSM_INN(LRK)
      if (Mul(LMK,LMI) /= JML) cycle
      IWDL = JUST(LRK,LRI)
      IWDR = JUST(LRK,LRJ)
      do MPL=1,MTYPE
        VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS11
        VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SS11
      end do
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
      ! SS(1-12) Ar(13)-Bl(32)-C"(21)-
      IWDL = JUST(LRK,LRI)
      IWDR = JUST(LRJ,LRK)
      do MPL=1,MTYPE
        VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS12
        VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SS12
      end do
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
      ! SS(1-13) Ar(23)-Bl(31)-C"(12)-
      IWDL = JUST(LRI,LRK)
      IWDR = JUST(LRK,LRJ)
      do MPL=1,MTYPE
        VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS13
        VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SS13
      end do
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
    end do
    !-------------------------------------------------------------------
  end do
end do

return

end subroutine SS_ARBL_ACT_BL_SGT0

subroutine ST_ARBL_ACT_BL_SGT0(LIN,LRA)
!***********************************************************************
! ST(2-3) Ar(13)-C'(22)-Bl(32)-
! ST(2-3) Ar(13)-Bl(32)-C'(22)-
!***********************************************************************

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, lsm_inn, mhlp, mtype, ngw2, ngw3, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w1, w1_st
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, jmlr, lmi, lmij, lmj, lmk, lri, lrj, lrk, mpl, ni
real(kind=wp) :: w1st3
integer(kind=iwp), external :: iwalk_ad

JMLR = Mul(JML,JMR)
ISMA = Mul(IML,IMR)
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    if (LMIJ /= JMLR) cycle

    W1ST3 = W1_ST(3)
    NI = mod(LRJ-LRI,2)
    if (NI == 0) W1ST3 = -W1ST3

    IJK = LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ)
    intpos = INTIND_IJKA(IJK)
    !-------------------------------------------------------------------
    ! ST(2-3) Ar(13)-C'(22)-Bl(32)-
    do LRK=LRI+1,LRJ-1
      LMK = LSM_INN(LRK)
      if (Mul(LMK,LMI) /= JML) cycle
      IWDL = JUST(LRK,LRI)
      IWDR = JUST(LRK,LRJ)
      do MPL=1,MTYPE
        VPLP_W0(MPL) = Zero
        VPLP_W1(MPL) = -VPLPNEW_W1(MPL)*W1ST3
      end do
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
    end do
    ! ST(2-3) Ar(13)-Bl(32)-C'(22)-
    do LRK=LRJ+1,NORB_DZ
      LMK = LSM_INN(LRK)
      if (Mul(LMK,LMI) /= JML) cycle
      IWDL = JUST(LRK,LRI)
      IWDR = JUST(LRJ,LRK)
      do MPL=1,MTYPE
        VPLP_W0(MPL) = Zero
        VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1ST3
      end do
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
    end do
  end do
end do

return

end subroutine ST_ARBL_ACT_BL_SGT0

subroutine ST_DRL_ACT_BL_SGT0(LIN,LRA)
! ST(2-7) Drl(12)-C"(22)-

use gugaci_global, only: ipae, ipael, jml, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w1, w1_st
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lmij, lmj, lri, lrj, mpl
integer(kind=iwp), external :: iwalk_ad

do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    if (LMIJ /= JML) cycle
    !-------------------------------------------------------------------
    ! ST(2-7) Drl(12)-C"(22)-
    IWDL = JUST(LRJ,LRI)
    IWDR = JUST(LRI,LRJ)
    do MPL=1,MHLP
      IWAL = LPNEW_LWEI(MPL)
      IWAR = LPNEW_RWEI(MPL)
      LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
      LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
    end do
    do MPL=1,MTYPE
      VPLP_W0(MPL) = Zero
      VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_ST(7)
    end do
    call Drl_BL_EXT_AR_NEW(LIN,LRI,LRA)
  end do
end do

return

end subroutine ST_DRL_ACT_BL_SGT0

subroutine TS_ARBL_ACT_BL_SGT0(LIN,LRA)
!=======================================================================
! TS(3) A&R-B^L-  ACT -B&L ............................................
! TS(3-3) Ar(23)-Bl(31)-C"(22)-

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, lsm_inn, mhlp, mtype, ngw2, ngw3, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w1, w1_ts
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, jmlr, lmi, lmij, lmj, lmk, lri, lrj, lrk, mpl, ni
real(kind=wp) :: w1ts3
integer(kind=iwp), external :: iwalk_ad

JMLR = Mul(JML,JMR)
ISMA = Mul(IML,IMR)
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    if (LMIJ /= JMLR) cycle
    IJK = LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ)
    intpos = INTIND_IJKA(IJK)
    !-------------------------------------------------------------------
    W1TS3 = W1_TS(3)
    NI = mod(LRJ-LRI,2)
    if (NI == 0) then
      W1TS3 = -W1TS3
    end if
    ! TS(3-3) Ar(23)-Bl(31)-C"(22)-
    do MPL=1,MTYPE
      VPLP_W0(MPL) = Zero
      VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1TS3
    end do
    do LRK=LRJ+1,NORB_DZ
      LMK = LSM_INN(LRK)
      if (Mul(LMK,LMI) /= JML) cycle
      IWDL = JUST(LRI,LRK)
      IWDR = JUST(LRK,LRJ)
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
    end do
  end do
end do

return

end subroutine TS_ARBL_ACT_BL_SGT0

subroutine STT_ARBL_ACT_BL_SGT1(LIN,LRA)
!=======================================================================
! STT(4) A&R-B^L-  ACT -B&L ............................................
! ST1(4-1) Ar(01)-Bl(31)-
! ST1(4-2) Ar(23)-Bl(31)-
! ST1(4-3) Ar(13)-C'(21)-Bl(31)-
! ST1(4-3) Ar(13)-Bl(31)-C"(21)-
! ST1(4-4) Ar(23)-C'(11)-Bl(31)-
! ST1(4-4) Ar(23)-Bl(31)-C"(11)-

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, lsm_inn, mhlp, mtype, ngw2, ngw3, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w1, w1_st1
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, lmi, lmij, lmj, lmk, lri, lrj, lrk, mpl
real(kind=wp) :: w1st1, w1st2, w1st3, w1st4
integer(kind=iwp), external :: iwalk_ad

ISMA = Mul(IML,IMR)
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    W1ST1 = W1_ST1(1)
    W1ST2 = W1_ST1(2)
    W1ST3 = W1_ST1(3)
    W1ST4 = W1_ST1(4)
    if (mod(LRJ-LRI,2) == 0) then
      W1ST1 = -W1ST1
      W1ST2 = -W1ST2
      W1ST3 = -W1ST3
      W1ST4 = -W1ST4
    end if
    IJK = LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ) !???
    intpos = INTIND_IJKA(IJK)                         !???
    ! ST1(4-1) Ar(01)-Bl(31)-
    ! ST1(4-2) Ar(23)-Bl(31)-
    if ((JML == 1) .and. (JMR == LMIJ)) then
      IWDL = JUST(LRI,LRI)
      IWDR = JUST(LRI,LRJ)
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      do MPL=1,MTYPE
        VPLP_W0(MPL) = Zero
        VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1ST1
      end do
      call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
    end if
    ! ST1(4-2) (11)Ar(23)-Bl(31)-
    do LRK=NORB_FRZ+1,LRI-1
      LMK = LSM_INN(LRK)
      if ((JML == Mul(LMK,LMI)) .and. (JMR == Mul(LMK,LMJ))) then
        IWDL = JUST(LRI,LRK)
        IWDR = JUST(LRK,LRJ)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = Zero
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1ST2
        end do
        call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
      end if
    end do
    ! ST1(4-3) Ar(13)-C'(21)-Bl(31)-
    do LRK=LRI+1,LRJ-1
      LMK = LSM_INN(LRK)
      if ((JML == Mul(LMI,LMK)) .and. (JMR == Mul(LMK,LMJ))) then
        IWDL = JUST(LRK,LRI)
        IWDR = JUST(LRK,LRJ)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = Zero
          VPLP_W1(MPL) = -VPLPNEW_W1(MPL)*W1ST3
        end do
        call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
        ! ST1(4-4) Ar(23)-C'(11)-Bl(31)-
        IWDL = JUST(LRI,LRK)
        IWDR = JUST(LRK,LRJ)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = Zero
          VPLP_W1(MPL) = -VPLPNEW_W1(MPL)*W1ST4
        end do
        call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
      end if
    end do
    ! ST1(4-3) Ar(13)-Bl(31)-C"(21)-
    do LRK=LRJ+1,NORB_DZ
      LMK = LSM_INN(LRK)
      if ((JML == Mul(LMI,LMK)) .and. (JMR == Mul(LMJ,LMK))) then
        IWDL = JUST(LRK,LRI)
        IWDR = JUST(LRJ,LRK)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = Zero
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1ST3
        end do
        call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
        ! ST1(4-4) Ar(23)-Bl(31)-C"(11)-
        IWDL = JUST(LRI,LRK)
        IWDR = JUST(LRJ,LRK)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = Zero
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1ST4
        end do
        call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
      end if
    end do
  end do
end do

return

end subroutine STT_ARBL_ACT_BL_SGT1

subroutine TTS_ARBL_ACT_BL_SGT1(LIN,LRA)
!=======================================================================
! TTS(5) A&R-B^L-  ACT -B&L ............................................
! T1S(5-1)   Ar(13)-Bl(10)-
! T1S(5-2)   (11)Ar(13)-Bl(32)-
! T1S(5-2)   Ar(13)-C'(11)-Bl(32)-
! T1S(5-3)   Ar(13)-Bl(31)-C"(12)-
! T1S(5-4)   Ar(13)-Bl(32)-C"(11)-

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, lsm_inn, mhlp, mtype, ngw2, ngw3, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w1, w1_t1s
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, lmi, lmij, lmj, lmk, lri, lrj, lrk, mpl
real(kind=wp) :: w1ts1, w1ts2, w1ts3, w1ts4
integer(kind=iwp), external :: iwalk_ad

ISMA = Mul(IML,IMR)
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    W1TS1 = W1_T1S(1)
    W1TS2 = W1_T1S(2)
    W1TS3 = W1_T1S(3)
    W1TS4 = W1_T1S(4)
    if (mod(LRJ-LRI,2) == 0) then
      W1TS1 = -W1TS1
      W1TS2 = -W1TS2
      W1TS3 = -W1TS3
      W1TS4 = -W1TS4
    end if
    IJK = LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ) !???
    intpos = INTIND_IJKA(IJK)                         !???
    ! T1S(5-1)   Ar(13)-Bl(10)-
    if ((JMR == 1) .and. (JML == LMIJ)) then
      IWDL = JUST(LRI,LRJ)
      IWDR = JUST(LRJ,LRJ)
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      do MPL=1,MTYPE
        VPLP_W0(MPL) = Zero
        VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1TS1
      end do
      call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
    end if
    ! T1S(5-2)   (11)Ar(13)-Bl(32)-
    do LRK=NORB_FRZ+1,LRI-1
      LMK = LSM_INN(LRK)
      if ((JML == Mul(LMK,LMI)) .and. (JMR == Mul(LMK,LMJ))) then
        IWDL = JUST(LRK,LRI)
        IWDR = JUST(LRJ,LRK)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = Zero
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1TS2
        end do
        call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
      end if
    end do
    ! T1S(5-2)   Ar(13)-C'(11)-Bl(32)-
    do LRK=LRI+1,LRJ-1
      LMK = LSM_INN(LRK)
      if ((JML == Mul(LMI,LMK)) .and. (JMR == Mul(LMK,LMJ))) then
        IWDL = JUST(LRI,LRK)
        IWDR = JUST(LRJ,LRK)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = Zero
          VPLP_W1(MPL) = -VPLPNEW_W1(MPL)*W1TS2
        end do
        call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
      end if
    end do
    ! T1S(5-3)   Ar(13)-Bl(31)-C"(12)-
    do LRK=LRJ+1,NORB_DZ
      LMK = LSM_INN(LRK)
      if ((JML == Mul(LMI,LMK)) .and. (JMR == Mul(LMJ,LMK))) then
        IWDL = JUST(LRI,LRK)
        IWDR = JUST(LRK,LRJ)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = Zero
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1TS3
        end do
        call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
        ! T1S(5-4)   Ar(13)-Bl(32)-C"(11)-
        IWDL = JUST(LRI,LRK)
        IWDR = JUST(LRJ,LRK)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = Zero
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1TS4
        end do
        call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
      end if
    end do
  end do
end do

return

end subroutine TTS_ARBL_ACT_BL_SGT1

subroutine TTS_Drl_ACT_BL_SGT1(LIN,LRA)
! T1S(5-5)   (11)Drl(12)-
! T1S(5-6)   Drl(12)-C"(12)-
! T1S(5-7)   Drl(12)-C"(11)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w1, w1_t1s
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lmj, lri, lrj, mpl
integer(kind=iwp), external :: iwalk_ad

do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    if ((JML /= Mul(LMI,LMJ)) .or. (JMR /= Mul(LMI,LMJ))) cycle
    ! T1S(5-5)   (11)Drl(12)-
    IWDL = JUST(LRI,LRJ)
    IWDR = JUST(LRJ,LRI)
    do MPL=1,MHLP
      IWAL = LPNEW_LWEI(MPL)
      IWAR = LPNEW_RWEI(MPL)
      LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
      LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
    end do
    do MPL=1,MTYPE
      VPLP_W0(MPL) = Zero
      VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_T1S(5)
    end do
    call Drl_BL_EXT_AR_NEW(LIN,LRJ,LRA)
    ! T1S(5-6)   Drl(11)-C"(12)-
    do MPL=1,MHLP
      IWAL = LPNEW_LWEI(MPL)
      IWAR = LPNEW_RWEI(MPL)
      LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
      LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
    end do
    do MPL=1,MTYPE
      VPLP_W0(MPL) = Zero
      VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_T1S(6)
    end do
    call Drl_BL_EXT_AR_NEW(LIN,LRI,LRA)
    ! T1S(5-7)   Drl(12)-C"(11)-
    IWDL = JUST(LRI,LRJ)
    IWDR = IWDL
    do MPL=1,MHLP
      IWAL = LPNEW_LWEI(MPL)
      IWAR = LPNEW_RWEI(MPL)
      LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
      LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
    end do
    do MPL=1,MTYPE
      VPLP_W0(MPL) = Zero
      VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_T1S(7)
    end do
    call Drl_BL_EXT_AR_NEW(LIN,LRI,LRA)
  end do
end do

return

end subroutine TTS_Drl_ACT_BL_SGT1

subroutine SDD_ABB_ACT_C_SGT0(LIN)
!SD1(8-5)    Ar(13)-Br(23)-BR(31)-
!SD1(8-6)    Ar(23)-Br(13)-BR(31)-
!SD1(8-7)    Ar(13)-Bl(31)-BL(23)-
!SD1(8-8)    Ar(23)-Bl(31)-BL(13)-

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, lsm_inn, mhlp, mtype, ngw2, ngw3, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, &
                         vplpnew_w1, w0_sd1, w1_sd1
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, lmi, lmij, lmj, lmk, lri, lrj, lrk, mpl, ni
real(kind=wp) :: w0sd5, w0sd6, w0sd7, w0sd8, w1sd5, w1sd6, w1sd7, w1sd8
integer(kind=iwp), external :: iwalk_ad

ISMA = Mul(IML,IMR)
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    do LRK=LRJ+1,NORB_DZ
      LMK = LSM_INN(LRK)
      W0SD5 = W0_SD1(5)
      W1SD5 = W1_SD1(5)
      W0SD6 = W0_SD1(6)
      W1SD6 = W1_SD1(6)
      W0SD7 = W0_SD1(7)
      W1SD7 = W1_SD1(7)
      W0SD8 = W0_SD1(8)
      W1SD8 = W1_SD1(8)
      NI = mod(LRJ-LRI+NORB_DZ-LRK,2)
      if (NI == 0) then
        W0SD5 = -W0SD5
        W1SD5 = -W1SD5
        W0SD6 = -W0SD6
        W1SD6 = -W1SD6
        W0SD7 = -W0SD7
        W1SD7 = -W1SD7
        W0SD8 = -W0SD8
        W1SD8 = -W1SD8
      end if
      IJK = LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRK-NORB_FRZ)
      INTPOS = INTIND_IJKA(IJK)

      if ((JML == LMIJ) .and. (JMR == LMK)) then
        ! SD1(8-5)    Ar(13)-Br(23)-BR(31)-
        IWDL = JUST(LRJ,LRI)
        IWDR = JUD(LRK)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SD5
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SD5
        end do
        call Ar_Br_Br_EXT_AR_NEW(LIN,INTPOS,ISMA)
        ! SD1(8-6)    Ar(23)-Br(13)-BR(31)-
        IWDL = JUST(LRI,LRJ)
        IWDR = JUD(LRK)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SD6
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SD6
        end do
        call Ar_Br_Br_EXT_AR_NEW(LIN,INTPOS,ISMA)
      end if
      ! SD1(8-7)    Ar(13)-Bl(31)-BL(23)-
      if ((JML == Mul(LMI,LMK)) .and. (JMR == LMJ)) then
        IWDL = JUST(LRK,LRI)
        IWDR = JUD(LRJ)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SD7
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SD7
        end do
        call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
        ! SD1(8-8)    Ar(23)-Bl(31)-BL(13)-
        IWDL = JUST(LRI,LRK)
        IWDR = JUD(LRJ)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SD8
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SD8
        end do
        call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
      end if
    end do
  end do
end do

return

end subroutine SDD_ABB_ACT_C_SGT0

subroutine SDD_DrlBl_ACT_C_SGT0(LIN)
! SD1(8-9)  Drl(33)-BL(01)-
! SD1(8-10) Drl(11)-BL(23)-
! SD1(8-11) (11)Drl(33)-BL(23)-
! SD1(8-11) Drl(33)-C"(11)-BL(23)-
! SD1(8-12) Drl(33)-BL(13)-C'(21)-
! SD1(8-13) Drl(33)-BL(23)-C'(11)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, &
                         mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_sd1, w1_sd1
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lmj, lmk, lri, lrj, lrk, mpl
real(kind=wp) :: w0sd10, w0sd11, w0sd12, w0sd13, w0sd9, w1sd10, w1sd11, w1sd12, w1sd13, w1sd9
integer(kind=iwp), external :: iwalk_ad

do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  if ((JML /= 1) .or. (JMR /= LMI)) cycle
  W0SD9 = W0_SD1(9)
  W1SD9 = W1_SD1(9)
  if (mod(NORB_DZ-LRI,2) == 1) then
    W0SD9 = -W0SD9
    W1SD9 = -W1SD9
  end if
  do LRK=1,LRI-1
    ! SD1(8-9)    Drl(33)-BL(01)-
    IWDL = JUST(LRI,LRI)
    IWDR = JUD(LRI)
    do MPL=1,MHLP
      IWAL = LPNEW_LWEI(MPL)
      IWAR = LPNEW_RWEI(MPL)
      LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
      LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
    end do
    do MPL=1,MTYPE
      VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SD9
      VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SD9
    end do
    call Drl_BL_EXT_AR_NEW(LIN,LRK,LRI)
  end do
end do
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    W0SD10 = W0_SD1(10)
    W1SD10 = W1_SD1(10)
    W0SD11 = W0_SD1(11)
    W1SD11 = W1_SD1(11)
    W0SD12 = W0_SD1(12)
    W1SD12 = W1_SD1(12)
    W0SD13 = W0_SD1(13)
    W1SD13 = W1_SD1(13)
    if (mod(NORB_DZ-LRJ,2) == 1) then
      W0SD10 = -W0SD10
      W1SD10 = -W1SD10
      W0SD11 = -W0SD11
      W1SD11 = -W1SD11
    end if
    if (mod(NORB_DZ-LRI,2) == 1) then
      W0SD12 = -W0SD12
      W1SD12 = -W1SD12
      W0SD13 = -W0SD13
      W1SD13 = -W1SD13
    end if
    ! SD1(8-10) Drl(11)-BL(23)-
    if ((JML == Mul(LMI,LMJ)) .and. (JMR == LMI)) then
      IWDL = JUST(LRJ,LRI)
      IWDR = JUD(LRI)
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
      end do
      do MPL=1,MTYPE
        VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SD10
        VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SD10
      end do
      call Drl_BL_EXT_AR_NEW(LIN,LRI,LRJ)
    end if
    do LRK=1,LRI-1
      LMK = LSM_INN(LRK)
      ! SD1(8-12) Drl(33)-BL(13)-C'(21)-
      if ((JML == Mul(LMI,LMJ)) .and. (JMR == LMJ)) then
        IWDL = JUST(LRJ,LRI)
        IWDR = JUD(LRJ)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = -VPLPNEW_W0(MPL)*W0SD12
          VPLP_W1(MPL) = -VPLPNEW_W1(MPL)*W1SD12
        end do
        call Drl_BL_EXT_AR_NEW(LIN,LRK,LRI)
        ! SD1(8-13) Drl(33)-BL(23)-C'(11)-
        IWDL = JUST(LRI,LRJ)
        IWDR = JUD(LRJ)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = -VPLPNEW_W0(MPL)*W0SD13
          VPLP_W1(MPL) = -VPLPNEW_W1(MPL)*W1SD13
        end do
        call Drl_BL_EXT_AR_NEW(LIN,LRK,LRI)
      end if
      if ((JML == Mul(LMI,LMJ)) .and. (JMR == LMI)) then
        ! SD1(8-11) Drl(33)-C"(11)-BL(23)-
        do MPL=1,MTYPE
          VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SD11
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SD11   !?????????
        end do
        IWDL = JUST(LRJ,LRI)
        IWDR = JUD(LRI)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
        end do
        call Drl_BL_EXT_AR_NEW(LIN,LRK,LRJ)
      end if
    end do
    do LRK=NORB_FRZ+1,LRI-1
      LMK = LSM_INN(LRK)
      ! SD1(8-11) (11)Drl(33)-BL(23)-
      if ((JML == Mul(LMK,LMJ)) .and. (JMR == LMK)) then
        IWDL = JUST(LRJ,LRK)
        IWDR = JUD(LRK)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SD11
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SD11
        end do
        call Drl_BL_EXT_AR_NEW(LIN,LRI,LRJ)
      end if
    end do
  end do
end do

return

end subroutine SDD_DrlBl_ACT_C_SGT0

subroutine SDD_DrrBR_ACT_C_SGT0(LIN)
! SD1(8-9)    Drr(03)-BR(31)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, &
                         mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, w0_sd1
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmj, lri, lrj, mpl
real(kind=wp) :: w0sd
integer(kind=iwp), external :: iwalk_ad

do LRI=NORB_FRZ+1,NORB_DZ
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    if ((JML /= 1) .or. (JMR /= LMJ)) cycle
    W0SD = W0_SD1(9)
    if (mod(NORB_DZ-LRJ,2) == 1) W0SD = -W0SD
    IWDL = JUST(LRI,LRI)
    IWDR = JUD(LRJ)
    do MPL=1,MHLP
      IWAL = LPNEW_LWEI(MPL)
      IWAR = LPNEW_RWEI(MPL)
      LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
      LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
    end do
    do MPL=1,MTYPE
      VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SD
      VPLP_W1(MPL) = Zero
    end do
    call Drr_BR_EXT_AR(LIN,LRI,LRJ)
  end do
end do

return

end subroutine SDD_DrrBR_ACT_C_SGT0

subroutine SDD_Ar_ACT_BrBR_SGT0(LIN,LRA)
! SD1(8-1)    Ar(01)-
! SD1(8-2)    Ar(23)-
! SD1(8-3)    Ar(13)-C'(21)-
! SD1(8-4)    Ar(23)-C'(11)-

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, lsm_inn, mhlp, mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_sd1, &
                         w1_sd1
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, lmi, lmij, lmj, lri, lrj, mpl, ni
real(kind=wp) :: w0sd1, w0sd2, w0sd3, w0sd4, w1sd1, w1sd2, w1sd3, w1sd4
integer(kind=iwp), external :: iwalk_ad

ISMA = Mul(IML,IMR)
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  W0SD1 = W0_SD1(1)
  W1SD1 = W1_SD1(1)
  W0SD2 = W0_SD1(2)
  W1SD2 = W1_SD1(2)
  W0SD3 = W0_SD1(3)
  W1SD3 = W1_SD1(3)
  W0SD4 = W0_SD1(4)
  W1SD4 = W1_SD1(4)
  NI = mod(NORB_DZ-LRI,2)
  if (NI == 1) then
    W0SD1 = -W0SD1
    W1SD1 = -W1SD1
    W0SD2 = -W0SD2
    W1SD2 = -W1SD2
    W0SD3 = -W0SD3
    W1SD3 = -W1SD3
    W0SD4 = -W0SD4
    W1SD4 = -W1SD4
  end if
  ! SD1(8-1)    Ar(01)-
  if ((JML == 1) .and. (JMR == LMI)) then
    IWDL = JUST(LRI,LRI)
    IWDR = JUD(LRI)
    do MPL=1,MHLP
      IWAL = LPNEW_LWEI(MPL)
      IWAR = LPNEW_RWEI(MPL)
      LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
      LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
    end do
    do MPL=1,MTYPE
      VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SD1
      VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SD1
    end do
    IJK = LRI-NORB_FRZ+LRA
    INTPOS = INTIND_IJKA(IJK)
    call Ar_Br_BR_EXT_AR_NEW(LIN,INTPOS,ISMA)
  end if
  do LRJ=NORB_FRZ+1,LRI-1
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    if (JML /= LMIJ) cycle
    if (JMR == LMJ) then
      ! SD1(8-2)    (11)Ar(23)-
      IWDL = JUST(LRI,LRJ)
      IWDR = JUD(LRJ)
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
      end do
      do MPL=1,MTYPE
        VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SD2
        VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SD2
      end do
      IJK = LRI-NORB_FRZ+LRA
      INTPOS = INTIND_IJKA(IJK)
      call Ar_Br_BR_EXT_AR_NEW(LIN,INTPOS,ISMA)
    end if
  end do
  ! SD1(8-3)    Ar(13)-C'(21)-
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    if (JML /= LMIJ) cycle
    if (JMR == LMJ) then
      IWDL = JUST(LRJ,LRI)
      IWDR = JUD(LRJ)
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
      end do
      do MPL=1,MTYPE
        VPLP_W0(MPL) = -VPLPNEW_W0(MPL)*W0SD3
        VPLP_W1(MPL) = -VPLPNEW_W1(MPL)*W1SD3
      end do
      IJK = LRI-NORB_FRZ+LRA
      INTPOS = INTIND_IJKA(IJK)
      call Ar_Br_BR_EXT_AR_NEW(LIN,INTPOS,ISMA)
      ! SD1(8-4)    Ar(23)-C'(11)-
      IWDL = JUST(LRI,LRJ)
      IWDR = JUD(LRJ)
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
      end do
      do MPL=1,MTYPE
        VPLP_W0(MPL) = -VPLPNEW_W0(MPL)*W0SD4
        VPLP_W1(MPL) = -VPLPNEW_W1(MPL)*W1SD4
      end do
      IJK = LRI-NORB_FRZ+LRA
      INTPOS = INTIND_IJKA(IJK)
      call Ar_Br_BR_EXT_AR_NEW(LIN,INTPOS,ISMA)
    end if
  end do
end do

return

end subroutine SDD_Ar_ACT_BrBR_SGT0

subroutine SDD_Ar_ACT_BlBL_SGT0(LIN,LRA)
! SD1(8-1)    Ar(01)-
! SD1(8-2)    Ar(23)-
! SD1(8-3)    Ar(13)-C'(21)-
! SD1(8-4)    Ar(23)-C'(11)-

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, lsm_inn, mhlp, mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_sd1, &
                         w1_sd1
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, lmi, lmij, lmj, lri, lrj, mpl, ni
real(kind=wp) :: w0sd1, w0sd2, w0sd3, w0sd4, w1sd1, w1sd2, w1sd3, w1sd4
integer(kind=iwp), external :: iwalk_ad

ISMA = Mul(IML,IMR)
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  W0SD1 = W0_SD1(1)
  W1SD1 = W1_SD1(1)
  W0SD2 = W0_SD1(2)
  W1SD2 = W1_SD1(2)
  W0SD3 = W0_SD1(3)
  W1SD3 = W1_SD1(3)
  W0SD4 = W0_SD1(4)
  W1SD4 = W1_SD1(4)
  NI = mod(NORB_DZ-LRI,2)
  if (NI == 1) then
    W0SD1 = -W0SD1
    W1SD1 = -W1SD1
    W0SD2 = -W0SD2
    W1SD2 = -W1SD2
    W0SD3 = -W0SD3
    W1SD3 = -W1SD3
    W0SD4 = -W0SD4
    W1SD4 = -W1SD4
  end if
  ! SD1(8-1)    Ar(01)-
  if ((JML == 1) .and. (JMR == LMI)) then
    IWDL = JUST(LRI,LRI)
    IWDR = JUD(LRI)
    do MPL=1,MHLP
      IWAL = LPNEW_LWEI(MPL)
      IWAR = LPNEW_RWEI(MPL)
      LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
      LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
    end do
    do MPL=1,MTYPE
      VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SD1
      VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SD1
    end do
    IJK = LRI-NORB_FRZ+LRA
    INTPOS = INTIND_IJKA(IJK)
    call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
  end if
  do LRJ=NORB_FRZ+1,LRI-1
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    if (JML /= LMIJ) cycle
    if (JMR == LMJ) then
      ! SD1(8-2)    (11)Ar(23)-
      IWDL = JUST(LRI,LRJ)
      IWDR = JUD(LRJ)
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
      end do
      do MPL=1,MTYPE
        VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SD2
        VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SD2
      end do
      IJK = LRI-NORB_FRZ+LRA
      INTPOS = INTIND_IJKA(IJK)
      call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
    end if
  end do
  ! SD1(8-3)    Ar(13)-C'(21)-
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    if (JML /= LMIJ) cycle
    if (JMR == LMJ) then
      IWDL = JUST(LRJ,LRI)
      IWDR = JUD(LRJ)
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
      end do
      do MPL=1,MTYPE
        VPLP_W0(MPL) = -VPLPNEW_W0(MPL)*W0SD3
        VPLP_W1(MPL) = -VPLPNEW_W1(MPL)*W1SD3
      end do
      IJK = LRI-NORB_FRZ+LRA
      INTPOS = INTIND_IJKA(IJK)
      call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
      ! SD1(8-4)    Ar(23)-C'(11)-
      IWDL = JUST(LRI,LRJ)
      IWDR = JUD(LRJ)
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
      end do
      do MPL=1,MTYPE
        VPLP_W0(MPL) = -VPLPNEW_W0(MPL)*W0SD4
        VPLP_W1(MPL) = -VPLPNEW_W1(MPL)*W1SD4
      end do
      IJK = LRI-NORB_FRZ+LRA
      INTPOS = INTIND_IJKA(IJK)
      call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
    end if
  end do
end do

return

end subroutine SDD_Ar_ACT_BlBL_SGT0

subroutine DDS_ABB_ACT_C_SGT0(LIN)
! D1S(9-2)    Ar(13)-Bl(31)-BR(32)-
! D1S(9-3)    Ar(13)-Bl(32)-BR(31)-

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, lsm_inn, mhlp, mtype, ngw2, ngw3, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w1, w1_d1s
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, lmi, lmj, lmk, lri, lrj, lrk, mpl, ni
real(kind=wp) :: w1ds2, w1ds3
integer(kind=iwp), external :: iwalk_ad

ISMA = Mul(IML,IMR)
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    do LRK=LRJ+1,NORB_DZ
      LMK = LSM_INN(LRK)
      W1DS2 = W1_D1S(2)
      W1DS3 = W1_D1S(3)
      NI = mod(LRJ-LRI+NORB_DZ-LRK,2)
      if (NI == 1) then
        W1DS2 = -W1DS2
        W1DS3 = -W1DS3
      end if
      IJK = LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRK-NORB_FRZ)
      INTPOS = INTIND_IJKA(IJK)
      ! D1S(9-2)    Ar(13)-Bl(31)-BR(32)-
      if ((JML == LMI) .and. (JMR == Mul(LMJ,LMK))) then
        IWDL = JUD(LRI)
        IWDR = JUST(LRK,LRJ)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = Zero
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1DS2
        end do
        call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
        ! D1S(9-3)    Ar(13)-Bl(32)-BR(31)-
        IWDL = JUD(LRI)
        IWDR = JUST(LRJ,LRK)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = Zero
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1DS3
        end do
        call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
      end if
    end do
  end do
end do

return

end subroutine DDS_ABB_ACT_C_SGT0

subroutine TTTT_Drl_ACT_BL_SGT1(LIN,LRA)
! T1T1(12-2)  (11)Drl(11)-
! T1T1(12-2)  Drl(11)-C"(11)-
! T1T1(12-3)  (11)Drl(33)-
! T1T1(12-3)  (11)Drl(33)-C"(11)-
! T1T1(12-3)  Drl(33)-C"(11)-C"(11)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_t1t1, w1_t1t1
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lmij, lmj, lri, lrj, lrk, mpl
integer(kind=iwp), external :: iwalk_ad

if (JML /= JMR) return
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    if (JML /= LMIJ) cycle
    ! T1T1(12-2)  (11)Drl(11)-
    ! T1T1(12-2)  Drl(11)-C"(11)-
    do MPL=1,MTYPE
      VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0_T1T1(2)
      VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_T1T1(2)
    end do
    IWDL = JUST(LRI,LRJ)
    IWDR = IWDL
    do MPL=1,MHLP
      IWAL = LPNEW_LWEI(MPL)
      IWAR = LPNEW_RWEI(MPL)
      LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
      LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
    end do
    call Drl_BL_EXT_AR_NEW(LIN,LRI,LRA)
    call Drl_BL_EXT_AR_NEW(LIN,LRJ,LRA)
    ! T1T1(12-3)  (11)Drl(33)-
    ! T1T1(12-3)  (11)Drl(33)-C"(11)-
    ! T1T1(12-3)  Drl(33)-C"(11)-C"(11)-
    do MPL=1,MTYPE
      VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0_T1T1(3)
      VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_T1T1(3)
    end do
    if (LRA > NORB_DZ) then
      call Drl_BL_SUM_AR_new(LIN,LRI,LRJ,LRA)
    else
      do LRK=1,NORB_DZ
        if (LRK == LRI) cycle
        if (LRK == LRJ) cycle
        call Drl_BL_EXT_AR_NEW(LIN,LRK,LRA)
      end do
    end if
  end do
end do

return

end subroutine TTTT_Drl_ACT_BL_SGT1

subroutine TTTT_ArBl_ACT_BL_SGT1(LIN,LRA)
! T1T1(12-1)  (11)Ar(13)-Bl(31)-
! T1T1(12-1)  Ar(13)-C'(11)-Bl(31)-
! T1T1(12-1)  Ar(13)-Bl(31)-C"(11)-

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, lsm_inn, mhlp, mtype, ngw2, ngw3, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, &
                         vplpnew_w1, w0_t1t1, w1_t1t1
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, lmi, lmj, lmk, lri, lrj, lrk, mpl
real(kind=wp) :: w0tt1, w1tt1
integer(kind=iwp), external :: iwalk_ad

ISMA = Mul(IML,IMR)
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    W0TT1 = W0_T1T1(1)
    W1TT1 = W1_T1T1(1)
    if (mod(LRJ-LRI,2) == 0) then
      W0TT1 = -W0TT1
      W1TT1 = -W1TT1
    end if
    do LRK=NORB_FRZ+1,LRI-1
      LMK = LSM_INN(LRK)
      ! T1T1(12-1)  (11)Ar(13)-Bl(31)-
      if ((JML == Mul(LMK,LMI)) .and. (JMR == Mul(LMK,LMJ))) then
        IWDL = JUST(LRK,LRI)
        IWDR = JUST(LRK,LRJ)
        IJK = LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ)
        INTPOS = INTIND_IJKA(IJK)
        do MPL=1,MTYPE
          VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0TT1
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1TT1
        end do
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
        end do
        call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
      end if
    end do
    do LRK=LRJ+1,NORB_DZ
      LMK = LSM_INN(LRK)
      if ((JML == Mul(LMI,LMk)) .and. (JMR == Mul(LMJ,LMK))) then
        ! T1T1(12-1)  Ar(13)-Bl(31)-C"(11)-
        IWDL = JUST(LRI,LRK)
        IWDR = JUST(LRJ,LRK)
        IJK = LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ)
        INTPOS = INTIND_IJKA(IJK)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0TT1
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1TT1
        end do
        call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
      end if
    end do
    do LRK=LRI+1,LRJ-1
      LMK = LSM_INN(LRK)
      if ((JML == Mul(LMI,LMK)) .and. (JMR == Mul(LMK,LMJ))) then
        ! T1T1(12-1)  Ar(13)-C'(11)-Bl(31)-
        IJK = LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ)
        INTPOS = INTIND_IJKA(IJK)
        IWDL = JUST(LRI,LRK)
        IWDR = JUST(LRK,LRJ)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = -VPLPNEW_W0(MPL)*W0TT1
          VPLP_W1(MPL) = -VPLPNEW_W1(MPL)*W1TT1
        end do
        call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
      end if
    end do
  end do
end do

return

end subroutine TTTT_ArBl_ACT_BL_SGT1

subroutine TTDD_DrlBL_ACT_C_SGT1(LIN)
! T1D1(15-4)  Drl(11)-BL(13)-
! T1D1(15-5)  (11)Drl(33)-BL(13)-
! T1D1(15-5)  Drl(33)-C"(11)-BL(13)-
! T1D1(15-5)  Drl(33)-BL(13)-C'(11)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, &
                         mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_t1d1, w1_t1d1
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lmj, lri, lrj, lrk, mpl
real(kind=wp) :: w0td4, w0td5, w1td4, w1td5
integer(kind=iwp), external :: iwalk_ad

do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    if ((JML == Mul(LMI,LMJ)) .and. (JMR == LMI)) then
      IWDL = JUST(LRI,LRJ)
      IWDR = JUD(LRI)
      W0TD4 = W0_T1D1(4)
      W1TD4 = W1_T1D1(4)
      W0TD5 = W0_T1D1(5)
      W1TD5 = W1_T1D1(5)
      if (mod(NORB_DZ-LRJ,2) == 1) then
        W0TD4 = -W0TD4
        W1TD4 = -W1TD4
        W0TD5 = -W0TD5
        W1TD5 = -W1TD5
      end if
      ! T1D1(15-4)  Drl(11)-BL(13)-
      ! T1D1(15-5)  (11)Drl(33)-BL(13)-
      do MPL=1,MTYPE
        VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0TD4
        VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1TD4
      end do
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
      end do
      call Drl_BL_EXT_AR_NEW(LIN,LRI,LRJ)

      do MPL=1,MTYPE
        VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0TD5
        VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1TD5
      end do
      do LRK=1,LRI-1
        ! T1D1(15-5)  Drl(33)-C"(11)-BL(13)-
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
        end do
        call Drl_BL_EXT_AR_NEW(LIN,LRK,LRJ)
      end do
      ! T1D1(15-5)  (11)Drl(33)-BL(13)-
      do LRK=LRI+1,LRJ-1
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
        end do
        call Drl_BL_EXT_AR_NEW(LIN,LRK,LRJ)
      end do
    end if
    W0TD5 = W0_T1D1(5)
    W1TD5 = W1_T1D1(5)
    if (mod(NORB_DZ-LRI,2) == 0) then
      W0TD5 = -W0TD5
      W1TD5 = -W1TD5
    end if
    if ((JML == Mul(LMI,LMJ)) .and. (JMR == LMJ)) then
      do MPL=1,MTYPE
        VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0TD5
        VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1TD5
      end do
      ! T1D1(15-5)  Drl(33)-BL(13)-C'(11)-
      do LRK=1,LRI-1
        IWDL = JUST(LRI,LRJ)
        IWDR = JUD(LRJ)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
        end do
        call Drl_BL_EXT_AR_NEW(LIN,LRK,LRI)
      end do
    end if
  end do
end do

return

end subroutine TTDD_DrlBL_ACT_C_SGT1

subroutine TTDD_ABB_ACT_C_SGT1(LIN)
! T1D1(15-2)  Ar(13)-Br(13)-BR(31)-
! T1D1(15-3)  Ar(13)-Bl(31)-Bl(13)-

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, lsm_inn, mhlp, mtype, ngw2, ngw3, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, &
                         vplpnew_w1, w0_t1d1, w1_t1d1
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, lmi, lmj, lmk, lri, lrj, lrk, mpl, ni
real(kind=wp) :: w0td3, w1td2, w1td3
integer(kind=iwp), external :: iwalk_ad

ISMA = Mul(IML,IMR)
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    do LRK=LRJ+1,NORB_DZ
      LMK = LSM_INN(LRK)
      W1TD2 = W1_T1D1(2)
      W0TD3 = W0_T1D1(3)
      W1TD3 = W1_T1D1(3)
      NI = mod(LRJ-LRI+NORB_DZ-LRK,2)
      if (NI == 0) then
        W1TD2 = -W1TD2
        W0TD3 = -W0TD3
        W1TD3 = -W1TD3
      end if
      IJK = LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRK-NORB_FRZ)
      INTPOS = INTIND_IJKA(IJK)
      if ((JML == Mul(LMI,LMJ)) .and. (JMR == LMK)) then
        ! T1D1(15-2)  Ar(13)-Br(13)-BR(31)-
        IWDL = JUST(LRI,LRJ)
        IWDR = JUD(LRK)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = Zero
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1TD2
        end do
        call Ar_Br_BR_EXT_AR_NEW(LIN,INTPOS,ISMA)
      end if
      if ((JML == Mul(LMI,LMK)) .and. (JMR == LMJ)) then
        ! T1D1(15-3)  Ar(13)-Bl(31)-Bl(13)-
        IWDL = JUST(LRI,LRK)
        IWDR = JUD(LRJ)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0TD3
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1TD3
        end do
        call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
      end if
    end do
  end do
end do

return

end subroutine TTDD_ABB_ACT_C_SGT1

subroutine TTDD_Ar_ACT_BrBR_SGT0(LIN,LRA)
! T1D1(15-1)  (11)Ar(13)-
! T1D1(15-1)  Ar(13)-C'(11)-

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, lsm_inn, mhlp, mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_t1d1, &
                         w1_t1d1
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, lmi, lmij, lmj, lri, lrj, mpl
real(kind=wp) :: w0td1, w1td1
integer(kind=iwp), external :: iwalk_ad

ISMA = Mul(IML,IMR)
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  W0TD1 = W0_T1D1(1)
  W1TD1 = W1_T1D1(1)
  if (mod(NORB_DZ-LRI,2) == 1) then
    W0TD1 = -W0TD1
    W1TD1 = -W1TD1
  end if
  do LRJ=NORB_FRZ+1,LRI-1
    ! T1D1(15-1)  (11)Ar(13)-
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    if ((JML == LMIJ) .and. (JMR == LMJ)) then
      IWDL = JUST(LRJ,LRI)
      IWDR = JUD(LRJ)
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
      end do
      do MPL=1,MTYPE
        VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0TD1
        VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1TD1
      end do
      IJK = LRI-NORB_FRZ+LRA
      INTPOS = INTIND_IJKA(IJK)
      call Ar_Br_Br_EXT_AR_NEW(LIN,INTPOS,ISMA)
    end if
  end do
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    ! T1D1(15-1)  Ar(13)-C'(11)-
    if ((JML == LMIJ) .and. (JMR == LMJ)) then
      IWDL = JUST(LRI,LRJ)
      IWDR = JUD(LRJ)
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
      end do
      do MPL=1,MTYPE
        VPLP_W0(MPL) = -VPLPNEW_W0(MPL)*W0TD1
        VPLP_W1(MPL) = -VPLPNEW_W1(MPL)*W1TD1
      end do
      IJK = LRI-NORB_FRZ+LRA
      INTPOS = INTIND_IJKA(IJK)
      call Ar_Br_Br_EXT_AR_NEW(LIN,INTPOS,ISMA)
    end if
  end do
end do

return

end subroutine TTDD_Ar_ACT_BrBR_SGT0

subroutine TTDD_Ar_ACT_BlBL_SGT0(LIN,LRA)
! T1D1(15-1)  (11)Ar(13)-
! T1D1(15-1)  Ar(13)-C'(11)-

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, lsm_inn, mhlp, mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_t1d1, &
                         w1_t1d1
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, lmi, lmij, lmj, lri, lrj, mpl
real(kind=wp) :: w0td1, w1td1
integer(kind=iwp), external :: iwalk_ad

ISMA = Mul(IML,IMR)
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  W0TD1 = W0_T1D1(1)
  W1TD1 = W1_T1D1(1)
  if (mod(NORB_DZ-LRI,2) == 1) then
    W0TD1 = -W0TD1
    W1TD1 = -W1TD1
  end if
  do LRJ=NORB_FRZ+1,LRI-1
    ! T1D1(15-1)  (11)Ar(13)-
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    if ((JML == LMIJ) .and. (JMR == LMJ)) then
      IWDL = JUST(LRJ,LRI)
      IWDR = JUD(LRJ)
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
      end do
      do MPL=1,MTYPE
        VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0TD1
        VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1TD1
      end do
      IJK = LRI-NORB_FRZ+LRA
      INTPOS = INTIND_IJKA(IJK)
      call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
    end if
  end do
  do LRJ=LRI+1,NORB_DZ
    ! T1D1(15-1)  Ar(13)-C'(11)-
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    if ((JML == LMIJ) .and. (JMR == LMJ)) then
      IWDL = JUST(LRI,LRJ)
      IWDR = JUD(LRJ)
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
      end do
      do MPL=1,MTYPE
        VPLP_W0(MPL) = -VPLPNEW_W0(MPL)*W0TD1
        VPLP_W1(MPL) = -VPLPNEW_W1(MPL)*W1TD1
      end do
      IJK = LRI-NORB_FRZ+LRA
      INTPOS = INTIND_IJKA(IJK)
      call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
    end if
  end do
end do

return

end subroutine TTDD_Ar_ACT_BlBL_SGT0

subroutine TTV_ArBr_ACT_C_SGT1(LIN,LRA)
! T1V(18) Ar(13)-Br(13)-

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, lsm_inn, mhlp, mtype, ngw2, ngw3, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w1, w1_t1v
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, lmi, lmij, lmj, lri, lrj, mpl
real(kind=wp) :: w1tv1
integer(kind=iwp), external :: iwalk_ad

isma = Mul(iml,imr)
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    if ((JML /= LMIJ) .or. (JMR /= 1)) cycle
    W1TV1 = W1_T1V
    if (mod(LRJ-LRI,2) == 0) W1TV1 = -W1TV1
    IJK = LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ)
    INTPOS = INTIND_IJKA(IJK)
    IWDL = JUST(LRI,LRJ)
    IWDR = 0
    do MPL=1,MHLP
      IWAL = LPNEW_LWEI(MPL)
      IWAR = LPNEW_RWEI(MPL)
      LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
      LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
    end do
    do MPL=1,MTYPE
      VPLP_W0(MPL) = Zero
      VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1TV1
    end do
    call Ar_Br_BR_EXT_AR_NEW(LIN,INTPOS,ISMA)
  end do
end do

return

end subroutine TTV_ArBr_ACT_C_SGT1

subroutine DDDD_ArBl_ACT_BL_SGT0(LIN,LRA)
! D1D1(20-1) Ar(13)-BL(31)-

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, jud, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, &
                         lsm_inn, mhlp, mtype, ngw2, ngw3, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_d1d1, &
                         w1_d1d1
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, lmi, lmj, lri, lrj, mpl
real(kind=wp) :: w0dd1, w1dd1
integer(kind=iwp), external :: iwalk_ad

ISMA = Mul(IML,IMR)
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    if ((JML /= LMI) .or. (JMR /= LMJ)) cycle
    W0DD1 = W0_D1D1(1)
    W1DD1 = W1_D1D1(1)
    if (mod(LRJ-LRI,2) == 0) then
      W0DD1 = -W0DD1
      W1DD1 = -W1DD1
    end if
    IJK = LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ)
    INTPOS = INTIND_IJKA(IJK)
    IWDL = JUD(LRI)
    IWDR = JUD(LRJ)
    do MPL=1,MHLP
      IWAL = LPNEW_LWEI(MPL)
      IWAR = LPNEW_RWEI(MPL)
      LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
      LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
    end do
    do MPL=1,MTYPE
      VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0DD1
      VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1DD1
    end do
    call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
  end do
end do

return

end subroutine DDDD_ArBl_ACT_BL_SGT0

subroutine DDDD_Drl_ACT_BL_SGT0(LIN,LRA)
! D1D1(20-2) Drl(11)-
! D1D1(20-3) (11)Drl(33)-
! D1D1(20-3) Drl(33)-C"(11)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, jud, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_d1d1, w1_d1d1
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lri, lrj, mpl
integer(kind=iwp), external :: iwalk_ad

if (JML /= JMR) return
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  if (JML /= LMI) cycle
  ! D1D1(20-2) Drl(11)-
  IWDL = JUD(LRI)
  IWDR = IWDL
  do MPL=1,MHLP
    IWAL = LPNEW_LWEI(MPL)
    IWAR = LPNEW_RWEI(MPL)
    LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
    LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
  end do
  do MPL=1,MTYPE
    VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0_D1D1(2)
    VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_D1D1(2)
  end do
  call Drl_BL_EXT_AR_NEW(LIN,LRI,LRA)
  do MPL=1,MTYPE
    VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0_D1D1(3)
    VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_D1D1(3)
  end do
  do LRJ=1,NORB_DZ
    if (LRJ == LRI) cycle
    ! D1D1(20-3) Drl(33)-C"(11)-
    ! D1D1(20-3) (11)Drl(33)-
    call Drl_BL_EXT_AR_NEW(LIN,LRJ,LRA)
  end do
end do

return

end subroutine DDDD_Drl_ACT_BL_SGT0

subroutine DD1_ArBl_ACT_BL_SGT0(LIN,LRA)
! DD1(21) Ar(23)-Bl(31)-

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, jud, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, &
                         lsm_inn, mhlp, mtype, ngw2, ngw3, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w1, w1_dd1
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, lmi, lmj, lri, lrj, mpl
real(kind=wp) :: w1dd1
integer(kind=iwp), external :: iwalk_ad

ISMA = Mul(IML,IMR)
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    if ((JML /= LMI) .or. (JMR /= LMJ)) cycle
    W1DD1 = W1_DD1
    if (mod(LRJ-LRI,2) == 0) then
      W1DD1 = -W1DD1
    end if
    IJK = LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ)
    INTPOS = INTIND_IJKA(IJK)
    IWDL = JUD(LRI)
    IWDR = JUD(LRJ)
    do MPL=1,MHLP
      IWAL = LPNEW_LWEI(MPL)
      IWAR = LPNEW_RWEI(MPL)
      LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
      LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
    end do
    do MPL=1,MTYPE
      VPLP_W0(MPL) = Zero
      VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1DD1
    end do
    call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
  end do
end do

return

end subroutine DD1_ArBl_ACT_BL_SGT0

subroutine D1D_ArBl_ACT_BL_SGT0(LIN,LRA)
! DD1(21) Ar(13)-Bl(32)-

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, jud, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, &
                         lsm_inn, mhlp, mtype, ngw2, ngw3, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w1, w1_d1d
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, lmi, lmj, lri, lrj, mpl
real(kind=wp) :: w1dd1
integer(kind=iwp), external :: iwalk_ad

ISMA = Mul(IML,IMR)
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    if ((JML /= LMI) .or. (JMR /= LMJ)) cycle
    W1DD1 = W1_D1D(1)
    if (mod(LRJ-LRI,2) == 0) then
      W1DD1 = -W1DD1
    end if
    IJK = LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ)
    INTPOS = INTIND_IJKA(IJK)
    IWDL = JUD(LRI)
    IWDR = JUD(LRJ)
    do MPL=1,MHLP
      IWAL = LPNEW_LWEI(MPL)
      IWAR = LPNEW_RWEI(MPL)
      LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
      LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
    end do
    do MPL=1,MTYPE
      VPLP_W0(MPL) = Zero
      VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1DD1
    end do
    call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
  end do
end do

return

end subroutine D1D_ArBl_ACT_BL_SGT0

subroutine D1D_Drl_ACT_BL_SGT0(LIN,LRA)
! D1D(22-2) Drl(12)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, jud, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w1, w1_d1d
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lri, mpl
integer(kind=iwp), external :: iwalk_ad

if (JML /= JMR) return
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  if (JML /= LMI) cycle
  ! D1D(22-2) Drl(12)-
  IWDL = JUD(LRI)
  IWDR = IWDL
  do MPL=1,MHLP
    IWAL = LPNEW_LWEI(MPL)
    IWAR = LPNEW_RWEI(MPL)
    LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
    LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
  end do
  do MPL=1,MTYPE
    VPLP_W0(MPL) = Zero
    VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_D1D(2)
  end do
  call Drl_BL_EXT_AR_NEW(LIN,LRI,LRA)
end do

return

end subroutine D1D_Drl_ACT_BL_SGT0

subroutine D1V_Drl_BL_ACT_C_SGT0(LIN)
! D1V(24-1)  Ar(13)-
! D1V(24-2)  Drl(33)-BL(13)-

use gugaci_global, only: iml, imr, ipae, ipael, jml, jmr, jpad, jpadl, jud, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, &
                         mhlp, mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_d1v
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: isma, iwal, iwar, iwdl, iwdr, lmi, lri, lrj, mpl, ni
real(kind=wp) :: w0
integer(kind=iwp), external :: iwalk_ad

ISMA = Mul(IML,IMR)
if (JMR /= 1) return
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  if (JML /= LMI) cycle
  NI = mod(NORB_DZ-LRI,2)
  ! D1V(24-1)  Ar(13)-
  W0 = W0_D1V(1)
  if (NI == 1) W0 = -W0_D1V(1)
  IWDL = JUD(LRI)
  IWDR = 0
  do MPL=1,MTYPE
    VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0
    VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W0
  end do
  do MPL=1,MHLP
    IWAL = LPNEW_LWEI(MPL)
    IWAR = LPNEW_RWEI(MPL)
    LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
    LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
  end do
  call Ar_DV_EXT_AR(51,isma,LRI,0)   !ar_dv
  ! D1V(24-2)  Drl(33)-BL(13)-
  W0 = W0_D1V(2)
  if (NI == 1) W0 = -W0_D1V(2)
  do MPL=1,MTYPE
    VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0
    VPLP_W1(MPL) = Zero
  end do
  do LRJ=1,LRI-1
    call Drl_BL_EXT_AR_NEW(LIN,LRJ,LRI)
  end do
end do

return

end subroutine D1V_Drl_BL_ACT_C_SGT0

subroutine D1V_Ar_ACT_BrBR_SGT0(LIN,LRA)
! D1V(24-1)  Ar(13)-

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, jud, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, &
                         lsm_inn, mhlp, mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_d1v, w1_d1v
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, lmi, lri, mpl
real(kind=wp) :: w0dv1, w1dv1
integer(kind=iwp), external :: iwalk_ad

ISMA = Mul(IML,IMR)
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  if ((JML /= LMI) .or. (JMR /= 1)) cycle
  W0DV1 = W0_D1V(1)
  W1DV1 = W1_D1V(1)
  if (mod(NORB_DZ-LRI,2) == 1) then
    W0DV1 = -W0_D1V(1)
    W1DV1 = -W1_D1V(1)
  end if
  IJK = LRI-NORB_FRZ+LRA
  INTPOS = INTIND_IJKA(IJK)
  IWDL = JUD(LRI)
  IWDR = 0
  do MPL=1,MHLP
    IWAL = LPNEW_LWEI(MPL)
    IWAR = LPNEW_RWEI(MPL)
    LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
    LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
  end do
  do MPL=1,MTYPE
    VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0DV1
    VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1DV1
  end do
  call Ar_Br_Br_EXT_AR_NEW(LIN,INTPOS,ISMA)
end do

return

end subroutine D1V_Ar_ACT_BrBR_SGT0

subroutine D1V_Ar_ACT_BlBL_SGT0(LIN,LRA)
! D1V(24-1)  Ar(13)-

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, jud, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, &
                         lsm_inn, mhlp, mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_d1v, w1_d1v
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, lmi, lri, mpl
real(kind=wp) :: w0dv1, w1dv1
integer(kind=iwp), external :: iwalk_ad

ISMA = Mul(IML,IMR)
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  if ((JML /= LMI) .or. (JMR /= 1)) cycle
  W0DV1 = W0_D1V(1)
  W1DV1 = W1_D1V(1)
  if (mod(NORB_DZ-LRI,2) == 1) then
    W0DV1 = -W0_D1V(1)
    W1DV1 = -W1_D1V(1)
  end if
  IJK = LRI-NORB_FRZ+LRA
  INTPOS = INTIND_IJKA(IJK)
  IWDL = JUD(LRI)
  IWDR = 0
  do MPL=1,MHLP
    IWAL = LPNEW_LWEI(MPL)
    IWAR = LPNEW_RWEI(MPL)
    LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
    LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
  end do
  do MPL=1,MTYPE
    VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0DV1
    VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1DV1
  end do
  call Ar_Bl_BL_EXT_AR_NEW(LIN,INTPOS,ISMA,1)
end do

return

end subroutine D1V_Ar_ACT_BlBL_SGT0

!=======================================================
! END OF DV, NEXT VD
!=======================================================

subroutine SS_DRL_ACT_BR_SGT0(LIN,LRA)
! SS(1-16) (11)-Drl(22)-
! SS(1-19) Drl(12)-C"(21)-
! SS(1-20) Drl(33)-C"(11)-C"(22)-
! SS(1-20) (11)Drl(33)-C"(22)-
! SS(1-20) (11)(22)Drl(33)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_ss, w1_ss
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lmij, lmj, lri, lrj, lrk, mpl
integer(kind=iwp), external :: iwalk_ad

if (JML /= JMR) return
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    if (LMIJ /= JML) cycle
    ! SS(1-16) (11)-Drl(22)-
    IWDL = JUST(LRJ,LRI)
    IWDR = IWDL
    do MPL=1,MHLP
      IWAL = LPNEW_LWEI(MPL)
      IWAR = LPNEW_RWEI(MPL)
      LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
      LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
    end do
    do MPL=1,MTYPE
      VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0_SS(16)
      VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_SS(16)
    end do
    call Drl_BR_EXT_AL_NEW(LIN,LRJ,LRA)
    ! SS(1-18) Drl(11)-C"(22)-
    do MPL=1,MTYPE
      VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0_SS(18)
      VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_SS(18)
    end do
    call Drl_BR_EXT_AL_NEW(LIN,LRI,LRA)
    ! SS(1-20) Drl(33)-C"(11)-C"(22)-
    ! SS(1-20) (11)Drl(33)-C"(22)-
    ! SS(1-20) (11)(22)Drl(33)-
    do MPL=1,MTYPE
      VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0_SS(20)
      VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_SS(20)
    end do
    if (LRA > NORB_DZ) then
      call Drl_BR_SUM_AL_new(LIN,LRI,LRJ,LRA)
    else
      do LRK=1,norb_dz
        if (LRK == LRI) cycle
        if (LRK == LRJ) cycle
        call Drl_BR_EXT_AL_NEW(LIN,LRK,LRA)
      end do
    end if
  end do
end do

return

end subroutine SS_DRL_ACT_BR_SGT0

subroutine SS_S_DRL_ACT_BR_SGT0(LIN,LRA)
! SS(1-19) Drl(12)-C"(21)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_ss, w1_ss
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lmij, lmj, lri, lrj, mpl
integer(kind=iwp), external :: iwalk_ad

if (JML /= JMR) return
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    if (LMIJ /= JML) cycle
    ! SS(1-19) Drl(12)-C"(21)-
    IWDL = JUST(LRJ,LRI)
    IWDR = JUST(LRI,LRJ)
    do MPL=1,MHLP
      IWAL = LPNEW_LWEI(MPL)
      IWAR = LPNEW_RWEI(MPL)
      LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
      LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
    end do
    do MPL=1,MTYPE
      VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0_SS(19)
      VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_SS(19)
    end do
    call Drl_BR_EXT_AL_NEW(LIN,LRI,LRA)
  end do
end do

return

end subroutine SS_S_DRL_ACT_BR_SGT0

subroutine SS_S_DRL_ACT_BL_SGT0(LIN,LRA)
! SS(1-19) Drl(12)-C"(21)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_ss, w1_ss
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lmij, lmj, lri, lrj, mpl
integer(kind=iwp), external :: iwalk_ad

if (JML /= JMR) return
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    if (LMIJ /= JML) cycle
    ! SS(1-19) Drl(12)-C"(21)-
    IWDL = JUST(LRJ,LRI)
    IWDR = JUST(LRI,LRJ)
    do MPL=1,MHLP
      IWAL = LPNEW_LWEI(MPL)
      IWAR = LPNEW_RWEI(MPL)
      LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
      LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
    end do
    do MPL=1,MTYPE
      VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0_SS(19)
      VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_SS(19)
    end do
    call Drl_BL_EXT_AR_NEW(LIN,LRI,LRA)
  end do
end do

return

end subroutine SS_S_DRL_ACT_BL_SGT0

subroutine SS_ARBL_ACT_BR_SGT0(LIN,LRA)
!=======================================================================
! SS(1-1)  Ar(01)-Bl(32)-
! SS(1-3)  Ar(13)-Bl(20)-
! SS(1-6)  (11)-Ar(23)-Bl(32)-
! SS(1-7)  Ar(13)-C'(21)-Bl(32)-
! SS(1-8)  Ar(13)-C'(22)-Bl(31)-
! SS(1-9)  Ar(23)-C'(11)-Bl(32)-
! SS(1-11) Ar(13)-Bl(31)-C"(22)-
! SS(1-12) Ar(13)-Bl(32)-C"(21)-
! SS(1-13) Ar(23)-Bl(31)-C"(12)-
!=======================================================================

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, lsm_inn, mhlp, mtype, ngw2, ngw3, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, &
                         vplpnew_w1, w0_ss, w1_ss
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, jmlr, lmi, lmij, lmj, lmk, lri, lrj, lrk, mpl, ni
real(kind=wp) :: w0ss1, w0ss11, w0ss12, w0ss13, w0ss3, w0ss6, w0ss7, w0ss8, w0ss9, w1ss1, w1ss11, w1ss12, w1ss13, w1ss3, w1ss6, &
                 w1ss7, w1ss8, w1ss9
integer(kind=iwp), external :: iwalk_ad

JMLR = Mul(JML,JMR)
ISMA = Mul(IML,IMR)
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    if (LMIJ /= JMLR) cycle
    IJK = LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ)
    intpos = INTIND_IJKA(IJK)
    !-------------------------------------------------------------------
    W0SS1 = W0_SS(1)
    W1SS1 = W1_SS(1)
    W0SS3 = W0_SS(3)
    W1SS3 = W1_SS(3)
    W0SS6 = W0_SS(6)
    W1SS6 = W1_SS(6)
    W0SS7 = W0_SS(7)
    W1SS7 = W1_SS(7)
    W0SS8 = W0_SS(8)
    W1SS8 = W1_SS(8)
    W0SS9 = W0_SS(9)
    W1SS9 = W1_SS(9)
    W0SS11 = W0_SS(11)
    W1SS11 = W1_SS(11)
    W0SS12 = W0_SS(12)
    W1SS12 = W1_SS(12)
    W0SS13 = W0_SS(13)
    W1SS13 = W1_SS(13)
    NI = mod(LRJ-LRI,2)
    if (NI == 0) then
      W0SS1 = -W0SS1
      W1SS1 = -W1SS1
      W0SS3 = -W0SS3
      W1SS3 = -W1SS3
      W0SS6 = -W0SS6
      W1SS6 = -W1SS6
      W0SS7 = -W0SS7
      W1SS7 = -W1SS7
      W0SS8 = -W0SS8
      W1SS8 = -W1SS8
      W0SS9 = -W0SS9
      W1SS9 = -W1SS9
      W0SS11 = -W0SS11
      W1SS11 = -W1SS11
      W0SS12 = -W0SS12
      W1SS12 = -W1SS12
      W0SS13 = -W0SS13
      W1SS13 = -W1SS13
    end if
    !-------------------------------------------------------------------
    ! SS(1-1)  Ar(01)-Bl(32)-
    if (JML == 1) then
      IWDL = JUST(LRI,LRI)
      IWDR = JUST(LRJ,LRI)
      do MPL=1,MTYPE
        VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS1
        VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SS1
      end do
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      call Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
    end if
    ! SS(1-3)  Ar(13)-Bl(20)-
    if (JMR == 1) then
      IWDL = JUST(LRJ,LRI)
      IWDR = JUST(LRJ,LRJ)
      do MPL=1,MTYPE
        VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS3
        VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SS3
      end do
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      call Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
    end if
    !-------------------------------------------------------------------
    ! SS(1-6)  (11)-Ar(23)-Bl(32)-
    do LRK=NORB_FRZ+1,LRI-1
      LMK = LSM_INN(LRK)
      if (Mul(LMK,LMI) /= JML) cycle
      IWDL = JUST(LRI,LRK)
      IWDR = JUST(LRJ,LRK)
      do MPL=1,MTYPE
        VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS6
        VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SS6
      end do
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      call Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
    end do
    !-------------------------------------------------------------------
    do LRK=LRI+1,LRJ-1
      LMK = LSM_INN(LRK)
      if (Mul(LMK,LMI) /= JML) cycle
      ! SS(1-7)  Ar(13)-C'(21)-Bl(32)-
      IWDL = JUST(LRK,LRI)
      IWDR = JUST(LRJ,LRK)
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      do MPL=1,MTYPE
        VPLP_W0(MPL) = -VPLPNEW_W0(MPL)*W0SS7
        VPLP_W1(MPL) = -VPLPNEW_W1(MPL)*W1SS7
      end do
      call Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
      ! SS(1-8)  Ar(13)-C'(22)-Bl(31)-
      IWDL = JUST(LRK,LRI)
      IWDR = JUST(LRK,LRJ)
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      do MPL=1,MTYPE
        VPLP_W0(MPL) = -VPLPNEW_W0(MPL)*W0SS8
        VPLP_W1(MPL) = -VPLPNEW_W1(MPL)*W1SS8
      end do
      call Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
      ! SS(1-9)  Ar(23)-C'(11)-Bl(32)-
      IWDL = JUST(LRI,LRK)
      IWDR = JUST(LRJ,LRK)
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      do MPL=1,MTYPE
        VPLP_W0(MPL) = -VPLPNEW_W0(MPL)*W0SS9
        VPLP_W1(MPL) = -VPLPNEW_W1(MPL)*W1SS9
      end do
      call Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
    end do
    !-------------------------------------------------------------------
    do LRK=LRJ+1,NORB_DZ
      LMK = LSM_INN(LRK)
      if (Mul(LMK,LMI) /= JML) cycle
      ! SS(1-11) Ar(13)-Bl(31)-C"(22)-
      IWDL = JUST(LRK,LRI)
      IWDR = JUST(LRK,LRJ)
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      do MPL=1,MTYPE
        VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS11
        VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SS11
      end do
      call Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
      ! SS(1-12) Ar(13)-Bl(32)-C"(21)-
      IWDL = JUST(LRK,LRI)
      IWDR = JUST(LRJ,LRK)
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      do MPL=1,MTYPE
        VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS12
        VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SS12
      end do
      call Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
      ! SS(1-13) Ar(23)-Bl(31)-C"(12)-
      IWDL = JUST(LRI,LRK)
      IWDR = JUST(LRK,LRJ)
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      do MPL=1,MTYPE
        VPLP_W0(MPL) = VPLPNEW_W0(MPL)*W0SS13
        VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1SS13
      end do
      call Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
    end do
    !-------------------------------------------------------------------
  end do
end do

return

end subroutine SS_ARBL_ACT_BR_SGT0

subroutine ST_ARBL_ACT_BR_SGT0(LIN,LRA)
! ST(2-3) Ar(13)-C'(22)-Bl(32)-
! ST(2-3) Ar(13)-Bl(32)-C'(22)-

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, lsm_inn, mhlp, mtype, ngw2, ngw3, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w1, w1_st
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, jmlr, lmi, lmij, lmj, lmk, lri, lrj, lrk, mpl, ni
real(kind=wp) :: w1st3
integer(kind=iwp), external :: iwalk_ad

JMLR = Mul(JML,JMR)
ISMA = Mul(IML,IMR)
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    if (LMIJ /= JMLR) cycle

    W1ST3 = W1_ST(3)
    NI = mod(LRJ-LRI,2)
    if (NI == 0) then
      W1ST3 = -W1ST3
    end if
    IJK = LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ)
    intpos = INTIND_IJKA(IJK)
    !-------------------------------------------------------------------
    ! ST(2-3) Ar(13)-C'(22)-Bl(32)-
    do LRK=LRI+1,LRJ-1
      LMK = LSM_INN(LRK)
      if (Mul(LMK,LMI) /= JML) cycle
      IWDL = JUST(LRK,LRI)
      IWDR = JUST(LRK,LRJ)
      do MPL=1,MTYPE
        VPLP_W0(MPL) = Zero
        VPLP_W1(MPL) = -VPLPNEW_W1(MPL)*W1ST3
      end do
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      call Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
    end do
    !-------------------------------------------------------------------
    ! ST(2-3) Ar(13)-Bl(32)-C'(22)-
    do LRK=LRJ+1,NORB_DZ
      LMK = LSM_INN(LRK)
      if (Mul(LMK,LMI) /= JML) cycle
      IWDL = JUST(LRK,LRI)
      IWDR = JUST(LRJ,LRK)
      do MPL=1,MTYPE
        VPLP_W0(MPL) = Zero
        VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1ST3
      end do
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      call Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
    end do
  end do
end do

return

end subroutine ST_ARBL_ACT_BR_SGT0

subroutine ST_DRL_ACT_BR_SGT0(LIN,LRA)
! ST(2-7) Drl(12)-C"(22)-

use gugaci_global, only: ipae, ipael, jml, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w1, w1_st
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lmij, lmj, lri, lrj, mpl
integer(kind=iwp), external :: iwalk_ad

do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    if (LMIJ /= JML) cycle
    !-------------------------------------------------------------------
    ! ST(2-7) Drl(12)-C"(22)-
    IWDL = JUST(LRJ,LRI)
    IWDR = JUST(LRI,LRJ)
    do MPL=1,MHLP
      IWAL = LPNEW_LWEI(MPL)
      IWAR = LPNEW_RWEI(MPL)
      LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
      LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
    end do
    do MPL=1,MTYPE
      VPLP_W0(MPL) = Zero
      VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_ST(7)
    end do
    call Drl_BR_EXT_AL_NEW(LIN,LRI,LRA)
  end do
end do

return

end subroutine ST_DRL_ACT_BR_SGT0

subroutine TS_ARBL_ACT_BR_SGT0(LIN,LRA)
!=======================================================================
! TS(3) A&R-B^L-  ACT -B&R ............................................
! TS(3-3) Ar(23)-Bl(31)-C"(22)-

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, lsm_inn, mhlp, mtype, ngw2, ngw3, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w1, w1_ts
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, lmi, lmik, lmj, lmjk, lmk, lri, lrj, lrk, mpl, ni
real(kind=wp) :: w1ts3
integer(kind=iwp), external :: iwalk_ad

ISMA = Mul(IML,IMR)
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    IJK = LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ)
    intpos = INTIND_IJKA(IJK)
    !-------------------------------------------------------------------
    W1TS3 = W1_TS(3)
    NI = mod(LRJ-LRI,2)
    if (NI == 0) then
      W1TS3 = -W1TS3
    end if
    !-------------------------------------------------------------------
    ! TS(3-3) Ar(23)-Bl(31)-C"(22)-
    do MPL=1,MTYPE
      VPLP_W0(MPL) = Zero
      VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1TS3
    end do
    do LRK=LRJ+1,NORB_DZ
      LMK = LSM_INN(LRK)
      LMIK = Mul(LMI,LMK)
      LMJK = Mul(LMJ,LMK)
      if ((LMIK /= JML) .or. (LMJK /= JMR)) cycle
      IWDL = JUST(LRI,LRK)
      IWDR = JUST(LRK,LRJ)
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      call Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
    end do
    !-------------------------------------------------------------------
  end do
end do

return

end subroutine TS_ARBL_ACT_BR_SGT0

subroutine STT_ARBL_ACT_BR_SGT1(LIN,LRA)
!=======================================================================
! STT(4) A&R-B^L-  ACT -BR ............................................
! ST1(4-1) Ar(01)-Bl(31)-
! ST1(4-2) Ar(23)-Bl(31)-
! ST1(4-3) Ar(13)-C'(21)-Bl(31)-
! ST1(4-3) Ar(13)-Bl(31)-C"(21)-
! ST1(4-4) Ar(23)-C'(11)-Bl(31)-
! ST1(4-4) Ar(23)-Bl(31)-C"(11)-

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, lsm_inn, mhlp, mtype, ngw2, ngw3, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w1, w1_st1
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, lmi, lmij, lmj, lmk, lri, lrj, lrk, mpl
real(kind=wp) :: w1st1, w1st2, w1st3, w1st4
integer(kind=iwp), external :: iwalk_ad

ISMA = Mul(IML,IMR)
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    W1ST1 = W1_ST1(1)
    W1ST2 = W1_ST1(2)
    W1ST3 = W1_ST1(3)
    W1ST4 = W1_ST1(4)
    if (mod(LRJ-LRI,2) == 0) then
      W1ST1 = -W1ST1
      W1ST2 = -W1ST2
      W1ST3 = -W1ST3
      W1ST4 = -W1ST4
    end if
    IJK = LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)+NGW3(LRA-NORB_FRZ) !???
    intpos = INTIND_IJKA(IJK)                         !???
    ! ST1(4-1) Ar(01)-Bl(31)-
    ! ST1(4-2) Ar(23)-Bl(31)-
    if ((JML == 1) .and. (JMR == LMIJ)) then
      IWDL = JUST(LRI,LRI)
      IWDR = JUST(LRI,LRJ)
      do MPL=1,MHLP
        IWAL = LPNEW_LWEI(MPL)
        IWAR = LPNEW_RWEI(MPL)
        LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      end do
      do MPL=1,MTYPE
        VPLP_W0(MPL) = Zero
        VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1ST1
      end do
      call Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
    end if
    ! ST1(4-2) (11)Ar(23)-Bl(31)-
    do LRK=NORB_FRZ+1,LRI-1
      LMK = LSM_INN(LRK)
      if ((JML == Mul(LMK,LMI)) .and. (JMR == Mul(LMK,LMJ))) then
        IWDL = JUST(LRI,LRK)
        IWDR = JUST(LRK,LRJ)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = Zero
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1ST2
        end do
        call Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
      end if
    end do
    ! ST1(4-3) Ar(13)-C'(21)-Bl(31)-
    do LRK=LRI+1,LRJ-1
      LMK = LSM_INN(LRK)
      if ((JML == Mul(LMI,LMK)) .and. (JMR == Mul(LMK,LMJ))) then
        IWDL = JUST(LRK,LRI)
        IWDR = JUST(LRK,LRJ)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = Zero
          VPLP_W1(MPL) = -VPLPNEW_W1(MPL)*W1ST3
        end do
        call Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
        ! ST1(4-4) Ar(23)-C'(11)-Bl(31)-
        IWDL = JUST(LRI,LRK)
        IWDR = JUST(LRK,LRJ)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = Zero
          VPLP_W1(MPL) = -VPLPNEW_W1(MPL)*W1ST4
        end do
        call Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
      end if
    end do
    ! ST1(4-3) Ar(13)-Bl(31)-C"(21)-
    do LRK=LRJ+1,NORB_DZ
      LMK = LSM_INN(LRK)
      if ((JML == Mul(LMI,LMK)) .and. (JMR == Mul(LMJ,LMK))) then
        IWDL = JUST(LRK,LRI)
        IWDR = JUST(LRJ,LRK)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = Zero
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1ST3
        end do
        call Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
        ! ST1(4-4) Ar(23)-Bl(31)-C"(11)-
        IWDL = JUST(LRI,LRK)
        IWDR = JUST(LRJ,LRK)
        do MPL=1,MHLP
          IWAL = LPNEW_LWEI(MPL)
          IWAR = LPNEW_RWEI(MPL)
          LP_LWEI(MPL) = IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
          LP_RWEI(MPL) = IWALK_AD(JPAD,Ipae,IWAR,IWDR)
        end do
        do MPL=1,MTYPE
          VPLP_W0(MPL) = Zero
          VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1ST4
        end do
        call Ar_Bl_BR_EXT_AL_NEW(LIN,INTPOS,ISMA,1)
      end if
    end do
  end do
end do

return

end subroutine STT_ARBL_ACT_BR_SGT1

subroutine TTS_Drl_ACT_BR_SGT1(LIN,LRA)
! T1S(5-5)   (11)Drl(12)-
! T1S(5-6)   Drl(12)-C"(12)-
! T1S(5-7)   Drl(12)-C"(11)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w1, w1_t1s
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lmj, lri, lrj, mpl
integer(kind=iwp), external :: iwalk_ad

do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    if ((JML /= Mul(LMI,LMJ)) .or. (JMR /= Mul(LMI,LMJ))) cycle
    ! T1S(5-5)   (11)Drl(12)-
    IWDL = JUST(LRI,LRJ)
    IWDR = JUST(LRJ,LRI)
    do MPL=1,MHLP
      IWAL = LPNEW_LWEI(MPL)
      IWAR = LPNEW_RWEI(MPL)
      LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
      LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
    end do
    do MPL=1,MTYPE
      VPLP_W0(MPL) = Zero
      VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_T1S(5)
    end do
    call Drl_BR_EXT_AL_NEW(LIN,LRJ,LRA)
    ! T1S(5-6)   Drl(11)-C"(12)-
    IWDL = JUST(LRI,LRJ)
    IWDR = JUST(LRJ,LRI)
    do MPL=1,MHLP
      IWAL = LPNEW_LWEI(MPL)
      IWAR = LPNEW_RWEI(MPL)
      LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
      LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
    end do
    do MPL=1,MTYPE
      VPLP_W0(MPL) = Zero
      VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_T1S(6)
    end do
    call Drl_BR_EXT_AL_NEW(LIN,LRI,LRA)
    ! T1S(5-7)   Drl(12)-C"(11)-
    IWDL = JUST(LRI,LRJ)
    IWDR = IWDL
    do MPL=1,MHLP
      IWAL = LPNEW_LWEI(MPL)
      IWAR = LPNEW_RWEI(MPL)
      LP_LWEI(MPL) = IWALK_AD(JPADL,IPAEL,IWAL,IWDL)
      LP_RWEI(MPL) = IWALK_AD(JPAD,IPAE,IWAR,IWDR)
    end do
    do MPL=1,MTYPE
      VPLP_W0(MPL) = Zero
      VPLP_W1(MPL) = VPLPNEW_W1(MPL)*W1_T1S(7)
    end do
    call Drl_BR_EXT_AL_NEW(LIN,LRI,LRA)
  end do
end do

return

end subroutine TTS_Drl_ACT_BR_SGT1

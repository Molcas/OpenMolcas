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

!=======================================================================
! have computed vlop0 and vlop1, link to extern space
! (ss, st, ts, tt)
! lin type of the act-ext space node
! lri index of ar
! lrj index of bl
!=======================================================================
subroutine arbl_act_c_link_ext_ab(lin,lri,lrj)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lri, lrj

if (lin == 1) call ar_bl_ext_ss(lri,lrj,1)
if (lin == 2) call ar_bl_ext_st(lri,lrj,1)
if (lin == 3) call ar_bl_ext_ts(lri,lrj,1)
if (lin == 11) call ar_bl_ext_tt(lri,lrj,1)
if (lin == 10) call ar_br_sv_ext_br_ar(lri,lrj)
if (lin == 17) call ar_br_tv_ext_br_ar(lri,lrj)

return

end subroutine arbl_act_c_link_ext_ab

!=======================================================================
! have computed vlop0 and vlop1, link to extern space
! (ss, st, ts, tt)
! lin type of the act-ext space node
! lri index of drl
!=======================================================================
subroutine drl_act_c_link_ext_ab(lin,lri)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lri

if (lin == 1) call drl_ss_ext(lri)
if (lin == 2) call drl_st_ext(lri)
if (lin == 3) call drl_ts_ext(lri)
if (lin == 11) call drl_tt_ext(lri)

return

end subroutine drl_act_c_link_ext_ab

subroutine drl_act_c_link_ext_ab_sum(lin,lri,lrj)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lri, lrj

if (lin == 1) call drl_ss_sum(lri,lrj)
if (lin == 11) call drl_tt_sum(lri,lrj)

return

end subroutine drl_act_c_link_ext_ab_sum

subroutine sd_ar_act_bl_sgt0(lin,lra)
! sd(6-3) a&r(13)c'(22)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, &
                         mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_sd
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, jmlr, lmi, lmk, lri, lrk, mpl, ni
real(kind=wp) :: w0sd3
integer(kind=iwp), external :: iwalk_ad

jmlr = Mul(jml,jmr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (lmi /= jmlr) cycle
  w0sd3 = w0_sd(3)
  ni = mod(norb_dz-lri,2)
  if (ni == 1) then
    w0sd3 = -w0sd3
  end if
  ! sd(6-3) a&r(13)c'(22)-
  do mpl=1,mtype
    vplp_w0(mpl) = -vplpnew_w0(mpl)*w0sd3
    vplp_w1(mpl) = -vplpnew_w1(mpl)*w0sd3
  end do
  do lrk=lri+1,norb_dz
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lrk,lri)
    iwdr = jud(lrk)
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    if (lin == 1) call ar_bl_ext_ss(lri,lra,1)
    if (lin == 2) call ar_bl_ext_st(lri,lra,1)
    if (lin == 3) call ar_bl_ext_ts(lri,lra,1)
    if (lin == 11) call ar_bl_ext_tt(lri,lra,1)
  end do
end do

return

end subroutine sd_ar_act_bl_sgt0

!=======================================================================
! extern space -bl-ar or -br-al (ss,st,ts,tt)
! double occupied space ar-bl-  (ss)
! active space -c"-
! spin great than 0
!=======================================================================
subroutine ss_arbl_act_c_ext_ab_sgt0(lin)
!-----------------------------------------------------------------------
! ss(1-1)  ar(01)-bl(32)-        act -c"-
! ss(1-3)  ar(13)-bl(20)-        act -c"-
! ss(1-6)  (11)-ar(23)-bl(32)-   act -c"-
! ss(1-7)  ar(13)-c'(21)-bl(32)- act -c"-
! ss(1-8)  ar(13)-c'(22)-bl(31)- act -c"-
! ss(1-9)  ar(23)-c'(11)-bl(32)- act -c"-
! ss(1-11) ar(13)-bl(31)-c"(22)- act -c"-
! ss(1-12) ar(13)-bl(32)-c"(21)- act -c"-
! ss(1-13) ar(23)-bl(31)-c"(12)- act -c"-
!-----------------------------------------------------------------------

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_ss, w1_ss
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lmij, lmj, lmk, lmki, lmkj, lri, lrj, lrk, mpl, ni
real(kind=wp) :: w0ss1, w0ss11, w0ss12, w0ss13, w0ss3, w0ss6, w0ss7, w0ss8, w0ss9, w1ss1, w1ss11, w1ss12, w1ss13, w1ss3, w1ss6, &
                 w1ss7, w1ss8, w1ss9
integer(kind=iwp), external :: iwalk_ad

do lri=norb_frz+1,norb_dz-1
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    w0ss1 = w0_ss(1)
    w1ss1 = w1_ss(1)
    w0ss3 = w0_ss(3)
    w1ss3 = w1_ss(3)
    w0ss6 = w0_ss(6)
    w1ss6 = w1_ss(6)
    w0ss7 = w0_ss(7)
    w1ss7 = w1_ss(7)
    w0ss8 = w0_ss(8)
    w1ss8 = w1_ss(8)
    w0ss9 = w0_ss(9)
    w1ss9 = w1_ss(9)
    w0ss11 = w0_ss(11)
    w1ss11 = w1_ss(11)
    w0ss12 = w0_ss(12)
    w1ss12 = w1_ss(12)
    w0ss13 = w0_ss(13)
    w1ss13 = w1_ss(13)
    ni = mod(lrj-lri,2)
    if (ni == 0) then
      w0ss1 = -w0ss1
      w1ss1 = -w1ss1
      w0ss3 = -w0ss3
      w1ss3 = -w1ss3
      w0ss6 = -w0ss6
      w1ss6 = -w1ss6
      w0ss7 = -w0ss7
      w1ss7 = -w1ss7
      w0ss8 = -w0ss8
      w1ss8 = -w1ss8
      w0ss9 = -w0ss9
      w1ss9 = -w1ss9
      w0ss11 = -w0ss11
      w1ss11 = -w1ss11
      w0ss12 = -w0ss12
      w1ss12 = -w1ss12
      w0ss13 = -w0ss13
      w1ss13 = -w1ss13
    end if
    if ((jml == 1) .and. (lmij == jmr)) then
      ! ss(1-1)  ar(01)-bl(32)-        act -c"-
      iwdl = just(lri,lri)
      iwdr = just(lrj,lri)
      do mpl=1,mhlp
        iwal = lpnew_lwei(mpl)
        iwar = lpnew_rwei(mpl)
        lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
        lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
      end do
      do mpl=1,mtype
        vplp_w0(mpl) = vplpnew_w0(mpl)*w0ss1
        vplp_w1(mpl) = vplpnew_w1(mpl)*w1ss1
      end do
      call arbl_act_c_link_ext_ab(lin,lri,lrj)
    end if
    if ((jmr == 1) .and. (lmij == jml)) then
      ! ss(1-3)  ar(13)-bl(20)-        act -c"-
      iwdl = just(lrj,lri)
      iwdr = just(lrj,lrj)
      do mpl=1,mhlp
        iwal = lpnew_lwei(mpl)
        iwar = lpnew_rwei(mpl)
        lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
        lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
      end do
      do mpl=1,mtype
        vplp_w0(mpl) = vplpnew_w0(mpl)*w0ss3
        vplp_w1(mpl) = vplpnew_w1(mpl)*w1ss3
      end do
      call arbl_act_c_link_ext_ab(lin,lri,lrj)
    end if
    ! ss(1-6)  (11)-ar(23)-bl(32)-   act -c"-
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0ss6
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1ss6
    end do
    do lrk=norb_frz+1,lri-1
      lmk = lsm_inn(lrk)
      lmki = Mul(lmi,lmk)
      lmkj = Mul(lmj,lmk)
      if ((lmki == jml) .and. (lmkj == jmr)) then
        iwdl = just(lri,lrk)
        iwdr = just(lrj,lrk)
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        call arbl_act_c_link_ext_ab(lin,lri,lrj)
      end if
    end do
    ! ss(1-7)  ar(13)-c'(21)-bl(32)- act -c"-
    do mpl=1,mtype
      vplp_w0(mpl) = -vplpnew_w0(mpl)*w0ss7
      vplp_w1(mpl) = -vplpnew_w1(mpl)*w1ss7
    end do
    do lrk=lri+1,lrj-1
      lmk = lsm_inn(lrk)
      lmki = Mul(lmk,lmi)
      lmkj = Mul(lmj,lmk)
      if ((lmki == jml) .and. (lmkj == jmr)) then
        iwdl = just(lrk,lri)
        iwdr = just(lrj,lrk)
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        call arbl_act_c_link_ext_ab(lin,lri,lrj)
      end if
    end do
    ! ss(1-8)  ar(13)-c'(22)-bl(31)- act -c"-
    do mpl=1,mtype
      vplp_w0(mpl) = -vplpnew_w0(mpl)*w0ss8
      vplp_w1(mpl) = -vplpnew_w1(mpl)*w1ss8
    end do
    do lrk=lri+1,lrj-1
      lmk = lsm_inn(lrk)
      lmki = Mul(lmk,lmi)
      lmkj = Mul(lmk,lmj)
      if ((lmki == jml) .and. (lmkj == jmr)) then
        iwdl = just(lrk,lri)
        iwdr = just(lrk,lrj)
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        call arbl_act_c_link_ext_ab(lin,lri,lrj)
      end if
    end do
    ! ss(1-9)  ar(23)-c'(11)-bl(32)- act -c"-
    do mpl=1,mtype
      vplp_w0(mpl) = -vplpnew_w0(mpl)*w0ss9
      vplp_w1(mpl) = -vplpnew_w1(mpl)*w1ss9
    end do
    do lrk=lri+1,lrj-1
      lmk = lsm_inn(lrk)
      lmki = Mul(lmk,lmi)
      lmkj = Mul(lmk,lmj)
      if ((lmki == jml) .and. (lmkj == jmr)) then
        iwdl = just(lri,lrk)
        iwdr = just(lrj,lrk)
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        call arbl_act_c_link_ext_ab(lin,lri,lrj)
      end if
    end do
    ! ss(1-11) ar(13)-bl(31)-c"(22)- act -c"-
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0ss11
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1ss11
    end do
    do lrk=lrj+1,norb_dz
      lmk = lsm_inn(lrk)
      lmki = Mul(lmk,lmi)
      lmkj = Mul(lmk,lmj)
      if ((lmki == jml) .and. (lmkj == jmr)) then
        iwdl = just(lrk,lri)
        iwdr = just(lrk,lrj)
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        call arbl_act_c_link_ext_ab(lin,lri,lrj)
      end if
    end do
    ! ss(1-12) ar(13)-bl(32)-c"(21)- act -c"-
    do mpl=1,mtype
      vplp_w0(mpl) = Zero
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1ss12
    end do
    do lrk=lrj+1,norb_dz
      lmk = lsm_inn(lrk)
      lmki = Mul(lmk,lmi)
      lmkj = Mul(lmj,lmk)
      if ((lmki == jml) .and. (lmkj == jmr)) then
        iwdl = just(lrk,lri)
        iwdr = just(lrj,lrk)
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        call arbl_act_c_link_ext_ab(lin,lri,lrj)
      end if
    end do
    ! ss(1-13) ar(23)-bl(31)-c"(12)- act -c"-
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0ss13
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1ss13
    end do
    do lrk=lrj+1,norb_dz
      lmk = lsm_inn(lrk)
      lmki = Mul(lmi,lmk)
      lmkj = Mul(lmk,lmj)
      if ((lmki == jml) .and. (lmkj == jmr)) then
        iwdl = just(lri,lrk)
        iwdr = just(lrk,lrj)
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        call arbl_act_c_link_ext_ab(lin,lri,lrj)
      end if
    end do
  end do
end do

return

end subroutine ss_arbl_act_c_ext_ab_sgt0

!=======================================================================
! extern space -bl-ar or -br-al (ss,st,ts,tt)
! double occupied space drl-  (ss)
! active space -c"-
! spin great than 0
!=======================================================================
subroutine ss_drl_act_c_ext_ab_sgt0(lin)
! ss(1-16) (11)-drl(22)-         act -c"-
! ss(1-18) drl(11)-c"(22)-       act -c"-
! ss(1-19) drl(12)-c"(21)-       act -c"-
! ss(1-20) (11)-drl(33)-c"(22)-  act -c"-
! ss(1-20) drl(33)-c"(11)-c"(22)-act -c"-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_ss, w1_ss
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lmij, lmj, lri, lrj, lrk, mpl
integer(kind=iwp), external :: iwalk_ad

do lri=norb_frz+1,norb_dz-1
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    if ((lmij /= jml) .or. (lmij /= jmr)) cycle
    iwdl = just(lrj,lri)
    iwdr = iwdl
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    ! ss(1-16) (11)-drl(22)-         act -c"-
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0_ss(16)
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1_ss(16)
    end do
    call drl_act_c_link_ext_ab(lin,lrj)
    ! ss(1-18) drl(11)-c"(22)-       act -c"-
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0_ss(18)
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1_ss(18)
    end do
    call drl_act_c_link_ext_ab(lin,lri)
    ! ss(1-20) (11)-drl(33)-c"(22)-  act -c"-
    ! ss(1-20) drl(33)-c"(11)-c"(22)-act -c"-
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0_ss(20)
      vplp_w1(mpl) = Zero
    end do
    do lrk=1,norb_dz
      if (lrk == lri) cycle
      if (lrk == lrj) cycle
      call drl_act_c_link_ext_ab(lin,lrk)
    end do
  end do
end do

return

end subroutine ss_drl_act_c_ext_ab_sgt0

!=======================================================================
! extern space -bl-ar or -br-al (ss,st,ts,tt)
! double occupied space drl-  (st)
! active space -c"-
! spin great than 0
!=======================================================================
subroutine st_drl_act_c_ext_ab_sgt0(lin)
! st(2-7) drl(12)-c"(22)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w1, w1_st
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lmij, lmj, lri, lrj, mpl
integer(kind=iwp), external :: iwalk_ad

if (jmr /= jml) return
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    if (jml /= lmij) cycle
    iwdl = just(lrj,lri)
    iwdr = just(lri,lrj)
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    do mpl=1,mtype
      vplp_w0(mpl) = Zero
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1_st(7)
    end do
    call drl_act_c_link_ext_ab(lin,lri)
  end do
end do

return

end subroutine st_drl_act_c_ext_ab_sgt0

!=======================================================================
! extern space -bl-ar or -br-al (ss,st,ts,tt)
! double occupied space ar-bl-  (st)
! active space -c"-
! spin great than 0
!=======================================================================
subroutine st_arbl_act_c_ext_ab_sgt0(lin)
! st(2-3) ar(13)-c'(22)-bl(32)-
! st(2-3) ar(13)-bl(32)-c'(22)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w1, w1_st
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lmj, lmk, lmki, lmkj, lri, lrj, lrk, mpl, ni
real(kind=wp) :: w1st3
integer(kind=iwp), external :: iwalk_ad

do lri=norb_frz+1,norb_dz-1
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    w1st3 = w1_st(3)
    ni = mod(lrj-lri,2)
    if (ni == 0) then
      w1st3 = -w1st3
    end if
    do lrk=lri+1,lrj-1
      lmk = lsm_inn(lrk)
      lmki = Mul(lmk,lmi)
      lmkj = Mul(lmk,lmj)
      if ((lmki == jml) .and. (lmkj == jmr)) then
        ! st(2-3) ar(13)-c'(22)-bl(32)-
        iwdl = just(lrk,lri)
        iwdr = just(lrk,lrj)
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        do mpl=1,mtype
          vplp_w0(mpl) = Zero
          vplp_w1(mpl) = -vplpnew_w1(mpl)*w1st3
        end do
        call arbl_act_c_link_ext_ab(lin,lri,lrj)
      end if
    end do

    do lrk=lrj+1,norb_dz
      lmk = lsm_inn(lrk)
      lmki = Mul(lmk,lmi)
      lmkj = Mul(lmk,lmj)
      if ((lmki == jml) .and. (lmkj == jmr)) then
        ! st(2-3) ar(13)-bl(32)-c'(22)-
        iwdl = just(lrk,lri)
        iwdr = just(lrj,lrk)
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        do mpl=1,mtype
          vplp_w0(mpl) = Zero
          vplp_w1(mpl) = vplpnew_w1(mpl)*w1st3
        end do
        call arbl_act_c_link_ext_ab(lin,lri,lrj)
      end if
    end do
  end do
end do

return

end subroutine st_arbl_act_c_ext_ab_sgt0

subroutine stt_arbl_act_c_ext_ab_sgt1(lin)
! st1(4-1) ar(01)-bl(31)-
! st1(4-2) (11)ar(23)-bl(31)-
! st1(4-3) ar(13)-c'(21)-bl(31)-
! st1(4-3) ar(13)-bl(31)-c"(21)-
! st1(4-4) ar(23)-c'(11)-bl(31)-
! st1(4-4) ar(23)-bl(31)-c"(11)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w1, w1_st1
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lmj, lmk, lri, lrj, lrk, mpl
real(kind=wp) :: w1st1, w1st2, w1st3, w1st4
integer(kind=iwp), external :: iwalk_ad

do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    w1st1 = w1_st1(1)
    w1st2 = w1_st1(2)
    w1st3 = w1_st1(3)
    w1st4 = w1_st1(4)
    if (mod(lrj-lri,2) == 0) then
      w1st1 = -w1st1
      w1st2 = -w1st2
      w1st3 = -w1st3
      w1st4 = -w1st4
    end if
    ! st1(4-1) ar(01)-bl(31)-
    if ((jml == 1) .and. (jmr == Mul(lmi,lmj))) then
      iwdl = just(lri,lri)
      iwdr = just(lri,lrj)
      do mpl=1,mhlp
        iwal = lpnew_lwei(mpl)
        iwar = lpnew_rwei(mpl)
        lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
        lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
      end do
      do mpl=1,mtype
        vplp_w0(mpl) = Zero
        vplp_w1(mpl) = vplpnew_w1(mpl)*w1st1
      end do
      call arbl_act_c_link_ext_ab(lin,lri,lrj)
    end if
    ! st1(4-2) (11)ar(23)-bl(31)-
    do lrk=norb_frz+1,lri-1
      lmk = lsm_inn(lrk)
      if ((jml /= Mul(lmk,lmi)) .or. (jmr /= Mul(lmk,lmj))) cycle
      iwdl = just(lri,lrk)
      iwdr = just(lrk,lrj)
      do mpl=1,mhlp
        iwal = lpnew_lwei(mpl)
        iwar = lpnew_rwei(mpl)
        lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
        lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
      end do
      do mpl=1,mtype
        vplp_w0(mpl) = Zero
        vplp_w1(mpl) = vplpnew_w1(mpl)*w1st2
      end do
      call arbl_act_c_link_ext_ab(lin,lri,lrj)
    end do
    ! st1(4-3) ar(13)-c'(21)-bl(31)-
    do lrk=lri+1,lrj-1
      lmk = lsm_inn(lrk)
      if ((jml == Mul(lmi,lmk)) .and. (jmr == Mul(lmk,lmj))) then
        iwdl = just(lrk,lri)
        iwdr = just(lrk,lrj)
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        do mpl=1,mtype
          vplp_w0(mpl) = Zero
          vplp_w1(mpl) = -vplpnew_w1(mpl)*w1st3
        end do
        call arbl_act_c_link_ext_ab(lin,lri,lrj)
        ! st1(4-4) ar(23)-c'(11)-bl(31)-
        iwdl = just(lri,lrk)
        iwdr = just(lrk,lrj)
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        do mpl=1,mtype
          vplp_w0(mpl) = Zero
          vplp_w1(mpl) = -vplpnew_w1(mpl)*w1st4
        end do
        call arbl_act_c_link_ext_ab(lin,lri,lrj)
      end if
    end do
    do lrk=lrj+1,norb_dz
      lmk = lsm_inn(lrk)
      ! st1(4-4) ar(23)-bl(31)-c"(11)-
      if ((jml == Mul(lmi,lmk)) .and. (jmr == Mul(lmj,lmk))) then
        iwdl = just(lri,lrk)
        iwdr = just(lrj,lrk)
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        do mpl=1,mtype
          vplp_w0(mpl) = Zero
          vplp_w1(mpl) = vplpnew_w1(mpl)*w1st4
        end do
        call arbl_act_c_link_ext_ab(lin,lri,lrj)
        ! st1(4-3) ar(13)-bl(31)-c"(21)-
        iwdl = just(lrk,lri)
        iwdr = just(lrj,lrk)
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        do mpl=1,mtype
          vplp_w0(mpl) = Zero
          vplp_w1(mpl) = vplpnew_w1(mpl)*w1st3
        end do
        call arbl_act_c_link_ext_ab(lin,lri,lrj)
      end if
    end do
  end do
end do

return

end subroutine stt_arbl_act_c_ext_ab_sgt1

subroutine tts_arbl_act_c_ext_ab_sgt1(lin)
! t1s(5-1)   ar(13)-bl(10)-
! t1s(5-2)   ar(13)-bl(32)-
! t1s(5-2)   ar(13)-c'(11)-bl(32)-
! t1s(5-3)   ar(13)-bl(31)-c"(12)-
! t1s(5-4)   ar(13)-bl(32)-c"(11)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w1, w1_t1s
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lmj, lmk, lri, lrj, lrk, mpl
real(kind=wp) :: w1ts1, w1ts2, w1ts3, w1ts4
integer(kind=iwp), external :: iwalk_ad

do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    w1ts1 = w1_t1s(1)
    w1ts2 = w1_t1s(2)
    w1ts3 = w1_t1s(3)
    w1ts4 = w1_t1s(4)
    if (mod(lrj-lri,2) == 0) then
      w1ts1 = -w1ts1
      w1ts2 = -w1ts2
      w1ts3 = -w1ts3
      w1ts4 = -w1ts4
    end if
    ! t1s(5-1)   ar(13)-bl(10)-
    if ((jml == Mul(lmi,lmj)) .and. (jmr == 1)) then
      iwdl = just(lri,lrj)
      iwdr = just(lrj,lrj)
      do mpl=1,mhlp
        iwal = lpnew_lwei(mpl)
        iwar = lpnew_rwei(mpl)
        lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
        lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
      end do
      do mpl=1,mtype
        vplp_w0(mpl) = Zero
        vplp_w1(mpl) = vplpnew_w1(mpl)*w1ts1
      end do
      call arbl_act_c_link_ext_ab(lin,lri,lrj)
    end if
    ! t1s(5-2)   (11)ar(13)-bl(32)-
    do mpl=1,mtype
      vplp_w0(mpl) = Zero
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1ts2
    end do
    do lrk=norb_frz+1,lri-1
      lmk = lsm_inn(lrk)
      if ((jml /= Mul(lmk,lmi)) .or. (jmr /= Mul(lmk,lmj))) cycle
      iwdl = just(lrk,lri)
      iwdr = just(lrj,lrk)
      do mpl=1,mhlp
        iwal = lpnew_lwei(mpl)
        iwar = lpnew_rwei(mpl)
        lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
        lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
      end do
      call arbl_act_c_link_ext_ab(lin,lri,lrj)
    end do
    ! t1s(5-2)   ar(13)-c'(11)-bl(32)-
    do mpl=1,mtype
      vplp_w0(mpl) = Zero
      vplp_w1(mpl) = -vplpnew_w1(mpl)*w1ts2
    end do
    do lrk=lri+1,lrj-1
      lmk = lsm_inn(lrk)
      if ((jml /= Mul(lmi,lmk)) .or. (jmr /= Mul(lmk,lmj))) cycle
      iwdl = just(lri,lrk)
      iwdr = just(lrj,lrk)
      do mpl=1,mhlp
        iwal = lpnew_lwei(mpl)
        iwar = lpnew_rwei(mpl)
        lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
        lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
      end do
      call arbl_act_c_link_ext_ab(lin,lri,lrj)
    end do
    ! t1s(5-3)   ar(13)-bl(31)-c"(12)-
    do mpl=1,mtype
      vplp_w0(mpl) = Zero
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1ts3
    end do
    do lrk=lrj+1,norb_dz
      lmk = lsm_inn(lrk)
      if ((jml /= Mul(lmi,lmk)) .or. (jmr /= Mul(lmj,lmk))) cycle
      iwdl = just(lri,lrk)
      iwdr = just(lrk,lrj)
      do mpl=1,mhlp
        iwal = lpnew_lwei(mpl)
        iwar = lpnew_rwei(mpl)
        lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
        lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
      end do
      call arbl_act_c_link_ext_ab(lin,lri,lrj)
    end do
    ! t1s(5-4)   ar(13)-bl(32)-c"(11)-
    do mpl=1,mtype
      vplp_w0(mpl) = Zero
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1ts4
    end do
    do lrk=lrj+1,norb_dz
      lmk = lsm_inn(lrk)
      if ((jml /= Mul(lmi,lmk)) .or. (jmr /= Mul(lmj,lmk))) cycle
      iwdl = just(lri,lrk)
      iwdr = just(lrj,lrk)
      do mpl=1,mhlp
        iwal = lpnew_lwei(mpl)
        iwar = lpnew_rwei(mpl)
        lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
        lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
      end do
      call arbl_act_c_link_ext_ab(lin,lri,lrj)
    end do
  end do
end do

return

end subroutine tts_arbl_act_c_ext_ab_sgt1

subroutine tts_drl_act_c_ext_ab_sgt1(lin)
! t1s(5-5)   (11)drl(12)-
! t1s(5-6)   drl(11)-c"(12)-
! t1s(5-7)   drl(12)-c"(11)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w1, w1_t1s
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lmij, lmj, lri, lrj, mpl
integer(kind=iwp), external :: iwalk_ad

if (jml /= jmr) return
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    if (jml /= lmij) cycle
    ! t1s(5-5)   (11)drl(12)-
    iwdl = just(lri,lrj)
    iwdr = just(lrj,lri)
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    do mpl=1,mtype
      vplp_w0(mpl) = Zero
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1_t1s(5)
    end do
    call drl_act_c_link_ext_ab(lin,lrj)
    ! t1s(5-6)   drl(11)-c"(12)-
    iwdl = just(lri,lrj)
    iwdr = just(lrj,lri)
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    do mpl=1,mtype
      vplp_w0(mpl) = Zero
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1_t1s(6)
    end do
    call drl_act_c_link_ext_ab(lin,lri)
    ! t1s(5-7)   drl(12)-c"(11)-
    iwdl = just(lri,lrj)
    iwdr = just(lri,lrj)
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    do mpl=1,mtype
      vplp_w0(mpl) = Zero
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1_t1s(7)
    end do
    call drl_act_c_link_ext_ab(lin,lri)
  end do
end do

return

end subroutine tts_drl_act_c_ext_ab_sgt1

subroutine sdd_ar_act_bl_sgt0(lin,lra)
! sd1(8-1)    ar(01)-
! sd1(8-2)    (11)ar(23)-
! sd1(8-3)    ar(13)-c'(21)-
! sd1(8-4)    ar(23)-c'(11)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, &
                         mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_sd1
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, jmlr, lmi, lmk, lri, lrk, mpl, ni
real(kind=wp) :: w0sd1, w0sd2, w0sd3, w0sd4
integer(kind=iwp), external :: iwalk_ad

jmlr = Mul(jml,jmr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (lmi /= jmlr) cycle
  w0sd1 = w0_sd1(1)
  w0sd2 = w0_sd1(2)
  w0sd3 = w0_sd1(3)
  w0sd4 = w0_sd1(4)
  ni = mod(norb_dz-lri,2)
  if (ni == 1) then
    w0sd1 = -w0sd1
    w0sd2 = -w0sd2
    w0sd3 = -w0sd3
    w0sd4 = -w0sd4
  end if
  if ((jml == 1) .and. (lmi == jmr)) then
    ! sd1(8-1)    ar(01)-
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd1
      vplp_w1(mpl) = vplpnew_w1(mpl)*w0sd1
    end do
    iwdl = just(lri,lri)
    iwdr = jud(lri)
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    if (lin == 1) call ar_bl_ext_ss(lri,lra,1)
    if (lin == 2) call ar_bl_ext_st(lri,lra,1)
    if (lin == 3) call ar_bl_ext_ts(lri,lra,1)
    if (lin == 11) call ar_bl_ext_tt(lri,lra,1)
  end if
  ! sd1(8-2)    (11)ar(23)-
  do mpl=1,mtype
    vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd2
    vplp_w1(mpl) = vplpnew_w1(mpl)*w0sd2
  end do
  do lrk=norb_frz+1,lri-1
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lri,lrk)
    iwdr = jud(lrk)
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    if (lin == 1) call ar_bl_ext_ss(lri,lra,1)
    if (lin == 2) call ar_bl_ext_st(lri,lra,1)
    if (lin == 3) call ar_bl_ext_ts(lri,lra,1)
    if (lin == 11) call ar_bl_ext_tt(lri,lra,1)

  end do
  ! sd1(8-3)    ar(13)-c'(21)-
  do mpl=1,mtype
    vplp_w0(mpl) = -vplpnew_w0(mpl)*w0sd3
    vplp_w1(mpl) = -vplpnew_w1(mpl)*w0sd3
  end do
  do lrk=lri+1,norb_dz
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lrk,lri)
    iwdr = jud(lrk)
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    if (lin == 1) call ar_bl_ext_ss(lri,lra,1)
    if (lin == 2) call ar_bl_ext_st(lri,lra,1)
    if (lin == 3) call ar_bl_ext_ts(lri,lra,1)
    if (lin == 11) call ar_bl_ext_tt(lri,lra,1)
  end do
  ! sd1(8-4)    ar(23)-c'(11)-
  do mpl=1,mtype
    vplp_w0(mpl) = -vplpnew_w0(mpl)*w0sd4
    vplp_w1(mpl) = -vplpnew_w1(mpl)*w0sd4
  end do
  do lrk=lri+1,norb_dz
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lri,lrk)
    iwdr = jud(lrk)
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    if (lin == 1) call ar_bl_ext_ss(lri,lra,1)
    if (lin == 2) call ar_bl_ext_st(lri,lra,1)
    if (lin == 3) call ar_bl_ext_ts(lri,lra,1)
    if (lin == 11) call ar_bl_ext_tt(lri,lra,1)
  end do
end do

return

end subroutine sdd_ar_act_bl_sgt0

subroutine tttt_arbl_act_c_ext_ab_sgt0(lin)
! t1t1(12-1)  ar(13)-bl(31)-
! t1t1(12-1)  ar(13)-c'(11)-bl(31)-
! t1t1(12-1)  ar(13)-bl(31)-c"(11)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_t1t1, w1_t1t1
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lmj, lmk, lmki, lmkj, lri, lrj, lrk, mpl, ni
real(kind=wp) :: w0tt1, w1tt1
integer(kind=iwp), external :: iwalk_ad

do lri=norb_frz+1,norb_dz-1
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    w0tt1 = w0_t1t1(1)
    w1tt1 = w1_t1t1(1)
    ni = mod(lrj-lri,2)
    if (ni == 0) then
      w0tt1 = -w0tt1
      w1tt1 = -w1tt1
    end if
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0tt1
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1tt1
    end do
    ! t1t1(12-1)  (11)ar(13)-bl(31)-
    do lrk=norb_frz+1,lri-1
      lmk = lsm_inn(lrk)
      lmki = Mul(lmk,lmi)
      lmkj = Mul(lmk,lmj)
      if ((lmki == jml) .and. (lmkj == jmr)) then
        iwdl = just(lrk,lri)
        iwdr = just(lrk,lrj)
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        call arbl_act_c_link_ext_ab(lin,lri,lrj)
      end if
    end do
    ! t1t1(12-1)  ar(13)-bl(31)-c"(11)-
    do lrk=lrj+1,norb_dz
      lmk = lsm_inn(lrk)
      lmki = Mul(lmk,lmi)
      lmkj = Mul(lmk,lmj)
      if ((lmki == jml) .and. (lmkj == jmr)) then
        iwdl = just(lri,lrk)
        iwdr = just(lrj,lrk)
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        call arbl_act_c_link_ext_ab(lin,lri,lrj)
      end if
    end do
    do mpl=1,mtype
      vplp_w0(mpl) = -vplp_w0(mpl)
      vplp_w1(mpl) = -vplp_w1(mpl)
    end do
    ! t1t1(12-1)  ar(13)-c'(11)-bl(31)-
    do lrk=lri+1,lrj-1
      lmk = lsm_inn(lrk)
      lmki = Mul(lmk,lmi)
      lmkj = Mul(lmk,lmj)
      if ((lmki == jml) .and. (lmkj == jmr)) then
        iwdl = just(lri,lrk)
        iwdr = just(lrk,lrj)
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        call arbl_act_c_link_ext_ab(lin,lri,lrj)
      end if
    end do
  end do
end do

return

end subroutine tttt_arbl_act_c_ext_ab_sgt0

subroutine tttt_drl_act_c_ext_ab_sgt0(lin)
! t1t1(12-2)  drl(11)-
! t1t1(12-2)  drl(11)-c"(11)-
! t1t1(12-3)  drl(33)-
! t1t1(12-3)  drl(33)-c"(11)-
! t1t1(12-3)  drl(33)-c"(11)-c"(11)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_t1t1, w1_t1t1
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lmij, lmj, lri, lrj, mpl
real(kind=wp) :: w0tt2, w0tt3, w1tt2
integer(kind=iwp), external :: iwalk_ad

w0tt2 = w0_t1t1(2)
w1tt2 = w1_t1t1(2)
w0tt3 = w0_t1t1(3)
do lri=norb_frz+1,norb_dz-1
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    if ((lmij /= jml) .or. (lmij /= jmr)) cycle
    ! t1t1(12-2)  drl(11)-
    ! t1t1(12-2)  drl(11)-c"(11)-
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0tt2
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1tt2
    end do
    iwdl = just(lri,lrj)
    iwdr = iwdl
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    call drl_act_c_link_ext_ab(lin,lri)
    call drl_act_c_link_ext_ab(lin,lrj)
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0tt3
      vplp_w1(mpl) = Zero
    end do
    ! t1t1(12-3)  drl(33)-
    ! t1t1(12-3)  drl(33)-c"(11)-
    ! t1t1(12-3)  drl(33)-c"(11)-c"(11)-
    !do lrk=1,norb_dz
    !  if (lrk == lri) cycle
    !  if (lrk == lrj) cycle
    !  call drl_ss_ext(lrk)
    !end do
    call drl_act_c_link_ext_ab_sum(lin,lri,lrj)
  end do
end do

return

end subroutine tttt_drl_act_c_ext_ab_sgt0

subroutine ttdd_ar_act_bl_sgt1(lin,lra)
! t1d1(15-1)  ar(13)-
! t1d1(15-1)  ar(13)-c'(11)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, &
                         mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_t1d1
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, jmlr, lmi, lmk, lri, lrk, mpl, ni
real(kind=wp) :: w0td1
integer(kind=iwp), external :: iwalk_ad

jmlr = Mul(jml,jmr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (lmi /= jmlr) cycle
  w0td1 = w0_t1d1(1)
  ni = mod(norb_dz-lri,2)
  if (ni == 1) w0td1 = -w0td1
  !---------------------------------------------------------------------
  ! t1d1(15-1)  (11)ar(13)-
  do lrk=norb_frz+1,lri-1
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lrk,lri)
    iwdr = jud(lrk)
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0td1
      vplp_w1(mpl) = vplpnew_w1(mpl)*w0td1
    end do
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    if (lin == 1) call ar_bl_ext_ss(lri,lra,1)
    if (lin == 2) call ar_bl_ext_st(lri,lra,1)
    if (lin == 3) call ar_bl_ext_ts(lri,lra,1)
    if (lin == 11) call ar_bl_ext_tt(lri,lra,1)
  end do
  !---------------------------------------------------------------------
  ! t1d1(15-1)  ar(13)-c'(11)-
  do lrk=lri+1,norb_dz
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lri,lrk)
    iwdr = jud(lrk)
    do mpl=1,mtype
      vplp_w0(mpl) = -vplpnew_w0(mpl)*w0td1
      vplp_w1(mpl) = -vplpnew_w1(mpl)*w0td1
    end do
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    if (lin == 1) call ar_bl_ext_ss(lri,lra,1)
    if (lin == 2) call ar_bl_ext_st(lri,lra,1)
    if (lin == 3) call ar_bl_ext_ts(lri,lra,1)
    if (lin == 11) call ar_bl_ext_tt(lri,lra,1)
  end do
end do

return

end subroutine ttdd_ar_act_bl_sgt1

subroutine d1d1_arbl_act_c_ext_ab_sgt0(lin)
! d1d1(20-1) ar(13)-bl(31)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, jud, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_d1d1, w1_d1d1
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lmj, lri, lrj, mpl, ni
real(kind=wp) :: w0dd1, w1dd1
integer(kind=iwp), external :: iwalk_ad

do lri=norb_frz+1,norb_dz-1
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    if ((jml /= lmi) .and. (jmr /= lmj)) cycle
    w0dd1 = w0_d1d1(1)
    w1dd1 = w1_d1d1(1)
    ni = mod(lrj-lri,2)
    if (ni == 0) then
      w0dd1 = -w0dd1
      w1dd1 = -w1dd1
    end if
    if ((lmi == jml) .and. (lmj == jmr)) then
      ! d1d1(20-1) ar(13)-bl(31)-
      do mpl=1,mtype
        vplp_w0(mpl) = vplpnew_w0(mpl)*w0dd1
        vplp_w1(mpl) = vplpnew_w1(mpl)*w1dd1
      end do
      iwdl = jud(lri)
      iwdr = jud(lrj)
      do mpl=1,mhlp
        iwal = lpnew_lwei(mpl)
        iwar = lpnew_rwei(mpl)
        lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
        lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
      end do
      call arbl_act_c_link_ext_ab(lin,lri,lrj)
    end if
  end do
end do

return

end subroutine d1d1_arbl_act_c_ext_ab_sgt0

subroutine d1d1_drl_act_c_ext_ab_sgt0(lin)
! d1d1(20-2) drl(11)-
! d1d1(20-3) drl(33)-
! d1d1(20-3) drl(33)-c"(11)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, jud, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_d1d1, w1_d1d1
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lri, lrk, mpl
real(kind=wp) :: w0dd2, w0dd3, w1dd2
integer(kind=iwp), external :: iwalk_ad

if (jml /= jmr) return
w0dd2 = w0_d1d1(2)
w1dd2 = w1_d1d1(2)
w0dd3 = w0_d1d1(3)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (lmi /= jml) cycle
  ! d1d1(20-2) drl(11)-
  do mpl=1,mtype
    vplp_w0(mpl) = vplpnew_w0(mpl)*w0dd2
    vplp_w1(mpl) = vplpnew_w1(mpl)*w1dd2
  end do
  iwdl = jud(lri)
  iwdr = iwdl
  do mpl=1,mhlp
    iwal = lpnew_lwei(mpl)
    iwar = lpnew_rwei(mpl)
    lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
    lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
  end do
  call drl_act_c_link_ext_ab(lin,lri)
  do mpl=1,mtype
    vplp_w0(mpl) = vplpnew_w0(mpl)*w0dd3
    vplp_w1(mpl) = Zero
  end do
  ! d1d1(20-3) drl(33)-
  ! d1d1(20-3) drl(33)-c"(11)-
  do lrk=1,norb_dz
    if (lrk == lri) cycle
    call drl_act_c_link_ext_ab(lin,lrk)
  end do
  !call drl_act_c_link_ext_ab_sum(lin,lri,0)
end do

return

end subroutine d1d1_drl_act_c_ext_ab_sgt0

subroutine dd1_arbl_act_c_ext_ab_sgt0(lin)
! dd1(21) ar(23)-bl(31)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, jud, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_dd1, w1_dd1
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, jmlr, lmi, lmij, lmj, lri, lrj, mpl, ni
real(kind=wp) :: w0dd1, w1dd1
integer(kind=iwp), external :: iwalk_ad

jmlr = Mul(jml,jmr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    if (lmij /= jmlr) cycle
    w0dd1 = w0_dd1
    w1dd1 = w1_dd1
    ni = mod(lrj-lri,2)
    if (ni == 0) then
      w0dd1 = -w0dd1
      w1dd1 = -w1dd1
    end if
    if ((lmi == jml) .and. (lmj == jmr)) then
      ! dd1(21)ar(23)-bl(31)-
      do mpl=1,mtype
        vplp_w0(mpl) = vplpnew_w0(mpl)*w0dd1
        vplp_w1(mpl) = vplpnew_w1(mpl)*w1dd1
      end do
      iwdl = jud(lri)
      iwdr = jud(lrj)
      do mpl=1,mhlp
        iwal = lpnew_lwei(mpl)
        iwar = lpnew_rwei(mpl)
        lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
        lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
      end do
      call arbl_act_c_link_ext_ab(lin,lri,lrj)
    end if
  end do
end do

return

end subroutine dd1_arbl_act_c_ext_ab_sgt0

subroutine d1d_arbl_act_c_ext_ab_sgt0(lin)
! d1d(22-1)   ar(13)-bl(32)-
! d1d(22-2)   drl(12)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, jud, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_d1d, w1_d1d
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, jmlr, lmi, lmij, lmj, lri, lrj, mpl, ni
real(kind=wp) :: w0dd1, w1dd1
integer(kind=iwp), external :: iwalk_ad

jmlr = Mul(jml,jmr)
do lri=norb_frz+1,norb_dz-1
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    if (lmij /= jmlr) cycle
    w0dd1 = w0_d1d(1)
    w1dd1 = w1_d1d(1)
    ni = mod(lrj-lri,2)
    if (ni == 0) then
      w0dd1 = -w0dd1
      w1dd1 = -w1dd1
    end if
    if ((lmi == jml) .and. (lmj == jmr)) then
      ! d1d(22-1)   ar(13)-bl(32)-
      do mpl=1,mtype
        vplp_w0(mpl) = vplpnew_w0(mpl)*w0dd1
        vplp_w1(mpl) = vplpnew_w1(mpl)*w1dd1
      end do
      iwdl = jud(lri)
      iwdr = jud(lrj)
      do mpl=1,mhlp
        iwal = lpnew_lwei(mpl)
        iwar = lpnew_rwei(mpl)
        lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
        lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
      end do
      call arbl_act_c_link_ext_ab(lin,lri,lrj)
    end if
  end do
end do

return

end subroutine d1d_arbl_act_c_ext_ab_sgt0

subroutine d1d_drl_act_c_ext_ab_sgt0(lin)
! d1d(22-2)   drl(12)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, jud, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_d1d, w1_d1d
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lri, mpl
integer(kind=iwp), external :: iwalk_ad

if (jmr /= jml) return
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (jml /= lmi) cycle
  do mpl=1,mtype
    vplp_w0(mpl) = vplpnew_w0(mpl)*w0_d1d(2)
    vplp_w1(mpl) = vplpnew_w1(mpl)*w1_d1d(2)
  end do
  iwdl = jud(lri)
  iwdr = jud(lri)
  do mpl=1,mhlp
    iwal = lpnew_lwei(mpl)
    iwar = lpnew_rwei(mpl)
    lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
    lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
  end do
  call drl_act_c_link_ext_ab(lin,lri)
end do

return

end subroutine d1d_drl_act_c_ext_ab_sgt0

subroutine d1v_ar_act_bl_ext_ab_sgt0(lin,lra)
! d1v(24-1)  ar(13)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, jud, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_d1v
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, jmlr, lmi, lri, mpl, ni
real(kind=wp) :: w0dv1
integer(kind=iwp), external :: iwalk_ad

jmlr = Mul(jml,jmr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (lmi /= jmlr) cycle
  w0dv1 = w0_d1v(1)
  ni = mod(norb_dz-lri,2)
  if (ni == 1) w0dv1 = -w0dv1
  ! d1v(24-1)  ar(13)-
  iwdl = jud(lri)
  iwdr = 0
  do mpl=1,mhlp
    iwal = lpnew_lwei(mpl)
    iwar = lpnew_rwei(mpl)
    lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
    lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
  end do
  do mpl=1,mtype
    vplp_w0(mpl) = vplpnew_w0(mpl)*w0dv1
    vplp_w1(mpl) = vplpnew_w1(mpl)*w0dv1
  end do
  call arbl_act_c_link_ext_ab(lin,lri,lra)
end do

return

end subroutine d1v_ar_act_bl_ext_ab_sgt0

subroutine ar_br_stv_link_ext_brar(lin,lri,lrj)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lri, lrj

if (lin == 17) call ar_br_tv_ext_br_ar(lri,lrj)
if (lin == 10) call ar_br_sv_ext_br_ar(lri,lrj)

end subroutine ar_br_stv_link_ext_brar

subroutine sd_ar_act_br_sgt0(lin,lra)
! sd(6-3) a&r(13)c'(22)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, &
                         mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_sd
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, jmlr, lmi, lmk, lri, lrk, mpl, ni
real(kind=wp) :: w0sd3
integer(kind=iwp), external :: iwalk_ad

jmlr = Mul(jml,jmr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (lmi /= jmlr) cycle
  w0sd3 = w0_sd(3)
  ni = mod(norb_dz-lri,2)
  if (ni == 1) w0sd3 = -w0sd3

  ! sd(6-3) a&r(13)c'(22)-
  do mpl=1,mtype
    vplp_w0(mpl) = -vplpnew_w0(mpl)*w0sd3
    vplp_w1(mpl) = -vplpnew_w1(mpl)*w0sd3
  end do
  do lrk=lri+1,norb_dz
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lrk,lri)
    iwdr = jud(lrk)
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    call ar_br_stv_link_ext_brar(lin,lri,lra)
  end do
end do

return

end subroutine sd_ar_act_br_sgt0

subroutine sdd_ar_act_br_sgt0(lin,lra)
! sd1(8-1)    ar(01)-
! sd1(8-2)    (11)ar(23)-
! sd1(8-3)    ar(13)-c'(21)-
! sd1(8-4)    ar(23)-c'(11)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, &
                         mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_sd1
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, jmlr, lmi, lmk, lri, lrk, mpl, ni
real(kind=wp) :: w0sd1, w0sd2, w0sd3, w0sd4
integer(kind=iwp), external :: iwalk_ad

jmlr = Mul(jml,jmr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (lmi /= jmlr) cycle
  w0sd1 = w0_sd1(1)
  w0sd2 = w0_sd1(2)
  w0sd3 = w0_sd1(3)
  w0sd4 = w0_sd1(4)
  ni = mod(norb_dz-lri,2)
  if (ni == 1) then
    w0sd1 = -w0sd1
    w0sd2 = -w0sd2
    w0sd3 = -w0sd3
    w0sd4 = -w0sd4
  end if
  if ((jml == 1) .and. (lmi == jmr)) then
    ! sd1(8-1)    ar(01)-
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd1
      vplp_w1(mpl) = vplpnew_w1(mpl)*w0sd1
    end do
    iwdl = just(lri,lri)
    iwdr = jud(lri)
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    call ar_br_stv_link_ext_brar(lin,lri,lra)
  end if
  ! sd1(8-2)    (11)ar(23)-
  do mpl=1,mtype
    vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd2
    vplp_w1(mpl) = vplpnew_w1(mpl)*w0sd2
  end do
  do lrk=norb_frz+1,lri-1
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lri,lrk)
    iwdr = jud(lrk)
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    call ar_br_stv_link_ext_brar(lin,lri,lra)
  end do
  ! sd1(8-3)    ar(13)-c'(21)-
  do mpl=1,mtype
    vplp_w0(mpl) = -vplpnew_w0(mpl)*w0sd3
    vplp_w1(mpl) = -vplpnew_w1(mpl)*w0sd3
  end do
  do lrk=lri+1,norb_dz
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lrk,lri)
    iwdr = jud(lrk)
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    call ar_br_stv_link_ext_brar(lin,lri,lra)
  end do
  ! sd1(8-4)    ar(23)-c'(11)-
  do mpl=1,mtype
    vplp_w0(mpl) = -vplpnew_w0(mpl)*w0sd4
    vplp_w1(mpl) = -vplpnew_w1(mpl)*w0sd4
  end do
  do lrk=lri+1,norb_dz
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lri,lrk)
    iwdr = jud(lrk)
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    call ar_br_stv_link_ext_brar(lin,lri,lra)
  end do
end do

return

end subroutine sdd_ar_act_br_sgt0

subroutine ttdd_ar_act_br_sgt1(lin,lra)
! t1d1(15-1)  ar(13)-
! t1d1(15-1)  ar(13)-c'(11)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, &
                         mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_t1d1
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, jmlr, lmi, lmk, lri, lrk, mpl, ni
real(kind=wp) :: w0td1
integer(kind=iwp), external :: iwalk_ad

jmlr = Mul(jml,jmr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (lmi /= jmlr) cycle
  w0td1 = w0_t1d1(1)
  ni = mod(norb_dz-lri,2)
  if (ni == 1) w0td1 = -w0td1
  !---------------------------------------------------------------------
  ! t1d1(15-1)  (11)ar(13)-
  do lrk=norb_frz+1,lri-1
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lrk,lri)
    iwdr = jud(lrk)
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0td1
      vplp_w1(mpl) = vplpnew_w1(mpl)*w0td1
    end do
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    call ar_br_stv_link_ext_brar(lin,lri,lra)
  end do
  !---------------------------------------------------------------------
  ! t1d1(15-1)  ar(13)-c'(11)-
  do lrk=lri+1,norb_dz
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lri,lrk)
    iwdr = jud(lrk)
    do mpl=1,mtype
      vplp_w0(mpl) = -vplpnew_w0(mpl)*w0td1
      vplp_w1(mpl) = -vplpnew_w1(mpl)*w0td1
    end do
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    call ar_br_stv_link_ext_brar(lin,lri,lra)
  end do
end do

return

end subroutine ttdd_ar_act_br_sgt1

subroutine ttv_arbr_act_c_stv_sgt1(lin)
! t1v(18) ar(13)-br(13)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w1, w1_t1v
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, jmlr, lmi, lmij, lmj, lri, lrj, mpl, ni
real(kind=wp) :: w1
integer(kind=iwp), external :: iwalk_ad

jmlr = Mul(jml,jmr)
do lri=norb_frz+1,norb_dz-1
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    if (lmij /= jmlr) cycle
    !-------------------------------------------------------------------
    iwdl = just(lri,lrj)
    iwdr = 0
    ni = mod(lrj-lri,2)
    w1 = w1_t1v
    if (ni == 0) w1 = -w1
    do mpl=1,mtype
      vplp_w0(mpl) = Zero
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1
    end do
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    call ar_br_stv_link_ext_brar(lin,lri,lrj)
  end do
end do

return

end subroutine ttv_arbr_act_c_stv_sgt1

subroutine ss_s_drl_act_c_ext_ab_sgt0(lin)
! ss(1-19) drl(12)-c"(21)-       act -c"-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_ss, w1_ss
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lmij, lmj, lri, lrj, mpl
integer(kind=iwp), external :: iwalk_ad

if (jml /= jmr) return
do lri=norb_frz+1,norb_dz-1
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    if (lmij /= jml) cycle
    ! ss(1-19) drl(12)-c"(21)-       act -c"-
    iwdl = just(lrj,lri)
    iwdr = just(lri,lrj)
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0_ss(19)
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1_ss(19)
    end do
    call drl_act_c_link_ext_ab(lin,lri)
  end do
end do

return

end subroutine ss_s_drl_act_c_ext_ab_sgt0

subroutine ds_arblbr_act_c1_sgt0(lin)
!=======================================================================
! ds(7-2) ar(23)-bl(31)-br(32)-
!=======================================================================

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, lsm_inn, mhlp, mtype, ngw2, ngw3, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w1, w1_ds
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, lmd, lmi, lmij, lmj, lrd, lri, lrj, mpl, ni
real(kind=wp) :: w1ds
integer(kind=iwp), external :: iwalk_ad

isma = Mul(iml,imr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    if (lmij /= jmr) cycle
    ! ds(7-2) ar(23)-bl(31)-br(32)-
    do lrd=norb_frz+1,lri-1
      lmd = lsm_inn(lrd)
      if (lmd /= jml) cycle
      ijk = lrd-norb_frz+ngw2(lri-norb_frz)+ngw3(lrj-norb_frz)
      intpos = intind_ijka(ijk)
      w1ds = w1_ds(2)
      ni = mod(norb_dz-lrj+lri-lrd,2)
      if (ni == 0) w1ds = -w1ds

      iwdr = just(lrj,lri)
      iwdl = jud(lrd)
      do mpl=1,mhlp
        iwal = lpnew_lwei(mpl)
        iwar = lpnew_rwei(mpl)
        lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
        lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
      end do
      do mpl=1,mtype
        vplp_w0(mpl) = Zero
        vplp_w1(mpl) = vplpnew_w1(mpl)*w1ds
      end do
      call ar_bl_br_ext_al_new(lin,intpos,isma,1)
    end do
  end do
end do

return

end subroutine ds_arblbr_act_c1_sgt0

subroutine sdd_ar_act_dlr_sgt0(lin,lra)
! sd1(8-1)    ar(01)-
! sd1(8-2)    (11)ar(23)-
! sd1(8-3)    ar(13)-c'(21)-
! sd1(8-4)    ar(23)-c'(11)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, &
                         mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_sd1, w1_sd1
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lmj, lri, lrj, mpl
real(kind=wp) :: w0sd1, w0sd2, w0sd3, w0sd4, w1sd1, w1sd2, w1sd3, w1sd4
integer(kind=iwp), external :: iwalk_ad

do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  w0sd1 = w0_sd1(1)
  w1sd1 = w1_sd1(1)
  w0sd2 = w0_sd1(2)
  w1sd2 = w1_sd1(2)
  w0sd3 = w0_sd1(3)
  w1sd3 = w1_sd1(3)
  w0sd4 = w0_sd1(4)
  w1sd4 = w1_sd1(4)
  if (mod(norb_dz-lri,2) == 1) then
    w0sd1 = -w0sd1
    w1sd1 = -w1sd1
    w0sd2 = -w0sd2
    w1sd2 = -w1sd2
    w0sd3 = -w0sd3
    w1sd3 = -w1sd3
    w0sd4 = -w0sd4
    w1sd4 = -w1sd4
  end if
  ! sd1(8-1)    ar(01)-
  if ((jml == 1) .and. (jmr == lmi)) then
    iwdl = just(lri,lri)
    iwdr = jud(lri)
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd1
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1sd1
    end do
    call ar_drl_ext_al_new(lin,lri,lra)
  end if
  ! sd1(8-2)    (11)ar(23)-
  do mpl=1,mtype
    vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd2
    vplp_w1(mpl) = vplpnew_w1(mpl)*w1sd2
  end do
  do lrj=norb_frz+1,lri-1
    lmj = lsm_inn(lrj)
    if ((jml /= Mul(lmi,lmj)) .or. (jmr /= lmj)) cycle
    iwdl = just(lri,lrj)
    iwdr = jud(lrj)
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    call ar_drl_ext_al_new(lin,lri,lra)
  end do
  ! sd1(8-3)    ar(13)-c'(21)-
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    if ((jml /= Mul(lmi,lmj)) .or. (jmr /= lmj)) cycle
    iwdl = just(lrj,lri)
    iwdr = jud(lrj)
    do mpl=1,mtype
      vplp_w0(mpl) = -vplpnew_w0(mpl)*w0sd3
      vplp_w1(mpl) = -vplpnew_w1(mpl)*w1sd3
    end do
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    call ar_drl_ext_al_new(lin,lri,lra)
    ! sd1(8-4)    ar(23)-c'(11)-
    iwdl = just(lri,lrj)
    iwdr = jud(lrj)
    do mpl=1,mtype
      vplp_w0(mpl) = -vplpnew_w0(mpl)*w0sd4
      vplp_w1(mpl) = -vplpnew_w1(mpl)*w1sd4
    end do
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    call ar_drl_ext_al_new(lin,lri,lra)
  end do
end do

return

end subroutine sdd_ar_act_dlr_sgt0

!subroutine tts_arbl_act_c_sgt1(lin,lra)
!!=======================================================================
!! tts(5) a&r-b^l-  act -b&l ............................................
!! t1s(5-1)   ar(13)-bl(10)-
!! t1s(5-2)   ar(13)-bl(32)-
!! t1s(5-2)   ar(13)-c'(11)-bl(32)-
!! t1s(5-3)   ar(13)-bl(31)-c"(12)-
!! t1s(5-4)   ar(13)-bl(32)-c"(11)-
!
!use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, &
!                         lpnew_rwei, lsm_inn, mhlp, mtype, ngw2, ngw3, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w1, w1_t1s
!use Symmetry_Info, only: Mul
!use Constants, only: Zero
!use Definitions, only: wp, iwp
!
!implicit none
!integer(kind=iwp), intent(in) :: lin, lra
!integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, lmi, lmij, lmj, lmk, lri, lrj, lrk, mpl
!real(kind=wp) :: w1ts1, w1ts2, w1ts3, w1ts4
!integer(kind=iwp), external :: iwalk_ad
!
!isma = Mul(iml,imr)
!do lri=norb_frz+1,norb_dz
!  lmi = lsm_inn(lri)
!  do lrj=lri+1,norb_dz
!    lmj = lsm_inn(lrj)
!    lmij = Mul(lmi,lmj)
!    w1ts1 = w1_t1s(1)
!    w1ts2 = w1_t1s(2)
!    w1ts3 = w1_t1s(3)
!    w1ts4 = w1_t1s(4)
!    if (mod(lrj-lri,2) == 0) then
!      w1ts1 = -w1ts1
!      w1ts2 = -w1ts2
!      w1ts3 = -w1ts3
!      w1ts4 = -w1ts4
!    end if
!    ijk = lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz) !???
!    intpos = intind_ijka(ijk)                         !???
!    ! t1s(5-1)   ar(13)-bl(10)-
!    if ((jmr == 1) .and. (jml == lmij)) then
!      iwdl = just(lri,lrj)
!      iwdr = just(lrj,lrj)
!      do mpl=1,mhlp
!        iwal = lpnew_lwei(mpl)
!        iwar = lpnew_rwei(mpl)
!        lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
!        lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
!      end do
!      do mpl=1,mtype
!        vplp_w0(mpl) = Zero
!        vplp_w1(mpl) = vplpnew_w1(mpl)*w1ts1
!      end do
!      call ar_bl_br_ext_al_new(lin,intpos,isma,1)
!    end if
!    ! t1s(5-2)   (11)ar(13)-bl(32)-
!    do lrk=norb_frz+1,lri-1
!      lmk = lsm_inn(lrk)
!      if ((jmr == Mul(lmk,lmi)) .and. (jml == Mul(lmk,lmj))) then
!        iwdl = just(lrk,lri)
!        iwdr = just(lrj,lrk)
!        do mpl=1,mhlp
!          iwal = lpnew_lwei(mpl)
!          iwar = lpnew_rwei(mpl)
!          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
!          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
!        end do
!        do mpl=1,mtype
!          vplp_w0(mpl) = Zero
!          vplp_w1(mpl) = vplpnew_w1(mpl)*w1ts2
!        end do
!        call ar_bl_br_ext_al_new(lin,intpos,isma,1)
!      end if
!    end do
!    ! t1s(5-2)   ar(13)-c'(11)-bl(32)-
!    do lrk=lri+1,lrj-1
!      lmk = lsm_inn(lrk)
!      if ((jmr == Mul(lmi,lmk)) .and. (jml == Mul(lmk,lmj))) then
!        iwdl = just(lri,lrk)
!        iwdr = just(lrj,lrk)
!        do mpl=1,mhlp
!          iwal = lpnew_lwei(mpl)
!          iwar = lpnew_rwei(mpl)
!          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
!          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
!        end do
!        do mpl=1,mtype
!          vplp_w0(mpl) = Zero
!          vplp_w1(mpl) = vplpnew_w1(mpl)*w1ts2
!        end do
!        call ar_bl_br_ext_al_new(lin,intpos,isma,1)
!      end if
!    end do
!    ! t1s(5-3)   ar(13)-bl(31)-c"(12)-
!    do lrk=lrj+1,norb_dz
!      lmk = lsm_inn(lrk)
!      if ((jmr == Mul(lmi,lmk)) .and. (jml == Mul(lmj,lmk))) then
!        iwdl = just(lri,lrk)
!        iwdr = just(lrk,lrj)
!        do mpl=1,mhlp
!          iwal = lpnew_lwei(mpl)
!          iwar = lpnew_rwei(mpl)
!          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
!          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
!        end do
!        do mpl=1,mtype
!          vplp_w0(mpl) = Zero
!          vplp_w1(mpl) = vplpnew_w1(mpl)*w1ts3
!        end do
!        call ar_bl_br_ext_al_new(lin,intpos,isma,1)
!        ! t1s(5-4)   ar(13)-bl(32)-c"(11)-
!        iwdl = just(lri,lrk)
!        iwdr = just(lrj,lrk)
!        do mpl=1,mhlp
!          iwal = lpnew_lwei(mpl)
!          iwar = lpnew_rwei(mpl)
!          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
!          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
!        end do
!        do mpl=1,mtype
!          vplp_w0(mpl) = Zero
!          vplp_w1(mpl) = vplpnew_w1(mpl)*w1ts4
!        end do
!        call ar_bl_br_ext_al_new(lin,intpos,isma,1)
!      end if
!    end do
!  end do
!end do
!
!return
!
!end subroutine tts_arbl_act_c_sgt1

subroutine tts_arbl_act_br_sgt1(lin,lra)
!=======================================================================
! tts(5) a&r-b^l-  act -b&l ............................................
! t1s(5-1)   ar(13)-bl(10)-
! t1s(5-2)   ar(13)-bl(32)-
! t1s(5-2)   ar(13)-c'(11)-bl(32)-
! t1s(5-3)   ar(13)-bl(31)-c"(12)-
! t1s(5-4)   ar(13)-bl(32)-c"(11)-

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

isma = Mul(iml,imr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    w1ts1 = w1_t1s(1)
    w1ts2 = w1_t1s(2)
    w1ts3 = w1_t1s(3)
    w1ts4 = w1_t1s(4)
    if (mod(lrj-lri,2) == 0) then
      w1ts1 = -w1ts1
      w1ts2 = -w1ts2
      w1ts3 = -w1ts3
      w1ts4 = -w1ts4
    end if
    ijk = lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz) !???
    intpos = intind_ijka(ijk)                         !???
    ! t1s(5-1)   ar(13)-bl(10)-
    if ((jmr == 1) .and. (jml == lmij)) then
      iwdl = just(lri,lrj)
      iwdr = just(lrj,lrj)
      do mpl=1,mhlp
        iwal = lpnew_lwei(mpl)
        iwar = lpnew_rwei(mpl)
        lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
        lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
      end do
      do mpl=1,mtype
        vplp_w0(mpl) = Zero
        vplp_w1(mpl) = vplpnew_w1(mpl)*w1ts1
      end do
      call ar_bl_br_ext_al_new(lin,intpos,isma,1)
    end if
    ! t1s(5-2)   (11)ar(13)-bl(32)-
    do lrk=norb_frz+1,lri-1
      lmk = lsm_inn(lrk)
      if ((jml == Mul(lmk,lmi)) .and. (jmr == Mul(lmk,lmj))) then
        iwdl = just(lrk,lri)
        iwdr = just(lrj,lrk)
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        do mpl=1,mtype
          vplp_w0(mpl) = Zero
          vplp_w1(mpl) = vplpnew_w1(mpl)*w1ts2
        end do
        call ar_bl_br_ext_al_new(lin,intpos,isma,1)
      end if
    end do
    ! t1s(5-2)   ar(13)-c'(11)-bl(32)-
    do lrk=lri+1,lrj-1
      lmk = lsm_inn(lrk)
      if ((jml == Mul(lmi,lmk)) .and. (jmr == Mul(lmk,lmj))) then
        iwdl = just(lri,lrk)
        iwdr = just(lrj,lrk)
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        do mpl=1,mtype
          vplp_w0(mpl) = Zero
          vplp_w1(mpl) = -vplpnew_w1(mpl)*w1ts2
        end do
        call ar_bl_br_ext_al_new(lin,intpos,isma,1)
      end if
    end do
    ! t1s(5-3)   ar(13)-bl(31)-c"(12)-
    do lrk=lrj+1,norb_dz
      lmk = lsm_inn(lrk)
      if ((jml == Mul(lmi,lmk)) .and. (jmr == Mul(lmj,lmk))) then
        iwdl = just(lri,lrk)
        iwdr = just(lrk,lrj)
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        do mpl=1,mtype
          vplp_w0(mpl) = Zero
          vplp_w1(mpl) = vplpnew_w1(mpl)*w1ts3
        end do
        call ar_bl_br_ext_al_new(lin,intpos,isma,1)
        ! t1s(5-4)   ar(13)-bl(32)-c"(11)-
        iwdl = just(lri,lrk)
        iwdr = just(lrj,lrk)
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        do mpl=1,mtype
          vplp_w0(mpl) = Zero
          vplp_w1(mpl) = vplpnew_w1(mpl)*w1ts4
        end do
        call ar_bl_br_ext_al_new(lin,intpos,isma,1)
      end if
    end do
  end do
end do

return

end subroutine tts_arbl_act_br_sgt1

subroutine sd_ar_act_dlr_sgt0(lin,lra)
! sd(6-3) a&r(13)c'(22)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, &
                         mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_sd
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, jmlr, lmi, lmk, lri, lrk, mpl, ni
real(kind=wp) :: w0sd3
integer(kind=iwp), external :: iwalk_ad

jmlr = Mul(jml,jmr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (lmi /= jmlr) cycle
  w0sd3 = w0_sd(3)
  ni = mod(norb_dz-lri,2)
  if (ni == 0) then
    w0sd3 = -w0sd3
  end if
  !---------------------------------------------------------------------
  ! sd(6-3) a&r(13)c'(22)-
  do mpl=1,mtype
    vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd3
    vplp_w1(mpl) = vplpnew_w1(mpl)*w0sd3
  end do
  do lrk=lri+1,norb_dz
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lrk,lri)
    iwdr = jud(lrk)
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    call ar_drl_ext_al_new(lin,lri,lra)
  end do
end do

return

end subroutine sd_ar_act_dlr_sgt0

subroutine sdd_ar_act_blbr_sgt0(lin,jk)
! sdd(8-1) a&r(01)-
! sdd(8-2) (11)a&(23)-
! sdd(8-3) a&r(13)c'(21)-
! sdd(8-4) a&r(23)c'(11)-

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, lsm_inn, mhlp, mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_sd1, &
                         w1_sd1
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, jk
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, lmi, lmk, lri, lrk, mpl
real(kind=wp) :: w0sd1, w0sd2, w0sd3, w0sd4, w1sd1, w1sd2, w1sd3, w1sd4
integer(kind=iwp), external :: iwalk_ad

isma = Mul(iml,imr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  w0sd1 = w0_sd1(1)
  w1sd1 = w1_sd1(1)
  w0sd2 = w0_sd1(2)
  w1sd2 = w1_sd1(2)
  w0sd3 = w0_sd1(3)
  w1sd3 = w1_sd1(3)
  w0sd4 = w0_sd1(4)
  w1sd4 = w1_sd1(4)
  if (mod(norb_dz-lri,2) == 1) then
    w0sd1 = -w0sd1
    w1sd1 = -w1sd1
    w0sd2 = -w0sd2
    w1sd2 = -w1sd2
    w0sd3 = -w0sd3
    w1sd3 = -w1sd3
    w0sd4 = -w0sd4
    w1sd4 = -w1sd4
  end if
  ijk = lri-norb_frz+jk
  intpos = intind_ijka(ijk)
  ! sdd(8-1) a&r(01)-
  if ((jml == 1) .and. (lmi == jmr)) then
    iwdl = just(lri,lri)
    iwdr = jud(lri)
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd1
      vplp_w1(mpl) = vplpnew_w1(mpl)*w0sd1
    end do
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    call ar_bl_br_ext_al_new(lin,intpos,isma,1)
  end if
  !---------------------------------------------------------------------
  ! sdd(8-2) c(11)a&(23)-
  do mpl=1,mtype
    vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd2
    vplp_w1(mpl) = vplpnew_w1(mpl)*w0sd2
  end do
  do lrk=norb_frz+1,lri-1
    lmk = lsm_inn(lrk)
    if ((lmk /= jmr) .or. (jml /= Mul(lmk,lmi))) cycle
    iwdl = just(lri,lrk)
    iwdr = jud(lrk)
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    call ar_bl_br_ext_al_new(lin,intpos,isma,1)
  end do
  !---------------------------------------------------------------------
  ! sdd(8-3) a&r(13)c'(21)-
  do lrk=lri+1,norb_dz
    lmk = lsm_inn(lrk)
    if ((lmk /= jmr) .or. (jml /= Mul(lmi,lmk))) cycle
    iwdl = just(lrk,lri)
    iwdr = jud(lrk)
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    do mpl=1,mtype
      vplp_w0(mpl) = -vplpnew_w0(mpl)*w0sd3
      vplp_w1(mpl) = -vplpnew_w1(mpl)*w0sd3
    end do
    call ar_bl_br_ext_al_new(lin,intpos,isma,1)
    ! sdd(8-4) a&r(23)c'(11)-
    iwdl = just(lri,lrk)
    iwdr = jud(lrk)
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    do mpl=1,mtype
      vplp_w0(mpl) = -vplpnew_w0(mpl)*w0sd4
      vplp_w1(mpl) = -vplpnew_w1(mpl)*w0sd4
    end do
    call ar_bl_br_ext_al_new(lin,intpos,isma,1)
  end do
end do

return

end subroutine sdd_ar_act_blbr_sgt0

subroutine sd_ar_act_blbr_sgt0(lin,jk)
! sd(6-3) a&r(13)c'(22)-

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, lsm_inn, mhlp, mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_sd
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, jk
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, jmlr, lmi, lmk, lri, lrk, mpl, ni
real(kind=wp) :: w0sd3
integer(kind=iwp), external :: iwalk_ad

jmlr = Mul(jml,jmr)
isma = Mul(iml,imr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (lmi /= jmlr) cycle
  ijk = lri-norb_frz+jk
  intpos = intind_ijka(ijk)
  w0sd3 = w0_sd(3)
  ni = mod(norb_dz-lri,2)
  if (ni == 0) then
    w0sd3 = -w0sd3
  end if
  !---------------------------------------------------------------------
  ! sd(6-3) a&r(13)c'(22)-
  do mpl=1,mtype
    vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd3
    vplp_w1(mpl) = vplpnew_w1(mpl)*w0sd3
  end do
  do lrk=lri+1,norb_dz
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lrk,lri)
    iwdr = jud(lrk)
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    call ar_bl_br_ext_al_new(lin,intpos,isma,1)
  end do
end do

return

end subroutine sd_ar_act_blbr_sgt0

subroutine sdd_ar_act_c_sd_ext_sgt0(lin)
! sd1(8-1)    ar(01)-
! sd1(8-2)    (11)ar(23)-
! sd1(8-3)    ar(13)-c'(21)-
! sd1(8-4)    ar(23)-c'(11)-

use gugaci_global, only: iml, imr, ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, &
                         lsm_inn, mhlp, mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_sd1, w1_sd1
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: isma, iwal, iwar, iwdl, iwdr, lmi, lmij, lmj, lri, lrj, mpl, ni
real(kind=wp) :: w0sd1, w0sd2, w0sd3, w0sd4, w1sd1, w1sd3, w1sd4
integer(kind=iwp), external :: iwalk_ad

isma = Mul(iml,imr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if ((jml /= 1) .or. (jmr /= lmi)) cycle
  w0sd1 = w0_sd1(1)
  w1sd1 = w1_sd1(1)
  if (mod(norb_dz-lri,2) == 1) then
    w0sd1 = -w0sd1
    w1sd1 = -w1sd1
  end if
  iwdl = just(lri,lri)
  iwdr = jud(lri)
  do mpl=1,mhlp
    iwal = lpnew_lwei(mpl)
    iwar = lpnew_rwei(mpl)
    lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
    lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
  end do
  do mpl=1,mtype
    vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd1
    vplp_w1(mpl) = vplpnew_w1(mpl)*w1sd1
  end do
  if (lin == 6) then
    call ar_sd_ext_ar(26,lri,lrj,isma)
    call ar_sd_ext_rest(lri)
  end if
  if (lin == 13) then
    call ar_td_ext_ar(26,lri,lrj,isma)
    call ar_td_ext_rest(lri)
  end if
  if (lin == 23) then
    call ar_dv_ext_ar(26,isma,lri,lrj)   !ar_dv
  end if
end do

do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    if (lmij /= jml) cycle

    if (lmi == jmr) then
      iwdl = just(lrj,lri)
      iwdr = jud(lri)
      w0sd2 = w0_sd1(2)
      !w1sd2 = w0_sd1(2)
      ni = mod(norb_dz-lrj,2)
      if (ni == 1) then
        !w1sd2 = -w1sd2
        w0sd2 = -w0sd2
      end if
      ! sd1(8-2)    (11)ar(23)-
      do mpl=1,mtype
        vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd2
        vplp_w1(mpl) = vplpnew_w1(mpl)*w0sd2
      end do
      do mpl=1,mhlp
        iwal = lpnew_lwei(mpl)
        iwar = lpnew_rwei(mpl)
        lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
        lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
      end do
      if (lin == 6) then
        call ar_sd_ext_ar(25,lrj,lri,isma)
        call ar_sd_ext_rest(lrj)
      end if
      if (lin == 13) then
        call ar_td_ext_ar(25,lrj,lri,isma)
        call ar_td_ext_rest(lrj)
      end if
      if (lin == 23) then
        call ar_dv_ext_ar(25,isma,lrj,lri)   !ar_dv
      end if
    end if
    ! sd1(8-3)    ar(13)-c'(21)-
    ! sd1(8-4)    ar(23)-c'(11)-
    if (lmj /= jmr) cycle
    w0sd3 = w0_sd1(3)
    w1sd3 = w1_sd1(3)
    w0sd4 = w0_sd1(4)
    w1sd4 = w1_sd1(4)
    ni = mod(norb_dz-lri,2)
    if (ni == 0) then
      w0sd3 = -w0sd3
      w1sd3 = -w1sd3
      w0sd4 = -w0sd4
      w1sd4 = -w1sd4
    end if
    if (lmj == jmr) then
      iwdl = just(lrj,lri)
      iwdr = jud(lrj)
      do mpl=1,mtype
        vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd3
        vplp_w1(mpl) = vplpnew_w1(mpl)*w1sd3
      end do
      do mpl=1,mhlp
        iwal = lpnew_lwei(mpl)
        iwar = lpnew_rwei(mpl)
        lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
        lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
      end do
      if (lin == 6) then
        call ar_sd_ext_ar(29,lri,lrj,isma)
        call ar_sd_ext_rest(lri)
      end if
      if (lin == 13) then
        call ar_td_ext_ar(29,lri,lrj,isma)
        call ar_td_ext_rest(lri)
      end if
      if (lin == 23) then
        call ar_dv_ext_ar(29,isma,lri,lrj)   !ar_dv
      end if
      iwdl = just(lri,lrj)
      iwdr = jud(lrj)
      do mpl=1,mtype
        vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd4
        vplp_w1(mpl) = vplpnew_w1(mpl)*w1sd4
      end do
      do mpl=1,mhlp
        iwal = lpnew_lwei(mpl)
        iwar = lpnew_rwei(mpl)
        lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
        lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
      end do
      if (lin == 6) then
        call ar_sd_ext_ar(57,lri,lrj,isma)
        call ar_sd_ext_rest(lri)
      end if
      if (lin == 13) then
        call ar_td_ext_ar(57,lri,lrj,isma)
        call ar_td_ext_rest(lri)
      end if
      if (lin == 23) then
        call ar_dv_ext_ar(57,isma,lri,lrj)   !ar_dv
      end if
    end if
  end do
end do

return

end subroutine sdd_ar_act_c_sd_ext_sgt0

subroutine ttdd_ar_act_c_ttdd_ext_sgt1(lin)
! t1d1(15-1)  ar(13)-
! t1d1(15-1)  ar(13)-c'(11)-

use gugaci_global, only: iml, imr, ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, &
                         lsm_inn, mhlp, mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_t1d1, w1_t1d1
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: isma, iwal, iwar, iwdl, iwdr, lmi, lmj, lri, lrj, mpl
real(kind=wp) :: w0td1, w1td1
integer(kind=iwp), external :: iwalk_ad

isma = Mul(iml,imr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=norb_frz+1,lri-1
    lmj = lsm_inn(lrj)
    if ((jml /= Mul(lmi,lmj)) .or. (jmr /= lmj)) cycle
    w0td1 = w0_t1d1(1)
    w1td1 = w1_t1d1(1)
    if (mod(norb_dz-lri,2) == 1) then
      w0td1 = -w0td1
      w1td1 = -w1td1
    end if
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0td1
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1td1
    end do
    ! t1d1(15-1)  (11)ar(13)-
    iwdl = just(lrj,lri)
    iwdr = jud(lrj)
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    if (lin == 6) then
      call ar_sd_ext_ar(43,lri,lrj,isma)
      call ar_sd_ext_rest(lri)
    end if
    if (lin == 13) then
      call ar_td_ext_ar(43,lri,lrj,isma)
      call ar_td_ext_rest(lri)
    end if
    if (lin == 23) then
      call ar_dv_ext_ar(43,isma,lri,lrj)
    end if
  end do
  ! t1d1(15-1)  ar(13)-c'(11)-
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    if ((jml /= Mul(lmi,lmj)) .or. (jmr /= lmj)) cycle
    w0td1 = w0_t1d1(1)
    w1td1 = w1_t1d1(1)
    if (mod(norb_dz-lri,2) == 1) then
      w0td1 = -w0td1
      w1td1 = -w1td1
    end if
    do mpl=1,mtype
      vplp_w0(mpl) = -vplpnew_w0(mpl)*w0td1
      vplp_w1(mpl) = -vplpnew_w1(mpl)*w1td1
    end do
    iwdl = just(lri,lrj)
    iwdr = jud(lrj)
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    if (lin == 6) then
      call ar_sd_ext_ar(46,lri,lrj,isma)
      call ar_sd_ext_rest(lri)
    end if
    if (lin == 13) then
      call ar_td_ext_ar(46,lri,lrj,isma)
      call ar_td_ext_rest(lri)
    end if
    if (lin == 23) then
      call ar_dv_ext_ar(46,isma,lri,lrj)
    end if
  end do
end do

return

end subroutine ttdd_ar_act_c_ttdd_ext_sgt1

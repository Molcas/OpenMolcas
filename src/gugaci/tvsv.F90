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

subroutine tv_drt_ci_new()

use gugaci_global, only: idisk_array, idisk_lp, iml, ipae, ipael, jml, jmr, jpad, jpadl, jpadlr, linelp, lpblock_tv
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: jptyl, jptyr, lpb

!iltype = 3
!irtype = 1
!w0_d48 = Zero
!w1_d48 = v_sqthree

idisk_lp = idisk_array(6)

do lpb=1,lpblock_tv
  call read_lp()
  ipael = iml+9
  ipae = 1
  call get_jpty(jpadlr,jptyl,jptyr)
  call get_jp(jptyl,jml,jpadl,1)
  call get_jp(jptyr,jmr,jpad,1)
  !jmlr = Mul(jml,jmr)
  if (linelp <= 12) then
    call tv_ext_head_in_act()
  else
    call tv_ext_head_in_dbl()
  end if
end do

return

end subroutine tv_drt_ci_new

subroutine tv_ext_head_in_dbl()

use gugaci_global, only: ipae, ipael, jb_sys, jml, jmr, jpad, jpadl, jpadlr, jud, just, linelp, logic_dh, lp_lwei, lp_rwei, &
                         lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, nlg1, nlg2, norb_dz, norb_frz, vplp_w0, vplp_w1, &
                         vplpnew_w0, vplpnew_w1, w0_dv, w0_sd, w0_sv, w0_td, w1_sv, w1_tv
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, jmlr, lmi, lmij, lmj, lmk, lpok, lra, lri, lrj, lrk, mpl, ni
real(kind=wp) :: w0, w0sd1, w0sd2, w0sd4, w0sv2, w0td1, w1, w1sv2
integer(kind=iwp), external :: iwalk_ad

logic_dh = .true.
lpok = jpadlr
jmlr = Mul(jml,jmr)
select case (lpok)
  case (6)
    !===================================================================
    ! sd(6-1) ar(02)-    act: -b&r-  ext tv: -b^r-a^r
    ! sd(6-2) c(22)-ar(13)-
    ! sd(6-3) ar(13)-c'(22)-
    ! sd(6-4) ar(23)-c'(12)-
    if (linelp /= 18) return
    lra = nlg1
    if (jb_sys > 0) call sd_ar_act_br_sgt0(17,lra)
    do lri=norb_frz+1,norb_dz
      lmi = lsm_inn(lri)
      if (lmi /= jmlr) cycle
      w0sd1 = w0_sd(1)
      w0sd2 = w0_sd(2)
      w0sd4 = w0_sd(4)
      ni = mod(norb_dz-lri,2)
      if (ni == 1) w0sd1 = -w0sd1
      if (ni == 1) w0sd2 = -w0sd2
      if (ni == 1) w0sd4 = -w0sd4
      ! sd(6-1) ar(02)-         act: -br-  ext tv: -br-ar
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
        call ar_br_tv_ext_br_ar(lri,lra)
      end if

      ! sd(6-2) c(22)-ar(13)-
      do mpl=1,mtype
        vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd2
        vplp_w1(mpl) = vplpnew_w1(mpl)*w0sd2
      end do
      do lrk=norb_frz+1,lri-1
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
        call ar_br_tv_ext_br_ar(lri,lra)
      end do
      ! sd(6-4)   ar(23)-c'(12)-
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
        call ar_br_tv_ext_br_ar(lri,lra)
      end do
    end do

  case (8)
    !===================================================================
    ! sd1(8) ar- act -br-
    if (linelp /= 18) return
    lra = nlg1
    call sdd_ar_act_br_sgt0(17,lra)

  case (10)
    !===================================================================
    ! sv(10-1) ar(13)br(23)  act -c"-  tv_ext -br-ar
    ! sv(10-2) ar(23)br(13)  act -c"-  tv_ext -br-ar
    ! sv(10-3) dr(03)        act -c"-  tv_ext -br-ar
    if ((linelp /= 14) .or. (nlg2 /= 2)) return
    if (jb_sys > 0) then
      call sv_arbr_act_c_ext_stv_sgt0(17)
    end if
    do lri=norb_frz+1,norb_dz
      lmi = lsm_inn(lri)
      do lrj=lri,norb_dz
        lmj = lsm_inn(lrj)
        lmij = Mul(lmi,lmj)
        if (lmij /= jmlr) cycle
        !---------------------------------------------------------------
        w0sv2 = w0_sv(2)
        w1sv2 = w1_sv(2)
        ni = mod(lrj-lri,2)
        if (ni == 0) then
          w0sv2 = -w0sv2
          w1sv2 = -w1sv2
        end if
        iwdl = just(lri,lrj)
        iwdr = 0
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        if (lri /= lrj) then
          do mpl=1,mtype
            vplp_w0(mpl) = vplpnew_w0(mpl)*w0sv2
            vplp_w1(mpl) = vplpnew_w1(mpl)*w1sv2
          end do
          call ar_br_tv_ext_br_ar(lri,lrj)
        else
          ! sv(10-3) dr(03)
          do mpl=1,mtype
            vplp_w0(mpl) = vplpnew_w0(mpl)*w0_sv(3)
            vplp_w1(mpl) = vplpnew_w1(mpl)*w1_sv(3)
          end do
          call drr_tv_ext_br_ar(lri)
        end if
      end do
    end do

  case (13)
    !===================================================================
    ! td(13-1) (22)-ar(23)-         act: -b&r-  ext tv: -br-ar
    ! td(13-2) ar(23)-c'(22)-
    if (linelp /= 18) return
    lra = nlg1
    do lri=norb_frz+1,norb_dz
      lmi = lsm_inn(lri)
      if (lmi /= jmlr) cycle
      w0td1 = w0_td(1)
      ni = mod(norb_dz-lri,2)
      if (ni == 1) w0td1 = -w0td1
      ! td(13-1) (22)-ar(23)-        act: -br-  ext tv: -br-ar
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
        call ar_br_tv_ext_br_ar(lri,lra)
      end do
      ! td(13-1)   ar(23)-c'(22)-
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
        call ar_br_tv_ext_br_ar(lri,lra)
      end do
    end do

  case (15)
    !===================================================================
    ! t1d1(15) ar- act -br-
    if (linelp /= 18) return
    lra = nlg1
    call ttdd_ar_act_br_sgt1(17,lra)

  case (17)
    !===================================================================
    ! tv(17) ar(23)br(23) act -c"- tv_ext -br-ar
    if ((linelp /= 14) .or. (nlg2 /= 2)) return
    do lri=norb_frz+1,norb_dz-1
      lmi = lsm_inn(lri)
      do lrj=lri+1,norb_dz
        lmj = lsm_inn(lrj)
        lmij = Mul(lmi,lmj)
        if (lmij /= jmlr) cycle
        !---------------------------------------------------------------
        iwdl = just(lri,lrj)
        iwdr = 0
        ni = mod(lrj-lri,2)
        w1 = w1_tv
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
        call ar_br_tv_ext_br_ar(lri,lrj)
      end do
    end do

  case (18)
    !===================================================================
    ! t1v(18) ar-br- act -c"-
    if ((linelp /= 14) .or. (nlg2 /= 2)) return
    call ttv_arbr_act_c_stv_sgt1(17)

  case (23)
    !===================================================================
    ! dv(23-1)  ar(23)-    act: -b&r-  ext tv: -br-ar
    if (linelp /= 18) return
    lra = nlg1
    do lri=norb_frz+1,norb_dz
      lmi = lsm_inn(lri)
      if (lmi /= jml) cycle
      iwdl = jud(lri)
      iwdr = 0
      ni = mod(norb_dz-lri,2)
      w0 = w0_dv(1)
      if (ni == 1) w0 = -w0
      do mpl=1,mtype
        vplp_w0(mpl) = vplpnew_w0(mpl)*w0
        vplp_w1(mpl) = vplpnew_w1(mpl)*w0
      end do
      do mpl=1,mhlp
        iwal = lpnew_lwei(mpl)
        iwar = lpnew_rwei(mpl)
        lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
        lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
      end do
      call ar_br_tv_ext_br_ar(lri,lra)
    end do

  case (24)
    !===================================================================
    ! d1v(24-1) ar- act -br-
    if (linelp /= 18) return
    lra = nlg1
    call d1v_ar_act_bl_ext_ab_sgt0(17,lra)

  case default
end select

return

end subroutine tv_ext_head_in_dbl

subroutine tv_ext_head_in_act()

use gugaci_global, only: linelp, logic_dh, nlg1, nlg2
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: lra, lrai, lraj

logic_dh = .false.

if (linelp == 3) then
  ! line=3 a&r--b&r<-->b^r-a^r
  lrai = nlg1
  lraj = nlg2
  call ar_br_tv_ext_br_ar(lrai,lraj)
end if
if (linelp == 8) then
  ! line=8 d&rr<-->b^r-a^r
  lra = nlg1
  call drr_tv_ext_br_ar(lra)
end if

return

end subroutine tv_ext_head_in_act

subroutine sv_drt_ci_new()

use gugaci_global, only: idisk_array, idisk_lp, iml, ipae, ipael, jml, jmr, jpad, jpadl, jpadlr, linelp, logic_g13, logic_g36a, &
                         lpblock_sv, v_sqtwo, w0g13a, w0g36a, w1g36a
use Constants, only: Zero, One
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: jptyl, jptyr, lpb

w0g36a = -One
w1g36a = Zero
w0g13a = -v_sqtwo

!iltype = 4
!irtype = 1
logic_g36a = .true.
logic_g13 = .false.

idisk_lp = idisk_array(10)

do lpb=1,lpblock_sv
  call read_lp()
  ipael = iml+17
  ipae = 1
  if (iml == 1) logic_g13 = .true.
  call get_jpty(jpadlr,jptyl,jptyr)
  call get_jp(jptyl,jml,jpadl,1)
  call get_jp(jptyr,jmr,jpad,1)
  !jmlr = Mul(jml,jmr)
  if (linelp <= 12) then
    call sv_ext_head_in_act()
  else
    call sv_ext_head_in_dbl()
  end if
end do

return

end subroutine sv_drt_ci_new

subroutine sv_ext_head_in_dbl()

use gugaci_global, only: ipae, ipael, jb_sys, jml, jmr, jpad, jpadl, jpadlr, jud, just, linelp, logic_dh, lp_lwei, lp_rwei, &
                         lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, nlg1, nlg2, norb_dz, norb_frz, vplp_w0, vplp_w1, &
                         vplpnew_w0, vplpnew_w1, w0_dv, w0_sd, w0_sv, w0_td, w1_sv, w1_tv
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, jmlr, lmi, lmij, lmj, lmk, lpok, lra, lri, lrj, lrk, mpl, ni
real(kind=wp) :: w0, w0sd1, w0sd2, w0sd4, w0sv2, w0td1, w1, w1sv2
integer(kind=iwp), external :: iwalk_ad

logic_dh = .true.
lpok = jpadlr
jmlr = Mul(jml,jmr)
select case (lpok)
  case (6)
    !===================================================================
    ! sd(6-1) ar(02)-    act: -b&r-  ext tv: -b^r-a^r
    ! sd(6-2) c(22)-ar(13)-
    ! sd(6-3) ar(13)-c'(22)-
    ! sd(6-4) ar(23)-c'(12)-
    if (linelp /= 18) return
    lra = nlg1
    if (jb_sys > 0) call sd_ar_act_br_sgt0(10,lra)
    do lri=norb_frz+1,norb_dz
      lmi = lsm_inn(lri)
      if (lmi /= jmlr) cycle
      w0sd1 = w0_sd(1)
      w0sd2 = w0_sd(2)
      w0sd4 = w0_sd(4)
      ni = mod(norb_dz-lri,2)
      if (ni == 1) w0sd1 = -w0sd1
      if (ni == 1) w0sd2 = -w0sd2
      if (ni == 1) w0sd4 = -w0sd4
      ! sd(6-1) ar(02)-         act: -br-  ext tv: -br-ar
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
        call ar_br_sv_ext_br_ar(lri,lra)
      end if

      ! sd(6-2) c(22)-ar(13)-
      do mpl=1,mtype
        vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd2
        vplp_w1(mpl) = vplpnew_w1(mpl)*w0sd2
      end do
      do lrk=norb_frz+1,lri-1
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
        call ar_br_sv_ext_br_ar(lri,lra)
      end do
      ! sd(6-4)   ar(23)-c'(12)-
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
        call ar_br_sv_ext_br_ar(lri,lra)
      end do
    end do

  case (8)
    !===================================================================
    ! sd1(8) ar- act -br-
    if (linelp /= 18) return
    lra = nlg1
    call sdd_ar_act_br_sgt0(10,lra)

  case (10)
    !===================================================================
    ! sv(10-1) ar(13)br(23)  act -c"-  tv_ext -br-ar
    ! sv(10-2) ar(23)br(13)  act -c"-  tv_ext -br-ar
    ! sv(10-1) dr(03)        act -c"-  tv_ext -br-ar
    if ((linelp /= 14) .or. (nlg2 /= 2)) return
    if (jb_sys > 0) then
      call sv_arbr_act_c_ext_stv_sgt0(10)
    end if
    do lri=norb_frz+1,norb_dz
      lmi = lsm_inn(lri)
      do lrj=lri,norb_dz
        lmj = lsm_inn(lrj)
        lmij = Mul(lmi,lmj)
        if (lmij /= jmlr) cycle
        !---------------------------------------------------------------
        w0sv2 = w0_sv(2)
        w1sv2 = w1_sv(2)
        ni = mod(lrj-lri,2)
        if (ni == 0) then
          w0sv2 = -w0sv2
          w1sv2 = -w1sv2
        end if
        iwdl = just(lri,lrj)
        iwdr = 0
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        if (lri /= lrj) then
          do mpl=1,mtype
            vplp_w0(mpl) = vplpnew_w0(mpl)*w0sv2
            vplp_w1(mpl) = vplpnew_w1(mpl)*w1sv2
          end do
          call ar_br_sv_ext_br_ar(lri,lrj)
        else
          ! sv(10-1) dr(03)
          do mpl=1,mtype
            vplp_w0(mpl) = vplpnew_w0(mpl)*w0_sv(3)
            vplp_w1(mpl) = vplpnew_w1(mpl)*w1_sv(3)
          end do
          call drr_sv_ext_br_ar(lri)
        end if
      end do
    end do

  case (13)
    !===================================================================
    ! td(13-1) (22)-ar(23)-         act: -b&r-  ext tv: -br-ar
    ! td(13-2) ar(23)-c'(22)-
    if (linelp /= 18) return
    lra = nlg1
    do lri=norb_frz+1,norb_dz
      lmi = lsm_inn(lri)
      if (lmi /= jmlr) cycle
      w0td1 = w0_td(1)
      ni = mod(norb_dz-lri,2)
      if (ni == 1) w0td1 = -w0td1
      ! td(13-1) (22)-ar(23)-        act: -br-  ext tv: -br-ar
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
        call ar_br_sv_ext_br_ar(lri,lra)
      end do
      ! td(13-1)   ar(23)-c'(22)-
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
        call ar_br_sv_ext_br_ar(lri,lra)
      end do
    end do

  case (15)
    !===================================================================
    ! t1d1(15) ar- act -br-
    if (linelp /= 18) return
    lra = nlg1
    call ttdd_ar_act_br_sgt1(10,lra)

  case (17)
    !===================================================================
    ! tv(17) ar(23)br(23) act -c"- tv_ext -br-ar
    if ((linelp /= 14) .or. (nlg2 /= 2)) return
    do lri=norb_frz+1,norb_dz-1
      lmi = lsm_inn(lri)
      do lrj=lri+1,norb_dz
        lmj = lsm_inn(lrj)
        lmij = Mul(lmi,lmj)
        if (lmij /= jmlr) cycle
        !---------------------------------------------------------------
        iwdl = just(lri,lrj)
        iwdr = 0
        ni = mod(lrj-lri,2)
        w1 = w1_tv
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
        call ar_br_sv_ext_br_ar(lri,lrj)
      end do
    end do

  case (18)
    !===================================================================
    ! t1v(18) ar-br- act -c"-
    if (linelp /= 14) return
    call ttv_arbr_act_c_stv_sgt1(10)

  case (23)
    !===================================================================
    ! dv(23-1)  ar(23)-    act: -b&r-  ext tv: -br-ar
    if (linelp /= 18) return
    lra = nlg1
    do lri=norb_frz+1,norb_dz
      lmi = lsm_inn(lri)
      if (lmi /= jml) cycle
      iwdl = jud(lri)
      iwdr = 0
      ni = mod(norb_dz-lri,2)
      w0 = w0_dv(1)
      if (ni == 1) w0 = -w0
      do mpl=1,mtype
        vplp_w0(mpl) = vplpnew_w0(mpl)*w0
        vplp_w1(mpl) = vplpnew_w1(mpl)*w0
      end do
      do mpl=1,mhlp
        iwal = lpnew_lwei(mpl)
        iwar = lpnew_rwei(mpl)
        lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
        lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
      end do
      call ar_br_sv_ext_br_ar(lri,lra)
    end do

  case (24)
    !===================================================================
    ! d1v(24-1) ar- act -br-
    if (linelp /= 18) return
    lra = nlg1
    call d1v_ar_act_bl_ext_ab_sgt0(10,lra)

  case default
end select

return

end subroutine sv_ext_head_in_dbl

subroutine sv_ext_head_in_act()

use gugaci_global, only: linelp, logic_dh, nlg1, nlg2
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: lra, lrai, lraj

logic_dh = .false.

if (linelp == 3) then
  ! line=3 a&r--b&r<-->b^r-a^r
  lrai = nlg1
  lraj = nlg2
  call ar_br_sv_ext_br_ar(lrai,lraj)
end if
if (linelp == 8) then
  ! line=8 d&rr<-->b^r-a^r
  lra = nlg1
  call drr_sv_ext_br_ar(lra)
end if

return

end subroutine sv_ext_head_in_act

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

subroutine dv_drt_ci_new()

use gugaci_global, only: idisk_array, idisk_lp, iml, imr, ipae, ipael, jml, jmr, jpad, jpadl, jpadlr, linelp, lpblock_dv
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: jptyl, jptyr, lpb

call external_space_plpmode_value_dv()

idisk_lp = idisk_array(2)
do lpb=1,lpblock_dv
  call read_lp()
  ipael = iml+1
  ipae = 1
  call get_jpty(jpadlr,jptyl,jptyr)
  call get_jp(jptyl,jml,jpadl,1)
  call get_jp(jptyr,jmr,jpad,1)
  !jmlr = Mul(jml,jmr)
  call gsd_determine_extarmode_paras(iml,imr,.false.)
  if (linelp <= 12) then
    call dv_ext_head_in_act()
  else
    call dv_ext_head_in_dbl()
  end if
end do

return

end subroutine dv_drt_ci_new

subroutine dv_ext_head_in_dbl()

use gugaci_global, only: iml, intind_ijka, ipae, ipael, jb_sys, jml, jmr, jpad, jpadl, jpadlr, jud, just, linelp, logic_dh, &
                         lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, ngw2, ngw3, nlg1, nlg2, norb_dz, &
                         norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_dv, w0_sd, w0_sv, w0_td, w0_vv, w1_sd, w1_sv, &
                         w1_td, w1_tv
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: ijk, imap_1, intpos, isma, iwal, iwar, iwdl, iwdr, jk, jmlr, lmd, lmi, lmij, lmj, lpok, lra, lrd, lri, lrj, &
                     lrk, mpl, ni
real(kind=wp) :: w0, w0sd1, w0sd11, w0sd12, w0sd14, w0sd16, w0sd2, w0sd4, w0sd5, w0sd8, w0sd9, w0sv2, w0td1, w0td2, w0td3, w0td4, &
                 w0td5, w1sd11, w1sd12, w1sd5, w1sd8, w1sd9, w1sv2, w1td2, w1td3, w1td4, w1tv
integer(kind=iwp), external :: iwalk_ad

logic_dh = .true.
isma = iml
lpok = jpadlr
jmlr = Mul(jml,jmr)
select case (lpok)
  case (1)
    !===================================================================
    ! ss(1)     act -b^l-
    if (linelp /= 15) return
    lra = nlg1
    if (nlg2 == 1) then
      call ss_drl_act_bl(23,lra)
      if (jb_sys > 0) then
        call ss_drl_act_bl_sgt0(23,lra)
      end if
    else
      call ss_arbl_act_bl(23,lra)
      if (jb_sys > 0) then
        call ss_arbl_act_bl_sgt0(23,lra)
        call ss_s_drl_act_bl_sgt0(23,lra)
      end if
    end if

  case (2)
    !===================================================================
    ! st(2)      act -b^l
    if ((linelp /= 15) .or. (nlg2 /= 2)) return
    lra = nlg1
    call st_arbl_act_bl(23,lra)
    if (jb_sys > 0) call st_arbl_act_bl_sgt0(23,lra)
    if (jml /= jmr) return
    call st_drl_act_bl(23,lra)
    if (jb_sys > 0) call st_drl_act_bl_sgt0(23,lra)

  case (3)
    !===================================================================
    ! ts(3)   act -b^l ............................................
    if (linelp /= 15) return
    lra = nlg1
    if (nlg2 /= 2) return
    call ts_arbl_act_bl(23,lra)
    if (jb_sys > 0) call ts_arbl_act_bl_sgt0(23,lra)

  case (4)
    !===================================================================
    ! stt(4) act-b^l ..............................................
    if (linelp /= 15) return
    lra = nlg1
    if (nlg2 == 1) return
    call stt_arbl_act_bl_sgt1(23,lra)

  case (5)
    !===================================================================
    ! tts(5) act-b^l ..............................................
    if ((linelp /= 15) .or. (nlg2 /= 2)) return
    lra = nlg1
    call tts_arbl_act_bl_sgt1(23,lra)
    call tts_drl_act_bl_sgt1(23,lra)

  case default ! (6)
    !===================================================================
    ! sd(6-1) a&r(02)-
    ! sd(6-2) c(22)a&(13)-
    ! sd(6-3) a&r(13)c'(22)-
    ! sd(6-4) a&r(23)c'(12)-
    ! sd(6-5) a&r(23)b&r(13)b^r(32)
    ! sd(6-6) a&r(13)b&r(23)b^r(32)
    ! sd(6-7) a&r(13)b&l(32)b^l(23)
    ! sd(6-8) a&r(23)b&l(32)b^l(13)
    ! sd(6-9) d&r&r(03)b^r(32)
    ! sd(6-10) d&r&l(12)b^l(23)
    ! sd(6-11) d&r&l(22)b^l(13)
    ! sd(6-12) d&r&l(33)b^l(02)
    ! sd(6-13) (22)d&r&l(33)b^l(13)
    ! sd(6-14) d&r&l(33)c"(22)b^l(13)
    ! sd(6-15) d&r&l(33)b^l(13)c'(22)
    ! sd(6-16) d&r&l(33)b^l(23)c'(12)
    if (linelp == 13) then
      if (jb_sys > 0) call sd_adb_act_c_ext_ar(23)
      do lri=norb_frz+1,norb_dz
        lmi = lsm_inn(lri)
        if (lmi /= jmlr) cycle
        w0sd1 = w0_sd(1)
        w0sd9 = w0_sd(9)
        w1sd9 = w1_sd(9)
        w0sd12 = w0_sd(12)
        w1sd12 = w1_sd(12)
        ni = mod(norb_dz-lri,2)
        if (ni == 1) w0sd1 = -w0sd1
        if (ni == 1) w0sd9 = -w0sd9
        if (ni == 1) w1sd9 = -w1sd9
        if (ni == 1) w0sd12 = -w0sd12
        if (ni == 1) w1sd12 = -w1sd12
        if ((jml == 1) .and. (lmi == jmr)) then
          iwdl = just(lri,lri)
          iwdr = jud(lri)
          do mpl=1,mhlp
            iwal = lpnew_lwei(mpl)
            iwar = lpnew_rwei(mpl)
            lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
          end do

          ! sd(6-1) a&r(02)-
          do mpl=1,mtype
            vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd1
            vplp_w1(mpl) = vplpnew_w1(mpl)*w0sd1
          end do
          call ar_dv_ext_ar(26,isma,lri,lrj)   !ar_dv

          ! sd(6-12) d&r&l(33)b^l(02)
          do mpl=1,mtype
            vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd12
            vplp_w1(mpl) = vplpnew_w1(mpl)*w1sd12
          end do
          do lrk=1,lri-1
            call drl_bl_ext_ar_new(23,lrk,lri)
          end do

          ! sd(6-9) d&r&r(03)b^r(32)
          do lrk=norb_frz+1,lri-1
            iwdl = just(lrk,lrk)
            iwdr = jud(lri)
            do mpl=1,mhlp
              iwal = lpnew_lwei(mpl)
              iwar = lpnew_rwei(mpl)
              lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
            end do
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd9
              vplp_w1(mpl) = vplpnew_w1(mpl)*w1sd9
            end do
            call drr_br_ext_ar(23,lrk,lri)
          end do
        end if
      end do

      do lri=norb_frz+1,norb_dz
        lmi = lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj = lsm_inn(lrj)
          lmij = Mul(lmi,lmj)
          if (lmij /= jml) cycle

          if (lmi == jmr) then
            iwdl = just(lri,lrj)
            iwdr = jud(lri)
            w0sd2 = w0_sd(2)
            w0sd11 = w0_sd(11)
            w1sd11 = w1_sd(11)
            w0sd14 = w0_sd(14)
            ni = mod(norb_dz-lrj,2)
            if (ni == 1) then
              w0sd2 = -w0sd2
              w0sd11 = -w0sd11
              w1sd11 = -w1sd11
              w0sd14 = -w0sd14
            end if
            ! sd(6-2) c(22)-a&r(13)-
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd2
            end do
            do mpl=1,mhlp
              iwal = lpnew_lwei(mpl)
              iwar = lpnew_rwei(mpl)
              lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
            end do
            call ar_dv_ext_ar(25,isma,lrj,lri) !ar_dv
            ! sd(6-11) d&r&l(22)b^l(13)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd11
              vplp_w1(mpl) = vplpnew_w1(mpl)*w1sd11
            end do
            call drl_bl_ext_ar_new(23,lri,lrj)
            ! sd(6-14) (22)d&r&l(33)b^l(13)
            ! sd(6-14) d&r&l(33)c"(22)b^l(13)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd14
              vplp_w1(mpl) = Zero
            end do
            do lrk=1,lrj-1
              if (lrk == lri) cycle
              call drl_bl_ext_ar_new(23,lrk,lrj)
            end do
          end if
          ! sd(6-4) a&r(23)-c'(12)-
          w0sd4 = w0_sd(4)
          w0sd16 = w0_sd(16)
          ni = mod(norb_dz-lri,2)
          if (ni == 0) then
            w0sd4 = -w0sd4
            w0sd16 = -w0sd16
          end if
          if (lmj == jmr) then
            iwdl = just(lri,lrj)
            iwdr = jud(lrj)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd4
            end do
            do mpl=1,mhlp
              iwal = lpnew_lwei(mpl)
              iwar = lpnew_rwei(mpl)
              lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
            end do
            call ar_dv_ext_ar(28,isma,lri,lrj)   !ar_dv
            ! sd(6-16) d&r&l(33)b^l(23)c'(12)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd16
              vplp_w1(mpl) = Zero
            end do
            do lrk=1,lri-1
              call drl_bl_ext_ar_new(23,lrk,lri)
            end do
          end if
          ! sd(6-5) a&r(23)b&r(13)b^r(32)
          do lrd=lrj+1,norb_dz
            lmd = lsm_inn(lrd)
            if (lmd /= jmr) cycle

            w0sd5 = w0_sd(5)
            w1sd5 = w1_sd(5)
            ni = mod(lrj-lri+norb_dz-lrd,2)
            if (ni == 0) then
              w0sd5 = -w0sd5
              w1sd5 = -w1sd5
            end if
            iwdl = just(lri,lrj)
            iwdr = jud(lrd)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd5
              vplp_w1(mpl) = vplpnew_w1(mpl)*w1sd5
            end do
            do mpl=1,mhlp
              iwal = lpnew_lwei(mpl)
              iwar = lpnew_rwei(mpl)
              lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
            end do
            ijk = lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lrd-norb_frz)
            intpos = intind_ijka(ijk)
            call ar_br_br_ext_ar_new(23,intpos,isma)
          end do
          ! sd(6-8) a&r(23)b&l(32)b^l(13)
          do lrd=lri+1,lrj-1
            lmd = lsm_inn(lrd)
            if (lmd /= jmr) cycle
            iwdl = just(lri,lrj)
            iwdr = jud(lrd)
            w0sd8 = w0_sd(8)
            w1sd8 = w1_sd(8)
            ni = mod(lrj-lri+norb_dz-lrd,2)
            if (ni == 0) then
              w0sd8 = -w0sd8
              w1sd8 = -w1sd8
            end if
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd8
              vplp_w1(mpl) = vplpnew_w1(mpl)*w1sd8
            end do
            do mpl=1,mhlp
              iwal = lpnew_lwei(mpl)
              iwar = lpnew_rwei(mpl)
              lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
            end do
            ijk = lri-norb_frz+ngw2(lrd-norb_frz)+ngw3(lrj-norb_frz)
            intpos = intind_ijka(ijk)
            call ar_bl_bl_ext_ar_new(23,intpos,isma,1)
          end do
        end do
      end do
    end if
    !===================================================================
    ! sd(6)    act -br-br-
    if (linelp == 20) then
      jk = nlg1
      call sd_ar_act_brbr(23,jk)
      if (jb_sys > 0) call sd_ar_act_brbr_sgt0(23,jk)
    end if
    ! sd(6)  act -bl-bl-
    if (linelp == 22) then
      jk = nlg1
      call sd_ar_act_blbl(23,jk)
      if (jb_sys > 0) call sd_ar_act_blbl_sgt0(23,jk)
    end if

  case (8)
    !===================================================================
    ! sdd(8) act-b^l -c'- -c"-..........................................
    if (linelp == 13) then
      lra = nlg1
      call sdd_drlbl_act_c_sgt0(23)
      call sdd_drrbr_act_c_sgt0(23)
      call sdd_abb_act_c_sgt0(23)
      call sdd_ar_act_c_sd_ext_sgt0(23)
    else
      lra = nlg1
      if (linelp == 20) call sdd_ar_act_brbr_sgt0(23,lra)
      if (linelp == 22) call sdd_ar_act_blbl_sgt0(23,lra)
    end if

  case (9)
    !===================================================================
    ! dds(9) act -c'- ..............................................
    if (linelp == 13) then
      lra = nlg1
      call dds_abb_act_c_sgt0(23)
    end if

  case (10)
    !===================================================================
    ! sv(10)     act -b^r-
    ! sv(10-1) ar(13)-br(23)-
    ! sv(10-2) ar(23)-br(13)-
    ! sv(10-3) drr(03)-
    if ((linelp /= 16) .or. (nlg2 /= 2)) return
    lra = nlg1
    if (jb_sys > 0) then
      call sv_arbr_act_br_sgt0(23,lra)
    end if
    do lri=norb_frz+1,norb_dz
      lmi = lsm_inn(lri)
      do lrj=lri,norb_dz
        lmj = lsm_inn(lrj)
        lmij = Mul(lmi,lmj)
        if (lmij /= jml) cycle
        w0sv2 = w0_sv(2)
        w1sv2 = w1_sv(2)
        ni = mod(lrj-lri,2)
        if (ni == 0) then
          w0sv2 = -w0sv2
          w1sv2 = -w1sv2
        end if
        !---------------------------------------------------------------
        iwdl = just(lri,lrj)
        iwdr = 0
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        if (lri /= lrj) then
          ! sv(10-2) ar(23)-br(13)-
          do mpl=1,mtype
            vplp_w0(mpl) = vplpnew_w0(mpl)*w0sv2
            vplp_w1(mpl) = vplpnew_w1(mpl)*w1sv2
          end do
          ijk = lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz)
          intpos = intind_ijka(ijk)
          call ar_br_br_ext_ar_new(23,intpos,isma)
        else
          ! sv(10-3) drr(03)-
          do mpl=1,mtype
            vplp_w0(mpl) = vplpnew_w0(mpl)*w0_sv(3)
            vplp_w1(mpl) = vplpnew_w1(mpl)*w1_sv(3)
          end do
          call drr_br_ext_ar(23,lri,lra)
        end if
      end do
    end do

  case (11)
    !===================================================================
    ! tt(11)   act -b^l
    if (linelp /= 15) return
    lra = nlg1
    if (nlg2 == 1) then
      call tt_drl_act_bl(23,lra)
    else
      call tt_arbl_act_bl(23,lra)
    end if

  case (12)
    !===================================================================
    ! tttt(12)   act -b^l
    if (linelp /= 15) return
    lra = nlg1
    if (nlg2 == 1) then
      call tttt_drl_act_bl_sgt1(23,lra)
    else
      call tttt_arbl_act_bl_sgt1(23,lra)
    end if

  case (13)
    !===================================================================
    ! td(13)                                  act -c'-
    ! td(13-1) (22)a&(23)
    ! td(13-1) a&(23)c'(22)
    ! td(13-2) a&(23)b&r(23)b^r(32)
    ! td(13-3) a&(23)b&l(32)b^l(23)
    ! td(13-4) d&r&l(22)b^l(23)
    ! td(13-5) (22)d&&l(33)b^l(23)
    ! td(13-5) d&rl(33)c"(22)b^l(23)
    ! td(13-5) d&rl(33)b^l(23)c'(22)
    if (linelp == 13) then
      do lri=norb_frz+1,norb_dz
        lmi = lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj = lsm_inn(lrj)
          lmij = Mul(lmi,lmj)
          if (lmij /= jml) cycle
          w0td1 = w0_td(1)
          w0td4 = w0_td(4)
          w1td4 = w1_td(4)
          w0td5 = w0_td(5)
          ni = mod(norb_dz-lrj,2)
          if (ni == 1) then
            w0td1 = -w0td1
            w0td4 = -w0td4
            w1td4 = -w1td4
            w0td5 = -w0td5
          end if
          iwdl = just(lri,lrj)
          !-------------------------------------------------------------
          ! td(13-1) (22)a&(23)
          if (lmi == jmr) then
            iwdr = jud(lri)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0td1
            end do
            do mpl=1,mhlp
              iwal = lpnew_lwei(mpl)
              iwar = lpnew_rwei(mpl)
              lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
            end do
            call ar_dv_ext_ar(43,isma,lrj,lri)
            ! td(13-4) d&r&l(22)b^l(23)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0td4
              vplp_w1(mpl) = vplpnew_w1(mpl)*w1td4
            end do
            call drl_bl_ext_ar_new(23,lri,lrj)
            ! td(13-5) (22)d&rl(33)b^l(23)
            ! td(13-5) d&rl(33)c"(22)b^l(23)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0td5
              vplp_w1(mpl) = Zero
            end do
            do lrk=1,lrj-1
              if (lrk == lri) cycle
              call drl_bl_ext_ar_new(23,lrk,lrj)
            end do
          end if
          !-------------------------------------------------------------
          ! td(13-1) a&(23)c'(22)
          if (lmj == jmr) then
            iwdr = jud(lrj)
            w0td1 = w0_td(1)
            w0td5 = w0_td(5)
            ni = mod(norb_dz-lri,2)
            if (ni == 0) then
              w0td1 = -w0td1
              w0td5 = -w0td5
            end if
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0td1
            end do
            do mpl=1,mhlp
              iwal = lpnew_lwei(mpl)
              iwar = lpnew_rwei(mpl)
              lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
            end do
            call ar_dv_ext_ar(46,isma,lri,lrj) !ar_dv
            ! td(13-5) d&rl(33)b^l(23)c'(22)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0td5
              vplp_w1(mpl) = Zero
            end do
            do lrk=1,lri-1
              call drl_bl_ext_ar_new(23,lrk,lri)
            end do
          end if
          !-------------------------------------------------------------
          ! td(13-2) a&(23)b&r(23)b^r(32)
          do lrd=lrj+1,norb_dz
            lmd = lsm_inn(lrd)
            if (lmd /= jmr) cycle
            iwdr = jud(lrd)
            ijk = lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lrd-norb_frz)
            intpos = intind_ijka(ijk)
            w0td2 = w0_td(2)
            w1td2 = w1_td(2)
            ni = mod(lrj-lri+norb_dz-lrd,2)
            if (ni == 0) then
              w0td2 = -w0td2
              w1td2 = -w1td2
            end if
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0td2
              vplp_w1(mpl) = vplpnew_w1(mpl)*w1td2
            end do
            do mpl=1,mhlp
              iwal = lpnew_lwei(mpl)
              iwar = lpnew_rwei(mpl)
              lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
            end do
            call ar_br_br_ext_ar_new(23,intpos,isma)
          end do
          !-------------------------------------------------------------
          ! td(13-3) a&(23)b&l(32)b^l(23)
          do lrd=lri+1,lrj-1
            lmd = lsm_inn(lrd)
            if (lmd /= jmr) cycle
            iwdr = jud(lrd)
            w0td3 = w0_td(3)
            w1td3 = w1_td(3)
            ni = mod(lrj-lri+norb_dz-lrd,2)
            if (ni == 0) then
              w0td3 = -w0td3
              w1td3 = -w1td3
            end if
            ijk = lri-norb_frz+ngw2(lrd-norb_frz)+ngw3(lrj-norb_frz)
            intpos = intind_ijka(ijk)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0td3
              vplp_w1(mpl) = vplpnew_w1(mpl)*w1td3
            end do
            do mpl=1,mhlp
              iwal = lpnew_lwei(mpl)
              iwar = lpnew_rwei(mpl)
              lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
            end do
            call ar_bl_bl_ext_ar_new(23,intpos,isma,1)
          end do
          !-------------------------------------------------------------
        end do
      end do
    end if
    !===================================================================
    ! td(13) act -br-br-
    if (linelp == 20) then
      jk = nlg1
      call td_ar_act_brbr(23,jk)
    end if
    ! td(13)  act -bl-bl-
    if (linelp == 22) then
      jk = nlg1
      call td_ar_act_blbl(23,jk)
    end if

  case (15)
    !===================================================================
    ! ttdd(15)    act -c'-   -br-br- -bl-bl
    lra = nlg1
    if (linelp == 13) then
      call ttdd_drlbl_act_c_sgt1(23)
      call ttdd_abb_act_c_sgt1(23)             ! ????
      call ttdd_ar_act_c_ttdd_ext_sgt1(23)
    end if
    if (linelp == 20) then
      call ttdd_ar_act_brbr_sgt0(23,lra)
    end if
    if (linelp == 22) then
      call ttdd_ar_act_blbl_sgt0(23,lra)
    end if

  case (17)
    !===================================================================
    ! tv(17)    act -b^r-
    if ((linelp /= 16) .or. (nlg2 /= 2)) return
    lra = nlg1
    do lri=norb_frz+1,norb_dz-1
      lmi = lsm_inn(lri)
      do lrj=lri+1,norb_dz
        lmj = lsm_inn(lrj)
        lmij = Mul(lmi,lmj)
        if (lmij /= jml) cycle

        w1tv = w1_tv
        ni = mod(lrj-lri,2)
        if (ni == 0) w1tv = -w1tv

        !---------------------------------------------------------------
        ! tv(17) ar(23)-br(23)-
        iwdl = just(lri,lrj)
        iwdr = 0
        do mpl=1,mtype
          vplp_w0(mpl) = Zero
          vplp_w1(mpl) = vplpnew_w1(mpl)*w1tv
        end do
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        ijk = lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz)
        intpos = intind_ijka(ijk)
        call ar_br_br_ext_ar_new(23,intpos,isma)
      end do
    end do

  case (18)
    !===================================================================
    ! ttv(19) act -b^r- ................................................
    if ((linelp /= 16) .or. (nlg2 /= 2)) return
    lra = nlg1
    call ttv_arbr_act_c_sgt1(23,lra)

  case (19)
    !===================================================================
    ! dd(19) act -b^l- .................................................
    if (linelp /= 15) return
    lra = nlg1
    if (nlg2 == 1) then
      call dd_drl_act_bl(23,lra)
    else
      call dd_arbl_act_bl(23,lra)
    end if

  case (20)
    !===================================================================
    ! dddd(19) act -b^l- -br- ..........................................
    lra = nlg1
    if (linelp == 15) then
      if (nlg2 == 2) then
        call dddd_arbl_act_bl_sgt0(23,lra)
      else
        call dddd_drl_act_bl_sgt0(23,lra)
      end if
    end if

  case (21)
    !===================================================================
    ! dd1(21) act -b^l-  ...............................................
    if ((linelp /= 15) .or. (nlg2 /= 2)) return
    lra = nlg1
    call dd1_arbl_act_bl_sgt0(23,lra)

  case (22)
    !===================================================================
    ! d1d(21) act -b^l-  ...............................................
    lra = nlg1
    if ((linelp == 15) .and. (nlg2 == 2)) then
      call d1d_arbl_act_bl_sgt0(23,lra)
      call d1d_drl_act_bl_sgt0(23,lra)
    end if

  case (23)
    !===================================================================
    ! dv(23) act -c'-..................................................
    ! dv(23-1) ar(23)-
    ! dv(23-2) drl(33)-bl(23)-
    if (linelp == 13) then
      imap_1 = 0
      do lrd=norb_frz+1,norb_dz
        lmd = lsm_inn(lrd)
        if (lmd /= jml) cycle
        ni = mod(norb_dz-lrd,2)

        ! dv(23-1) ar(23)-
        w0 = w0_dv(1)
        if (ni == 1) w0 = -w0_dv(1)
        iwdl = jud(lrd)
        iwdr = 0
        imap_1 = imap_1+1
        do mpl=1,mtype
          vplp_w0(mpl) = vplpnew_w0(mpl)*w0
        end do
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        call ar_dv_ext_ar(51,isma,lrd,0)   !ar_dv
        ! dv(23-2) drl(33)-bl(23)-
        w0 = w0_dv(2)
        if (ni == 1) w0 = -w0_dv(2)
        do mpl=1,mtype
          vplp_w0(mpl) = vplpnew_w0(mpl)*w0
          vplp_w1(mpl) = Zero
        end do
        do lrk=1,lrd-1
          call drl_bl_ext_ar_new(23,lrk,lrd)
        end do
      end do
    end if
    !===================================================================
    ! dv(23) act -br-br-................................................
    if (linelp == 20) then
      jk = nlg1
      do lrd=norb_frz+1,norb_dz
        lmd = lsm_inn(lrd)
        if (lmd /= jml) cycle
        ni = mod(norb_dz-lrd,2)
        w0 = w0_dv(1)
        if (ni == 1) w0 = -w0_dv(1)
        !..............................................................
        ! dv(23-1) ar(23)-
        iwdl = jud(lrd)
        iwdr = 0
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
        ijk = lrd-norb_frz+jk
        intpos = intind_ijka(ijk)
        call ar_br_br_ext_ar_new(23,intpos,isma)
      end do
    end if
    ! dv(23) act -bl-bl-................................................
    if (linelp /= 22) return
    jk = nlg1
    do lrd=norb_frz+1,norb_dz
      lmd = lsm_inn(lrd)
      if (lmd /= jml) cycle
      ni = mod(norb_dz-lrd,2)
      w0 = w0_dv(1)
      if (ni == 1) w0 = -w0_dv(1)
      ! dv(23-1) ar(23)-
      iwdl = jud(lrd)
      iwdr = 0
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
      ijk = lrd-norb_frz+jk
      intpos = intind_ijka(ijk)
      call ar_bl_bl_ext_ar_new(23,intpos,isma,1)
    end do

  case (24)
    !===================================================================
    ! d1v(24) act -bl-bl -br-br- -c'- ..................................
    lra = nlg1
    if (linelp == 13) then
      call d1v_drl_bl_act_c_sgt0(23)    ! drl-bl- ar-
    end if
    if (linelp == 22) then
      call d1v_ar_act_blbl_sgt0(23,lra)
    end if
    if (linelp == 20) then
      call d1v_ar_act_brbr_sgt0(23,lra)
    end if

  case (25)
    !===================================================================
    ! vv(25) act -b^l- .................................................
    if (linelp /= 15) return
    lra = nlg1
    if (nlg2 == 1) then
      iwdl = 0
      iwdr = 0
      do mpl=1,mtype
        vplp_w0(mpl) = vplpnew_w0(mpl)*w0_vv
        vplp_w1(mpl) = Zero
      end do
      do mpl=1,mhlp
        iwal = lpnew_lwei(mpl)
        iwar = lpnew_rwei(mpl)
        lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
        lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
      end do
      call drl_bl_sum_ar_new(23,0,0,lra)
      !do lrk=1,norb_dz
      !  call drl_bl_ext_ar_new(23,lrk,lra)
      !end do
    end if

  case (7,14,16,26)
end select

return

end subroutine dv_ext_head_in_dbl

subroutine dv_ext_head_in_act()

use gugaci_global, only: iml, linelp, logic_dh, nlg1, nlg2
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: intpos, isma, lri, lrj

logic_dh = .false.
lri = nlg1
lrj = nlg2
intpos = nlg1
isma = iml

select case (linelp)
  case default ! (1)
    ! line=1 a&r<-->a^r
    call ar_dv_ext_ar(51,isma,lri,lrj)      ! ar_dv_ext_ar(100_???
  case (4)
    ! line=4 a&r--b&r--b^r<-->a^r
    call ar_br_br_ext_ar_new(23,intpos,isma)
  case (7)
    ! line=7 a&r--b&l--b^l<-->a^r
    call ar_bl_bl_ext_ar_new(23,intpos,isma,1)
  case (10)
    ! line=10 d&rr--b^r<-->a^r
    call drr_br_ext_ar(23,lri,lrj)
  case (12)
    ! line=12 d&rl--b^l<-->a^r
    call drl_bl_ext_ar_new(23,lri,lrj)
  case (2:3,5:6,8:9,11)
end select

return

end subroutine dv_ext_head_in_act

subroutine sd_drt_ci_new_den()

use gugaci_global, only: idisk_array, idisk_lp, iml, imr, ipae, ipael, jml, jmr, jpad, jpadl, jpadlr, linelp, lpblock_sd
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: jptyl, jptyr, lpb

!write(u6,*) '  sd_wyb'

!logic_sd = .true.
call external_space_plpmode_value_sd()

idisk_lp = idisk_array(11)

do lpb=1,lpblock_sd
  call read_lp()
  ipael = iml+17
  ipae = imr+1
  call get_jpty(jpadlr,jptyl,jptyr)
  call get_jp(jptyl,jml,jpadl,1)
  call get_jp(jptyr,jmr,jpad,1)
  !jmlr = Mul(jml,jmr)
  call gsd_determine_extarmode_paras(iml,imr,.true.)
  if (linelp <= 12) then
    call sd_ext_head_in_act()
  else
    call sd_ext_head_in_dbl()
  end if
end do

return

end subroutine sd_drt_ci_new_den

subroutine sd_drt_ci_new()

use gugaci_global, only: idisk_array, idisk_lp, iml, imr, ipae, ipael, jml, jmr, jpad, jpadl, jpadlr, linelp, lpblock_sd
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: jptyl, jptyr, lpb

!write(u6,*) '  sd_wyb'

!logic_sd = .true.
call external_space_plpmode_value_sd()

idisk_lp = idisk_array(11)

do lpb=1,lpblock_sd
  call read_lp()
  ipael = iml+17
  ipae = imr+1
  call get_jpty(jpadlr,jptyl,jptyr)
  call get_jp(jptyl,jml,jpadl,1)
  call get_jp(jptyr,jmr,jpad,1)
  !jmlr = Mul(jml,jmr)
  call gsd_determine_extarmode_paras(iml,imr,.true.)
  if (linelp <= 12) then
    call sd_ext_head_in_act()
  else
    call sd_ext_head_in_dbl()
  end if
end do

return

end subroutine sd_drt_ci_new

subroutine sd_ext_head_in_dbl()

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jb_sys, jml, jmr, jpad, jpadl, jpadlr, jud, just, linelp, logic_dh, &
                         lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, ngw2, ngw3, nlg1, nlg2, norb_dz, &
                         norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_dv, w0_sd, w0_sv, w0_td, w0_vv, w1_sd, w1_sv, &
                         w1_td, w1_tv
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: ijk, imap_1, intpos, isma, iwal, iwar, iwdl, iwdr, jk, jmlr, lmd, lmi, lmij, lmj, lpok, lra, lrd, lri, lrj, &
                     lrk, mpl, ni
real(kind=wp) :: w0, w0sd1, w0sd11, w0sd12, w0sd14, w0sd16, w0sd2, w0sd4, w0sd5, w0sd8, w0sd9, w0sv2, w0td1, w0td2, w0td3, w0td4, &
                 w0td5, w1sd11, w1sd12, w1sd5, w1sd8, w1sd9, w1sv2, w1td2, w1td3, w1td4, w1tv
integer(kind=iwp), external :: iwalk_ad

logic_dh = .true.
isma = Mul(iml,imr)
jmlr = Mul(jml,jmr)
lpok = jpadlr
select case (lpok)
  case (1)
    !===================================================================
    ! ss(1)     act -b^l-
    if (linelp /= 15) return
    lra = nlg1
    if (nlg2 == 1) then
      call ss_drl_act_bl(6,lra)
      if (jb_sys > 0) call ss_drl_act_bl_sgt0(6,lra)
    else
      call ss_arbl_act_bl(6,lra)
      if (jb_sys > 0) then
        call ss_arbl_act_bl_sgt0(6,lra)
        call ss_s_drl_act_bl_sgt0(6,lra)
      end if
    end if

  case (2)
    !===================================================================
    ! st(2)      act -b^l
    if ((linelp /= 15) .or. (nlg2 /= 2)) return
    lra = nlg1
    call st_arbl_act_bl(6,lra)
    if (jml == jmr) then
      call st_drl_act_bl(6,lra)
    end if
    if (jb_sys > 0) then
      call st_arbl_act_bl_sgt0(6,lra)
      if (jml == jmr) then
        call st_drl_act_bl_sgt0(6,lra)
      end if
    end if

  case (3)
    !===================================================================
    ! ts(3)   act -b^l ............................................
    if (linelp /= 15) return
    lra = nlg1
    if (nlg2 == 2) then
      call ts_arbl_act_bl(6,lra)
    end if

  case (4)
    !===================================================================
    ! stt(4) act-b^l ..............................................
    if (linelp /= 15) return
    lra = nlg1
    if (nlg2 == 1) return
    call stt_arbl_act_bl_sgt1(6,lra)

  case (5)
    !===================================================================
    ! tts(5) act-b^l ..............................................
    if ((linelp /= 15) .or. (nlg2 /= 2)) return
    lra = nlg1
    call tts_arbl_act_bl_sgt1(6,lra)
    call tts_drl_act_bl_sgt1(6,lra)

  case default ! (6)
    !===================================================================
    ! sd(6-1) a&r(02)-
    ! sd(6-2) c(22)a&(13)-
    ! sd(6-3) a&r(13)c'(22)-
    ! sd(6-4) a&r(23)c'(12)-
    ! sd(6-5) a&r(23)b&r(13)b^r(32)
    ! sd(6-6) a&r(13)b&r(23)b^r(32)
    ! sd(6-7) a&r(13)b&l(32)b^l(23)
    ! sd(6-8) a&r(23)b&l(32)b^l(13)
    ! sd(6-9) d&r&r(03)b^r(32)
    ! sd(6-10) d&r&l(12)b^l(23)
    ! sd(6-11) d&r&l(22)b^l(13)
    ! sd(6-12) d&r&l(33)b^l(02)
    ! sd(6-13) (22)d&r&l(33)b^l(13)
    ! sd(6-14) d&r&l(33)c"(22)b^l(13)
    ! sd(6-15) d&r&l(33)b^l(13)c'(22)
    ! sd(6-16) d&r&l(33)b^l(23)c'(12)
    if (linelp == 13) then
      if (jb_sys > 0) then
        call sd_adb_act_c_ext_ar(6)
      end if
      do lri=norb_frz+1,norb_dz
        lmi = lsm_inn(lri)
        if (lmi /= jmlr) cycle
        w0sd1 = w0_sd(1)
        w0sd9 = w0_sd(9)
        w1sd9 = w1_sd(9)
        w0sd12 = w0_sd(12)
        w1sd12 = w1_sd(12)
        ni = mod(norb_dz-lri,2)
        if (ni == 1) w0sd1 = -w0sd1
        if (ni == 1) w0sd9 = -w0sd9
        if (ni == 1) w1sd9 = -w1sd9
        if (ni == 1) w0sd12 = -w0sd12
        if (ni == 1) w1sd12 = -w1sd12
        if ((jml == 1) .and. (lmi == jmr)) then
          iwdl = just(lri,lri)
          iwdr = jud(lri)
          do mpl=1,mhlp
            iwal = lpnew_lwei(mpl)
            iwar = lpnew_rwei(mpl)
            lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
          end do

          ! sd(6-1) a&r(02)-
          do mpl=1,mtype
            vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd1
            vplp_w1(mpl) = vplpnew_w1(mpl)*w0sd1
          end do
          call ar_sd_ext_ar(26,lri,lrj,isma)
          call ar_sd_ext_rest(lri)

          ! sd(6-12) d&r&l(33)b^l(02)
          do mpl=1,mtype
            vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd12
            vplp_w1(mpl) = vplpnew_w1(mpl)*w1sd12
          end do
          do lrk=1,lri-1
            call drl_bl_ext_ar_new(6,lrk,lri)
          end do
          ! sd(6-9) d&r&r(03)b^r(32)
          do lrk=norb_frz+1,lri-1
            iwdl = just(lrk,lrk)
            iwdr = jud(lri)
            do mpl=1,mhlp
              iwal = lpnew_lwei(mpl)
              iwar = lpnew_rwei(mpl)
              lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
            end do
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd9
              vplp_w1(mpl) = vplpnew_w1(mpl)*w1sd9
            end do
            call drr_br_ext_ar(6,lrk,lri)
          end do
        end if
      end do

      do lri=norb_frz+1,norb_dz
        lmi = lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj = lsm_inn(lrj)
          lmij = Mul(lmi,lmj)
          if (lmij /= jml) cycle

          if (lmi == jmr) then
            iwdl = just(lri,lrj)
            iwdr = jud(lri)
            w0sd2 = w0_sd(2)
            w0sd11 = w0_sd(11)
            w1sd11 = w1_sd(11)
            w0sd14 = w0_sd(14)
            ni = mod(norb_dz-lrj,2)
            if (ni == 1) then
              w0sd2 = -w0sd2
              w0sd11 = -w0sd11
              w1sd11 = -w1sd11
              w0sd14 = -w0sd14
            end if
            ! sd(6-2) c(22)-a&r(13)-
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd2
            end do
            do mpl=1,mhlp
              iwal = lpnew_lwei(mpl)
              iwar = lpnew_rwei(mpl)
              lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
            end do
            call ar_sd_ext_ar(25,lrj,lri,isma)
            call ar_sd_ext_rest(lrj)
            ! sd(6-11) d&r&l(22)b^l(13)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd11
              vplp_w1(mpl) = vplpnew_w1(mpl)*w1sd11
            end do
            call drl_bl_ext_ar_new(6,lri,lrj)
            ! sd(6-14) (22)d&r&l(33)b^l(13)
            ! sd(6-14) d&r&l(33)c"(22)b^l(13)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd14
              vplp_w1(mpl) = Zero
            end do
            do lrk=1,lrj-1
              if (lrk == lri) cycle
              call drl_bl_ext_ar_new(6,lrk,lrj)
            end do
          end if
          ! sd(6-4) a&r(23)-c'(12)-
          w0sd4 = w0_sd(4)
          w0sd16 = w0_sd(16)
          ni = mod(norb_dz-lri,2)
          if (ni == 0) then
            w0sd4 = -w0sd4
            w0sd16 = -w0sd16
          end if
          if (lmj == jmr) then
            iwdl = just(lri,lrj)
            iwdr = jud(lrj)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd4
            end do
            do mpl=1,mhlp
              iwal = lpnew_lwei(mpl)
              iwar = lpnew_rwei(mpl)
              lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
            end do
            call ar_sd_ext_ar(28,lri,lrj,isma)
            call ar_sd_ext_rest(lri)
            ! sd(6-16) d&r&l(33)b^l(23)c'(12)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd16
              vplp_w1(mpl) = Zero
            end do
            do lrk=1,lri-1
              call drl_bl_ext_ar_new(6,lrk,lri)
            end do
          end if
          ! sd(6-5) a&r(23)b&r(13)b^r(32)
          do lrd=lrj+1,norb_dz
            lmd = lsm_inn(lrd)
            if (lmd /= jmr) cycle

            w0sd5 = w0_sd(5)
            w1sd5 = w1_sd(5)
            ni = mod(lrj-lri+norb_dz-lrd,2)
            if (ni == 0) then
              w0sd5 = -w0sd5
              w1sd5 = -w1sd5
            end if
            iwdl = just(lri,lrj)
            iwdr = jud(lrd)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd5
              vplp_w1(mpl) = vplpnew_w1(mpl)*w1sd5
            end do
            do mpl=1,mhlp
              iwal = lpnew_lwei(mpl)
              iwar = lpnew_rwei(mpl)
              lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
            end do
            ijk = lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lrd-norb_frz)
            intpos = intind_ijka(ijk)
            call ar_br_br_ext_ar_new(6,intpos,isma)
          end do
          ! sd(6-8) a&r(23)b&l(32)b^l(13)
          do lrd=lri+1,lrj-1
            lmd = lsm_inn(lrd)
            if (lmd /= jmr) cycle
            iwdl = just(lri,lrj)
            iwdr = jud(lrd)
            w0sd8 = w0_sd(8)
            w1sd8 = w1_sd(8)
            ni = mod(lrj-lri+norb_dz-lrd,2)
            if (ni == 0) then
              w0sd8 = -w0sd8
              w1sd8 = -w1sd8
            end if
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd8
              vplp_w1(mpl) = vplpnew_w1(mpl)*w1sd8
            end do
            do mpl=1,mhlp
              iwal = lpnew_lwei(mpl)
              iwar = lpnew_rwei(mpl)
              lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
            end do
            ijk = lri-norb_frz+ngw2(lrd-norb_frz)+ngw3(lrj-norb_frz)
            intpos = intind_ijka(ijk)
            call ar_bl_bl_ext_ar_new(6,intpos,isma,1)
          end do
        end do
      end do
    end if
    !===================================================================
    ! sd(6)    act -br-br-
    if (linelp == 20) then
      jk = nlg1
      call sd_ar_act_brbr(6,jk)
      if (jb_sys > 0) call sd_ar_act_brbr_sgt0(6,jk)
    end if
    ! sd(6)  act -bl-bl-
    if (linelp == 22) then
      jk = nlg1
      call sd_ar_act_blbl(6,jk)
      if (jb_sys > 0) call sd_ar_act_blbl_sgt0(6,jk)
    end if

  case (8)
    !===================================================================
    ! sdd(8) act-b^l -c'- -c"-..........................................
    lra = nlg1
    if (linelp == 13) then
      call sdd_drlbl_act_c_sgt0(6)
      call sdd_drrbr_act_c_sgt0(6)
      call sdd_abb_act_c_sgt0(6)
      call sdd_ar_act_c_sd_ext_sgt0(6)
    else
      lra = nlg1
      if (linelp == 20) call sdd_ar_act_brbr_sgt0(6,lra)
      if (linelp == 22) call sdd_ar_act_blbl_sgt0(6,lra)
    end if

  !case (9)
  !  !==================================================================
  !  ! dds(9) act -c'- ..............................................
  !  if (linelp == 13) then
  !    lra = nlg1
  !    call dds_abb_act_c_sgt0(6)
  !  end if

  case (10)
    !===================================================================
    ! sv(10)     act -b^r-
    ! sv(10-1) ar(13)-br(23)-
    ! sv(10-2) ar(23)-br(13)-
    ! sv(10-3) drr(03)-
    if ((linelp /= 16) .or. (nlg2 /= 2)) return
    lra = nlg1
    if (jb_sys > 0) then
      call sv_arbr_act_br_sgt0(6,lra)
    end if
    do lri=norb_frz+1,norb_dz
      lmi = lsm_inn(lri)
      do lrj=lri,norb_dz
        lmj = lsm_inn(lrj)
        lmij = Mul(lmi,lmj)
        if (lmij /= jml) cycle
        w0sv2 = w0_sv(2)
        w1sv2 = w1_sv(2)
        ni = mod(lrj-lri,2)
        if (ni == 0) then
          w0sv2 = -w0sv2
          w1sv2 = -w1sv2
        end if
        !---------------------------------------------------------------
        iwdl = just(lri,lrj)
        iwdr = 0
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        if (lri /= lrj) then
          ! sv(10-2) ar(23)-br(13)-
          do mpl=1,mtype
            vplp_w0(mpl) = vplpnew_w0(mpl)*w0sv2
            vplp_w1(mpl) = vplpnew_w1(mpl)*w1sv2
          end do
          ijk = lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz)
          intpos = intind_ijka(ijk)
          call ar_br_br_ext_ar_new(6,intpos,isma)
        else
          ! sv(10-3) drr(03)-
          do mpl=1,mtype
            vplp_w0(mpl) = vplpnew_w0(mpl)*w0_sv(3)
            vplp_w1(mpl) = vplpnew_w1(mpl)*w1_sv(3)
          end do
          call drr_br_ext_ar(6,lri,lra)
        end if
      end do
    end do

  case (11)
    !===================================================================
    ! tt(11)   act -b^l
    if (linelp /= 15) return
    lra = nlg1
    if (nlg2 == 1) then
      call tt_drl_act_bl(6,lra)
    else
      call tt_arbl_act_bl(6,lra)
    end if

  case (12)
    !===================================================================
    ! tttt(12)   act -b^l
    if (linelp /= 15) return
    lra = nlg1
    if (nlg2 == 1) then
      call tttt_drl_act_bl_sgt1(6,lra)
    else
      call tttt_arbl_act_bl_sgt1(6,lra)
    end if

  case (13)
    !===================================================================
    ! td(13)                                  act -c'-
    ! td(13-1) (22)a&(23)
    ! td(13-1) a&(23)c'(22)
    ! td(13-2) a&(23)b&r(23)b^r(32)
    ! td(13-3) a&(23)b&l(32)b^l(23)
    ! td(13-4) d&r&l(22)b^l(23)
    ! td(13-5) (22)d&&l(33)b^l(23)
    ! td(13-5) d&rl(33)c"(22)b^l(23)
    ! td(13-5) d&rl(33)b^l(23)c'(22)
    if (linelp == 13) then
      do lri=norb_frz+1,norb_dz
        lmi = lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj = lsm_inn(lrj)
          lmij = Mul(lmi,lmj)
          if (lmij /= jml) cycle
          w0td1 = w0_td(1)
          w0td4 = w0_td(4)
          w1td4 = w1_td(4)
          w0td5 = w0_td(5)
          ni = mod(norb_dz-lrj,2)
          if (ni == 1) then
            w0td1 = -w0td1
            w0td4 = -w0td4
            w1td4 = -w1td4
            w0td5 = -w0td5
          end if
          iwdl = just(lri,lrj)
          !-------------------------------------------------------------
          ! td(13-1) (22)a&(23)
          if (lmi == jmr) then
            iwdr = jud(lri)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0td1
            end do
            do mpl=1,mhlp
              iwal = lpnew_lwei(mpl)
              iwar = lpnew_rwei(mpl)
              lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
            end do
            call ar_sd_ext_ar(43,lrj,lri,isma)
            call ar_sd_ext_rest(lrj)
            ! td(13-4) d&r&l(22)b^l(23)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0td4
              vplp_w1(mpl) = vplpnew_w1(mpl)*w1td4
            end do
            call drl_bl_ext_ar_new(6,lri,lrj)
            ! td(13-5) (22)d&rl(33)b^l(23)
            ! td(13-5) d&rl(33)c"(22)b^l(23)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0td5
              vplp_w1(mpl) = Zero
            end do
            do lrk=1,lrj-1
              if (lrk == lri) cycle
              call drl_bl_ext_ar_new(6,lrk,lrj)
            end do
          end if
          !-------------------------------------------------------------
          ! td(13-1) a&(23)c'(22)
          if (lmj == jmr) then
            iwdr = jud(lrj)
            w0td1 = w0_td(1)
            w0td5 = w0_td(5)
            ni = mod(norb_dz-lri,2)
            if (ni == 0) then
              w0td1 = -w0td1
              w0td5 = -w0td5
            end if
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0td1
            end do
            do mpl=1,mhlp
              iwal = lpnew_lwei(mpl)
              iwar = lpnew_rwei(mpl)
              lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
            end do
            call ar_sd_ext_ar(46,lri,lrj,isma)
            call ar_sd_ext_rest(lri)
            ! td(13-5) d&rl(33)b^l(23)c'(22)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0td5
              vplp_w1(mpl) = Zero
            end do
            do lrk=1,lri-1
              call drl_bl_ext_ar_new(6,lrk,lri)
            end do
          end if
          !-------------------------------------------------------------
          ! td(13-2) a&(23)b&r(23)b^r(32)
          do lrd=lrj+1,norb_dz
            lmd = lsm_inn(lrd)
            if (lmd /= jmr) cycle
            iwdr = jud(lrd)
            ijk = lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lrd-norb_frz)
            intpos = intind_ijka(ijk)
            w0td2 = w0_td(2)
            w1td2 = w1_td(2)
            ni = mod(lrj-lri+norb_dz-lrd,2)
            if (ni == 0) then
              w0td2 = -w0td2
              w1td2 = -w1td2
            end if
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0td2
              vplp_w1(mpl) = vplpnew_w1(mpl)*w1td2
            end do
            do mpl=1,mhlp
              iwal = lpnew_lwei(mpl)
              iwar = lpnew_rwei(mpl)
              lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
            end do
            call ar_br_br_ext_ar_new(6,intpos,isma)
          end do
          !-------------------------------------------------------------
          ! td(13-3) a&(23)b&l(32)b^l(23)
          do lrd=lri+1,lrj-1
            lmd = lsm_inn(lrd)
            if (lmd /= jmr) cycle
            iwdr = jud(lrd)
            w0td3 = w0_td(3)
            w1td3 = w1_td(3)
            ni = mod(lrj-lri+norb_dz-lrd,2)
            if (ni == 0) then
              w0td3 = -w0td3
              w1td3 = -w1td3
            end if
            ijk = lri-norb_frz+ngw2(lrd-norb_frz)+ngw3(lrj-norb_frz)
            intpos = intind_ijka(ijk)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0td3
              vplp_w1(mpl) = vplpnew_w1(mpl)*w1td3
            end do
            do mpl=1,mhlp
              iwal = lpnew_lwei(mpl)
              iwar = lpnew_rwei(mpl)
              lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
            end do
            call ar_bl_bl_ext_ar_new(6,intpos,isma,1)
          end do
          !-------------------------------------------------------------
        end do
      end do
    end if
    !===================================================================
    ! td(13) act -br-br-
    if (linelp == 20) then
      jk = nlg1
      call td_ar_act_brbr(6,jk)
    end if
    ! td(13)  act -bl-bl-
    if (linelp == 22) then
      jk = nlg1
      call td_ar_act_blbl(6,jk)
    end if

  case (15)
    !===================================================================
    ! ttdd(15)    act -c'-   -br-br- -bl-bl
    lra = nlg1
    if (linelp == 13) then                  ! no complete ar-
      !if (nlg2 == 1) then
      call ttdd_drlbl_act_c_sgt1(6)
      call ttdd_abb_act_c_sgt1(6)           ! ????
      call ttdd_ar_act_c_ttdd_ext_sgt1(6)
      !end if
    end if
    if (linelp == 20) then
      call ttdd_ar_act_brbr_sgt0(6,lra)
    end if
    if (linelp == 22) then
      call ttdd_ar_act_blbl_sgt0(6,lra)
    end if

  case (17)
    !===================================================================
    ! tv(17)    act -b^r-
    if ((linelp /= 16) .or. (nlg2 /= 2)) return
    lra = nlg1
    do lri=norb_frz+1,norb_dz-1
      lmi = lsm_inn(lri)
      do lrj=lri+1,norb_dz
        lmj = lsm_inn(lrj)
        lmij = Mul(lmi,lmj)
        if (lmij /= jml) cycle

        w1tv = w1_tv
        ni = mod(lrj-lri,2)
        if (ni == 0) w1tv = -w1tv

        !---------------------------------------------------------------
        ! tv(17) ar(23)-br(23)-
        iwdl = just(lri,lrj)
        iwdr = 0
        do mpl=1,mtype
          vplp_w0(mpl) = Zero
          vplp_w1(mpl) = vplpnew_w1(mpl)*w1tv
        end do
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        ijk = lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz)
        intpos = intind_ijka(ijk)
        call ar_br_br_ext_ar_new(6,intpos,isma)
      end do
    end do

  case (18)
    !===================================================================
    ! ttv(19) act -b^r- ................................................
    if ((linelp /= 16) .or. (nlg2 /= 2)) return
    lra = nlg1
    call ttv_arbr_act_c_sgt1(6,lra)

  case (19)
    !===================================================================
    ! dd(19) act -b^l- .................................................
    if (linelp /= 15) return
    lra = nlg1
    if (nlg2 == 1) then
      call dd_drl_act_bl(6,lra)
    else
      call dd_arbl_act_bl(6,lra)
    end if

  case (20)
    !===================================================================
    ! dddd(19) act -b^l- -br- ..........................................
    lra = nlg1
    if (linelp == 15) then
      if (nlg2 == 2) then
        call dddd_arbl_act_bl_sgt0(6,lra)
      else
        call dddd_drl_act_bl_sgt0(6,lra)
      end if
    end if

  case (21)
    !===================================================================
    ! dd1(21) act -b^l-  ...............................................
    if ((linelp /= 15) .or. (nlg2 /= 2)) return
    lra = nlg1
    call dd1_arbl_act_bl_sgt0(6,lra)

  case (22)
    !===================================================================
    ! d1d(22) act -b^l-  ...............................................
    lra = nlg1
    if ((linelp == 15) .and. (nlg2 == 2)) then
      call d1d_arbl_act_bl_sgt0(6,lra)
      call d1d_drl_act_bl_sgt0(6,lra)
    end if

  case (23)
    !===================================================================
    ! dv(23) act -c'-..................................................
    ! dv(23-1) ar(23)-
    ! dv(23-2) drl(33)-bl(23)-
    if (linelp == 13) then
      imap_1 = 0
      do lrd=norb_frz+1,norb_dz
        lmd = lsm_inn(lrd)
        if (lmd /= jml) cycle
        ni = mod(norb_dz-lrd,2)

        ! dv(23-1) ar(23)-
        w0 = w0_dv(1)
        if (ni == 1) w0 = -w0_dv(1)
        iwdl = jud(lrd)
        iwdr = 0
        imap_1 = imap_1+1
        do mpl=1,mtype
          vplp_w0(mpl) = vplpnew_w0(mpl)*w0
        end do
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        call ar_sd_ext_ar(51,lrd,0,isma)
        call ar_sd_ext_rest(lrd)
        ! dv(23-2) drl(33)-bl(23)-
        w0 = w0_dv(2)
        if (ni == 1) w0 = -w0_dv(2)
        do mpl=1,mtype
          vplp_w0(mpl) = vplpnew_w0(mpl)*w0
          vplp_w1(mpl) = Zero
        end do
        do lrk=1,lrd-1
          call drl_bl_ext_ar_new(6,lrk,lrd)
        end do
      end do
    end if
    !===================================================================
    ! dv(23) act -br-br-................................................
    if (linelp == 20) then
      jk = nlg1
      do lrd=norb_frz+1,norb_dz
        lmd = lsm_inn(lrd)
        if (lmd /= jml) cycle
        ni = mod(norb_dz-lrd,2)
        w0 = w0_dv(1)
        if (ni == 1) w0 = -w0_dv(1)
        !...............................................................
        ! dv(23-1) ar(23)-
        iwdl = jud(lrd)
        iwdr = 0
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
        ijk = lrd-norb_frz+jk
        intpos = intind_ijka(ijk)
        call ar_br_br_ext_ar_new(6,intpos,isma)
      end do
    end if
    ! dv(23) act -bl-bl-................................................
    if (linelp /= 22) return
    jk = nlg1
    do lrd=norb_frz+1,norb_dz
      lmd = lsm_inn(lrd)
      if (lmd /= jml) cycle
      ni = mod(norb_dz-lrd,2)
      w0 = w0_dv(1)
      if (ni == 1) w0 = -w0_dv(1)
      ! dv(23-1) ar(23)-
      iwdl = jud(lrd)
      iwdr = 0
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
      ijk = lrd-norb_frz+jk
      intpos = intind_ijka(ijk)
      call ar_bl_bl_ext_ar_new(6,intpos,isma,1)
    end do

  case (24)
    !===================================================================
    ! d1v(24) act -c'-..................................................
    ! d1v(24-1) ar(13)-
    ! d1v(24-2) drl(33)-bl(13)-
    if (linelp == 13) then
      call d1v_ar_act_c_ext(6)
    end if
    if (linelp == 20) then
      jk = nlg1
      call d1v_ar_act_brbr_ext(6,jk)
    end if
    if (linelp == 22) then
      jk = nlg1
      call d1v_ar_act_blbl_ext(6,jk)
    end if

  case (25)
    !===================================================================
    ! vv(25) act -b^l- .................................................
    if (linelp /= 15) return
    lra = nlg1
    if (nlg2 == 1) then
      iwdl = 0
      iwdr = 0
      do mpl=1,mtype
        vplp_w0(mpl) = vplpnew_w0(mpl)*w0_vv
        vplp_w1(mpl) = Zero
      end do
      do mpl=1,mhlp
        iwal = lpnew_lwei(mpl)
        iwar = lpnew_rwei(mpl)
        lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
        lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
      end do
      call drl_bl_sum_ar_new(6,0,0,lra)
      !do lrk=1,norb_dz
      !  call drl_bl_ext_ar_new(6,lrk,lra)
      !end do
    end if

  case (7,9,14,16,26)
end select

return

end subroutine sd_ext_head_in_dbl

subroutine sd_ext_head_in_act()

use gugaci_global, only: iml, imr, linelp, logic_dh, nlg1, nlg2
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: intpos, isma, lri, lrj

logic_dh = .false.
lri = nlg1
lrj = nlg2
intpos = nlg1
isma = Mul(iml,imr)

select case (linelp)
  case default ! (1)
    ! line=1 a&r<-->a^r
    call ar_sd_ext_ar(100,lri,lrj,isma)
    call ar_sd_ext_rest(lri)
  case (4)
    ! line=4 a&r--b&r--b^r<-->a^r
    call ar_br_br_ext_ar_new(6,intpos,isma)
  case (7)
    ! line=7 a&r--b&l--b^l<-->a^r
    call ar_bl_bl_ext_ar_new(6,intpos,isma,1)
  case (10)
    ! line=10 d&rr--b^r<-->a^r
    call drr_br_ext_ar(6,lri,lrj)
  case (12)
    ! line=12 d&rl--b^l<-->a^r
    call drl_bl_ext_ar_new(6,lri,lrj)
  case (2:3,5:6,8:9,11)
end select

return

end subroutine sd_ext_head_in_act

subroutine sd_ext_space_w01plp_value()

use gugaci_global, only: w0_sdplp, w0g25a, w0g26a, w0g27, w0g28a, w0g29, w0g30, w0g31, w0g32, w0plp25, w0plp26, w0plp27, w0plp28, &
                         w0plp29, w0plp30, w0plp31, w0plp32, w1g27, w1g31, w1g32, w1plp27, w1plp31, w1plp32

implicit none

w0plp25 = w0_sdplp*w0g25a
w0plp26 = w0_sdplp*w0g26a
w0plp27 = w0_sdplp*w0g27
w1plp27 = w0_sdplp*w1g27
w0plp28 = w0_sdplp*w0g28a
w0plp29 = w0_sdplp*w0g29
w0plp30 = w0_sdplp*w0g30
w0plp31 = w0_sdplp*w0g31
w1plp31 = w0_sdplp*w1g31
w0plp32 = w0_sdplp*w0g32
w1plp32 = w0_sdplp*w1g32

end subroutine sd_ext_space_w01plp_value

subroutine gsd_samesym_aaa(lri,isma)

use gugaci_global, only: ibsm_ext, icnt_base, iesm_ext, intind_iabc, intind_iaqq, iwt_orb_ext, m_jd, nabc, ngw2, ngw3, norb_ext, &
                         norb_number, value_lpext, vint_ci, w0g28a, w0plp27, w0plp28, w0plp31, w0plp32, w1plp27, w1plp31, w1plp32
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lri, isma
integer(kind=iwp) :: ia, ia0, iabc, iabc0, iaqq, iasta, ib, ibend, ibsta, ic, ilwei, intpos, iposint, jcoffset, lrc

ia0 = (lri-1)*norb_ext
iabc0 = (lri-1)*nabc

ic = m_jd
lrc = norb_number(ic)
jcoffset = lrc*2-2

iasta = ibsm_ext(isma)
ibend = iesm_ext(isma)
ibsta = iasta+1

ilwei = icnt_base+iwt_orb_ext(iasta,ibsta)
do ib=ibsta,ic-1
  do ia=iasta,ib-1
    ! g31   type_12 arbrb^ra^r
    iabc = iabc0+ia+ngw2(ib)+ngw3(ic)
    iposint = intind_iabc(iabc)
    value_lpext(ilwei) = vint_ci(iposint+1)*w0plp31+vint_ci(iposint+2)*w1plp31
    ilwei = ilwei+1
  end do
end do

ib = ic
ilwei = icnt_base+iwt_orb_ext(iasta,ib)
do ia=iasta,ib-1
  ! g28     cw-ar                    330
  iaqq = ia0+ia
  intpos = intind_iaqq(iaqq)
  iposint = intpos+jcoffset
  value_lpext(ilwei) = (vint_ci(iposint)/w0g28a+vint_ci(iposint+1))*w0plp28
  ilwei = ilwei+1
end do

ia = ic
do ib=ic+1,ibend
  ! g25,g27: bl(20)-drl(11)         220
  iaqq = ia0+ib
  intpos = intind_iaqq(iaqq)
  iposint = intpos+jcoffset
  ilwei = icnt_base+iwt_orb_ext(ic,ib)
  value_lpext(ilwei) = vint_ci(iposint)*w0plp27-vint_ci(iposint+1)*w1plp27
end do

do ib=ic+1,ibend
  ilwei = icnt_base+iwt_orb_ext(iasta,ib)
  do ia=iasta,ic-1
    ! g32a   type g12
    iabc = iabc0+ia+ngw2(ic)+ngw3(ib)
    iposint = intind_iabc(iabc)
    value_lpext(ilwei) = vint_ci(iposint+2)*w0plp32-vint_ci(iposint)*w1plp32 !severe_new_er
    ilwei = ilwei+1
  end do
end do

do ib=ic+2,ibend
  ilwei = icnt_base+iwt_orb_ext(ic+1,ib)
  do ia=ic+1,ib-1
    iabc = iabc0+ic+ngw2(ia)+ngw3(ib)
    iposint = intind_iabc(iabc)
    ! g32b  type 11
    value_lpext(ilwei) = vint_ci(iposint+1)*w0plp32-vint_ci(iposint)*w1plp32 !severe_new_
    ilwei = ilwei+1
  end do
end do

end subroutine gsd_samesym_aaa

subroutine gsd_diffsamesym_abb(lri,isma,ismb)

use gugaci_global, only: ibsm_ext, icnt_base, iesm_ext, intind_iabc, intind_iaqq, iwt_orb_ext, m_jc, m_jd, nabc, ngw2, ngw3, &
                         norb_ext, norb_number, value_lpext, vint_ci, w0g28a, w0plp28, w0plp31, w0plp32, w1plp31, w1plp32
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lri, isma, ismb
integer(kind=iwp) :: ia, ia0, iabc, iabc0, iaend, iaqq, iasta, ib, ibend, ibsta, ic, ilwei, intpos, iposint, jb, jcoffset, lrc

ia0 = (lri-1)*norb_ext
iabc0 = (lri-1)*nabc

ic = m_jd
lrc = norb_number(ic)
jcoffset = lrc*2-2

iasta = ibsm_ext(isma)
iaend = iesm_ext(isma)
ibsta = ibsm_ext(ismb)
ibend = iesm_ext(ismb)

ilwei = icnt_base+iwt_orb_ext(iasta,ibsta)
do ib=ibsta,ic-1
  do ia=iasta,iaend
    ! g31
    ! g31   type_12 arbrb^ra^r
    iabc = iabc0+ia+ngw2(ib)+ngw3(ic)
    iposint = intind_iabc(iabc)
    value_lpext(ilwei) = vint_ci(iposint+1)*w0plp31+vint_ci(iposint+2)*w1plp31 !severe_new_
    ilwei = ilwei+1
  end do
end do

ilwei = icnt_base+iwt_orb_ext(iasta,ic+1)
jb = m_jc
do ib=ic+1,ibend
  jb = jb+1
  do ia=iasta,iaend
    ! g32a
    iabc = iabc0+ia+ngw2(ic)+ngw3(ib)
    iposint = intind_iabc(iabc)
    value_lpext(ilwei) = vint_ci(iposint+2)*w0plp32-vint_ci(iposint)*w1plp32 !severe_new_err
    ilwei = ilwei+1
  end do
end do

ib = ic
ilwei = icnt_base+iwt_orb_ext(iasta,ib)
do ia=iasta,iaend      !ib-1      !severe_error_1020
  ! g28
  iaqq = ia0+ia
  intpos = intind_iaqq(iaqq)
  iposint = intpos+jcoffset
  !value_lpext(ilwei) = (vint_ci(iposint2)+
  value_lpext(ilwei) = w0plp28*(vint_ci(iposint)/w0g28a+vint_ci(iposint+1))
  ilwei = ilwei+1
end do

end subroutine gsd_diffsamesym_abb

subroutine gsd_diffsamesym_aab(lri,isma,ismb)

use gugaci_global, only: ibsm_ext, icnt_base, iesm_ext, intind_iabc, intind_iaqq, iwt_orb_ext, m_jd, nabc, ngw2, ngw3, norb_ext, &
                         norb_number, value_lpext, vint_ci, w0plp27, w0plp32, w1plp27, w1plp32
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lri, isma, ismb
integer(kind=iwp) :: ia, ia0, iabc, iabc0, iaend, iaqq, iasta, ib, ibend, ibsta, ic, ilwei, intpos, iposint, jcoffset, lrc

ia0 = (lri-1)*norb_ext
iabc0 = (lri-1)*nabc

ic = m_jd
lrc = norb_number(ic)
jcoffset = lrc*2-2

ibsta = ibsm_ext(ismb)
ibend = iesm_ext(ismb)
iasta = ibsm_ext(isma)
iaend = iesm_ext(isma)

! if (ic-1 > iasta) then
do ib=ibsta,ibend
  ilwei = icnt_base+iwt_orb_ext(iasta,ib)
  do ia=iasta,ic-1
    ! g32a
    iabc = iabc0+ia+ngw2(ic)+ngw3(ib)
    iposint = intind_iabc(iabc)
    value_lpext(ilwei) = vint_ci(iposint+2)*w0plp32-vint_ci(iposint)*w1plp32 !severe_new_
    ilwei = ilwei+1
  end do
end do

do ib=ibsta,ibend
  ilwei = icnt_base+iwt_orb_ext(ic+1,ib)
  do ia=ic+1,iaend      !ib-1      !severe_error_1020
    iabc = iabc0+ic+ngw2(ia)+ngw3(ib)
    iposint = intind_iabc(iabc)
    ! g32b
    value_lpext(ilwei) = vint_ci(iposint+1)*w0plp32-vint_ci(iposint)*w1plp32 !severe_new_err
    ilwei = ilwei+1
  end do
end do

ia = ic
do ib=ibsta,ibend
  ! g25
  iaqq = ia0+ib
  intpos = intind_iaqq(iaqq)
  iposint = intpos+jcoffset
  ilwei = icnt_base+iwt_orb_ext(ic,ib)
  value_lpext(ilwei) = vint_ci(iposint)*w0plp27-vint_ci(iposint+1)*w1plp27
end do

end subroutine gsd_diffsamesym_aab

subroutine gsd_arlp_s1(lri)

use gugaci_global, only: icnt_base, intind_iaqq, isegdownwei, m_jd, norb_ext, norb_number, value_lpext, vint_ci, w0plp26, w0plp29, &
                         w0plp30
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lri
integer(kind=iwp) :: ia0, iaqq, ic, ilwei, intoffset, intpos, iposint, is1orb, lrk

ia0 = (lri-1)*norb_ext
ic = m_jd

iaqq = ia0+ic
intpos = intind_iaqq(iaqq)

ilwei = icnt_base+isegdownwei-norb_ext+1
do is1orb=1,ic-1
  ! g30 -b^rd^rr
  lrk = norb_number(is1orb)
  intoffset = (lrk-1)*2
  iposint = intpos+intoffset
  value_lpext(ilwei) = vint_ci(iposint)*w0plp30
  ilwei = ilwei+1
end do

! g26 -a^r     610
lrk = norb_number(ic)
intoffset = (lrk-1)*2
iposint = intpos+intoffset
value_lpext(ilwei) = vint_ci(iposint)*w0plp26
ilwei = ilwei+1

! g29 -dl^ra^l
do is1orb=ic+1,norb_ext
  lrk = norb_number(is1orb)
  intoffset = (lrk-1)*2
  iposint = intpos+intoffset
  value_lpext(ilwei) = vint_ci(iposint)*w0plp29
  ilwei = ilwei+1
end do

end subroutine gsd_arlp_s1

subroutine td_drt_ci_new()

use gugaci_global, only: idisk_array, idisk_lp, iml, imr, ipae, ipael, jml, jmr, jpad, jpadl, jpadlr, linelp, lpblock_td, w0_sdplp
use Constants, only: One
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: jptyl, jptyr, lpb

!write(u6,*) '  td_wyb'

!logic_sd = .true.
call external_space_plpmode_value_td()
w0_sdplp = One
call sd_ext_space_w01plp_value()

idisk_lp = idisk_array(7)
do lpb=1,lpblock_td
  call read_lp()
  ipael = iml+9
  ipae = imr+1
  call get_jpty(jpadlr,jptyl,jptyr)
  call get_jp(jptyl,jml,jpadl,1)
  call get_jp(jptyr,jmr,jpad,1)
  !jmlr = Mul(jml,jmr)
  call gsd_determine_extarmode_paras(iml,imr,.false.)
  if (linelp <= 12) then
    call td_ext_head_in_act()
  else
    call td_ext_head_in_dbl()
  end if
end do

return

end subroutine td_drt_ci_new

subroutine td_ext_head_in_dbl()

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jb_sys, jml, jmr, jpad, jpadl, jpadlr, jud, just, linelp, logic_dh, &
                         lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, ngw2, ngw3, nlg1, nlg2, norb_dz, &
                         norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_dv, w0_sd, w0_sv, w0_td, w0_vv, w1_sd, w1_sv, &
                         w1_td, w1_tv
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: ijk, imap_1, intpos, isma, iwal, iwar, iwdl, iwdr, jk, jmlr, lmd, lmi, lmij, lmj, lpok, lra, lrd, lri, lrj, &
                     lrk, mpl, ni
real(kind=wp) :: w0, w0sd1, w0sd11, w0sd12, w0sd14, w0sd16, w0sd2, w0sd4, w0sd5, w0sd8, w0sd9, w0sv2, w0td1, w0td2, w0td3, w0td4, &
                 w0td5, w1sd11, w1sd12, w1sd5, w1sd8, w1sd9, w1sv2, w1td2, w1td3, w1td4, w1tv
integer(kind=iwp), external :: iwalk_ad

logic_dh = .true.
isma = Mul(iml,imr)
jmlr = Mul(jml,jmr)
lpok = jpadlr
select case (lpok)
  case (1)
    !===================================================================
    ! ss(1)     act -b^l-
    if (linelp /= 15) return
    lra = nlg1
    if (nlg2 == 1) then
      call ss_drl_act_bl(13,lra)
      if (jb_sys > 0) call ss_drl_act_bl_sgt0(13,lra)
    else
      call ss_arbl_act_bl(13,lra)
      if (jb_sys > 0) then
        call ss_arbl_act_bl_sgt0(13,lra)
        call ss_s_drl_act_bl_sgt0(13,lra)
      end if
    end if

  case (2)
    !===================================================================
    ! st(2)      act -b^l
    if ((linelp /= 15) .or. (nlg2 /= 2)) return
    lra = nlg1
    call st_arbl_act_bl(13,lra)
    if (jb_sys > 0) call st_arbl_act_bl_sgt0(13,lra)
    if (jml == jmr) then
      call st_drl_act_bl(13,lra)
      if (jb_sys > 0) call st_drl_act_bl_sgt0(13,lra)
    end if

  case (3)
    !===================================================================
    ! ts(3)   act -b^l ............................................
    if (linelp /= 15) return
    lra = nlg1
    if (nlg2 == 2) then
      call ts_arbl_act_bl(13,lra)
      if (jb_sys > 0) call ts_arbl_act_bl_sgt0(13,lra)
    end if

  case (4)
    !===================================================================
    ! stt(4) act-b^l ..............................................
    if (linelp /= 15) return
    lra = nlg1
    if (nlg2 == 1) return
    call stt_arbl_act_bl_sgt1(13,lra)

  case (5)
    !===================================================================
    ! tts(5) act-b^l ..............................................
    if ((linelp /= 15) .or. (nlg2 /= 2)) return
    lra = nlg1
    call tts_arbl_act_bl_sgt1(13,lra)
    call tts_drl_act_bl_sgt1(13,lra)

  case default ! (6)
    !===================================================================
    ! sd(6-1) a&r(02)-
    ! sd(6-2) c(22)a&(13)-
    ! sd(6-3) a&r(13)c'(22)-
    ! sd(6-4) a&r(23)c'(12)-
    ! sd(6-5) a&r(23)b&r(13)b^r(32)
    ! sd(6-6) a&r(13)b&r(23)b^r(32)
    ! sd(6-7) a&r(13)b&l(32)b^l(23)
    ! sd(6-8) a&r(23)b&l(32)b^l(13)
    ! sd(6-9) d&r&r(03)b^r(32)
    ! sd(6-10) d&r&l(12)b^l(23)
    ! sd(6-11) d&r&l(22)b^l(13)
    ! sd(6-12) d&r&l(33)b^l(02)
    ! sd(6-13) (22)d&r&l(33)b^l(13)
    ! sd(6-14) d&r&l(33)c"(22)b^l(13)
    ! sd(6-15) d&r&l(33)b^l(13)c'(22)
    ! sd(6-16) d&r&l(33)b^l(23)c'(12)
    if (linelp == 13) then
      if (jb_sys > 0) then
        call sd_adb_act_c_ext_ar(13)
      end if
      do lri=norb_frz+1,norb_dz
        lmi = lsm_inn(lri)
        if (lmi /= jmlr) cycle
        w0sd1 = w0_sd(1)
        w0sd9 = w0_sd(9)
        w1sd9 = w1_sd(9)
        w0sd12 = w0_sd(12)
        w1sd12 = w1_sd(12)
        ni = mod(norb_dz-lri,2)
        if (ni == 1) w0sd1 = -w0sd1
        if (ni == 1) w0sd9 = -w0sd9
        if (ni == 1) w1sd9 = -w1sd9
        if (ni == 1) w0sd12 = -w0sd12
        if (ni == 1) w1sd12 = -w1sd12
        if ((jml == 1) .and. (lmi == jmr)) then
          iwdl = just(lri,lri)
          iwdr = jud(lri)
          do mpl=1,mhlp
            iwal = lpnew_lwei(mpl)
            iwar = lpnew_rwei(mpl)
            lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
          end do

          ! sd(6-1) a&r(02)-
          do mpl=1,mtype
            vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd1
            vplp_w1(mpl) = vplpnew_w1(mpl)*w0sd1
          end do
          call ar_td_ext_ar(26,lri,lrj,isma)
          call ar_td_ext_rest(lri)

          ! sd(6-12) d&r&l(33)b^l(02)
          do mpl=1,mtype
            vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd12
            vplp_w1(mpl) = vplpnew_w1(mpl)*w1sd12
          end do
          do lrk=1,lri-1
            call drl_bl_ext_ar_new(13,lrk,lri)
          end do

          ! sd(6-9) d&r&r(03)b^r(32)
          do lrk=norb_frz+1,lri-1
            iwdl = just(lrk,lrk)
            iwdr = jud(lri)
            do mpl=1,mhlp
              iwal = lpnew_lwei(mpl)
              iwar = lpnew_rwei(mpl)
              lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
            end do
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd9
              vplp_w1(mpl) = vplpnew_w1(mpl)*w1sd9
            end do
            call drr_br_ext_ar(13,lrk,lri)
          end do
        end if
      end do

      do lri=norb_frz+1,norb_dz
        lmi = lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj = lsm_inn(lrj)
          lmij = Mul(lmi,lmj)
          if (lmij /= jml) cycle

          if (lmi == jmr) then
            iwdl = just(lri,lrj)
            iwdr = jud(lri)
            w0sd2 = w0_sd(2)
            w0sd11 = w0_sd(11)
            w1sd11 = w1_sd(11)
            w0sd14 = w0_sd(14)
            ni = mod(norb_dz-lrj,2)
            if (ni == 1) then
              w0sd2 = -w0sd2
              w0sd11 = -w0sd11
              w1sd11 = -w1sd11
              w0sd14 = -w0sd14
            end if
            ! sd(6-2) c(22)-a&r(13)-
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd2
            end do
            do mpl=1,mhlp
              iwal = lpnew_lwei(mpl)
              iwar = lpnew_rwei(mpl)
              lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
            end do
            call ar_td_ext_ar(25,lrj,lri,isma)
            call ar_td_ext_rest(lrj)
            ! sd(6-11) d&r&l(22)b^l(13)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd11
              vplp_w1(mpl) = vplpnew_w1(mpl)*w1sd11
            end do
            call drl_bl_ext_ar_new(13,lri,lrj)
            ! sd(6-14) (22)d&r&l(33)b^l(13)
            ! sd(6-14) d&r&l(33)c"(22)b^l(13)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd14
              vplp_w1(mpl) = Zero
            end do
            do lrk=1,lrj-1
              if (lrk == lri) cycle
              call drl_bl_ext_ar_new(13,lrk,lrj)
            end do
          end if
          ! sd(6-4) a&r(23)-c'(12)-
          w0sd4 = w0_sd(4)
          w0sd16 = w0_sd(16)
          ni = mod(norb_dz-lri,2)
          if (ni == 0) then
            w0sd4 = -w0sd4
            w0sd16 = -w0sd16
          end if
          if (lmj == jmr) then
            iwdl = just(lri,lrj)
            iwdr = jud(lrj)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd4
              !vplp_w1(mpl) = vplpnew_w1(mpl)*w0sd4
            end do
            do mpl=1,mhlp
              iwal = lpnew_lwei(mpl)
              iwar = lpnew_rwei(mpl)
              lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
            end do
            call ar_td_ext_ar(28,lri,lrj,isma)
            call ar_td_ext_rest(lri)
            ! sd(6-16) d&r&l(33)b^l(23)c'(12)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd16
              vplp_w1(mpl) = Zero
            end do
            do lrk=1,lri-1
              call drl_bl_ext_ar_new(13,lrk,lri)
            end do
          end if
          ! sd(6-5) a&r(23)b&r(13)b^r(32)
          do lrd=lrj+1,norb_dz
            lmd = lsm_inn(lrd)
            if (lmd /= jmr) cycle

            w0sd5 = w0_sd(5)
            w1sd5 = w1_sd(5)
            ni = mod(lrj-lri+norb_dz-lrd,2)
            if (ni == 0) then
              w0sd5 = -w0sd5
              w1sd5 = -w1sd5
            end if
            iwdl = just(lri,lrj)
            iwdr = jud(lrd)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd5
              vplp_w1(mpl) = vplpnew_w1(mpl)*w1sd5
            end do
            do mpl=1,mhlp
              iwal = lpnew_lwei(mpl)
              iwar = lpnew_rwei(mpl)
              lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
            end do
            ijk = lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lrd-norb_frz)
            intpos = intind_ijka(ijk)
            call ar_br_br_ext_ar_new(13,intpos,isma)
          end do
          ! sd(6-8) a&r(23)b&l(32)b^l(13)
          do lrd=lri+1,lrj-1
            lmd = lsm_inn(lrd)
            if (lmd /= jmr) cycle
            iwdl = just(lri,lrj)
            iwdr = jud(lrd)
            w0sd8 = w0_sd(8)
            w1sd8 = w1_sd(8)
            ni = mod(lrj-lri+norb_dz-lrd,2)
            if (ni == 0) then
              w0sd8 = -w0sd8
              w1sd8 = -w1sd8
            end if
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0sd8
              vplp_w1(mpl) = vplpnew_w1(mpl)*w1sd8
            end do
            do mpl=1,mhlp
              iwal = lpnew_lwei(mpl)
              iwar = lpnew_rwei(mpl)
              lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
            end do
            ijk = lri-norb_frz+ngw2(lrd-norb_frz)+ngw3(lrj-norb_frz)
            intpos = intind_ijka(ijk)
            call ar_bl_bl_ext_ar_new(13,intpos,isma,1)
          end do
        end do
      end do
    end if
    !===================================================================
    ! sd(6)    act -br-br-
    if (linelp == 20) then
      jk = nlg1
      call sd_ar_act_brbr(13,jk)
      if (jb_sys > 0) call sd_ar_act_brbr_sgt0(13,jk)
    end if
    ! sd(6)  act -bl-bl-
    if (linelp == 22) then
      jk = nlg1
      call sd_ar_act_blbl(13,jk)
      if (jb_sys > 0) call sd_ar_act_blbl_sgt0(13,jk)
    end if

  case (8)
    !===================================================================
    ! sdd(8) act-b^l -c'- -c"-..........................................
    if (linelp == 13) then
      lra = nlg1
      call sdd_drlbl_act_c_sgt0(13)
      call sdd_drrbr_act_c_sgt0(13)
      call sdd_abb_act_c_sgt0(13)
      call sdd_ar_act_c_sd_ext_sgt0(13)
    else
      lra = nlg1
      if (linelp == 20) call sdd_ar_act_brbr_sgt0(13,lra)
      if (linelp == 22) call sdd_ar_act_blbl_sgt0(13,lra)
    end if

  case (9)
    !===================================================================
    ! dds(9) act -c'- ..............................................
    if (linelp == 13) then
      lra = nlg1
      call dds_abb_act_c_sgt0(13)
    end if

  case (10)
    !===================================================================
    ! sv(10)     act -b^r-
    ! sv(10-1) ar(13)-br(23)-
    ! sv(10-2) ar(23)-br(13)-
    ! sv(10-3) drr(03)-
    if ((linelp /= 16) .or. (nlg2 /= 2)) return
    lra = nlg1
    if (jb_sys > 0) then
      call sv_arbr_act_br_sgt0(13,lra)
    end if
    do lri=norb_frz+1,norb_dz
      lmi = lsm_inn(lri)
      do lrj=lri,norb_dz
        lmj = lsm_inn(lrj)
        lmij = Mul(lmi,lmj)
        if (lmij /= jml) cycle
        w0sv2 = w0_sv(2)
        w1sv2 = w1_sv(2)
        ni = mod(lrj-lri,2)
        if (ni == 0) then
          w0sv2 = -w0sv2
          w1sv2 = -w1sv2
        end if
        !---------------------------------------------------------------
        iwdl = just(lri,lrj)
        iwdr = 0
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        if (lri /= lrj) then
          ! sv(10-2) ar(23)-br(13)-
          do mpl=1,mtype
            vplp_w0(mpl) = vplpnew_w0(mpl)*w0sv2
            vplp_w1(mpl) = vplpnew_w1(mpl)*w1sv2
          end do
          ijk = lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz)
          intpos = intind_ijka(ijk)
          call ar_br_br_ext_ar_new(13,intpos,isma)
        else
          ! sv(10-3) drr(03)-
          do mpl=1,mtype
            vplp_w0(mpl) = vplpnew_w0(mpl)*w0_sv(3)
            vplp_w1(mpl) = vplpnew_w1(mpl)*w1_sv(3)
          end do
          call drr_br_ext_ar(13,lri,lra)
        end if
      end do
    end do

  case (11)
    !===================================================================
    ! tt(11)   act -b^l
    if (linelp /= 15) return
    lra = nlg1
    if (nlg2 == 1) then
      call tt_drl_act_bl(13,lra)
    else
      call tt_arbl_act_bl(13,lra)
    end if

  case (12)
    !===================================================================
    ! tttt(12)   act -b^l
    if (linelp /= 15) return
    lra = nlg1
    if (nlg2 == 1) then
      call tttt_drl_act_bl_sgt1(13,lra)
    else
      call tttt_arbl_act_bl_sgt1(13,lra)
    end if

  case (13)
    !===================================================================
    ! td(13)                                  act -c'-
    ! td(13-1) (22)a&(23)
    ! td(13-1) a&(23)c'(22)
    ! td(13-2) a&(23)b&r(23)b^r(32)
    ! td(13-3) a&(23)b&l(32)b^l(23)
    ! td(13-4) d&r&l(22)b^l(23)
    ! td(13-5) (22)d&&l(33)b^l(23)
    ! td(13-5) d&rl(33)c"(22)b^l(23)
    ! td(13-5) d&rl(33)b^l(23)c'(22)
    if (linelp == 13) then
      do lri=norb_frz+1,norb_dz
        lmi = lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj = lsm_inn(lrj)
          lmij = Mul(lmi,lmj)
          if (lmij /= jml) cycle
          w0td1 = w0_td(1)
          w0td4 = w0_td(4)
          w1td4 = w1_td(4)
          w0td5 = w0_td(5)
          ni = mod(norb_dz-lrj,2)
          if (ni == 1) then
            w0td1 = -w0td1
            w0td4 = -w0td4
            w1td4 = -w1td4
            w0td5 = -w0td5
          end if
          iwdl = just(lri,lrj)
          !-------------------------------------------------------------
          ! td(13-1) (22)a&(23)
          if (lmi == jmr) then
            iwdr = jud(lri)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0td1
            end do
            do mpl=1,mhlp
              iwal = lpnew_lwei(mpl)
              iwar = lpnew_rwei(mpl)
              lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
            end do
            call ar_td_ext_ar(43,lrj,lri,isma)
            call ar_td_ext_rest(lrj)
            ! td(13-4) d&r&l(22)b^l(23)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0td4
              vplp_w1(mpl) = vplpnew_w1(mpl)*w1td4
            end do
            call drl_bl_ext_ar_new(13,lri,lrj)
            ! td(13-5) (22)d&rl(33)b^l(23)
            ! td(13-5) d&rl(33)c"(22)b^l(23)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0td5
              vplp_w1(mpl) = Zero
            end do
            do lrk=1,lrj-1
              if (lrk == lri) cycle
              call drl_bl_ext_ar_new(13,lrk,lrj)
            end do
          end if
          !-------------------------------------------------------------
          ! td(13-1) a&(23)c'(22)
          if (lmj == jmr) then
            iwdr = jud(lrj)
            w0td1 = w0_td(1)
            w0td5 = w0_td(5)
            ni = mod(norb_dz-lri,2)
            if (ni == 0) then
              w0td1 = -w0td1
              w0td5 = -w0td5
            end if
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0td1
            end do
            do mpl=1,mhlp
              iwal = lpnew_lwei(mpl)
              iwar = lpnew_rwei(mpl)
              lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
            end do
            call ar_td_ext_ar(46,lri,lrj,isma)
            call ar_td_ext_rest(lri)
            ! td(13-5) d&rl(33)b^l(23)c'(22)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0td5
              vplp_w1(mpl) = Zero
            end do
            do lrk=1,lri-1
              call drl_bl_ext_ar_new(13,lrk,lri)
            end do
          end if
          !-------------------------------------------------------------
          ! td(13-2) a&(23)b&r(23)b^r(32)
          do lrd=lrj+1,norb_dz
            lmd = lsm_inn(lrd)
            if (lmd /= jmr) cycle
            iwdr = jud(lrd)
            ijk = lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lrd-norb_frz)
            intpos = intind_ijka(ijk)
            w0td2 = w0_td(2)
            w1td2 = w1_td(2)
            ni = mod(lrj-lri+norb_dz-lrd,2)
            if (ni == 0) then
              w0td2 = -w0td2
              w1td2 = -w1td2
            end if
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0td2
              vplp_w1(mpl) = vplpnew_w1(mpl)*w1td2
            end do
            do mpl=1,mhlp
              iwal = lpnew_lwei(mpl)
              iwar = lpnew_rwei(mpl)
              lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
            end do
            call ar_br_br_ext_ar_new(13,intpos,isma)
          end do
          !-------------------------------------------------------------
          ! td(13-3) a&(23)b&l(32)b^l(23)
          do lrd=lri+1,lrj-1
            lmd = lsm_inn(lrd)
            if (lmd /= jmr) cycle
            iwdr = jud(lrd)
            w0td3 = w0_td(3)
            w1td3 = w1_td(3)
            ni = mod(lrj-lri+norb_dz-lrd,2)
            if (ni == 0) then
              w0td3 = -w0td3
              w1td3 = -w1td3
            end if
            ijk = lri-norb_frz+ngw2(lrd-norb_frz)+ngw3(lrj-norb_frz)
            intpos = intind_ijka(ijk)
            do mpl=1,mtype
              vplp_w0(mpl) = vplpnew_w0(mpl)*w0td3
              vplp_w1(mpl) = vplpnew_w1(mpl)*w1td3
            end do
            do mpl=1,mhlp
              iwal = lpnew_lwei(mpl)
              iwar = lpnew_rwei(mpl)
              lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
            end do
            call ar_bl_bl_ext_ar_new(13,intpos,isma,1)
          end do
          !-------------------------------------------------------------
        end do
      end do
    end if
    !===================================================================
    ! td(13) act -br-br-
    if (linelp == 20) then
      jk = nlg1
      call td_ar_act_brbr(13,jk)
    end if
    ! td(13)  act -bl-bl-
    if (linelp == 22) then
      jk = nlg1
      call td_ar_act_blbl(13,jk)
    end if

  case (15)
    !===================================================================
    ! ttdd(15)    act -c'-   -br-br- -bl-bl
    lra = nlg1
    if (linelp == 13) then                  ! no complete ar-
      !if (nlg2 == 1) then
      call ttdd_drlbl_act_c_sgt1(13)
      call ttdd_abb_act_c_sgt1(13)             ! ????
      call ttdd_ar_act_c_ttdd_ext_sgt1(13)
      !end if
    end if
    if (linelp == 20) then
      call ttdd_ar_act_brbr_sgt0(13,lra)
    end if
    if (linelp == 22) then
      call ttdd_ar_act_blbl_sgt0(13,lra)
    end if

  case (17)
    !===================================================================
    ! tv(17)    act -b^r-
    if ((linelp /= 16) .or. (nlg2 /= 2)) return
    lra = nlg1
    do lri=norb_frz+1,norb_dz-1
      lmi = lsm_inn(lri)
      do lrj=lri+1,norb_dz
        lmj = lsm_inn(lrj)
        lmij = Mul(lmi,lmj)
        if (lmij /= jml) cycle

        w1tv = w1_tv
        ni = mod(lrj-lri,2)
        if (ni == 0) w1tv = -w1tv

        !---------------------------------------------------------------
        ! tv(17) ar(23)-br(23)-
        iwdl = just(lri,lrj)
        iwdr = 0
        do mpl=1,mtype
          vplp_w0(mpl) = Zero
          vplp_w1(mpl) = vplpnew_w1(mpl)*w1tv
        end do
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        ijk = lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz)
        intpos = intind_ijka(ijk)
        call ar_br_br_ext_ar_new(13,intpos,isma)
      end do
    end do

  case (18)
    !===================================================================
    ! t1v(18) arbr act -c"- ext -br-ar -ar  ! no complete
    if ((nlg2 /= 2) .or. (linelp /= 16)) return
    lra = nlg1
    call ttv_arbr_act_c_sgt1(13,lra)

  case (19)
    !===================================================================
    ! dd(19) act -b^l- .................................................
    if (linelp /= 15) return
    lra = nlg1
    if (nlg2 == 1) then
      call dd_drl_act_bl(13,lra)
    else
      call dd_arbl_act_bl(13,lra)
    end if

  case (20)
    !===================================================================
    ! dddd(20) act -b^l-  ..............................................
    if (linelp /= 15) return
    lra = nlg1
    if (nlg2 == 1) then
      call dddd_drl_act_bl_sgt0(13,lra)
    else
      call dddd_arbl_act_bl_sgt0(13,lra)
    end if

  case (21)
    !===================================================================
    ! dd1(21) act -b^l-  ...............................................
    if ((linelp /= 15) .or. (nlg2 /= 2)) return
    lra = nlg1
    call dd1_arbl_act_bl_sgt0(13,lra)

  case (22)
    !===================================================================
    ! d1d(22) act -b^l-  ...............................................
    lra = nlg1
    if ((linelp == 15) .and. (nlg2 == 2)) then
      call d1d_arbl_act_bl_sgt0(13,lra)
      call d1d_drl_act_bl_sgt0(13,lra)
    end if

  case (23)
    !===================================================================
    ! dv(23) act -c'-..................................................
    ! dv(23-1) ar(23)-
    ! dv(23-2) drl(33)-bl(23)-
    if (linelp == 13) then
      imap_1 = 0
      do lrd=norb_frz+1,norb_dz
        lmd = lsm_inn(lrd)
        if (lmd /= jml) cycle
        ni = mod(norb_dz-lrd,2)

        ! dv(23-1) ar(23)-
        w0 = w0_dv(1)
        if (ni == 1) w0 = -w0_dv(1)
        iwdl = jud(lrd)
        iwdr = 0
        imap_1 = imap_1+1
        do mpl=1,mtype
          vplp_w0(mpl) = vplpnew_w0(mpl)*w0
        end do
        do mpl=1,mhlp
          iwal = lpnew_lwei(mpl)
          iwar = lpnew_rwei(mpl)
          lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
        end do
        call ar_td_ext_ar(51,lrd,0,isma)
        call ar_td_ext_rest(lrd)
        ! dv(23-2) drl(33)-bl(23)-
        w0 = w0_dv(2)
        if (ni == 1) w0 = -w0_dv(2)
        do mpl=1,mtype
          vplp_w0(mpl) = vplpnew_w0(mpl)*w0
          vplp_w1(mpl) = Zero
        end do
        do lrk=1,lrd-1
          call drl_bl_ext_ar_new(13,lrk,lrd)
        end do
      end do
    end if
    !===================================================================
    ! dv(23) act -br-br-................................................
    if (linelp == 20) then
      jk = nlg1
      do lrd=norb_frz+1,norb_dz
        lmd = lsm_inn(lrd)
        if (lmd /= jml) cycle
        ni = mod(norb_dz-lrd,2)
        w0 = w0_dv(1)
        if (ni == 1) w0 = -w0_dv(1)
        !...............................................................
        ! dv(23-1) ar(23)-
        iwdl = jud(lrd)
        iwdr = 0
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
        ijk = lrd-norb_frz+jk
        intpos = intind_ijka(ijk)
        call ar_br_br_ext_ar_new(13,intpos,isma)
      end do
    end if
    ! dv(23) act -bl-bl-................................................
    if (linelp /= 22) return
    jk = nlg1
    do lrd=norb_frz+1,norb_dz
      lmd = lsm_inn(lrd)
      if (lmd /= jml) cycle
      ni = mod(norb_dz-lrd,2)
      w0 = w0_dv(1)
      if (ni == 1) w0 = -w0_dv(1)
      ! dv(23-1) ar(23)-
      iwdl = jud(lrd)
      iwdr = 0
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
      ijk = lrd-norb_frz+jk
      intpos = intind_ijka(ijk)
      call ar_bl_bl_ext_ar_new(13,intpos,isma,1)
    end do

  case (24)
    !===================================================================
    ! d1v(24) act -c'-..................................................
    ! d1v(24-1) ar(13)-
    ! d1v(24-2) drl(33)-bl(13)-
    if (linelp == 13) then
      call d1v_ar_act_c_ext(13)
    end if
    if (linelp == 20) then
      jk = nlg1
      call d1v_ar_act_brbr_ext(13,jk)
    end if
    if (linelp == 22) then
      jk = nlg1
      call d1v_ar_act_blbl_ext(13,jk)
    end if

  case (25)
    !===================================================================
    ! vv(25) act -b^l- .................................................
    if (linelp /= 15) return
    lra = nlg1
    if (nlg2 == 1) then
      iwdl = 0
      iwdr = 0
      do mpl=1,mtype
        vplp_w0(mpl) = vplpnew_w0(mpl)*w0_vv
        vplp_w1(mpl) = Zero
      end do
      do mpl=1,mhlp
        iwal = lpnew_lwei(mpl)
        iwar = lpnew_rwei(mpl)
        lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
        lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
      end do
      call drl_bl_sum_ar_new(13,0,0,lra)
      !do lrk=1,norb_dz
      !  call drl_bl_ext_ar_new(13,lrk,lra)
      !end do
    end if

  case (7,14,16,26)
end select

return

end subroutine td_ext_head_in_dbl

subroutine td_ext_head_in_act()

use gugaci_global, only: iml, imr, linelp, logic_dh, nlg1, nlg2
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: intpos, isma, lri, lrj

logic_dh = .false.
lri = nlg1
lrj = nlg2
intpos = nlg1
isma = Mul(iml,imr)

select case (linelp)
  case default ! (1)
    ! line=1 a&r<-->a^r
    call ar_td_ext_ar(100,lri,lrj,isma)
    call ar_td_ext_rest(lri)
  case (4)
    ! line=4 a&r--b&r--b^r<-->a^r
    call ar_br_br_ext_ar_new(13,intpos,isma)
  case (7)
    ! line=7 a&r--b&l--b^l<-->a^r
    call ar_bl_bl_ext_ar_new(13,intpos,isma,1)
  case (10)
    ! line=10 d&rr--b^r<-->a^r
    call drr_br_ext_ar(13,lri,lrj)
  case (12)
    ! line=12 d&rl--b^l<-->a^r
    call drl_bl_ext_ar_new(13,lri,lrj)
  case (2:3,5:6,8:9,11)
end select

return

end subroutine td_ext_head_in_act

subroutine gsd_ext_sequence(iltype,ilsm,irsm,lri)

use gugaci_global, only: ibsm_ext, icano_nnend, icano_nnsta, icnt_base, iesm_ext, iseg_downwei, isegdownwei, m_jc, m_jd, ng_sm
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iltype, ilsm, irsm, lri
integer(kind=iwp) :: ic, icano_nn, icend, icsta, ilnodedownwei, indl, isma, ismb, ismnoded, ismnodes

ismnodes = ilsm
ismnoded = irsm
indl = 0 !?
if (iltype == 2) indl = 1+ismnodes
if (iltype == 3) indl = 9+ismnodes
if (iltype == 4) indl = 17+ismnodes
ilnodedownwei = iseg_downwei(indl)
isegdownwei = ilnodedownwei
icano_nnsta = 1
icnt_base = 0
icsta = ibsm_ext(ismnoded)
icend = iesm_ext(ismnoded)
m_jc = 0

do ic=icsta,icend
  m_jd = ic
  m_jc = ic-icsta+1
  icano_nn = m_jc
  icano_nnend = icano_nn
  do ismb=1,ismnoded-1
    isma = Mul(ismnodes,ismb)
    if (isma > ismb) cycle
    call g31_diffsym(lri,isma,ismb)
  end do

  ismb = ismnoded
  isma = Mul(ismnodes,ismb)
  if (isma == ismb) then
    call gsd_samesym_aaa(lri,isma)
  else if (isma < ismb) then
    call gsd_diffsamesym_abb(lri,isma,ismb)
  end if

  do ismb=ismnoded+1,ng_sm
    isma = Mul(ismnodes,ismb)
    if (isma > ismb) cycle
    if (ismnoded > isma) then
      call g32a_diffsym(lri,isma,ismb)
    else if (ismnoded == isma) then
      call gsd_diffsamesym_aab(lri,isma,ismb)
    else
      call g32b_diffsym(lri,isma,ismb)
    end if
  end do

  if ((ismnodes == 1) .and. (iltype == 4)) then
    call gsd_arlp_s1(lri)
  end if
  icnt_base = icnt_base+ilnodedownwei

end do

end subroutine gsd_ext_sequence

subroutine ar_sd_ext_rest(lri)

use gugaci_global, only: icnt_base, ihy, ihyl, ilsegdownwei, iml, imr, ipae, ipael, irsegdownwei, iseg_downwei, jpad, jpad_upwei, &
                         jpadl, jphy, logic_dh, logic_grad, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, mtype, ndim, nstaval, &
                         nvalue, value_lpext, value_lpext1, vplp_w0, vplpnew_w0, w0_sdplp
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lri
integer(kind=iwp) :: ihypos, iiext, ilpend, ilpsta, ilw, in_, iplp, irw, iw0, iwal, iwal0, iwar, iwar0, iwd, iwuplwei
real(kind=wp) :: w0_old, w0multi
integer(kind=iwp), external :: iwalk_ad

iwuplwei = jpad_upwei(jpadl)
ilsegdownwei = iseg_downwei(ipael)
irsegdownwei = iseg_downwei(ipae)

if (logic_grad) then
  call gsd_ext_sequence_g(4,iml,imr,lri)

  w0_old = One
  do iw0=1,mtype
    w0_sdplp = vplpnew_w0(iw0)
    if (logic_dh) w0_sdplp = vplp_w0(iw0)
    w0multi = w0_sdplp/w0_old
    w0_old = w0_sdplp
    do iiext=1,icnt_base
      value_lpext(iiext) = value_lpext(iiext)*w0multi
      value_lpext1(iiext) = value_lpext1(iiext)*w0multi
    end do
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend
      if (logic_dh) then              !lp_head is in dbl_space
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call complete_sd_ar_ext_loop_g(ilw,irw,ilsegdownwei)
      else
        ihypos = jphy(iplp)
        ndim = ihy(ihypos)
        iwal0 = lpnew_lwei(iplp)
        iwar0 = lpnew_rwei(iplp)
        do in_=1,ndim
          iwal = iwal0+ihyl(ihypos+in_)
          iwar = iwar0+ihy(ihypos+in_)
          do iwd=0,iwuplwei-1
            ilw = iwalk_ad(jpadl,ipael,iwal,iwd)
            irw = iwalk_ad(jpad,ipae,iwar,iwd)
            call complete_sd_ar_ext_loop_g(ilw,irw,ilsegdownwei)
          end do
        end do
      end if
    end do
  end do

else
  call gsd_ext_sequence(4,iml,imr,lri)

  w0_old = One
  do iw0=1,mtype
    w0_sdplp = vplpnew_w0(iw0)
    if (logic_dh) w0_sdplp = vplp_w0(iw0)
    w0multi = w0_sdplp/w0_old
    w0_old = w0_sdplp
    do iiext=1,icnt_base
      value_lpext(iiext) = value_lpext(iiext)*w0multi
    end do
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend
      if (logic_dh) then              !lp_head is in dbl_space
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call complete_sd_ar_ext_loop(ilw,irw,ilsegdownwei)
      else                            !lp_head is in act_s
        ihypos = jphy(iplp)
        ndim = ihy(ihypos)
        iwal0 = lpnew_lwei(iplp)
        iwar0 = lpnew_rwei(iplp)
        do in_=1,ndim
          iwal = iwal0+ihyl(ihypos+in_)
          iwar = iwar0+ihy(ihypos+in_)
          do iwd=0,iwuplwei-1
            ilw = iwalk_ad(jpadl,ipael,iwal,iwd)
            irw = iwalk_ad(jpad,ipae,iwar,iwd)
            call complete_sd_ar_ext_loop(ilw,irw,ilsegdownwei)
          end do
        end do
      end if
    end do
  end do
end if

return

end subroutine ar_sd_ext_rest

subroutine ar_td_ext_rest(lri)

use gugaci_global, only: icnt_base, ihy, ihyl, ilsegdownwei, iml, imr, ipae, ipael, irsegdownwei, iseg_downwei, jpad, jpad_upwei, &
                         jpadl, jphy, logic_dh, logic_grad, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, mtype, ndim, nstaval, &
                         nvalue, value_lpext, value_lpext1, vplp_w0, vplpnew_w0, w0_sdplp
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lri
integer(kind=iwp) :: ihypos, iiext, ilpend, ilpsta, ilw, in_, iplp, irw, iw0, iwal, iwal0, iwar, iwar0, iwd, iwuplwei
real(kind=wp) :: w0_old, w0multi
integer(kind=iwp), external :: iwalk_ad

iwuplwei = jpad_upwei(jpadl)
ilsegdownwei = iseg_downwei(ipael)
irsegdownwei = iseg_downwei(ipae)

if (logic_grad) then
  call gsd_ext_sequence_g(3,iml,imr,lri)

  w0_old = One
  do iw0=1,mtype
    w0_sdplp = vplpnew_w0(iw0)
    if (logic_dh) w0_sdplp = vplp_w0(iw0)
    w0multi = w0_sdplp/w0_old
    w0_old = w0_sdplp
    do iiext=1,icnt_base
      value_lpext(iiext) = value_lpext(iiext)*w0multi
      value_lpext1(iiext) = value_lpext1(iiext)*w0multi
    end do
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend
      if (logic_dh) then              !lp_head is in dbl_space
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call complete_sd_ar_ext_loop_g(ilw,irw,ilsegdownwei)
      else                            !lp_head is in act_spa
        ihypos = jphy(iplp)
        ndim = ihy(ihypos)
        iwal0 = lpnew_lwei(iplp)
        iwar0 = lpnew_rwei(iplp)
        do in_=1,ndim
          iwal = iwal0+ihyl(ihypos+in_)
          iwar = iwar0+ihy(ihypos+in_)
          do iwd=0,iwuplwei-1
            ilw = iwalk_ad(jpadl,ipael,iwal,iwd)
            irw = iwalk_ad(jpad,ipae,iwar,iwd)
            call complete_sd_ar_ext_loop_g(ilw,irw,ilsegdownwei)
          end do
        end do
      end if
    end do
  end do

else

  call gsd_ext_sequence(3,iml,imr,lri)

  w0_old = One
  do iw0=1,mtype
    w0_sdplp = vplpnew_w0(iw0)
    if (logic_dh) w0_sdplp = vplp_w0(iw0)
    w0multi = w0_sdplp/w0_old
    w0_old = w0_sdplp
    do iiext=1,icnt_base
      value_lpext(iiext) = value_lpext(iiext)*w0multi
    end do
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend
      if (logic_dh) then              !lp_head is in dbl_space
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call complete_sd_ar_ext_loop(ilw,irw,ilsegdownwei)
      else                            !lp_head is in act_spa
        ihypos = jphy(iplp)
        ndim = ihy(ihypos)
        iwal0 = lpnew_lwei(iplp)
        iwar0 = lpnew_rwei(iplp)
        do in_=1,ndim
          iwal = iwal0+ihyl(ihypos+in_)
          iwar = iwar0+ihy(ihypos+in_)
          do iwd=0,iwuplwei-1
            ilw = iwalk_ad(jpadl,ipael,iwal,iwd)
            irw = iwalk_ad(jpad,ipae,iwar,iwd)
            call complete_sd_ar_ext_loop(ilw,irw,ilsegdownwei)
          end do
        end do
      end if
    end do
  end do
end if

return

end subroutine ar_td_ext_rest

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

subroutine vd_drt_ci_new()

use gugaci_global, only: idisk_array, idisk_lp, iml, imr, ipae, ipael, jml, jmr, jpad, jpadl, jpadlr, linelp, lpblock_vd
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: jptyl, jptyr, lpb

call external_space_plpmode_value_vd()

idisk_lp = idisk_array(1)
do lpb=1,lpblock_vd
  call read_lp()
  ipael = 1
  ipae = imr+1
  call get_jpty(jpadlr,jptyl,jptyr)
  call get_jp(jptyl,jml,jpadl,1)
  call get_jp(jptyr,jmr,jpad,1)
  call gsd_determine_extarmode_paras(iml,imr,.false.)
  if (linelp <= 12) then
    call vd_ext_head_in_act()
  else
    call vd_ext_head_in_dbl()
  end if
end do

return

end subroutine vd_drt_ci_new

subroutine vd_ext_head_in_dbl()

use gugaci_global, only: ipae, ipael, jb_sys, jml, jmr, jpad, jpadl, jpadlr, linelp, logic_dh, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, mhlp, mtype, nlg1, nlg2, vplp_w0, vplp_w1, vplpnew_w0, w0_vv
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, jk, lpok, lra, mpl
integer(kind=iwp), external :: iwalk_ad

logic_dh = .true.
lpok = jpadlr
!jmlr = Mul(jml,jmr)
select case (lpok)
  case default ! (1)
    !===================================================================
    ! ss(1) 16: act -b^r-
    if (linelp /= 16) return
    lra = nlg1
    if ((nlg2 == 1) .and. (jml == jmr)) then
      call ss_drl_act_br(26,lra)
      if (jb_sys > 0) call ss_drl_act_br_sgt0(26,lra)
    end if
    if (nlg2 == 2) then
      call ss_arbl_act_br(26,lra)
      if (jb_sys > 0) then
        call ss_s_drl_act_br_sgt0(26,lra)
        call ss_arbl_act_br_sgt0(26,lra)
      end if
    end if

  case (2)
    !===================================================================
    ! st(2)  18: act -b^r-
    if ((linelp /= 16) .or. (nlg2 /= 2)) return
    lra = nlg1
    call st_arbl_act_br(26,lra)
    if (jb_sys > 0) call st_arbl_act_br_sgt0(26,lra)
    if (jml == jmr) then
      call st_drl_act_br(26,lra)
      if (jb_sys > 0) call st_drl_act_br_sgt0(26,lra)
    end if

  case (3)
    !===================================================================
    ! ts(3) a&r-b^l-  act -b^r .........................................
    if ((linelp /= 16) .or. (nlg2 /= 2)) return
    lra = nlg1
    call ts_arbl_act_br(26,lra)
    if (jb_sys > 0) call ts_arbl_act_br_sgt0(26,lra)

  case (4)
    !===================================================================
    ! stt(4) ar-bl- act-br
    if ((linelp /= 16) .or. (nlg2 /= 2)) return
    lra = nlg1
    call stt_arbl_act_br_sgt1(26,lra)

  case (5)
    !===================================================================
    ! tts(5) ar-bl- act-br;    drl- act -br
    if ((linelp /= 16) .or. (nlg2 /= 2)) return
    lra = nlg1
    call tts_drl_act_br_sgt1(26,lra)
    call tts_arbl_act_br_sgt1(26,lra)

  case (6)
    !===================================================================
    ! sd(6)  linelp=19  act -d&l^r-
    ! sd(6)  linelp=21  act -b&l-b^r-
    if (linelp == 19) then
      lra = nlg1
      call sd_ar_act_dlr(26,lra)
      if (jb_sys > 0) call sd_ar_act_dlr_sgt0(26,lra)
    end if
    if (linelp == 21) then
      lra = nlg1
      call sd_ar_act_blbr(26,lra)
      if (jb_sys > 0) call sd_ar_act_blbr_sgt0(26,lra)
    end if

  case (7)
    !===================================================================
    ! ds(7) act -c'-
    if (linelp /= 13) return
    if (jmr == 1) call ds_ardlr_act_c1(26)
    call ds_arblbr_act_c1(26)
    if (jb_sys > 0) then
      call ds_arblbr_act_c1_sgt0(26)
    end if

  case (8)
    !===================================================================
    ! sdd(8) ar- act -drl-
    if (linelp == 19) then
      lra = nlg1
      call sdd_ar_act_dlr_sgt0(26,lra)
    end if
    if (linelp == 21) then
      lra = nlg1
      call sdd_ar_act_blbr_sgt0(26,lra)
    end if

  case (9)
    !===================================================================
    ! dds(9) ar-bl-br ar-dlr act -c'-
    if (linelp /= 13) return
    if (jmr == 1) call dds_ardlr_act_c_sgt0(26)
    lra = nlg1
    call dds_drlbr_act_c_sgt0(26)
    call dds_arblbr_act_c_sgt0(26)

  case (11)
    !===================================================================
    ! tt(11) act -b^r
    if (linelp /= 16) return
    lra = nlg1
    if ((nlg2 == 1) .and. (jml == jmr)) then
      call tt_drl_act_br(26,lra)
    end if
    if (nlg2 == 2) then
      call tt_arbl_act_br(26,lra)
    end if

  case (12)
    !===================================================================
    ! tttt(12) act -b^r
    if (linelp /= 16) return
    lra = nlg1
    if ((nlg2 == 1) .and. (jml == jmr)) then
      call tttt_drl_act_br_sgt1(26,lra)
    end if
    if (nlg2 == 2) then
      call tttt_arbl_act_br_sgt1(26,lra)
    end if

  case (13)
    !===================================================================
    ! td(13) 19: act -d&l^r-  21: act -b&l-b^r-
    jk = nlg1
    lra = nlg1
    if (linelp == 19) then
      call td_ar_act_dlr(26,lra)
    end if
    if (linelp == 21) then
      call td_ar_act_blbr(26,jk)
    end if

  case (14)
    !===================================================================
    ! dt(14) 13: act -c'-
    if (linelp /= 13) return
    call dt_arblbr_act_c1(26)

  case (15)
    !===================================================================
    ! ttdd(15) 19: act -d&l^r-  21: act -b&l-b^r-
    lra = nlg1
    if (linelp == 19) then
      call ttdd_ar_act_dlr_sgt1(26,lra)
    end if
    if (linelp == 21) then
      call ttdd_ar_act_blbr_sgt1(26,lra)
    end if

  case (16)
    !===================================================================
    ! ddtt(16) 13: act -c'-
    if (linelp /= 13) return
    call ddtt_arblbr_act_c1_sgt1(26)

  case (19)
    !===================================================================
    ! dd(19) act -br- ..................................................
    if (linelp /= 16) return
    lra = nlg1
    if ((nlg2 == 1) .and. (jml == jmr)) then
      call dd_drl_act_br(26,lra)
    end if
    if (nlg2 == 2) then
      call dd_arbl_act_br(26,lra)
    end if

  case (20)
    !===================================================================
    ! dddd(20) act -br- ................................................
    if (linelp /= 16) return
    lra = nlg1
    if ((nlg2 == 1) .and. (jml == jmr)) then
      call dddd_drl_act_br_sgt0(26,lra)
    end if
    if (nlg2 == 2) then
      call dddd_arbl_act_br_sgt0(26,lra)
    end if

  case (21)
    !===================================================================
    ! dd1(21) 16: ar-bl- act -br- ......................................
    if ((linelp /= 16) .or. (nlg2 /= 2)) return
    lra = nlg1
    call dd1_arbl_act_br_sgt0(26,lra)

  case (22)
    !===================================================================
    ! d1d(22) 16: ar-bl- drl act -br-  .................................
    if ((linelp /= 16) .or. (nlg2 /= 2)) return
    lra = nlg1
    if (jml == jmr) then
      call d1d_drl_act_br_sgt0(26,lra)
    end if
    call d1d_arbl_act_br_sgt0(26,lra)

  case (23)
    !===================================================================
    ! dv(23)  19: act -d&l^r-  21: act -b&l-b^r-
    jk = nlg1
    lra = nlg1
    if (linelp == 19) then
      call dv_ar_act_dlr(26,lra)
    end if
    if (linelp == 21) then
      call dv_ar_act_blbr(26,jk)
    end if

  case (24)
    !===================================================================
    ! d1v(24)  19: act -d&l^r-  21: act -b&l-b^r-
    lra = nlg1
    if (linelp == 19) then
      call ddv_ar_act_dlr_sgt0(26,lra)
    end if
    if (linelp == 21) then
      call ddv_ar_act_blbr_sgt0(26,lra)
    end if

  case (25)
    !===================================================================
    ! vv(25)  d&r^l-  act -b^r
    if ((linelp == 16) .and. (nlg2 == 1)) then
      iwdl = 0
      iwdr = 0
      lra = nlg1
      !lra = nlg2
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
      call drl_br_sum_al_new(26,0,0,lra)
      !do lrk=1,norb_dz
      !  call drl_br_ext_al_new(26,lrk,lra)
      !end do
    end if

  case (10,17:18,26)
end select

return

end subroutine vd_ext_head_in_dbl

subroutine vd_ext_head_in_act()

use gugaci_global, only: imr, linelp, logic_dh, nlg1, nlg2
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: intpos, isma, lri, lrj

logic_dh = .false.

lri = nlg1
lrj = nlg2
intpos = nlg1
isma = imr
select case (linelp)
  case (2)
    ! linelp=2 a&r--d&l^r<-->a^l
    call ar_drl_ext_al_new(26,lri,lrj)
  case (6)
    ! linelp=6 a&r--b&l--b^r<-->a^l
    call ar_bl_br_ext_al_new(26,intpos,isma,1)
  case (11)
    ! linelp=11 d&rl--b^r<-->a^l
    call drl_br_ext_al_new(26,lri,lrj)   ! start
end select

return

end subroutine vd_ext_head_in_act

subroutine ds_drt_ci_new()

use gugaci_global, only: idisk_array, idisk_lp, iml, imr, ipae, ipael, jml, jmr, jpad, jpadl, jpadlr, linelp, lpblock_ds
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: jptyl, jptyr, lpb

!write(u6,*) '  sd_wyb'

call external_space_plpmode_value_ds()

idisk_lp = idisk_array(5)
do lpb=1,lpblock_ds
  call read_lp()
  ipael = iml+1
  ipae = imr+17
  call get_jpty(jpadlr,jptyl,jptyr)
  call get_jp(jptyl,jml,jpadl,1)
  call get_jp(jptyr,jmr,jpad,1)
  !jmlr = Mul(jml,jmr)
  call gsd_determine_extarmode_paras(imr,iml,.true.)
  if (linelp <= 12) then
    call ds_ext_head_in_act()
  else
    call ds_ext_head_in_dbl()
  end if
end do

return

end subroutine ds_drt_ci_new

subroutine ds_ext_head_in_dbl()

use gugaci_global, only: ipae, ipael, jb_sys, jml, jmr, jpad, jpadl, jpadlr, linelp, logic_dh, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, mhlp, mtype, nlg1, nlg2, vplp_w0, vplp_w1, vplpnew_w0, w0_vv
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, jk, lpok, lra, mpl
integer(kind=iwp), external :: iwalk_ad

logic_dh = .true.
!jmlr = Mul(jml,jmr)
lpok = jpadlr
select case (lpok)
  case default ! (1)
    !===================================================================
    ! ss(1) 16: act -b^r-
    if (linelp /= 16) return
    lra = nlg1
    if ((nlg2 == 1) .and. (jml == jmr)) then
      call ss_drl_act_br(7,lra)
      if (jb_sys > 0) call ss_drl_act_br_sgt0(7,lra)
    end if
    if (nlg2 == 2) then
      call ss_arbl_act_br(7,lra)
      if (jb_sys > 0) then
        call ss_s_drl_act_br_sgt0(7,lra)
        call ss_arbl_act_br_sgt0(7,lra)
      end if
    end if

  case (2)
    !===================================================================
    ! st(2)  18: act -b^r-
    if ((linelp /= 16) .or. (nlg2 /= 2)) return
    lra = nlg1
    call st_arbl_act_br(7,lra)
    if (jb_sys > 0) call st_arbl_act_br_sgt0(7,lra)
    if (jml == jmr) then
      call st_drl_act_br(7,lra)
      if (jb_sys > 0) call st_drl_act_br_sgt0(7,lra)
    end if

  case (3)
    !===================================================================
    ! ts(3) a&r-b^l-  act -b^r .........................................
    if ((linelp /= 16) .or. (nlg2 /= 2)) return
    lra = nlg1
    call ts_arbl_act_br(7,lra)
    if (jb_sys > 0) then
      call ts_arbl_act_br_sgt0(7,lra)
    end if

  case (4)
    !===================================================================
    ! st1(4) a&r-b^l-  act -b^r ........................................
    if ((linelp /= 16) .or. (nlg2 /= 2)) return
    lra = nlg1
    call stt_arbl_act_br_sgt1(7,lra)

  case (5)
    !===================================================================
    ! tts(5) ar-bl- act-br;    drl- act -br
    if ((linelp /= 16) .or. (nlg2 /= 2)) return
    lra = nlg1
    call tts_drl_act_br_sgt1(7,lra)
    call tts_arbl_act_br_sgt1(7,lra)

  case (6)
    !===================================================================
    ! sd(6)  linelp=19  act -d&l^r-  linelp=21  act -b&l-b^r-
    if (linelp == 19) then
      lra = nlg1
      call sd_ar_act_dlr(7,lra)
      if (jb_sys > 0) call sd_ar_act_dlr_sgt0(7,lra)
    end if
    if (linelp == 21) then
      jk = nlg1
      call sd_ar_act_blbr(7,jk)
      if (jb_sys > 0) call sd_ar_act_blbr_sgt0(7,jk)
    end if

  case (7)
    !===================================================================
    ! ds(7) act -c'-
    if (linelp /= 13) return
    call ds_arblbr_act_c1(7)
    if (jmr == 1) call ds_ardlr_act_c1(7)
    if (jb_sys > 0) then
      call ds_arblbr_act_c1_sgt0(7)
    end if

  case (8)
    !===================================================================
    ! sdd(8) ar- act -drl-
    if (linelp == 19) then
      lra = nlg1
      call sdd_ar_act_dlr_sgt0(7,lra)
    end if
    if (linelp == 21) then
      lra = nlg1
      call sdd_ar_act_blbr_sgt0(7,lra)
    end if

  case (9)
    !===================================================================
    ! dds(9) ar-bl-br ar-dlr act -c'-
    if (linelp /= 13) return
    if (jmr == 1) call dds_ardlr_act_c_sgt0(7)
    lra = nlg1
    call dds_drlbr_act_c_sgt0(7)
    call dds_arblbr_act_c_sgt0(7)

  case (11)
    !===================================================================
    ! tt(11) act -b^r
    if (linelp /= 16) return
    lra = nlg1
    if ((nlg2 == 1) .and. (jml == jmr)) then
      call tt_drl_act_br(7,lra)
    end if
    if (nlg2 == 2) then
      call tt_arbl_act_br(7,lra)
    end if

  case (12)
    !===================================================================
    ! tttt(12) act -b^r
    if (linelp /= 16) return
    lra = nlg1
    if ((nlg2 == 1) .and. (jml == jmr)) then
      call tttt_drl_act_br_sgt1(7,lra)
    end if
    if (nlg2 == 2) then
      call tttt_arbl_act_br_sgt1(7,lra)
    end if

  case (13)
    !===================================================================
    ! td(13) 19: act -d&l^r-  21: act -b&l-b^r-
    jk = nlg1
    lra = nlg1
    if (linelp == 19) then
      call td_ar_act_dlr(7,lra)
    end if
    if (linelp == 21) then
      call td_ar_act_blbr(7,jk)
    end if

  case (14)
    !===================================================================
    ! dt(14) 13: act -c'-
    if (linelp /= 13) return
    call dt_arblbr_act_c1(7)

  case (15)
    !===================================================================
    ! ttdd(15) 19: act -d&l^r-  21: act -b&l-b^r-
    lra = nlg1
    if (linelp == 19) then
      call ttdd_ar_act_dlr_sgt1(7,lra)
    end if
    if (linelp == 21) then
      call ttdd_ar_act_blbr_sgt1(7,lra)
    end if

  case (16)
    !===================================================================
    ! ddtt(16) 13: act -c'-
    if (linelp /= 13) return
    call ddtt_arblbr_act_c1_sgt1(7)

  case (19)
    !===================================================================
    ! dd(19) act -br- ..................................................
    if (linelp /= 16) return
    lra = nlg1
    if ((nlg2 == 1) .and. (jml == jmr)) then
      call dd_drl_act_br(7,lra)
    end if
    if (nlg2 == 2) then
      call dd_arbl_act_br(7,lra)
    end if

  case (20)
    !===================================================================
    ! dddd(20) act -br- ................................................
    if (linelp /= 16) return
    lra = nlg1
    if ((nlg2 == 1) .and. (jml == jmr)) then
      call dddd_drl_act_br_sgt0(7,lra)
    end if
    if (nlg2 == 2) then
      call dddd_arbl_act_br_sgt0(7,lra)
    end if

  case (21)
    !===================================================================
    ! dd1(21) 16: ar-bl- act -br- ......................................
    if ((linelp /= 16) .or. (nlg2 /= 2)) return
    lra = nlg1
    call dd1_arbl_act_br_sgt0(7,lra)

  case (22)
    !===================================================================
    ! d1d(22) 16: ar-bl- drl act -br-  .................................
    if ((linelp /= 16) .or. (nlg2 /= 2)) return
    lra = nlg1
    if (jml == jmr) then
      call d1d_drl_act_br_sgt0(7,lra)
    end if
    call d1d_arbl_act_br_sgt0(7,lra)

  case (23)
    !===================================================================
    ! dv(23)  19: act -d&l^r-  21: act -b&l-b^r-
    jk = nlg1
    lra = nlg1
    if (linelp == 19) then
      call dv_ar_act_dlr(7,lra)
    end if
    if (linelp == 21) then
      call dv_ar_act_blbr(7,jk)
    end if

  case (24)
    !===================================================================
    ! d1v(24)  19: act -d&l^r-  21: act -b&l-b^r-
    lra = nlg1
    if (linelp == 19) then
      call ddv_ar_act_dlr_sgt0(7,lra)
    end if
    if (linelp == 21) then
      call ddv_ar_act_blbr_sgt0(7,lra)
    end if

  case (25)
    !===================================================================
    ! vv(25)  d&r^l-  act -b^r
    if ((linelp == 16) .and. (nlg2 == 1)) then
      iwdl = 0
      iwdr = 0
      lra = nlg1
      !lra = nlg2
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
      call drl_br_sum_al_new(7,0,0,lra)
      !do lrk=1,norb_dz
      !  call drl_br_ext_al_new(7,lrk,lra)
      !end do
    end if

  case (10,17:18,26)
end select

return

end subroutine ds_ext_head_in_dbl

subroutine ds_ext_head_in_act()

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
  case default ! (2)
    ! linelp=2 a&r--d&l^r<-->a^l
    call ar_drl_ext_al_new(7,lri,lrj)
  case (6)
    ! linelp=6 a&r--b&l--b^r<-->a^l
    call ar_bl_br_ext_al_new(7,intpos,isma,1)
  case (11)
    ! linelp=11 d&rl--b^r<-->a^l
    call drl_br_ext_al_new(7,lri,lrj)   ! start
  case (1,3:5,7:10,12)
end select

return

end subroutine ds_ext_head_in_act

subroutine dt_drt_ci_new()

use gugaci_global, only: idisk_array, idisk_lp, iml, imr, ipae, ipael, jml, jmr, jpad, jpadl, jpadlr, linelp, lpblock_dt
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: jptyl, jptyr, lpb

!write(u6,*) '  td_wyb'
call external_space_plpmode_value_dt()

idisk_lp = idisk_array(4)
do lpb=1,lpblock_dt
  call read_lp()
  ipael = iml+1
  ipae = imr+9
  call get_jpty(jpadlr,jptyl,jptyr)
  call get_jp(jptyl,jml,jpadl,1)
  call get_jp(jptyr,jmr,jpad,1)
  !jmlr = Mul(jml,jmr)
  call gsd_determine_extarmode_paras(imr,iml,.false.)
  if (linelp <= 12) then
    call dt_ext_head_in_act()
  else
    call dt_ext_head_in_dbl()
  end if
end do

return

end subroutine dt_drt_ci_new

subroutine dt_ext_head_in_dbl()

use gugaci_global, only: ipae, ipael, jb_sys, jml, jmr, jpad, jpadl, jpadlr, linelp, logic_dh, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, mhlp, mtype, nlg1, nlg2, vplp_w0, vplp_w1, vplpnew_w0, w0_vv
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, jk, lpok, lra, mpl
integer(kind=iwp), external :: iwalk_ad

logic_dh = .true.
!jmlr = Mul(jml,jmr)
lpok = jpadlr
select case (lpok)
  case default ! (1)
    !===================================================================
    ! ss(1) 16: act -b^r-
    if (linelp /= 16) return
    lra = nlg1
    if ((nlg2 == 1) .and. (jml == jmr)) then
      call ss_drl_act_br(14,lra)
      if (jb_sys > 0) call ss_drl_act_br_sgt0(14,lra)
    end if
    if (nlg2 == 2) then
      call ss_arbl_act_br(14,lra)
      if (jb_sys > 0) then
        call ss_s_drl_act_br_sgt0(14,lra)
        call ss_arbl_act_br_sgt0(14,lra)
      end if
    end if

  case (2)
    !===================================================================
    ! st(2)  16: act -b^r-
    if ((linelp /= 16) .or. (nlg2 /= 2)) return
    lra = nlg1
    call st_arbl_act_br(14,lra)
    if (jb_sys > 0) call st_arbl_act_br_sgt0(14,lra)
    if (jml == jmr) then
      call st_drl_act_br(14,lra)
      if (jb_sys > 0) call st_drl_act_br_sgt0(14,lra)
    end if

  case (3)
    !===================================================================
    ! ts(3) a&r-b^l-  act -b^r .........................................
    if ((linelp /= 16) .or. (nlg2 /= 2)) return
    lra = nlg1
    call ts_arbl_act_br(14,lra)
    if (jb_sys > 0) then
      call ts_arbl_act_br_sgt0(14,lra)
    end if

  case (4)
    !===================================================================
    ! st1(4) a&r-b^l-  act -b^r ........................................
    if ((linelp /= 16) .or. (nlg2 /= 2)) return
    lra = nlg1
    call stt_arbl_act_br_sgt1(14,lra)

  case (5)
    !===================================================================
    ! tts(5) ar-bl- act-br;    drl- act -br
    if ((linelp /= 16) .or. (nlg2 /= 2)) return
    lra = nlg1
    call tts_drl_act_br_sgt1(14,lra)
    call tts_arbl_act_br_sgt1(14,lra)

  case (6)
    !===================================================================
    ! sd(6)  linelp=19  act -d&l^r-  linelp=21  act -b&l-b^r-
    if (linelp == 19) then
      lra = nlg1
      call sd_ar_act_dlr(14,lra)
      if (jb_sys > 0) call sd_ar_act_dlr_sgt0(14,lra)
    end if
    if (linelp == 21) then
      jk = nlg1
      call sd_ar_act_blbr(14,jk)
      if (jb_sys > 0) call sd_ar_act_blbr_sgt0(14,jk)
    end if

  case (7)
    !===================================================================
    ! ds(7) act -c'-
    if (linelp /= 13) return
    if (jmr == 1) call ds_ardlr_act_c1(14)
    call ds_arblbr_act_c1(14)
    if (jb_sys > 0) then
      call ds_arblbr_act_c1_sgt0(14)
    end if

  case (8)
    !===================================================================
    ! sdd(8) ar- act -drl-
    if (linelp == 19) then
      lra = nlg1
      call sdd_ar_act_dlr_sgt0(14,lra)
    end if
    if (linelp == 21) then
      lra = nlg1
      call sdd_ar_act_blbr_sgt0(14,lra)
    end if

  case (9)
    !===================================================================
    ! dds(9) ar-bl-br ar-dlr act -c'-
    if (linelp /= 13) return
    if (jmr == 1) call dds_ardlr_act_c_sgt0(14)
    lra = nlg1
    call dds_drlbr_act_c_sgt0(14)
    call dds_arblbr_act_c_sgt0(14)

  case (11)
    !===================================================================
    ! tt(11) act -b^r
    if (linelp /= 16) return
    lra = nlg1
    if ((nlg2 == 1) .and. (jml == jmr)) then
      call tt_drl_act_br(14,lra)
    end if
    if (nlg2 == 2) then
      call tt_arbl_act_br(14,lra)
    end if

  case (12)
    !===================================================================
    ! tttt(12) act -b^r
    if (linelp /= 16) return
    lra = nlg1
    if ((nlg2 == 1) .and. (jml == jmr)) then
      call tttt_drl_act_br_sgt1(14,lra)
    end if
    if (nlg2 == 2) then
      call tttt_arbl_act_br_sgt1(14,lra)
    end if

  case (13)
    !===================================================================
    ! td(13) 19: act -d&l^r-  21: act -b&l-b^r-
    jk = nlg1
    lra = nlg1
    if (linelp == 19) then
      call td_ar_act_dlr(14,lra)
    end if
    if (linelp == 21) then
      call td_ar_act_blbr(14,jk)
    end if

  case (14)
    !===================================================================
    ! dt(14) 13: act -c'-
    if (linelp /= 13) return
    call dt_arblbr_act_c1(14)

  case (15)
    !===================================================================
    ! ttdd(15) 19: act -d&l^r-  21: act -b&l-b^r-
    lra = nlg1
    if (linelp == 19) then
      call ttdd_ar_act_dlr_sgt1(14,lra)
    end if
    if (linelp == 21) then
      call ttdd_ar_act_blbr_sgt1(14,lra)
    end if

  case (16)
    !===================================================================
    ! ddtt(16) 13: act -c'-
    if (linelp /= 13) return
    call ddtt_arblbr_act_c1_sgt1(14)

  case (19)
    !===================================================================
    ! dd(19) act -br- ..................................................
    if (linelp /= 16) return
    lra = nlg1
    if ((nlg2 == 1) .and. (jml == jmr)) then
      call dd_drl_act_br(14,lra)
    end if
    if (nlg2 == 2) then
      call dd_arbl_act_br(14,lra)
    end if

  case (20)
    !===================================================================
    ! dddd(20) act -br- ................................................
    if (linelp /= 16) return
    lra = nlg1
    if (nlg2 == 1) then
      call dddd_drl_act_br_sgt0(14,lra)
    end if
    if (nlg2 == 2) then
      call dddd_arbl_act_br_sgt0(14,lra)
    end if

  case (21)
    !===================================================================
    ! dd1(21) 16: ar-bl- act -br- ......................................
    if ((linelp /= 16) .or. (nlg2 /= 2)) return
    lra = nlg1
    call dd1_arbl_act_br_sgt0(14,lra)

  case (22)
    !===================================================================
    ! d1d(22) 16: ar-bl- drl act -br-  .................................
    if ((linelp /= 16) .or. (nlg2 /= 2)) return
    lra = nlg1
    if (jml == jmr) then
      call d1d_drl_act_br_sgt0(14,lra)
    end if
    call d1d_arbl_act_br_sgt0(14,lra)

  case (23)
    !===================================================================
    ! dv(23)  19: act -d&l^r-  21: act -b&l-b^r-
    jk = nlg1
    lra = nlg1
    if (linelp == 19) then
      call dv_ar_act_dlr(14,lra)
    end if
    if (linelp == 21) then
      call dv_ar_act_blbr(14,jk)
    end if

  case (24)
    !===================================================================
    ! d1v(24)  19: act -d&l^r-  21: act -b&l-b^r-
    lra = nlg1
    if (linelp == 19) then
      call ddv_ar_act_dlr_sgt0(14,lra)
    end if
    if (linelp == 21) then
      call ddv_ar_act_blbr_sgt0(14,lra)
    end if

  case (25)
    !===================================================================
    ! vv(25)  d&r^l-  act -b^r
    if ((linelp == 16) .and. (nlg2 == 1)) then
      iwdl = 0
      iwdr = 0
      lra = nlg1
      !lra = nlg2
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
      call drl_br_sum_al_new(14,0,0,lra)
      !do lrk=1,norb_dz
      !  call drl_br_ext_al_new(14,lrk,lra)
      !end do
    end if

  case (10,17:18,26)
end select

return

end subroutine dt_ext_head_in_dbl

subroutine dt_ext_head_in_act()

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
  case default ! (2)
    ! linelp=2 a&r--d&l^r<-->a^l
    call ar_drl_ext_al_new(14,lri,lrj)
  case (6)
    ! linelp=6 a&r--b&l--b^r<-->a^l
    call ar_bl_br_ext_al_new(14,intpos,isma,1)
  case (11)
    ! linelp=11 d&rl--b^r<-->a^l
    call drl_br_ext_al_new(14,lri,lrj)   ! start
  case (1,3:5,7:10,12)
end select

return

end subroutine dt_ext_head_in_act

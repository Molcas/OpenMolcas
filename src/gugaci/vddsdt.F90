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

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
#include "lpdisk.fh"

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

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

logic_dh = .true.
lpok = jpadlr
!jmlr = mul_tab(jml,jmr)
goto(101,102,103,104,105,106,107,108,109,10,111,112,113,114,115,116,10,10,119,120,121,122,123,124,125,10),lpok
!=======================================================================
! ss(1) 16: act -b^r-
101 continue
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
return
!=======================================================================
! st(2)  18: act -b^r-
102 continue
if ((linelp /= 16) .or. (nlg2 /= 2)) return
lra = nlg1
call st_arbl_act_br(26,lra)
if (jb_sys > 0) call st_arbl_act_br_sgt0(26,lra)
if (jml == jmr) then
  call st_drl_act_br(26,lra)
  if (jb_sys > 0) call st_drl_act_br_sgt0(26,lra)
end if
return
!=======================================================================
! ts(3) a&r-b^l-  act -b^r ............................................
103 continue
if ((linelp /= 16) .or. (nlg2 /= 2)) return
lra = nlg1
call ts_arbl_act_br(26,lra)
if (jb_sys > 0) call ts_arbl_act_br_sgt0(26,lra)
return
!=======================================================================
! stt(4) ar-bl- act-br
104 continue
if ((linelp /= 16) .or. (nlg2 /= 2)) return
lra = nlg1
call stt_arbl_act_br_sgt1(26,lra)
return

!=======================================================================
! tts(5) ar-bl- act-br;    drl- act -br
105 continue
if ((linelp /= 16) .or. (nlg2 /= 2)) return
lra = nlg1
call tts_drl_act_br_sgt1(26,lra)
call tts_arbl_act_br_sgt1(26,lra)
return
!=======================================================================
! sd(6)  linelp=19  act -d&l^r-
! sd(6)  linelp=21  act -b&l-b^r-
106 continue
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
return
!=======================================================================
! ds(7) act -c'-
107 continue
if (linelp /= 13) return
if (jmr == 1) call ds_ardlr_act_c1(26)
call ds_arblbr_act_c1(26)
if (jb_sys > 0) then
  call ds_arblbr_act_c1_sgt0(26)
end if
return
!=======================================================================
! sdd(8) ar- act -drl-
108 continue
if (linelp == 19) then
  lra = nlg1
  call sdd_ar_act_dlr_sgt0(26,lra)
end if
if (linelp == 21) then
  lra = nlg1
  call sdd_ar_act_blbr_sgt0(26,lra)
end if
return
!=======================================================================
! dds(9) ar-bl-br ar-dlr act -c'-
109 continue
if (linelp /= 13) return
if (jmr == 1) call dds_ardlr_act_c_sgt0(26)
lra = nlg1
call dds_drlbr_act_c_sgt0(26)
call dds_arblbr_act_c_sgt0(26)
return
!=======================================================================
! tt(11) act -b^r
111 continue
if (linelp /= 16) return
lra = nlg1
if ((nlg2 == 1) .and. (jml == jmr)) then
  call tt_drl_act_br(26,lra)
end if
if (nlg2 == 2) then
  call tt_arbl_act_br(26,lra)
end if
return
!=======================================================================
! tttt(12) act -b^r
112 continue
if (linelp /= 16) return
lra = nlg1
if ((nlg2 == 1) .and. (jml == jmr)) then
  call tttt_drl_act_br_sgt1(26,lra)
end if
if (nlg2 == 2) then
  call tttt_arbl_act_br_sgt1(26,lra)
end if
return
!=======================================================================
! td(13) 19: act -d&l^r-  21: act -b&l-b^r-
113 continue
jk = nlg1
lra = nlg1
if (linelp == 19) then
  call td_ar_act_dlr(26,lra)
end if
if (linelp == 21) then
  call td_ar_act_blbr(26,jk)
end if
return
!=======================================================================
! dt(14) 13: act -c'-
114 continue
if (linelp /= 13) return
call dt_arblbr_act_c1(26)
return
!=======================================================================
! ttdd(15) 19: act -d&l^r-  21: act -b&l-b^r-
115 continue
lra = nlg1
if (linelp == 19) then
  call ttdd_ar_act_dlr_sgt1(26,lra)
end if
if (linelp == 21) then
  call ttdd_ar_act_blbr_sgt1(26,lra)
end if
return
!=======================================================================
! ddtt(16) 13: act -c'-
116 continue
if (linelp /= 13) return
call ddtt_arblbr_act_c1_sgt1(26)
return
!=======================================================================
! dd(19) act -br- ....................................................
119 continue
if (linelp /= 16) return
lra = nlg1
if ((nlg2 == 1) .and. (jml == jmr)) then
  call dd_drl_act_br(26,lra)
end if
if (nlg2 == 2) then
  call dd_arbl_act_br(26,lra)
end if
return
!=======================================================================
! dddd(20) act -br- ....................................................
120 continue
if (linelp /= 16) return
lra = nlg1
if ((nlg2 == 1) .and. (jml == jmr)) then
  call dddd_drl_act_br_sgt0(26,lra)
end if
if (nlg2 == 2) then
  call dddd_arbl_act_br_sgt0(26,lra)
end if
return
!=======================================================================
! dd1(21) 16: ar-bl- act -br- ...........................................
121 continue
if ((linelp /= 16) .or. (nlg2 /= 2)) return
lra = nlg1
call dd1_arbl_act_br_sgt0(26,lra)
return
!=======================================================================
! d1d(22) 16: ar-bl- drl act -br-  ......................................
122 continue
if ((linelp /= 16) .or. (nlg2 /= 2)) return
lra = nlg1
if (jml == jmr) then
  call d1d_drl_act_br_sgt0(26,lra)
end if
call d1d_arbl_act_br_sgt0(26,lra)
return
!=======================================================================
! dv(23)  19: act -d&l^r-  21: act -b&l-b^r-
123 continue
jk = nlg1
lra = nlg1
if (linelp == 19) then
  call dv_ar_act_dlr(26,lra)
end if
if (linelp == 21) then
  call dv_ar_act_blbr(26,jk)
end if
return
!=======================================================================
! d1v(24)  19: act -d&l^r-  21: act -b&l-b^r-
124 continue
lra = nlg1
if (linelp == 19) then
  call ddv_ar_act_dlr_sgt0(26,lra)
end if
if (linelp == 21) then
  call ddv_ar_act_blbr_sgt0(26,lra)
end if
return
!=======================================================================
! vv(25)  d&r^l-  act -b^r
125 continue
if ((linelp == 16) .and. (nlg2 == 1)) then
  iwdl = 0
  iwdr = 0
  lra = nlg1
  !lra = nlg2
  do mpl=1,mtype
    vplp_w0(mpl) = vplpnew_w0(mpl)*w0_vv
    vplp_w1(mpl) = 0.d0
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
return
!=======================================================================

10 continue
return

end subroutine vd_ext_head_in_dbl

subroutine vd_ext_head_in_act()

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

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

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
#include "lpdisk.fh"

!write(6,*) '  sd_wyb'

call external_space_plpmode_value_ds()

idisk_lp = idisk_array(5)
do lpb=1,lpblock_ds
  call read_lp()
  ipael = iml+1
  ipae = imr+17
  call get_jpty(jpadlr,jptyl,jptyr)
  call get_jp(jptyl,jml,jpadl,1)
  call get_jp(jptyr,jmr,jpad,1)
  !jmlr = mul_tab(jml,jmr)
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

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

logic_dh = .true.
!jmlr = mul_tab(jml,jmr)
lpok = jpadlr
!goto(101,102,103,104,10,106,107,108,109,10,111,112,113,114,115,116,10,10,119,120,121,122,123,124,125,10),lpok
goto(101,102,103,104,105,106,107,108,109,10,111,112,113,114,115,116,10,10,119,120,121,122,123,124,125,10),lpok
!=======================================================================
! ss(1) 16: act -b^r-
101 continue
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
return
!=======================================================================
! st(2)  18: act -b^r-
102 continue
if ((linelp /= 16) .or. (nlg2 /= 2)) return
lra = nlg1
call st_arbl_act_br(7,lra)
if (jb_sys > 0) call st_arbl_act_br_sgt0(7,lra)
if (jml == jmr) then
  call st_drl_act_br(7,lra)
  if (jb_sys > 0) call st_drl_act_br_sgt0(7,lra)
end if
return
!=======================================================================
! ts(3) a&r-b^l-  act -b^r ............................................
103 continue
if ((linelp /= 16) .or. (nlg2 /= 2)) return
lra = nlg1
call ts_arbl_act_br(7,lra)
if (jb_sys > 0) then
  call ts_arbl_act_br_sgt0(7,lra)
end if
return
!=======================================================================
! st1(4) a&r-b^l-  act -b^r ............................................
104 continue
if ((linelp /= 16) .or. (nlg2 /= 2)) return
lra = nlg1
call stt_arbl_act_br_sgt1(7,lra)
return
!=======================================================================
! tts(5) ar-bl- act-br;    drl- act -br
105 continue
if ((linelp /= 16) .or. (nlg2 /= 2)) return
lra = nlg1
call tts_drl_act_br_sgt1(7,lra)
call tts_arbl_act_br_sgt1(7,lra)
return
!=======================================================================
! sd(6)  linelp=19  act -d&l^r-  linelp=21  act -b&l-b^r-
106 continue
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
return
!=======================================================================
! ds(7) act -c'-
107 continue
if (linelp /= 13) return
call ds_arblbr_act_c1(7)
if (jmr == 1) call ds_ardlr_act_c1(7)
if (jb_sys > 0) then
  call ds_arblbr_act_c1_sgt0(7)
end if
return
!=======================================================================
! sdd(8) ar- act -drl-
108 continue
if (linelp == 19) then
  lra = nlg1
  call sdd_ar_act_dlr_sgt0(7,lra)
end if
if (linelp == 21) then
  lra = nlg1
  call sdd_ar_act_blbr_sgt0(7,lra)
end if
return
!=======================================================================
! dds(9) ar-bl-br ar-dlr act -c'-
109 continue
if (linelp /= 13) return
if (jmr == 1) call dds_ardlr_act_c_sgt0(7)
lra = nlg1
call dds_drlbr_act_c_sgt0(7)
call dds_arblbr_act_c_sgt0(7)
return
!=======================================================================
! tt(11) act -b^r
111 continue
if (linelp /= 16) return
lra = nlg1
if ((nlg2 == 1) .and. (jml == jmr)) then
  call tt_drl_act_br(7,lra)
end if
if (nlg2 == 2) then
  call tt_arbl_act_br(7,lra)
end if
return
!=======================================================================
! tttt(12) act -b^r
112 continue
if (linelp /= 16) return
lra = nlg1
if ((nlg2 == 1) .and. (jml == jmr)) then
  call tttt_drl_act_br_sgt1(7,lra)
end if
if (nlg2 == 2) then
  call tttt_arbl_act_br_sgt1(7,lra)
end if
return
!=======================================================================
! td(13) 19: act -d&l^r-  21: act -b&l-b^r-
113 continue
jk = nlg1
lra = nlg1
if (linelp == 19) then
  call td_ar_act_dlr(7,lra)
end if
if (linelp == 21) then
  call td_ar_act_blbr(7,jk)
end if
return
!=======================================================================
! dt(14) 13: act -c'-
114 continue
if (linelp /= 13) return
call dt_arblbr_act_c1(7)
return
!=======================================================================
! ttdd(15) 19: act -d&l^r-  21: act -b&l-b^r-
115 continue
lra = nlg1
if (linelp == 19) then
  call ttdd_ar_act_dlr_sgt1(7,lra)
end if
if (linelp == 21) then
  call ttdd_ar_act_blbr_sgt1(7,lra)
end if
return
!=======================================================================
! ddtt(16) 13: act -c'-
116 continue
if (linelp /= 13) return
call ddtt_arblbr_act_c1_sgt1(7)
return
!=======================================================================
! dd(19) act -br- ....................................................
119 continue
if (linelp /= 16) return
lra = nlg1
if ((nlg2 == 1) .and. (jml == jmr)) then
  call dd_drl_act_br(7,lra)
end if
if (nlg2 == 2) then
  call dd_arbl_act_br(7,lra)
end if
return
!=======================================================================
! dddd(20) act -br- ....................................................
120 continue
if (linelp /= 16) return
lra = nlg1
if ((nlg2 == 1) .and. (jml == jmr)) then
  call dddd_drl_act_br_sgt0(7,lra)
end if
if (nlg2 == 2) then
  call dddd_arbl_act_br_sgt0(7,lra)
end if
return
!=======================================================================
! dd1(21) 16: ar-bl- act -br- ...........................................
121 continue
if ((linelp /= 16) .or. (nlg2 /= 2)) return
lra = nlg1
call dd1_arbl_act_br_sgt0(7,lra)
return
!=======================================================================
! d1d(22) 16: ar-bl- drl act -br-  ......................................
122 continue
if ((linelp /= 16) .or. (nlg2 /= 2)) return
lra = nlg1
if (jml == jmr) then
  call d1d_drl_act_br_sgt0(7,lra)
end if
call d1d_arbl_act_br_sgt0(7,lra)
return
!=======================================================================
! dv(23)  19: act -d&l^r-  21: act -b&l-b^r-
123 continue
jk = nlg1
lra = nlg1
if (linelp == 19) then
  call dv_ar_act_dlr(7,lra)
end if
if (linelp == 21) then
  call dv_ar_act_blbr(7,jk)
end if
return
!=======================================================================
! d1v(24)  19: act -d&l^r-  21: act -b&l-b^r-
124 continue
lra = nlg1
if (linelp == 19) then
  call ddv_ar_act_dlr_sgt0(7,lra)
end if
if (linelp == 21) then
  call ddv_ar_act_blbr_sgt0(7,lra)
end if
return
!=======================================================================
! vv(25)  d&r^l-  act -b^r
125 continue
if ((linelp == 16) .and. (nlg2 == 1)) then
  iwdl = 0
  iwdr = 0
  lra = nlg1
  !lra = nlg2
  do mpl=1,mtype
    vplp_w0(mpl) = vplpnew_w0(mpl)*w0_vv
    vplp_w1(mpl) = 0.d0
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

10 continue
return

end subroutine ds_ext_head_in_dbl

subroutine ds_ext_head_in_act()

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

logic_dh = .false.
lri = nlg1
lrj = nlg2
intpos = nlg1
isma = mul_tab(iml,imr)
goto(1,2,3,4,5,6,7,8,9,10,11,12),linelp
! linelp=2 a&r--d&l^r<-->a^l
2 continue
call ar_drl_ext_al_new(7,lri,lrj)
goto 100
! linelp=6 a&r--b&l--b^r<-->a^l
6 continue
call ar_bl_br_ext_al_new(7,intpos,isma,1)
goto 100
! linelp=11 d&rl--b^r<-->a^l
11 continue
call drl_br_ext_al_new(7,lri,lrj)   ! start
goto 100
1 continue
goto 100
3 continue
goto 100
4 continue
goto 100
5 continue
goto 100
7 continue
goto 100
8 continue
goto 100
9 continue
goto 100
10 continue
goto 100
12 continue
goto 100

100 continue
return

end subroutine ds_ext_head_in_act

subroutine dt_drt_ci_new()

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
#include "lpdisk.fh"

!write(6,*) '  td_wyb'
call external_space_plpmode_value_dt()

idisk_lp = idisk_array(4)
do lpb=1,lpblock_dt
  call read_lp()
  ipael = iml+1
  ipae = imr+9
  call get_jpty(jpadlr,jptyl,jptyr)
  call get_jp(jptyl,jml,jpadl,1)
  call get_jp(jptyr,jmr,jpad,1)
  !jmlr = mul_tab(jml,jmr)
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

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

logic_dh = .true.
!jmlr = mul_tab(jml,jmr)
lpok = jpadlr
!goto(101,102,103,10,10,106,107,10,10,10,111,10,113,114,10,10,10,10,119,10,10,10,123,10,125,10),lpok
goto(101,102,103,104,105,106,107,108,109,10,111,112,113,114,115,116,10,10,119,120,121,122,123,124,125,10),lpok
!=======================================================================
! ss(1) 16: act -b^r-
101 continue
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
return
!=======================================================================
! st(2)  16: act -b^r-
102 continue
if ((linelp /= 16) .or. (nlg2 /= 2)) return
lra = nlg1
call st_arbl_act_br(14,lra)
if (jb_sys > 0) call st_arbl_act_br_sgt0(14,lra)
if (jml == jmr) then
  call st_drl_act_br(14,lra)
  if (jb_sys > 0) call st_drl_act_br_sgt0(14,lra)
end if
return
!=======================================================================
! ts(3) a&r-b^l-  act -b^r ............................................
103 continue
if ((linelp /= 16) .or. (nlg2 /= 2)) return
lra = nlg1
call ts_arbl_act_br(14,lra)
if (jb_sys > 0) then
  call ts_arbl_act_br_sgt0(14,lra)
end if
return
!=======================================================================
! st1(4) a&r-b^l-  act -b^r ............................................
104 continue
if ((linelp /= 16) .or. (nlg2 /= 2)) return
lra = nlg1
call stt_arbl_act_br_sgt1(14,lra)
return
!=======================================================================
! tts(5) ar-bl- act-br;    drl- act -br
105 continue
if ((linelp /= 16) .or. (nlg2 /= 2)) return
lra = nlg1
call tts_drl_act_br_sgt1(14,lra)
call tts_arbl_act_br_sgt1(14,lra)
return
!=======================================================================
! sd(6)  linelp=19  act -d&l^r-  linelp=21  act -b&l-b^r-
106 continue
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
return
!=======================================================================
! ds(7) act -c'-
107 continue
if (linelp /= 13) return
if (jmr == 1) call ds_ardlr_act_c1(14)
call ds_arblbr_act_c1(14)
if (jb_sys > 0) then
  call ds_arblbr_act_c1_sgt0(14)
end if
return
!=======================================================================
! sdd(8) ar- act -drl-
108 continue
if (linelp == 19) then
  lra = nlg1
  call sdd_ar_act_dlr_sgt0(14,lra)
end if
if (linelp == 21) then
  lra = nlg1
  call sdd_ar_act_blbr_sgt0(14,lra)
end if
return
!=======================================================================
! dds(9) ar-bl-br ar-dlr act -c'-
109 continue
if (linelp /= 13) return
if (jmr == 1) call dds_ardlr_act_c_sgt0(14)
lra = nlg1
call dds_drlbr_act_c_sgt0(14)
call dds_arblbr_act_c_sgt0(14)
return
!=======================================================================
! tt(11) act -b^r
111 continue
if (linelp /= 16) return
lra = nlg1
if ((nlg2 == 1) .and. (jml == jmr)) then
  call tt_drl_act_br(14,lra)
end if
if (nlg2 == 2) then
  call tt_arbl_act_br(14,lra)
end if
return
!=======================================================================
! tttt(12) act -b^r
112 continue
if (linelp /= 16) return
lra = nlg1
if ((nlg2 == 1) .and. (jml == jmr)) then
  call tttt_drl_act_br_sgt1(14,lra)
end if
if (nlg2 == 2) then
  call tttt_arbl_act_br_sgt1(14,lra)
end if
return
!=======================================================================
! td(13) 19: act -d&l^r-  21: act -b&l-b^r-
113 continue
jk = nlg1
lra = nlg1
if (linelp == 19) then
  call td_ar_act_dlr(14,lra)
end if
if (linelp == 21) then
  call td_ar_act_blbr(14,jk)
end if
return
!=======================================================================
! dt(14) 13: act -c'-
114 continue
if (linelp /= 13) return
call dt_arblbr_act_c1(14)
return
!=======================================================================
! ttdd(15) 19: act -d&l^r-  21: act -b&l-b^r-
115 continue
lra = nlg1
if (linelp == 19) then
  call ttdd_ar_act_dlr_sgt1(14,lra)
end if
if (linelp == 21) then
  call ttdd_ar_act_blbr_sgt1(14,lra)
end if
return
!=======================================================================
! ddtt(16) 13: act -c'-
116 continue
if (linelp /= 13) return
call ddtt_arblbr_act_c1_sgt1(14)
return
!=======================================================================
! dd(19) act -br- ....................................................
119 continue
if (linelp /= 16) return
lra = nlg1
if ((nlg2 == 1) .and. (jml == jmr)) then
  call dd_drl_act_br(14,lra)
end if
if (nlg2 == 2) then
  call dd_arbl_act_br(14,lra)
end if
return
!=======================================================================
! dddd(20) act -br- ....................................................
120 continue
if (linelp /= 16) return
lra = nlg1
if (nlg2 == 1) then
  call dddd_drl_act_br_sgt0(14,lra)
end if
if (nlg2 == 2) then
  call dddd_arbl_act_br_sgt0(14,lra)
end if
return
!=======================================================================
! dd1(21) 16: ar-bl- act -br- ..........................................
121 continue
if ((linelp /= 16) .or. (nlg2 /= 2)) return
lra = nlg1
call dd1_arbl_act_br_sgt0(14,lra)
return
!=======================================================================
! d1d(22) 16: ar-bl- drl act -br-  .....................................
122 continue
if ((linelp /= 16) .or. (nlg2 /= 2)) return
lra = nlg1
if (jml == jmr) then
  call d1d_drl_act_br_sgt0(14,lra)
end if
call d1d_arbl_act_br_sgt0(14,lra)
return
!=======================================================================
! dv(23)  19: act -d&l^r-  21: act -b&l-b^r-
123 continue
jk = nlg1
lra = nlg1
if (linelp == 19) then
  call dv_ar_act_dlr(14,lra)
end if
if (linelp == 21) then
  call dv_ar_act_blbr(14,jk)
end if
return
!=======================================================================
! d1v(24)  19: act -d&l^r-  21: act -b&l-b^r-
124 continue
lra = nlg1
if (linelp == 19) then
  call ddv_ar_act_dlr_sgt0(14,lra)
end if
if (linelp == 21) then
  call ddv_ar_act_blbr_sgt0(14,lra)
end if
return
!=======================================================================
! vv(25)  d&r^l-  act -b^r
125 continue
if ((linelp == 16) .and. (nlg2 == 1)) then
  iwdl = 0
  iwdr = 0
  lra = nlg1
  !lra = nlg2
  do mpl=1,mtype
    vplp_w0(mpl) = vplpnew_w0(mpl)*w0_vv
    vplp_w1(mpl) = 0.d0
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

10 continue
return

end subroutine dt_ext_head_in_dbl

subroutine dt_ext_head_in_act()

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

logic_dh = .false.
lri = nlg1
lrj = nlg2
intpos = nlg1
isma = mul_tab(iml,imr)
goto(1,2,3,4,5,6,7,8,9,10,11,12),linelp
! linelp=2 a&r--d&l^r<-->a^l
2 continue
call ar_drl_ext_al_new(14,lri,lrj)
goto 100
! linelp=6 a&r--b&l--b^r<-->a^l
6 continue
call ar_bl_br_ext_al_new(14,intpos,isma,1)
goto 100
! linelp=11 d&rl--b^r<-->a^l
11 continue
call drl_br_ext_al_new(14,lri,lrj)   ! start
goto 100
1 continue
goto 100
3 continue
goto 100
4 continue
goto 100
5 continue
goto 100
7 continue
goto 100
8 continue
goto 100
9 continue
goto 100
10 continue
goto 100
12 continue
goto 100

100 continue
return

end subroutine dt_ext_head_in_act

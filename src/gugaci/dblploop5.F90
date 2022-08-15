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

subroutine ds_ardlr_act_c1(lin)
!=======================================================================
! ds(7-1) ar(23)-drl(30)-
!=======================================================================

use gugaci_global, only: ipae, ipael, jml, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, w0_ds
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmd, lrd, lri, mpl, ni
real(kind=wp) :: w0ds1
integer(kind=iwp), external :: iwalk_ad

do lri=norb_frz+1,norb_dz
  do lrd=norb_frz+1,lri-1
    lmd = lsm_inn(lrd)
    if (lmd /= jml) cycle
    iwdr = just(lri,lri)
    iwdl = jud(lrd)
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    w0ds1 = w0_ds(1)
    ni = mod(norb_dz-lrd,2)
    if (ni == 0) w0ds1 = -w0ds1
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0ds1
      vplp_w1(mpl) = Zero
    end do
    call ar_drl_ext_al_new(lin,lrd,lri)
  end do
end do

return

end subroutine ds_ardlr_act_c1

subroutine ds_arblbr_act_c1(lin)
!=======================================================================
! ds(7-3) ar(23)-bl(32)-br(31)-
!=======================================================================

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, lsm_inn, mhlp, mtype, ngw2, ngw3, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, &
                         vplpnew_w1, w0_ds, w1_ds
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, lmd, lmi, lmij, lmj, lrd, lri, lrj, mpl, ni
real(kind=wp) :: w0ds, w1ds
integer(kind=iwp), external :: iwalk_ad

isma = Mul(iml,imr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    if (lmij /= jmr) cycle
    ! ds(7-3) ar(23)-bl(32)-br(31)-
    do lrd=norb_frz+1,lri-1
      lmd = lsm_inn(lrd)
      if (lmd /= jml) cycle
      ijk = lrd-norb_frz+ngw2(lri-norb_frz)+ngw3(lrj-norb_frz)
      intpos = intind_ijka(ijk)
      w0ds = w0_ds(3)
      w1ds = w1_ds(3)
      ni = mod(norb_dz-lrj+lri-lrd,2)
      if (ni == 0) w0ds = -w0ds
      if (ni == 0) w1ds = -w1ds

      iwdr = just(lri,lrj)
      iwdl = jud(lrd)
      do mpl=1,mhlp
        iwal = lpnew_lwei(mpl)
        iwar = lpnew_rwei(mpl)
        lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
        lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
      end do
      do mpl=1,mtype
        vplp_w0(mpl) = vplpnew_w0(mpl)*w0ds
        vplp_w1(mpl) = vplpnew_w1(mpl)*w1ds
      end do
      call ar_bl_br_ext_al_new(lin,intpos,isma,1)
    end do
  end do
end do

return

end subroutine ds_arblbr_act_c1

subroutine dt_arblbr_act_c1(lin)
!=======================================================================
! dt(14) ar(23)-bl(32)-br(32)-
!=======================================================================

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, lsm_inn, mhlp, mtype, ngw2, ngw3, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, &
                         vplpnew_w1, w0_dt, w1_dt
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, lmd, lmi, lmij, lmj, lrd, lri, lrj, mpl, ni
real(kind=wp) :: w0dt, w1dt
integer(kind=iwp), external :: iwalk_ad

isma = Mul(iml,imr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    if (lmij /= jmr) cycle
    ! dt(14) ar(23)-bl(32)-br(32)-
    do lrd=norb_frz+1,lri-1
      lmd = lsm_inn(lrd)
      if (lmd /= jml) cycle
      ijk = lrd-norb_frz+ngw2(lri-norb_frz)+ngw3(lrj-norb_frz)
      intpos = intind_ijka(ijk)
      w0dt = w0_dt
      w1dt = w1_dt
      ni = mod(norb_dz-lrj+lri-lrd,2)
      if (ni == 0) w0dt = -w0dt
      if (ni == 0) w1dt = -w1dt

      iwdr = just(lri,lrj)
      iwdl = jud(lrd)
      do mpl=1,mhlp
        iwal = lpnew_lwei(mpl)
        iwar = lpnew_rwei(mpl)
        lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
        lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
      end do
      do mpl=1,mtype
        vplp_w0(mpl) = vplpnew_w0(mpl)*w0dt
        vplp_w1(mpl) = vplpnew_w1(mpl)*w1dt
      end do
      call ar_bl_br_ext_al_new(lin,intpos,isma,1)
    end do
  end do
end do

return

end subroutine dt_arblbr_act_c1

subroutine dv_ar_act_blbr(lin,jk)
!=======================================================================
! dv(23-1) ar(23)-
!=======================================================================

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jpad, jpadl, jud, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, &
                         lsm_inn, mhlp, mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_dv
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, jk
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, lmi, lri, mpl, ni
real(kind=wp) :: w0dv1
integer(kind=iwp), external :: iwalk_ad

isma = Mul(iml,imr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (lmi /= jml) cycle
  ijk = lri-norb_frz+jk
  intpos = intind_ijka(ijk)

  w0dv1 = w0_dv(1)
  ni = mod(norb_dz-lri,2)
  if (ni == 1) w0dv1 = -w0dv1
  iwdl = jud(lri)
  iwdr = 0
  do mpl=1,mtype
    vplp_w0(mpl) = vplpnew_w0(mpl)*w0dv1
    vplp_w1(mpl) = vplpnew_w1(mpl)*w0dv1
  end do
  do mpl=1,mhlp
    iwal = lpnew_lwei(mpl)
    iwar = lpnew_rwei(mpl)
    lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
    lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
  end do
  call ar_bl_br_ext_al_new(lin,intpos,isma,1)
end do

return

end subroutine dv_ar_act_blbr

!subroutine dv_ar_act_brbr(lin,lrai,lraj)
!!=======================================================================
!! dv(23-1) ar(23)-
!!=======================================================================
!
!use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jpad, jpadl, jud, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, &
!                         lsm_inn, mhlp, mtype, ngw2, ngw3, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_dv
!use Symmetry_Info, only: Mul
!use Definitions, only: wp, iwp
!
!implicit none
!integer(kind=iwp), intent(in) :: lin, lrai, lraj
!integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, lmi, lri, mpl, ni
!real(kind=wp) :: w0dv1
!integer(kind=iwp), external :: iwalk_ad
!
!isma = Mul(iml,imr)
!do lri=norb_frz+1,norb_dz
!  lmi = lsm_inn(lri)
!  if (lmi /= jml) cycle
!  ijk = lri-norb_frz+ngw2(lrai-norb_frz)+ngw3(lraj-norb_frz)
!  intpos = intind_ijka(ijk)
!
!  w0dv1 = w0_dv(1)
!  ni = mod(norb_dz-lri,2)
!  if (ni == 1) w0dv1 = -w0dv1
!  iwdl = jud(lri)
!  iwdr = 0
!  do mpl=1,mtype
!    vplp_w0(mpl) = vplpnew_w0(mpl)*w0dv1
!    vplp_w1(mpl) = vplpnew_w1(mpl)*w0dv1
!  end do
!  do mpl=1,mhlp
!    iwal = lpnew_lwei(mpl)
!    iwar = lpnew_rwei(mpl)
!    lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
!    lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
!  end do
!  call ar_br_br_ext_ar_new(lin,intpos,isma)
!end do
!
!return
!
!end subroutine dv_ar_act_brbr

subroutine dv_ar_act_dlr(lin,lra)
!=======================================================================
! dv(23-1) ar(23)-
!=======================================================================

use gugaci_global, only: ipae, ipael, jml, jpad, jpadl, jud, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_dv
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lri, mpl, ni
real(kind=wp) :: w0dv1
integer(kind=iwp), external :: iwalk_ad

do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (lmi /= jml) cycle

  w0dv1 = w0_dv(1)
  ni = mod(norb_dz-lri,2)
  if (ni == 1) w0dv1 = -w0dv1
  iwdl = jud(lri)
  iwdr = 0
  do mpl=1,mtype
    vplp_w0(mpl) = vplpnew_w0(mpl)*w0dv1
    vplp_w1(mpl) = vplpnew_w1(mpl)*w0dv1
  end do
  do mpl=1,mhlp
    iwal = lpnew_lwei(mpl)
    iwar = lpnew_rwei(mpl)
    lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
    lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
  end do
  call ar_drl_ext_al_new(lin,lri,lra)
end do

return

end subroutine dv_ar_act_dlr

subroutine tt_arbl_act_bl(lin,lra)
!=======================================================================
! tt(11-1) (22)ar(23)-bl(32)-
! tt(11-1) ar(23)-c'(22)-bl(32)-
! tt(11-1) ar(23)-bl(32)-c"(22)-
!=======================================================================

use gugaci_global, only: iml, imr, intind_ijka, jml, jmr, lsm_inn, ngw2, ngw3, norb_dz, norb_frz
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, jmlr, lmi, lmij, lmj, lri, lrj, nk

jmlr = Mul(jml,jmr)
isma = Mul(iml,imr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    if (lmij /= jmlr) cycle
    ijk = lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz)
    intpos = intind_ijka(ijk)
    call tt1_ext(lri,lrj,nk,1)
    call ar_bl_bl_ext_ar_new(lin,intpos,isma,nk)
    call tt1_ext(lri,lrj,nk,-1)
    call ar_bl_bl_ext_ar_new(lin,intpos,isma,nk)
  end do
end do

return

end subroutine tt_arbl_act_bl

subroutine tt_drl_act_bl(lin,lra)
! tt(11-2) (22)drl(22)-
! tt(11-2) drl(22)-c"(22)-
! tt(11-3) (22)(22)drl(33)-
! tt(11-3) (22)drl(33)-c"(22)-
! tt(11-3) drl(33)-c"(22)-c"(22)-

use gugaci_global, only: ipae, ipael, jml, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_tt, w1_tt
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lmij, lmj, lri, lrj, lrk, mpl
integer(kind=iwp), external :: iwalk_ad

do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    if (lmij /= jml) cycle
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0_tt(2)
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1_tt(2)
    end do
    ! tt(11-2) (22)drl(22)-
    ! tt(11-2) drl(22)-c"(22)-
    iwdl = just(lri,lrj)
    iwdr = iwdl
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    call drl_bl_ext_ar_new(lin,lri,lra)
    call drl_bl_ext_ar_new(lin,lrj,lra)
    ! tt(11-3) drl(33)-c"(22)-c"(22)-
    ! tt(11-3) (22)drl(33)-c"(22)-
    ! tt(11-3) (22)(22)drl(33)-
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0_tt(3)
      vplp_w1(mpl) = Zero
    end do
    if (lra > norb_dz) then
      call drl_bl_sum_ar_new(lin,lri,lrj,lra)
    else
      do lrk=1,norb_dz
        if (lrk == lri) cycle
        if (lrk == lrj) cycle
        call drl_bl_ext_ar_new(lin,lrk,lra)
      end do
    end if
    !do lrk=1,norb_dz
    !  if (lrk == lri) cycle
    !  if (lrk == lrj) cycle
    !  call drl_bl_ext_ar_new(lin,lrk,lra)
    !end do
  end do
end do

return

end subroutine tt_drl_act_bl

subroutine tt_arbl_act_br(lin,lra)
!=======================================================================
! tt(11-1) (22)ar(23)-bl(32)-
! tt(11-1) ar(23)-c'(22)-bl(32)-
! tt(11-1) ar(23)-bl(32)-c"(22)-
!=======================================================================

use gugaci_global, only: iml, imr, intind_ijka, jml, jmr, lsm_inn, ngw2, ngw3, norb_dz, norb_frz
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, jmlr, lmi, lmij, lmj, lri, lrj, nk

jmlr = Mul(jml,jmr)
isma = Mul(iml,imr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    if (lmij /= jmlr) cycle
    ijk = lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz)
    intpos = intind_ijka(ijk)
    call tt1_ext(lri,lrj,nk,1)
    call ar_bl_br_ext_al_new(lin,intpos,isma,nk)
    call tt1_ext(lri,lrj,nk,-1)
    call ar_bl_br_ext_al_new(lin,intpos,isma,nk)
  end do
end do

return

end subroutine tt_arbl_act_br

subroutine tt_drl_act_br(lin,lra)
! tt(11-2) (22)drl(22)-
! tt(11-2) drl(22)-c"(22)-
! tt(11-3) (22)(22)drl(33)-
! tt(11-3) (22)drl(33)-c"(22)-
! tt(11-3) drl(33)-c"(22)-c"(22)-

use gugaci_global, only: ipae, ipael, jml, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_tt, w1_tt
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lmij, lmj, lri, lrj, lrk, mpl
integer(kind=iwp), external :: iwalk_ad

do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    if (lmij /= jml) cycle
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0_tt(2)
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1_tt(2)
    end do
    ! tt(11-2) (22)drl(22)-
    ! tt(11-2) drl(22)-c"(22)-
    iwdl = just(lri,lrj)
    iwdr = iwdl
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    call drl_br_ext_al_new(lin,lri,lra)
    call drl_br_ext_al_new(lin,lrj,lra)
    ! tt(11-3) drl(33)-c"(22)-c"(22)-
    ! tt(11-3) (22)drl(33)-c"(22)-
    ! tt(11-3) (22)(22)drl(33)-
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0_tt(3)
      vplp_w1(mpl) = Zero
    end do
    if (lra > norb_dz) then
      call drl_br_sum_al_new(lin,lri,lrj,lra)
    else
      do lrk=1,norb_dz
        if (lrk == lri) cycle
        if (lrk == lrj) cycle
        call drl_br_ext_al_new(lin,lrk,lra)
      end do
    end if
    !do lrk=1,norb_dz
    !  if (lrk == lri) cycle
    !  if (lrk == lrj) cycle
    !  call drl_br_ext_al_new(lin,lrk,lra)
    !end do
  end do
end do

return

end subroutine tt_drl_act_br

subroutine st_drl_act_bl(lin,lra)
! st(2-5) (22)drl(12)-
! st(2-6) drl(22)-c"(12)-
! st(2-7) drl(12)-c"(22)-

use gugaci_global, only: ipae, ipael, jml, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w1, w1_st
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lmij, lmj, lri, lrj, mpl
integer(kind=iwp), external :: iwalk_ad

do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    if (lmij /= jml) cycle
    !-------------------------------------------------------------------
    ! st(2-5) (22)drl(12)-
    iwdl = just(lri,lrj)
    iwdr = iwdl
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    do mpl=1,mtype
      vplp_w0(mpl) = Zero
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1_st(5)
    end do
    call drl_bl_ext_ar_new(lin,lrj,lra)
    ! st(2-6) drl(22)-c"(12)-
    do mpl=1,mtype
      vplp_w0(mpl) = Zero
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1_st(6)
    end do
    call drl_bl_ext_ar_new(lin,lri,lra)
  end do
end do

return

end subroutine st_drl_act_bl

subroutine st_arbl_act_bl(lin,lra)

use gugaci_global, only: iml, imr, intind_ijka, ngw2, ngw3, norb_dz, norb_frz
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, lri, lrj, nk

isma = Mul(iml,imr)
do lri=norb_frz+1,norb_dz-1
  do lrj=lri+1,norb_dz
    ijk = lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz)
    intpos = intind_ijka(ijk)
    call st1_ext(lri,lrj,nk)
    if (nk /= 0) call ar_bl_bl_ext_ar_new(lin,intpos,isma,1)
    call st2_ext(lri,lrj,nk)
    if (nk /= 0) call ar_bl_bl_ext_ar_new(lin,intpos,isma,nk)
    call st4_ext(lri,lrj,nk,1)
    if (nk /= 0) call ar_bl_bl_ext_ar_new(lin,intpos,isma,nk)
    call st4_ext(lri,lrj,nk,-1)
    if (nk /= 0) call ar_bl_bl_ext_ar_new(lin,intpos,isma,nk)
  end do
end do

return

end subroutine st_arbl_act_bl

subroutine st_drl_act_br(lin,lra)
! st(2-5) (22)drl(12)-
! st(2-6) drl(22)-c"(12)-
! st(2-7) drl(12)-c"(22)-

use gugaci_global, only: ipae, ipael, jml, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w1, w1_st
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lmij, lmj, lri, lrj, mpl
integer(kind=iwp), external :: iwalk_ad

do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    if (lmij /= jml) cycle
    !-------------------------------------------------------------------
    ! st(2-5) (22)drl(12)-
    iwdl = just(lri,lrj)
    iwdr = iwdl
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    do mpl=1,mtype
      vplp_w0(mpl) = Zero
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1_st(5)
    end do
    call drl_br_ext_al_new(lin,lrj,lra)
    ! st(2-6) drl(22)-c"(12)-
    do mpl=1,mtype
      vplp_w0(mpl) = Zero
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1_st(6)
    end do
    call drl_br_ext_al_new(lin,lri,lra)
  end do
end do

return

end subroutine st_drl_act_br

subroutine st_arbl_act_br(lin,lra)

use gugaci_global, only: iml, imr, intind_ijka, ngw2, ngw3, norb_dz, norb_frz
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, lri, lrj, nk

isma = Mul(iml,imr)
do lri=norb_frz+1,norb_dz-1
  do lrj=lri+1,norb_dz
    ijk = lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz)
    intpos = intind_ijka(ijk)
    call st1_ext(lri,lrj,nk)
    if (nk /= 0) call ar_bl_br_ext_al_new(lin,intpos,isma,1)
    call st2_ext(lri,lrj,nk)
    if (nk /= 0) call ar_bl_br_ext_al_new(lin,intpos,isma,nk)
    call st4_ext(lri,lrj,nk,1)
    if (nk /= 0) call ar_bl_br_ext_al_new(lin,intpos,isma,nk)
    call st4_ext(lri,lrj,nk,-1)
    if (nk /= 0) call ar_bl_br_ext_al_new(lin,intpos,isma,nk)
  end do
end do

return

end subroutine st_arbl_act_br

subroutine ss_arbl_act_br(lin,lra)
!=======================================================================
! ss(1-3)  ar(13)-bl(20)-
! ss(1-6)  (11)-ar(23)-bl(32)-
! ss(1-7)  ar(13)-c'(21)-bl(32)-
! ss(1-8)  ar(13)-c'(22)-bl(31)-
! ss(1-9)  ar(23)-c'(11)-bl(32)-
! ss(1-11) ar(13)-bl(31)-c"(22)-
! ss(1-12) ar(13)-bl(32)-c"(21)-
! ss(1-13) ar(23)-bl(31)-c"(12)-
!=======================================================================

use gugaci_global, only: iml, imr, intind_ijka, ngw2, ngw3, norb_dz, norb_frz
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, lri, lrj, nk

isma = Mul(iml,imr)
do lri=norb_frz+1,norb_dz
  do lrj=lri+1,norb_dz
    ijk = lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz)
    intpos = intind_ijka(ijk)
    !-------------------------------------------------------------------
    call ss2_ext(lri,lrj,nk)
    if (nk /= 0) call ar_bl_br_ext_al_new(lin,intpos,isma,1)
    call ss4_ext(lri,lrj,nk)
    if (nk /= 0) call ar_bl_br_ext_al_new(lin,intpos,isma,1)
    call ss5_ext(lri,lrj,nk)
    if (nk /= 0) call ar_bl_br_ext_al_new(lin,intpos,isma,nk)
    call ss10_ext(lri,lrj,nk)
    if (nk /= 0) call ar_bl_br_ext_al_new(lin,intpos,isma,nk)
    call ss14_ext(lri,lrj,nk)
    if (nk /= 0) call ar_bl_br_ext_al_new(lin,intpos,isma,nk)
    !-------------------------------------------------------------------
  end do
end do

return

end subroutine ss_arbl_act_br

subroutine ss_drl_act_br(lin,lra)

use gugaci_global, only: ipae, ipael, jml, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_ss, w1_ss
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lmij, lmj, lri, lrj, lrk, mpl
integer(kind=iwp), external :: iwalk_ad

do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (jml == 1) then
    ! ss(1-20) drl(33)-c"(00)-
    iwdl = just(lri,lri)
    iwdr = iwdl
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0_ss(20)
      vplp_w1(mpl) = Zero
    end do
    if (lri > norb_dz) then
      call drl_br_sum_al_new(lin,lri,0,lri)
    else
      do lrk=1,norb_dz
        if (lrk == lri) cycle
        call drl_br_ext_al_new(lin,lrk,lra)
      end do
    end if
  end if
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    if (lmij /= jml) cycle
    ! ss(1-15) (22)-drl(11)-
    iwdl = just(lri,lrj)
    iwdr = iwdl
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0_ss(15)
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1_ss(15)
    end do
    call drl_br_ext_al_new(lin,lrj,lra)
    ! ss(1-17) drl(22)-c"(11)-
    iwdl = just(lri,lrj)
    iwdr = iwdl
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0_ss(17)
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1_ss(17)
    end do
    call drl_br_ext_al_new(lin,lri,lra)
    ! ss(1-20) drl(33)-c"(22)-c"(11)-
    ! ss(1-20) (22)drl(33)-c"(11)-
    ! ss(1-20) (22)(11)drl(33)-
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0_ss(20)
      vplp_w1(mpl) = Zero
    end do
    if (lra > norb_dz) then
      call drl_br_sum_al_new(lin,lri,lrj,lra)
    else
      do lrk=1,norb_dz
        if (lrk == lri) cycle
        if (lrk == lrj) cycle
        call drl_br_ext_al_new(lin,lrk,lra)
      end do
    end if
  end do
end do

return

end subroutine ss_drl_act_br

subroutine ts_arbl_act_br(lin,lra)
!=======================================================================
! ts(3) a&r-b^l-  act -b&l ............................................

use gugaci_global, only: iml, imr, intind_ijka, ngw2, ngw3, norb_dz, norb_frz
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, lri, lrj, nk

isma = Mul(iml,imr)
do lri=norb_frz+1,norb_dz
  do lrj=lri+1,norb_dz
    ijk = lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz)
    intpos = intind_ijka(ijk)
    call ts1_ext(lri,lrj,nk)
    if (nk /= 0) call ar_bl_br_ext_al_new(lin,intpos,isma,1)
    call ts2_ext(lri,lrj,nk,1)
    if (nk /= 0) call ar_bl_br_ext_al_new(lin,intpos,isma,nk)
    call ts2_ext(lri,lrj,nk,-1)
    if (nk /= 0) call ar_bl_br_ext_al_new(lin,intpos,isma,nk)
    call ts4_ext(lri,lrj,nk)
    if (nk /= 0) call ar_bl_br_ext_al_new(lin,intpos,isma,nk)
  end do
end do

return

end subroutine ts_arbl_act_br

subroutine ts_arbl_act_bl(lin,lra)
!=======================================================================
! ts(3) a&r-b^l-  act -b&l ............................................

use gugaci_global, only: iml, imr, intind_ijka, ngw2, ngw3, norb_dz, norb_frz
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, lri, lrj, nk

isma = Mul(iml,imr)
do lri=norb_frz+1,norb_dz
  do lrj=lri+1,norb_dz
    ijk = lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz)
    intpos = intind_ijka(ijk)
    call ts1_ext(lri,lrj,nk)
    if (nk /= 0) call ar_bl_bl_ext_ar_new(lin,intpos,isma,1)
    call ts2_ext(lri,lrj,nk,1)
    if (nk /= 0) call ar_bl_bl_ext_ar_new(lin,intpos,isma,nk)
    call ts2_ext(lri,lrj,nk,-1)
    if (nk /= 0) call ar_bl_bl_ext_ar_new(lin,intpos,isma,nk)
    call ts4_ext(lri,lrj,nk)
    if (nk /= 0) call ar_bl_bl_ext_ar_new(lin,intpos,isma,nk)
  end do
end do

return

end subroutine ts_arbl_act_bl

subroutine ss_drl_act_bl(lin,lra)
!=======================================================================
! ss(1)    act -bl-
! ss(1-16) (11)-drl(22)-
! ss(1-17) drl(22)-c"(11)-
! ss(1-18) drl(11)-c"(22)-
! ss(1-19) drl(12)-c"(21)-
! ss(1-20) drl(33)-c"(11)-c"(22)-
! ss(1-20) (11)drl(33)-c"(22)-
! ss(1-20) (11)(22)drl(33)-
!=======================================================================

use gugaci_global, only: ipae, ipael, jml, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_ss, w1_ss
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lmi, lmij, lmj, lri, lrj, lrk, mpl
integer(kind=iwp), external :: iwalk_ad

do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (jml == 1) then
    ! ss(1-20) drl(33)-c"(00)-
    iwdl = just(lri,lri)
    iwdr = iwdl
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0_ss(20)
      vplp_w1(mpl) = Zero
    end do
    if (lra > norb_dz) then
      call drl_bl_sum_ar_new(lin,lri,0,lra)
    else
      do lrk=1,norb_dz
        if (lrk == lri) cycle
        call drl_bl_ext_ar_new(lin,lrk,lra)
      end do
    end if
    !do lrk=1,norb_dz
    !  if (lrk == lri) cycle
    !  call drl_bl_ext_ar_new(lin,lrk,lra)
    !end do
  end if
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    if (lmij /= jml) cycle
    iwdl = just(lri,lrj)
    iwdr = iwdl
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    ! ss(1-15) (22)-drl(11)-
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0_ss(15)
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1_ss(15)
    end do
    call drl_bl_ext_ar_new(lin,lrj,lra)
    ! ss(1-17) drl(22)-c"(11)-
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0_ss(17)
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1_ss(17)
    end do
    call drl_bl_ext_ar_new(lin,lri,lra)
    ! ss(1-20) drl(33)-c"(22)-c"(11)-
    ! ss(1-20) (22)drl(33)-c"(11)-
    ! ss(1-20) (22)(11)-drl(33)
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0_ss(20)
      vplp_w1(mpl) = Zero
    end do
    if (lra > norb_dz) then
      call drl_bl_sum_ar_new(lin,lri,lrj,lra)
    else
      do lrk=1,norb_dz
        if (lrk == lri) cycle
        if (lrk == lrj) cycle
        call drl_bl_ext_ar_new(lin,lrk,lra)
      end do
    end if
    !do lrk=1,norb_dz
    !  if (lrk == lri) cycle
    !  if (lrk == lrj) cycle
    !  call drl_bl_ext_ar_new(lin,lrk,lra)
    !end do
  end do
end do

return

end subroutine ss_drl_act_bl

subroutine ss_arbl_act_bl(lin,lra)
!=======================================================================
! ss(1)    act -bl-
! ss(1-3)  ar(13)-bl(20)-
! ss(1-6)  (11)-ar(23)-bl(32)-
! ss(1-7)  ar(13)-c'(21)-bl(32)-
! ss(1-8)  ar(13)-c'(22)-bl(31)-
! ss(1-9)  ar(23)-c'(11)-bl(32)-
! ss(1-11) ar(13)-bl(31)-c"(22)-
! ss(1-12) ar(13)-bl(32)-c"(21)-
! ss(1-13) ar(23)-bl(31)-c"(12)-
!=======================================================================

use gugaci_global, only: iml, imr, intind_ijka, ngw2, ngw3, norb_dz, norb_frz
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, lri, lrj, nk

isma = Mul(iml,imr)
do lri=norb_frz+1,norb_dz
  do lrj=lri+1,norb_dz
    ijk = lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz)
    intpos = intind_ijka(ijk)
    call ss2_ext(lri,lrj,nk)
    if (nk /= 0) call ar_bl_bl_ext_ar_new(lin,intpos,isma,1)
    call ss4_ext(lri,lrj,nk)
    if (nk /= 0) call ar_bl_bl_ext_ar_new(lin,intpos,isma,1)
    call ss5_ext(lri,lrj,nk)
    if (nk /= 0) call ar_bl_bl_ext_ar_new(lin,intpos,isma,nk)
    call ss10_ext(lri,lrj,nk)
    if (nk /= 0) call ar_bl_bl_ext_ar_new(lin,intpos,isma,nk)
    call ss14_ext(lri,lrj,nk)
    if (nk /= 0) call ar_bl_bl_ext_ar_new(lin,intpos,isma,nk)
  end do
end do

return

end subroutine ss_arbl_act_bl

subroutine sd_ar_act_bl(lin,lra)
! sd(6-1) a&r(02)-
! sd(6-2) (22)a&(13)-
! sd(6-3) a&r(13)c'(22)-
! sd(6-4) a&r(23)c'(12)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, &
                         mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_sd
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, jmlr, lmi, lmk, lri, lrk, mpl, ni
real(kind=wp) :: w0sd1, w0sd2, w0sd4
integer(kind=iwp), external :: iwalk_ad

jmlr = Mul(jml,jmr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (lmi /= jmlr) cycle
  w0sd1 = w0_sd(1)
  w0sd2 = w0_sd(2)
  w0sd4 = w0_sd(4)
  ni = mod(norb_dz-lri,2)
  if (ni == 1) then
    w0sd1 = -w0sd1
    w0sd2 = -w0sd2
    w0sd4 = -w0sd4
  end if
  if ((jml == 1) .and. (lmi == jmr)) then
    ! sd(6-1) a&r(02)-
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
  ! sd(6-2) (22)a&(13)-
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
    if (lin == 1) call ar_bl_ext_ss(lri,lra,1)
    if (lin == 2) call ar_bl_ext_st(lri,lra,1)
    if (lin == 3) call ar_bl_ext_ts(lri,lra,1)
    if (lin == 11) call ar_bl_ext_tt(lri,lra,1)
  end do
  ! sd(6-4) a&r(23)-c'(12)-
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

end subroutine sd_ar_act_bl

subroutine sd_ar_act_blbl(lin,jk)
! sd(6-1) a&r(02)-
! sd(6-2) (22)a&(13)-
! sd(6-3) a&r(13)c'(22)-
! sd(6-4) a&r(23)c'(12)-

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, lsm_inn, mhlp, mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_sd
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, jk
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, jmlr, lmi, lmk, lri, lrk, mpl, ni
real(kind=wp) :: w0sd1, w0sd2, w0sd4
integer(kind=iwp), external :: iwalk_ad

jmlr = Mul(jml,jmr)
isma = Mul(iml,imr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (lmi /= jmlr) cycle
  ijk = lri-norb_frz+jk
  intpos = intind_ijka(ijk)
  w0sd1 = w0_sd(1)
  w0sd2 = w0_sd(2)
  w0sd4 = w0_sd(4)
  ni = mod(norb_dz-lri,2)
  if (ni == 1) then
    w0sd1 = -w0sd1
    w0sd2 = -w0sd2
    w0sd4 = -w0sd4
  end if
  ! sd(6-1) a&r(02)-
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
    call ar_bl_bl_ext_ar_new(lin,intpos,isma,1)
  end if
  !---------------------------------------------------------------------
  ! sd(6-2) c(22)a&(13)-
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
    call ar_bl_bl_ext_ar_new(lin,intpos,isma,1)
  end do
  !---------------------------------------------------------------------
  ! sd(6-4) a&r(23)c'(12)-
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
    call ar_bl_bl_ext_ar_new(lin,intpos,isma,1)
  end do
end do

return

end subroutine sd_ar_act_blbl

subroutine sd_ar_act_brbr(lin,jk)
! sd(6-1) a&r(02)-
! sd(6-2) (22)a&(13)-
! sd(6-3) a&r(13)c'(22)-
! sd(6-4) a&r(23)c'(12)-

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, lsm_inn, mhlp, mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_sd
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, jk
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, jmlr, lmi, lmk, lri, lrk, mpl, ni
real(kind=wp) :: w0sd1, w0sd2, w0sd4
integer(kind=iwp), external :: iwalk_ad

jmlr = Mul(jml,jmr)
isma = Mul(iml,imr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (lmi /= jmlr) cycle
  ijk = lri-norb_frz+jk
  intpos = intind_ijka(ijk)
  w0sd1 = w0_sd(1)
  w0sd2 = w0_sd(2)
  w0sd4 = w0_sd(4)
  ni = mod(norb_dz-lri,2)
  if (ni == 1) then
    w0sd1 = -w0sd1
    w0sd2 = -w0sd2
    w0sd4 = -w0sd4
  end if
  ! sd(6-1) a&r(02)-
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
    call ar_br_br_ext_ar_new(lin,intpos,isma)
  end if
  !---------------------------------------------------------------------
  ! sd(6-2) c(22)a&(13)-
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
    call ar_br_br_ext_ar_new(lin,intpos,isma)
  end do
  !---------------------------------------------------------------------
  ! sd(6-4) a&r(23)c'(12)-
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
    call ar_br_br_ext_ar_new(lin,intpos,isma)
  end do
end do

return

end subroutine sd_ar_act_brbr

subroutine sd_ar_act_blbr(lin,jk)
! sd(6-1) a&r(02)-
! sd(6-2) (22)a&(13)-
! sd(6-3) a&r(13)c'(22)-
! sd(6-4) a&r(23)c'(12)-

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, lsm_inn, mhlp, mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_sd
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, jk
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, jmlr, lmi, lmk, lri, lrk, mpl, ni
real(kind=wp) :: w0sd1, w0sd2, w0sd4
integer(kind=iwp), external :: iwalk_ad

jmlr = Mul(jml,jmr)
isma = Mul(iml,imr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (lmi /= jmlr) cycle
  ijk = lri-norb_frz+jk
  intpos = intind_ijka(ijk)
  w0sd1 = w0_sd(1)
  w0sd2 = w0_sd(2)
  w0sd4 = w0_sd(4)
  ni = mod(norb_dz-lri,2)
  if (ni == 1) then
    w0sd1 = -w0sd1
    w0sd2 = -w0sd2
    w0sd4 = -w0sd4
  end if
  ! sd(6-1) a&r(02)-
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
  ! sd(6-2) c(22)a&(13)-
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
    call ar_bl_br_ext_al_new(lin,intpos,isma,1)
  end do
  !---------------------------------------------------------------------
  ! sd(6-4) a&r(23)c'(12)-
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
    call ar_bl_br_ext_al_new(lin,intpos,isma,1)
  end do
end do

return

end subroutine sd_ar_act_blbr

subroutine sd_ar_act_dlr(lin,lra)
! sd(6-1) a&r(02)-
! sd(6-2) (22)a&(13)-
! sd(6-3) a&r(13)c'(22)-
! sd(6-4) a&r(23)c'(12)-

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, &
                         mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_sd
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, jmlr, lmi, lmk, lri, lrk, mpl, ni
real(kind=wp) :: w0sd1, w0sd2, w0sd4
integer(kind=iwp), external :: iwalk_ad

jmlr = Mul(jml,jmr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (lmi /= jmlr) cycle
  w0sd1 = w0_sd(1)
  w0sd2 = w0_sd(2)
  w0sd4 = w0_sd(4)
  ni = mod(norb_dz-lri,2)
  if (ni == 1) then
    w0sd1 = -w0sd1
    w0sd2 = -w0sd2
    w0sd4 = -w0sd4
  end if
  ! sd(6-1) a&r(02)-
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
    call ar_drl_ext_al_new(lin,lri,lra)
  end if
  !---------------------------------------------------------------------
  ! sd(6-2) c(22)a&(13)-
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
    call ar_drl_ext_al_new(lin,lri,lra)
  end do
  !---------------------------------------------------------------------
  ! sd(6-4) a&r(23)c'(12)-
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
    call ar_drl_ext_al_new(lin,lri,lra)
  end do
end do

return

end subroutine sd_ar_act_dlr

subroutine td_ar_act_bl(lin,lra)
! td(13-1) (22)a(23)
! td(13-1) a(23)c'(22)

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, &
                         mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_td
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
  w0td1 = w0_td(1)
  ni = mod(norb_dz-lri,2)
  if (ni == 1) w0td1 = -w0td1
  !---------------------------------------------------------------------
  ! td(13-1) (22)a(23)
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
  ! td(13-1) a(23)c'(22)
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

end subroutine td_ar_act_bl

subroutine td_ar_act_brbr(lin,jk)
! td(13-1) (22)a(23)
! td(13-1) a(23)c'(22)

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, lsm_inn, mhlp, mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_td
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, jk
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, jmlr, lmi, lmk, lri, lrk, mpl, ni
real(kind=wp) :: w0td1
integer(kind=iwp), external :: iwalk_ad

jmlr = Mul(jml,jmr)
isma = Mul(iml,imr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (lmi /= jmlr) cycle
  w0td1 = w0_td(1)
  ni = mod(norb_dz-lri,2)
  if (ni == 1) w0td1 = -w0td1
  ijk = lri-norb_frz+jk
  intpos = intind_ijka(ijk)
  !---------------------------------------------------------------------
  ! td(13-1) (22)a(23)
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
    call ar_br_br_ext_ar_new(lin,intpos,isma)
  end do
  !---------------------------------------------------------------------
  ! td(13-1) a(23)c'(22)
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
    call ar_br_br_ext_ar_new(lin,intpos,isma)
  end do
end do

return

end subroutine td_ar_act_brbr

subroutine td_ar_act_blbl(lin,jk)
! td(13-1) (22)a(23)
! td(13-1) a(23)c'(22)

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, lsm_inn, mhlp, mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_td
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, jk
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, jmlr, lmi, lmk, lri, lrk, mpl, ni
real(kind=wp) :: w0td1
integer(kind=iwp), external :: iwalk_ad

jmlr = Mul(jml,jmr)
isma = Mul(iml,imr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (lmi /= jmlr) cycle
  w0td1 = w0_td(1)
  ni = mod(norb_dz-lri,2)
  if (ni == 1) w0td1 = -w0td1
  ijk = lri-norb_frz+jk
  intpos = intind_ijka(ijk)
  !---------------------------------------------------------------------
  ! td(13-1) (22)a(23)
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
    call ar_bl_bl_ext_ar_new(lin,intpos,isma,1)
  end do
  !---------------------------------------------------------------------
  ! td(13-1) a(23)c'(22)
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
    call ar_bl_bl_ext_ar_new(lin,intpos,isma,1)
  end do
end do

return

end subroutine td_ar_act_blbl

subroutine td_ar_act_blbr(lin,jk)
! td(13-1) (22)a(23)
! td(13-1) a(23)c'(22)

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, &
                         lpnew_rwei, lsm_inn, mhlp, mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_td
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, jk
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, jmlr, lmi, lmk, lri, lrk, mpl, ni
real(kind=wp) :: w0td1
integer(kind=iwp), external :: iwalk_ad

jmlr = Mul(jml,jmr)
isma = Mul(iml,imr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (lmi /= jmlr) cycle
  w0td1 = w0_td(1)
  ni = mod(norb_dz-lri,2)
  if (ni == 1) w0td1 = -w0td1
  ijk = lri-norb_frz+jk
  intpos = intind_ijka(ijk)
  !---------------------------------------------------------------------
  ! td(13-1) (22)a(23)
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
    call ar_bl_br_ext_al_new(lin,intpos,isma,1)
  end do
  !---------------------------------------------------------------------
  ! td(13-1) a(23)c'(22)
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
    call ar_bl_br_ext_al_new(lin,intpos,isma,1)
  end do
end do

return

end subroutine td_ar_act_blbr

subroutine td_ar_act_dlr(lin,lra)
! td(13-1) (22)a(23)
! td(13-1) a(23)c'(22)

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, jud, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, &
                         mtype, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_td
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, jmlr, lmi, lmk, lri, lrk, mpl, ni
real(kind=wp) :: w0td1
integer(kind=iwp), external :: iwalk_ad

jmlr = Mul(jml,jmr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (lmi /= jmlr) cycle
  w0td1 = w0_td(1)
  ni = mod(norb_dz-lri,2)
  if (ni == 1) w0td1 = -w0td1
  !---------------------------------------------------------------------
  ! td(13-1) (22)a(23)
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
    call ar_drl_ext_al_new(lin,lri,lra)
  end do
  !---------------------------------------------------------------------
  ! td(13-1) a(23)c'(22)
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
    call ar_drl_ext_al_new(lin,lri,lra)
  end do
end do

return

end subroutine td_ar_act_dlr

subroutine dd_drl_act_bl(lin,lra)
! dd(19-2) drl(22)-
! dd(19-3) drl(33)-c"(22)-
! dd(19-3) (22)drl(33)-

use gugaci_global, only: ipae, ipael, jml, jpad, jpadl, jud, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_dd, w1_dd
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lml, lrk, lrl, mpl
real(kind=wp) :: w0, w1
integer(kind=iwp), external :: iwalk_ad

do lrl=norb_frz+1,norb_dz
  lml = lsm_inn(lrl)
  if (lml /= jml) cycle
  !---------------------------------------------------------------------
  iwdl = jud(lrl)
  iwdr = iwdl
  do mpl=1,mhlp
    iwal = lpnew_lwei(mpl)
    iwar = lpnew_rwei(mpl)
    lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
    lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
  end do
  w0 = w0_dd(2)
  w1 = w1_dd(2)
  do mpl=1,mtype
    vplp_w0(mpl) = vplpnew_w0(mpl)*w0
    vplp_w1(mpl) = vplpnew_w1(mpl)*w1
  end do
  call drl_bl_ext_ar_new(lin,lrl,lra)

  w0 = w0_dd(3)
  w1 = Zero
  do mpl=1,mtype
    vplp_w0(mpl) = vplpnew_w0(mpl)*w0
    vplp_w1(mpl) = vplpnew_w1(mpl)*w1
  end do
  if (lra > norb_dz) then
    call drl_bl_sum_ar_new(lin,lrl,0,lra)
  else
    do lrk=1,norb_dz
      if (lrk == lrl) cycle
      call drl_bl_ext_ar_new(lin,lrk,lra)
    end do
  end if
  !do lrk=1,norb_dz
  !  if (lrk == lrl) cycle
  !  call drl_bl_ext_ar_new(lin,lrk,lra)
  !end do
end do

return

end subroutine dd_drl_act_bl

subroutine dd_drl_act_br(lin,lra)
! dd(19-2) drl(22)-
! dd(19-3) drl(33)-c"(22)-
! dd(19-3) (22)drl(33)-

use gugaci_global, only: ipae, ipael, jml, jpad, jpadl, jud, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_dd, w1_dd
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, lml, lrk, lrl, mpl
real(kind=wp) :: w0, w1
integer(kind=iwp), external :: iwalk_ad

do lrl=norb_frz+1,norb_dz
  lml = lsm_inn(lrl)
  if (lml /= jml) cycle
  !---------------------------------------------------------------------
  iwdl = jud(lrl)
  iwdr = iwdl
  do mpl=1,mhlp
    iwal = lpnew_lwei(mpl)
    iwar = lpnew_rwei(mpl)
    lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
    lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
  end do
  w0 = w0_dd(2)
  w1 = w1_dd(2)
  do mpl=1,mtype
    vplp_w0(mpl) = vplpnew_w0(mpl)*w0
    vplp_w1(mpl) = vplpnew_w1(mpl)*w1
  end do
  call drl_br_ext_al_new(lin,lrl,lra)

  w0 = w0_dd(3)
  w1 = Zero
  do mpl=1,mtype
    vplp_w0(mpl) = vplpnew_w0(mpl)*w0
    vplp_w1(mpl) = vplpnew_w1(mpl)*w1
  end do
  do lrk=1,norb_dz
    if (lrk == lrl) cycle
    call drl_br_ext_al_new(lin,lrk,lra)
  end do
end do

return

end subroutine dd_drl_act_br

subroutine dd_arbl_act_bl(lin,lra)
! dd(19-1) a&r(23)-b&l(32)-

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, jud, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, &
                         lsm_inn, mhlp, mtype, ngw2, ngw3, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_dd, &
                         w1_dd
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, lml, lmr, lrl, lrr, mpl, ni
real(kind=wp) :: w0dd1, w1dd1
integer(kind=iwp), external :: iwalk_ad

isma = Mul(iml,imr)
do lrl=norb_frz+1,norb_dz-1
  lml = lsm_inn(lrl)
  if (lml /= jml) cycle
  do lrr=lrl+1,norb_dz
    lmr = lsm_inn(lrr)
    if (lmr /= jmr) cycle
    ijk = lrl-norb_frz+ngw2(lrr-norb_frz)+ngw3(lra-norb_frz)
    intpos = intind_ijka(ijk)
    w0dd1 = w0_dd(1)
    w1dd1 = w1_dd(1)
    ni = mod(lrr-lrl,2)
    if (ni == 0) w0dd1 = -w0_dd(1)
    if (ni == 0) w1dd1 = -w1_dd(1)
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0dd1
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1dd1
    end do

    iwdl = jud(lrl)
    iwdr = jud(lrr)
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    call ar_bl_bl_ext_ar_new(lin,intpos,isma,1)
    !-------------------------------------------------------------------
  end do
end do

return

end subroutine dd_arbl_act_bl

subroutine dd_arbl_act_br(lin,lra)
! dd(19-1) a&r(23)-b&l(32)-

use gugaci_global, only: iml, imr, intind_ijka, ipae, ipael, jml, jmr, jpad, jpadl, jud, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, &
                         lsm_inn, mhlp, mtype, ngw2, ngw3, norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_dd, w1_dd
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: ijk, intpos, isma, iwal, iwar, iwdl, iwdr, lml, lmr, lrl, lrr, mpl, ni
real(kind=wp) :: w0dd1, w1dd1
integer(kind=iwp), external :: iwalk_ad

isma = Mul(iml,imr)
do lrl=norb_frz+1,norb_dz-1
  lml = lsm_inn(lrl)
  if (lml /= jml) cycle
  do lrr=lrl+1,norb_dz
    lmr = lsm_inn(lrr)
    if (lmr /= jmr) cycle
    ijk = lrl-norb_frz+ngw2(lrr-norb_frz)+ngw3(lra-norb_frz)
    intpos = intind_ijka(ijk)
    w0dd1 = w0_dd(1)
    w1dd1 = w1_dd(1)
    ni = mod(lrr-lrl,2)
    if (ni == 0) w0dd1 = -w0_dd(1)
    if (ni == 0) w1dd1 = -w1_dd(1)
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0dd1
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1dd1
    end do

    iwdl = jud(lrl)
    iwdr = jud(lrr)
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    call ar_bl_br_ext_al_new(lin,intpos,isma,1)
    !-------------------------------------------------------------------
  end do
end do

return

end subroutine dd_arbl_act_br

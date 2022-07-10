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

subroutine drl_bl_ext_ar_new(lin,lrk,lri)

use gugaci_global, only: ihy, ihyl, ilsegdownwei, ipae, ipael, irsegdownwei, iseg_downwei, jpad, jpad_upwei, jpadl, jphy, jphyl, &
                         log_prod, logic_dh, logic_grad, lp_lwei, lp_rwei, lpnew_head, lpnew_lwei, lpnew_rwei, lsm_inn, mtype, &
                         ndim, nstaval, nvalue, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_sdplp, w1_sdplp
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lrk, lri
integer(kind=iwp) :: ihypos, ihyposl, ihyposr, ilpend, ilpsta, ilw, in_, iplp, irw, isma, iw0, iwal, iwal0, iwar, iwar0, iwd, &
                     iwuplwei, lphead
integer(kind=iwp), external :: iwalk_ad

isma = lsm_inn(lri)
iwuplwei = jpad_upwei(jpadl)
ilsegdownwei = iseg_downwei(ipael)
irsegdownwei = iseg_downwei(ipae)
do iw0=1,mtype
  w0_sdplp = vplpnew_w0(iw0)
  w1_sdplp = vplpnew_w1(iw0)
  if (logic_dh) w0_sdplp = vplp_w0(iw0)
  if (logic_dh) w1_sdplp = vplp_w1(iw0)
  if (logic_grad) then
    call lp9_drlbl_ext_calcuvalue_g(lri,lrk,isma)
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend
      if (logic_dh) then
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        if (lin == 6) call gsd_sequence_extspace_g(ilw,irw)
        if (lin == 13) call gtd_sequence_extspace_g(ilw,irw)
        if (lin == 23) call gdv_sequence_extspace_g(ilw,irw)
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
            if (lin == 6) call gsd_sequence_extspace_g(ilw,irw)
            if (lin == 13) call gtd_sequence_extspace_g(ilw,irw)
            if (lin == 23) call gdv_sequence_extspace_g(ilw,irw)
          end do
        end do
      end if
    end do
  else
    call lp9_drlbl_ext_calcuvalue_wyb(lri,lrk,isma)
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend
      if (logic_dh) then
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        if (lin == 6) call gsd_sequence_extspace(ilw,irw)
        if (lin == 13) call gtd_sequence_extspace(ilw,irw)
        if (lin == 23) call gdv_sequence_extspace(ilw,irw)
      else
        if (log_prod == 3) then
          lphead = lpnew_head(iplp)
          ihyposl = jphyl(lphead)
          ihyposr = jphy(lphead)
          ndim = ihyl(ihyposl)
        else
          ihyposl = jphy(iplp)
          ihyposr = ihyposl
          ndim = ihy(ihyposl)
        end if
        iwal0 = lpnew_lwei(iplp)
        iwar0 = lpnew_rwei(iplp)
        do in_=1,ndim
          iwal = iwal0+ihyl(ihyposl+in_)
          iwar = iwar0+ihy(ihyposr+in_)
          do iwd=0,iwuplwei-1
            ilw = iwalk_ad(jpadl,ipael,iwal,iwd)
            irw = iwalk_ad(jpad,ipae,iwar,iwd)
            if (lin == 6) call gsd_sequence_extspace(ilw,irw)
            if (lin == 13) call gtd_sequence_extspace(ilw,irw)
            if (lin == 23) call gdv_sequence_extspace(ilw,irw)
          end do
        end do
      end if
    end do

  end if
end do

return

end subroutine drl_bl_ext_ar_new

subroutine drl_br_ext_al_new(lin,lrk,lri)   !to be revised

use gugaci_global, only: ihy, ihyl, ilsegdownwei, ipae, ipael, irsegdownwei, iseg_downwei, jpad, jpad_upwei, jpadl, jphy, jphyl, &
                         log_prod, logic_dh, logic_grad, lp_lwei, lp_rwei, lpnew_head, lpnew_lwei, lpnew_rwei, lsm_inn, mtype, &
                         ndim, nstaval, nvalue, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_sdplp
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lrk, lri
integer(kind=iwp) :: ihypos, ihyposl, ihyposr, ilpend, ilpsta, ilw, in_, iplp, irw, isma, iw0, iwal, iwal0, iwar, iwar0, iwd, &
                     iwuplwei, lphead, nlp_value
real(kind=wp), parameter :: crl = 1.0e-8_wp
integer(kind=iwp), external :: iwalk_ad

iwuplwei = jpad_upwei(jpad)
ilsegdownwei = iseg_downwei(ipae)
irsegdownwei = iseg_downwei(ipael)
isma = lsm_inn(lri)

do iw0=1,mtype
  w0_sdplp = vplpnew_w0(iw0)-vplpnew_w1(iw0)
  if (logic_dh) w0_sdplp = vplp_w0(iw0)-vplp_w1(iw0)
  if (abs(w0_sdplp) < crl) cycle
  if (logic_grad) then
    call lp678_ext_calcuvalue_g(lri,lrk,isma,nlp_value)
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend
      if (logic_dh) then
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        if (lin == 7) call gsd_sequence_extspace_g(irw,ilw)
        if (lin == 14) call gtd_sequence_extspace_g(irw,ilw)
        if (lin == 26) call gdv_sequence_extspace_g(irw,ilw)
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
            if (lin == 7) call gsd_sequence_extspace_g(irw,ilw)
            if (lin == 14) call gtd_sequence_extspace_g(irw,ilw)
            if (lin == 26) call gdv_sequence_extspace_g(irw,ilw)
          end do
        end do
      end if
    end do

  else
    call lp678_ext_wyb_calcuvalue(lri,lrk,isma,nlp_value)
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend
      if (logic_dh) then
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        if (lin == 7) call gsd_sequence_extspace(irw,ilw)
        if (lin == 14) call gtd_sequence_extspace(irw,ilw)
        if (lin == 26) call gdv_sequence_extspace(irw,ilw)
      else
        if (log_prod == 3) then
          lphead = lpnew_head(iplp)
          ihyposl = jphyl(lphead)
          ihyposr = jphy(lphead)
          ndim = ihyl(ihyposl)
        else
          ihyposl = jphy(iplp)
          ihyposr = ihyposl
          ndim = ihy(ihyposl)
        end if
        iwal0 = lpnew_lwei(iplp)
        iwar0 = lpnew_rwei(iplp)
        do in_=1,ndim
          iwal = iwal0+ihyl(ihyposl+in_)
          iwar = iwar0+ihy(ihyposr+in_)
          do iwd=0,iwuplwei-1
            ilw = iwalk_ad(jpadl,ipael,iwal,iwd)
            irw = iwalk_ad(jpad,ipae,iwar,iwd)
            if (lin == 7) call gsd_sequence_extspace(irw,ilw)
            if (lin == 14) call gtd_sequence_extspace(irw,ilw)
            if (lin == 26) call gdv_sequence_extspace(irw,ilw)
          end do
        end do
      end if
    end do
  end if
end do

return

end subroutine drl_br_ext_al_new

subroutine ar_bl_br_ext_al_new(lin,intentry,isma,nk)

use gugaci_global, only: ihy, ihyl, ilsegdownwei, ipae, ipael, irsegdownwei, iseg_downwei, jpad, jpad_upwei, jpadl, jphy, jphyl, &
                         log_prod, logic_dh, logic_grad, lp_lwei, lp_rwei, lpnew_head, lpnew_lwei, lpnew_rwei, mtype, ndim, &
                         nstaval, nvalue, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_sdplp, w1_sdplp
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, intentry, isma, nk
integer(kind=iwp) :: ihypos, ihyposl, ihyposr, ilpend, ilpsta, ilw, in_, iplp, irw, iw0, iwal, iwal0, iwar, iwar0, iwd, iwuplwei, &
                     lphead, nlp_value
integer(kind=iwp), external :: iwalk_ad

!if ((lin == 14) .and. (intentry == 15195) .and. (isma == 1) .and. (jpad == 20) .and. (jpad == jpadl)) write(u6,*) 'bbs_tmp'

iwuplwei = jpad_upwei(jpad)
ilsegdownwei = iseg_downwei(ipae)
irsegdownwei = iseg_downwei(ipael)
do iw0=1,mtype
  w0_sdplp = vplpnew_w0(iw0)
  w1_sdplp = vplpnew_w1(iw0)
  if (logic_dh) w0_sdplp = vplp_w0(iw0)
  if (logic_dh) w1_sdplp = vplp_w1(iw0)

  if (logic_grad) then
    call lp11_arblbr_ext_calcuvalue_g(intentry,isma,nlp_value)
    ilpsta = nstaval(iw0)*nk+1
    ilpend = (nstaval(iw0)+nvalue(iw0))*nk
    do iplp=ilpsta,ilpend
      if (logic_dh) then                     !lp_head is in dbl_spa
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        if (lin == 7) call gsd_sequence_extspace_g(irw,ilw)
        if (lin == 14) call gtd_sequence_extspace_g(irw,ilw)
        if (lin == 26) call gdv_sequence_extspace_g(irw,ilw)
      else                                   !lp_head is in act_s
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
            if (lin == 7) call gsd_sequence_extspace_g(irw,ilw)
            if (lin == 14) call gtd_sequence_extspace_g(irw,ilw)
            if (lin == 26) call gdv_sequence_extspace_g(irw,ilw)
          end do
        end do
      end if
    end do
  else
    call lp11_arblbr_ext_calcuvalue(intentry,isma,nlp_value)
    ilpsta = nstaval(iw0)*nk+1
    ilpend = (nstaval(iw0)+nvalue(iw0))*nk
    do iplp=ilpsta,ilpend
      if (logic_dh) then                     !lp_head is in dbl_spa
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        if (lin == 7) call gsd_sequence_extspace(irw,ilw)
        if (lin == 14) call gtd_sequence_extspace(irw,ilw)
        if (lin == 26) call gdv_sequence_extspace(irw,ilw)
      else                                   !lp_head is in act_s
        if (log_prod == 3) then
          lphead = lpnew_head(iplp)
          ihyposl = jphyl(lphead)
          ihyposr = jphy(lphead)
          ndim = ihyl(ihyposl)
        else
          ihyposl = jphy(iplp)
          ihyposr = ihyposl
          ndim = ihy(ihyposl)
        end if

        iwal0 = lpnew_lwei(iplp)
        iwar0 = lpnew_rwei(iplp)
        do in_=1,ndim
          iwal = iwal0+ihyl(ihyposl+in_)
          iwar = iwar0+ihy(ihyposr+in_)
          do iwd=0,iwuplwei-1
            ilw = iwalk_ad(jpadl,ipael,iwal,iwd)
            irw = iwalk_ad(jpad,ipae,iwar,iwd)
            if (lin == 7) call gsd_sequence_extspace(irw,ilw)
            if (lin == 14) call gtd_sequence_extspace(irw,ilw)
            if (lin == 26) call gdv_sequence_extspace(irw,ilw)
          end do
        end do
      end if
    end do
  end if
end do

return

end subroutine ar_bl_br_ext_al_new

subroutine ar_drl_ext_al_new(lin,lri,lrk)

use gugaci_global, only: ihy, ihyl, ilsegdownwei, ipae, ipael, irsegdownwei, iseg_downwei, jpad, jpad_upwei, jpadl, jphy, jphyl, &
                         log_prod, logic_dh, logic_grad, lp_lwei, lp_rwei, lpnew_head, lpnew_lwei, lpnew_rwei, lsm_inn, mtype, &
                         ndim, nstaval, nvalue, value_lpext, value_lpext1, vplp_w0, vplpnew_w0, w0_sdplp
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lri, lrk
integer(kind=iwp) :: ihypos, ihyposl, ihyposr, iiext, ilpend, ilpsta, ilw, in_, iplp, irw, isma, iw0, iwal, iwal0, iwar, iwar0, &
                     iwd, iwuplwei, lphead, nlp_value
real(kind=wp) :: w0_sdold, w0multi
real(kind=wp), parameter :: crl = 1.0e-8_wp
integer(kind=iwp), external :: iwalk_ad

isma = lsm_inn(lri)
iwuplwei = jpad_upwei(jpad)
ilsegdownwei = iseg_downwei(ipae)
irsegdownwei = iseg_downwei(ipael)
w0_sdplp = vplpnew_w0(1)
if (logic_dh) w0_sdplp = vplp_w0(1)

if (logic_grad) then
  call lp678_ext_calcuvalue_g(lri,lrk,isma,nlp_value)
  w0_sdold = w0_sdplp
  do iw0=1,mtype
    if (iw0 /= 1) then
      w0_sdplp = vplpnew_w0(iw0)
      if (logic_dh) w0_sdplp = vplp_w0(iw0)
      if (abs(w0_sdplp) < crl) cycle
      w0multi = w0_sdplp/w0_sdold
      w0_sdold = w0_sdplp
      do iiext=1,nlp_value
        value_lpext(iiext) = value_lpext(iiext)*w0multi
        value_lpext1(iiext) = value_lpext1(iiext)*w0multi
      end do
    end if
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend
      if (logic_dh) then
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        if (lin == 7) call gsd_sequence_extspace_g(irw,ilw)
        if (lin == 14) call gtd_sequence_extspace_g(irw,ilw)
        if (lin == 26) call gdv_sequence_extspace_g(irw,ilw)
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
            if (lin == 7) call gsd_sequence_extspace_g(irw,ilw)
            if (lin == 14) call gtd_sequence_extspace_g(irw,ilw)
            if (lin == 26) call gdv_sequence_extspace_g(irw,ilw)
          end do
        end do
      end if
    end do
  end do

else
  call lp678_ext_wyb_calcuvalue(lri,lrk,isma,nlp_value)
  w0_sdold = w0_sdplp
  do iw0=1,mtype
    if (iw0 /= 1) then
      w0_sdplp = vplpnew_w0(iw0)
      if (logic_dh) w0_sdplp = vplp_w0(iw0)
      if (abs(w0_sdplp) < crl) cycle
      w0multi = w0_sdplp/w0_sdold
      w0_sdold = w0_sdplp
      do iiext=1,nlp_value
        value_lpext(iiext) = value_lpext(iiext)*w0multi
      end do
    end if
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend
      if (logic_dh) then
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        if (lin == 7) call gsd_sequence_extspace(irw,ilw)
        if (lin == 14) call gtd_sequence_extspace(irw,ilw)
        if (lin == 26) call gdv_sequence_extspace(irw,ilw)
      else
        if (log_prod == 3) then
          lphead = lpnew_head(iplp)
          ihyposl = jphyl(lphead)
          ihyposr = jphy(lphead)
          ndim = ihyl(ihyposl)
        else
          ihyposl = jphy(iplp)
          ihyposr = ihyposl
          ndim = ihy(ihyposl)
        end if

        iwal0 = lpnew_lwei(iplp)
        iwar0 = lpnew_rwei(iplp)
        do in_=1,ndim
          iwal = iwal0+ihyl(ihyposl+in_)
          iwar = iwar0+ihy(ihyposr+in_)
          do iwd=0,iwuplwei-1
            ilw = iwalk_ad(jpadl,ipael,iwal,iwd)
            irw = iwalk_ad(jpad,ipae,iwar,iwd)
            if (lin == 7) call gsd_sequence_extspace(irw,ilw)
            if (lin == 14) call gtd_sequence_extspace(irw,ilw)
            if (lin == 26) call gdv_sequence_extspace(irw,ilw)
          end do
        end do
      end if
    end do
  end do
end if

return

end subroutine ar_drl_ext_al_new

subroutine drr_br_ext_ar(lin,lrk,lri)

use gugaci_global, only: ihy, ihyl, ilsegdownwei, ipae, ipael, irsegdownwei, iseg_downwei, jpad, jpad_upwei, jpadl, jphy, jphyl, &
                         log_prod, logic_dh, logic_grad, lp_lwei, lp_rwei, lpnew_head, lpnew_lwei, lpnew_rwei, lsm_inn, mtype, &
                         ndim, nstaval, nvalue, vplp_w0, vplpnew_w0, w0_sdplp
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lrk, lri
integer(kind=iwp) :: ihypos, ihyposl, ihyposr, ilpend, ilpsta, ilw, in_, iplp, irw, isma, iw0, iwal, iwal0, iwar, iwar0, iwd, &
                     iwuplwei, lphead, nlp_value
real(kind=wp), parameter :: crl = 1.0e-8_wp
integer(kind=iwp), external :: iwalk_ad

isma = lsm_inn(lri)
iwuplwei = jpad_upwei(jpadl)
ilsegdownwei = iseg_downwei(ipael)
irsegdownwei = iseg_downwei(ipae)

do iw0=1,mtype
  w0_sdplp = vplpnew_w0(iw0)
  if (logic_dh) w0_sdplp = vplp_w0(iw0)
  if (abs(w0_sdplp) < crl) cycle
  if (logic_grad) then
    call lp678_ext_calcuvalue_g(lri,lrk,isma,nlp_value)
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend
      if (logic_dh) then
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        if (lin == 6) call gsd_sequence_extspace_g(ilw,irw)
        if (lin == 13) call gtd_sequence_extspace_g(ilw,irw)
        if (lin == 23) call gdv_sequence_extspace_g(ilw,irw)
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
            if (lin == 6) call gsd_sequence_extspace_g(ilw,irw)
            if (lin == 13) call gtd_sequence_extspace_g(ilw,irw)
            if (lin == 23) call gdv_sequence_extspace_g(ilw,irw)
          end do
        end do
      end if
    end do

  else
    call lp678_ext_wyb_calcuvalue(lri,lrk,isma,nlp_value)
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend
      if (logic_dh) then
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        if (lin == 6) call gsd_sequence_extspace(ilw,irw)
        if (lin == 13) call gtd_sequence_extspace(ilw,irw)
        if (lin == 23) call gdv_sequence_extspace(ilw,irw)
      else
        if (log_prod == 3) then
          lphead = lpnew_head(iplp)
          ihyposl = jphyl(lphead)
          ihyposr = jphy(lphead)
          ndim = ihyl(ihyposl)
        else
          ihyposl = jphy(iplp)
          ihyposr = ihyposl
          ndim = ihy(ihyposl)
        end if

        iwal0 = lpnew_lwei(iplp)
        iwar0 = lpnew_rwei(iplp)
        do in_=1,ndim
          iwal = iwal0+ihyl(ihyposl+in_)
          iwar = iwar0+ihy(ihyposr+in_)
          do iwd=0,iwuplwei-1
            ilw = iwalk_ad(jpadl,ipael,iwal,iwd)
            irw = iwalk_ad(jpad,ipae,iwar,iwd)
            if (lin == 6) call gsd_sequence_extspace(ilw,irw)
            if (lin == 13) call gtd_sequence_extspace(ilw,irw)
            if (lin == 23) call gdv_sequence_extspace(ilw,irw)
          end do
        end do
      end if
    end do
  end if
end do

return

end subroutine drr_br_ext_ar

subroutine ar_br_br_ext_ar_new(lin,intentry,isma)

use gugaci_global, only: ihy, ihyl, ilsegdownwei, ipae, ipael, irsegdownwei, iseg_downwei, jpad, jpad_upwei, jpadl, jphy, jphyl, &
                         log_prod, logic_dh, logic_grad, lp_lwei, lp_rwei, lpnew_head, lpnew_lwei, lpnew_rwei, mtype, ndim, &
                         nstaval, nvalue, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_sdplp, w1_sdplp
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, intentry, isma
integer(kind=iwp) :: ihypos, ihyposl, ihyposr, ilpend, ilpsta, ilw, in_, iplp, irw, iw0, iwal, iwal0, iwar, iwar0, iwd, iwuplwei, &
                     lphead, nlp_value
integer(kind=iwp), external :: iwalk_ad

iwuplwei = jpad_upwei(jpadl)
ilsegdownwei = iseg_downwei(ipael)
irsegdownwei = iseg_downwei(ipae)
do iw0=1,mtype
  w0_sdplp = vplpnew_w0(iw0)
  w1_sdplp = vplpnew_w1(iw0)
  if (logic_dh) w0_sdplp = vplp_w0(iw0)
  if (logic_dh) w1_sdplp = vplp_w1(iw0)
  if (logic_grad) then
    call lp10_arbrbr_ext_calcuvalue_g(intentry,isma,nlp_value)
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend
      if (logic_dh) then
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        if (lin == 6) call gsd_sequence_extspace_g(ilw,irw)
        if (lin == 13) call gtd_sequence_extspace_g(ilw,irw)
        if (lin == 23) call gdv_sequence_extspace_g(ilw,irw)
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
            if (lin == 6) call gsd_sequence_extspace_g(ilw,irw)
            if (lin == 13) call gtd_sequence_extspace_g(ilw,irw)
            if (lin == 23) call gdv_sequence_extspace_g(ilw,irw)
          end do
        end do
      end if
    end do

  else
    call lp10_arbrbr_ext_calcuvalue(intentry,isma,nlp_value)
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend
      if (logic_dh) then
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        if (lin == 6) call gsd_sequence_extspace(ilw,irw)
        if (lin == 13) call gtd_sequence_extspace(ilw,irw)
        if (lin == 23) call gdv_sequence_extspace(ilw,irw)
      else
        if (log_prod == 3) then
          lphead = lpnew_head(iplp)
          ihyposl = jphyl(lphead)
          ihyposr = jphy(lphead)
          ndim = ihyl(ihyposl)
        else
          ihyposl = jphy(iplp)
          ihyposr = ihyposl
          ndim = ihy(ihyposl)
        end if

        iwal0 = lpnew_lwei(iplp)
        iwar0 = lpnew_rwei(iplp)
        do in_=1,ndim
          iwal = iwal0+ihyl(ihyposl+in_)
          iwar = iwar0+ihy(ihyposr+in_)
          do iwd=0,iwuplwei-1
            ilw = iwalk_ad(jpadl,ipael,iwal,iwd)
            irw = iwalk_ad(jpad,ipae,iwar,iwd)
            if (lin == 6) call gsd_sequence_extspace(ilw,irw)
            if (lin == 13) call gtd_sequence_extspace(ilw,irw)
            if (lin == 23) call gdv_sequence_extspace(ilw,irw)
          end do
        end do
      end if
    end do
  end if
end do

return

end subroutine ar_br_br_ext_ar_new

subroutine ar_bl_bl_ext_ar_new(lin,intentry,isma,nk)

use gugaci_global, only: ihy, ihyl, ilsegdownwei, ipae, ipael, irsegdownwei, iseg_downwei, jpad, jpad_upwei, jpadl, jphy, jphyl, &
                         log_prod, logic_dh, logic_grad, lp_lwei, lp_rwei, lpnew_head, lpnew_lwei, lpnew_rwei, mtype, ndim, &
                         nstaval, nvalue, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_sdplp, w1_sdplp
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, intentry, isma, nk
integer(kind=iwp) :: ihypos, ihyposl, ihyposr, ilpend, ilpsta, ilw, in_, iplp, irw, iw0, iwal, iwal0, iwar, iwar0, iwd, iwuplwei, &
                     lphead, nlp_value
integer(kind=iwp), external :: iwalk_ad

iwuplwei = jpad_upwei(jpadl)
ilsegdownwei = iseg_downwei(ipael)
irsegdownwei = iseg_downwei(ipae)
do iw0=1,mtype
  w0_sdplp = vplpnew_w0(iw0)
  w1_sdplp = vplpnew_w1(iw0)
  if (logic_dh) w0_sdplp = vplp_w0(iw0)
  if (logic_dh) w1_sdplp = vplp_w1(iw0)
  if (logic_grad) then
    call lp12_arblbl_ext_calcuvalue_g(intentry,isma,nlp_value)
    ilpsta = nstaval(iw0)*nk+1
    ilpend = (nstaval(iw0)+nvalue(iw0))*nk
    do iplp=ilpsta,ilpend
      if (logic_dh) then                     !lp_head is in dbl_spa
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        if (lin == 6) call gsd_sequence_extspace_g(ilw,irw)
        if (lin == 13) call gtd_sequence_extspace_g(ilw,irw)
        if (lin == 23) call gdv_sequence_extspace_g(ilw,irw)
      else                                   !lp_head is in act_s
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
            if (lin == 6) call gsd_sequence_extspace_g(ilw,irw)
            if (lin == 13) call gtd_sequence_extspace_g(ilw,irw)
            if (lin == 23) call gdv_sequence_extspace_g(ilw,irw)
          end do
        end do
      end if
    end do

  else
    call lp12_arblbl_ext_calcuvalue(intentry,isma,nlp_value)
    ilpsta = nstaval(iw0)*nk+1
    ilpend = (nstaval(iw0)+nvalue(iw0))*nk
    do iplp=ilpsta,ilpend
      if (logic_dh) then                     !lp_head is in dbl_spa
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        if (lin == 6) call gsd_sequence_extspace(ilw,irw)
        if (lin == 13) call gtd_sequence_extspace(ilw,irw)
        if (lin == 23) call gdv_sequence_extspace(ilw,irw)
      else                                   !lp_head is in act_s
        if (log_prod == 3) then
          lphead = lpnew_head(iplp)
          ihyposl = jphyl(lphead)
          ihyposr = jphy(lphead)
          ndim = ihyl(ihyposl)
        else
          ihyposl = jphy(iplp)
          ihyposr = ihyposl
          ndim = ihy(ihyposl)
        end if

        iwal0 = lpnew_lwei(iplp)
        iwar0 = lpnew_rwei(iplp)
        do in_=1,ndim
          iwal = iwal0+ihyl(ihyposl+in_)
          iwar = iwar0+ihy(ihyposr+in_)
          do iwd=0,iwuplwei-1
            ilw = iwalk_ad(jpadl,ipael,iwal,iwd)
            irw = iwalk_ad(jpad,ipae,iwar,iwd)
            if (lin == 6) call gsd_sequence_extspace(ilw,irw)
            if (lin == 13) call gtd_sequence_extspace(ilw,irw)
            if (lin == 23) call gdv_sequence_extspace(ilw,irw)
          end do
        end do
      end if
    end do
  end if
end do

return

end subroutine ar_bl_bl_ext_ar_new

subroutine drl_br_sum_al_new(lin,lrp,lrq,lri)   !to be revised

use gugaci_global, only: ilsegdownwei, ipae, ipael, irsegdownwei, iseg_downwei, logic_grad, lp_lwei, lp_rwei, lsm_inn, mtype, &
                         norb_dz, nstaval, nvalue, value_lpext, value_lpext1, vplp_w0, w0_sdplp
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lrp, lrq, lri
integer(kind=iwp) :: iiext, ilpend, ilpsta, ilw, iplp, irw, isma, iw0, lrk, nlp_value
real(kind=wp) :: w0_sdold, w0multi
real(kind=wp), parameter :: crl = 1.0e-8_wp

if (logic_grad) then
  do lrk=1,norb_dz
    if (lrp == lrk) cycle
    if (lrq == lrk) cycle
    ilsegdownwei = iseg_downwei(ipae)
    irsegdownwei = iseg_downwei(ipael)
    isma = lsm_inn(lri)
    w0_sdplp = vplp_w0(1)
    call lp8_drlbr_sum_calcuvalue_g(lri,lrk,isma,nlp_value)
    w0_sdold = w0_sdplp
    do iw0=1,mtype
      if (iw0 /= 1) then
        w0_sdplp = vplp_w0(iw0)
        if (abs(w0_sdplp) < crl) cycle
        w0multi = w0_sdplp/w0_sdold
        w0_sdold = w0_sdplp
        do iiext=1,nlp_value
          value_lpext(iiext) = value_lpext(iiext)*w0multi
          value_lpext1(iiext) = value_lpext1(iiext)*w0multi
        end do
      end if
      ilpsta = nstaval(iw0)+1
      ilpend = nstaval(iw0)+nvalue(iw0)
      do iplp=ilpsta,ilpend
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        if (lin == 7) call gsd_sequence_extspace_g(irw,ilw)
        if (lin == 14) call gtd_sequence_extspace_g(irw,ilw)
        if (lin == 26) call gdv_sequence_extspace_g(irw,ilw)
      end do
    end do
  end do

else
  ilsegdownwei = iseg_downwei(ipae)
  irsegdownwei = iseg_downwei(ipael)
  isma = lsm_inn(lri)
  w0_sdplp = vplp_w0(1)
  call lp8_drlbr_sum_calcuvalue_wyb(lri,lrp,lrq,isma,nlp_value)
  w0_sdold = w0_sdplp
  do iw0=1,mtype
    if (iw0 /= 1) then
      w0_sdplp = vplp_w0(iw0)
      if (abs(w0_sdplp) < crl) cycle
      w0multi = w0_sdplp/w0_sdold
      w0_sdold = w0_sdplp
      do iiext=1,nlp_value
        value_lpext(iiext) = value_lpext(iiext)*w0multi
      end do
    end if
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend
      ilw = lp_lwei(iplp)
      irw = lp_rwei(iplp)
      if (lin == 7) call gsd_sequence_extspace(irw,ilw)
      if (lin == 14) call gtd_sequence_extspace(irw,ilw)
      if (lin == 26) call gdv_sequence_extspace(irw,ilw)
    end do
  end do
end if

return

end subroutine drl_br_sum_al_new

subroutine drl_bl_sum_ar_new(lin,lrp,lrq,lri)

use gugaci_global, only: ilsegdownwei, ipae, ipael, irsegdownwei, iseg_downwei, logic_grad, lp_lwei, lp_rwei, lsm_inn, mtype, &
                         norb_dz, nstaval, nvalue, value_lpext, value_lpext1, vplp_w0, w0_sdplp
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lrp, lrq, lri
integer(kind=iwp) :: iiext, ilpend, ilpsta, ilw, iplp, irw, isma, iw0, lrk, nlp_value
real(kind=wp) :: w0_sdold, w0multi
real(kind=wp), parameter :: crl = 1.0e-8_wp

if (logic_grad) then
  do lrk=1,norb_dz
    if (lrp == lrk) cycle
    if (lrq == lrk) cycle
    isma = lsm_inn(lri)
    ilsegdownwei = iseg_downwei(ipael)
    irsegdownwei = iseg_downwei(ipae)
    w0_sdplp = vplp_w0(1)
    call lp9_drlbl_sum_calcuvalue_g(lri,lrk,isma,nlp_value)
    w0_sdold = w0_sdplp
    do iw0=1,mtype
      if (iw0 /= 1) then
        w0_sdplp = vplp_w0(iw0)
        if (abs(w0_sdplp) < crl) cycle
        w0multi = w0_sdplp/w0_sdold
        w0_sdold = w0_sdplp
        do iiext=1,nlp_value
          value_lpext(iiext) = value_lpext(iiext)*w0multi
          value_lpext1(iiext) = value_lpext1(iiext)*w0multi
        end do
      end if
      ilpsta = nstaval(iw0)+1
      ilpend = nstaval(iw0)+nvalue(iw0)
      do iplp=ilpsta,ilpend
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        if (lin == 6) call gsd_sequence_extspace_g(ilw,irw)
        if (lin == 13) call gtd_sequence_extspace_g(ilw,irw)
        if (lin == 23) call gdv_sequence_extspace_g(ilw,irw)
      end do
    end do
  end do

else
  isma = lsm_inn(lri)
  ilsegdownwei = iseg_downwei(ipael)
  irsegdownwei = iseg_downwei(ipae)
  w0_sdplp = vplp_w0(1)

  call lp9_drlbl_sum_calcuvalue_wyb(lri,lrp,lrq,isma,nlp_value)
  w0_sdold = w0_sdplp
  do iw0=1,mtype
    if (iw0 /= 1) then
      w0_sdplp = vplp_w0(iw0)
      if (abs(w0_sdplp) < crl) cycle
      w0multi = w0_sdplp/w0_sdold
      w0_sdold = w0_sdplp
      do iiext=1,nlp_value
        value_lpext(iiext) = value_lpext(iiext)*w0multi
      end do
    end if
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend
      ilw = lp_lwei(iplp)
      irw = lp_rwei(iplp)
      if (lin == 6) call gsd_sequence_extspace(ilw,irw)
      if (lin == 13) call gtd_sequence_extspace(ilw,irw)
      if (lin == 23) call gdv_sequence_extspace(ilw,irw)
    end do
  end do
end if

return

end subroutine drl_bl_sum_ar_new

subroutine ar_bl_ext_ss(lri,lrj,nk)

use gugaci_global, only: ihy, ihyl, ilsegdownwei, iml, imr, ipae, ipael, irsegdownwei, iseg_downwei, jpad, jpad_upwei, jpadl, &
                         jphy, logic_dh, logic_g13, logic_g1415, logic_g2g4b, logic_g34b, logic_g35b, logic_g36b, logic_grad, &
                         lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, mtype, ndim, nstaval, nvalue, value_lpext, value_lpext1, &
                         vplp_w0, vplpnew_w0, w0_plp, w1_plp
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lri, lrj, nk
integer(kind=iwp) :: ihypos, iiext, ilpend, ilpsta, ilw, imlr, in_, iplp, irw, iw0, iwal, iwal0, iwar, iwar0, iwd, iwuplwei, &
                     nlp_value
real(kind=wp) :: w0_old, w0multi
real(kind=wp), parameter :: crl = 1.0e-8_wp
integer(kind=iwp), external :: iwalk_ad

iwuplwei = jpad_upwei(jpadl)
ilsegdownwei = iseg_downwei(ipael)
irsegdownwei = iseg_downwei(ipae)
imlr = Mul(iml,imr)
if (imlr == 1) then
  logic_g1415 = .true.
  if (iml == 1) logic_g13 = .true.
  if (iml == 1) logic_g2g4b = .true.
  logic_g36b = .true.
  logic_g35b = .true.
  logic_g34b = .true.
end if

w0_plp = vplpnew_w0(1)
w1_plp = Zero
if (logic_dh) w0_plp = vplp_w0(1)

if (logic_grad) then
  call lp_arbl_ext_st_calcuvalue_g(lri,lrj,nlp_value)
  w0_old = w0_plp
  do iw0=1,mtype
    if (iw0 /= 1) then
      w0_plp = vplpnew_w0(iw0)
      if (logic_dh) w0_plp = vplp_w0(iw0)
      if (abs(w0_plp) < crl) cycle
      w0multi = w0_plp/w0_old
      w0_old = w0_plp
      do iiext=1,nlp_value
        value_lpext(iiext) = value_lpext(iiext)*w0multi
        value_lpext1(iiext) = value_lpext1(iiext)*w0multi
      end do
    end if

    ilpsta = nstaval(iw0)*nk+1
    ilpend = (nstaval(iw0)+nvalue(iw0))*nk
    do iplp=ilpsta,ilpend
      if (logic_dh) then                    !lp_head is in dbl_space
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call inn_ext_ss_loop_unpack_g(ilw,irw)
      else                                  !lp_head is in act_spa
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
            call inn_ext_ss_loop_unpack_g(ilw,irw)
          end do
        end do
      end if
    end do
  end do
else
  call lp_arbl_ext_st_calcuvalue(lri,lrj,nlp_value)
  w0_old = w0_plp
  do iw0=1,mtype
    if (iw0 /= 1) then
      w0_plp = vplpnew_w0(iw0)
      if (logic_dh) w0_plp = vplp_w0(iw0)
      if (abs(w0_plp) < crl) cycle
      w0multi = w0_plp/w0_old
      w0_old = w0_plp
      do iiext=1,nlp_value
        value_lpext(iiext) = value_lpext(iiext)*w0multi
      end do
    end if

    ilpsta = nstaval(iw0)*nk+1
    ilpend = (nstaval(iw0)+nvalue(iw0))*nk
    do iplp=ilpsta,ilpend
      if (logic_dh) then                    !lp_head is in dbl_space
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call inn_ext_ss_loop_unpack(ilw,irw)
      else                                  !lp_head is in act_spa
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
            call inn_ext_ss_loop_unpack(ilw,irw)
          end do
        end do
      end if
    end do
  end do
end if

return

end subroutine ar_bl_ext_ss

subroutine ar_bl_ext_st(lri,lrj,nk)     !w0=0

use gugaci_global, only: ihy, ihyl, ilsegdownwei, ipae, ipael, irsegdownwei, iseg_downwei, jpad, jpad_upwei, jpadl, jphy, &
                         logic_dh, logic_grad, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, mtype, ndim, nstaval, nvalue, &
                         value_lpext, value_lpext1, vplp_w1, vplpnew_w1, w0_plp, w1_plp
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lri, lrj, nk
integer(kind=iwp) :: ihypos, iiext, ilpend, ilpsta, ilw, in_, iplp, irw, iw0, iwal, iwal0, iwar, iwar0, iwd, iwuplwei, nlp_value
real(kind=wp) :: w0multi, w1_old
real(kind=wp), parameter :: crl = 1.0e-8_wp
integer(kind=iwp), external :: iwalk_ad

iwuplwei = jpad_upwei(jpadl)
ilsegdownwei = iseg_downwei(ipael)
irsegdownwei = iseg_downwei(ipae)
w0_plp = Zero
w1_plp = vplpnew_w1(1)
if (logic_dh) w1_plp = vplp_w1(1)

if (logic_grad) then
  call lp_arbl_ext_st_calcuvalue_g(lri,lrj,nlp_value)
  w1_old = w1_plp
  do iw0=1,mtype
    if (iw0 /= 1) then
      w1_plp = vplpnew_w1(iw0)
      if (logic_dh) w1_plp = vplp_w1(iw0)
      if (abs(w1_plp) < crl) cycle
      w0multi = w1_plp/w1_old
      w1_old = w1_plp
      do iiext=1,nlp_value
        value_lpext(iiext) = value_lpext(iiext)*w0multi
        value_lpext1(iiext) = value_lpext1(iiext)*w0multi
      end do
    end if

    ilpsta = nstaval(iw0)*nk+1
    ilpend = (nstaval(iw0)+nvalue(iw0))*nk
    do iplp=ilpsta,ilpend
      if (logic_dh) then                    !lp_head is in dbl_space
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call inn_ext_st_loop_unpack_g(ilw,irw)
      else                                  !lp_head is in act_spa
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
            call inn_ext_st_loop_unpack_g(ilw,irw)
          end do
        end do
      end if
    end do
  end do
else
  call lp_arbl_ext_st_calcuvalue(lri,lrj,nlp_value)
  w1_old = w1_plp
  do iw0=1,mtype
    if (iw0 /= 1) then
      w1_plp = vplpnew_w1(iw0)
      if (logic_dh) w1_plp = vplp_w1(iw0)
      if (abs(w1_plp) < crl) cycle
      w0multi = w1_plp/w1_old
      w1_old = w1_plp
      do iiext=1,nlp_value
        value_lpext(iiext) = value_lpext(iiext)*w0multi
      end do
    end if

    ilpsta = nstaval(iw0)*nk+1
    ilpend = (nstaval(iw0)+nvalue(iw0))*nk
    do iplp=ilpsta,ilpend
      if (logic_dh) then                    !lp_head is in dbl_space
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call inn_ext_st_loop_unpack(ilw,irw)
      else                                  !lp_head is in act_spa
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
            call inn_ext_st_loop_unpack(ilw,irw)
          end do
        end do
      end if
    end do
  end do
end if

return

end subroutine ar_bl_ext_st

subroutine ar_bl_ext_ts(lri,lrj,nk)

use gugaci_global, only: ihy, ihyl, ilsegdownwei, ipae, ipael, irsegdownwei, iseg_downwei, jpad, jpad_upwei, jpadl, jphy, &
                         logic_dh, logic_grad, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, mtype, ndim, nstaval, nvalue, &
                         value_lpext, value_lpext1, vplp_w1, vplpnew_w1, w0_plp, w1_plp
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lri, lrj, nk
integer(kind=iwp) :: ihypos, iiext, ilpend, ilpsta, ilw, in_, iplp, irw, iw0, iwal, iwal0, iwar, iwar0, iwd, iwuplwei, nlp_value
real(kind=wp) :: w0multi, w1_old
real(kind=wp), parameter :: crl = 1.0e-8_wp
integer(kind=iwp), external :: iwalk_ad

iwuplwei = jpad_upwei(jpadl)
ilsegdownwei = iseg_downwei(ipael)
irsegdownwei = iseg_downwei(ipae)
w0_plp = Zero
w1_plp = vplpnew_w1(1)
if (logic_dh) w1_plp = vplp_w1(1)

if (logic_grad) then
  call lp_arbl_ext_st_calcuvalue_g(lri,lrj,nlp_value)
  w1_old = w1_plp
  do iw0=1,mtype
    if (iw0 /= 1) then
      w1_plp = vplpnew_w1(iw0)
      if (logic_dh) w1_plp = vplp_w1(iw0)
      if (abs(w1_plp) < crl) cycle
      w0multi = w1_plp/w1_old
      w1_old = w1_plp
      do iiext=1,nlp_value
        value_lpext(iiext) = value_lpext(iiext)*w0multi
        value_lpext1(iiext) = value_lpext1(iiext)*w0multi
      end do
    end if

    ilpsta = nstaval(iw0)*nk+1
    ilpend = (nstaval(iw0)+nvalue(iw0))*nk
    do iplp=ilpsta,ilpend
      if (logic_dh) then                    !lp_head is in dbl_space
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call inn_ext_ts_loop_unpack_g(ilw,irw)
      else                                  !lp_head is in act_spa
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
            call inn_ext_ts_loop_unpack_g(ilw,irw)
          end do
        end do
      end if
    end do
  end do
else
  call lp_arbl_ext_st_calcuvalue(lri,lrj,nlp_value)
  w1_old = w1_plp
  do iw0=1,mtype
    if (iw0 /= 1) then
      w1_plp = vplpnew_w1(iw0)
      if (logic_dh) w1_plp = vplp_w1(iw0)
      if (abs(w1_plp) < crl) cycle
      w0multi = w1_plp/w1_old
      w1_old = w1_plp
      do iiext=1,nlp_value
        value_lpext(iiext) = value_lpext(iiext)*w0multi
      end do
    end if

    ilpsta = nstaval(iw0)*nk+1
    ilpend = (nstaval(iw0)+nvalue(iw0))*nk
    do iplp=ilpsta,ilpend
      if (logic_dh) then                    !lp_head is in dbl_space
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call inn_ext_ts_loop_unpack(ilw,irw)
      else                                  !lp_head is in act_spa
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
            call inn_ext_ts_loop_unpack(ilw,irw)
          end do
        end do
      end if
    end do
  end do
end if

return

end subroutine ar_bl_ext_ts

subroutine ar_bl_ext_tt(lri,lrj,nk)

use gugaci_global, only: ihy, ihyl, ilsegdownwei, iml, imr, ipae, ipael, irsegdownwei, iseg_downwei, jpad, jpad_upwei, jpadl, &
                         jphy, logic_dh, logic_g1415, logic_g34b, logic_g35b, logic_g36b, logic_grad, lp_lwei, lp_rwei, &
                         lpnew_lwei, lpnew_rwei, mtype, ndim, nstaval, nvalue, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_plp, &
                         w1_plp
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lri, lrj, nk
integer(kind=iwp) :: ihypos, ilpend, ilpsta, ilw, imlr, in_, iplp, irw, iw0, iwal, iwal0, iwar, iwar0, iwd, iwuplwei, nlp_value
integer(kind=iwp), external :: iwalk_ad

iwuplwei = jpad_upwei(jpadl)
ilsegdownwei = iseg_downwei(ipael)
irsegdownwei = iseg_downwei(ipae)
imlr = Mul(iml,imr)
if (imlr == 1) then
  logic_g1415 = .true.
  logic_g36b = .true.
  logic_g35b = .true.
  logic_g34b = .true.
end if
do iw0=1,mtype
  w0_plp = vplpnew_w0(iw0)
  w1_plp = vplpnew_w1(iw0)
  if (logic_dh) w0_plp = vplp_w0(iw0)
  if (logic_dh) w1_plp = vplp_w1(iw0)

  if (logic_grad) then
    call lp_arbl_ext_st_calcuvalue_g(lri,lrj,nlp_value)
    ilpsta = nstaval(iw0)*nk+1
    ilpend = (nstaval(iw0)+nvalue(iw0))*nk
    do iplp=ilpsta,ilpend
      if (logic_dh) then                    !lp_head is in dbl_space
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call inn_ext_tt_loop_unpack_g(ilw,irw)
      else                                  !lp_head is in act_spa
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
            call inn_ext_tt_loop_unpack_g(ilw,irw)
          end do
        end do
      end if
    end do
  else
    call lp_arbl_ext_st_calcuvalue(lri,lrj,nlp_value)
    ilpsta = nstaval(iw0)*nk+1
    ilpend = (nstaval(iw0)+nvalue(iw0))*nk
    do iplp=ilpsta,ilpend
      if (logic_dh) then                    !lp_head is in dbl_space
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call inn_ext_tt_loop_unpack(ilw,irw)
      else                                  !lp_head is in act_spa
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
            call inn_ext_tt_loop_unpack(ilw,irw)
          end do
        end do
      end if
    end do

  end if
end do

return

end subroutine ar_bl_ext_tt

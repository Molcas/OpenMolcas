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

! bbs act-ext
subroutine ar_bl_dd_ext(lri,lrj,nk)

use gugaci_global, only: ihy, ihyl, ildownwei_segdd, ilsegdownwei, iml, imr, ipae, ipael, irdownwei_segdd, irsegdownwei, &
                         iseg_downwei, jpad, jpad_upwei, jpadl, jphy, logic_dh, logic_g49b, logic_grad, lp_lwei, lp_rwei, &
                         lpnew_lwei, lpnew_rwei, mtype, ndim, nstaval, nvalue, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_plp, &
                         w1_plp
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lri, lrj, nk
integer(kind=iwp) :: ihypos, ilpend, ilpsta, ilw, iml0, imr0, in_, iplp, irw, iw0, iwal, iwal0, iwar, iwar0, iwd, iwuplwei, &
                     nlp_value
integer(kind=iwp), external :: iwalk_ad

!write(nf2,*) 'ar_bl_dd_ext'
logic_g49b = .true.
iwuplwei = jpad_upwei(jpadl)
ilsegdownwei = iseg_downwei(ipael)
irsegdownwei = iseg_downwei(ipae)
ildownwei_segdd = iseg_downwei(ipael)
irdownwei_segdd = iseg_downwei(ipae)
iml0 = iml
imr0 = imr
do iw0=1,mtype
  w0_plp = vplpnew_w0(iw0)
  w1_plp = vplpnew_w1(iw0)
  if (logic_dh) w0_plp = vplp_w0(iw0)
  if (logic_dh) w1_plp = vplp_w1(iw0)

  if (logic_grad) then
    call lp_arbl_ext_dd_calcuvalue_g(lri,lrj,iml0,imr0,nlp_value)
    ilpsta = nstaval(iw0)*nk+1
    ilpend = (nstaval(iw0)+nvalue(iw0))*nk
    do iplp=ilpsta,ilpend
      if (logic_dh) then                    !lp_head is in dbl_space
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call inn_ext_dd_loop_unpack_g(ilw,irw)
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
            call inn_ext_dd_loop_unpack_g(ilw,irw)
          end do
        end do
      end if
    end do
  else
    call lp_arbl_ext_dd_calcuvalue(lri,lrj,iml0,imr0,nlp_value)
    ilpsta = nstaval(iw0)*nk+1
    ilpend = (nstaval(iw0)+nvalue(iw0))*nk
    do iplp=ilpsta,ilpend
      if (logic_dh) then                    !lp_head is in dbl_space
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call inn_ext_dd_loop_unpack(ilw,irw)
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
            call inn_ext_dd_loop_unpack(ilw,irw)
          end do
        end do
      end if
    end do
  end if
end do

return

end subroutine ar_bl_dd_ext

subroutine drl_dd_ext(lri)

use gugaci_global, only: ihy, ihyl, ildownwei_segdd, iml, ipae, ipael, irdownwei_segdd, iseg_downwei, jpad, jpad_upwei, jpadl, &
                         jphy, logic_dh, logic_g49b, logic_grad, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, mtype, ndim, nstaval, &
                         nvalue, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_plp, w1_plp
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lri
integer(kind=iwp) :: ihypos, ilpend, ilpsta, ilw, in_, iplp, irw, iw0, iwal, iwal0, iwar, iwar0, iwd, iwuplwei, nlp_value
integer(kind=iwp), external :: iwalk_ad

iwuplwei = jpad_upwei(jpadl)
ildownwei_segdd = iseg_downwei(ipael)
irdownwei_segdd = iseg_downwei(ipae)

do iw0=1,mtype
  w0_plp = vplpnew_w0(iw0)
  w1_plp = vplpnew_w1(iw0)
  if (logic_dh) w0_plp = vplp_w0(iw0)
  if (logic_dh) w1_plp = vplp_w1(iw0)
  ilpsta = nstaval(iw0)+1
  ilpend = nstaval(iw0)+nvalue(iw0)
  do iplp=ilpsta,ilpend
    ilw = lp_lwei(iplp)
    irw = lp_rwei(iplp)
    logic_g49b = .false.
    if (ilw /= irw) logic_g49b = .true.

    if (logic_grad) then
      !call lp_drl_ext_dd_calcuvalue_g(lri,iml,nlp_value)
      if (logic_dh) then                      !lp_head is in dbl_spa
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        logic_g49b = .false.
        if (ilw /= irw) logic_g49b = .true.
        call lp_drl_ext_dd_calcuvalue_g(lri,iml,nlp_value)
        call inn_ext_dd_loop_unpack_g(ilw,irw)
      else                                    !lp_head is in act_spa
        ihypos = jphy(iplp)
        ndim = ihy(ihypos)
        iwal0 = lpnew_lwei(iplp)
        iwar0 = lpnew_rwei(iplp)
        logic_g49b = .false.
        if (iwal0 /= iwar0) logic_g49b = .true.
        call lp_drl_ext_dd_calcuvalue_g(lri,iml,nlp_value)
        do in_=1,ndim
          iwal = iwal0+ihyl(ihypos+in_)
          iwar = iwar0+ihy(ihypos+in_)
          do iwd=0,iwuplwei-1
            ilw = iwalk_ad(jpadl,ipael,iwal,iwd)
            irw = iwalk_ad(jpad,ipae,iwar,iwd)
            call inn_ext_dd_loop_unpack_g(ilw,irw)
          end do
        end do
      end if
    else
      !call lp_drl_ext_dd_calcuvalue_wyb(lri,iml,nlp_value)
      if (logic_dh) then                      !lp_head is in dbl_spa
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        logic_g49b = .false.
        if (ilw /= irw) logic_g49b = .true.
        call lp_drl_ext_dd_calcuvalue_wyb(lri,iml,nlp_value)
        call inn_ext_dd_loop_unpack(ilw,irw)
      else                                    !lp_head is in act_spa
        ihypos = jphy(iplp)
        ndim = ihy(ihypos)
        iwal0 = lpnew_lwei(iplp)
        iwar0 = lpnew_rwei(iplp)
        logic_g49b = .false.
        if (iwal0 /= iwar0) logic_g49b = .true.
        call lp_drl_ext_dd_calcuvalue_wyb(lri,iml,nlp_value)
        do in_=1,ndim
          iwal = iwal0+ihyl(ihypos+in_)
          iwar = iwar0+ihy(ihypos+in_)
          do iwd=0,iwuplwei-1
            ilw = iwalk_ad(jpadl,ipael,iwal,iwd)
            irw = iwalk_ad(jpad,ipae,iwar,iwd)
            call inn_ext_dd_loop_unpack(ilw,irw)
          end do
        end do
      end if
    end if
  end do
end do

return

end subroutine drl_dd_ext

subroutine drl_ss_ext(lri)

use gugaci_global, only: ihy, ihyl, ildownwei_segdd, ilsegdownwei, ipae, ipael, irdownwei_segdd, irsegdownwei, iseg_downwei, jpad, &
                         jpad_upwei, jpadl, jphy, logic_dh, logic_g1415, logic_g2g4b, logic_g34b, logic_g35b, logic_g36b, &
                         logic_grad, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, mtype, ndim, nstaval, nvalue, value_lpext, &
                         value_lpext1, vplp_w0, vplpnew_w0, w0_plp
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lri
integer(kind=iwp) :: ihypos, iiext, ilpend, ilpsta, ilw, in_, iplp, irw, iw0, iwal, iwal0, iwar, iwar0, iwd, iwuplwei, nlp_value
real(kind=wp) :: w0_old, w0multi
real(kind=wp), parameter :: crl = 1.0e-8_wp
integer(kind=iwp), external :: iwalk_ad

logic_g1415 = .false.
logic_g2g4b = .false.
logic_g36b = .false.
logic_g35b = .false.
logic_g34b = .false.
iwuplwei = jpad_upwei(jpadl)
ilsegdownwei = iseg_downwei(ipael)
irsegdownwei = iseg_downwei(ipae)
ildownwei_segdd = iseg_downwei(ipael)
irdownwei_segdd = iseg_downwei(ipae)

w0_plp = vplpnew_w0(1)
if (logic_dh) w0_plp = vplp_w0(1)

if (logic_grad) then
  call lp_drl_ext_ss_calcuvalue_g(lri,nlp_value)
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

    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend
      if (logic_dh) then                    !lp_head is in dbl_space
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call inn_ext_ss_drl_loop_unpack_g(ilw,irw)
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
            call inn_ext_ss_drl_loop_unpack_g(ilw,irw)
          end do
        end do
      end if
    end do
  end do

else
  call lp_drl_ext_ss_calcuvalue(lri,nlp_value)
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

    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend
      if (logic_dh) then                    !lp_head is in dbl_space
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call inn_ext_ss_drl_loop_unpack(ilw,irw)
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
            call inn_ext_ss_drl_loop_unpack(ilw,irw)
          end do
        end do
      end if
    end do
  end do
end if

return

end subroutine drl_ss_ext

subroutine drl_ss_sum(lri,lrj)

use gugaci_global, only: ildownwei_segdd, ilsegdownwei, ipae, ipael, irdownwei_segdd, irsegdownwei, iseg_downwei, logic_g1415, &
                         logic_g2g4b, logic_g34b, logic_g35b, logic_g36b, logic_grad, lp_lwei, lp_rwei, mtype, norb_dz, nstaval, &
                         nvalue, value_lpext, value_lpext1, vplp_w0, w0_plp
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lri, lrj
integer(kind=iwp) :: iiext, ilpend, ilpsta, ilw, iplp, irw, iw0, lrk, nlp_value
real(kind=wp) :: w0_old, w0multi
real(kind=wp), parameter :: crl = 1.0e-8_wp

logic_g1415 = .false.
logic_g2g4b = .false.
logic_g36b = .false.
logic_g35b = .false.
logic_g34b = .false.

if (logic_grad) then
  do lrk=1,norb_dz
    if (lri == lrk) cycle
    if (lrj == lrk) cycle
    ilsegdownwei = iseg_downwei(ipael)
    irsegdownwei = iseg_downwei(ipae)
    ildownwei_segdd = iseg_downwei(ipael)
    irdownwei_segdd = iseg_downwei(ipae)
    w0_plp = vplp_w0(1)
    call lp_drl_ext_ss_calcuvalue_g(lrk,nlp_value)
    w0_old = w0_plp
    do iw0=1,mtype
      if (iw0 /= 1) then
        w0_plp = vplp_w0(iw0)
        if (abs(w0_plp) < crl) cycle
        w0multi = w0_plp/w0_old
        w0_old = w0_plp
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
        call inn_ext_ss_drl_loop_unpack_g(ilw,irw)
      end do
    end do
  end do

else
  ilsegdownwei = iseg_downwei(ipael)
  irsegdownwei = iseg_downwei(ipae)
  ildownwei_segdd = iseg_downwei(ipael)
  irdownwei_segdd = iseg_downwei(ipae)
  w0_plp = vplp_w0(1)
  call lp_drl_sum_ss_calcuvalue(lri,lrj,nlp_value)
  w0_old = w0_plp
  do iw0=1,mtype
    if (iw0 /= 1) then
      w0_plp = vplp_w0(iw0)
      if (abs(w0_plp) < crl) cycle
      w0multi = w0_plp/w0_old
      w0_old = w0_plp
      do iiext=1,nlp_value
        value_lpext(iiext) = value_lpext(iiext)*w0multi
      end do
    end if
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend
      ilw = lp_lwei(iplp)
      irw = lp_rwei(iplp)
      call inn_ext_ss_drl_loop_unpack(ilw,irw)
    end do
  end do
end if

return

end subroutine drl_ss_sum

subroutine drl_st_ext(lri)

use gugaci_global, only: ihy, ihyl, ildownwei_segdd, ilsegdownwei, ipae, ipael, irdownwei_segdd, irsegdownwei, iseg_downwei, jpad, &
                         jpad_upwei, jpadl, jphy, logic_dh, logic_grad, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, mtype, ndim, &
                         nstaval, nvalue, value_lpext, value_lpext1, vplp_w1, vplpnew_w1, w1_plp
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lri
integer(kind=iwp) :: ihypos, iiext, ilpend, ilpsta, ilw, in_, iplp, irw, iw0, iwal, iwal0, iwar, iwar0, iwd, iwuplwei, nlp_value
real(kind=wp) :: w1_old, w1multi
real(kind=wp), parameter :: crl = 1.0e-8_wp
integer(kind=iwp), external :: iwalk_ad

iwuplwei = jpad_upwei(jpadl)
ilsegdownwei = iseg_downwei(ipael)
irsegdownwei = iseg_downwei(ipae)
ildownwei_segdd = iseg_downwei(ipael)
irdownwei_segdd = iseg_downwei(ipae)
!logic_g1415 = .false.
!logic_g34b = .false.
!logic_g35b = .false.
!logic_g36b = .false.

w1_plp = vplpnew_w1(1)
if (logic_dh) w1_plp = vplp_w1(1)

if (logic_grad) then
  call lp_drl_ext_st_calcuvalue_g(lri,nlp_value)
  w1_old = w1_plp
  do iw0=1,mtype
    if (iw0 /= 1) then
      w1_plp = vplpnew_w1(iw0)
      if (logic_dh) w1_plp = vplp_w1(iw0)
      if (abs(w1_plp) < crl) cycle
      w1multi = w1_plp/w1_old
      w1_old = w1_plp

      do iiext=1,nlp_value
        value_lpext(iiext) = value_lpext(iiext)*w1multi
        value_lpext1(iiext) = value_lpext1(iiext)*w1multi
      end do
    end if
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend

      if (logic_dh) then                    !lp_head is in dbl_space
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call inn_ext_st_drl_loop_unpack_g(ilw,irw)
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
            call inn_ext_st_drl_loop_unpack_g(ilw,irw)
          end do
        end do
      end if
    end do
  end do

else
  call lp_drl_ext_st_calcuvalue(lri,nlp_value)
  w1_old = w1_plp
  do iw0=1,mtype
    if (iw0 /= 1) then
      w1_plp = vplpnew_w1(iw0)
      if (logic_dh) w1_plp = vplp_w1(iw0)
      if (abs(w1_plp) < crl) cycle
      w1multi = w1_plp/w1_old
      w1_old = w1_plp
      !call calcu_drl2_value_wyb(1,lri,lrj)
      !call lp_drl_ext_st_calcuvalue(lri,nlp_value)
      do iiext=1,nlp_value
        value_lpext(iiext) = value_lpext(iiext)*w1multi
      end do
    end if
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend
      !write(u6,*) '      calcuvalue_w',intentry

      if (logic_dh) then                    !lp_head is in dbl_space
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call inn_ext_st_drl_loop_unpack(ilw,irw)
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
            call inn_ext_st_drl_loop_unpack(ilw,irw)
          end do
        end do
      end if
    end do
  end do

end if

return

end subroutine drl_st_ext

subroutine drl_tt_ext(lri)

use gugaci_global, only: ihy, ihyl, ildownwei_segdd, ilsegdownwei, ipae, ipael, irdownwei_segdd, irsegdownwei, iseg_downwei, jpad, &
                         jpad_upwei, jpadl, jphy, logic_dh, logic_grad, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, mtype, ndim, &
                         nstaval, nvalue, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_plp, w1_plp
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lri
integer(kind=iwp) :: ihypos, ilpend, ilpsta, ilw, in_, iplp, irw, iw0, iwal, iwal0, iwar, iwar0, iwd, iwuplwei, n1415, nlp_value
integer(kind=iwp), external :: iwalk_ad

!logic_g1415 = .false.
!logic_g2g4b = .false.
!logic_g36b = .false.
!logic_g35b = .false.
!logic_g34b = .false.
iwuplwei = jpad_upwei(jpadl)
ilsegdownwei = iseg_downwei(ipael)
irsegdownwei = iseg_downwei(ipae)
ildownwei_segdd = iseg_downwei(ipael)
irdownwei_segdd = iseg_downwei(ipae)

do iw0=1,mtype
  w0_plp = vplpnew_w0(iw0)
  w1_plp = vplpnew_w1(iw0)
  if (logic_dh) w0_plp = vplp_w0(iw0)
  if (logic_dh) w1_plp = vplp_w1(iw0)

  if (logic_grad) then
    call lp_drl_ext_tt_calcuvalue_g(lri,n1415,nlp_value)
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend

      if (logic_dh) then                    !lp_head is in dbl_space
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call inn_ext_tt_drl_loop_unpack_g(ilw,irw,n1415)
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
            call inn_ext_tt_drl_loop_unpack_g(ilw,irw,n1415)
          end do
        end do
      end if
    end do
  else
    call lp_drl_ext_tt_calcuvalue(lri,n1415,nlp_value)
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend

      if (logic_dh) then                    !lp_head is in dbl_space
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call inn_ext_tt_drl_loop_unpack(ilw,irw,n1415)
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
            call inn_ext_tt_drl_loop_unpack(ilw,irw,n1415)
          end do
        end do
      end if
    end do
  end if
end do

return

end subroutine drl_tt_ext

subroutine drl_tt_sum(lri,lrj)

use gugaci_global, only: ildownwei_segdd, ilsegdownwei, ipae, ipael, irdownwei_segdd, irsegdownwei, iseg_downwei, logic_grad, &
                         lp_lwei, lp_rwei, mtype, norb_dz, nstaval, nvalue, value_lpext, value_lpext1, vplp_w0, w0_plp
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lri, lrj
integer(kind=iwp) :: iiext, ilpend, ilpsta, ilw, iplp, irw, iw0, lrk, n1415, nlp_value
real(kind=wp) :: w0_old, w0multi
real(kind=wp), parameter :: crl = 1.0e-8_wp

!logic_g1415 = .false.
!logic_g2g4b = .false.
!logic_g36b = .false.
!logic_g35b = .false.
!logic_g34b = .false.
!write(nf2,*) 'drl_tt_sum'
if (logic_grad) then
  do lrk=1,norb_dz
    if (lri == lrk) cycle
    if (lrj == lrk) cycle
    ilsegdownwei = iseg_downwei(ipael)
    irsegdownwei = iseg_downwei(ipae)
    ildownwei_segdd = iseg_downwei(ipael)
    irdownwei_segdd = iseg_downwei(ipae)
    w0_plp = vplp_w0(1)
    call lp_drl_sum_tt_calcuvalue_g(lrk,n1415,nlp_value)
    w0_old = w0_plp
    do iw0=1,mtype
      if (iw0 /= 1) then
        w0_plp = vplp_w0(iw0)
        if (abs(w0_plp) < crl) cycle
        w0multi = w0_plp/w0_old
        w0_old = w0_plp
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
        call inn_ext_tt_drl_loop_unpack_g(ilw,irw,n1415)
      end do
    end do
  end do

else
  ilsegdownwei = iseg_downwei(ipael)
  irsegdownwei = iseg_downwei(ipae)
  ildownwei_segdd = iseg_downwei(ipael)
  irdownwei_segdd = iseg_downwei(ipae)
  w0_plp = vplp_w0(1)
  call lp_drl_sum_tt_calcuvalue(lri,lrj,n1415,nlp_value)
  w0_old = w0_plp
  do iw0=1,mtype
    if (iw0 /= 1) then
      w0_plp = vplp_w0(iw0)
      if (abs(w0_plp) < crl) cycle
      w0multi = w0_plp/w0_old
      w0_old = w0_plp
      do iiext=1,nlp_value
        value_lpext(iiext) = value_lpext(iiext)*w0multi
      end do
      !call lp_drl_sum_tt_calcuvalue(lri,lrj,n1415,nlp_value)
    end if
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend
      ilw = lp_lwei(iplp)
      irw = lp_rwei(iplp)
      call inn_ext_tt_drl_loop_unpack(ilw,irw,n1415)
    end do
  end do
end if

return

end subroutine drl_tt_sum

subroutine drl_ts_ext(lri)

use gugaci_global, only: ihy, ihyl, ildownwei_segdd, ilsegdownwei, ipae, ipael, irdownwei_segdd, irsegdownwei, iseg_downwei, jpad, &
                         jpad_upwei, jpadl, jphy, logic_dh, logic_grad, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, mtype, ndim, &
                         nstaval, nvalue, value_lpext, value_lpext1, vplp_w1, vplpnew_w1, w1_plp
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lri
integer(kind=iwp) :: ihypos, iiext, ilpend, ilpsta, ilw, in_, iplp, irw, iw0, iwal, iwal0, iwar, iwar0, iwd, iwuplwei, nlp_value
real(kind=wp) :: w1_old, w1multi
real(kind=wp), parameter :: crl = 1.0e-8_wp
integer(kind=iwp), external :: iwalk_ad

!logic_g1415 = .false.
!logic_g2g4b = .false.
!logic_g36b = .false.
!logic_g35b = .false.
!logic_g34b = .false.
ildownwei_segdd = iseg_downwei(ipael)
irdownwei_segdd = iseg_downwei(ipae)
iwuplwei = jpad_upwei(jpadl)
ilsegdownwei = iseg_downwei(ipael)
irsegdownwei = iseg_downwei(ipae)
w1_plp = vplpnew_w1(1)
if (logic_dh) w1_plp = vplp_w1(1)

if (logic_grad) then
  call lp_drl_ext_ts_calcuvalue_g(lri,nlp_value)
  w1_old = w1_plp
  do iw0=1,mtype
    if (iw0 /= 1) then
      w1_plp = vplpnew_w1(iw0)
      if (logic_dh) w1_plp = vplp_w1(iw0)
      if (abs(w1_plp) < crl) cycle
      w1multi = w1_plp/w1_old
      w1_old = w1_plp
      do iiext=1,nlp_value
        value_lpext(iiext) = value_lpext(iiext)*w1multi
        value_lpext1(iiext) = value_lpext1(iiext)*w1multi
      end do
    end if
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend

      if (logic_dh) then                    !lp_head is in dbl_space
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call inn_ext_ts_drl_loop_unpack_g(ilw,irw)
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
            call inn_ext_ts_drl_loop_unpack_g(ilw,irw)
          end do
        end do
      end if
    end do
  end do

else
  call lp_drl_ext_ts_calcuvalue(lri,nlp_value)
  w1_old = w1_plp
  do iw0=1,mtype
    if (iw0 /= 1) then
      w1_plp = vplpnew_w1(iw0)
      if (logic_dh) w1_plp = vplp_w1(iw0)
      if (abs(w1_plp) < crl) cycle
      w1multi = w1_plp/w1_old
      w1_old = w1_plp
      do iiext=1,nlp_value
        value_lpext(iiext) = value_lpext(iiext)*w1multi
      end do
    end if
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend

      if (logic_dh) then                    !lp_head is in dbl_space
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call inn_ext_ts_drl_loop_unpack(ilw,irw)
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
            call inn_ext_ts_drl_loop_unpack(ilw,irw)
          end do
        end do
      end if
    end do
  end do

end if

return

end subroutine drl_ts_ext

subroutine ar_br_tv_ext_br_ar(lri,lrj)

use gugaci_global, only: ihy, ihyl, ilsegdownwei, ipae, ipael, irsegdownwei, iseg_downwei, jpad, jpad_upwei, jpadl, jphy, jphyl, &
                         log_prod, logic_dh, logic_grad, lp_lwei, lp_rwei, lpnew_head, lpnew_lwei, lpnew_rwei, mtype, ndim, &
                         nstaval, nvalue, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_plp, w0g36a, w1_plp, w1g36a
use Constants, only: Zero, One
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lri, lrj
integer(kind=iwp) :: ihypos, ihyposl, ihyposr, ilpend, ilpsta, ilw, in_, iplp, irw, iw0, iwal, iwal0, iwar, iwar0, iwd, iwuplwei, &
                     lphead, nlp_value
integer(kind=iwp), external :: iwalk_ad

ilsegdownwei = iseg_downwei(ipael)
irsegdownwei = iseg_downwei(ipae)
iwuplwei = jpad_upwei(jpadl)
w0g36a = Zero
w1g36a = -One
do iw0=1,mtype
  w0_plp = vplpnew_w0(iw0)
  w1_plp = vplpnew_w1(iw0)
  if (logic_dh) w0_plp = vplp_w0(iw0)
  if (logic_dh) w1_plp = vplp_w1(iw0)

  if (logic_grad) then
    call lp_arbr_ext_svtv_calcuvalue_g(lri,lrj,nlp_value)
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend
      if (logic_dh) then
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call inn_ext_sv_loop_unpack_g(ilw,irw)
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
            call inn_ext_sv_loop_unpack_g(ilw,irw)
          end do
        end do
      end if
    end do
  else
    call lp_arbr_ext_svtv_calcuvalue_wyb(lri,lrj,nlp_value)
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend
      if (logic_dh) then
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call inn_ext_sv_loop_unpack(ilw,irw)
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
            call inn_ext_sv_loop_unpack(ilw,irw)
          end do
        end do
      end if
    end do
  end if
end do

return

end subroutine ar_br_tv_ext_br_ar

subroutine drr_sv_ext_br_ar(lri)

use gugaci_global, only: ihy, ihyl, ilsegdownwei, ipae, ipael, irsegdownwei, iseg_downwei, jpad, jpad_upwei, jpadl, jphy, jphyl, &
                         log_prod, logic_dh, logic_grad, lp_lwei, lp_rwei, lpnew_head, lpnew_lwei, lpnew_rwei, mtype, ndim, &
                         norb_dz, nstaval, nvalue, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_plp, w1_plp
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lri
integer(kind=iwp) :: ihypos, ihyposl, ihyposr, ilpend, ilpsta, ilw, in_, iplp, irw, iw0, iwal, iwal0, iwar, iwar0, iwd, iwuplwei, &
                     lphead, ndorb, nlp_value
integer(kind=iwp), external :: iwalk_ad

ndorb = norb_dz-lri
ilsegdownwei = iseg_downwei(ipael)
irsegdownwei = iseg_downwei(ipae)
iwuplwei = jpad_upwei(jpadl)
do iw0=1,mtype
  w0_plp = vplpnew_w0(iw0)
  w1_plp = vplpnew_w1(iw0)
  if (ndorb >= 0) w0_plp = vplp_w0(iw0)
  if (ndorb >= 0) w1_plp = vplp_w1(iw0)

  if (logic_grad) then
    call lp_drr_ext_svtv_calcuvalue_g(lri,nlp_value)
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend
      if (logic_dh) then
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call inn_ext_sv_loop_unpack_g(ilw,irw)
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
            call inn_ext_sv_loop_unpack_g(ilw,irw)
          end do
        end do
      end if
    end do
  else
    call lp_drr_ext_svtv_calcuvalue_wyb(lri,nlp_value)
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend
      if (logic_dh) then
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call inn_ext_sv_loop_unpack(ilw,irw)
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
            call inn_ext_sv_loop_unpack(ilw,irw)
          end do
        end do
      end if
    end do
  end if
end do

return

end subroutine drr_sv_ext_br_ar

subroutine drr_tv_ext_br_ar(lri)

use gugaci_global, only: ihy, ihyl, ilsegdownwei, ipae, ipael, irsegdownwei, iseg_downwei, jpad, jpad_upwei, jpadl, jphy, jphyl, &
                         log_prod, logic_dh, logic_grad, lp_lwei, lp_rwei, lpnew_head, lpnew_lwei, lpnew_rwei, mtype, ndim, &
                         norb_dz, nstaval, nvalue, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_plp, w1_plp
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lri
integer(kind=iwp) :: ihypos, ihyposl, ihyposr, ilpend, ilpsta, ilw, in_, iplp, irw, iw0, iwal, iwal0, iwar, iwar0, iwd, iwuplwei, &
                     lphead, ndorb, nlp_value
integer(kind=iwp), external :: iwalk_ad

ndorb = norb_dz-lri
ilsegdownwei = iseg_downwei(ipael)
irsegdownwei = iseg_downwei(ipae)
iwuplwei = jpad_upwei(jpadl)
do iw0=1,mtype
  w0_plp = vplpnew_w0(iw0)
  w1_plp = vplpnew_w1(iw0)
  if (ndorb >= 0) w0_plp = vplp_w0(iw0)
  if (ndorb >= 0) w1_plp = vplp_w1(iw0)

  if (logic_grad) then
    call lp_drr_ext_svtv_calcuvalue_g(lri,nlp_value)
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend
      if (logic_dh) then
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call inn_ext_sv_loop_unpack_g(ilw,irw)
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
            call inn_ext_sv_loop_unpack_g(ilw,irw)
          end do
        end do
      end if
    end do
  else
    call lp_drr_ext_svtv_calcuvalue_wyb(lri,nlp_value)
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend
      if (logic_dh) then
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call inn_ext_sv_loop_unpack(ilw,irw)
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
            call inn_ext_sv_loop_unpack(ilw,irw)
          end do
        end do
      end if
    end do
  end if
end do

return

end subroutine drr_tv_ext_br_ar

subroutine ar_dv_ext_ar(idtu,isma,lri,lrj)

use gugaci_global, only: ihy, ihyl, ilsegdownwei, ipae, ipael, irsegdownwei, iseg_downwei, jpad, jpad_upwei, jpadl, jphy, jphyl, &
                         log_prod, logic_dh, logic_grad, lp_lwei, lp_rwei, lpnew_coe, lpnew_head, lpnew_lwei, lpnew_rwei, mtype, &
                         ndim, norb_dz, norb_inn, nstaval, nvalue, vplp_w0, vplpnew_w0, w0_sdplp
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: idtu, isma, lri, lrj
integer(kind=iwp) :: ihypos, ihyposl, ihyposr, ilpend, ilpsta, ilw, in_, iplp, irw, iw0, iwal, iwal0, iwar, iwar0, iwd, iwuplwei, &
                     lphead, nlp_value, norb, nvalue1
integer(kind=iwp), allocatable :: lpcoe(:)
integer(kind=iwp), external :: iwalk_ad

iwuplwei = jpad_upwei(jpadl)
ilsegdownwei = iseg_downwei(ipael)
irsegdownwei = iseg_downwei(ipae)
call mma_allocate(lpcoe,[norb_dz+1,norb_inn],label='lpcoe')
do iw0=1,mtype
  w0_sdplp = vplpnew_w0(iw0)
  if (logic_dh) w0_sdplp = vplp_w0(iw0)
  ilpsta = nstaval(iw0)+1
  ilpend = nstaval(iw0)+nvalue(iw0)
  do iplp=ilpsta,ilpend
    !ipcoe = lp_arpos(iplp)
    do norb=norb_dz+1,norb_inn
      lpcoe(norb) = lpnew_coe(norb,iplp)  !01.12.25
    end do

    if (logic_grad) then
      call lp_ar_coe_calcuvalue_g(idtu,isma,lri,lrj,nlp_value,lpcoe,nvalue1)
      if (logic_dh) then                  !lp_head is in dbl_space
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call gdv_sequence_extspace1_g(ilw,irw,nvalue1)
      else                                !lp_head is in act_spa
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
            call gdv_sequence_extspace1_g(ilw,irw,nvalue1)
          end do
        end do
      end if
    else
      call lp_ar_coe_calcuvalue_wyb(idtu,isma,lri,lrj,nlp_value,lpcoe)
      if (logic_dh) then                  !lp_head is in dbl_space
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call gdv_sequence_extspace(ilw,irw)
      else                                !lp_head is in act_space
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
            call gdv_sequence_extspace(ilw,irw)
          end do
        end do
      end if
    end if
  end do
end do
call mma_deallocate(lpcoe)

return

end subroutine ar_dv_ext_ar

subroutine ar_sd_ext_ar(idtu,lri,lrj,isma)

use gugaci_global, only: ihy, ihyl, ilsegdownwei, ipae, ipael, irsegdownwei, iseg_downwei, jpad, jpad_upwei, jpadl, jphy, &
                         logic_dh, logic_grad, lp_lwei, lp_rwei, lpnew_coe, lpnew_lwei, lpnew_rwei, mtype, ndim, norb_dz, &
                         norb_inn, nstaval, nvalue, vplp_w0, vplpnew_w0, w0_sdplp
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: idtu, lri, lrj, isma
integer(kind=iwp) :: ihypos, ilpend, ilpsta, ilw, in_, iplp, irw, iw0, iwal, iwal0, iwar, iwar0, iwd, iwuplwei, nlp_value, norb, &
                     nvalue1
integer(kind=iwp), allocatable :: lpcoe(:)
integer(kind=iwp), external :: iwalk_ad

iwuplwei = jpad_upwei(jpadl)
ilsegdownwei = iseg_downwei(ipael)
irsegdownwei = iseg_downwei(ipae)
call mma_allocate(lpcoe,[norb_dz+1,norb_inn],label='lpcoe')
do iw0=1,mtype
  w0_sdplp = vplpnew_w0(iw0)
  if (logic_dh) w0_sdplp = vplp_w0(iw0)
  ilpsta = nstaval(iw0)+1
  ilpend = nstaval(iw0)+nvalue(iw0)
  do iplp=ilpsta,ilpend

    do norb=norb_dz+1,norb_inn
      lpcoe(norb) = lpnew_coe(norb,iplp)
    end do

    if (logic_grad) then
      call lp_ar_coe_calcuvalue_g(idtu,isma,lri,lrj,nlp_value,lpcoe,nvalue1)
      if (logic_dh) then              !lp_head is in dbl_space
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call gsd_sequence_extspace1_g(ilw,irw,nvalue1)
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
            call gsd_sequence_extspace1_g(ilw,irw,nvalue1)
          end do
        end do
      end if
    else
      call lp_ar_coe_calcuvalue_wyb(idtu,isma,lri,lrj,nlp_value,lpcoe)
      if (logic_dh) then              !lp_head is in dbl_space
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call gsd_sequence_extspace(ilw,irw)
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
            call gsd_sequence_extspace(ilw,irw)
          end do
        end do
      end if
    end if
  end do
end do
call mma_deallocate(lpcoe)

return

end subroutine ar_sd_ext_ar

subroutine ar_td_ext_ar(idtu,lri,lrj,isma)

use gugaci_global, only: ihy, ihyl, ilsegdownwei, ipae, ipael, irsegdownwei, iseg_downwei, jpad, jpad_upwei, jpadl, jphy, &
                         logic_dh, logic_grad, lp_lwei, lp_rwei, lpnew_coe, lpnew_lwei, lpnew_rwei, mtype, ndim, norb_dz, &
                         norb_inn, nstaval, nvalue, vplp_w0, vplpnew_w0, w0_sdplp
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: idtu, lri, lrj, isma
integer(kind=iwp) :: ihypos, ilpend, ilpsta, ilw, in_, iplp, irw, iw0, iwal, iwal0, iwar, iwar0, iwd, iwuplwei, nlp_value, norb, &
                     nvalue1
integer(kind=iwp), allocatable :: lpcoe(:)
integer(kind=iwp), external :: iwalk_ad

iwuplwei = jpad_upwei(jpadl)
ilsegdownwei = iseg_downwei(ipael)
irsegdownwei = iseg_downwei(ipae)
call mma_allocate(lpcoe,[norb_dz+1,norb_inn],label='lpcoe')
do iw0=1,mtype
  w0_sdplp = vplpnew_w0(iw0)
  if (logic_dh) w0_sdplp = vplp_w0(iw0)
  ilpsta = nstaval(iw0)+1
  ilpend = nstaval(iw0)+nvalue(iw0)
  do iplp=ilpsta,ilpend
    do norb=norb_dz+1,norb_inn
      lpcoe(norb) = lpnew_coe(norb,iplp)
    end do

    if (logic_grad) then
      call lp_ar_coe_calcuvalue_g(idtu,isma,lri,lrj,nlp_value,lpcoe,nvalue1)
      if (logic_dh) then              !lp_head is in dbl_space
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)

        call gtd_sequence_extspace1_g(ilw,irw,nvalue1)
      else                            !lp_head is in act_space
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
            call gtd_sequence_extspace1_g(ilw,irw,nvalue1)
          end do
        end do
      end if
    else
      call lp_ar_coe_calcuvalue_wyb(idtu,isma,lri,lrj,nlp_value,lpcoe)
      if (logic_dh) then             !lp_head is in dbl_space
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call gtd_sequence_extspace(ilw,irw)
      else                           !lp_head is in act_space
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
            call gtd_sequence_extspace(ilw,irw)
          end do
        end do
      end if
    end if
  end do
end do
call mma_deallocate(lpcoe)

return

end subroutine ar_td_ext_ar

subroutine ar_br_sv_ext_br_ar(lri,lrj)

use gugaci_global, only: ihy, ihyl, ilsegdownwei, ipae, ipael, irsegdownwei, iseg_downwei, jpad, jpad_upwei, jpadl, jphy, jphyl, &
                         log_prod, logic_dh, logic_grad, lp_lwei, lp_rwei, lpnew_head, lpnew_lwei, lpnew_rwei, mtype, ndim, &
                         nstaval, nvalue, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_plp, w1_plp
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lri, lrj
integer(kind=iwp) :: ihypos, ihyposl, ihyposr, ilpend, ilpsta, ilw, in_, iplp, irw, iw0, iwal, iwal0, iwar, iwar0, iwd, iwuplwei, &
                     lphead, nlp_value
integer(kind=iwp), external :: iwalk_ad

ilsegdownwei = iseg_downwei(ipael)
irsegdownwei = iseg_downwei(ipae)
iwuplwei = jpad_upwei(jpadl)
do iw0=1,mtype
  w0_plp = vplpnew_w0(iw0)
  w1_plp = vplpnew_w1(iw0)
  if (logic_dh) w0_plp = vplp_w0(iw0)
  if (logic_dh) w1_plp = vplp_w1(iw0)

  if (logic_grad) then
    call lp_arbr_ext_svtv_calcuvalue_g(lri,lrj,nlp_value)
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend
      if (logic_dh) then
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call inn_ext_sv_loop_unpack_g(ilw,irw)
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
            call inn_ext_sv_loop_unpack_g(ilw,irw)
          end do
        end do
      end if
    end do
  else
    call lp_arbr_ext_svtv_calcuvalue_wyb(lri,lrj,nlp_value)
    ilpsta = nstaval(iw0)+1
    ilpend = nstaval(iw0)+nvalue(iw0)
    do iplp=ilpsta,ilpend
      if (logic_dh) then
        ilw = lp_lwei(iplp)
        irw = lp_rwei(iplp)
        call inn_ext_sv_loop_unpack(ilw,irw)
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
            call inn_ext_sv_loop_unpack(ilw,irw)
          end do
        end do
      end if
    end do
  end if
end do

return

end subroutine ar_br_sv_ext_br_ar

subroutine logicg_dd(ilnodesm,irnodesm)

use gugaci_global, only: logic_g49a, logic_g49b, logic_g50
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: ilnodesm, irnodesm

logic_g50 = .false.
logic_g49a = .false.
logic_g49b = .false.
if (ilnodesm < irnodesm) then
  logic_g49a = .true.
else if (ilnodesm == irnodesm) then
  logic_g49a = .true.
  logic_g49b = .true.
  logic_g50 = .true.
else
  logic_g49b = .true.
end if

end subroutine logicg_dd

!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2009, Yubin Wang                                       *
!               2009, Bingbing Suo                                     *
!***********************************************************************
! Sep. 28, 2009 -BSuo- Firstly write by YBWang, revised by BSuo
! Active space partial loops

subroutine link_b1_at_given_orb(mh)       !b^l:lstep<rstep

use gugaci_global, only: istep_occ, iy, iyl, jb, jj_sub, jjl_sub, lp_head, lp_ltail, lp_lwei, lp_rtail, lp_rwei, lpnew_head, &
                         lpnew_ltail, lpnew_lwei, lpnew_rtail, lpnew_rwei, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(inout) :: mh
integer(kind=iwp) :: iactploop, idb, idocc, ilc, ilstep, ind0, irc, irstep, jbr, lphead, lpltail, lplwei, lplwei0, lpnew, &
                     lpnextltail, lpnextrtail, lprtail, lprwei, lprwei0, ni
real(kind=wp) :: w, w0, w1, ww
integer(kind=iwp), parameter :: ismb1(8) = [3,8,34,35,40,44,66,76]

lpnew = 0
do iactploop=1,mh
  lphead = lp_head(iactploop)
  lpltail = lp_ltail(iactploop)
  lprtail = lp_rtail(iactploop)
  lplwei0 = lp_lwei(iactploop)
  lprwei0 = lp_rwei(iactploop)
  w0 = vplp_w0(iactploop)
  w1 = vplp_w1(iactploop)
  idb = jb(lprtail)-jb(lpltail)
  do ilstep=1,4
    ilc = istep_occ(ilstep)
    lpnextltail = jjl_sub(ilstep,lpltail)
    if (lpnextltail == 0) cycle
    do irstep=ilstep+1,4
      irc = istep_occ(irstep)
      idocc = abs(ilc-irc)
      if (idocc /= 1) cycle
      lpnextrtail = jj_sub(irstep,lprtail)
      if (lpnextrtail == 0) cycle
      ind0 = (idb+2)*16+(ilstep-1)*4+irstep
      jbr = jb(lprtail)
      do ni=1,8
        if (ind0 /= ismb1(ni)) cycle
        call segmidb1(w,ww,ni,jbr)
        lpnew = lpnew+1
        lpnew_head(lpnew) = lphead
        lpnew_ltail(lpnew) = lpnextltail
        lpnew_rtail(lpnew) = lpnextrtail
        lplwei = lplwei0
        lprwei = lprwei0
        if (ilstep /= 1) lplwei = lplwei+iyl(ilstep,lpltail)
        if (irstep /= 1) lprwei = lprwei+iy(irstep,lprtail)
        lpnew_lwei(lpnew) = lplwei
        lpnew_rwei(lpnew) = lprwei
        vplpnew_w0(lpnew) = w0*w
        vplpnew_w1(lpnew) = w1*ww
        if ((vplpnew_w0(lpnew) == 0) .and. (vplpnew_w1(lpnew) == 0)) then
          lpnew = lpnew-1
        end if
      end do
    end do
  end do
end do
mh = lpnew
call change_vplp_pointer_arrays()

return

end subroutine link_b1_at_given_orb

subroutine link_b2_at_given_orb(mh)        !b^r:lstep>rstep

use gugaci_global, only: istep_occ, iy, iyl, jb, jj_sub, jjl_sub, logic_br, lp_head, lp_ltail, lp_lwei, lp_rtail, lp_rwei, &
                         lpnew_head, lpnew_ltail, lpnew_lwei, lpnew_rtail, lpnew_rwei, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(inout) :: mh
integer(kind=iwp) :: iactploop, idb, idocc, ilc, ilstep, ind0, irc, irstep, jbr, lphead, lpltail, lplwei, lplwei0, lpnew, &
                     lpnextltail, lpnextrtail, lprtail, lprwei, lprwei0, ni
real(kind=wp) :: w, w0, w1, ww
logical(kind=iwp) :: logic_plbr
integer(kind=iwp), parameter :: ismb2(8) = [5,15,37,41,46,47,73,78]

lpnew = 0
do iactploop=1,mh
  logic_plbr = logic_br(iactploop)
  if (.not. logic_plbr) cycle
  lphead = lp_head(iactploop)
  lpltail = lp_ltail(iactploop)
  lprtail = lp_rtail(iactploop)
  lplwei0 = lp_lwei(iactploop)
  lprwei0 = lp_rwei(iactploop)
  w0 = vplp_w0(iactploop)
  w1 = vplp_w1(iactploop)
  idb = jb(lprtail)-jb(lpltail)
  do ilstep=1,4
    ilc = istep_occ(ilstep)
    lpnextltail = jjl_sub(ilstep,lpltail)
    if (lpnextltail == 0) cycle
    do irstep=1,ilstep-1
      irc = istep_occ(irstep)
      idocc = abs(ilc-irc)
      if (idocc /= 1) cycle
      lpnextrtail = jj_sub(irstep,lprtail)
      if (lpnextrtail == 0) cycle
      ind0 = (idb+2)*16+(ilstep-1)*4+irstep
      jbr = jb(lprtail)
      do ni=1,8
        if (ind0 /= ismb2(ni)) cycle
        call segmidb2(w,ww,ni,jbr)
        lpnew = lpnew+1
        !if (lpnew == 103) write(u6,*) 'bbs_tmp'
        lpnew_head(lpnew) = lphead
        lpnew_ltail(lpnew) = lpnextltail
        lpnew_rtail(lpnew) = lpnextrtail
        lplwei = lplwei0
        lprwei = lprwei0
        if (ilstep /= 1) lplwei = lplwei+iyl(ilstep,lpltail)
        if (irstep /= 1) lprwei = lprwei+iy(irstep,lprtail)
        lpnew_lwei(lpnew) = lplwei
        lpnew_rwei(lpnew) = lprwei
        vplpnew_w0(lpnew) = w0*w
        vplpnew_w1(lpnew) = w1*ww
        if ((vplpnew_w0(lpnew) == 0) .and. (vplpnew_w1(lpnew) == 0)) then
          lpnew = lpnew-1
        end if
      end do
    end do
  end do
end do
mh = lpnew
call change_vplp_pointer_arrays()

return

end subroutine link_b2_at_given_orb

subroutine link_b3_at_given_orb(mh)      !b&l:lstep>rstep

use gugaci_global, only: istep_occ, iy, iyl, jb, jj_sub, jjl_sub, lp_head, lp_ltail, lp_lwei, lp_rtail, lp_rwei, lpnew_head, &
                         lpnew_ltail, lpnew_lwei, lpnew_rtail, lpnew_rwei, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(inout) :: mh
integer(kind=iwp) :: iactploop, idb, idocc, ilc, ilstep, ind0, irc, irstep, jbr, lphead, lpltail, lplwei, lplwei0, lpnew, &
                     lpnextltail, lpnextrtail, lprtail, lprwei, lprwei0, ni
real(kind=wp) :: w, w0, w1, ww
integer(kind=iwp), parameter :: ismb3(8) = [21,25,30,31,53,57,62,63]

lpnew = 0
do iactploop=1,mh
  lphead = lp_head(iactploop)
  lpltail = lp_ltail(iactploop)
  lprtail = lp_rtail(iactploop)
  lplwei0 = lp_lwei(iactploop)
  lprwei0 = lp_rwei(iactploop)
  w0 = vplp_w0(iactploop)
  w1 = vplp_w1(iactploop)
  idb = jb(lprtail)-jb(lpltail)
  do ilstep=1,4
    ilc = istep_occ(ilstep)
    lpnextltail = jjl_sub(ilstep,lpltail)
    if (lpnextltail == 0) cycle
    do irstep=1,ilstep-1
      irc = istep_occ(irstep)
      idocc = abs(ilc-irc)
      if (idocc /= 1) cycle
      lpnextrtail = jj_sub(irstep,lprtail)
      if (lpnextrtail == 0) cycle
      ind0 = (idb+2)*16+(ilstep-1)*4+irstep
      jbr = jb(lprtail)
      do ni=1,8
        if (ind0 /= ismb3(ni)) cycle
        call segmidb3(w,ww,ni,jbr)
        lpnew = lpnew+1
        lpnew_head(lpnew) = lphead
        lpnew_ltail(lpnew) = lpnextltail
        lpnew_rtail(lpnew) = lpnextrtail
        lplwei = lplwei0
        lprwei = lprwei0
        if (ilstep /= 1) lplwei = lplwei+iyl(ilstep,lpltail)
        if (irstep /= 1) lprwei = lprwei+iy(irstep,lprtail)
        lpnew_lwei(lpnew) = lplwei
        lpnew_rwei(lpnew) = lprwei
        vplpnew_w0(lpnew) = w0*w
        vplpnew_w1(lpnew) = w1*ww
        if ((vplpnew_w0(lpnew) == 0) .and. (vplpnew_w1(lpnew) == 0)) then
          lpnew = lpnew-1
        end if
      end do
    end do
  end do
end do
mh = lpnew
call change_vplp_pointer_arrays()

return

end subroutine link_b3_at_given_orb

subroutine link_b4_at_given_orb(mh)         !b&r:lstep<rstep

use gugaci_global, only: istep_occ, iy, iyl, jb, jj_sub, jjl_sub, lp_head, lp_ltail, lp_lwei, lp_rtail, lp_rwei, lpnew_head, &
                         lpnew_ltail, lpnew_lwei, lpnew_rtail, lpnew_rwei, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(inout) :: mh
integer(kind=iwp) :: iactploop, idb, idocc, ilc, ilstep, ind0, irc, irstep, jbr, lphead, lpltail, lplwei, lplwei0, lpnew, &
                     lpnextltail, lpnextrtail, lprtail, lprwei, lprwei0, ni
real(kind=wp) :: w, w0, w1, ww
integer(kind=iwp), parameter :: ismb4(8) = [18,19,24,28,50,51,56,60]

lpnew = 0
do iactploop=1,mh
  lphead = lp_head(iactploop)
  lpltail = lp_ltail(iactploop)
  lprtail = lp_rtail(iactploop)
  lplwei0 = lp_lwei(iactploop)
  lprwei0 = lp_rwei(iactploop)
  w0 = vplp_w0(iactploop)
  w1 = vplp_w1(iactploop)
  idb = jb(lprtail)-jb(lpltail)
  do ilstep=1,4
    ilc = istep_occ(ilstep)
    lpnextltail = jjl_sub(ilstep,lpltail)
    if (lpnextltail == 0) cycle
    do irstep=ilstep+1,4
      irc = istep_occ(irstep)
      idocc = abs(ilc-irc)
      if (idocc /= 1) cycle
      lpnextrtail = jj_sub(irstep,lprtail)
      if (lpnextrtail == 0) cycle
      ind0 = (idb+2)*16+(ilstep-1)*4+irstep
      jbr = jb(lprtail)
      do ni=1,8
        if (ind0 /= ismb4(ni)) cycle
        call segmidb4(w,ww,ni,jbr)
        lpnew = lpnew+1
        lpnew_head(lpnew) = lphead
        lpnew_ltail(lpnew) = lpnextltail
        lpnew_rtail(lpnew) = lpnextrtail
        lplwei = lplwei0
        lprwei = lprwei0
        if (ilstep /= 1) lplwei = lplwei+iyl(ilstep,lpltail)
        if (irstep /= 1) lprwei = lprwei+iy(irstep,lprtail)
        lpnew_lwei(lpnew) = lplwei
        lpnew_rwei(lpnew) = lprwei
        vplpnew_w0(lpnew) = w0*w
        vplpnew_w1(lpnew) = w1*ww
        if ((vplpnew_w0(lpnew) == 0) .and. (vplpnew_w1(lpnew) == 0)) then
          lpnew = lpnew-1
        end if
      end do
    end do
  end do
end do
mh = lpnew
call change_vplp_pointer_arrays()

return

end subroutine link_b4_at_given_orb

subroutine link_d10_at_given_orb(mh)

use gugaci_global, only: iy, iyl, jb, jj_sub, jjl_sub, lp_head, lp_ltail, lp_lwei, lp_rtail, lp_rwei, lpnew_head, lpnew_ltail, &
                         lpnew_lwei, lpnew_rtail, lpnew_rwei, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(inout) :: mh
integer(kind=iwp) :: iactploop, idb, ilstep, ind0, irstep, jbr, lphead, lpltail, lplwei, lplwei0, lpnew, lpnextltail, lpnextrtail, &
                     lprtail, lprwei, lprwei0, ni
real(kind=wp) :: w, w0, w1, ww
integer(kind=iwp), parameter :: ismd10(2) = [29,61]

lpnew = 0
do iactploop=1,mh
  lphead = lp_head(iactploop)
  lpltail = lp_ltail(iactploop)
  lprtail = lp_rtail(iactploop)
  lplwei0 = lp_lwei(iactploop)
  lprwei0 = lp_rwei(iactploop)
  w0 = vplp_w0(iactploop)
  w1 = vplp_w1(iactploop)
  idb = jb(lprtail)-jb(lpltail)
  ilstep = 4
  lpnextltail = jjl_sub(ilstep,lpltail)
  if (lpnextltail == 0) cycle
  irstep = 1
  lpnextrtail = jj_sub(irstep,lprtail)
  if (lpnextrtail == 0) cycle
  ind0 = (idb+2)*16+(ilstep-1)*4+irstep
  jbr = jb(lprtail)
  do ni=1,2
    if (ind0 /= ismd10(ni)) cycle
    call segmidd10(w,ww,ni,jbr)
    lpnew = lpnew+1
    lpnew_head(lpnew) = lphead
    lpnew_ltail(lpnew) = lpnextltail
    lpnew_rtail(lpnew) = lpnextrtail
    lplwei = lplwei0
    lprwei = lprwei0
    if (ilstep /= 1) lplwei = lplwei+iyl(ilstep,lpltail)
    if (irstep /= 1) lprwei = lprwei+iy(irstep,lprtail)
    lpnew_lwei(lpnew) = lplwei
    lpnew_rwei(lpnew) = lprwei
    vplpnew_w0(lpnew) = w0*w
    vplpnew_w1(lpnew) = w1*ww
    if ((vplpnew_w0(lpnew) == 0) .and. (vplpnew_w1(lpnew) == 0)) then
      lpnew = lpnew-1
    end if
  end do
end do
mh = lpnew
call change_vplp_pointer_arrays()

return

end subroutine link_d10_at_given_orb

subroutine link_c2_to_given_orb(mh,lrsta,lrend)

use gugaci_global, only: istep_occ, iy, iyl, jb, jj_sub, jjl_sub, jpad, jpadl, logic_br, logic_newbr, lp_head, lp_ltail, lp_lwei, &
                         lp_rtail, lp_rwei, lpnew_head, lpnew_ltail, lpnew_lwei, lpnew_rtail, lpnew_rwei, norb_act, norb_dz, &
                         vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(inout) :: mh
integer(kind=iwp), intent(in) :: lrsta, lrend
integer(kind=iwp) :: iactploop, idb, ilc, ilstep, ind0, iorb, irc, irstep, jbr, lphead, lpltail, lplwei, lplwei0, lpnew, &
                     lpnextltail, lpnextrtail, lprtail, lprwei, lprwei0, ni
real(kind=wp) :: w, w0, w1, ww
logical(kind=iwp) :: logic_plbr
integer(kind=iwp), parameter :: ismc2(16) = [1,6,7,11,16,33,38,39,42,43,48,65,70,74,75,80]

if ((norb_act == 0) .or. (lrsta == norb_dz+1)) then
  mh = 1
  lp_head(mh) = 0
  lp_ltail(mh) = jpadl
  lp_rtail(mh) = jpad
  lp_lwei(mh) = 0
  lp_rwei(mh) = 0
  vplp_w0(mh) = One
  vplp_w1(mh) = One
  !logic_br(mh) = .true.
end if
if (norb_act == 0) return
do iorb=lrsta,lrend
  lpnew = 0
  do iactploop=1,mh
    lphead = lp_head(iactploop)
    lpltail = lp_ltail(iactploop)
    lprtail = lp_rtail(iactploop)
    lplwei0 = lp_lwei(iactploop)
    lprwei0 = lp_rwei(iactploop)
    w0 = vplp_w0(iactploop)
    w1 = vplp_w1(iactploop)
    logic_plbr = logic_br(iactploop)
    idb = jb(lprtail)-jb(lpltail)
    do ilstep=1,4
      ilc = istep_occ(ilstep)
      lpnextltail = jjl_sub(ilstep,lpltail)
      if (lpnextltail == 0) cycle
      do irstep=1,4
        irc = istep_occ(irstep)
        if (ilc /= irc) cycle
        lpnextrtail = jj_sub(irstep,lprtail)
        if (lpnextrtail == 0) cycle
        jbr = jb(lprtail)
        ind0 = (idb+2)*16+(ilstep-1)*4+irstep
        if ((ilstep == 3) .and. (irstep == 2) .and. (.not. logic_plbr)) cycle
        do ni=1,16
          if (ind0 /= ismc2(ni)) cycle
          call segmidc2(w,ww,ni,jbr)
          lpnew = lpnew+1
          !if (lpnew == 57) write(u6,*) 'bbs_tmp'
          lpnew_head(lpnew) = lphead
          lpnew_ltail(lpnew) = lpnextltail
          lpnew_rtail(lpnew) = lpnextrtail
          lplwei = lplwei0
          lprwei = lprwei0
          if (ilstep /= 1) lplwei = lplwei+iyl(ilstep,lpltail)
          if (irstep /= 1) lprwei = lprwei+iy(irstep,lprtail)
          lpnew_lwei(lpnew) = lplwei
          lpnew_rwei(lpnew) = lprwei
          logic_newbr(lpnew) = logic_plbr
          if (irstep > ilstep) logic_newbr(lpnew) = .true.
          vplpnew_w0(lpnew) = w0*w
          vplpnew_w1(lpnew) = w1*ww
          if ((vplpnew_w0(lpnew) == 0) .and. (vplpnew_w1(lpnew) == 0)) then
            lpnew = lpnew-1
            cycle
          end if
        end do
      end do
    end do
  end do
  mh = lpnew
  call change_vplp_pointer_arrays()
  call change_br_pointer_arrays()
end do

end subroutine link_c2_to_given_orb

subroutine link_c1_to_given_orb(mh,lrsta,lrend)

use gugaci_global, only: istep_occ, iy, iyl, jb, jj_sub, jjl_sub, jpad, jpadl, lp_head, lp_ltail, lp_lwei, lp_rtail, lp_rwei, &
                         lpnew_head, lpnew_ltail, lpnew_lwei, lpnew_rtail, lpnew_rwei, norb_act, norb_dz, vplp_w0, vplp_w1, &
                         vplpnew_w0, vplpnew_w1
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(inout) :: mh
integer(kind=iwp), intent(in) :: lrsta, lrend
integer(kind=iwp) :: iactploop, idb, ilc, ilstep, ind0, iorb, irc, irstep, jbr, lphead, lpltail, lplwei, lplwei0, lpnew, &
                     lpnextltail, lpnextrtail, lprtail, lprwei, lprwei0, ni
real(kind=wp) :: w, w0, w1, ww
integer(kind=iwp), parameter :: ismc1(10) = [17,22,23,27,32,49,54,58,59,64]

if ((norb_act == 0) .or. (lrsta == norb_dz+1)) then
  mh = 1
  lp_head(mh) = 0
  lp_ltail(mh) = jpadl
  lp_rtail(mh) = jpad
  lp_lwei(mh) = 0
  lp_rwei(mh) = 0
  vplp_w0(mh) = One
  vplp_w1(mh) = One
end if
if (norb_act == 0) return
do iorb=lrsta,lrend
  lpnew = 0
  do iactploop=1,mh
    lphead = lp_head(iactploop)
    lpltail = lp_ltail(iactploop)
    lprtail = lp_rtail(iactploop)
    lplwei0 = lp_lwei(iactploop)
    lprwei0 = lp_rwei(iactploop)
    w0 = vplp_w0(iactploop)
    w1 = vplp_w1(iactploop)
    jbr = jb(lprtail)
    idb = jb(lprtail)-jb(lpltail)
    do ilstep=1,4
      ilc = istep_occ(ilstep)
      lpnextltail = jjl_sub(ilstep,lpltail)
      if (lpnextltail == 0) cycle
      do irstep=1,4
        irc = istep_occ(irstep)
        if (ilc /= irc) cycle
        lpnextrtail = jj_sub(irstep,lprtail)
        if (lpnextrtail == 0) cycle
        ind0 = (idb+2)*16+(ilstep-1)*4+irstep
        do ni=1,10
          if (ind0 /= ismc1(ni)) cycle
          call segmidc1(w,ww,ni,jbr)
          lpnew = lpnew+1
          !if (lpnew == 15) write(u6,*) 'bbs_tmp'
          lpnew_head(lpnew) = lphead
          lpnew_ltail(lpnew) = lpnextltail
          lpnew_rtail(lpnew) = lpnextrtail
          lplwei = lplwei0
          lprwei = lprwei0
          if (ilstep /= 1) lplwei = lplwei+iyl(ilstep,lpltail)
          if (irstep /= 1) lprwei = lprwei+iy(irstep,lprtail)
          lpnew_lwei(lpnew) = lplwei
          lpnew_rwei(lpnew) = lprwei
          vplpnew_w0(lpnew) = w0*w
          vplpnew_w1(lpnew) = w1*ww
          if ((vplpnew_w0(lpnew) /= 0) .or. (vplpnew_w1(lpnew) /= 0)) exit
          lpnew = lpnew-1
        end do
      end do
    end do
  end do
  mh = lpnew
  call change_vplp_pointer_arrays()
end do

end subroutine link_c1_to_given_orb

subroutine link_c1_to_given_orb_coe(mh,lrsta,lrend)

use gugaci_global, only: istep_occ, iy, iyl, jb, jj_sub, jjl_sub, jpad, jpadl, lp_coe, lp_head, lp_ltail, lp_lwei, lp_rtail, &
                         lp_rwei, lpnew_coe, lpnew_head, lpnew_ltail, lpnew_lwei, lpnew_rtail, lpnew_rwei, norb_act, norb_dz, &
                         vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(inout) :: mh
integer(kind=iwp), intent(in) :: lrsta, lrend
integer(kind=iwp) :: iactploop, idb, ilc, ilstep, ind0, iorb, irc, irstep, jbl, jbr, lphead, lpltail, lplwei, lplwei0, lpnew, &
                     lpnextltail, lpnextrtail, lprtail, lprwei, lprwei0, lr, ni
real(kind=wp) :: w, w0, w1, ww
integer(kind=iwp), parameter :: ismc1(10) = [17,22,23,27,32,49,54,58,59,64]
integer(kind=iwp), external :: k_coe

if ((norb_act == 0) .or. (lrsta == norb_dz+1)) then
  mh = 1
  lp_head(mh) = 0
  lp_ltail(mh) = jpadl
  lp_rtail(mh) = jpad
  lp_lwei(mh) = 0
  lp_rwei(mh) = 0
  vplp_w0(mh) = One
  vplp_w1(mh) = One
end if
if (norb_act == 0) return
do iorb=lrsta,lrend
  lpnew = 0
  do iactploop=1,mh
    lphead = lp_head(iactploop)
    lpltail = lp_ltail(iactploop)
    lprtail = lp_rtail(iactploop)
    lplwei0 = lp_lwei(iactploop)
    lprwei0 = lp_rwei(iactploop)
    w0 = vplp_w0(iactploop)
    w1 = vplp_w1(iactploop)
    jbl = jb(lpltail)
    jbr = jb(lprtail)
    idb = jbr-jbl
    !iposib = jb(lprtail)*80
    do ilstep=1,4
      ilc = istep_occ(ilstep)
      lpnextltail = jjl_sub(ilstep,lpltail)
      if (lpnextltail == 0) cycle
      do irstep=1,4
        irc = istep_occ(irstep)
        if (ilc /= irc) cycle
        lpnextrtail = jj_sub(irstep,lprtail)
        if (lpnextrtail == 0) cycle
        ind0 = (idb+2)*16+(ilstep-1)*4+irstep
        do ni=1,10
          if (ind0 /= ismc1(ni)) cycle
          call segmidc1(w,ww,ni,jbr)
          lpnew = lpnew+1
          lpnew_head(lpnew) = lphead
          lpnew_ltail(lpnew) = lpnextltail
          lpnew_rtail(lpnew) = lpnextrtail
          lplwei = lplwei0
          lprwei = lprwei0
          if (ilstep /= 1) lplwei = lplwei+iyl(ilstep,lpltail)
          if (irstep /= 1) lprwei = lprwei+iy(irstep,lprtail)
          lpnew_lwei(lpnew) = lplwei
          lpnew_rwei(lpnew) = lprwei
          vplpnew_w0(lpnew) = w0*w
          vplpnew_w1(lpnew) = w1*ww
          if ((vplpnew_w0(lpnew) /= 0) .or. (vplpnew_w1(lpnew) /= 0)) then
            do lr=norb_dz+1,iorb-1
              lpnew_coe(lr,lpnew) = lp_coe(lr,iactploop)
            end do
            lpnew_coe(iorb,lpnew) = k_coe(jbl,jbr,ilstep,irstep)
            exit
          end if
          lpnew = lpnew-1
        end do
      end do
    end do
  end do
  mh = lpnew
  call change_vplp_pointer_arrays()
  call change_coe_pointer_arrays()
end do

end subroutine link_c1_to_given_orb_coe

subroutine head_drr_at_given_orb(mh,lri)

use gugaci_global, only: iy, iyl, jj_sub, jjl_sub, lp_head, lp_ltail, lp_lwei, lp_rtail, lp_rwei, no, vplp_w0, vplp_w1
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: mh
integer(kind=iwp), intent(in) :: lri
integer(kind=iwp) :: iactploop, ind0, jdl, jdr, jl, jp, jpend, jpsta, jr
real(kind=wp) :: w, ww
integer(kind=iwp), parameter :: ishdrr = 36

iactploop = 0
jpsta = no(lri-1)+1
jpend = no(lri)
do jp=jpsta,jpend
  if ((iy(1,jp) == 0) .or. (iyl(1,jp) == 0)) cycle
  jdl = 1
  jl = jjl_sub(jdl,jp)        ! no cycle
  if (jl == 0) cycle
  jdr = 4
  jr = jj_sub(jdr,jp)
  if (jr == 0) cycle
  ind0 = 32+(jdl-1)*4+jdr
  !start '^'
  if (ind0 /= ishdrr) cycle
  call stermhd5(w,ww)
  iactploop = iactploop+1
  lp_head(iactploop) = jp
  lp_ltail(iactploop) = jl
  lp_rtail(iactploop) = jr
  vplp_w0(iactploop) = w
  vplp_w1(iactploop) = ww
  lp_lwei(iactploop) = 0
  lp_rwei(iactploop) = 0
  if (jdl /= 1) lp_lwei(iactploop) = iyl(jdl,jp)
  if (jdr /= 1) lp_rwei(iactploop) = iy(jdr,jp)
end do
mh = iactploop

return

end subroutine head_drr_at_given_orb

subroutine head_drl_at_given_orb(mh,lri)

use gugaci_global, only: istep_occ, iy, iyl, jb, jj_sub, jjl_sub, logic_br, lp_head, lp_ltail, lp_lwei, lp_rtail, lp_rwei, no, &
                         vplp_w0, vplp_w1
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: mh
integer(kind=iwp), intent(in) :: lri
integer(kind=iwp) :: iactploop, ind0, jbr, jdl, jdr, jl, jp, jpend, jpsta, jr, ndl, ndr, ni
real(kind=wp) :: w, ww
integer(kind=iwp), parameter :: ishd1(4) = [38,39,43,48]

jpsta = no(lri-1)+1
jpend = no(lri)
iactploop = 0
do jp=jpsta,jpend
  if ((iy(1,jp) == 0) .or. (iyl(1,jp) == 0)) cycle
  do jdl=2,4
    jl = jjl_sub(jdl,jp)
    if (jl == 0) cycle
    ndl = istep_occ(jdl)
    do jdr=jdl,4
      ndr = istep_occ(jdr)
      if (ndr /= ndl) cycle
      jr = jj_sub(jdr,jp)
      if (jr == 0) cycle
      jbr = jb(jp)
      ind0 = 32+(jdl-1)*4+jdr
      !start '^'
      do ni=1,4
        if (ind0 /= ishd1(ni)) cycle
        call stermhd1(w,ww,ni,jbr)
        iactploop = iactploop+1
        lp_head(iactploop) = jp
        lp_ltail(iactploop) = jl
        lp_rtail(iactploop) = jr
        vplp_w0(iactploop) = w
        vplp_w1(iactploop) = ww
        lp_lwei(iactploop) = 0
        lp_rwei(iactploop) = 0
        if (jdl /= 1) lp_lwei(iactploop) = iyl(jdl,jp)
        if (jdr /= 1) lp_rwei(iactploop) = iy(jdr,jp)
        logic_br(iactploop) = .false.
        if (jdr > jdl) logic_br(iactploop) = .true.
      end do
    end do
  end do
end do
mh = iactploop

return

end subroutine head_drl_at_given_orb

subroutine head_ar_at_given_orb(mh,lri)

use gugaci_global, only: istep_occ, iy, iyl, jb, jj_sub, jjl_sub, lp_coe, lp_head, lp_ltail, lp_lwei, lp_rtail, lp_rwei, no, &
                         vplp_w0, vplp_w1
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: mh
integer(kind=iwp), intent(in) :: lri
integer(kind=iwp) :: iactploop, ind0, jbr, jd_occ, jdl, jdl_occ, jdr, jdr_occ, jl, jp, jpend, jpsta, jr, ni
real(kind=wp) :: w, ww
logical(kind=iwp) :: found
integer(kind=iwp), parameter :: isha4(4) = [34,35,40,44]

iactploop = 0
jpsta = no(lri-1)+1
jpend = no(lri)
do jp=jpsta,jpend
  if ((iy(1,jp) == 0) .or. (iyl(1,jp) == 0)) cycle
  do jdl=1,3
    jl = jjl_sub(jdl,jp)
    if (jl == 0) cycle
    jdl_occ = istep_occ(jdl)
    do jdr=jdl+1,4
      jr = jj_sub(jdr,jp)
      if (jr == 0) cycle
      jdr_occ = istep_occ(jdr)
      jd_occ = jdr_occ-jdl_occ
      if (jd_occ /= 1) cycle
      jbr = jb(jp)
      ind0 = 32+(jdl-1)*4+jdr
      !start '^'
      found = .false.
      do ni=1,4
        if (ind0 == isha4(ni)) then
          call stermha4(w,ww,ni,jbr)
          found = .true.
          exit
        end if
      end do
      if (.not. found) cycle
      iactploop = iactploop+1
      !if (iactploop == 4) write(u6,*) 'bbs_tmp'
      lp_head(iactploop) = jp
      lp_ltail(iactploop) = jl
      lp_rtail(iactploop) = jr
      lp_coe(lri,iactploop) = 0
      if (jdr_occ == 2) lp_coe(lri,iactploop) = 100
      vplp_w0(iactploop) = w
      vplp_w1(iactploop) = ww
      lp_lwei(iactploop) = 0
      lp_rwei(iactploop) = 0
      if (jdl /= 1) lp_lwei(iactploop) = iyl(jdl,jp)
      if (jdr /= 1) lp_rwei(iactploop) = iy(jdr,jp)
    end do
  end do
end do
mh = iactploop

return

end subroutine head_ar_at_given_orb

subroutine tail_ar_at_given_orb(mh)       !a^r:lstep>rstep

use gugaci_global, only: istep_occ, iy, iyl, ja, jb, jj_sub, jjl_sub, jm, lp_head, lp_ltail, lp_lwei, lp_rtail, lp_rwei, &
                         lpnew_head, lpnew_ltail, lpnew_lwei, lpnew_rtail, lpnew_rwei, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(inout) :: mh
integer(kind=iwp) :: iactploop, idb, idocc, ilc, ilstep, ind0, irc, irstep, jbr, lphead, lpltail, lplwei, lplwei0, lpnew, &
                     lpnextltail, lpnextrtail, lprtail, lprwei, lprwei0, ni
real(kind=wp) :: w, w0, w1, ww
integer(kind=iwp), parameter :: isla2(4) = [21,31,57,62]

lpnew = 0
do iactploop=1,mh
  lphead = lp_head(iactploop)
  lpltail = lp_ltail(iactploop)
  lprtail = lp_rtail(iactploop)
  lplwei0 = lp_lwei(iactploop)
  lprwei0 = lp_rwei(iactploop)
  w0 = vplp_w0(iactploop)
  w1 = vplp_w1(iactploop)
  idb = jb(lprtail)-jb(lpltail)
  do ilstep=1,4
    ilc = istep_occ(ilstep)
    lpnextltail = jjl_sub(ilstep,lpltail)
    if (lpnextltail == 0) cycle
    do irstep=1,ilstep-1
      irc = istep_occ(irstep)
      idocc = abs(ilc-irc)
      if (idocc /= 1) cycle
      lpnextrtail = jj_sub(irstep,lprtail)
      if (lpnextrtail == 0) cycle

      if (ja(lpnextltail) /= ja(lpnextrtail)) cycle
      if (jb(lpnextltail) /= jb(lpnextrtail)) cycle
      if (jm(lpnextltail) /= jm(lpnextrtail)) cycle

      ind0 = (idb+2)*16+(ilstep-1)*4+irstep
      jbr = jb(lprtail)
      do ni=1,4
        if (ind0 /= isla2(ni)) cycle
        call stermla2(w,ww,ni,jbr)
        lpnew = lpnew+1
        lpnew_head(lpnew) = lphead
        lpnew_ltail(lpnew) = lpnextltail
        lpnew_rtail(lpnew) = lpnextrtail
        lplwei = lplwei0
        lprwei = lprwei0
        if (ilstep /= 1) lplwei = lplwei+iyl(ilstep,lpltail)
        if (irstep /= 1) lprwei = lprwei+iy(irstep,lprtail)
        lpnew_lwei(lpnew) = lplwei
        lpnew_rwei(lpnew) = lprwei
        vplpnew_w0(lpnew) = w0*w
        vplpnew_w1(lpnew) = w1*ww
        if ((vplpnew_w0(lpnew) == 0) .and. (vplpnew_w1(lpnew) == 0)) then
          lpnew = lpnew-1
        end if
      end do
    end do
  end do
end do
mh = lpnew
!call change_vplp_pointer_arrays()

return

end subroutine tail_ar_at_given_orb

subroutine tail_al_at_given_orb(mh)       !a^l:lstep<rstep

use gugaci_global, only: istep_occ, iy, iyl, ja, jb, jj_sub, jjl_sub, jm, lp_head, lp_ltail, lp_lwei, lp_rtail, lp_rwei, &
                         lpnew_head, lpnew_ltail, lpnew_lwei, lpnew_rtail, lpnew_rwei, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(inout) :: mh
integer(kind=iwp) :: iactploop, idb, idocc, ilc, ilstep, ind0, irc, irstep, jbr, lphead, lpltail, lplwei, lplwei0, lpnew, &
                     lpnextltail, lpnextrtail, lprtail, lprwei, lprwei0, ni
real(kind=wp) :: w, w0, w1, ww
integer(kind=iwp), parameter :: isla1(4) = [19,24,50,60]

lpnew = 0
do iactploop=1,mh
  lphead = lp_head(iactploop)
  lpltail = lp_ltail(iactploop)
  lprtail = lp_rtail(iactploop)
  lplwei0 = lp_lwei(iactploop)
  lprwei0 = lp_rwei(iactploop)
  w0 = vplp_w0(iactploop)
  w1 = vplp_w1(iactploop)
  idb = jb(lprtail)-jb(lpltail)
  do ilstep=1,4
    ilc = istep_occ(ilstep)
    lpnextltail = jjl_sub(ilstep,lpltail)
    if (lpnextltail == 0) cycle
    do irstep=ilstep+1,4
      irc = istep_occ(irstep)
      idocc = abs(ilc-irc)
      if (idocc /= 1) cycle
      lpnextrtail = jj_sub(irstep,lprtail)
      if (lpnextrtail == 0) cycle

      if (ja(lpnextltail) /= ja(lpnextrtail)) cycle
      if (jb(lpnextltail) /= jb(lpnextrtail)) cycle
      if (jm(lpnextltail) /= jm(lpnextrtail)) cycle

      ind0 = (idb+2)*16+(ilstep-1)*4+irstep
      jbr = jb(lprtail)
      do ni=1,4
        if (ind0 /= isla1(ni)) cycle
        call stermla1(w,ww,ni,jbr)
        lpnew = lpnew+1
        lpnew_head(lpnew) = lphead
        lpnew_ltail(lpnew) = lpnextltail
        lpnew_rtail(lpnew) = lpnextrtail
        lplwei = lplwei0
        lprwei = lprwei0
        if (ilstep /= 1) lplwei = lplwei+iyl(ilstep,lpltail)
        if (irstep /= 1) lprwei = lprwei+iy(irstep,lprtail)
        lpnew_lwei(lpnew) = lplwei
        lpnew_rwei(lpnew) = lprwei
        vplpnew_w0(lpnew) = w0*w
        vplpnew_w1(lpnew) = w1*ww
        if ((vplpnew_w0(lpnew) == 0) .and. (vplpnew_w1(lpnew) == 0)) then
          lpnew = lpnew-1
        end if
      end do
    end do
  end do
end do
mh = lpnew
!call change_vplp_pointer_arrays()

return

end subroutine tail_al_at_given_orb

subroutine tail_drr_at_given_orb(mh)       !d^r^r(3,0):locc=

use gugaci_global, only: iy, iyl, ja, jb, jj_sub, jjl_sub, jm, lp_head, lp_ltail, lp_lwei, lp_rtail, lp_rwei, lpnew_head, &
                         lpnew_ltail, lpnew_lwei, lpnew_rtail, lpnew_rwei, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(inout) :: mh
integer(kind=iwp) :: iactploop, idb, ilstep, ind0, irstep, lphead, lpltail, lplwei, lplwei0, lpnew, lpnextltail, lpnextrtail, &
                     lprtail, lprwei, lprwei0
real(kind=wp) :: w, w0, w1, ww
integer(kind=iwp), parameter :: isld6 = 45

lpnew = 0
do iactploop=1,mh
  lphead = lp_head(iactploop)
  lpltail = lp_ltail(iactploop)
  lprtail = lp_rtail(iactploop)
  lplwei0 = lp_lwei(iactploop)
  lprwei0 = lp_rwei(iactploop)
  w0 = vplp_w0(iactploop)
  w1 = vplp_w1(iactploop)
  idb = jb(lprtail)-jb(lpltail)
  ilstep = 4
  lpnextltail = jjl_sub(ilstep,lpltail)
  if (lpnextltail == 0) cycle
  irstep = 1
  lpnextrtail = jj_sub(irstep,lprtail)
  if (lpnextrtail == 0) cycle

  if (ja(lpnextltail) /= ja(lpnextrtail)) cycle
  if (jb(lpnextltail) /= jb(lpnextrtail)) cycle
  if (jm(lpnextltail) /= jm(lpnextrtail)) cycle

  ind0 = (idb+2)*16+(ilstep-1)*4+irstep
  if (ind0 /= isld6) cycle
  call stermld6(w,ww)
  lpnew = lpnew+1
  lpnew_head(lpnew) = lphead
  lpnew_ltail(lpnew) = lpnextltail
  lpnew_rtail(lpnew) = lpnextrtail
  lplwei = lplwei0
  lprwei = lprwei0
  if (ilstep /= 1) lplwei = lplwei+iyl(ilstep,lpltail)
  if (irstep /= 1) lprwei = lprwei+iy(irstep,lprtail)
  lpnew_lwei(lpnew) = lplwei
  lpnew_rwei(lpnew) = lprwei
  vplpnew_w0(lpnew) = w0*w
  vplpnew_w1(lpnew) = w1*ww
  if ((vplpnew_w0(lpnew) == 0) .and. (vplpnew_w1(lpnew) == 0)) then
    lpnew = lpnew-1
  end if
end do
mh = lpnew
!call change_vplp_pointer_arrays()

return

end subroutine tail_drr_at_given_orb

subroutine tail_drl_at_given_orb(mh)       !d^r^l:locc=rocc

use gugaci_global, only: istep_occ, iy, iyl, ja, jb, jj_sub, jjl_sub, jm, logic_br, lp_head, lp_ltail, lp_lwei, lp_rtail, lp_rwei, &
                         lpnew_head, lpnew_ltail, lpnew_lwei, lpnew_rtail, lpnew_rwei, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(inout) :: mh
logical(kind=iwp) :: logic_plbr
integer(kind=iwp) :: iactploop, idb, idocc, ilc, ilstep, ind0, irc, irstart, irstep, jbr, lphead, lpltail, lplwei, lplwei0, lpnew, &
                     lpnextltail, lpnextrtail, lprtail, lprwei, lprwei0, ni
real(kind=wp) :: w, w0, w1, ww
integer(kind=iwp), parameter :: isld2(5) = [7,38,43,48,74]

lpnew = 0
do iactploop=1,mh
  lphead = lp_head(iactploop)
  lpltail = lp_ltail(iactploop)
  lprtail = lp_rtail(iactploop)
  lplwei0 = lp_lwei(iactploop)
  lprwei0 = lp_rwei(iactploop)
  w0 = vplp_w0(iactploop)
  w1 = vplp_w1(iactploop)
  idb = jb(lprtail)-jb(lpltail)
  logic_plbr = logic_br(iactploop)
  do ilstep=1,4
    ilc = istep_occ(ilstep)
    lpnextltail = jjl_sub(ilstep,lpltail)
    if (lpnextltail == 0) cycle
    irstart = 1
    if (.not. logic_plbr) irstart = ilstep+1
    do irstep=irstart,4
      irc = istep_occ(irstep)
      idocc = abs(ilc-irc)
      if (idocc /= 0) cycle
      lpnextrtail = jj_sub(irstep,lprtail)
      if (lpnextrtail == 0) cycle

      if (ja(lpnextltail) /= ja(lpnextrtail)) cycle
      if (jb(lpnextltail) /= jb(lpnextrtail)) cycle
      if (jm(lpnextltail) /= jm(lpnextrtail)) cycle

      ind0 = (idb+2)*16+(ilstep-1)*4+irstep
      jbr = jb(lprtail)
      do ni=1,5
        if (ind0 /= isld2(ni)) cycle
        call stermld2(w,ww,ni,jbr)
        lpnew = lpnew+1
        lpnew_head(lpnew) = lphead
        lpnew_ltail(lpnew) = lpnextltail
        lpnew_rtail(lpnew) = lpnextrtail
        lplwei = lplwei0
        lprwei = lprwei0
        if (ilstep /= 1) lplwei = lplwei+iyl(ilstep,lpltail)
        if (irstep /= 1) lprwei = lprwei+iy(irstep,lprtail)
        lpnew_lwei(lpnew) = lplwei
        lpnew_rwei(lpnew) = lprwei
        vplpnew_w0(lpnew) = w0*w
        vplpnew_w1(lpnew) = w1*ww
        if ((vplpnew_w0(lpnew) == 0) .and. (vplpnew_w1(lpnew) == 0)) then
          lpnew = lpnew-1
        end if
      end do
    end do
  end do
end do
mh = lpnew
!call change_vplp_pointer_arrays()

return

end subroutine tail_drl_at_given_orb

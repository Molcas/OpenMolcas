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

subroutine link_b1_at_given_orb(mh,lract)       !b^l:lstep<rstep

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
dimension ismb1(8)
data ismb1/3,8,34,35,40,44,66,76/

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
! Avoid unused argument warnings
if (.false.) call Unused_integer(lract)

end subroutine link_b1_at_given_orb

subroutine link_b2_at_given_orb(mh,lract)        !b^r:lstep>rstep

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
dimension ismb2(8)
data ismb2/5,15,37,41,46,47,73,78/
logical logic_plbr

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
        !if (lpnew == 103) write(6,*) 'bbs_tmp'
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
! Avoid unused argument warnings
if (.false.) call Unused_integer(lract)

end subroutine link_b2_at_given_orb

subroutine link_b3_at_given_orb(mh,lract)      !b&l:lstep>rstep

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
dimension ismb3(8)
data ismb3/21,25,30,31,53,57,62,63/

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
! Avoid unused argument warnings
if (.false.) call Unused_integer(lract)

end subroutine link_b3_at_given_orb

subroutine link_b4_at_given_orb(mh,lract)         !b&r:lstep<rstep

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
dimension ismb4(8)
data ismb4/18,19,24,28,50,51,56,60/

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
! Avoid unused argument warnings
if (.false.) call Unused_integer(lract)

end subroutine link_b4_at_given_orb

subroutine link_d10_at_given_orb(mh,lract)

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
dimension ismd10(2)
data ismd10/29,61/

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
! Avoid unused argument warnings
if (.false.) call Unused_integer(lract)

end subroutine link_d10_at_given_orb

subroutine link_c2_to_given_orb(mh,lrsta,lrend)

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
dimension ismc2(16)
logical logic_plbr
data ismc2/1,6,7,11,16,33,38,39,42,43,48,65,70,74,75,80/

if ((norb_act == 0) .or. (lrsta == norb_dz+1)) then
  mh = 1
  lp_head(mh) = 0
  lp_ltail(mh) = jpadl
  lp_rtail(mh) = jpad
  lp_lwei(mh) = 0
  lp_rwei(mh) = 0
  vplp_w0(mh) = 1.0d0
  vplp_w1(mh) = 1.0d0
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
          !if (lpnew == 57) write(6,*) 'bbs_tmp'
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

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
dimension ismc1(10)
data ismc1/17,22,23,27,32,49,54,58,59,64/

if ((norb_act == 0) .or. (lrsta == norb_dz+1)) then
  mh = 1
  lp_head(mh) = 0
  lp_ltail(mh) = jpadl
  lp_rtail(mh) = jpad
  lp_lwei(mh) = 0
  lp_rwei(mh) = 0
  vplp_w0(mh) = 1.0d0
  vplp_w1(mh) = 1.0d0
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
          !if (lpnew == 15) write(6,*) 'bbs_tmp'
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

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
dimension ismc1(10)
data ismc1/17,22,23,27,32,49,54,58,59,64/

if ((norb_act == 0) .or. (lrsta == norb_dz+1)) then
  mh = 1
  lp_head(mh) = 0
  lp_ltail(mh) = jpadl
  lp_rtail(mh) = jpad
  lp_lwei(mh) = 0
  lp_rwei(mh) = 0
  vplp_w0(mh) = 1.0d0
  vplp_w1(mh) = 1.0d0
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

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
data ishdrr/36/

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

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
dimension ishd1(4)
data ishd1/38,39,43,48/

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

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
dimension isha4(4)
logical found
data isha4/34,35,40,44/

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
      !if (iactploop == 4) write(6,*) 'bbs_tmp'
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

subroutine tail_ar_at_given_orb(mh,lract)       !a^r:lstep>rstep

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
dimension isla2(4)
data isla2/21,31,57,62/

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
! Avoid unused argument warnings
if (.false.) call Unused_integer(lract)

end subroutine tail_ar_at_given_orb

subroutine tail_al_at_given_orb(mh,lract)       !a^l:lstep<rstep

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
dimension isla1(4)
data isla1/19,24,50,60/

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
! Avoid unused argument warnings
if (.false.) call Unused_integer(lract)

end subroutine tail_al_at_given_orb

subroutine tail_drr_at_given_orb(mh,lract)       !d^r^r(3,0):locc=

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
data isld6/45/

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
! Avoid unused argument warnings
if (.false.) call Unused_integer(lract)

end subroutine tail_drr_at_given_orb

subroutine tail_drl_at_given_orb(mh,lract)       !d^r^l:locc=rocc

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
dimension isld2(5)
data isld2/7,38,43,48,74/
logical logic_plbr

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
! Avoid unused argument warnings
if (.false.) call Unused_integer(lract)

end subroutine tail_drl_at_given_orb

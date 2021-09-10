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

! inner space loop calculation,
! loops in dbl space, dbl-act space, act space
!***********************************************************************
!   ar      =1 (+a^r)     drr     =2 (+d^rr)   drl     =3 (+d^rl)
!   arbr    =4 (+d^rr)    arbl    =5 (+d^rl)   ard_l^r =6 (+a^l)
!   drrb^r  =7 (+a^r)     drlb^r  =8 (+a^l)    drlb^l  =9 (+a^r)
!   arbrb^r =10 (+a^r)    arblb^r =11 (+a^l)   arblb^l =12 (+a^r)
!   drl        =13 (*)
!***********************************************************************

subroutine inner_space_loop()

implicit real*8(a-h,o-z)

!wsc0 = c_time()
call dbl_space_loop()
!wsc1 = c_time()
call act_space_cloop()
call act_space_ploop()
!wsc2 = c_time()
!write(6,'(2x,2(a5,f12.4))') 'dbl',wsc1-wsc0,'act',wsc2-wsc1

return

end subroutine inner_space_loop

subroutine act_space_cloop()        ! one sub_drt

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"

if (norb_act == 0) return
do ipae_=1,25
  ipae = ipae_ ! ipae is in common block, is this necessary?
  jpae = nu_ae(ipae)
  if (jpae == 0) cycle
  do jpad_=1,mxnode
    jpad = jpad_ ! jpad is in common block, is this necessary?
    if (nu_ad(jpad) == 0) cycle
    call seg_drt()
    if (ndim == 0) cycle
    call copy_to_drtl()
    call cloop_in_act()
  end do
end do

return

end subroutine act_space_cloop

subroutine act_space_ploop()        ! two sub_drt with same ipae

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"

if (norb_act == 0) return
do ipae_=1,25
  ipae = ipae_ ! ipae is in common block, is this necessary?
  jpae = nu_ae(ipae)
  if (jpae == 0) cycle
  do jpadl_=1,mxnode                               ! jpadl
    jpadl = jpadl_ ! jpadl is in common block, is this necessary?
    if (nu_ad(jpadl) == 0) cycle
    jpad = jpadl
    call seg_drt()
    if (ndim == 0) cycle
    call copy_to_drtl()
    do jpad_=1,mxnode                               !jpadr
      jpad = jpad_ ! jpad is in common block, is this necessary?
      if (nu_ad(jpad) == 0) cycle
      call seg_drt()
      if (ndim == 0) cycle
      !if ((ipae == 18) .and. (jpadl == 2) .and. (jpad == 1)) write(6,*)
      call ploop_in_act()
    end do
  end do
end do

return

end subroutine act_space_ploop

subroutine cloop_in_act()

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

do lrai=norb_dz+1,norb_inn-1
  lmi = lsm_inn(lrai)
  do lraj=lrai+1,norb_inn
    lmj = lsm_inn(lraj)
    lmij = mul_tab(lmi,lmj)
    !-------------------------------------------------------------------
    ! line=8 d&r&r--d^r^r
    call head_drr_at_given_orb(mh,lrai)
    logic_br(1:mh) = .true.
    call link_c2_to_given_orb(mh,lrai+1,lraj-1)
    call tail_drr_at_given_orb(mh,lraj)
    !write(6,'(6i6)') 8,mh,lrai,lraj,0,0
    if (mh /= 0) call act_cloop(8,mh,lrai,lraj,0,0)
    !-------------------------------------------------------------------
    ! line=9 d&r&l--d^r^l
    call head_drl_at_given_orb(mh,lrai)
    call link_c2_to_given_orb(mh,lrai+1,lraj-1)
    call tail_drl_at_given_orb(mh,lraj)
    !write(6,'(6i6)') 9,mh,lrai,lraj,0,0
    if (mh /= 0) call act_cloop(9,mh,lrai,lraj,0,0)
    lsmij = mul_tab(lmi,lmj)
    !-------------------------------------------------------------------
    if (lsmij /= 1) goto 400
    !-------------------------------------------------------------------
    ! line=1 a&r--a^r
    call head_ar_at_given_orb(mh,lrai)
    call link_c1_to_given_orb_coe(mh,lrai+1,lraj-1)
    call tail_ar_at_given_orb_coe(mh,lraj)
!          write(6,'(6i6)')1,mh,lrai,lraj,0,0
    if (mh /= 0) call act_cloop(1,mh,lrai,lraj,0,0)
    !call save_clp(1,mh,lra,0)
    !-------------------------------------------------------------------
    ! line=2 a&r-d^r&l-a^l
    do lrak=lrai+1,lraj-1
      call head_ar_at_given_orb(mh,lrai)
      call link_c1_to_given_orb(mh,lrai+1,lrak-1)
      call link_d10_at_given_orb(mh,lrak)
      call link_c1_to_given_orb(mh,lrak+1,lraj-1)
      call tail_al_at_given_orb(mh,lraj)
      !write(6,'(6i6)') 2,mh,lrai,lraj,lrak,0
      if (mh /= 0) call act_cloop(2,mh,lrai,lraj,lrak,0)
    end do
      !-----------------------------------------------------------------
      ! line=3 a&r-b&r-d^r^r
    do lrak=lraj+1,norb_inn
      call head_ar_at_given_orb(mh,lrai)
      call link_c1_to_given_orb(mh,lrai+1,lraj-1)
      call link_b4_at_given_orb(mh,lraj)
      logic_br(1:mh) = .true.
      call link_c2_to_given_orb(mh,lraj+1,lrak-1)
      call tail_drr_at_given_orb(mh,lrak)
      !write(6,'(6i6)') 3,mh,lrai,lraj,lrak,0
      if (mh /= 0) call act_cloop(3,mh,lrai,lraj,lrak,0)
      !-----------------------------------------------------------------
      ! line=5 a&r-b&l-d^r^l
      call head_ar_at_given_orb(mh,lrai)
      call link_c1_to_given_orb(mh,lrai+1,lraj-1)
      call link_b3_at_given_orb(mh,lraj)
      logic_br(1:mh) = .true.
      call link_c2_to_given_orb(mh,lraj+1,lrak-1)
      call tail_drl_at_given_orb(mh,lrak)
      !write(6,'(6i6)') 5,mh,lrai,lraj,lrak,0
      if (mh /= 0) call act_cloop(5,mh,lrai,lraj,lrak,0)
    end do
    !-------------------------------------------------------------------
    do lrak=norb_dz+1,lrai-1
      ! line=10 d&rr--b^r--a^r
      call head_drr_at_given_orb(mh,lrak)
      logic_br(1:mh) = .true.
      call link_c2_to_given_orb(mh,lrak+1,lrai-1)
      call link_b2_at_given_orb(mh,lrai)
      call link_c1_to_given_orb(mh,lrai+1,lraj-1)
      call tail_ar_at_given_orb(mh,lraj)
      !write(6,'(6i6)') 10,mh,lrai,lraj,lrak,0
      if (mh /= 0) call act_cloop(10,mh,lrai,lraj,lrak,0)
      !-----------------------------------------------------------------
      ! line=11 d&r&l-b^r-a^l
      call head_drl_at_given_orb(mh,lrak)
      call link_c2_to_given_orb(mh,lrak+1,lrai-1)
      call link_b2_at_given_orb(mh,lra)
      call link_c1_to_given_orb(mh,lrai+1,lraj-1)
      call tail_al_at_given_orb(mh,lraj)
      !write(6,'(6i6)') 11,mh,lrai,lraj,lrak,0
      if (mh /= 0) call act_cloop(11,mh,lrai,lraj,lrak,0)
      !-----------------------------------------------------------------
      ! line=12 d&r&l-b^l-a^r
      call head_drl_at_given_orb(mh,lrak)
      call link_c2_to_given_orb(mh,lrak+1,lrai-1)
      call link_b1_at_given_orb(mh,lrai)
      call link_c1_to_given_orb(mh,lrai+1,lraj-1)
      call tail_ar_at_given_orb(mh,lraj)
      !write(6,'(6i6)') 12,mh,lrai,lraj,lrak,0
      if (mh /= 0) call act_cloop(12,mh,lrai,lraj,lrak,0)
    end do
400 continue
    if (lraj > norb_inn-2) cycle
    do lrak=lraj+1,norb_inn
      lmk = lsm_inn(lrak)
      lmk = mul_tab(lmij,lmk)
      do lral=lrak+1,norb_inn
        lml = lsm_inn(lral)
        lml = mul_tab(lmk,lml)
        if (lml /= 1) cycle
        ! line=4  a&r--b&r--b^r--a^r
        call head_ar_at_given_orb(mh,lrai)
        call link_c1_to_given_orb(mh,lrai+1,lraj-1)
        call link_b4_at_given_orb(mh,lraj)
        logic_br(1:mh) = .true.
        call link_c2_to_given_orb(mh,lraj+1,lrak-1)
        call link_b2_at_given_orb(mh,lrak)
        call link_c1_to_given_orb(mh,lrak+1,lral-1)
        call tail_ar_at_given_orb(mh,lral)
        !write(6,'(6i6)') 4,mh,lrai,lral,lraj,lrak
        if (mh /= 0) call act_cloop(4,mh,lrai,lral,lraj,lrak)
        !---------------------------------------------------------------
        ! line=6  a&r--b&l--b^r--a^l
        call head_ar_at_given_orb(mh,lrai)
        call link_c1_to_given_orb(mh,lrai+1,lraj-1)
        call link_b3_at_given_orb(mh,lraj)
        logic_br(1:mh) = .true.
        call link_c2_to_given_orb(mh,lraj+1,lrak-1)
        call link_b2_at_given_orb(mh,lrak)
        call link_c1_to_given_orb(mh,lrak+1,lral-1)
        call tail_al_at_given_orb(mh,lral)
        !write(6,'(6i6)') 6,mh,lrai,lral,lraj,lrak
        if (mh /= 0) call act_cloop(6,mh,lrai,lral,lraj,lrak)
        !---------------------------------------------------------------
        ! line=7 a&r--b&l--b^l--a^r
        call head_ar_at_given_orb(mh,lrai)
        call link_c1_to_given_orb(mh,lrai+1,lraj-1)
        call link_b3_at_given_orb(mh,lraj)
        logic_br(1:mh) = .true.
        call link_c2_to_given_orb(mh,lraj+1,lrak-1)
        call link_b1_at_given_orb(mh,lrak)
        call link_c1_to_given_orb(mh,lrak+1,lral-1)
        call tail_ar_at_given_orb(mh,lral)
        !write(6,'(6i6)') 7,mh,lrai,lral,lraj,lrak
        if (mh /= 0) call act_cloop(7,mh,lrai,lral,lraj,lrak)
        !---------------------------------------------------------------
      end do
    end do
  end do
end do

return

end subroutine cloop_in_act

subroutine ploop_in_act()

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

!=======================================================================
do lrai=norb_dz+1,norb_inn
  ! line=25 -c"-d^r^r
  logic_br(1) = .true.
  call link_c2_to_given_orb(mh,norb_dz+1,lrai-1)
  call tail_drr_at_given_orb(mh,lrai)
  if (mh /= 0) call lp_act_tail(25,mh,0,lrai)
  !---------------------------------------------------------------------
  ! line=26 -c"-d^r^l
  logic_br(1) = .true.
  call link_c2_to_given_orb(mh,norb_dz+1,lrai-1)
  call tail_drl_at_given_orb(mh,lrai)
  if (mh /= 0) call lp_act_tail(26,mh,0,lrai)
  !=====================================================================
  ! line=23 -c'-a^l
  call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
  call tail_al_at_given_orb(mh,lrai)
  if (mh /= 0) call lp_act_tail(23,mh,0,0)
  !---------------------------------------------------------------------
  ! line=24 -c'-a^r
  call link_c1_to_given_orb_coe(mh,norb_dz+1,lrai-1)
  call tail_ar_at_given_orb_coe(mh,lrai)
  if (mh /= 0) call lp_act_tail(24,mh,0,0)
  !=====================================================================
  do lrak=lrai+1,norb_inn
    !-------------------------------------------------------------------
    ! line=30 -c'-b&r-d^r^r
    call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
    call link_b4_at_given_orb(mh,lrai)
    logic_br(1:mh) = .true.
    call link_c2_to_given_orb(mh,lrai+1,lrak-1)
    call tail_drr_at_given_orb(mh,lrak)
    if (mh /= 0) call lp_act_tail(30,mh,lrai,0)
    !-------------------------------------------------------------------
    ! line=31 -c'-b&l-d^r^l
    call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
    call link_b3_at_given_orb(mh,lrai)
    logic_br(1:mh) = .true.
    call link_c2_to_given_orb(mh,lrai+1,lrak-1)
    call tail_drl_at_given_orb(mh,lrak)
    if (mh /= 0) call lp_act_tail(31,mh,lrai,0)
  end do
  !=====================================================================
  do lrak=norb_dz+1,lrai-1
    !-------------------------------------------------------------------
    ! line=35 -c'-d^r&l-a^l
    call link_c1_to_given_orb(mh,norb_dz+1,lrak-1)
    call link_d10_at_given_orb(mh,lrak)
    call link_c1_to_given_orb(mh,lrak+1,lrai-1)
    call tail_al_at_given_orb(mh,lrai)
    if (mh /= 0) call lp_act_tail(35,mh,lrak,lrak)
  end do
  !=====================================================================
  do lraj=lrai+1,norb_inn
    !-------------------------------------------------------------------
    ! line=27 -c"-b^r-a^r
    call link_c2_to_given_orb(mh,norb_dz+1,lrai-1)
    call link_b2_at_given_orb(mh,lrai)
    call link_c1_to_given_orb(mh,lrai+1,lraj-1)
    call tail_ar_at_given_orb(mh,lraj)
    if (mh /= 0) call lp_act_tail(27,mh,0,lrai)
    !-------------------------------------------------------------------
    ! line=28 -c"-b^l-a^r
    call link_c2_to_given_orb(mh,norb_dz+1,lrai-1)
    call link_b1_at_given_orb(mh,lrai)
    call link_c1_to_given_orb(mh,lrai+1,lraj-1)
    call tail_ar_at_given_orb(mh,lraj)
    if (mh /= 0) call lp_act_tail(28,mh,0,lrai)
    !-------------------------------------------------------------------
    ! line=29 -c"-b^r-a^l
    call link_c2_to_given_orb(mh,norb_dz+1,lrai-1)
    call link_b2_at_given_orb(mh,lrai)
    call link_c1_to_given_orb(mh,lrai+1,lraj-1)
    call tail_al_at_given_orb(mh,lraj)
    if (mh /= 0) call lp_act_tail(29,mh,0,lrai)
    !-------------------------------------------------------------------
    do lral=lraj+1,norb_inn
      !-----------------------------------------------------------------
      ! line=32 -c'-b&r-c"-b^r-a^r
      call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
      call link_b4_at_given_orb(mh,lrai)
      logic_br(1:mh) = .true.
      call link_c2_to_given_orb(mh,lrai+1,lraj-1)
      call link_b2_at_given_orb(mh,lraj)
      call link_c1_to_given_orb(mh,lraj+1,lral-1)
      call tail_ar_at_given_orb(mh,lral)
      if (mh /= 0) call lp_act_tail(32,mh,lrai,lraj)
      !-----------------------------------------------------------------
      ! line=33 -c'-b&l-c"-b^r-a^l
      call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
      call link_b3_at_given_orb(mh,lrai)
      logic_br(1:mh) = .true.
      call link_c2_to_given_orb(mh,lrai+1,lraj-1)
      call link_b2_at_given_orb(mh,lraj)
      call link_c1_to_given_orb(mh,lraj+1,lral-1)
      call tail_al_at_given_orb(mh,lral)
      if (mh /= 0) call lp_act_tail(33,mh,lrai,lraj)
      !-----------------------------------------------------------------
      ! line=34 -c'-b&l-c"-b^l-a^r
      call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
      call link_b3_at_given_orb(mh,lrai)
      logic_br(1:mh) = .true.
      call link_c2_to_given_orb(mh,lrai+1,lraj-1)
      call link_b1_at_given_orb(mh,lraj)
      call link_c1_to_given_orb(mh,lraj+1,lral-1)
      call tail_ar_at_given_orb(mh,lral)
      if (mh /= 0) call lp_act_tail(34,mh,lrai,lraj)
    end do
  end do
  !---------------------------------------------------------------------
end do

return

end subroutine ploop_in_act

subroutine lp_act_tail(lin,mh,lrg0,lrs0)

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
dimension lpcoe(norb_dz+1:norb_inn)
#include "onepl.fh"

line = lin
lrg = lrg0
lrs = lrs0
jph = 0
do mhlp_=1,mh
  mhlp = mhlp_ ! mhlp is in common block, is this necessary?
  jpel = lpnew_ltail(mhlp)
  jper = lpnew_rtail(mhlp)
  jwl = lpnew_lwei(mhlp)
  jwr = lpnew_rwei(mhlp)
  w0 = vplpnew_w0(mhlp)
  w1 = vplpnew_w1(mhlp)
  if (line == 24) then
    do iorb=norb_dz+1,norb_inn
      lpcoe(iorb) = lpnew_coe(iorb,mhlp)
    end do
  end if
  call dbl_head_act_tail(lpcoe)
end do

return

end subroutine lp_act_tail

subroutine tail_ar_at_given_orb_coe(mh,lract)       !a^r:lstep>rst

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
dimension isla2(4)
data isla2/21,31,57,62/

iorb = lract
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
        do lr=norb_dz+1,iorb-1
          lpnew_coe(lr,lpnew) = lp_coe(lr,iactploop)
        end do
        kcoe = 0
        if (ilstep == 4) kcoe = 100
        lpnew_coe(iorb,lpnew) = kcoe
      end do
    end do
  end do
end do

mh = lpnew
!call change_vplp_pointer_arrays()

end subroutine tail_ar_at_given_orb_coe

subroutine act_cloop(lin,mh,lr0,lr,lrg0,lrs0)

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
dimension lpcoe(norb_dz+1:norb_inn)
#include "onepl.fh"

line = lin
lrg = lrg0
lrs = lrs0
do mhlp_=1,mh
  mhlp = mhlp_ ! mhlp is in common block, is this necessary?
  jph = lpnew_head(mhlp)
  jpel = lpnew_ltail(mhlp)
  jper = lpnew_rtail(mhlp)
  jwl = lpnew_lwei(mhlp)
  jwr = lpnew_rwei(mhlp)
  vlop0 = vplpnew_w0(mhlp)
  vlop1 = vplpnew_w1(mhlp)

  goto(31,32,21,13,22,11,12,51,52,41,42,43),line
  !---------------------------------------------------------------------
31 continue
  do l=norb_dz+1,lr
    lpcoe(l) = lpnew_coe(l,mhlp)
  end do
  wl = voint(lr0,lr)
  do l=lr0,lr
    list = list3(lr0,lr,l)
    kcoe = lpcoe(l)
    call neoc(kcoe,nocc,tcoe)
    wl = wl+nocc*(vint_ci(list+1)+tcoe*vint_ci(list))
  end do
  wl = wl*vlop0
  goto 500
32 continue
  list = list3(lr0,lr,lrg)
  wl = vlop0*vint_ci(list)
  goto 500
  !---------------------------------------------------------------------
11 continue
  list = list4(lr0,lrg,lrs,lr)
  wl = vlop0*(vint_ci(list)-2*vint_ci(list+1))-vlop1*vint_ci(list)
  goto 500
  !---------------------------------------------------------------------
12 continue
  list = list4(lr0,lrg,lrs,lr)
  wl = vlop0*(vint_ci(list+2)-2.0d0*vint_ci(list+1))-vlop1*vint_ci(list+2)
  goto 500
  !---------------------------------------------------------------------
13 continue
  list = list4(lr0,lrg,lrs,lr)
  wl = vlop0*(vint_ci(list+2)+vint_ci(list))-vlop1*(vint_ci(list+2)-vint_ci(list))
  goto 500
  !---------------------------------------------------------------------
21 continue
  !list = list3(lr0,lrg,lrg)
  list = list3(lr0,lr,lrg)
  wl = (vlop0+vlop1)*vint_ci(list)
  goto 500
  !---------------------------------------------------------------------
22 continue
  !list = list3(lr0,lrg,lr)
  list = list3(lr0,lr,lrg)
  wl = vlop0*(vint_ci(list)-2.0d0*vint_ci(list+1))-vlop1*vint_ci(list)
  goto 500
  !---------------------------------------------------------------------
51 continue
  wl = vlop0*voint(lr,lr0)*0.5d0
  goto 500
  !---------------------------------------------------------------------
52 continue
  wl = (vlop0-vlop1)*voint(lr,lr0)
  goto 500
  !---------------------------------------------------------------------
41 continue
  !list = list3(lrg,lr,lr0)
  list = list3(lr0,lr,lrg)
  wl = (vlop0+vlop1)*vint_ci(list)
  goto 500
  !---------------------------------------------------------------------
42 continue
  !list = list3(lrg,lr,lr0)
  list = list3(lr0,lr,lrg)
  wl = (vlop0-vlop1)*vint_ci(list)
  goto 500
  !---------------------------------------------------------------------
43 continue
  !list = list3(lrg,lr,lr0)
  list = list3(lr0,lr,lrg)
  wl = vlop0*(vint_ci(list)-2.0d0*vint_ci(list+1))-vlop1*vint_ci(list)
  !---------------------------------------------------------------------
500 continue
  call prodab(2,jph,jpel,jwl,jwr,0,wl,jper)
end do

return

end subroutine act_cloop

subroutine dbl_head_act_tail_0(lpcoe)

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
dimension lpcoe(norb_dz+1:norb_inn)
#include "onepl.fh"

lra = kk(jpel)-1
jml = mod((jpadl-1),8)
jmr = mod((jpad-1),8)
itypadl = (jpadl-1)/8+2
itypadr = (jpad-1)/8+2
if (jml == 0) then
  jml = 8
  itypadl = itypadl-1
end if
if (jmr == 0) then
  jmr = 8
  itypadr = itypadr-1
end if
if (jpadl == 1) itypadl = 1
if (jpad == 1) itypadr = 1
if (jpadl == 1) jml = ns_sm
if (jpad == 1) jmr = ns_sm
jml = mul_tab(jml,ns_sm)
jmr = mul_tab(jmr,ns_sm)
jmlr = mul_tab(jml,jmr)
lpok = map_jplr(itypadl,itypadr)
if (lpok == 0) return
! 23   24 25 26 27 28 29 30 31 32 33 34 35
goto(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300),line-22
!line=23:-a^l<-->ds(7),dds(9),dt(14),ddtt(16)
100 continue
goto(10,10,10,10,10,10,107,10,109,10,10,10,10,114,10,116,10,10,10,10,10,10,10,10,10,10),lpok
! ds(7-1) ar(23)-drl(30)-
107 continue
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (jmr /= 1) goto 106
  iwdr = just(lri,lri)
  do lrd=norb_frz+1,lri-1
    lmd = lsm_inn(lrd)
    if (lmd /= jml) cycle
    iwdl = jud(lrd)
    w0ds1 = w0_ds(1)
    ni = mod(norb_dz-lrd,2)
    if (ni == 0) w0ds1 = -w0ds1
    vlop0 = w0*w0ds1
    list = list3(lrd,lra,lri)
    wl = vlop0*vint_ci(list)          !   3.2
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end do
106 continue
  do lrj=lri+1,norb_dz
    ! ds(7-3) ar(23)-bl(32)-br(31)-
    lmj = lsm_inn(lrj)
    lmij = mul_tab(lmi,lmj)
    if (lmij /= jmr) cycle
    do lrd=norb_frz+1,lri-1
      iwdr = just(lri,lrj)
      lmd = lsm_inn(lrd)
      if (lmd /= jml) cycle
      iwdl = jud(lrd)
      w0ds3 = w0_ds(3)
      w1ds3 = w1_ds(3)
      ni = mod(norb_dz-lrj+lri-lrd,2)
      if (ni == 0) w0ds3 = -w0ds3
      if (ni == 0) w1ds3 = -w1ds3
      vlop0 = w0*w0ds3
      vlop1 = w1*w1ds3
      list = list4(lrd,lri,lrj,lra)
      wl = (vlop0-vlop1)*vint_ci(list)-2*vlop0*vint_ci(list+1)            !1.1
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
      if (jb_sys > 0) then
        ! ds(7-2) ar(23)-bl(31)-br(32)-         the symmetry problem
        iwdr = just(lrj,lri)
        w0ds3 = w0_ds(2)
        w1ds3 = w1_ds(2)
        ni = mod(norb_dz-lrj+lri-lrd,2)
        if (ni == 0) w0ds3 = -w0ds3
        if (ni == 0) w1ds3 = -w1ds3
        vlop0 = w0*w0ds3
        vlop1 = w1*w1ds3
        list = list4(lrd,lri,lrj,lra)
        wl = (vlop0-vlop1)*vint_ci(list)-2*vlop0*vint_ci(list+1)            !1.1
        call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
      end if
    end do
  end do
end do
return

109 continue
do lri=norb_frz+1,norb_dz
  ! d1s(9-1) ar(13)-drl(30)-
  lmi = lsm_inn(lri)
  do lrd=norb_frz+1,lri-1
    lmd = lsm_inn(lrd)
    if ((lmd == jml) .and. (jmr == 1)) then
      iwdr = just(lri,lri)
      iwdl = jud(lrd)
      w0ds1 = w0_d1s(1)
      ni = mod(norb_dz-lri+lri-lrd,2)
      if (ni == 0) w0ds1 = -w0ds1
      vlop0 = w0*w0ds1
      list = list3(lrd,lra,lri)
      wl = vlop0*vint_ci(list)          !   3.2
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end if
    ! d1s(9-4) drl(12)-br(31)-
    if ((jml == lmd) .and. (jmr == mul_tab(lmd,lmi))) then
      iwdr = just(lrd,lri)
      iwdl = jud(lrd)
      w1ds = w1_d1s(4)
      if (mod(norb_dz-lri,2) == 1) w1ds = -w1ds
      vlop1 = w1*w1ds
      list = list3(lri,lra,lrd)
      wl = -vlop1*vint_ci(list)
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end if
  end do
end do
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    ! d1s(9-3) ar(13)-bl(32)-br(31)-
    lmj = lsm_inn(lrj)
    lmij = mul_tab(lmi,lmj)
    if (lmij /= jmr) cycle
    do lrd=norb_frz+1,lri-1
      iwdr = just(lri,lrj)
      lmd = lsm_inn(lrd)
      if (lmd /= jml) cycle
      iwdl = jud(lrd)
      w0ds3 = w0_d1s(3)
      w1ds3 = w1_d1s(3)
      ni = mod(norb_dz-lrj+lri-lrd,2)
      if (ni == 0) w0ds3 = -w0ds3
      if (ni == 0) w1ds3 = -w1ds3
      vlop0 = w0*w0ds3
      vlop1 = w1*w1ds3
      list = list4(lrd,lri,lrj,lra)
      wl = (vlop0-vlop1)*vint_ci(list)-2*vlop0*vint_ci(list+1)            !1.1
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
      if (jb_sys > 0) then
        ! d1s(9-2)   ar(13)-bl(31)-br(32)-   the symmetry problem
        iwdr = just(lrj,lri)
        w0ds3 = w0_d1s(2)
        w1ds3 = w1_d1s(2)
        ni = mod(norb_dz-lrj+lri-lrd,2)
        if (ni == 0) w0ds3 = -w0ds3
        if (ni == 0) w1ds3 = -w1ds3
        vlop0 = w0*w0ds3
        vlop1 = w1*w1ds3
        list = list4(lrd,lri,lrj,lra)
        wl = (vlop0-vlop1)*vint_ci(list)-2*vlop0*vint_ci(list+1)            !1.1
        call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
      end if
    end do
  end do
end do
return

! dt(14) ar(23)-bl(32)-br(32)-
114 continue
do lri=norb_frz+1,norb_dz-1
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = mul_tab(lmi,lmj)
    if (lmij /= jmr) cycle
    iwdr = just(lri,lrj)
    do lrd=norb_frz+1,lri-1
      lmd = lsm_inn(lrd)
      lmd = mul_tab(lmd,1)
      if (lmd /= jml) cycle
      iwdl = jud(lrd)
      vlop0 = w0*w0_dt
      vlop1 = w1*w1_dt
      ni = mod(lri-lrd+norb_dz-lrj,2)
      if (ni == 0) then
        vlop0 = -vlop0
        vlop1 = -vlop1
      end if
      list = list4(lrd,lri,lrj,lra)
      wl = (vlop0-vlop1)*vint_ci(list)-2*vlop0*vint_ci(list+1) !1.1
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end do
  end do
end do
return

! d1t1(16)  ar(13)-bl(31)-br(31)-
116 continue
do lri=norb_frz+1,norb_dz-1
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = mul_tab(lmi,lmj)
    if (lmij /= jmr) cycle
    iwdr = just(lri,lrj)
    do lrd=norb_frz+1,lri-1
      lmd = lsm_inn(lrd)
      lmd = mul_tab(lmd,1)
      if (lmd /= jml) cycle
      iwdl = jud(lrd)
      vlop0 = w0*w0_d1t1
      vlop1 = w1*w1_d1t1
      ni = mod(lri-lrd+norb_dz-lrj,2)
      if (ni == 0) then
        vlop0 = -vlop0
        vlop1 = -vlop1
      end if
      list = list4(lrd,lri,lrj,lra)
      wl = (vlop0-vlop1)*vint_ci(list)-2*vlop0*vint_ci(list+1) !1.1
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end do
  end do
end do
return
! line=24:-a^r<-->sd(6),sdd(8),td(13),ttdd(15),dv(23),ddv(24)
200 continue
goto(10,10,10,10,10,206,10,208,10,10,10,10,213,10,215,10,10,10,10,10,10,10,223,224,10,10),lpok

206 continue
call sd_head_dbl_tail_act(lra,lpcoe)
return

208 continue
call sdd_head_dbl_tail_act(lra,lpcoe)
return

!call td_head_dbl_tail_act(lra,lpcoe)
! td(13-1) (22)a&(23)
! td(13-1) a&(23)c'(22)
! td(13-5) (22)d&&l(33)b^l(23)
213 continue
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (lmi /= jmlr) cycle
  w0td1 = w0_td(1)
  w0td4 = w0_td(4)
  w0td5 = w0_td(5)
  ni = mod(norb_dz-lri,2)
  if (ni == 1) w0td1 = -w0td1
  if (ni == 1) w0td4 = -w0td4
  if (ni == 1) w0td5 = -w0td5

  ! td(13-1) a&(23)c'(22)
  do lrd=lri+1,norb_dz
    lmd = lsm_inn(lrd)
    if (lmd /= jmr) cycle
    iwdl = just(lri,lrd)
    iwdr = jud(lrd)
    vlop0 = -w0*w0td1
    list = list3(lri,lra,lri)
    wl = voint(lri,lra)+vint_ci(list)           !310,act_coe,610,7
    list = list3(lri,lra,lrd)
    wl = wl+vint_ci(list+1)
    do lr=lri+1,norb_dz
      if (lr == lrd) cycle
      list = list3(lri,lra,lr)
      wl = wl+2*vint_ci(list+1)-vint_ci(list)       !310:neoc=2,coe=
    end do
    do lrk=norb_dz+1,lra
      list = list3(lri,lra,lrk)
      kcoe = lpcoe(lrk)
      call neoc(kcoe,nocc,tcoe)
      wl = wl+nocc*(vint_ci(list+1)+tcoe*vint_ci(list))
    end do
    wl = wl*vlop0
    ! td(13-5) d&rl(33)b^l(23)c'(22)
    vlop0 = -w0*w0td5
    do lrk=1,lri-1
      list = list3(lri,lra,lrk)
      wl = wl-vlop0*(2*vint_ci(list+1)-vint_ci(list))
    end do
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end do
  !---------------------------------------------------------------------
  do lrd=norb_frz+1,lri-1
    lmd = lsm_inn(lrd)
    if (lmd /= jmr) cycle
    iwdl = just(lrd,lri)
    iwdr = jud(lrd)
    ! td(13-1) (22)a&(23)
    vlop0 = w0*w0td1
    list = list3(lri,lra,lri)
    wl = vlop0*(voint(lri,lra)+vint_ci(list))             !310,act_c
    do lr=lri+1,norb_dz
      list = list3(lri,lra,lr)
      wl = wl+vlop0*(2*vint_ci(list+1)-vint_ci(list)) !  310:neoc=2,
    end do
    do lrk=norb_dz+1,lra
      list = list3(lri,lra,lrk)
      kcoe = lpcoe(lrk)
      call neoc(kcoe,nocc,tcoe)
      wl = wl+vlop0*nocc*(vint_ci(list+1)+tcoe*vint_ci(list))
    end do
    !wl = wl*vlop0
    ! td(13-4) d&r&l(22)b^l(23)
    vlop0 = w0*w0td4
    vlop1 = w1*w0td4
    list = list3(lri,lra,lrd)
    wl = wl+(vlop0-vlop1)*vint_ci(list)-2*vlop0*vint_ci(list+1)
    ! td(13-5) d&rl(33)c"(22)b^l(23)
    vlop0 = w0*w0td5
    do lrk=1,lri-1
      if (lrk == lrd) cycle
      list = list3(lri,lra,lrk)
      wl = wl+vlop0*(vint_ci(list)-2*vint_ci(list+1))      !4.3
    end do
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end do
end do
do lri=norb_frz+1,norb_dz-1
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = mul_tab(lmi,lmj)
    if (lmij /= jml) cycle
    iwdl = just(lri,lrj)

    ! td(13-2) a&(23)b&r(23)b^r(32)
    do lrd=lrj+1,norb_dz
      lmd = lsm_inn(lrd)
      if (lmd /= jmr) cycle
      w0td2 = w0_td(2)
      w1td2 = w1_td(2)
      ni = mod(lrj-lri+norb_dz-lrd,2)
      if (ni == 0) w0td2 = -w0td2
      if (ni == 0) w1td2 = -w1td2

      iwdr = jud(lrd)
      vlop0 = w0*w0td2
      vlop1 = w1*w1td2
      list = list4(lri,lrj,lrd,lra)
      wl = vlop0*(vint_ci(list+2)+vint_ci(list))-vlop1*(vint_ci(list+2)-vint_ci(list)) !1.3
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end do
    ! td(13-3) a&(23)b&l(32)b^l(23)
    do lrd=lri+1,lrj-1
      lmd = lsm_inn(lrd)
      if (lmd /= jmr) cycle
      iwdr = jud(lrd)
      w0td3 = w0_td(3)
      w1td3 = w1_td(3)
      ni = mod(lrd-lri+norb_dz-lrj,2)
      if (ni == 0) w0td3 = -w0td3
      if (ni == 0) w1td3 = -w1td3
      vlop0 = w0*w0td3                !d6-8
      vlop1 = w1*w1td3
      list = list4(lri,lrd,lrj,lra)
      wl = vlop0*(vint_ci(list+2)-2*vint_ci(list+1))-vlop1*vint_ci(list+2) !1.2
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end do
  end do
end do
return

215 continue
call ttdd_head_dbl_tail_act(lra,lpcoe)
return

! dv(23-1) ar(23)-
! dv(23-2) drl(33)-bl(23)-
223 continue
iwdr = 0
do lrd=norb_frz+1,norb_dz
  imd = lsm_inn(lrd)
  if (imd /= jml) cycle
  iwdl = jud(lrd)
  w0dv1 = w0_dv(1)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv1 = -w0dv1
  vlop0 = w0*w0dv1                !d23-1
  vlop1 = w1*w0dv1
  !*********************************************************************
  lr0 = lrd
  lr = kk(jpel)-1
  list = list3(lr0,lr,lr0)
  wl = vlop0*(voint(lr0,lr)+vint_ci(list))       !310+710
  do l=lr0+1,norb_dz
    list = list3(lr0,lr,l)
    nocc = 2
    tcoe = -0.5d0
    wl = wl+nocc*vlop0*(vint_ci(list+1)+tcoe*vint_ci(list))  !dbl_
  end do
  do l=norb_dz+1,lr
    list = list3(lr0,lr,l)
    kcoe = lpcoe(l)
    call neoc(kcoe,nocc,tcoe)
    wl = wl+nocc*vlop0*(vint_ci(list+1)+tcoe*vint_ci(list))   !act_c
  end do
  wl_430 = 0.d0
  w0dv2 = w0_dv(2)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv2 = -w0dv2
  do lrk=1,lrd-1
    list = list3(lr0,lr,lrk)
    vlop0 = w0*w0dv2
    wl_430 = wl_430+vlop0*(vint_ci(list)-2*vint_ci(list+1))
  end do
  wl = wl+wl_430
  !*********************************************************************
  call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
end do
return
! d1v(24-1) ar(13)-
! d1v(24-2) drl(33)-bl(13)-
224 continue
iwdr = 0
do lrd=norb_frz+1,norb_dz
  imd = lsm_inn(lrd)
  if (imd /= jml) cycle
  iwdl = jud(lrd)
  w0dv1 = w0_d1v(1)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv1 = -w0dv1
  vlop0 = w0*w0dv1                !d24-1
  vlop1 = w1*w0dv1
  !*********************************************************************
  lr0 = lrd
  lr = kk(jpel)-1
  list = list3(lr0,lr,lr0)
  wl = vlop0*(voint(lr0,lr)+vint_ci(list))       !310+710
  do l=lr0+1,norb_dz
    list = list3(lr0,lr,l)
    nocc = 2
    tcoe = -0.5d0
    wl = wl+nocc*vlop0*(vint_ci(list+1)+tcoe*vint_ci(list))  !dbl_
  end do
  do l=norb_dz+1,lr
    list = list3(lr0,lr,l)
    kcoe = lpcoe(l)
    call neoc(kcoe,nocc,tcoe)
    wl = wl+nocc*vlop0*(vint_ci(list+1)+tcoe*vint_ci(list))   !act_c
  end do
  wl_430 = 0.d0
  w0dv2 = w0_d1v(2)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv2 = -w0dv2
  do lrk=1,lrd-1
    list = list3(lr0,lr,lrk)
    vlop0 = w0*w0dv2
    wl_430 = wl_430+vlop0*(vint_ci(list)-2*vint_ci(list+1))
  end do
  wl = wl+wl_430
  !*********************************************************************
  call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
end do
return
! line=25:-d^r^r<-->sv(10),tv(17),ttv(18)
300 continue
goto(10,10,10,10,10,10,10,10,10,310,10,10,10,10,10,10,317,318,10,10,10,10,10,10,10,10),lpok

310 continue
call sv_head_dbl_tail_act(lra)
return

! tv(17) ar(23)-br(23)-
317 continue
iwdr = 0
do lri=norb_frz+1,norb_dz
  imi = lsm_inn(lri)
  do lrj=lri,norb_dz
    imj = lsm_inn(lrj)
    imij = mul_tab(imi,imj)
    if (imij /= jml) cycle
    iwdl = just(lri,lrj)
    vlop1 = w1*w1_tv             !d17 vlop0=0
    list = list3(lri,lrj,lra)
    wl = vlop1*vint_ci(list)        !2.1                  !!!!!
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end do
end do
return

! t1v(18) ar(13)-br(13)-
318 continue
iwdr = 0
do lri=norb_frz+1,norb_dz
  imi = lsm_inn(lri)
  do lrj=lri,norb_dz
    imj = lsm_inn(lrj)
    imij = mul_tab(imi,imj)
    if (imij /= jml) cycle
    iwdl = just(lri,lrj)
    vlop1 = w1*w1_t1v             !d18 vlop0=0
    list = list3(lri,lrj,lra)
    wl = vlop1*vint_ci(list)        !2.1                  !!!!!
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end do
end do
return   ! tmp for spin=0

! line=26:-d^r^l<-->ss(1),st(2),ts(3),stt(4),tts(5),tt(11),tttt(12),dd(19
400 continue
goto(401,402,403,404,405,10,10,10,10,10,411,412,10,10,10,10,10,10,419,420,421,422,10,10,425,10),lpok
401 continue
call ss_head_dbl_tail_act(lra)
return

402 continue
call st_head_dbl_tail_act(lra)
return
!=======================================================================
403 continue
call ts_head_dbl_tail_act(lra)
return

404 continue
call stt_head_dbl_tail_act(lra)
return

405 continue
call tts_head_dbl_tail_act(lra)
return

411 continue
call tt_head_dbl_tail_act(lra)
return

412 continue
call tttt_head_dbl_tail_act(lra)
return

419 continue
call dd_head_dbl_tail_act(lra)
return

420 continue
call dddd_head_dbl_tail_act(lra)
return

421 continue
call dd1_head_dbl_tail_act(lra)
return

422 continue
call d1d_head_dbl_tail_act(lra)
return

! vv(25) drl(33)-
425 continue
if (jwl == jwr) return
vlop0 = w0*w0_vv             !d25
wl = 0.d0
iwdl = 0
iwdr = 0
do lri=1,norb_dz
  wl = wl+vlop0*voint(lri,lra)
end do
call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
return

! sv(10),tv(17),ttv(18)
! line=27:-b^r-a^r<-->sv(10),tv(17),ttv(18)
500 continue
goto(10,10,10,10,10,10,10,10,10,510,10,10,10,10,10,10,517,518,10,10,10,10,10,10,10,10),lpok

510 continue
call sv_head_dbl_tail_act(lra)
return

! tv(17) ar(23)-br(23)-
517 continue
iwdr = 0
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = mul_tab(lmi,lmj)
    if (lmij /= jml) cycle
    w1tv = w1_tv
    if (mod(lrj-lri,2) == 0) w1tv = -w1tv
    iwdl = just(lri,lrj)
    vlop1 = w1*w1tv             !d17
    list = list4(lri,lrj,lrs,lra)
    wl = vlop1*(vint_ci(list)-vint_ci(list+2)) !1.3 vlop0=0      !!!!!
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end do
end do
return
! t1v(18) ar(13)-br(13)-
518 continue
iwdr = 0
do lri=norb_frz+1,norb_dz-1
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = mul_tab(lmi,lmj)
    if (lmij /= jml) cycle
    w1tv = w1_t1v
    if (mod(lrj-lri,2) == 0) w1tv = -w1tv
    iwdl = just(lri,lrj)
    vlop1 = w1*w1tv             !d17
    list = list4(lri,lrj,lrs,lra)
    wl = vlop1*(vint_ci(list)-vint_ci(list+2)) !1.3 vlop0=0      !!!!!
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end do
end do
return   ! tmp for spin=0

! line=28:-b^l-a^r<-->ss(1),st(2),ts(3),stt(4),tts(5),tt(11),tttt(12),dd(
600 continue
goto(601,602,603,604,605,10,10,10,10,10,611,612,10,10,10,10,10,10,619,620,621,622,10,10,625,10),lpok
601 continue
call ss_head_dbl_tail_act(lra)
return

602 continue
call st_head_dbl_tail_act(lra)
return

603 continue
call ts_head_dbl_tail_act(lra)
return

604 continue
call stt_head_dbl_tail_act(lra)
return

605 continue
call tts_head_dbl_tail_act(lra)
return

611 continue
call tt_head_dbl_tail_act(lra)
return

612 continue
call tttt_head_dbl_tail_act(lra)
return

619 continue
call dd_head_dbl_tail_act(lra)
return

620 continue
call dddd_head_dbl_tail_act(lra)
return

621 continue
call dd1_head_dbl_tail_act(lra)
return

622 continue
call d1d_head_dbl_tail_act(lra)
return

! vv(25) drl(33)-
625 continue
vlop0 = w0*w0_vv             !d25
wl = 0.d0
iwdl = 0
iwdr = 0
do lrk=1,norb_dz
  list = list3(lrs,lra,lrk)
  wl = wl+vlop0*(vint_ci(list)-2*vint_ci(list+1))   !4.3 vlop1=0
end do
call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
return

! line=29:-b^r-a^l<-->ss(1),st(2),ts(3),stt(4),tts(5),tt(11),tttt(12),dd(
700 continue
goto(701,702,703,704,705,10,10,10,10,10,711,712,10,10,10,10,10,10,719,720,721,722,10,10,725,10),lpok
701 continue
call ss_head_dbl_tail_act(lra)
return

702 continue
call st_head_dbl_tail_act(lra)
return

703 continue
call ts_head_dbl_tail_act(lra)
return

704 continue
call stt_head_dbl_tail_act(lra)
return

705 continue
call tts_head_dbl_tail_act(lra)
return

711 continue
call tt_head_dbl_tail_act(lra)
return

712 continue
call tttt_head_dbl_tail_act(lra)
return

719 continue
call dd_head_dbl_tail_act(lra)
return

720 continue
call dddd_head_dbl_tail_act(lra)
return

721 continue
call dd1_head_dbl_tail_act(lra)
return

722 continue
call d1d_head_dbl_tail_act(lra)
return

! vv(25) drl(33)-
725 continue
if (jwl >= jwr) return
vlop0 = w0*w0_vv             !d25
wl = 0.d0
iwdl = 0
iwdr = 0
do lrk=1,norb_dz
  list = list3(lrs,lra,lrk)
  wl = vlop0*vint_ci(list)      !4.2  vlop1=0         !!!!!
end do
call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
return

! line=30:-b&r-d^r^r<-->sd(6),sdd(8),td(13),ttdd(15),dv(23),ddv(24)
800 continue
goto(10,10,10,10,10,806,10,808,10,10,10,10,813,10,815,10,10,10,10,10,10,10,823,824,10,10),lpok
! sd(6-3) a&r(13)c'(22)-
806 continue
call dbl_sd_act_comp(3,lra)
return
! sd(6-3) tmp for spin=0
!return
808 continue
call dbl_sdd_act_comp(3,lra)
return
813 continue
call dbl_td_act_comp(3,lra)
return

815 continue
call dbl_ttdd_act_comp(3,lra)
return
! dv(23-1) ar(23)-
! dv(23-2) drl(33)-bl(23)-
823 continue
iwdr = 0
do lrd=norb_frz+1,norb_dz
  imd = lsm_inn(lrd)
  if (imd /= jml) cycle
  iwdl = jud(lrd)
  w0dv1 = w0_dv(1)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv1 = -w0dv1
  vlop0 = w0*w0dv1                   !d23-1
  vlop1 = w1*w0dv1
  list = list3(lrd,lrg,lra)
  wl = (vlop0+vlop1)*vint_ci(list)        !2.1       !!!!!
  call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
end do
return
! d1v(24-1)  ar(13)-
824 continue
iwdr = 0
do lrd=norb_frz+1,norb_dz
  imd = lsm_inn(lrd)
  if (imd /= jml) cycle
  iwdl = jud(lrd)
  w0dv1 = w0_d1v(1)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv1 = -w0dv1
  vlop0 = w0*w0dv1                   !d23-1
  vlop1 = w1*w0dv1
  list = list3(lrd,lrg,lra)
  wl = (vlop0+vlop1)*vint_ci(list)        !2.1       !!!!!
  call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
end do
return

! line=31:-b&l-d^r^l<-->sd(6),sdd(8),td(13),ttdd(15),dv(23),ddv(24)
900 continue
goto(10,10,10,10,10,906,10,908,10,10,10,10,913,10,915,10,10,10,10,10,10,10,923,924,10,10),lpok
906 continue
call dbl_sd_act_comp(5,lra)
return

908 continue
call dbl_sdd_act_comp(5,lra)
return
913 continue
call dbl_td_act_comp(5,lra)
return

915 continue
call dbl_ttdd_act_comp(5,lra)
return
! dv(23-1) ar(23)-
! dv(23-2) drl(33)-bl(23)-
923 continue
iwdr = 0
do lrd=norb_frz+1,norb_dz
  imd = lsm_inn(lrd)
  if (imd /= jml) cycle
  iwdl = jud(lrd)
  w0dv1 = w0_dv(1)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv1 = -w0dv1
  vlop0 = w0*w0dv1                !d23-1
  vlop1 = w1*w0dv1
  list = list3(lrd,lrg,lra)
  wl = vlop0*(vint_ci(list)-2*vint_ci(list+1))-vlop1*(vint_ci(list)) !2.2          !!!!!
  call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
end do
return
! d1v(24-1)  ar(13)-
924 continue
iwdr = 0
do lrd=norb_frz+1,norb_dz
  imd = lsm_inn(lrd)
  if (imd /= jml) cycle
  iwdl = jud(lrd)
  w0dv1 = w0_d1v(1)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv1 = -w0dv1
  vlop0 = w0*w0dv1                !d23-1
  vlop1 = w1*w0dv1
  list = list3(lrd,lrg,lra)
  wl = vlop0*(vint_ci(list)-2*vint_ci(list+1))-vlop1*(vint_ci(list)) !2.2          !!!!!
  call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
end do
return

! line=32:-b&r-b^r-a^r<-->sd(6),sdd(8),td(13),ttdd(15),dv(23),ddv(24)
1000 continue
goto(10,10,10,10,10,1006,10,1008,10,10,10,10,1013,10,1015,10,10,10,10,10,10,10,1023,1024,10,10),lpok
1006 continue
call dbl_sd_act_comp(4,lra)
return
1008 continue
call dbl_sdd_act_comp(4,lra)
return
1013 continue
call dbl_td_act_comp(4,lra)
return
1015 continue
call dbl_ttdd_act_comp(4,lra)
return
! dv(23-1) ar(23)-
! dv(23-2) drl(33)-bl(23)-
1023 continue
iwdr = 0
do lrd=norb_frz+1,norb_dz
  imd = lsm_inn(lrd)
  if (imd /= jml) cycle
  iwdl = jud(lrd)
  w0dv1 = w0_dv(1)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv1 = -w0dv1
  vlop0 = w0*w0dv1                !d23-1
  vlop1 = w1*w0dv1
  list = list4(lrd,lrg,lrs,lra)
  wl = vlop0*(vint_ci(list+2)+vint_ci(list))-vlop1*(vint_ci(list+2)-vint_ci(list)) !1.3        !!!!!
  call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
end do
return
! d1v(24-1)  ar(13)-
1024 continue
iwdr = 0
do lrd=norb_frz+1,norb_dz
  imd = lsm_inn(lrd)
  if (imd /= jml) cycle
  iwdl = jud(lrd)
  w0dv1 = w0_d1v(1)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv1 = -w0dv1
  vlop0 = w0*w0dv1                !d23-1
  vlop1 = w1*w0dv1
  list = list4(lrd,lrg,lrs,lra)
  wl = vlop0*(vint_ci(list+2)+vint_ci(list))-vlop1*(vint_ci(list+2)-vint_ci(list)) !1.3        !!!!!
  call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
end do
return
! line=33:-b&l-b^r-a^l<-->sd(6),sdd(8),td(13),ttdd(15),dv(23),ddv(24)
1100 continue
goto(10,10,10,10,10,1106,10,1108,10,10,10,10,1113,10,1115,10,10,10,10,10,10,10,1123,1124,10,10),lpok
1106 continue
call dbl_sd_act_comp(6,lra)
return
1108 continue
call dbl_sdd_act_comp(6,lra)
return
1113 continue
call dbl_td_act_comp(6,lra)
return
1115 continue
call dbl_ttdd_act_comp(6,lra)
return

! dv(23-1) ar(23)-
! dv(23-2) drl(33)-bl(23)-
1123 continue
iwdr = 0
do lrd=norb_frz+1,norb_dz
  imd = lsm_inn(lrd)
  if (imd /= jml) cycle
  iwdl = jud(lrd)
  w0dv1 = w0_dv(1)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv1 = -w0dv1
  vlop0 = w0*w0dv1                !d23-1
  vlop1 = w1*w0dv1
  list = list4(lrd,lrg,lrs,lra)
  wl = vlop0*(vint_ci(list)-2*vint_ci(list+1))-vlop1*vint_ci(list) !1.1          !!!!!
  call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
end do
return
! d1v(24-1)  ar(13)-
1124 continue
iwdr = 0
do lrd=norb_frz+1,norb_dz
  imd = lsm_inn(lrd)
  if (imd /= jml) cycle
  iwdl = jud(lrd)
  w0dv1 = w0_d1v(1)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv1 = -w0dv1
  vlop0 = w0*w0dv1                !d23-1
  vlop1 = w1*w0dv1
  list = list4(lrd,lrg,lrs,lra)
  wl = vlop0*(vint_ci(list)-2*vint_ci(list+1))-vlop1*vint_ci(list) !1.1          !!!!!
  call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
end do
return
! line=34:-b&l-b^l-a^r<-->sd(6),sdd(8),td(13),ttdd(15),dv(23),ddv(24)
1200 continue
goto(10,10,10,10,10,1206,10,1208,10,10,10,10,1213,10,1215,10,10,10,10,10,10,10,1223,1224,10,10),lpok
1206 continue
call dbl_sd_act_comp(7,lra)
return
1208 continue
call dbl_sdd_act_comp(7,lra)
return
1213 continue
call dbl_td_act_comp(7,lra)
return
1215 continue
call dbl_ttdd_act_comp(7,lra)
return
! dv(23-1) ar(23)-
! dv(23-2) drl(33)-bl(23)-
1223 continue
iwdr = 0
do lrd=norb_frz+1,norb_dz
  imd = lsm_inn(lrd)
  if (imd /= jml) cycle
  iwdl = jud(lrd)
  w0dv1 = w0_dv(1)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv1 = -w0dv1
  vlop0 = w0*w0dv1                !d23-1
  vlop1 = w1*w0dv1
  list = list4(lrd,lrg,lrs,lra)
  wl = vlop0*(vint_ci(list+2)-2.0d0*vint_ci(list+1))-vlop1*vint_ci(list+2) !1.2      !!!!!
  call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
end do
return
! d1v(24-1)  ar(13)-
1224 continue
iwdr = 0
do lrd=norb_frz+1,norb_dz
  imd = lsm_inn(lrd)
  if (imd /= jml) cycle
  iwdl = jud(lrd)
  w0dv1 = w0_d1v(1)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv1 = -w0dv1
  vlop0 = w0*w0dv1                !d23-1
  vlop1 = w1*w0dv1
  list = list4(lrd,lrg,lrs,lra)
  wl = vlop0*(vint_ci(list+2)-2.0d0*vint_ci(list+1))-vlop1*vint_ci(list+2) !1.2      !!!!!
  call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
end do
return
! line=35:-d&r^l-a^l<-->sd(6),sdd(8),td(13),ttdd(15),dv(23),ddv(24)
1300 continue
goto(10,10,10,10,10,1306,10,1308,10,10,10,10,1313,10,1315,10,10,10,10,10,10,10,1323,1324,10,10),lpok
1306 continue
call dbl_sd_act_comp(2,lra)
return
1308 continue
call dbl_sdd_act_comp(2,lra)
return
1313 continue
call dbl_td_act_comp(2,lra)
return
1315 continue
call dbl_ttdd_act_comp(2,lra)
return
! dv(23-1) ar(23)-
! dv(23-2) drl(33)-bl(23)-
1323 continue
iwdr = 0
do lrd=norb_frz+1,norb_dz
  imd = lsm_inn(lrd)
  if (imd /= jml) cycle
  iwdl = jud(lrd)
  w0dv1 = w0_dv(1)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv1 = -w0dv1
  vlop0 = w0*w0dv1                !d23-1
  list = list3(lrd,lra,lrg)
  wl = vlop0*vint_ci(list)          !3.2         !!!!!
  call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
end do
return
! d1v(24-1)  ar(13)-
1324 continue
iwdr = 0
do lrd=norb_frz+1,norb_dz
  imd = lsm_inn(lrd)
  if (imd /= jml) cycle
  iwdl = jud(lrd)
  w0dv1 = w0_d1v(1)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv1 = -w0dv1
  vlop0 = w0*w0dv1                !d23-1
  list = list3(lrd,lra,lrg)
  wl = vlop0*vint_ci(list)          !3.2         !!!!!
  call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
end do
return

10 continue
return

end subroutine dbl_head_act_tail_0

subroutine dbl_td_act_comp(lin,lra)

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "onepl.fh"

! td(13-1) (22)a&(23)
! td(13-1) a&(23)c'(22)
! td(13-2) a&(23)b&r(23)b^r(32)
! td(13-3) a&(23)b&l(32)b^l(23)
! td(13-4) d&r&l(22)b^l(23)
! td(13-5) (22)d&&l(33)b^l(23)
! td(13-5) d&rl(33)c"(22)b^l(23)
! td(13-5) d&rl(33)b^l(23)c'(22)

jmlr = mul_tab(jml,jmr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (lmi /= jmlr) cycle
  w0td1 = w0_td(1)
  ni = mod(norb_dz-lri,2)
  if (ni == 1) then
    w0td1 = -w0td1
  end if

  ! td(13-1) a&(23)c'(22)
  vlop0 = -w0*w0td1
  vlop1 = -w1*w0td1
  call comp_loop(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl)
  do lrk=lri+1,norb_dz
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lri,lrk)
    iwdr = jud(lrk)
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end do
  ! td(13-1) (22)a&(23)
  vlop0 = w0*w0td1
  vlop1 = w1*w0td1
  call comp_loop(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl)
  do lrk=norb_frz+1,lri-1
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lrk,lri)
    iwdr = jud(lrk)
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end do
  !---------------------------------------------------------------------
end do

return

end subroutine dbl_td_act_comp

subroutine dbl_ttdd_act_comp(lin,lra)

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "onepl.fh"

! t1d1(15-1)  ar(13)-
! t1d1(15-1)  ar(13)-c'(11)-

jmlr = mul_tab(jml,jmr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (lmi /= jmlr) cycle
  w0td1 = w0_t1d1(1)
  ni = mod(norb_dz-lri,2)
  if (ni == 1) then
    w0td1 = -w0td1
  end if
  ! t1d1(15-1)  ar(13)-c'(11)-
  vlop0 = -w0*w0td1
  vlop1 = -w1*w0td1
  call comp_loop(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl)
  do lrk=lri+1,norb_dz
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lri,lrk)
    iwdr = jud(lrk)
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end do
  ! t1d1(15-1)  (11)ar(13)-
  vlop0 = w0*w0td1
  vlop1 = w1*w0td1
  call comp_loop(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl)
  do lrk=norb_frz+1,lri-1
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lrk,lri)
    iwdr = jud(lrk)
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end do
  !---------------------------------------------------------------------
end do

return

end subroutine dbl_ttdd_act_comp

subroutine comp_loop(line,lr0,lrg,lrs,lr,vlop0,vlop1,wl)

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"

goto(31,32,21,13,22,11,12,51,52,41,42,43),line
!-----------------------------------------------------------------------
31 continue
goto 500
!do l=norb_dz+1,lr-1
! lpcoe(l) = lp_coe(l,mpl)
!end do
!lpcoe(lr) = kcoe
!wl = voint(lr0,lr)
!do l=lr0,lr
!  list = list3(lr0,lr,l)
!  kcoe = lpcoe(l)
!  call neoc(kcoe,nocc,tcoe)
!  wl = wl+nocc*(vint_ci(list+1)+tcoe*vint_ci(list))
!end do
!wl = wl*vlop0
32 continue
list = list3(lr0,lr,lrs)   ! lrg
wl = vlop0*vint_ci(list)
goto 500
!-----------------------------------------------------------------------
11 continue
list = list4(lr0,lrg,lrs,lr)
wl = vlop0*(vint_ci(list)-2*vint_ci(list+1))-vlop1*vint_ci(list)
goto 500
!-----------------------------------------------------------------------
12 continue
list = list4(lr0,lrg,lrs,lr)
wl = vlop0*(vint_ci(list+2)-2.0d0*vint_ci(list+1))-vlop1*vint_ci(list+2)
goto 500
!-----------------------------------------------------------------------
13 continue
list = list4(lr0,lrg,lrs,lr)
wl = vlop0*(vint_ci(list+2)+vint_ci(list))-vlop1*(vint_ci(list+2)-vint_ci(list))
goto 500
!-----------------------------------------------------------------------
21 continue
list = list3(lr0,lrg,lr)
wl = (vlop0+vlop1)*vint_ci(list)
goto 500
!-----------------------------------------------------------------------
22 continue
list = list3(lr0,lrg,lr)
wl = vlop0*(vint_ci(list)-2.0d0*vint_ci(list+1))-vlop1*vint_ci(list)
goto 500
!-----------------------------------------------------------------------
51 continue
wl = vlop0*voint(lr,lr0)*0.5d0
goto 500
!-----------------------------------------------------------------------
52 continue
wl = (vlop0-vlop1)*voint(lr,lr0)
goto 500
!-----------------------------------------------------------------------
41 continue
list = list3(lrs,lr,lr0)      ! lrg
wl = (vlop0+vlop1)*vint_ci(list)
goto 500
!----------------------------------------------------------------------
42 continue
list = list3(lrs,lr,lr0)      ! lrg
wl = (vlop0-vlop1)*vint_ci(list)
goto 500
!-----------------------------------------------------------------------
43 continue
list = list3(lrs,lr,lr0)      ! lrg
wl = vlop0*(vint_ci(list)-2.0d0*vint_ci(list+1))-vlop1*vint_ci(list)
!-----------------------------------------------------------------------
500 continue

return

end subroutine comp_loop

subroutine dbl_sd_act_comp(lin,lra)

! sd(6-1) a&r(02)-
! sd(6-2) (22)a&(13)-
! sd(6-3) a&r(13)c'(22)-
! sd(6-4) a&r(23)c'(12)-

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "onepl.fh"

jmlr = mul_tab(jml,jmr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (lmi /= jmlr) cycle
  w0sd1 = w0_sd(1)
  w0sd2 = w0_sd(2)
  w0sd3 = w0_sd(3)
  w0sd4 = w0_sd(4)
  ni = mod(norb_dz-lri,2)
  if (ni == 1) then
    w0sd1 = -w0sd1
    w0sd2 = -w0sd2
    w0sd3 = -w0sd3
    w0sd4 = -w0sd4
  end if
  if ((jml == 1) .and. (lmi == jmr)) then
    ! sd(6-1) a&r(02)-
    iwdl = just(lri,lri)
    iwdr = jud(lri)
    vlop0 = w0*w0sd1
    vlop1 = w1*w0sd1
    call comp_loop(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl)
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end if
  ! sd(6-2) (22)a&(13)-
  vlop0 = w0*w0sd2
  vlop1 = w1*w0sd2
  call comp_loop(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl)
  do lrk=norb_frz+1,lri-1
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lrk,lri)
    iwdr = jud(lrk)
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end do
  ! sd(6-4) a&r(23)-c'(12)-
  vlop0 = -w0*w0sd4
  vlop1 = -w1*w0sd4
  call comp_loop(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl)
  do lrk=lri+1,norb_dz
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lri,lrk)
    iwdr = jud(lrk)
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end do
  ! sd(6-3) a&r(13)c'(22)-
  if (jb_sys > 0) then
    vlop0 = -w0*w0sd3
    vlop1 = -w1*w0sd3
    call comp_loop(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl)
    do lrk=lri+1,norb_dz
      lmk = lsm_inn(lrk)
      if (lmk /= jmr) cycle
      iwdl = just(lrk,lri)
      iwdr = jud(lrk)
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end do
  end if
end do

return

end subroutine dbl_sd_act_comp

subroutine dbl_sdd_act_comp(lin,lra)

!sd1(8-1)    ar(01)-
!sd1(8-2)    (11)ar(23)-
!sd1(8-3)    ar(13)-c'(21)-
!sd1(8-4)    ar(23)-c'(11)-

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "onepl.fh"

jmlr = mul_tab(jml,jmr)
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
    iwdl = just(lri,lri)
    iwdr = jud(lri)
    vlop0 = w0*w0sd1
    vlop1 = w1*w0sd1
    call comp_loop(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl)
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end if
  ! sd1(8-2)    (11)ar(23)-
  vlop0 = w0*w0sd2
  vlop1 = w1*w0sd2
  call comp_loop(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl)
  do lrk=norb_frz+1,lri-1
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lri,lrk)
    iwdr = jud(lrk)
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end do
  ! sd1(8-4)    ar(23)-c'(11)-
  vlop0 = -w0*w0sd4
  vlop1 = -w1*w0sd4
  call comp_loop(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl)
  do lrk=lri+1,norb_dz
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lri,lrk)
    iwdr = jud(lrk)
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end do
  ! sd1(8-3)    ar(13)-c'(21)-
  if (jb_sys > 0) then
    vlop0 = -w0*w0sd3
    vlop1 = -w1*w0sd3
    call comp_loop(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl)
    do lrk=lri+1,norb_dz
      lmk = lsm_inn(lrk)
      if (lmk /= jmr) cycle
      iwdl = just(lrk,lri)
      iwdr = jud(lrk)
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end do
  end if
end do

return

end subroutine dbl_sdd_act_comp

subroutine ext_space_loop()

#include "drt_h.fh"
#include "gext_sequence.fh"

ism = 0
do inx=18,25
  ism = ism+1
  if (nu_ae(inx) == 0) cycle
  isegsta = iseg_sta(inx)
  isegupwei = iseg_upwei(inx)
  isegdownwei = iseg_downwei(inx)
  call g_ss_ext_sequence(ism,4)
end do
ism = 0
do inx=10,17
  ism = ism+1
  if (nu_ae(inx) == 0) cycle
  isegsta = iseg_sta(inx)
  isegupwei = iseg_upwei(inx)
  isegdownwei = iseg_downwei(inx)
  call g_tt_ext_sequence(ism)
end do
ism = 0
do inx=2,9
  ism = ism+1
  if (nu_ae(inx) == 0) cycle
  isegsta = iseg_sta(inx)
  isegupwei = iseg_upwei(inx)
  isegdownwei = iseg_downwei(inx)
  call g_dd_ext_sequence(ism)
end do

return

end subroutine ext_space_loop

subroutine g_tt_ext_sequence(ism)

#include "drt_h.fh"
#include "intsort_h.fh"
#include "gext_sequence.fh"

icano_nnsta = 2
icnt_base = 0
do ismd=1,ng_sm
  ismc = mul_tab(ism,ismd)
  if (ismc > ismd) cycle
  id_sta = ibsm_ext(ismd)
  idsta = id_sta
  idend = iesm_ext(ismd)
  ic_sta = ibsm_ext(ismc)
  icend = iesm_ext(ismc)
  if (ismd == ismc) idsta = idsta+1
  do id=idsta,idend
    m_jd = id-id_sta+1
    do ic=ic_sta,min(icend,id-1)
      m_jc = ic-ic_sta+1
      icano_nn = iwt_orb_ext(ic,id)
      if (icnt_base+icano_nn-1 > max_tmpvalue) then
        call complete_ext_loop()
        icnt_base = 0
        icano_nnsta = icano_nn
      end if
      icano_nnend = icano_nn
      do ismb=1,ismd-1
        isma = mul_tab(ism,ismb)
        if (isma > ismb) cycle
        if (ismc > ismb) then
          call g12_t_diffsym(isma,ismb,ismc)
        else if (ismc > isma) then
          call g11a_t_diffsym(isma,ismb,ismc)
        else
          call g11b_t_diffsym(isma,ismb,ismc)
        end if
      end do
      if (ism == 1) then
        isma = ismd
        call g1112_t_symaaaa(isma,ic,id)
      else
        isma = mul_tab(ism,ismd)
        call g11a11b_t_symaacc(isma,ismd,ic,id)
      end if
      call g36_t_ext(ismc,ic,id)
      !write(6,*) '  g36',value_lpext(1)
      call g5_t_ext(ismd,ic,id)
      !write(6,*) '  g5' ,value_lpext(1)
      if (ism == 1) call g9_t_ext(ismd,ic,id)
      !write(6,*) '  g9',value_lpext(1)
      icnt_base = icnt_base+icano_nn-1
    end do
  end do
end do
call complete_ext_loop()

end subroutine g_tt_ext_sequence

subroutine dbl_head_act_tail(lpcoe)

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
dimension lpcoe(norb_dz+1:norb_inn)
#include "onepl.fh"

lra = kk(jpel)-1
jml = mod((jpadl-1),8)
jmr = mod((jpad-1),8)
itypadl = (jpadl-1)/8+2
itypadr = (jpad-1)/8+2
if (jml == 0) then
  jml = 8
  itypadl = itypadl-1
end if
if (jmr == 0) then
  jmr = 8
  itypadr = itypadr-1
end if
if (jpadl == 1) itypadl = 1
if (jpad == 1) itypadr = 1
if (jpadl == 1) jml = ns_sm
if (jpad == 1) jmr = ns_sm
jml = mul_tab(jml,ns_sm)
jmr = mul_tab(jmr,ns_sm)
jmlr = mul_tab(jml,jmr)
lpok = map_jplr(itypadl,itypadr)
if (lpok == 0) return
! 23   24 25 26 27 28 29 30 31 32 33 34 35
goto(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300),line-22
! line=23:-a^l<-->ds(7),dds(9),dt(14),ddtt(16)
100 continue
goto(10,10,10,10,10,10,107,10,109,10,10,10,10,114,10,116,10,10,10,10,10,10,10,10,10,10),lpok
! ds(7-2) ar(23)-bl(31)-br(32)-

! ds(7-1) ar(23)-drl(30)-
107 continue
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (jmr /= 1) goto 106
  iwdr = just(lri,lri)
  do lrd=norb_frz+1,lri-1
    lmd = lsm_inn(lrd)
    if (lmd /= jml) cycle
    iwdl = jud(lrd)
    w0ds1 = w0_ds(1)
    ni = mod(norb_dz-lri+lri-lrd,2)
    if (ni == 0) w0ds1 = -w0ds1
    vlop0 = w0*w0ds1
    list = list3(lrd,lra,lri)
    wl = vlop0*vint_ci(list)          !   3.2
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end do
106 continue
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = mul_tab(lmi,lmj)
    if (lmij /= jmr) cycle
    do lrd=norb_frz+1,lri-1
      lmd = lsm_inn(lrd)
      if (lmd /= jml) cycle
      list = list4(lrd,lri,lrj,lra)
      iwdl = jud(lrd)
      w0ds2 = w0_ds(2)
      w1ds2 = w1_ds(2)
      w0ds3 = w0_ds(3)
      w1ds3 = w1_ds(3)
      ni = mod(norb_dz-lrj+lri-lrd,2)
      if (ni == 0) then
        w0ds2 = -w0ds2
        w1ds2 = -w1ds2
        w0ds3 = -w0ds3
        w1ds3 = -w1ds3
      end if
      ! ds(7-3) ar(23)-bl(32)-br(31)-
      iwdr = just(lri,lrj)
      vlop0 = w0*w0ds3
      vlop1 = w1*w1ds3
      wl = (vlop0-vlop1)*vint_ci(list)-2*vlop0*vint_ci(list+1)            !1.1
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
      if (jb_sys > 0) then
        ! ds(7-2) ar(23)-bl(31)-br(32)-         the symmetry problem
        iwdr = just(lrj,lri)
        vlop0 = w0*w0ds2
        vlop1 = w1*w1ds2
        wl = (vlop0-vlop1)*vint_ci(list)-2*vlop0*vint_ci(list+1)            !1.1
        call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
      end if
    end do
  end do
end do
goto 10
109 continue
do lri=norb_frz+1,norb_dz
  ! d1s(9-1) ar(13)-drl(30)-
  lmi = lsm_inn(lri)
  do lrd=norb_frz+1,lri-1
    lmd = lsm_inn(lrd)
    if ((lmd == jml) .and. (jmr == 1)) then
      iwdr = just(lri,lri)
      iwdl = jud(lrd)
      w0ds1 = w0_d1s(1)
      ni = mod(norb_dz-lri+lri-lrd,2)
      if (ni == 0) w0ds1 = -w0ds1
      vlop0 = w0*w0ds1
      list = list3(lrd,lra,lri)
      wl = vlop0*vint_ci(list)          !   3.2
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end if
    ! d1s(9-4) drl(12)-br(31)-
    if ((jml == lmd) .and. (jmr == mul_tab(lmd,lmi))) then
      iwdr = just(lrd,lri)
      iwdl = jud(lrd)
      w1ds = w1_d1s(4)
      if (mod(norb_dz-lri,2) == 1) w1ds = -w1ds
      vlop1 = w1*w1ds
      list = list3(lri,lra,lrd)
      wl = -vlop1*vint_ci(list)
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end if
  end do
end do
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    ! d1s(9-3) ar(13)-bl(32)-br(31)-
    lmj = lsm_inn(lrj)
    lmij = mul_tab(lmi,lmj)
    if (lmij /= jmr) cycle
    do lrd=norb_frz+1,lri-1
      iwdr = just(lri,lrj)
      lmd = lsm_inn(lrd)
      if (lmd /= jml) cycle
      iwdl = jud(lrd)
      w0ds3 = w0_d1s(3)
      w1ds3 = w1_d1s(3)
      ni = mod(norb_dz-lrj+lri-lrd,2)
      if (ni == 0) w0ds3 = -w0ds3
      if (ni == 0) w1ds3 = -w1ds3
      vlop0 = w0*w0ds3
      vlop1 = w1*w1ds3
      list = list4(lrd,lri,lrj,lra)
      wl = (vlop0-vlop1)*vint_ci(list)-2*vlop0*vint_ci(list+1)            !1.1
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
      if (jb_sys > 0) then
        ! d1s(9-2)   ar(13)-bl(31)-br(32)-   the symmetry problem
        iwdr = just(lrj,lri)
        w0ds3 = w0_d1s(2)
        w1ds3 = w1_d1s(2)
        ni = mod(norb_dz-lrj+lri-lrd,2)
        if (ni == 0) w0ds3 = -w0ds3
        if (ni == 0) w1ds3 = -w1ds3
        vlop0 = w0*w0ds3
        vlop1 = w1*w1ds3
        list = list4(lrd,lri,lrj,lra)
        wl = (vlop0-vlop1)*vint_ci(list)-2*vlop0*vint_ci(list+1)            !1.1
        call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
      end if
    end do
  end do
end do
return

! dt(14) ar(23)-bl(32)-br(32)-
114 continue
do lri=norb_frz+1,norb_dz-1
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = mul_tab(lmi,lmj)
    if (lmij /= jmr) cycle
    iwdr = just(lri,lrj)
    do lrd=norb_frz+1,lri-1
      lmd = lsm_inn(lrd)
      if (lmd /= jml) cycle
      iwdl = jud(lrd)
      vlop0 = w0*w0_dt
      vlop1 = w1*w1_dt
      ni = mod(lri-lrd+norb_dz-lrj,2)
      if (ni == 0) then
        vlop0 = -vlop0
        vlop1 = -vlop1
      end if
      list = list4(lrd,lri,lrj,lra)
      wl = (vlop0-vlop1)*vint_ci(list)-2*vlop0*vint_ci(list+1) !1.1
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end do
  end do
end do
return
! d1t1(16)  ar(13)-bl(31)-br(31)-
116 continue
do lri=norb_frz+1,norb_dz-1
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = mul_tab(lmi,lmj)
    if (lmij /= jmr) cycle
    iwdr = just(lri,lrj)
    do lrd=norb_frz+1,lri-1
      lmd = lsm_inn(lrd)
      lmd = mul_tab(lmd,1)
      if (lmd /= jml) cycle
      iwdl = jud(lrd)
      vlop0 = w0*w0_d1t1
      vlop1 = w1*w1_d1t1
      ni = mod(lri-lrd+norb_dz-lrj,2)
      if (ni == 0) then
        vlop0 = -vlop0
        vlop1 = -vlop1
      end if
      list = list4(lrd,lri,lrj,lra)
      wl = (vlop0-vlop1)*vint_ci(list)-2*vlop0*vint_ci(list+1) !1.1
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end do
  end do
end do
return

! line=24:-a^r<-->sd(6),sdd(8),td(13),ttdd(15),dv(23),ddv(24)
200 continue
goto(10,10,10,10,10,206,10,208,10,10,10,10,213,10,215,10,10,10,10,10,10,10,223,224,10,10),lpok
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
! sd(6-1) a&r(02)-

206 continue
call sd_head_dbl_tail_act(lra,lpcoe)
return

208 continue
call sdd_head_dbl_tail_act(lra,lpcoe)
return
! td(13-1) (22)a&(23)
! td(13-1) a&(23)c'(22)
! td(13-5) (22)d&&l(33)b^l(23)
213 continue
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (lmi /= jmlr) cycle
  w0td1 = w0_td(1)
  w0td4 = w0_td(4)
  w0td5 = w0_td(5)
  ni = mod(norb_dz-lri,2)
  if (ni == 1) w0td1 = -w0td1
  if (ni == 1) w0td4 = -w0td4
  if (ni == 1) w0td5 = -w0td5

  ! td(13-1) a&(23)c'(22)
  do lrd=lri+1,norb_dz
    lmd = lsm_inn(lrd)
    if (lmd /= jmr) cycle
    iwdl = just(lri,lrd)
    iwdr = jud(lrd)
    vlop0 = -w0*w0td1
    list = list3(lri,lra,lri)
    wl = voint(lri,lra)+vint_ci(list)         !310,act_coe,610,7
    list = list3(lri,lra,lrd)
    wl = wl+vint_ci(list+1)
    do lr=lri+1,norb_dz
      if (lr == lrd) cycle
      list = list3(lri,lra,lr)
      wl = wl+2*vint_ci(list+1)-vint_ci(list)       !310:neoc=2,coe=
    end do
    do lrk=norb_dz+1,lra
      list = list3(lri,lra,lrk)
      kcoe = lpcoe(lrk)
      call neoc(kcoe,nocc,tcoe)
      wl = wl+nocc*(vint_ci(list+1)+tcoe*vint_ci(list))
    end do
    wl = wl*vlop0
    ! td(13-5) d&rl(33)b^l(23)c'(22)
    vlop0 = -w0*w0td5
    do lrk=1,lri-1
      list = list3(lri,lra,lrk)
      wl = wl-vlop0*(2*vint_ci(list+1)-vint_ci(list))
    end do
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end do
  !---------------------------------------------------------------------
  do lrd=norb_frz+1,lri-1
    lmd = lsm_inn(lrd)
    if (lmd /= jmr) cycle
    iwdl = just(lrd,lri)
    iwdr = jud(lrd)
    ! td(13-1) (22)a&(23)
    vlop0 = w0*w0td1
    list = list3(lri,lra,lri)
    wl = vlop0*(voint(lri,lra)+vint_ci(list))             !310,act_c
    do lr=lri+1,norb_dz
      list = list3(lri,lra,lr)
      wl = wl+vlop0*(2*vint_ci(list+1)-vint_ci(list)) !  310:neoc=2,
    end do
    do lrk=norb_dz+1,lra
      list = list3(lri,lra,lrk)
      kcoe = lpcoe(lrk)
      call neoc(kcoe,nocc,tcoe)
      wl = wl+vlop0*nocc*(vint_ci(list+1)+tcoe*vint_ci(list))
    end do
    !wl = wl*vlop0
    ! td(13-4) d&r&l(22)b^l(23)
    vlop0 = w0*w0td4
    vlop1 = w1*w0td4
    list = list3(lri,lra,lrd)
    wl = wl+(vlop0-vlop1)*vint_ci(list)-2*vlop0*vint_ci(list+1)
    ! td(13-5) d&rl(33)c"(22)b^l(23)
    vlop0 = w0*w0td5
    do lrk=1,lri-1
      if (lrk == lrd) cycle
      list = list3(lri,lra,lrk)
      wl = wl+vlop0*(vint_ci(list)-2*vint_ci(list+1))      !4.3
    end do
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end do
end do
do lri=norb_frz+1,norb_dz-1
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = mul_tab(lmi,lmj)
    if (lmij /= jml) cycle
    iwdl = just(lri,lrj)

    ! td(13-2) a&(23)b&r(23)b^r(32)
    do lrd=lrj+1,norb_dz
      lmd = lsm_inn(lrd)
      if (lmd /= jmr) cycle
      w0td2 = w0_td(2)
      w1td2 = w1_td(2)
      ni = mod(lrj-lri+norb_dz-lrd,2)
      if (ni == 0) w0td2 = -w0td2
      if (ni == 0) w1td2 = -w1td2

      iwdr = jud(lrd)
      vlop0 = w0*w0td2
      vlop1 = w1*w1td2
      list = list4(lri,lrj,lrd,lra)
      wl = vlop0*(vint_ci(list+2)+vint_ci(list))-vlop1*(vint_ci(list+2)-vint_ci(list)) !1.3
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end do
    ! td(13-3) a&(23)b&l(32)b^l(23)
    do lrd=lri+1,lrj-1
      lmd = lsm_inn(lrd)
      if (lmd /= jmr) cycle
      iwdr = jud(lrd)
      w0td3 = w0_td(3)
      w1td3 = w1_td(3)
      ni = mod(lrd-lri+norb_dz-lrj,2)
      if (ni == 0) w0td3 = -w0td3
      if (ni == 0) w1td3 = -w1td3
      vlop0 = w0*w0td3                !d6-8
      vlop1 = w1*w1td3
      list = list4(lri,lrd,lrj,lra)
      wl = vlop0*(vint_ci(list+2)-2*vint_ci(list+1))-vlop1*vint_ci(list+2) !1.2
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end do
  end do
end do
goto 10
215 continue
call ttdd_head_dbl_tail_act(lra,lpcoe)
return

! dv(23-1) ar(23)-
! dv(23-2) drl(33)-bl(23)-
223 continue
iwdr = 0
do lrd=norb_frz+1,norb_dz
  imd = lsm_inn(lrd)
  if (imd /= jml) cycle
  iwdl = jud(lrd)
  w0dv1 = w0_dv(1)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv1 = -w0dv1
  vlop0 = w0*w0dv1                !d23-1
  vlop1 = w1*w0dv1
  !*********************************************************************
  lr0 = lrd
  lr = kk(jpel)-1
  list = list3(lr0,lr,lr0)
  wl = vlop0*(voint(lr0,lr)+vint_ci(list))       !310+710
  do l=lr0+1,norb_dz
    list = list3(lr0,lr,l)
    nocc = 2
    tcoe = -0.5d0
    wl = wl+nocc*vlop0*(vint_ci(list+1)+tcoe*vint_ci(list))  !dbl_
  end do
  do l=norb_dz+1,lr
    list = list3(lr0,lr,l)
    kcoe = lpcoe(l)
    call neoc(kcoe,nocc,tcoe)
    wl = wl+nocc*vlop0*(vint_ci(list+1)+tcoe*vint_ci(list))   !act_c
  end do
  wl_430 = 0.d0
  w0dv2 = w0_dv(2)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv2 = -w0dv2
  do lrk=1,lrd-1
    list = list3(lr0,lr,lrk)
    vlop0 = w0*w0dv2
    wl_430 = wl_430+vlop0*(vint_ci(list)-2*vint_ci(list+1))
  end do
  wl = wl+wl_430
  !*********************************************************************
  call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
end do
return
! d1v(24-1) ar(13)-
! d1v(24-2) drl(33)-bl(13)-
224 continue
iwdr = 0
do lrd=norb_frz+1,norb_dz
  imd = lsm_inn(lrd)
  if (imd /= jml) cycle
  iwdl = jud(lrd)
  w0dv1 = w0_d1v(1)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv1 = -w0dv1
  vlop0 = w0*w0dv1                !d24-1
  vlop1 = w1*w0dv1
  !*********************************************************************
  lr0 = lrd
  lr = kk(jpel)-1
  list = list3(lr0,lr,lr0)
  wl = vlop0*(voint(lr0,lr)+vint_ci(list))       !310+710
  do l=lr0+1,norb_dz
    list = list3(lr0,lr,l)
    nocc = 2
    tcoe = -0.5d0
    wl = wl+nocc*vlop0*(vint_ci(list+1)+tcoe*vint_ci(list))  !dbl_
  end do
  do l=norb_dz+1,lr
    list = list3(lr0,lr,l)
    kcoe = lpcoe(l)
    call neoc(kcoe,nocc,tcoe)
    wl = wl+nocc*vlop0*(vint_ci(list+1)+tcoe*vint_ci(list))   !act_c
  end do
  wl_430 = 0.d0
  w0dv2 = w0_d1v(2)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv2 = -w0dv2
  do lrk=1,lrd-1
    list = list3(lr0,lr,lrk)
    vlop0 = w0*w0dv2
    wl_430 = wl_430+vlop0*(vint_ci(list)-2*vint_ci(list+1))
  end do
  wl = wl+wl_430
  !*********************************************************************
  call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
end do
return

! line=25:-d^r^r<-->sv(10),tv(17),ttv(18)
300 continue
goto(10,10,10,10,10,10,10,10,10,310,10,10,10,10,10,10,317,318,10,10,10,10,10,10,10,10),lpok
310 continue
call sv_head_dbl_tail_act(lra)
return
! tv(17) ar(23)-br(23)-
317 continue
iwdr = 0
do lri=norb_frz+1,norb_dz
  imi = lsm_inn(lri)
  do lrj=lri,norb_dz
    imj = lsm_inn(lrj)
    imij = mul_tab(imi,imj)
    if (imij /= jml) cycle
    iwdl = just(lri,lrj)
    vlop1 = w1*w1_tv             !d17 vlop0=0
    list = list3(lri,lrj,lra)
    wl = vlop1*vint_ci(list)        !2.1                  !!!!!
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end do
end do
return

! t1v(18) ar(13)-br(13)-
318 continue
iwdr = 0
do lri=norb_frz+1,norb_dz
  imi = lsm_inn(lri)
  do lrj=lri,norb_dz
    imj = lsm_inn(lrj)
    imij = mul_tab(imi,imj)
    if (imij /= jml) cycle
    iwdl = just(lri,lrj)
    vlop1 = w1*w1_t1v             !d18 vlop0=0
    list = list3(lri,lrj,lra)
    wl = vlop1*vint_ci(list)        !2.1                  !!!!!
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end do
end do
return

! line=26:-d^r^l<-->ss(1),st(2),ts(3),stt(4),tts(5),tt(11),tttt(12),dd(19
400 continue
goto(401,402,403,404,405,10,10,10,10,10,411,412,10,10,10,10,10,10,419,420,421,422,10,10,425,10),lpok
401 continue
call ss_head_dbl_tail_act(lra)
return

402 continue
call st_head_dbl_tail_act(lra)
return
!=======================================================================
403 continue
call ts_head_dbl_tail_act(lra)
return

404 continue
call stt_head_dbl_tail_act(lra)
return

405 continue
call tts_head_dbl_tail_act(lra)
return

411 continue
call tt_head_dbl_tail_act(lra)
return

412 continue
call tttt_head_dbl_tail_act(lra)
return

419 continue
call dd_head_dbl_tail_act(lra)
return

420 continue
call dddd_head_dbl_tail_act(lra)
return

421 continue
call dd1_head_dbl_tail_act(lra)
return

422 continue
call d1d_head_dbl_tail_act(lra)
return

! vv(25) drl(33)-
425 continue
if (jwl == jwr) return
vlop0 = w0*w0_vv             !d25
wl = 0.d0
iwdl = 0
iwdr = 0
do lri=1,norb_dz
  wl = wl+vlop0*voint(lri,lra)
end do
call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
return

! line=27:-b^r-a^r<-->sv(10),tv(17),ttv(18)
500 continue
goto(10,10,10,10,10,10,10,10,10,510,10,10,10,10,10,10,517,518,10,10,10,10,10,10,10,10),lpok
510 continue
call sv_head_dbl_tail_act(lra)
return
! tv(17) ar(23)-br(23)-
517 continue
iwdr = 0
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = mul_tab(lmi,lmj)
    if (lmij /= jml) cycle
    w1tv = w1_tv
    if (mod(lrj-lri,2) == 0) w1tv = -w1tv
    iwdl = just(lri,lrj)
    vlop1 = w1*w1tv             !d17
    list = list4(lri,lrj,lrs,lra)
    wl = vlop1*(vint_ci(list)-vint_ci(list+2)) !1.3 vlop0=0      !!!!!
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end do
end do
return
! t1v(18) ar(13)-br(13)-
518 continue
iwdr = 0
do lri=norb_frz+1,norb_dz-1
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = mul_tab(lmi,lmj)
    if (lmij /= jml) cycle
    w1tv = w1_t1v
    if (mod(lrj-lri,2) == 0) w1tv = -w1tv
    iwdl = just(lri,lrj)
    vlop1 = w1*w1tv             !d17
    list = list4(lri,lrj,lrs,lra)
    wl = vlop1*(vint_ci(list)-vint_ci(list+2)) !1.3 vlop0=0      !!!!!
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end do
end do
return

! line=28:-b^l-a^r<-->ss(1),st(2),ts(3),stt(4),tts(5),tt(11),tttt(12),dd(
600 continue
goto(601,602,603,604,605,10,10,10,10,10,611,612,10,10,10,10,10,10,619,620,621,622,10,10,625,10),lpok
601 continue
call ss_head_dbl_tail_act(lra)
return

602 continue
call st_head_dbl_tail_act(lra)
return

603 continue
call ts_head_dbl_tail_act(lra)
return

604 continue
call stt_head_dbl_tail_act(lra)
return

605 continue
call tts_head_dbl_tail_act(lra)
return

611 continue
call tt_head_dbl_tail_act(lra)
return

612 continue
call tttt_head_dbl_tail_act(lra)
return

619 continue
call dd_head_dbl_tail_act(lra)
return

620 continue
call dddd_head_dbl_tail_act(lra)
return

621 continue
call dd1_head_dbl_tail_act(lra)
return

622 continue
call d1d_head_dbl_tail_act(lra)
return

! vv(25) drl(33)-
625 continue
vlop0 = w0*w0_vv             !d25
wl = 0.d0
iwdl = 0
iwdr = 0
do lrk=1,norb_dz
  list = list3(lrs,lra,lrk)
  wl = wl+vlop0*(vint_ci(list)-2*vint_ci(list+1))   !4.3 vlop1=0
end do
call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
return

! line=29:-b^r-a^l<-->ss(1),st(2),ts(3),stt(4),tts(5),tt(11),tttt(12),dd(
700 continue
goto(701,702,703,704,705,10,10,10,10,10,711,712,10,10,10,10,10,10,719,720,721,722,10,10,725,10),lpok
701 continue
call ss_head_dbl_tail_act(lra)
return

702 continue
call st_head_dbl_tail_act(lra)
return

703 continue
call ts_head_dbl_tail_act(lra)
return

704 continue
call stt_head_dbl_tail_act(lra)
return

705 continue
call tts_head_dbl_tail_act(lra)
return

711 continue
call tt_head_dbl_tail_act(lra)
return

712 continue
call tttt_head_dbl_tail_act(lra)
return

719 continue
call dd_head_dbl_tail_act(lra)
return

720 continue
call dddd_head_dbl_tail_act(lra)
return

721 continue
call dd1_head_dbl_tail_act(lra)
return

722 continue
call d1d_head_dbl_tail_act(lra)
return

! vv(25) drl(33)-
725 continue
if (jwl >= jwr) return
vlop0 = w0*w0_vv             !d25
wl = 0.d0
iwdl = 0
iwdr = 0
do lrk=1,norb_dz
  list = list3(lrs,lra,lrk)
  wl = vlop0*vint_ci(list)      !4.2  vlop1=0         !!!!!
end do
call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
return

! line=30:-b&r-d^r^r<-->sd(6),sdd(8),td(13),ttdd(15),dv(23),ddv(24)
800 continue
goto(10,10,10,10,10,806,10,808,10,10,10,10,813,10,815,10,10,10,10,10,10,10,823,824,10,10),lpok
! sd(6-3) a&r(13)c'(22)-
806 continue
call dbl_sd_act_comp(3,lra)
return
! sd(6-3) tmp for spin=0
! return
808 continue
call dbl_sdd_act_comp(3,lra)
return
813 continue
call dbl_td_act_comp(3,lra)
return

815 continue
call dbl_ttdd_act_comp(3,lra)
return
! dv(23-1) ar(23)-
! dv(23-2) drl(33)-bl(23)-
823 continue
iwdr = 0
do lrd=norb_frz+1,norb_dz
  imd = lsm_inn(lrd)
  if (imd /= jml) cycle
  iwdl = jud(lrd)
  w0dv1 = w0_dv(1)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv1 = -w0dv1
  vlop0 = w0*w0dv1                   !d23-1
  vlop1 = w1*w0dv1
  list = list3(lrd,lrg,lra)
  wl = (vlop0+vlop1)*vint_ci(list)        !2.1       !!!!!
  call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
end do
return
! d1v(24-1)  ar(13)-
824 continue
iwdr = 0
do lrd=norb_frz+1,norb_dz
  imd = lsm_inn(lrd)
  if (imd /= jml) cycle
  iwdl = jud(lrd)
  w0dv1 = w0_d1v(1)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv1 = -w0dv1
  vlop0 = w0*w0dv1                   !d23-1
  vlop1 = w1*w0dv1
  list = list3(lrd,lrg,lra)
  wl = (vlop0+vlop1)*vint_ci(list)        !2.1
  call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
end do
return

! line=31:-b&l-d^r^l<-->sd(6),sdd(8),td(13),ttdd(15),dv(23),ddv(24)
900 continue
goto(10,10,10,10,10,906,10,908,10,10,10,10,913,10,915,10,10,10,10,10,10,10,923,924,10,10),lpok
906 continue
call dbl_sd_act_comp(5,lra)
return

908 continue
call dbl_sdd_act_comp(5,lra)
return
913 continue
call dbl_td_act_comp(5,lra)
return

915 continue
call dbl_ttdd_act_comp(5,lra)
return
! dv(23-1) ar(23)-
! dv(23-1) ar(23)-
! dv(23-2) drl(33)-bl(23)-
923 continue
iwdr = 0
do lrd=norb_frz+1,norb_dz
  imd = lsm_inn(lrd)
  if (imd /= jml) cycle
  iwdl = jud(lrd)
  w0dv1 = w0_dv(1)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv1 = -w0dv1
  vlop0 = w0*w0dv1                !d23-1
  vlop1 = w1*w0dv1
  list = list3(lrd,lrg,lra)
  wl = vlop0*(vint_ci(list)-2*vint_ci(list+1))-vlop1*(vint_ci(list)) !2.2          !!!!!
  call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
end do
return
! d1v(24-1)  ar(13)-
924 continue
iwdr = 0
do lrd=norb_frz+1,norb_dz
  imd = lsm_inn(lrd)
  if (imd /= jml) cycle
  iwdl = jud(lrd)
  w0dv1 = w0_d1v(1)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv1 = -w0dv1
  vlop0 = w0*w0dv1                !d23-1
  vlop1 = w1*w0dv1
  list = list3(lrd,lrg,lra)
  wl = vlop0*(vint_ci(list)-2*vint_ci(list+1))-vlop1*(vint_ci(list)) !2.2          !!!!!
  call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
end do
return

! line=32:-b&r-b^r-a^r<-->sd(6),sdd(8),td(13),ttdd(15),dv(23),ddv(24)
1000 continue
goto(10,10,10,10,10,1006,10,1008,10,10,10,10,1013,10,1015,10,10,10,10,10,10,10,1023,1024,10,10),lpok
1006 continue
call dbl_sd_act_comp(4,lra)
return
1008 continue
call dbl_sdd_act_comp(4,lra)
return
1013 continue
call dbl_td_act_comp(4,lra)
return
1015 continue
call dbl_ttdd_act_comp(4,lra)
return
! dv(23-1) ar(23)-
! dv(23-2) drl(33)-bl(23)-
1023 continue
iwdr = 0
do lrd=norb_frz+1,norb_dz
  imd = lsm_inn(lrd)
  if (imd /= jml) cycle
  iwdl = jud(lrd)
  w0dv1 = w0_dv(1)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv1 = -w0dv1
  vlop0 = w0*w0dv1                !d23-1
  vlop1 = w1*w0dv1
  list = list4(lrd,lrg,lrs,lra)
  wl = vlop0*(vint_ci(list+2)+vint_ci(list))-vlop1*(vint_ci(list+2)-vint_ci(list)) !1.3        !!!!!
  call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
end do
return
! d1v(24-1)  ar(13)-
1024 continue
iwdr = 0
do lrd=norb_frz+1,norb_dz
  imd = lsm_inn(lrd)
  if (imd /= jml) cycle
  iwdl = jud(lrd)
  w0dv1 = w0_d1v(1)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv1 = -w0dv1
  vlop0 = w0*w0dv1                !d23-1
  vlop1 = w1*w0dv1
  list = list4(lrd,lrg,lrs,lra)
  wl = vlop0*(vint_ci(list+2)+vint_ci(list))-vlop1*(vint_ci(list+2)-vint_ci(list)) !1.3        !!!!!
  call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
end do
return

! line=33:-b&l-b^r-a^l<-->sd(6),sdd(8),td(13),ttdd(15),dv(23),ddv(24)
1100 continue
goto(10,10,10,10,10,1106,10,1108,10,10,10,10,1113,10,1115,10,10,10,10,10,10,10,1123,1124,10,10),lpok
1106 continue
call dbl_sd_act_comp(6,lra)
return
1108 continue
call dbl_sdd_act_comp(6,lra)
return
1113 continue
call dbl_td_act_comp(6,lra)
return
1115 continue
call dbl_ttdd_act_comp(6,lra)
return

! dv(23-1) ar(23)-
! dv(23-2) drl(33)-bl(23)-
1123 continue
iwdr = 0
do lrd=norb_frz+1,norb_dz
  imd = lsm_inn(lrd)
  if (imd /= jml) cycle
  iwdl = jud(lrd)
  w0dv1 = w0_dv(1)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv1 = -w0dv1
  vlop0 = w0*w0dv1                !d23-1
  vlop1 = w1*w0dv1
  list = list4(lrd,lrg,lrs,lra)
  wl = vlop0*(vint_ci(list)-2*vint_ci(list+1))-vlop1*vint_ci(list) !1.1          !!!!!
  call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
end do
return
! d1v(24-1)  ar(13)-
1124 continue
iwdr = 0
do lrd=norb_frz+1,norb_dz
  imd = lsm_inn(lrd)
  if (imd /= jml) cycle
  iwdl = jud(lrd)
  w0dv1 = w0_d1v(1)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv1 = -w0dv1
  vlop0 = w0*w0dv1                !d23-1
  vlop1 = w1*w0dv1
  list = list4(lrd,lrg,lrs,lra)
  wl = vlop0*(vint_ci(list)-2*vint_ci(list+1))-vlop1*vint_ci(list) !1.1          !!!!!
  call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
end do
return

! line=34:-b&l-b^l-a^r<-->sd(6),sdd(8),td(13),ttdd(15),dv(23),ddv(24)
1200 continue
goto(10,10,10,10,10,1206,10,1208,10,10,10,10,1213,10,1215,10,10,10,10,10,10,10,1223,1224,10,10),lpok
1206 continue
call dbl_sd_act_comp(7,lra)
return
1208 continue
call dbl_sdd_act_comp(7,lra)
return
1213 continue
call dbl_td_act_comp(7,lra)
return
1215 continue
call dbl_ttdd_act_comp(7,lra)
return
! dv(23-1) ar(23)-
! dv(23-2) drl(33)-bl(23)-
1223 continue
iwdr = 0
do lrd=norb_frz+1,norb_dz
  imd = lsm_inn(lrd)
  if (imd /= jml) cycle
  iwdl = jud(lrd)
  w0dv1 = w0_dv(1)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv1 = -w0dv1
  vlop0 = w0*w0dv1                !d23-1
  vlop1 = w1*w0dv1
  list = list4(lrd,lrg,lrs,lra)
  wl = vlop0*(vint_ci(list+2)-2.0d0*vint_ci(list+1))-vlop1*vint_ci(list+2) !1.2      !!!!!
  call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
end do
return
! d1v(24-1)  ar(13)-
1224 continue
iwdr = 0
do lrd=norb_frz+1,norb_dz
  imd = lsm_inn(lrd)
  if (imd /= jml) cycle
  iwdl = jud(lrd)
  w0dv1 = w0_d1v(1)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv1 = -w0dv1
  vlop0 = w0*w0dv1                !d23-1
  vlop1 = w1*w0dv1
  list = list4(lrd,lrg,lrs,lra)
  wl = vlop0*(vint_ci(list+2)-2.0d0*vint_ci(list+1))-vlop1*vint_ci(list+2) !1.2      !!!!!
  call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
end do
return

! line=35:-d&r^l-a^l<-->sd(6),sdd(8),td(13),ttdd(15),dv(23),ddv(24)
1300 continue
goto(10,10,10,10,10,1306,10,1308,10,10,10,10,1313,10,1315,10,10,10,10,10,10,10,1323,1324,10,10),lpok
1306 continue
call dbl_sd_act_comp(2,lra)
return
1308 continue
call dbl_sdd_act_comp(2,lra)
return
1313 continue
call dbl_td_act_comp(2,lra)
return
1315 continue
call dbl_ttdd_act_comp(2,lra)
return
! dv(23-1) ar(23)-
! dv(23-2) drl(33)-bl(23)-
1323 continue
iwdr = 0
do lrd=norb_frz+1,norb_dz
  imd = lsm_inn(lrd)
  if (imd /= jml) cycle
  iwdl = jud(lrd)
  w0dv1 = w0_dv(1)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv1 = -w0dv1
  vlop0 = w0*w0dv1                !d23-1
  list = list3(lrd,lra,lrs)
  wl = vlop0*vint_ci(list)          !3.2         !!!!!
  call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
end do
return
! d1v(24-1)  ar(13)-
1324 continue
iwdr = 0
do lrd=norb_frz+1,norb_dz
  imd = lsm_inn(lrd)
  if (imd /= jml) cycle
  iwdl = jud(lrd)
  w0dv1 = w0_d1v(1)
  ni = mod(norb_dz-lrd,2)
  if (ni == 1) w0dv1 = -w0dv1
  vlop0 = w0*w0dv1                !d23-1
  list = list3(lrd,lra,lrg)
  wl = vlop0*vint_ci(list)          !3.2         !!!!!
  call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
end do
return

10 continue
return

end subroutine dbl_head_act_tail

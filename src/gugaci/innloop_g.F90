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

!***********************************************************************
!   ar      =1 (+a^r)     drr     =2 (+d^rr)   drl     =3 (+d^rl)
!   arbr    =4 (+d^rr)    arbl    =5 (+d^rl)   ard_l^r =6 (+a^l)
!   drrb^r  =7 (+a^r)     drlb^r  =8 (+a^l)    drlb^l  =9 (+a^r)
!   arbrb^r =10 (+a^r)    arblb^r =11 (+a^l)   arblb^l =12 (+a^r)
!   drl     =13 (*)
!***********************************************************************

subroutine inner_space_loop_g()

implicit none

!wsc0 = c_time()
call dbl_space_loop_g()
!wsc1 = c_time()
call act_space_cloop_g()
call act_space_ploop_g()
!wsc2 = c_time()
!write(u6,'(2x,2(a5,f12.4))') 'dbl',wsc1-wsc0,'act',wsc2-wsc1

return

end subroutine inner_space_loop_g

subroutine act_space_cloop_g()        ! one sub_drt

use gugaci_global, only: ipae, jpad, jpae, mxnode, ndim, norb_act, nu_ad, nu_ae
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ipae_, jpad_

if (norb_act == 0) return
do ipae_=1,25
  ipae = ipae_ ! ipae is in global module, is this necessary?
  jpae = nu_ae(ipae)
  if (jpae == 0) cycle
  do jpad_=1,mxnode
    jpad = jpad_ ! jpad is in global module, is this necessary?
    if (nu_ad(jpad) == 0) cycle
    call seg_drt()
    if (ndim == 0) cycle
    call copy_to_drtl()
    call cloop_in_act_g()
  end do
end do

return

end subroutine act_space_cloop_g

subroutine act_space_ploop_g()        ! two sub_drt with same ipae

use gugaci_global, only: ipae, jpad, jpadl, jpae, mxnode, ndim, norb_act, nu_ad, nu_ae
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ipae_, jpad_, jpadl_

if (norb_act == 0) return
do ipae_=1,25
  ipae = ipae_ ! ipae is in global module, is this necessary?
  jpae = nu_ae(ipae)
  if (jpae == 0) cycle
  do jpadl_=1,mxnode                               ! jpadl
    jpadl = jpadl_ ! jpadl is in global module, is this necessary?
    if (nu_ad(jpadl) == 0) cycle
    jpad = jpadl
    call seg_drt()
    if (ndim == 0) cycle
    call copy_to_drtl()
    do jpad_=1,mxnode                               !jpadr
      jpad = jpad_ ! jpad is in global module, is this necessary?
      if (nu_ad(jpad) == 0) cycle
      call seg_drt()
      if (ndim == 0) cycle
      !if ((ipae == 18) .and. (jpadl == 2) .and. (jpad == 1)) write(u6,*)
      call ploop_in_act_g()
    end do
  end do
end do

return

end subroutine act_space_ploop_g

subroutine cloop_in_act_g()

use gugaci_global, only: logic_br, lsm_inn, norb_dz, norb_inn
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: lmi, lmij, lmj, lmk, lml, lrai, lraj, lrak, lral, lsmij, mh

do lrai=norb_dz+1,norb_inn-1
  lmi = lsm_inn(lrai)
  do lraj=lrai+1,norb_inn
    lmj = lsm_inn(lraj)
    lmij = Mul(lmi,lmj)
    !-------------------------------------------------------------------
    ! line=8 d&r&r--d^r^r
    call head_drr_at_given_orb(mh,lrai)
    logic_br(1:mh) = .true.
    call link_c2_to_given_orb(mh,lrai+1,lraj-1)
    call tail_drr_at_given_orb(mh)
    !write(u6,'(6i6)') 8,mh,lrai,lraj,0,0
    if (mh /= 0) call act_cloop_g(8,mh,lrai,lraj,0,0)
    !-------------------------------------------------------------------
    ! line=9 d&r&l--d^r^l
    call head_drl_at_given_orb(mh,lrai)
    call link_c2_to_given_orb(mh,lrai+1,lraj-1)
    call tail_drl_at_given_orb(mh)
    !write(u6,'(6i6)') 9,mh,lrai,lraj,0,0
    if (mh /= 0) call act_cloop_g(9,mh,lrai,lraj,0,0)
    lsmij = Mul(lmi,lmj)
    !-------------------------------------------------------------------
    if (lsmij == 1) then
      !-----------------------------------------------------------------
      ! line=1 a&r--a^r
      call head_ar_at_given_orb(mh,lrai)
      call link_c1_to_given_orb_coe(mh,lrai+1,lraj-1)
      call tail_ar_at_given_orb_coe(mh,lraj)
      !write(u6,'(6i6)') 1,mh,lrai,lraj,0,0
      if (mh /= 0) call act_cloop_g(1,mh,lrai,lraj,0,0)
      !call save_clp(1,mh,lra,0)
      !-----------------------------------------------------------------
      ! line=2 a&r-d^r&l-a^l
      do lrak=lrai+1,lraj-1
        call head_ar_at_given_orb(mh,lrai)
        call link_c1_to_given_orb(mh,lrai+1,lrak-1)
        call link_d10_at_given_orb(mh)
        call link_c1_to_given_orb(mh,lrak+1,lraj-1)
        call tail_al_at_given_orb(mh)
        !write(u6,'(6i6)') 2,mh,lrai,lraj,lrak,0
        if (mh /= 0) call act_cloop_g(2,mh,lrai,lraj,lrak,0)
      end do
      !-----------------------------------------------------------------
      ! line=3 a&r-b&r-d^r^r
      do lrak=lraj+1,norb_inn
        call head_ar_at_given_orb(mh,lrai)
        call link_c1_to_given_orb(mh,lrai+1,lraj-1)
        call link_b4_at_given_orb(mh)
        logic_br(1:mh) = .true.
        call link_c2_to_given_orb(mh,lraj+1,lrak-1)
        call tail_drr_at_given_orb(mh)
        !write(u6,'(6i6)') 3,mh,lrai,lraj,lrak,0
        if (mh /= 0) call act_cloop_g(3,mh,lrai,lraj,lrak,0)
        !---------------------------------------------------------------
        ! line=5 a&r-b&l-d^r^l
        call head_ar_at_given_orb(mh,lrai)
        call link_c1_to_given_orb(mh,lrai+1,lraj-1)
        call link_b3_at_given_orb(mh)
        logic_br(1:mh) = .true.
        call link_c2_to_given_orb(mh,lraj+1,lrak-1)
        call tail_drl_at_given_orb(mh)
        !write(u6,'(6i6)') 5,mh,lrai,lraj,lrak,0
        if (mh /= 0) call act_cloop_g(5,mh,lrai,lraj,lrak,0)
      end do
      !-----------------------------------------------------------------
      do lrak=norb_dz+1,lrai-1
        ! line=10 d&rr--b^r--a^r
        call head_drr_at_given_orb(mh,lrak)
        logic_br(1:mh) = .true.
        call link_c2_to_given_orb(mh,lrak+1,lrai-1)
        call link_b2_at_given_orb(mh)
        call link_c1_to_given_orb(mh,lrai+1,lraj-1)
        call tail_ar_at_given_orb(mh)
        !write(u6,'(6i6)') 10,mh,lrai,lraj,lrak,0
        if (mh /= 0) call act_cloop_g(10,mh,lrai,lraj,lrak,0)
        !---------------------------------------------------------------
        ! line=11 d&r&l-b^r-a^l
        call head_drl_at_given_orb(mh,lrak)
        call link_c2_to_given_orb(mh,lrak+1,lrai-1)
        call link_b2_at_given_orb(mh)
        call link_c1_to_given_orb(mh,lrai+1,lraj-1)
        call tail_al_at_given_orb(mh)
        !write(u6,'(6i6)') 11,mh,lrai,lraj,lrak,0
        if (mh /= 0) call act_cloop_g(11,mh,lrai,lraj,lrak,0)
        !---------------------------------------------------------------
        ! line=12 d&r&l-b^l-a^r
        call head_drl_at_given_orb(mh,lrak)
        call link_c2_to_given_orb(mh,lrak+1,lrai-1)
        call link_b1_at_given_orb(mh)
        call link_c1_to_given_orb(mh,lrai+1,lraj-1)
        call tail_ar_at_given_orb(mh)
        !write(u6,'(6i6)') 12,mh,lrai,lraj,lrak,0
        if (mh /= 0) call act_cloop_g(12,mh,lrai,lraj,lrak,0)
      end do
    end if
    if (lraj > norb_inn-2) cycle
    do lrak=lraj+1,norb_inn
      lmk = lsm_inn(lrak)
      lmk = Mul(lmij,lmk)
      do lral=lrak+1,norb_inn
        lml = lsm_inn(lral)
        lml = Mul(lmk,lml)
        if (lml /= 1) cycle
        ! line=4  a&r--b&r--b^r--a^r
        call head_ar_at_given_orb(mh,lrai)
        call link_c1_to_given_orb(mh,lrai+1,lraj-1)
        call link_b4_at_given_orb(mh)
        logic_br(1:mh) = .true.
        call link_c2_to_given_orb(mh,lraj+1,lrak-1)
        call link_b2_at_given_orb(mh)
        call link_c1_to_given_orb(mh,lrak+1,lral-1)
        call tail_ar_at_given_orb(mh)
        !  write(u6,'(6i6)') 4,mh,lrai,lral,lraj,lrak
        if (mh /= 0) call act_cloop_g(4,mh,lrai,lral,lraj,lrak)
        !---------------------------------------------------------------
        ! line=6  a&r--b&l--b^r--a^l
        call head_ar_at_given_orb(mh,lrai)
        call link_c1_to_given_orb(mh,lrai+1,lraj-1)
        call link_b3_at_given_orb(mh)
        logic_br(1:mh) = .true.
        call link_c2_to_given_orb(mh,lraj+1,lrak-1)
        call link_b2_at_given_orb(mh)
        call link_c1_to_given_orb(mh,lrak+1,lral-1)
        call tail_al_at_given_orb(mh)
        !  write(u6,'(6i6)') 6,mh,lrai,lral,lraj,lrak
        if (mh /= 0) call act_cloop_g(6,mh,lrai,lral,lraj,lrak)
        !---------------------------------------------------------------
        ! line=7 a&r--b&l--b^l--a^r
        call head_ar_at_given_orb(mh,lrai)
        call link_c1_to_given_orb(mh,lrai+1,lraj-1)
        call link_b3_at_given_orb(mh)
        logic_br(1:mh) = .true.
        call link_c2_to_given_orb(mh,lraj+1,lrak-1)
        call link_b1_at_given_orb(mh)
        call link_c1_to_given_orb(mh,lrak+1,lral-1)
        call tail_ar_at_given_orb(mh)
        !  write(u6,'(6i6)') 7,mh,lrai,lral,lraj,lrak
        if (mh /= 0) call act_cloop_g(7,mh,lrai,lral,lraj,lrak)
        !---------------------------------------------------------------
      end do
    end do
  end do
end do

return

end subroutine cloop_in_act_g

subroutine ploop_in_act_g()

use gugaci_global, only: logic_br, norb_dz, norb_inn
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: lrai, lraj, lrak, lral, mh

do lrai=norb_dz+1,norb_inn
  ! line=25 -c"-d^r^r
  logic_br(1) = .true.
  call link_c2_to_given_orb(mh,norb_dz+1,lrai-1)
  call tail_drr_at_given_orb(mh)
  if (mh /= 0) call lp_act_tail_g(25,mh,0,lrai)
  !---------------------------------------------------------------------
  ! line=26 -c"-d^r^l
  logic_br(1) = .true.
  call link_c2_to_given_orb(mh,norb_dz+1,lrai-1)
  call tail_drl_at_given_orb(mh)
  if (mh /= 0) call lp_act_tail_g(26,mh,0,lrai)
  !=====================================================================
  ! line=23 -c'-a^l
  call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
  call tail_al_at_given_orb(mh)
  if (mh /= 0) call lp_act_tail_g(23,mh,0,0)
  !---------------------------------------------------------------------
  ! line=24 -c'-a^r
  call link_c1_to_given_orb_coe(mh,norb_dz+1,lrai-1)
  call tail_ar_at_given_orb_coe(mh,lrai)
  if (mh /= 0) call lp_act_tail_g(24,mh,0,0)
  !=====================================================================
  do lrak=lrai+1,norb_inn
    !-------------------------------------------------------------------
    ! line=30 -c'-b&r-d^r^r
    call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
    call link_b4_at_given_orb(mh)
    logic_br(1:mh) = .true.
    call link_c2_to_given_orb(mh,lrai+1,lrak-1)
    call tail_drr_at_given_orb(mh)
    if (mh /= 0) call lp_act_tail_g(30,mh,lrai,0)
    !-------------------------------------------------------------------
    ! line=31 -c'-b&l-d^r^l
    call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
    call link_b3_at_given_orb(mh)
    logic_br(1:mh) = .true.
    call link_c2_to_given_orb(mh,lrai+1,lrak-1)
    call tail_drl_at_given_orb(mh)
    if (mh /= 0) call lp_act_tail_g(31,mh,lrai,0)
  end do
  !=====================================================================
  do lrak=norb_dz+1,lrai-1
    !-------------------------------------------------------------------
    ! line=35 -c'-d^r&l-a^l
    call link_c1_to_given_orb(mh,norb_dz+1,lrak-1)
    call link_d10_at_given_orb(mh)
    call link_c1_to_given_orb(mh,lrak+1,lrai-1)
    call tail_al_at_given_orb(mh)
    if (mh /= 0) call lp_act_tail_g(35,mh,lrak,lrak)
  end do
  !=====================================================================
  do lraj=lrai+1,norb_inn
    !-------------------------------------------------------------------
    ! line=27 -c"-b^r-a^r
    call link_c2_to_given_orb(mh,norb_dz+1,lrai-1)
    call link_b2_at_given_orb(mh)
    call link_c1_to_given_orb(mh,lrai+1,lraj-1)
    call tail_ar_at_given_orb(mh)
    if (mh /= 0) call lp_act_tail_g(27,mh,0,lrai)
    !-------------------------------------------------------------------
    ! line=28 -c"-b^l-a^r
    call link_c2_to_given_orb(mh,norb_dz+1,lrai-1)
    call link_b1_at_given_orb(mh)
    call link_c1_to_given_orb(mh,lrai+1,lraj-1)
    call tail_ar_at_given_orb(mh)
    if (mh /= 0) call lp_act_tail_g(28,mh,0,lrai)
    !-------------------------------------------------------------------
    ! line=29 -c"-b^r-a^l
    call link_c2_to_given_orb(mh,norb_dz+1,lrai-1)
    call link_b2_at_given_orb(mh)
    call link_c1_to_given_orb(mh,lrai+1,lraj-1)
    call tail_al_at_given_orb(mh)
    if (mh /= 0) call lp_act_tail_g(29,mh,0,lrai)
    !-------------------------------------------------------------------
    do lral=lraj+1,norb_inn
      !-----------------------------------------------------------------
      ! line=32 -c'-b&r-c"-b^r-a^r
      call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
      call link_b4_at_given_orb(mh)
      logic_br(1:mh) = .true.
      call link_c2_to_given_orb(mh,lrai+1,lraj-1)
      call link_b2_at_given_orb(mh)
      call link_c1_to_given_orb(mh,lraj+1,lral-1)
      call tail_ar_at_given_orb(mh)
      if (mh /= 0) call lp_act_tail_g(32,mh,lrai,lraj)
      !-----------------------------------------------------------------
      ! line=33 -c'-b&l-c"-b^r-a^l
      call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
      call link_b3_at_given_orb(mh)
      logic_br(1:mh) = .true.
      call link_c2_to_given_orb(mh,lrai+1,lraj-1)
      call link_b2_at_given_orb(mh)
      call link_c1_to_given_orb(mh,lraj+1,lral-1)
      call tail_al_at_given_orb(mh)
      if (mh /= 0) call lp_act_tail_g(33,mh,lrai,lraj)
      !-----------------------------------------------------------------
      ! line=34 -c'-b&l-c"-b^l-a^r
      call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
      call link_b3_at_given_orb(mh)
      logic_br(1:mh) = .true.
      call link_c2_to_given_orb(mh,lrai+1,lraj-1)
      call link_b1_at_given_orb(mh)
      call link_c1_to_given_orb(mh,lraj+1,lral-1)
      call tail_ar_at_given_orb(mh)
      if (mh /= 0) call lp_act_tail_g(34,mh,lrai,lraj)
    end do
  end do
  !---------------------------------------------------------------------
end do

return

end subroutine ploop_in_act_g

subroutine lp_act_tail_g(lin,mh,lrg0,lrs0)

use gugaci_global, only: jpel, jper, jph_, jwl, jwr, line, lpnew_coe, lpnew_ltail, lpnew_lwei, lpnew_rtail, lpnew_rwei, lrg, lrs, &
                         mhlp, norb_dz, norb_inn, vplpnew_w0, vplpnew_w1, w0, w1
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lin, mh, lrg0, lrs0
integer(kind=iwp) :: iorb, mhlp_
integer(kind=iwp), allocatable :: lpcoe(:)

call mma_allocate(lpcoe,[norb_dz+1,norb_inn],label='lpcoe')
line = lin
lrg = lrg0
lrs = lrs0
jph_ = 0
do mhlp_=1,mh
  mhlp = mhlp_ ! mhlp is in global module, is this necessary?
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
  call dbl_head_act_tail_g(lpcoe)
end do
call mma_deallocate(lpcoe)

return

end subroutine lp_act_tail_g

subroutine act_cloop_g(lin,mh,lr0,lr,lrg0,lrs0)

use gugaci_global, only: jpel, jper, jph_, jwl, jwr, line, lpnew_coe, lpnew_head, lpnew_ltail, lpnew_lwei, lpnew_rtail, &
                         lpnew_rwei, lrg, lrs, mhlp, norb_dz, norb_inn, vplpnew_w0, vplpnew_w1
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, mh, lr0, lr, lrg0, lrs0
integer(kind=iwp) :: kcoe, l, mhlp_, nocc, nxo
real(kind=wp) :: tcoe, vlop0, vlop1, wl
integer(kind=iwp), allocatable :: lpcoe(:)

line = lin
lrg = lrg0
lrs = lrs0
do mhlp_=1,mh
  mhlp = mhlp_ ! mhlp is in global module, is this necessary?
  jph_ = lpnew_head(mhlp)
  jpel = lpnew_ltail(mhlp)
  jper = lpnew_rtail(mhlp)
  jwl = lpnew_lwei(mhlp)
  jwr = lpnew_rwei(mhlp)
  vlop0 = vplpnew_w0(mhlp)
  vlop1 = vplpnew_w1(mhlp)

  select case (line)
    case default ! (1)
      !-----------------------------------------------------------------
      call mma_allocate(lpcoe,[norb_dz+1,norb_inn],label='lpcoe')
      do l=norb_dz+1,lr
        lpcoe(l) = lpnew_coe(l,mhlp)
      end do
      !wl = voint(lr0,lr)
      wl = vlop0
      !write(nf2,*) 'ar-ar'
      call prodab_1(2,jph_,jpel,jwl,jwr,0,wl,jper,lr0,lr)
      do l=lr0,lr
        !list = list3(lr0,lr,l)
        kcoe = lpcoe(l)
        !write(nf2,*) 'arw-c-arw'
        call neoc(kcoe,nocc,tcoe)
        !wl = wl+nocc*(vint_ci(list+1)+tcoe*vint_ci(list))
        wl = vlop0*nocc
        call trans_ijkl_intpos(lr,lr0,l,l,nxo)
        !write(nf2,'(4i4)') lr,lr0,l,l
        call prodab_2(2,jph_,jpel,jwl,jwr,0,wl,jper,nxo)
        wl = vlop0*nocc*tcoe
        call trans_ijkl_intpos(lr,l,lr0,l,nxo)
        call prodab_2(2,jph_,jpel,jwl,jwr,0,wl,jper,nxo)

      end do
      !wl = wl*vlop0
      call mma_deallocate(lpcoe)

    case (2)
      !-----------------------------------------------------------------
      wl = vlop0
      !write(nf2,*) 'ar-drl-al'
      !list = list3(lr0,lr,lrg)
      !wl = vlop0*vint_ci(list)
      call trans_ijkl_intpos(lr,lrg,lr0,lrg,nxo)
      call prodab_2(2,jph_,jpel,jwl,jwr,0,wl,jper,nxo)

    case (3)
      !-----------------------------------------------------------------
      wl = vlop0+vlop1
      !write(nf2,*) 'ar-br-drr'
      !list = list3(lr0,lr,lrg)
      !wl = (vlop0+vlop1)*vint_ci(list)

      call trans_ijkl_intpos(lr,lrg,lr0,lrg,nxo)
      call prodab_2(2,jph_,jpel,jwl,jwr,0,wl,jper,nxo)

    case (4)
      !-----------------------------------------------------------------
      wl = vlop0+vlop1
      !write(nf2,*) 'ar-br-br-ar'
      !list = list4(lr0,lrg,lrs,lr)
      !wl = vlop0*(vint_ci(list+2)+vint_ci(list))-vlop1*(vint_ci(list+2)-vint_ci(list))

      call trans_ijkl_intpos(lr,lrg,lrs,lr0,nxo)
      call prodab_2(2,jph_,jpel,jwl,jwr,0,wl,jper,nxo)
      wl = vlop0-vlop1
      call trans_ijkl_intpos(lr,lr0,lrg,lrs,nxo)
      call prodab_2(2,jph_,jpel,jwl,jwr,0,wl,jper,nxo)

    case (5)
      !-----------------------------------------------------------------
      wl = vlop0-vlop1
      !write(nf2,*) 'ar-bl-drl'
      !list = list3(lr0,lr,lrg)
      !wl = vlop0*(vint_ci(list)-Two*vint_ci(list+1))-vlop1*vint_ci(list)

      call trans_ijkl_intpos(lr,lrg,lr0,lrg,nxo)
      call prodab_2(2,jph_,jpel,jwl,jwr,0,wl,jper,nxo)
      wl = -vlop0*Two
      call trans_ijkl_intpos(lr,lr0,lrg,lrg,nxo)
      call prodab_2(2,jph_,jpel,jwl,jwr,0,wl,jper,nxo)

    case (6)
      !-----------------------------------------------------------------
      wl = vlop0-vlop1
      !list = list4(lr0,lrg,lrs,lr)
      !wl = vlop0*(vint_ci(list)-Two*vint_ci(list+1))-vlop1*vint_ci(list)
      !vint_ci(numb) = vfutei(la,lc,lb,ld)
      !vint_ci(numb+1) = vfutei(la,lb,lc,ld)
      !vint_ci(numb+2) = vfutei(la,ld,lc,lb)
      !write(nf2,*) 'ar-bl-br-al'
      call trans_ijkl_intpos(lr,lrg,lrs,lr0,nxo)
      call prodab_2(2,jph_,jpel,jwl,jwr,0,wl,jper,nxo)
      wl = -vlop0*Two
      call trans_ijkl_intpos(lr,lrs,lrg,lr0,nxo)
      call prodab_2(2,jph_,jpel,jwl,jwr,0,wl,jper,nxo)

    case (7)
      !-----------------------------------------------------------------
      wl = vlop0-vlop1
      !list = list4(lr0,lrg,lrs,lr)
      !wl = vlop0*(vint_ci(list+2)-Two*vint_ci(list+1))-vlop1*vint_ci(list+2)
      !write(nf2,*) 'ar-bl-bl-ar'
      call trans_ijkl_intpos(lr,lr0,lrg,lrs,nxo)
      call prodab_2(2,jph_,jpel,jwl,jwr,0,wl,jper,nxo)
      wl = -vlop0*Two
      call trans_ijkl_intpos(lr,lrs,lrg,lr0,nxo)
      call prodab_2(2,jph_,jpel,jwl,jwr,0,wl,jper,nxo)

    case (8)
      !-----------------------------------------------------------------
      ! drr-drr
      ! wl=vlop0 but not Half*vlop0 is based on that the non-diagonal
      ! just uses the non-triangle <ci|h|cj> which designates that i > j.
      wl = vlop0

      !wl = vlop0*voint(lr,lr0)*Half
      call trans_ijkl_intpos(lr,lr0,lr,lr0,nxo)
      !write(nf2,*) 'drr-drr'
      call prodab_2(2,jph_,jpel,jwl,jwr,0,wl,jper,nxo)

    case (9)
      !-----------------------------------------------------------------
      wl = (vlop0-vlop1)*Two

      ! drl-drl
      !wl = (vlop0-vlop1)*voint(lr,lr0)
      call trans_ijkl_intpos(lr,lr0,lr,lr0,nxo)
      !write(nf2,'(a7,2i4,f18.10)') 'drl-drl',lr,lr0,wl
      call prodab_2(2,jph_,jpel,jwl,jwr,0,wl,jper,nxo)

    case (10)
      !-----------------------------------------------------------------
      wl = vlop0+vlop1
      !write(nf2,*) 'drr-br-ar'
      !list = list3(lr0,lr,lrg)
      !wl = (vlop0+vlop1)*vint_ci(list)

      call trans_ijkl_intpos(lr,lrg,lr0,lrg,nxo)
      call prodab_2(2,jph_,jpel,jwl,jwr,0,wl,jper,nxo)

    case (11)
      !-----------------------------------------------------------------
      wl = vlop0-vlop1
      !write(nf2,*) 'drl-br-al'
      !list = list3(lr0,lr,lrg)
      !wl = (vlop0-vlop1)*vint_ci(list)

      call trans_ijkl_intpos(lr,lrg,lr0,lrg,nxo)
      call prodab_2(2,jph_,jpel,jwl,jwr,0,wl,jper,nxo)

    case (12)
      !-----------------------------------------------------------------
      wl = vlop0-vlop1
      !write(nf2,*) 'drl-bl-ar'
      !list = list3(lr0,lr,lrg)
      !wl = vlop0*(vint_ci(list)-Two*vint_ci(list+1))-vlop1*vint_ci(list)

      call trans_ijkl_intpos(lr,lrg,lr0,lrg,nxo)
      call prodab_2(2,jph_,jpel,jwl,jwr,0,wl,jper,nxo)
      wl = -vlop0*Two
      call trans_ijkl_intpos(lr,lr0,lrg,lrg,nxo)
      call prodab_2(2,jph_,jpel,jwl,jwr,0,wl,jper,nxo)
  end select
end do

return

end subroutine act_cloop_g

subroutine dbl_td_act_comp_g(lin,lra)

use gugaci_global, only: jml, jmr, jpel, jper, jud, just, jwl, jwr, lrg, lrs, lsm_inn, norb_dz, norb_frz, w0, w0_td, w1
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwdl, iwdr, jmlr, list0, list1, lmi, lmk, lri, lrk, ni
real(kind=wp) :: vlop0, vlop1, w0td1, wl0, wl1

! td(13-1) (22)a&(23)
! td(13-1) a&(23)c'(22)
! td(13-2) a&(23)b&r(23)b^r(32)
! td(13-3) a&(23)b&l(32)b^l(23)
! td(13-4) d&r&l(22)b^l(23)
! td(13-5) (22)d&&l(33)b^l(23)
! td(13-5) d&rl(33)c"(22)b^l(23)
! td(13-5) d&rl(33)b^l(23)c'(22)

jmlr = Mul(jml,jmr)
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
  call comp_loop_g(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl0,list0,wl1,list1)
  do lrk=lri+1,norb_dz
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lri,lrk)
    iwdr = jud(lrk)
    !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    !call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl2,jper,list2,list3)
    call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl0,jper,list0)
    if (list1 /= 0) call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl1,jper,list1)
  end do
  ! td(13-1) (22)a&(23)
  vlop0 = w0*w0td1
  vlop1 = w1*w0td1
  call comp_loop_g(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl0,list0,wl1,list1)
  do lrk=norb_frz+1,lri-1
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lrk,lri)
    iwdr = jud(lrk)
    !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    !call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl2,jper,list2,list3)
    call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl0,jper,list0)
    if (list1 /= 0) call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl1,jper,list1)
  end do
  !---------------------------------------------------------------------
end do

return

end subroutine dbl_td_act_comp_g

subroutine dbl_ttdd_act_comp_g(lin,lra)

use gugaci_global, only: jml, jmr, jpel, jper, jud, just, jwl, jwr, lrg, lrs, lsm_inn, norb_dz, norb_frz, w0, w0_t1d1, w1
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwdl, iwdr, jmlr, list0, list1, lmi, lmk, lri, lrk, ni
real(kind=wp) :: vlop0, vlop1, w0td1, wl0, wl1

! t1d1(15-1)  ar(13)-
! t1d1(15-1)  ar(13)-c'(11)-

jmlr = Mul(jml,jmr)
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
  call comp_loop_g(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl0,list0,wl1,list1)
  do lrk=lri+1,norb_dz
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lri,lrk)
    iwdr = jud(lrk)
    !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    !call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl2,jper,list2,list3)
    call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl0,jper,list0)
    if (list1 /= 0) call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl1,jper,list1)
  end do
  ! t1d1(15-1)  (11)ar(13)-
  vlop0 = w0*w0td1
  vlop1 = w1*w0td1
  call comp_loop_g(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl0,list0,wl1,list1)
  do lrk=norb_frz+1,lri-1
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lrk,lri)
    iwdr = jud(lrk)
    !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    !call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl2,jper,list2,list3)
    call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl0,jper,list0)
    if (list1 /= 0) call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl1,jper,list1)
  end do
  !---------------------------------------------------------------------
end do

return

end subroutine dbl_ttdd_act_comp_g

subroutine comp_loop_g(line,lr0,lrg,lrs,lr,vlop0,vlop1,wl0,list0,wl1,list1)

use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: line, lr0, lrg, lrs, lr
real(kind=wp), intent(in) :: vlop0, vlop1
integer(kind=iwp), intent(out) :: list0, list1
real(kind=wp), intent(out) :: wl0, wl1
integer(kind=iwp) :: nxo

wl0 = Zero
list0 = 0
wl1 = Zero
list1 = 0
select case (line)
  case default ! (1)
    !-------------------------------------------------------------------
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
  case (2)
    !-------------------------------------------------------------------
    wl0 = vlop0
    !lis = =list3(lr0,lr,lrs)   ! lrg
    !wl = vlop0*vint_ci(list)
    !write(nf2,*) 'ar-cw-ar'
    call trans_ijkl_intpos(lr,lrs,lr0,lrs,nxo)
    list0 = nxo

  case (3)
    !-------------------------------------------------------------------
    wl0 = vlop0+vlop1
    !list = list3(lr0,lrg,lr)
    !wl = (vlop0+vlop1)*vint_ci(list)

    call trans_ijkl_intpos(lrg,lr,lr0,lr,nxo)
    list0 = nxo

  case (4)
    !-----------------------------------------------------------
    wl0 = vlop0-vlop1
    !list = list4(lr0,lrg,lrs,lr)
    !wl = vlop0*(vint_ci(list+2)+vint_ci(list))-vlop1*(vint_ci(list+2)-vint_ci(list))

    call trans_ijkl_intpos(lr,lr0,lrg,lrs,nxo)
    list0 = nxo

    wl1 = vlop0+vlop1
    call trans_ijkl_intpos(lr,lrg,lrs,lr0,nxo)
    list1 = nxo

  case (5)
    !-------------------------------------------------------------------
    wl0 = vlop0-vlop1
    !list = list3(lr0,lrg,lr)
    !wl = vlop0*(vint_ci(list)-Two*vint_ci(list+1))-vlop1*vint_ci(list)

    call trans_ijkl_intpos(lrg,lr,lr0,lr,nxo)
    list0 = nxo

    wl1 = -Two*vlop0
    call trans_ijkl_intpos(lrg,lr0,lr,lr,nxo)
    list1 = nxo

  case (6)
    !-------------------------------------------------------------------
    wl0 = vlop0-vlop1
    !list = list4(lr0,lrg,lrs,lr)
    !wl = vlop0*(vint_ci(list)-Two*vint_ci(list+1))-vlop1*vint_ci(list)

    call trans_ijkl_intpos(lr,lrg,lrs,lr0,nxo)
    list0 = nxo
    wl1 = -Two*vlop0
    call trans_ijkl_intpos(lr,lrs,lrg,lr0,nxo)
    list1 = nxo

  case (7)
    !-----------------------------------------------------------
    wl0 = vlop0-vlop1
    !list = list4(lr0,lrg,lrs,lr)
    !wl = vlop0*(vint_ci(list+2)-Two*vint_ci(list+1))-vlop1*vint_ci(list+2)

    call trans_ijkl_intpos(lr,lr0,lrg,lrs,nxo)
    list0 = nxo
    wl1 = -Two*vlop0
    call trans_ijkl_intpos(lr,lrs,lrg,lr0,nxo)
    list1 = nxo

  case (8)
    !-------------------------------------------------------------------
    ! drr-drr
    ! w0lp=Two*w0lp but not One*w0lp is based on that the non-diagonal
    ! just uses the non-triangle <ci|h|cj> which designates that i > j.
    wl0 = vlop0
    !wl = vlop0*voint(lr,lr0)*Half

    call trans_ijkl_intpos(lr,lr0,lr,lr0,nxo)
    list0 = nxo

  case (9)
    !-------------------------------------------------------------------
    wl0 = (vlop0-vlop1)*Two
    !wl = (vlop0-vlop1)*voint(lr,lr0)

    call trans_ijkl_intpos(lr,lr0,lr,lr0,nxo)
    list0 = nxo

  case (10)
    !-------------------------------------------------------------------
    wl0 = vlop0+vlop1
    !list = list3(lrs,lr,lr0)      ! lrg
    !wl = (vlop0+vlop1)*vint_ci(list)

    call trans_ijkl_intpos(lr,lr0,lrs,lr0,nxo)
    list0 = nxo

  case (11)
    !-------------------------------------------------------------------
    wl0 = vlop0-vlop1
    !list = list3(lrs,lr,lr0)      ! lrg
    !wl = (vlop0-vlop1)*vint_ci(list)

    call trans_ijkl_intpos(lr,lr0,lrs,lr0,nxo)
    list0 = nxo

  case (12)
    !-------------------------------------------------------------------
    wl0 = vlop0-vlop1
    !list = list3(lrs,lr,lr0)      ! lrg
    !wl = vlop0*(vint_ci(list)-Two*vint_ci(list+1))-vlop1*vint_ci(list)

    call trans_ijkl_intpos(lr,lr0,lrs,lr0,nxo)
    list0 = nxo

    wl1 = -Two*vlop0
    call trans_ijkl_intpos(lr,lrs,lr0,lr0,nxo)
    list1 = nxo
end select

return

end subroutine comp_loop_g

subroutine dbl_sd_act_comp_g(lin,lra)
!sd(6-1) a&r(02)-
!sd(6-2) (22)a&(13)-
!sd(6-3) a&r(13)c'(22)-
!sd(6-4) a&r(23)c'(12)-

use gugaci_global, only: jb_sys, jml, jmr, jpel, jper, jud, just, jwl, jwr, lrg, lrs, lsm_inn, norb_dz, norb_frz, w0, w0_sd, w1
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwdl, iwdr, jmlr, list0, list1, lmi, lmk, lri, lrk, ni
real(kind=wp) :: vlop0, vlop1, w0sd1, w0sd2, w0sd3, w0sd4, wl0, wl1

jmlr = Mul(jml,jmr)
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
    call comp_loop_g(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl0,list0,wl1,list1)
    !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    !call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl2,jper,list2,list3)
    call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl0,jper,list0)
    if (list1 /= 0) call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl1,jper,list1)
  end if
  ! sd(6-2) (22)a&(13)-
  vlop0 = w0*w0sd2
  vlop1 = w1*w0sd2
  call comp_loop_g(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl0,list0,wl1,list1)
  do lrk=norb_frz+1,lri-1
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lrk,lri)
    iwdr = jud(lrk)
    !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    !call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl2,jper,list2,list3)
    call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl0,jper,list0)
    if (list1 /= 0) call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl1,jper,list1)
  end do
  ! sd(6-4) a&r(23)-c'(12)-
  vlop0 = -w0*w0sd4
  vlop1 = -w1*w0sd4
  call comp_loop_g(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl0,list0,wl1,list1)
  do lrk=lri+1,norb_dz
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lri,lrk)
    iwdr = jud(lrk)
    !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    !call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl2,jper,list2,list3)
    call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl0,jper,list0)
    if (list1 /= 0) call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl1,jper,list1)
  end do
  ! sd(6-3) a&r(13)c'(22)-
  if (jb_sys > 0) then
    vlop0 = -w0*w0sd3
    vlop1 = -w1*w0sd3
    call comp_loop_g(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl0,list0,wl1,list1)
    do lrk=lri+1,norb_dz
      lmk = lsm_inn(lrk)
      if (lmk /= jmr) cycle
      iwdl = just(lrk,lri)
      iwdr = jud(lrk)
      !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
      !call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl2,jper,list2,list3)
      call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl0,jper,list0)
      if (list1 /= 0) call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl1,jper,list1)
    end do
  end if
end do

return

end subroutine dbl_sd_act_comp_g

subroutine dbl_sdd_act_comp_g(lin,lra)
! sd1(8-1)    ar(01)-
! sd1(8-2)    (11)ar(23)-
! sd1(8-3)    ar(13)-c'(21)-
! sd1(8-4)    ar(23)-c'(11)-

use gugaci_global, only: jb_sys, jml, jmr, jpel, jper, jud, just, jwl, jwr, lrg, lrs, lsm_inn, norb_dz, norb_frz, w0, w0_sd1, w1
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin, lra
integer(kind=iwp) :: iwdl, iwdr, jmlr, list0, list1, lmi, lmk, lri, lrk, ni
real(kind=wp) :: vlop0, vlop1, w0sd1, w0sd2, w0sd3, w0sd4, wl0, wl1

jmlr = Mul(jml,jmr)
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
    call comp_loop_g(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl0,list0,wl1,list1)
    !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    !call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl2,jper,list2,list3)
    call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl0,jper,list0)
    if (list1 /= 0) call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl1,jper,list1)
  end if
  ! sd1(8-2)    (11)ar(23)-
  vlop0 = w0*w0sd2
  vlop1 = w1*w0sd2
  call comp_loop_g(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl0,list0,wl1,list1)
  do lrk=norb_frz+1,lri-1
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lri,lrk)
    iwdr = jud(lrk)
    !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    !call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl2,jper,list2,list3)
    call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl0,jper,list0)
    if (list1 /= 0) call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl1,jper,list1)
  end do
  ! sd1(8-4)    ar(23)-c'(11)-
  vlop0 = -w0*w0sd4
  vlop1 = -w1*w0sd4
  call comp_loop_g(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl0,list0,wl1,list1)
  do lrk=lri+1,norb_dz
    lmk = lsm_inn(lrk)
    if (lmk /= jmr) cycle
    iwdl = just(lri,lrk)
    iwdr = jud(lrk)
    !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    !call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl2,jper,list2,list3)
    call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl0,jper,list0)
    if (list1 /= 0) call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl1,jper,list1)
  end do
  ! sd1(8-3)    ar(13)-c'(21)-
  if (jb_sys > 0) then
    vlop0 = -w0*w0sd3
    vlop1 = -w1*w0sd3
    call comp_loop_g(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl0,list0,wl1,list1)
    do lrk=lri+1,norb_dz
      lmk = lsm_inn(lrk)
      if (lmk /= jmr) cycle
      iwdl = just(lrk,lri)
      iwdr = jud(lrk)
      !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
      !call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl2,jper,list2,list3)
      call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl0,jper,list0)
      if (list1 /= 0) call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl1,jper,list1)
    end do
  end if
end do

return

end subroutine dbl_sdd_act_comp_g

subroutine ext_space_loop_g()

use gugaci_global, only: iseg_downwei, iseg_sta, iseg_upwei, isegdownwei, isegsta, isegupwei, nu_ae
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: inx, ism

ism = 0
do inx=18,25
  ism = ism+1
  if (nu_ae(inx) == 0) cycle
  isegsta = iseg_sta(inx)
  isegupwei = iseg_upwei(inx)
  isegdownwei = iseg_downwei(inx)
  call g_ss_ext_sequence_g(ism,4)
end do
ism = 0
do inx=10,17
  ism = ism+1
  if (nu_ae(inx) == 0) cycle
  isegsta = iseg_sta(inx)
  isegupwei = iseg_upwei(inx)
  isegdownwei = iseg_downwei(inx)
  call g_tt_ext_sequence_g(ism)
end do
ism = 0
do inx=2,9
  ism = ism+1
  if (nu_ae(inx) == 0) cycle
  isegsta = iseg_sta(inx)
  isegupwei = iseg_upwei(inx)
  isegdownwei = iseg_downwei(inx)
  call g_dd_ext_sequence_g(ism)
end do

return

end subroutine ext_space_loop_g

subroutine dbl_head_act_tail_g(lpcoe)

use gugaci_global, only: jb_sys, jml, jmr, jpad, jpadl, jpel, jper, jud, just, jwl, jwr, kk, line, lrg, lrs, lsm_inn, map_jplr, &
                         norb_dz, norb_frz, norb_inn, ns_sm, w0, w0_d1s, w0_d1t1, w0_d1v, w0_ds, w0_dt, w0_dv, w0_td, w0_vv, w1, &
                         w1_d1s, w1_d1t1, w1_ds, w1_dt, w1_t1v, w1_td, w1_tv
use Symmetry_Info, only: Mul
use Constants, only: Zero, Two, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lpcoe(norb_dz+1:norb_inn)
integer(kind=iwp) :: imd, imi, imij, imj, itypadl, itypadr, iwdl, iwdr, jmlr, kcoe, l, lmd, lmi, lmij, lmj, lpok, lr, lr0, lra, &
                     lrd, lri, lrj, lrk, ni, nocc, nxo
real(kind=wp) :: tcoe, vlop0, vlop1, w0ds1, w0ds2, w0ds3, w0dv1, w0dv2, w0td1, w0td2, w0td3, w0td4, w0td5, w1ds, w1ds2, w1ds3, &
                 w1td2, w1td3, w1tv, wl

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
jml = Mul(jml,ns_sm)
jmr = Mul(jmr,ns_sm)
jmlr = Mul(jml,jmr)
lpok = map_jplr(itypadl,itypadr)
if (lpok == 0) return

select case (line)
  case default ! (23)
    ! line=23:-a^l<-->ds(7),dds(9),dt(14),ddtt(16)
    select case (lpok)
      case default ! (7)
        ! ds(7-2) ar(23)-bl(31)-br(32)-

        ! ds(7-1) ar(23)-drl(30)-
        do lri=norb_frz+1,norb_dz
          lmi = lsm_inn(lri)
          if (jmr == 1) then
            iwdr = just(lri,lri)
            do lrd=norb_frz+1,lri-1
              lmd = lsm_inn(lrd)
              if (lmd /= jml) cycle
              iwdl = jud(lrd)
              w0ds1 = w0_ds(1)
              ni = mod(norb_dz-lri+lri-lrd,2)
              if (ni == 0) w0ds1 = -w0ds1
              vlop0 = w0*w0ds1
              !list = list3(lrd,lra,lri)
              call trans_ijkl_intpos(lra,lri,lrd,lri,nxo)
              !wl = vlop0*vint_ci(list)
              wl = vlop0
              !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            end do
          end if
          do lrj=lri+1,norb_dz
            lmj = lsm_inn(lrj)
            lmij = Mul(lmi,lmj)
            if (lmij /= jmr) cycle
            do lrd=norb_frz+1,lri-1
              lmd = lsm_inn(lrd)
              if (lmd /= jml) cycle
              !list = list4(lrd,lri,lrj,lra)
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
              !wl = (vlop0-vlop1)*vint_ci(list)-Two*vlop0*vint_ci(list+1)            !1.1
              wl = vlop0-vlop1
              call trans_ijkl_intpos(lra,lri,lrj,lrd,nxo)

              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
              wl = -Two*vlop0
              call trans_ijkl_intpos(lra,lrj,lri,lrd,nxo)

              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
              if (jb_sys > 0) then
                ! ds(7-2) ar(23)-bl(31)-br(32)-         the symmetry problem
                iwdr = just(lrj,lri)
                vlop0 = w0*w0ds2
                vlop1 = w1*w1ds2
                !wl = (vlop0-vlop1)*vint_ci(list)-Two*vlop0*vint_ci(list+1)            !1.1
                wl = vlop0-vlop1
                call trans_ijkl_intpos(lra,lri,lrj,lrd,nxo)

                call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
                wl = -Two*vlop0
                call trans_ijkl_intpos(lra,lrj,lri,lrd,nxo)

                call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

              end if
            end do
          end do
        end do

      case (9)
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
              !list = list3(lrd,lra,lri)
              call trans_ijkl_intpos(lra,lri,lrd,lri,nxo)

              !wl = vlop0*vint_ci(list)          !   3.2
              wl = vlop0
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            end if
            ! d1s(9-4) drl(12)-br(31)-
            if ((jml == lmd) .and. (jmr == Mul(lmd,lmi))) then
              iwdr = just(lrd,lri)
              iwdl = jud(lrd)
              w1ds = w1_d1s(4)
              if (mod(norb_dz-lri,2) == 1) w1ds = -w1ds
              vlop1 = w1*w1ds
              !list = list3(lri,lra,lrd)
              call trans_ijkl_intpos(lra,lrd,lri,lrd,nxo)

              !wl = -vlop1*vint_ci(list)
              wl = -vlop1
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            end if
          end do
        end do
        do lri=norb_frz+1,norb_dz
          lmi = lsm_inn(lri)
          do lrj=lri+1,norb_dz
            ! d1s(9-3) ar(13)-bl(32)-br(31)-
            lmj = lsm_inn(lrj)
            lmij = Mul(lmi,lmj)
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
              !list = list4(lrd,lri,lrj,lra)
              !wl = (vlop0-vlop1)*vint_ci(list)-Two*vlop0*vint_ci(list+1)            !1.1
              wl = vlop0-vlop1
              call trans_ijkl_intpos(lra,lri,lrj,lrd,nxo)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
              wl = -Two*vlop0
              call trans_ijkl_intpos(lra,lrj,lri,lrd,nxo)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

              if (jb_sys > 0) then
                !d1s(9-2)   ar(13)-bl(31)-br(32)-   the symmetry problem
                iwdr = just(lrj,lri)
                w0ds3 = w0_d1s(2)
                w1ds3 = w1_d1s(2)
                ni = mod(norb_dz-lrj+lri-lrd,2)
                if (ni == 0) w0ds3 = -w0ds3
                if (ni == 0) w1ds3 = -w1ds3
                vlop0 = w0*w0ds3
                vlop1 = w1*w1ds3
                !list = list4(lrd,lri,lrj,lra)
                !wl = (vlop0-vlop1)*vint_ci(list)-Two*vlop0*vint_ci(list+1)            !1.1
                wl = vlop0-vlop1
                call trans_ijkl_intpos(lra,lri,lrj,lrd,nxo)
                call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
                wl = -Two*vlop0
                call trans_ijkl_intpos(lra,lrj,lri,lrd,nxo)
                call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

              end if
            end do
          end do
        end do

      case (14)
        ! dt(14) ar(23)-bl(32)-br(32)-
        do lri=norb_frz+1,norb_dz-1
          lmi = lsm_inn(lri)
          do lrj=lri+1,norb_dz
            lmj = lsm_inn(lrj)
            lmij = Mul(lmi,lmj)
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
              !list = list4(lrd,lri,lrj,lra)
              !wl = (vlop0-vlop1)*vint_ci(list)-Two*vlop0*vint_ci(list+1) !1.1
              wl = vlop0-vlop1
              call trans_ijkl_intpos(lra,lri,lrj,lrd,nxo)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
              wl = -Two*vlop0
              call trans_ijkl_intpos(lra,lrj,lri,lrd,nxo)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

            end do
          end do
        end do

      case (16)
        ! d1t1(16)  ar(13)-bl(31)-br(31)-
        do lri=norb_frz+1,norb_dz-1
          lmi = lsm_inn(lri)
          do lrj=lri+1,norb_dz
            lmj = lsm_inn(lrj)
            lmij = Mul(lmi,lmj)
            if (lmij /= jmr) cycle
            iwdr = just(lri,lrj)
            do lrd=norb_frz+1,lri-1
              lmd = lsm_inn(lrd)
              lmd = Mul(lmd,1)
              if (lmd /= jml) cycle
              iwdl = jud(lrd)
              vlop0 = w0*w0_d1t1
              vlop1 = w1*w1_d1t1
              ni = mod(lri-lrd+norb_dz-lrj,2)
              if (ni == 0) then
                vlop0 = -vlop0
                vlop1 = -vlop1
              end if
              !list = list4(lrd,lri,lrj,lra)
              !wl = (vlop0-vlop1)*vint_ci(list)-Two*vlop0*vint_ci(list+1) !1.1
              wl = vlop0-vlop1
              call trans_ijkl_intpos(lra,lri,lrj,lrd,nxo)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
              wl = -Two*vlop0
              call trans_ijkl_intpos(lra,lrj,lri,lrd,nxo)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

            end do
          end do
        end do

      case (1:6,8,10:13,15,17:26)
    end select

  case (24)
    ! line=24:-a^r<-->sd(6),sdd(8),td(13),ttdd(15),dv(23),ddv(24)
    select case (lpok)
      case default ! (6)
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

        call sd_head_dbl_tail_act_g(lra,lpcoe)

      case (8)
        call sdd_head_dbl_tail_act_g(lra,lpcoe)

      case (13)
        ! td(13-1) (22)a&(23)
        ! td(13-1) a&(23)c'(22)
        ! td(13-5) (22)d&&l(33)b^l(23)
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
            !list = list3(lri,lra,lri)
            !wl = voint(lri,lra)+vint_ci(list)             !310,act_coe,610,
            wl = vlop0
            call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,lri,lra)
            call trans_ijkl_intpos(lra,lri,lri,lri,nxo)
            call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

            !list = list3(lri,lra,lrd)
            !wl = wl+vint_ci(list+1)                          !310 c'(22) co
            call trans_ijkl_intpos(lra,lri,lrd,lrd,nxo)
            call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

            do lr=lri+1,norb_dz
              if (lr == lrd) cycle
              !list = list3(lri,lra,lr)
              !wl = wl+Two*vint_ci(list+1)-vint_ci(list)       !310:neoc=2,coe
              wl = Two*vlop0
              call trans_ijkl_intpos(lra,lri,lr,lr,nxo)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
              wl = -vlop0
              call trans_ijkl_intpos(lra,lr,lri,lr,nxo)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

            end do
            do lrk=norb_dz+1,lra
              !list = list3(lri,lra,lrk)
              kcoe = lpcoe(lrk)
              call neoc(kcoe,nocc,tcoe)
              !wl = wl+nocc*(vint_ci(list+1)+tcoe*vint_ci(list))
              wl = vlop0*nocc
              call trans_ijkl_intpos(lra,lri,lrk,lrk,nxo)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
              wl = vlop0*tcoe*nocc
              call trans_ijkl_intpos(lra,lrk,lri,lrk,nxo)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

            end do
            !wl = wl*vlop0
            ! td(13-5) d&rl(33)b^l(23)c'(22)
            vlop0 = -w0*w0td5
            do lrk=1,lri-1
              !list = list3(lri,lra,lrk)
              !wl = wl-vlop0*(Two*vint_ci(list+1)-vint_ci(list))
              wl = -vlop0*2
              call trans_ijkl_intpos(lra,lri,lrk,lrk,nxo)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
              wl = vlop0
              call trans_ijkl_intpos(lra,lrk,lri,lrk,nxo)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

            end do
            !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
          end do
          !---------------------------------------------------------------------
          do lrd=norb_frz+1,lri-1
            lmd = lsm_inn(lrd)
            if (lmd /= jmr) cycle
            iwdl = just(lrd,lri)
            iwdr = jud(lrd)
            ! td(13-1) (22)a&(23)
            vlop0 = w0*w0td1
            !list = list3(lri,lra,lri)
            !wl = vlop0*(voint(lri,lra)+vint_ci(list))     !310,act_coe,610,
            wl = vlop0
            call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,lri,lra)
            call trans_ijkl_intpos(lra,lri,lri,lri,nxo)
            call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

            do lr=lri+1,norb_dz
              !list = list3(lri,lra,lr)
              !wl = wl+vlop0*(Two*vint_ci(list+1)-vint_ci(list)) !  310:neoc=2
              wl = vlop0*2
              call trans_ijkl_intpos(lra,lri,lr,lr,nxo)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
              wl = -vlop0
              call trans_ijkl_intpos(lra,lr,lri,lr,nxo)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

            end do
            do lrk=norb_dz+1,lra
              !list = list3(lri,lra,lrk)
              kcoe = lpcoe(lrk)
              call neoc(kcoe,nocc,tcoe)
              !wl = wl+vlop0*nocc*(vint_ci(list+1)+tcoe*vint_ci(list))
              wl = vlop0*nocc
              call trans_ijkl_intpos(lra,lri,lrk,lrk,nxo)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
              wl = vlop0*tcoe*nocc
              call trans_ijkl_intpos(lra,lrk,lri,lrk,nxo)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

            end do
            !wl = wl*vlop0
            ! td(13-4) d&r&l(22)b^l(23)
            vlop0 = w0*w0td4
            vlop1 = w1*w0td4
            !list = list3(lri,lra,lrd)
            !wl = wl+(vlop0-vlop1)*vint_ci(list)-Two*vlop0*vint_ci(list+1)
            wl = -Two*vlop0
            call trans_ijkl_intpos(lra,lri,lrd,lrd,nxo)
            call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            wl = vlop0-vlop1
            call trans_ijkl_intpos(lra,lrd,lri,lrd,nxo)
            call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

            ! td(13-5) d&rl(33)c"(22)b^l(23)
            vlop0 = w0*w0td5
            do lrk=1,lri-1
              if (lrk == lrd) cycle
              !list = list3(lri,lra,lrk)
              !wl = wl+vlop0*(vint_ci(list)-Two*vint_ci(list+1))      !4.3
              wl = -Two*vlop0
              call trans_ijkl_intpos(lra,lri,lrk,lrk,nxo)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
              wl = vlop0
              call trans_ijkl_intpos(lra,lrk,lri,lrk,nxo)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

            end do
            !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
          end do
        end do

        do lri=norb_frz+1,norb_dz-1
          lmi = lsm_inn(lri)
          do lrj=lri+1,norb_dz
            lmj = lsm_inn(lrj)
            lmij = Mul(lmi,lmj)
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
              !list = list4(lri,lrj,lrd,lra)
              !wl = vlop0*(vint_ci(list+2)+vint_ci(list))-vlop1*(vint_ci(list+2)-vint_ci(list))  !1.3
              wl = vlop0-vlop1
              call trans_ijkl_intpos(lra,lri,lrj,lrd,nxo)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
              wl = vlop0+vlop1
              call trans_ijkl_intpos(lra,lrj,lrd,lri,nxo)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

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
              !list = list4(lri,lrd,lrj,lra)
              !wl = vlop0*(vint_ci(list+2)-Two*vint_ci(list+1))-vlop1*vint_ci(list+2)      !1.2

              wl = vlop0-vlop1
              call trans_ijkl_intpos(lra,lri,lrd,lrj,nxo)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
              wl = -vlop0*2
              call trans_ijkl_intpos(lra,lrj,lrd,lri,nxo)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

            end do
          end do
        end do

      case (15)
        call ttdd_head_dbl_tail_act_g(lra,lpcoe)

      case (23)
        ! dv(23-1) ar(23)-
        ! dv(23-2) drl(33)-bl(23)-
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
          !list = list3(lr0,lr,lr0)
          !wl = vlop0*(voint(lr0,lr)+vint_ci(list))       !310+710
          wl = vlop0

          call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,lr0,lr)
          call trans_ijkl_intpos(lr,lr0,lr0,lr0,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
          do l=lr0+1,norb_dz
            !list=list3(lr0,lr,l)
            nocc = 2
            tcoe = -Half
            !wl = wl+nocc*vlop0*(vint_ci(list+1)+tcoe*vint_ci(list))  !dbl_
            wl = vlop0*nocc
            call trans_ijkl_intpos(lr,lr0,l,l,nxo)
            call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            wl = vlop0*tcoe*nocc
            call trans_ijkl_intpos(lr,l,lr0,l,nxo)
            call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

          end do
          do l=norb_dz+1,lr
            !list = list3(lr0,lr,l)
            kcoe = lpcoe(l)
            call neoc(kcoe,nocc,tcoe)
            !wl = wl+nocc*vlop0*(vint_ci(list+1)+tcoe*vint_ci(list))   !act_c
            wl = vlop0*nocc
            call trans_ijkl_intpos(lr,lr0,l,l,nxo)
            call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            wl = vlop0*tcoe*nocc
            call trans_ijkl_intpos(lr,l,lr0,l,nxo)
            call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

          end do
          !wl_430 = Zero
          w0dv2 = w0_dv(2)
          ni = mod(norb_dz-lrd,2)
          if (ni == 1) w0dv2 = -w0dv2
          do lrk=1,lrd-1
            !list = list3(lr0,lr,lrk)
            vlop0 = w0*w0dv2
            !wl_430 = wl_430+vlop0*(vint_ci(list)-Two*vint_ci(list+1))
            wl = vlop0
            call trans_ijkl_intpos(lr,lrk,lr0,lrk,nxo)
            call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            wl = -vlop0*2
            call trans_ijkl_intpos(lr,lr0,lrk,lrk,nxo)
            call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

          end do
          !wl = wl+wl_430
          !*********************************************************************
          !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
        end do

      case (24)
        ! d1v(24-1) ar(13)-
        ! d1v(24-2) drl(33)-bl(13)-
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
          !list = list3(lr0,lr,lr0)
          !wl = vlop0*(voint(lr0,lr)+vint_ci(list))       !310+710
          wl = vlop0
          call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,lr0,lr)
          call trans_ijkl_intpos(lr,lr0,lr0,lr0,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

          do l=lr0+1,norb_dz
            !list=list3(lr0,lr,l)
            nocc = 2
            tcoe = -Half
            !wl = wl+nocc*vlop0*(vint_ci(list+1)+tcoe*vint_ci(list))  !dbl_
            wl = vlop0*nocc
            call trans_ijkl_intpos(lr,lr0,l,l,nxo)
            call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            wl = vlop0*tcoe*nocc
            call trans_ijkl_intpos(lr,l,lr0,l,nxo)
            call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

          end do
          do l=norb_dz+1,lr
            !list = list3(lr0,lr,l)
            kcoe = lpcoe(l)
            call neoc(kcoe,nocc,tcoe)
            !wl = wl+nocc*vlop0*(vint_ci(list+1)+tcoe*vint_ci(list))   !act_c
            wl = vlop0*nocc
            call trans_ijkl_intpos(lr,lr0,l,l,nxo)
            call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            wl = vlop0*tcoe*nocc
            call trans_ijkl_intpos(lr,l,lr0,l,nxo)
            call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

          end do
          !wl_430 = Zero
          w0dv2 = w0_d1v(2)
          ni = mod(norb_dz-lrd,2)
          if (ni == 1) w0dv2 = -w0dv2
          do lrk=1,lrd-1
            !list = list3(lr0,lr,lrk)
            vlop0 = w0*w0dv2
            !wl_430 = wl_430+vlop0*(vint_ci(list)-Two*vint_ci(list+1))
            wl = vlop0
            call trans_ijkl_intpos(lr,lrk,lr0,lrk,nxo)
            call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            wl = -vlop0*2
            call trans_ijkl_intpos(lr,lr0,lrk,lrk,nxo)
            call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

          end do
          !wl = wl+wl_430
          !*********************************************************************
          !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
        end do

      case (1:5,7,9:12,14,16:22,25:26)
    end select

  case (25)
    ! line=25:-d^r^r<-->sv(10),tv(17),ttv(18)
    select case (lpok)
      case default ! (10)
        call sv_head_dbl_tail_act_g(lra)

      case (17)
        ! tv(17) ar(23)-br(23)-
        iwdr = 0
        do lri=norb_frz+1,norb_dz
          imi = lsm_inn(lri)
          do lrj=lri,norb_dz
            imj = lsm_inn(lrj)
            imij = Mul(imi,imj)
            if (imij /= jml) cycle
            iwdl = just(lri,lrj)
            vlop1 = w1*w1_tv             !d17 vlop0=0
            !list = list3(lri,lrj,lra)
            !wl = vlop1*vint_ci(list)        !2.1                  !!!!!
            wl = vlop1
            call trans_ijkl_intpos(lrj,lra,lri,lra,nxo)
            call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
          end do
        end do

      case (18)
        ! t1v(18) ar(13)-br(13)-
        iwdr = 0
        do lri=norb_frz+1,norb_dz
          imi = lsm_inn(lri)
          do lrj=lri,norb_dz
            imj = lsm_inn(lrj)
            imij = Mul(imi,imj)
            if (imij /= jml) cycle
            iwdl = just(lri,lrj)
            vlop1 = w1*w1_t1v             !d18 vlop0=0
            !list = list3(lri,lrj,lra)
            !wl = vlop1*vint_ci(list)        !2.1                  !!!!!
            wl = vlop1
            call trans_ijkl_intpos(lrj,lra,lri,lra,nxo)
            call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

          end do
        end do

      case (1:9,11:16,19:26)
    end select

  case (26)
    ! line=26:-d^r^l<-->ss(1),st(2),ts(3),stt(4),tts(5),tt(11),tttt(12),dd(19
    select case (lpok)
      case default ! (1)
        call ss_head_dbl_tail_act_g(lra)

      case (2)
        call st_head_dbl_tail_act_g(lra)

      case (3)
        call ts_head_dbl_tail_act_g(lra)

      case (4)
        call stt_head_dbl_tail_act_g(lra)

      case (5)
        call tts_head_dbl_tail_act_g(lra)

      case (11)
        call tt_head_dbl_tail_act_g(lra)

      case (12)
        call tttt_head_dbl_tail_act_g(lra)

      case (19)
        call dd_head_dbl_tail_act_g(lra)

      case (20)
        call dddd_head_dbl_tail_act_g(lra)

      case (21)
        call dd1_head_dbl_tail_act_g(lra)

      case (22)
        call d1d_head_dbl_tail_act_g(lra)

      case (25)
        ! vv(25) drl(33)-
        if (jwl == jwr) return
        vlop0 = w0*w0_vv             !d25
        wl = Zero
        iwdl = 0
        iwdr = 0
        do lri=1,norb_dz
          !wl = wl+vlop0*voint(lri,lra)
          wl = vlop0
          call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,lri,lra)
        end do
        !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)

      case (6:10,13:18,23:24,26)
    end select

  case (27)
    ! line=27:-b^r-a^r<-->sv(10),tv(17),ttv(18)
    select case (lpok)
      case default ! (10)
        call sv_head_dbl_tail_act_g(lra)

      case (17)
        ! tv(17) ar(23)-br(23)-
        iwdr = 0
        do lri=norb_frz+1,norb_dz
          lmi = lsm_inn(lri)
          do lrj=lri+1,norb_dz
            lmj = lsm_inn(lrj)
            lmij = Mul(lmi,lmj)
            if (lmij /= jml) cycle
            w1tv = w1_tv
            if (mod(lrj-lri,2) == 0) w1tv = -w1tv
            iwdl = just(lri,lrj)
            vlop1 = w1*w1tv             !d17
            !list = list4(lri,lrj,lrs,lra)
            !wl = vlop1*(vint_ci(list)-vint_ci(list+2)) !1.3 vlop0=0      !!!!
            !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
            wl = vlop1
            call trans_ijkl_intpos(lra,lrj,lrs,lri,nxo)
            call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            wl = -vlop1
            call trans_ijkl_intpos(lra,lri,lrj,lrs,nxo)
            call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

          end do
        end do

      case (18)
        ! t1v(18) ar(13)-br(13)-
        iwdr = 0
        do lri=norb_frz+1,norb_dz-1
          lmi = lsm_inn(lri)
          do lrj=lri+1,norb_dz
            lmj = lsm_inn(lrj)
            lmij = Mul(lmi,lmj)
            if (lmij /= jml) cycle
            w1tv = w1_t1v
            if (mod(lrj-lri,2) == 0) w1tv = -w1tv
            iwdl = just(lri,lrj)
            vlop1 = w1*w1tv             !d17
            !list = list4(lri,lrj,lrs,lra)
            !wl = vlop1*(vint_ci(list)-vint_ci(list+2)) !1.3 vlop0=0      !!!!
            !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
            wl = vlop1
            call trans_ijkl_intpos(lra,lrj,lrs,lri,nxo)
            call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            wl = -vlop1
            call trans_ijkl_intpos(lra,lri,lrj,lrs,nxo)
            call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
          end do
        end do

      case (1:9,11:16,19:26)
    end select

  case (28)
    ! line=28:-b^l-a^r<-->ss(1),st(2),ts(3),stt(4),tts(5),tt(11),tttt(12),dd(
    select case (lpok)
      case default ! (1)
        call ss_head_dbl_tail_act_g(lra)

      case (2)
        call st_head_dbl_tail_act_g(lra)

      case (3)
        call ts_head_dbl_tail_act_g(lra)

      case (4)
        call stt_head_dbl_tail_act_g(lra)

      case (5)
        call tts_head_dbl_tail_act_g(lra)

      case (11)
        call tt_head_dbl_tail_act_g(lra)

      case (12)
        call tttt_head_dbl_tail_act_g(lra)

      case (19)
        call dd_head_dbl_tail_act_g(lra)

      case (20)
        call dddd_head_dbl_tail_act_g(lra)

      case (21)
        call dd1_head_dbl_tail_act_g(lra)

      case (22)
        call d1d_head_dbl_tail_act_g(lra)

      case (25)
        ! vv(25) drl(33)-
        vlop0 = w0*w0_vv             !d25
        wl = 0.d0
        iwdl = 0
        iwdr = 0
        do lrk=1,norb_dz
          !list = list3(lrs,lra,lrk)
          !wl = wl+vlop0*(vint_ci(list)-Two*vint_ci(list+1))   !4.3 vlop1=0
          wl = vlop0
          call trans_ijkl_intpos(lra,lrk,lrs,lrk,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
          wl = -vlop0*2
          call trans_ijkl_intpos(lra,lrs,lrk,lrk,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

        end do
        !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)

      case (6:10,13:18,23:24,26)
    end select

  case (29)
    ! line=29:-b^r-a^l<-->ss(1),st(2),ts(3),stt(4),tts(5),tt(11),tttt(12),dd(
    select case (lpok)
      case default ! (1)
        call ss_head_dbl_tail_act_g(lra)

      case (2)
        call st_head_dbl_tail_act_g(lra)

      case (3)
        call ts_head_dbl_tail_act_g(lra)

      case (4)
        call stt_head_dbl_tail_act_g(lra)

      case (5)
        call tts_head_dbl_tail_act_g(lra)

      case (11)
        call tt_head_dbl_tail_act_g(lra)

      case (12)
        call tttt_head_dbl_tail_act_g(lra)

      case (19)
        call dd_head_dbl_tail_act_g(lra)

      case (20)
        call dddd_head_dbl_tail_act_g(lra)

      case (21)
        call dd1_head_dbl_tail_act_g(lra)

      case (22)
        call d1d_head_dbl_tail_act_g(lra)

      case (25)
        ! vv(25) drl(33)-
        if (jwl >= jwr) return
        vlop0 = w0*w0_vv             !d25
        wl = Zero
        iwdl = 0
        iwdr = 0
        do lrk=1,norb_dz
          !list = list3(lrs,lra,lrk)
          !wl = vlop0*vint_ci(list)      !4.2  vlop1=0         !!!!!
          wl = vlop0
          call trans_ijkl_intpos(lra,lrk,lrs,lrk,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
        end do
        !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)

      case (6:10,13:18,23:24,26)
    end select

  case (30)
    ! line=30:-b&r-d^r^r<-->sd(6),sdd(8),td(13),ttdd(15),dv(23),ddv(24)
    select case (lpok)
      case default ! (6)
        ! sd(6-3) a&r(13)c'(22)-
        call dbl_sd_act_comp_g(3,lra)

      case (8)
        ! sd(6-3) tmp for spin=0
        !return
        call dbl_sdd_act_comp_g(3,lra)

      case (13)
        call dbl_td_act_comp_g(3,lra)

      case (15)
        call dbl_ttdd_act_comp_g(3,lra)

      case (23)
        ! dv(23-1) ar(23)-
        ! dv(23-2) drl(33)-bl(23)-
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
          !list = list3(lrd,lrg,lra)
          !wl = (vlop0+vlop1)*vint_ci(list)        !2.1       !!!!!
          !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
          wl = vlop0+vlop1
          call trans_ijkl_intpos(lrg,lra,lrd,lra,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
        end do

      case (24)
        ! d1v(24-1)  ar(13)-
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
          !list = list3(lrd,lrg,lra)
          !wl = (vlop0+vlop1)*vint_ci(list)        !2.1
          !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
          wl = vlop0+vlop1
          call trans_ijkl_intpos(lrg,lra,lrd,lra,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
        end do

      case (1:5,7,9:12,14,16:22,25:26)
    end select

  case (31)
    ! line=31:-b&l-d^r^l<-->sd(6),sdd(8),td(13),ttdd(15),dv(23),ddv(24)
    select case (lpok)
      case default ! (6)
        call dbl_sd_act_comp_g(5,lra)

      case (8)
        call dbl_sdd_act_comp_g(5,lra)

      case (13)
        call dbl_td_act_comp_g(5,lra)

      case (15)
        call dbl_ttdd_act_comp_g(5,lra)

      case (23)
        ! dv(23-1) ar(23)-
        ! dv(23-1) ar(23)-
        ! dv(23-2) drl(33)-bl(23)-
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
          !list = list3(lrd,lrg,lra)
          !wl = vlop0*(vint_ci(list)-Two*vint_ci(list+1))-vlop1*(vint_ci(list))   !2.2          !!!!
          !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
          wl = vlop0-vlop1
          call trans_ijkl_intpos(lrg,lra,lrd,lra,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
          wl = -vlop0*2
          call trans_ijkl_intpos(lrg,lrd,lra,lra,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
        end do

      case (24)
        ! d1v(24-1)  ar(13)-
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
          !list = list3(lrd,lrg,lra)
          !wl = vlop0*(vint_ci(list)-Two*vint_ci(list+1))-vlop1*(vint_ci(list))   !2.2          !!!!
          !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
          wl = vlop0-vlop1
          call trans_ijkl_intpos(lrg,lra,lrd,lra,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
          wl = -vlop0*2
          call trans_ijkl_intpos(lrg,lrd,lra,lra,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

        end do

      case (1:5,7,9:12,14,16:22,25:26)
    end select

  case (32)
    ! line=32:-b&r-b^r-a^r<-->sd(6),sdd(8),td(13),ttdd(15),dv(23),ddv(24)
    select case (lpok)
      case default ! (6)
        call dbl_sd_act_comp_g(4,lra)

      case (8)
        call dbl_sdd_act_comp_g(4,lra)

      case (13)
        call dbl_td_act_comp_g(4,lra)

      case (15)
        call dbl_ttdd_act_comp_g(4,lra)

      case (23)
        ! dv(23-1) ar(23)-
        ! dv(23-2) drl(33)-bl(23)-
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
          !list = list4(lrd,lrg,lrs,lra)
          !wl = vlop0*(vint_ci(list+2)+vint_ci(list))-vlop1*(vint_ci(list+2)-vint_ci(list))  !1.3        !!!!!
          !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
          wl = vlop0+vlop1
          call trans_ijkl_intpos(lra,lrg,lrs,lrd,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
          wl = vlop0-vlop1
          call trans_ijkl_intpos(lra,lrd,lrg,lrs,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

        end do

      case (24)
        ! d1v(24-1)  ar(13)-
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
          !list = list4(lrd,lrg,lrs,lra)
          !wl = vlop0*(vint_ci(list+2)+vint_ci(list))-vlop1*(vint_ci(list+2)-vint_ci(list))  !1.3        !!!!!
          !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
          wl = vlop0+vlop1
          call trans_ijkl_intpos(lra,lrg,lrs,lrd,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
          wl = vlop0-vlop1
          call trans_ijkl_intpos(lra,lrd,lrg,lrs,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
        end do

      case (1:5,7,9:12,14,16:22,25:26)
    end select

  case (33)
    ! line=33:-b&l-b^r-a^l<-->sd(6),sdd(8),td(13),ttdd(15),dv(23),ddv(24)
    select case (lpok)
      case default ! (6)
        call dbl_sd_act_comp_g(6,lra)

      case (8)
        call dbl_sdd_act_comp_g(6,lra)

      case (13)
        call dbl_td_act_comp_g(6,lra)

      case (15)
        call dbl_ttdd_act_comp_g(6,lra)

      case (23)
        ! dv(23-1) ar(23)-
        ! dv(23-2) drl(33)-bl(23)-
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
          !list = list4(lrd,lrg,lrs,lra)
          !wl = vlop0*(vint_ci(list)-Two*vint_ci(list+1))-vlop1*vint_ci(list)  !1.1          !!!!!
          !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
          wl = vlop0-vlop1
          call trans_ijkl_intpos(lra,lrg,lrs,lrd,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
          wl = -vlop0*2
          call trans_ijkl_intpos(lra,lrs,lrg,lrd,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
        end do

      case (24)
        ! d1v(24-1)  ar(13)-
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
          !list = list4(lrd,lrg,lrs,lra)
          !wl = vlop0*(vint_ci(list)-Two*vint_ci(list+1))-vlop1*vint_ci(list)  !1.1          !!!!!
          !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
          wl = vlop0-vlop1
          call trans_ijkl_intpos(lra,lrg,lrs,lrd,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
          wl = -vlop0*2
          call trans_ijkl_intpos(lra,lrs,lrg,lrd,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
        end do

      case (1:5,7,9:12,14,16:22,25:26)
    end select

  case (34)
    ! line=34:-b&l-b^l-a^r<-->sd(6),sdd(8),td(13),ttdd(15),dv(23),ddv(24)
    select case (lpok)
      case default ! (6)
        call dbl_sd_act_comp_g(7,lra)

      case (8)
        call dbl_sdd_act_comp_g(7,lra)

      case (13)
        call dbl_td_act_comp_g(7,lra)

      case (15)
        call dbl_ttdd_act_comp_g(7,lra)

      case (23)
        ! dv(23-1) ar(23)-
        ! dv(23-2) drl(33)-bl(23)-
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
          !list = list4(lrd,lrg,lrs,lra)
          !wl = vlop0*(vint_ci(list+2)-Two*vint_ci(list+1))-vlop1*vint_ci(list+2) !1.2      !!!!
          !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
          wl = vlop0-vlop1
          call trans_ijkl_intpos(lra,lrd,lrg,lrs,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
          wl = -vlop0*2
          call trans_ijkl_intpos(lra,lrs,lrg,lrd,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

        end do

      case (24)
        ! d1v(24-1)  ar(13)-
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
          !list = list4(lrd,lrg,lrs,lra)
          !wl = vlop0*(vint_ci(list+2)-Two*vint_ci(list+1))-vlop1*vint_ci(list+2) !1.2      !!!!
          !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
          call trans_ijkl_intpos(lra,lrd,lrg,lrs,nxo)
          wl = vlop0-vlop1
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
          wl = -vlop0*2
          call trans_ijkl_intpos(lra,lrs,lrg,lrd,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
        end do

      case (1:5,7,9:12,14,16:22,25:26)
    end select

  case (35)
    ! line=35:-d&r^l-a^l<-->sd(6),sdd(8),td(13),ttdd(15),dv(23),ddv(24)
    select case (lpok)
      case default ! (6)
        call dbl_sd_act_comp_g(2,lra)

      case (8)
        call dbl_sdd_act_comp_g(2,lra)

      case (13)
        call dbl_td_act_comp_g(2,lra)

      case (15)
        call dbl_ttdd_act_comp_g(2,lra)

      case (23)
        ! dv(23-1) ar(23)-
        ! dv(23-2) drl(33)-bl(23)-
        iwdr = 0
        do lrd=norb_frz+1,norb_dz
          imd = lsm_inn(lrd)
          if (imd /= jml) cycle
          iwdl = jud(lrd)
          w0dv1 = w0_dv(1)
          ni = mod(norb_dz-lrd,2)
          if (ni == 1) w0dv1 = -w0dv1
          vlop0 = w0*w0dv1                !d23-1
          !list = list3(lrd,lra,lrs)
          !wl = vlop0*vint_ci(list)          !3.2         !!!!!
          !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
          wl = vlop0
          call trans_ijkl_intpos(lra,lrs,lrd,lrs,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
        end do

      case (24)
        ! d1v(24-1)  ar(13)-
        iwdr = 0
        do lrd=norb_frz+1,norb_dz
          imd = lsm_inn(lrd)
          if (imd /= jml) cycle
          iwdl = jud(lrd)
          w0dv1 = w0_d1v(1)
          ni = mod(norb_dz-lrd,2)
          if (ni == 1) w0dv1 = -w0dv1
          vlop0 = w0*w0dv1                !d23-1
          !list = list3(lrd,lra,lrg)
          !wl = vlop0*vint_ci(list)          !3.2         !!!!!
          !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
          wl = vlop0
          call trans_ijkl_intpos(lra,lrg,lrd,lrg,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
        end do

      case (1:5,7,9:12,14,16:22,25:26)
    end select
end select

return

end subroutine dbl_head_act_tail_g

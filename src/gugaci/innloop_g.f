************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
!****************************************************************
!   ar      =1 (+a^r)     drr     =2 (+d^rr)   drl     =3 (+d^rl)
!   arbr    =4 (+d^rr)    arbl    =5 (+d^rl)   ard_l^r =6 (+a^l)
!   drrb^r  =7 (+a^r)     drlb^r  =8 (+a^l)    drlb^l  =9 (+a^r)
!   arbrb^r =10 (+a^r)    arblb^r =11 (+a^l)   arblb^l =12 (+a^r)
!   drl        =13 (*)
!****************************************************************

      subroutine inner_space_loop_g()
      implicit REAL*8 (a-h,o-z)
c      wsc0=c_time()
      call dbl_space_loop_g()
c      wsc1=c_time()
      call act_space_cloop_g()
      call act_space_ploop_g()
c      wsc2=c_time()
c      write(6,'(2x,2(a5,f12.4))')'dbl',wsc1-wsc0,'act',wsc2-wsc1
      return
      end

      subroutine act_space_cloop_g()        ! one sub_drt
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      if(norb_act.eq.0) return
      do ipae_=1,25
        ipae=ipae_ ! ipae is in common block, is this necessary?
        jpae=nu_ae(ipae)
        if(jpae.eq.0) cycle
        do jpad_=1,mxnode
          jpad=jpad_ ! jpad is in common block, is this necessary?
          if(nu_ad(jpad).eq.0) cycle
          call seg_drt()
          if(ndim .eq. 0) cycle
          call copy_to_drtl()
          call cloop_in_act_g()
        enddo
      enddo
      return
      end

      subroutine act_space_ploop_g()        ! two sub_drt with same ipae
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"

      if(norb_act.eq.0) return
      do ipae_=1,25
        ipae=ipae_ ! ipae is in common block, is this necessary?
        jpae=nu_ae(ipae)
        if(jpae.eq.0) cycle
        do jpadl_=1,mxnode                               ! jpadl
          jpadl=jpadl_ ! jpadl is in common block, is this necessary?
          if(nu_ad(jpadl).eq.0) cycle
            jpad=jpadl
            call seg_drt()
            if(ndim .eq. 0) cycle
            call copy_to_drtl()
          do jpad_=1,mxnode                               !jpadr
            jpad=jpad_ ! jpad is in common block, is this necessary?
            if(nu_ad(jpad).eq.0) cycle
            call seg_drt()
            if(ndim .eq. 0) cycle
c        if( ipae.eq.18.and. jpadl.eq.2.and.jpad.eq.1) then
c      write(6,*)
c     endif
           call ploop_in_act_g()
          enddo
        enddo
      enddo
      return
      end

      subroutine cloop_in_act_g()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      do lrai=norb_dz+1,norb_inn-1
        lmi=lsm_inn(lrai)
        do lraj=lrai+1,norb_inn
          lmj=lsm_inn(lraj)
          lmij=mul_tab(lmi,lmj)
c-----------------------------------------------------------
!line=8 d&r&r--d^r^r
          call head_drr_at_given_orb(mh,lrai)
          logic_br(1:mh)=.true.
          call link_c2_to_given_orb(mh,lrai+1,lraj-1)
          call tail_drr_at_given_orb(mh,lraj)
c          write(6,'(6i6)')8,mh,lrai,lraj,0,0
          if(mh.ne.0) call act_cloop_g(8,mh,lrai,lraj,0,0)
c-----------------------------------------------------------
!line=9 d&r&l--d^r^l
          call head_drl_at_given_orb(mh,lrai)
          call link_c2_to_given_orb(mh,lrai+1,lraj-1)
          call tail_drl_at_given_orb(mh,lraj)
c          write(6,'(6i6)')9,mh,lrai,lraj,0,0
          if(mh.ne.0) call act_cloop_g(9,mh,lrai,lraj,0,0)
          lsmij=mul_tab(lmi,lmj)
c-----------------------------------------------------------
          if(lsmij.ne.1) goto 400
c-----------------------------------------------------------
!line=1 a&r--a^r
          call head_ar_at_given_orb(mh,lrai)
          call link_c1_to_given_orb_coe(mh,lrai+1,lraj-1)
          call tail_ar_at_given_orb_coe(mh,lraj)
c          write(6,'(6i6)')1,mh,lrai,lraj,0,0
          if(mh.ne.0) call act_cloop_g(1,mh,lrai,lraj,0,0)
c          call save_clp(1,mh,lra,0)
c-----------------------------------------------------------
!line=2 a&r-d^r&l-a^l
          do lrak=lrai+1,lraj-1
            call head_ar_at_given_orb(mh,lrai)
            call link_c1_to_given_orb(mh,lrai+1,lrak-1)
            call link_d10_at_given_orb(mh,lrak)
            call link_c1_to_given_orb(mh,lrak+1,lraj-1)
            call tail_al_at_given_orb(mh,lraj)
c          write(6,'(6i6)')2,mh,lrai,lraj,lrak,0
            if(mh.ne.0) call act_cloop_g(2,mh,lrai,lraj,lrak,0)
          enddo
c-----------------------------------------------------------
!line=3 a&r-b&r-d^r^r
          do lrak=lraj+1,norb_inn
           call head_ar_at_given_orb(mh,lrai)
            call link_c1_to_given_orb(mh,lrai+1,lraj-1)
            call link_b4_at_given_orb(mh,lraj)
            logic_br(1:mh)=.true.
            call link_c2_to_given_orb(mh,lraj+1,lrak-1)
            call tail_drr_at_given_orb(mh,lrak)
c          write(6,'(6i6)')3,mh,lrai,lraj,lrak,0
            if(mh.ne.0) call act_cloop_g(3,mh,lrai,lraj,lrak,0)
c-----------------------------------------------------------
!line=5 a&r-b&l-d^r^l
            call head_ar_at_given_orb(mh,lrai)
            call link_c1_to_given_orb(mh,lrai+1,lraj-1)
            call link_b3_at_given_orb(mh,lraj)
            logic_br(1:mh)=.true.
            call link_c2_to_given_orb(mh,lraj+1,lrak-1)
            call tail_drl_at_given_orb(mh,lrak)
c          write(6,'(6i6)')5,mh,lrai,lraj,lrak,0
            if(mh.ne.0) call act_cloop_g(5,mh,lrai,lraj,lrak,0)
          enddo
c-----------------------------------------------------------
          do lrak=norb_dz+1,lrai-1
!line=10 d&rr--b^r--a^r
            call head_drr_at_given_orb(mh,lrak)
            logic_br(1:mh)=.true.
            call link_c2_to_given_orb(mh,lrak+1,lrai-1)
            call link_b2_at_given_orb(mh,lrai)
            call link_c1_to_given_orb(mh,lrai+1,lraj-1)
            call tail_ar_at_given_orb(mh,lraj)
c          write(6,'(6i6)')10,mh,lrai,lraj,lrak,0
            if(mh.ne.0) call act_cloop_g(10,mh,lrai,lraj,lrak,0)
c-----------------------------------------------------------
!line=11 d&r&l-b^r-a^l
           call head_drl_at_given_orb(mh,lrak)
            call link_c2_to_given_orb(mh,lrak+1,lrai-1)
            call link_b2_at_given_orb(mh,lra)
            call link_c1_to_given_orb(mh,lrai+1,lraj-1)
            call tail_al_at_given_orb(mh,lraj)
c          write(6,'(6i6)')11,mh,lrai,lraj,lrak,0
            if(mh.ne.0) call act_cloop_g(11,mh,lrai,lraj,lrak,0)
c-----------------------------------------------------------
!line=12 d&r&l-b^l-a^r
           call head_drl_at_given_orb(mh,lrak)
            call link_c2_to_given_orb(mh,lrak+1,lrai-1)
            call link_b1_at_given_orb(mh,lrai)
            call link_c1_to_given_orb(mh,lrai+1,lraj-1)
            call tail_ar_at_given_orb(mh,lraj)
c          write(6,'(6i6)')12,mh,lrai,lraj,lrak,0
            if(mh.ne.0) call act_cloop_g(12,mh,lrai,lraj,lrak,0)
          enddo
400       if(lraj.gt.norb_inn-2) cycle
         do lrak=lraj+1,norb_inn
            lmk=lsm_inn(lrak)
            lmk=mul_tab(lmij,lmk)
           do lral=lrak+1,norb_inn
              lml=lsm_inn(lral)
              lml=mul_tab(lmk,lml)
            if(lml.ne.1) cycle
!line=4  a&r--b&r--b^r--a^r
              call head_ar_at_given_orb(mh,lrai)
              call link_c1_to_given_orb(mh,lrai+1,lraj-1)
              call link_b4_at_given_orb(mh,lraj)
              logic_br(1:mh)=.true.
              call link_c2_to_given_orb(mh,lraj+1,lrak-1)
              call link_b2_at_given_orb(mh,lrak)
              call link_c1_to_given_orb(mh,lrak+1,lral-1)
              call tail_ar_at_given_orb(mh,lral)
c          write(6,'(6i6)')4,mh,lrai,lral,lraj,lrak
              if(mh.ne.0) call act_cloop_g(4,mh,lrai,lral,lraj,lrak)
c-----------------------------------------------------------
!line=6  a&r--b&l--b^r--a^l
              call head_ar_at_given_orb(mh,lrai)
              call link_c1_to_given_orb(mh,lrai+1,lraj-1)
              call link_b3_at_given_orb(mh,lraj)
              logic_br(1:mh)=.true.
              call link_c2_to_given_orb(mh,lraj+1,lrak-1)
              call link_b2_at_given_orb(mh,lrak)
              call link_c1_to_given_orb(mh,lrak+1,lral-1)
              call tail_al_at_given_orb(mh,lral)
c          write(6,'(6i6)')6,mh,lrai,lral,lraj,lrak
              if(mh.ne.0) call act_cloop_g(6,mh,lrai,lral,lraj,lrak)
c-----------------------------------------------------------
!line=7 a&r--b&l--b^l--a^r
              call head_ar_at_given_orb(mh,lrai)
              call link_c1_to_given_orb(mh,lrai+1,lraj-1)
              call link_b3_at_given_orb(mh,lraj)
              logic_br(1:mh)=.true.
              call link_c2_to_given_orb(mh,lraj+1,lrak-1)
              call link_b1_at_given_orb(mh,lrak)
              call link_c1_to_given_orb(mh,lrak+1,lral-1)
              call tail_ar_at_given_orb(mh,lral)
c          write(6,'(6i6)')7,mh,lrai,lral,lraj,lrak
              if(mh.ne.0) call act_cloop_g(7,mh,lrai,lral,lraj,lrak)
c-----------------------------------------------------------
            enddo
          enddo
        enddo
      enddo
      return
      end

      subroutine ploop_in_act_g()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
c===========================================================
      do lrai=norb_dz+1,norb_inn
!line=25 -c"-d^r^r
        logic_br(1)=.true.
        call link_c2_to_given_orb(mh,norb_dz+1,lrai-1)
        call tail_drr_at_given_orb(mh,lrai)
        if(mh.ne.0)  call lp_act_tail_g(25,mh,0,lrai)
c-----------------------------------------------------------
!line=26 -c"-d^r^l
        logic_br(1)=.true.
        call link_c2_to_given_orb(mh,norb_dz+1,lrai-1)
        call tail_drl_at_given_orb(mh,lrai)
        if(mh.ne.0)  call lp_act_tail_g(26,mh,0,lrai)
c===========================================================
!line=23 -c'-a^l
        call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
        call tail_al_at_given_orb(mh,lrai)
        if(mh.ne.0)  call lp_act_tail_g(23,mh,0,0)
c-----------------------------------------------------------
!line=24 -c'-a^r
        call link_c1_to_given_orb_coe(mh,norb_dz+1,lrai-1)
        call tail_ar_at_given_orb_coe(mh,lrai)
        if(mh.ne.0)  call lp_act_tail_g(24,mh,0,0)
c===========================================================
        do lrak=lrai+1,norb_inn
c-----------------------------------------------------------
!line=30 -c'-b&r-d^r^r
          call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
          call link_b4_at_given_orb(mh,lrai)
          logic_br(1:mh)=.true.
          call link_c2_to_given_orb(mh,lrai+1,lrak-1)
          call tail_drr_at_given_orb(mh,lrak)
          if(mh.ne.0)  call lp_act_tail_g(30,mh,lrai,0)
c-----------------------------------------------------------
!line=31 -c'-b&l-d^r^l
          call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
          call link_b3_at_given_orb(mh,lrai)
          logic_br(1:mh)=.true.
          call link_c2_to_given_orb(mh,lrai+1,lrak-1)
          call tail_drl_at_given_orb(mh,lrak)
          if(mh.ne.0)  call lp_act_tail_g(31,mh,lrai,0)
        enddo
c===========================================================
        do lrak=norb_dz+1,lrai-1
c-----------------------------------------------------------
!line=35 -c'-d^r&l-a^l
          call link_c1_to_given_orb(mh,norb_dz+1,lrak-1)
          call link_d10_at_given_orb(mh,lrak)
          call link_c1_to_given_orb(mh,lrak+1,lrai-1)
          call tail_al_at_given_orb(mh,lrai)
          if(mh.ne.0)  call lp_act_tail_g(35,mh,lrak,lrak)
         enddo
c===========================================================
        do lraj=lrai+1,norb_inn
c-----------------------------------------------------------
!line=27 -c"-b^r-a^r
          call link_c2_to_given_orb(mh,norb_dz+1,lrai-1)
          call link_b2_at_given_orb(mh,lrai)
          call link_c1_to_given_orb(mh,lrai+1,lraj-1)
          call tail_ar_at_given_orb(mh,lraj)
          if(mh.ne.0)  call lp_act_tail_g(27,mh,0,lrai)
c-----------------------------------------------------------
!line=28 -c"-b^l-a^r
          call link_c2_to_given_orb(mh,norb_dz+1,lrai-1)
          call link_b1_at_given_orb(mh,lrai)
          call link_c1_to_given_orb(mh,lrai+1,lraj-1)
          call tail_ar_at_given_orb(mh,lraj)
          if(mh.ne.0)  call lp_act_tail_g(28,mh,0,lrai)
c-----------------------------------------------------------
!line=29 -c"-b^r-a^l
          call link_c2_to_given_orb(mh,norb_dz+1,lrai-1)
          call link_b2_at_given_orb(mh,lrai)
          call link_c1_to_given_orb(mh,lrai+1,lraj-1)
          call tail_al_at_given_orb(mh,lraj)
          if(mh.ne.0)  call lp_act_tail_g(29,mh,0,lrai)
c-----------------------------------------------------------
          do lral=lraj+1,norb_inn
c-----------------------------------------------------------
!line=32 -c'-b&r-c"-b^r-a^r
            call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
            call link_b4_at_given_orb(mh,lrai)
            logic_br(1:mh)=.true.
            call link_c2_to_given_orb(mh,lrai+1,lraj-1)
            call link_b2_at_given_orb(mh,lraj)
            call link_c1_to_given_orb(mh,lraj+1,lral-1)
            call tail_ar_at_given_orb(mh,lral)
           if(mh.ne.0)  call lp_act_tail_g(32,mh,lrai,lraj)
c-----------------------------------------------------------
!line=33 -c'-b&l-c"-b^r-a^l
            call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
            call link_b3_at_given_orb(mh,lrai)
            logic_br(1:mh)=.true.
            call link_c2_to_given_orb(mh,lrai+1,lraj-1)
            call link_b2_at_given_orb(mh,lraj)
            call link_c1_to_given_orb(mh,lraj+1,lral-1)
            call tail_al_at_given_orb(mh,lral)
           if(mh.ne.0)  call lp_act_tail_g(33,mh,lrai,lraj)
c-----------------------------------------------------------
!line=34 -c'-b&l-c"-b^l-a^r
            call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
            call link_b3_at_given_orb(mh,lrai)
            logic_br(1:mh)=.true.
            call link_c2_to_given_orb(mh,lrai+1,lraj-1)
            call link_b1_at_given_orb(mh,lraj)
            call link_c1_to_given_orb(mh,lraj+1,lral-1)
            call tail_ar_at_given_orb(mh,lral)
            if(mh.ne.0)  call lp_act_tail_g(34,mh,lrai,lraj)
          enddo
        enddo
c-----------------------------------------------------------
      enddo
      return
      end

      subroutine lp_act_tail_g(lin,mh,lrg0,lrs0)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      dimension lpcoe(norb_dz+1:norb_inn)
      common/onepl/line,jph,jpel,jper,lrg,lrs,jwl,jwr,w0,w1
      line=lin
      lrg=lrg0
      lrs=lrs0
      jph=0
      do mhlp_=1,mh
        mhlp=mhlp_ ! mhlp is in common block, is this necessary?
        jpel=lpnew_ltail(mhlp)
        jper=lpnew_rtail(mhlp)
        jwl =lpnew_lwei(mhlp)
        jwr =lpnew_rwei(mhlp)
          w0  =vplpnew_w0(mhlp)
          w1  =vplpnew_w1(mhlp)
        if(line.eq.24) then
          do iorb=norb_dz+1,norb_inn
            lpcoe(iorb)=lpnew_coe(iorb,mhlp)
          enddo
        endif
        call dbl_head_act_tail_g(lpcoe)
      enddo
      return
      end

      subroutine act_cloop_g(lin,mh,lr0,lr,lrg0,lrs0)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      dimension lpcoe(norb_dz+1:norb_inn)
      common/onepl/line,jph,jpel,jper,lrg,lrs,jwl,jwr,w0,w1
      parameter (two=2.0d+00,half=0.5d+00)
      line=lin
      lrg=lrg0
      lrs=lrs0
      do mhlp_=1,mh
        mhlp=mhlp_ ! mhlp is in common block, is this necessary?
        jph =lpnew_head(mhlp)
        jpel=lpnew_ltail(mhlp)
        jper=lpnew_rtail(mhlp)
        jwl =lpnew_lwei(mhlp)
        jwr =lpnew_rwei(mhlp)
        vlop0 =vplpnew_w0(mhlp)
        vlop1 =vplpnew_w1(mhlp)

        goto (31,32,21,13,22,11,12,51,52,41,42,43),line
c-----------------------------------------------------------
31          do l=norb_dz+1,lr
              lpcoe(l)=lpnew_coe(l,mhlp)
            enddo
c            wl=voint(lr0,lr)
            wl=vlop0
c            write(nf2,*) 'ar-ar'
          call prodab_1(2,jph,jpel,jwl,jwr,0,wl,jper,lr0,lr)
            do l=lr0,lr
c              list=list3(lr0,lr,l)
              kcoe=lpcoe(l)
c            write(nf2,*) 'arw-c-arw'
              call neoc(kcoe,nocc,tcoe)
c              wl=wl+nocc*(vint_ci(list+1)+tcoe*vint_ci(list))
              wl=vlop0*nocc
          call trans_ijkl_intpos(lr,lr0,l,l,nxo)
c          write(nf2,'(4i4)') lr,lr0,l,l
          call prodab_2(2,jph,jpel,jwl,jwr,0,wl,jper,nxo)
             wl=vlop0*nocc*tcoe
          call trans_ijkl_intpos(lr,l,lr0,l,nxo)
          call prodab_2(2,jph,jpel,jwl,jwr,0,wl,jper,nxo)

            enddo
c            wl=wl*vlop0
            goto 500
c-----------------------------------------------------------
32         wl= vlop0
c         write(nf2,*) 'ar-drl-al'
c        list=list3(lr0,lr,lrg)
c            wl=vlop0*vint_ci(list)
          call trans_ijkl_intpos(lr,lrg,lr0,lrg,nxo)
          call prodab_2(2,jph,jpel,jwl,jwr,0,wl,jper,nxo)
          goto 500
c-----------------------------------------------------------
11        wl=vlop0-vlop1
c           list=list4(lr0,lrg,lrs,lr)
c             wl=vlop0*(vint_ci(list)-2*vint_ci(list+1))
c     :             -vlop1*vint_ci(list)
c        vint_ci(numb)=vfutei(la,lc,lb,ld)
c        vint_ci(numb+1)=vfutei(la,lb,lc,ld)
c        vint_ci(numb+2)=vfutei(la,ld,lc,lb)
c          write(nf2,*)  'ar-bl-br-al'
          call trans_ijkl_intpos(lr,lrg,lrs,lr0,nxo)
          call prodab_2(2,jph,jpel,jwl,jwr,0,wl,jper,nxo)
               wl=-vlop0*2.0d0
          call trans_ijkl_intpos(lr,lrs,lrg,lr0,nxo)
          call prodab_2(2,jph,jpel,jwl,jwr,0,wl,jper,nxo)
          goto 500

c-----------------------------------------------------------
12        wl=vlop0-vlop1
c      list=list4(lr0,lrg,lrs,lr)
c             wl=vlop0*(vint_ci(list+2)-2.0d0*vint_ci(list+1))
c     :           -vlop1*vint_ci(list+2)
c          write(nf2,*) 'ar-bl-bl-ar'
          call trans_ijkl_intpos(lr,lr0,lrg,lrs,nxo)
          call prodab_2(2,jph,jpel,jwl,jwr,0,wl,jper,nxo)
               wl=-vlop0*2.0d0
          call trans_ijkl_intpos(lr,lrs,lrg,lr0,nxo)
          call prodab_2(2,jph,jpel,jwl,jwr,0,wl,jper,nxo)
          goto 500

c-----------------------------------------------------------
13       wl=vlop0+vlop1
c         write(nf2,*) 'ar-br-br-ar'
c           list=list4(lr0,lrg,lrs,lr)
c             wl=vlop0*(vint_ci(list+2)+vint_ci(list))
c     :             -vlop1*(vint_ci(list+2)-vint_ci(list))

          call trans_ijkl_intpos(lr,lrg,lrs,lr0,nxo)
          call prodab_2(2,jph,jpel,jwl,jwr,0,wl,jper,nxo)
              wl=vlop0-vlop1
          call trans_ijkl_intpos(lr,lr0,lrg,lrs,nxo)
          call prodab_2(2,jph,jpel,jwl,jwr,0,wl,jper,nxo)
          goto 500
c-----------------------------------------------------------
21       wl=vlop0+vlop1
c         write(nf2,*)  'ar-br-drr'
c         list=list3(lr0,lr,lrg)
c             wl=(vlop0+vlop1)*vint_ci(list)

          call trans_ijkl_intpos(lr,lrg,lr0,lrg,nxo)
          call prodab_2(2,jph,jpel,jwl,jwr,0,wl,jper,nxo)
          goto 500
c-----------------------------------------------------------
22       wl=vlop0-vlop1
c         write(nf2,*) 'ar-bl-drl'
c          list=list3(lr0,lr,lrg)
c             wl=vlop0*(vint_ci(list)-2.0d0*vint_ci(list+1))
c     :           -vlop1*vint_ci(list)

          call trans_ijkl_intpos(lr,lrg,lr0,lrg,nxo)
          call prodab_2(2,jph,jpel,jwl,jwr,0,wl,jper,nxo)
               wl=-vlop0*2.0d0
          call trans_ijkl_intpos(lr,lr0,lrg,lrg,nxo)
          call prodab_2(2,jph,jpel,jwl,jwr,0,wl,jper,nxo)
          goto 500
c-----------------------------------------------------------
c  drr-drr
c  wl=vlop0 but not 0.5d0*vlop0 is based on that the non-diagonal
c  just uses the non-triangle <ci|h|cj> which designates that i > j.
c
51           wl=vlop0

c          wl=vlop0*voint(lr,lr0)*0.5d0
          call trans_ijkl_intpos(lr,lr0,lr,lr0,nxo)
c          write(nf2,*) 'drr-drr'
          call prodab_2(2,jph,jpel,jwl,jwr,0,wl,jper,nxo)
          goto 500
c-----------------------------------------------------------
52        wl=(vlop0-vlop1)*2.0d0

c drl-drl
c          wl=(vlop0-vlop1)*voint(lr,lr0)
          call trans_ijkl_intpos(lr,lr0,lr,lr0,nxo)
c          write(nf2,'(a7,2i4,f18.10)') 'drl-drl',lr,lr0,wl
          call prodab_2(2,jph,jpel,jwl,jwr,0,wl,jper,nxo)
          goto 500
c-----------------------------------------------------------
41        wl=vlop0+vlop1
c          write(nf2,*)  'drr-br-ar'
c           list=list3(lr0,lr,lrg)
c              wl=(vlop0+vlop1)*vint_ci(list)

          call trans_ijkl_intpos(lr,lrg,lr0,lrg,nxo)
          call prodab_2(2,jph,jpel,jwl,jwr,0,wl,jper,nxo)
          goto 500
c-----------------------------------------------------------
42        wl=vlop0-vlop1
c          write(nf2,*) 'drl-br-al'
c          list=list3(lr0,lr,lrg)
c              wl=(vlop0-vlop1)*vint_ci(list)

          call trans_ijkl_intpos(lr,lrg,lr0,lrg,nxo)
          call prodab_2(2,jph,jpel,jwl,jwr,0,wl,jper,nxo)
          goto 500
c-----------------------------------------------------------
43       wl=vlop0-vlop1
c         write(nf2,*) 'drl-bl-ar'
c           list=list3(lr0,lr,lrg)
c              wl=vlop0*(vint_ci(list)-2.0d0*vint_ci(list+1))
c     :           -vlop1*vint_ci(list)

          call trans_ijkl_intpos(lr,lrg,lr0,lrg,nxo)
          call prodab_2(2,jph,jpel,jwl,jwr,0,wl,jper,nxo)
               wl=-vlop0*2.0d0
          call trans_ijkl_intpos(lr,lr0,lrg,lrg,nxo)
          call prodab_2(2,jph,jpel,jwl,jwr,0,wl,jper,nxo)
500     continue
c-----------------------------------------------------------
      enddo
      return
      end

      subroutine dbl_td_act_comp_g(lin,lra)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      common/onepl/line,jph,jpel,jper,lrg,lrs,jwl,jwr,w0,w1
!td(13-1) (22)a&(23)
!td(13-1) a&(23)c'(22)
!td(13-2) a&(23)b&r(23)b^r(32)
!td(13-3) a&(23)b&l(32)b^l(23)
!td(13-4) d&r&l(22)b^l(23)
!td(13-5) (22)d&&l(33)b^l(23)
!td(13-5) d&rl(33)c"(22)b^l(23)
!td(13-5) d&rl(33)b^l(23)c'(22)

      jmlr=mul_tab(jml,jmr)
        do lri=norb_frz+1,norb_dz
          lmi=lsm_inn(lri)
          if(lmi.ne.jmlr) cycle
          w0td1=w0_td(1)
          ni=mod(norb_dz-lri,2)
          if(ni.eq.1) then
            w0td1=-w0td1
          endif

!td(13-1) a&(23)c'(22)
          vlop0=-w0*w0td1
          vlop1=-w1*w0td1
          call comp_loop_g(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl0,list0,
     :                     wl1,list1)
          do lrk=lri+1,norb_dz
            lmk=lsm_inn(lrk)
            if(lmk.ne.jmr) cycle
            iwdl=just(lri,lrk)   !
            iwdr=jud(lrk)
c            call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
c           call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl2,jper,list2,list3)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl0,jper,list0)
      if(list1.ne.0)
     :        call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl1,jper,list1)
          enddo
!td(13-1) (22)a&(23)
          vlop0=w0*w0td1
          vlop1=w1*w0td1
          call comp_loop_g(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl0,list0,
     :                     wl1,list1)
          do lrk=norb_frz+1,lri-1
            lmk=lsm_inn(lrk)
            if(lmk.ne.jmr) cycle
            iwdl=just(lrk,lri)   !
            iwdr=jud(lrk)
c            call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
c           call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl2,jper,list2,list3)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl0,jper,list0)
      if(list1.ne.0)
     :        call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl1,jper,list1)
          enddo
!-------------------------------------------------------------------
        enddo
      return
        end

      subroutine dbl_ttdd_act_comp_g(lin,lra)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      common/onepl/line,jph,jpel,jper,lrg,lrs,jwl,jwr,w0,w1
!t1d1(15-1)  ar(13)-
!t1d1(15-1)  ar(13)-c'(11)-

      jmlr=mul_tab(jml,jmr)
        do lri=norb_frz+1,norb_dz
          lmi=lsm_inn(lri)
          if(lmi.ne.jmlr) cycle
          w0td1=w0_t1d1(1)
          ni=mod(norb_dz-lri,2)
          if(ni.eq.1) then
            w0td1=-w0td1
          endif
!t1d1(15-1)  ar(13)-c'(11)-
          vlop0=-w0*w0td1
          vlop1=-w1*w0td1
          call comp_loop_g(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl0,list0,
     :                     wl1,list1)
          do lrk=lri+1,norb_dz
            lmk=lsm_inn(lrk)
            if(lmk.ne.jmr) cycle
            iwdl=just(lri,lrk)   !
            iwdr=jud(lrk)
c            call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
c           call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl2,jper,list2,list3)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl0,jper,list0)
      if(list1.ne.0)
     :        call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl1,jper,list1)
          enddo
!t1d1(15-1)  (11)ar(13)-
          vlop0=w0*w0td1
          vlop1=w1*w0td1
          call comp_loop_g(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl0,list0,
     :                     wl1,list1)
          do lrk=norb_frz+1,lri-1
            lmk=lsm_inn(lrk)
            if(lmk.ne.jmr) cycle
            iwdl=just(lrk,lri)   !
            iwdr=jud(lrk)
c            call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
c           call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl2,jper,list2,list3)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl0,jper,list0)
      if(list1.ne.0)
     :        call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl1,jper,list1)
          enddo
!-------------------------------------------------------------------
        enddo
      return
        end

      subroutine comp_loop_g(line,lr0,lrg,lrs,lr,vlop0,vlop1,wl0,list0,
     :                     wl1,list1)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
        wl0=0.0d0
        list0=0
        wl1=0.0d0
        list1=0
      goto (31,32,21,13,22,11,12,51,52,41,42,43),line
c-----------------------------------------------------------
31    goto 500
c         do l=norb_dz+1,lr-1
c                lpcoe(l)=lp_coe(l,mpl)
c             enddo
c            lpcoe(lr)=kcoe
c             wl=voint(lr0,lr)
c             do l=lr0,lr
c               list=list3(lr0,lr,l)
c               kcoe=lpcoe(l)
c               call neoc(kcoe,nocc,tcoe)
c               wl=wl+nocc*(vint_ci(list+1)+tcoe*vint_ci(list))
c             enddo
c             wl=wl*vlop0
32                wl0=vlop0
c     list=list3(lr0,lr,lrs)   ! lrg
c             wl=vlop0*vint_ci(list)
c         write(nf2,*) 'ar-cw-ar'
         call trans_ijkl_intpos(lr,lrs,lr0,lrs,nxo)
             list0=nxo
             goto 500
c-----------------------------------------------------------
11       wl0=vlop0-vlop1
c          list=list4(lr0,lrg,lrs,lr)
c             wl=vlop0*(vint_ci(list)-2*vint_ci(list+1))
c     :             -vlop1*vint_ci(list)

         call trans_ijkl_intpos(lr,lrg,lrs,lr0,nxo)
             list0=nxo
         wl1=-2*vlop0
         call trans_ijkl_intpos(lr,lrs,lrg,lr0,nxo)
             list1=nxo
             goto 500
c-----------------------------------------------------------
12       wl0=vlop0-vlop1
c           list=list4(lr0,lrg,lrs,lr)
c             wl=vlop0*(vint_ci(list+2)-2.0d0*vint_ci(list+1))
c     :           -vlop1*vint_ci(list+2)

         call trans_ijkl_intpos(lr,lr0,lrg,lrs,nxo)
             list0=nxo
              wl1=-2*vlop0
         call trans_ijkl_intpos(lr,lrs,lrg,lr0,nxo)
             list1=nxo
             goto 500
c-----------------------------------------------------------
13       wl0=vlop0-vlop1
c     list=list4(lr0,lrg,lrs,lr)
c             wl=vlop0*(vint_ci(list+2)+vint_ci(list))
c     :             -vlop1*(vint_ci(list+2)-vint_ci(list))

         call trans_ijkl_intpos(lr,lr0,lrg,lrs,nxo)
             list0=nxo

              wl1=vlop0+vlop1
         call trans_ijkl_intpos(lr,lrg,lrs,lr0,nxo)
             list1=nxo

             goto 500
c-----------------------------------------------------------
21       wl0=vlop0+vlop1
c          list=list3(lr0,lrg,lr)
c             wl=(vlop0+vlop1)*vint_ci(list)

         call trans_ijkl_intpos(lrg,lr,lr0,lr,nxo)
             list0=nxo

             goto 500
c-----------------------------------------------------------
22       wl0=vlop0-vlop1
c    list=list3(lr0,lrg,lr)
c             wl=vlop0*(vint_ci(list)-2.0d0*vint_ci(list+1))
c     :           -vlop1*vint_ci(list)

         call trans_ijkl_intpos(lrg,lr,lr0,lr,nxo)
             list0=nxo

              wl1=-2.0d0*vlop0
         call trans_ijkl_intpos(lrg,lr0,lr,lr,nxo)
             list1=nxo

             goto 500
c-----------------------------------------------------------
c  drr-drr
c  w0lp=2.0d0*w0lp but not 1.0d0*w0lp is based on that the non-diagonal
c  just uses the non-triangle <ci|h|cj> which designates that i > j.
c
51        wl0=vlop0
c          wl=vlop0*voint(lr,lr0)*0.5d0

         call trans_ijkl_intpos(lr,lr0,lr,lr0,nxo)
             list0=nxo
              goto 500
c-----------------------------------------------------------
52       wl0=(vlop0-vlop1)*2.0d0
c         wl=(vlop0-vlop1)*voint(lr,lr0)

         call trans_ijkl_intpos(lr,lr0,lr,lr0,nxo)
             list0=nxo

              goto 500
c-----------------------------------------------------------
41       wl0=vlop0+vlop1
c        list=list3(lrs,lr,lr0)      ! lrg
c              wl=(vlop0+vlop1)*vint_ci(list)

         call trans_ijkl_intpos(lr,lr0,lrs,lr0,nxo)
             list0=nxo

              goto 500
c-----------------------------------------------------------
42       wl0=vlop0-vlop1
c       list=list3(lrs,lr,lr0)      ! lrg
c              wl=(vlop0-vlop1)*vint_ci(list)

         call trans_ijkl_intpos(lr,lr0,lrs,lr0,nxo)
             list0=nxo

              goto 500
c-----------------------------------------------------------
43       wl0=vlop0-vlop1
c      list=list3(lrs,lr,lr0)      ! lrg
c              wl=vlop0*(vint_ci(list)-2.0d0*vint_ci(list+1))
c     :           -vlop1*vint_ci(list)

         call trans_ijkl_intpos(lr,lr0,lrs,lr0,nxo)
             list0=nxo

              wl1=-2*vlop0
         call trans_ijkl_intpos(lr,lrs,lr0,lr0,nxo)
             list1=nxo

c-----------------------------------------------------------
500   return
      end


      subroutine dbl_sd_act_comp_g(lin,lra)
!sd(6-1) a&r(02)-
!sd(6-2) (22)a&(13)-
!sd(6-3) a&r(13)c'(22)-
!sd(6-4) a&r(23)c'(12)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      common/onepl/line,jph,jpel,jper,lrg,lrs,jwl,jwr,w0,w1
      jmlr=mul_tab(jml,jmr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jmlr) cycle
        w0sd1 =w0_sd(1)
        w0sd2 =w0_sd(2)
          w0sd3 =w0_sd(3)
        w0sd4 =w0_sd(4)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1) then
          w0sd1 =-w0sd1
          w0sd2 =-w0sd2
          w0sd3 =-w0sd3
          w0sd4 =-w0sd4
        endif
        if(jml.eq.1.and.lmi.eq.jmr) then
!sd(6-1) a&r(02)-
          iwdl=just(lri,lri)
          iwdr=jud(lri)
          vlop0=w0*w0sd1
          vlop1=w1*w0sd1
          call comp_loop_g(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl0,list0,
     :                     wl1,list1)
c          call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
c           call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl2,jper,list2,list3)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl0,jper,list0)
      if(list1.ne.0)
     :        call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl1,jper,list1)
        endif
!sd(6-2) (22)a&(13)-
        vlop0=w0*w0sd2
        vlop1=w1*w0sd2
        call comp_loop_g(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl0,list0,
     :                     wl1,list1)
        do lrk=norb_frz+1,lri-1
          lmk=lsm_inn(lrk)
          if(lmk.ne.jmr) cycle
          iwdl=just(lrk,lri)
          iwdr=jud(lrk)
c          call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
c           call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl2,jper,list2,list3)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl0,jper,list0)
      if(list1.ne.0)
     :        call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl1,jper,list1)
        enddo
!sd(6-4) a&r(23)-c'(12)-
        vlop0=-w0*w0sd4
        vlop1=-w1*w0sd4
        call comp_loop_g(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl0,list0,
     :                     wl1,list1)
        do lrk=lri+1,norb_dz
          lmk=lsm_inn(lrk)
          if(lmk.ne.jmr) cycle
          iwdl=just(lri,lrk)
          iwdr=jud(lrk)
c          call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
c           call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl2,jper,list2,list3)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl0,jper,list0)
      if(list1.ne.0)
     :        call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl1,jper,list1)
        enddo
!sd(6-3) a&r(13)c'(22)-
        if(jb_sys.gt.0) then
          vlop0=-w0*w0sd3
          vlop1=-w1*w0sd3
          call comp_loop_g(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl0,list0,
     :                     wl1,list1)
          do lrk=lri+1,norb_dz
            lmk=lsm_inn(lrk)
            if(lmk.ne.jmr) cycle
            iwdl=just(lrk,lri)
            iwdr=jud(lrk)
c            call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
c           call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl2,jper,list2,list3)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl0,jper,list0)
      if(list1.ne.0)
     :        call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl1,jper,list1)
          enddo
        endif
      enddo
      return
      end

      subroutine dbl_sdd_act_comp_g(lin,lra)
!sd1(8-1)    ar(01)-
!sd1(8-2)    (11)ar(23)-
!sd1(8-3)    ar(13)-c'(21)-
!sd1(8-4)    ar(23)-c'(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      common/onepl/line,jph,jpel,jper,lrg,lrs,jwl,jwr,w0,w1
      jmlr=mul_tab(jml,jmr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jmlr) cycle
        w0sd1 =w0_sd1(1)
        w0sd2 =w0_sd1(2)
          w0sd3 =w0_sd1(3)
        w0sd4 =w0_sd1(4)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1) then
          w0sd1 =-w0sd1
          w0sd2 =-w0sd2
            w0sd3 =-w0sd3
          w0sd4 =-w0sd4
        endif
        if(jml.eq.1.and.lmi.eq.jmr) then
!sd1(8-1)    ar(01)-
          iwdl=just(lri,lri)
          iwdr=jud(lri)
          vlop0=w0*w0sd1
          vlop1=w1*w0sd1
          call comp_loop_g(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl0,list0,
     :                     wl1,list1)
c          call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
c           call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl2,jper,list2,list3)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl0,jper,list0)
      if(list1.ne.0)
     :        call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl1,jper,list1)
        endif
!sd1(8-2)    (11)ar(23)-
        vlop0=w0*w0sd2
        vlop1=w1*w0sd2
        call comp_loop_g(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl0,list0,
     :                     wl1,list1)
        do lrk=norb_frz+1,lri-1
          lmk=lsm_inn(lrk)
          if(lmk.ne.jmr) cycle
          iwdl=just(lri,lrk)
          iwdr=jud(lrk)
c          call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
c           call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl2,jper,list2,list3)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl0,jper,list0)
      if(list1.ne.0)
     :        call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl1,jper,list1)
        enddo
!sd1(8-4)    ar(23)-c'(11)-
        vlop0=-w0*w0sd4
        vlop1=-w1*w0sd4
        call comp_loop_g(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl0,list0,
     :                     wl1,list1)
        do lrk=lri+1,norb_dz
          lmk=lsm_inn(lrk)
          if(lmk.ne.jmr) cycle
          iwdl=just(lri,lrk)
          iwdr=jud(lrk)
c          call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
c           call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl2,jper,list2,list3)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl0,jper,list0)
      if(list1.ne.0)
     :        call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl1,jper,list1)
        enddo
!sd1(8-3)    ar(13)-c'(21)-
        if(jb_sys.gt.0) then
          vlop0=-w0*w0sd3
          vlop1=-w1*w0sd3
          call comp_loop_g(lin,lri,lrg,lrs,lra,vlop0,vlop1,wl0,list0,
     :                     wl1,list1)
          do lrk=lri+1,norb_dz
            lmk=lsm_inn(lrk)
            if(lmk.ne.jmr) cycle
            iwdl=just(lrk,lri)
            iwdr=jud(lrk)
c            call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
c           call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl2,jper,list2,list3)
              call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl0,jper,list0)
      if(list1.ne.0)
     :        call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl1,jper,list1)
          enddo
        endif
      enddo
      return
      end


      subroutine ext_space_loop_g()
#include "drt_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *      m_jc,m_jd, isegsta,isegupwei,isegdownwei

      ism=0
      do inx=18,25
          ism=ism+1
            if(nu_ae(inx).eq.0) cycle
            isegsta=iseg_sta(inx)
            isegupwei=iseg_upwei(inx)
            isegdownwei=iseg_downwei(inx)
          call g_ss_ext_sequence_g(ism,4)
      enddo
      ism=0
      do inx=10,17
          ism=ism+1
            if(nu_ae(inx).eq.0) cycle
            isegsta=iseg_sta(inx)
            isegupwei=iseg_upwei(inx)
            isegdownwei=iseg_downwei(inx)
             call g_tt_ext_sequence_g(ism)
      enddo
      ism=0
      do inx=2,9
          ism=ism+1
            if(nu_ae(inx).eq.0) cycle
            isegsta=iseg_sta(inx)
            isegupwei=iseg_upwei(inx)
            isegdownwei=iseg_downwei(inx)
             call g_dd_ext_sequence_g(ism)
      enddo
      return
      end

      subroutine dbl_head_act_tail_g(lpcoe)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      dimension lpcoe(norb_dz+1:norb_inn)
      common/onepl/line,jph,jpel,jper,lrg,lrs,jwl,jwr,w0,w1
      lra=kk(jpel)-1
      jml=mod((jpadl-1),8)
      jmr=mod((jpad-1),8)
      itypadl=(jpadl-1)/8+2
      itypadr=(jpad-1)/8+2
      if(jml.eq.0) then
        jml=8
        itypadl=itypadl-1
      endif
      if(jmr.eq.0) then
        jmr=8
        itypadr=itypadr-1
      endif
      if(jpadl.eq.1) itypadl=1
      if(jpad.eq.1)  itypadr=1
      if(jpadl.eq.1) jml=ns_sm
      if(jpad.eq.1)  jmr=ns_sm
      jml=mul_tab(jml,ns_sm)
      jmr=mul_tab(jmr,ns_sm)
      jmlr=mul_tab(jml,jmr)
      lpok=map_jplr(itypadl,itypadr)
      if(lpok.eq.0) return

!          23   24 25 26 27 28 29 30 31 32 33 34 35
      goto(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300),
     :    line-22
!line=23:-a^l<-->ds(7),dds(9),dt(14),ddtt(16)
100   goto(10,10,10,10,10,10,107,10,109,10,10,10,
     :     10,114,10,116,10,10,10,10,10,10,10,10,10,10),lpok
!ds(7-2) ar(23)-bl(31)-br(32)-

!ds(7-1) ar(23)-drl(30)-

107   do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(jmr.ne.1) goto 106
        iwdr=just(lri,lri)
        do lrd=norb_frz+1,lri-1
          lmd=lsm_inn(lrd)
          if(lmd.ne.jml) cycle
          iwdl=jud(lrd)
          w0ds1 =w0_ds(1)
          ni=mod(norb_dz-lri+lri-lrd,2)
          if(ni.eq.0) w0ds1 =-w0ds1
          vlop0=w0*w0ds1
c          list=list3(lrd,lra,lri)
         call trans_ijkl_intpos(lra,lri,lrd,lri,nxo)
c        wl=vlop0*vint_ci(list)
          wl=vlop0
c         call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
        enddo
106     do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jmr) cycle
          do lrd=norb_frz+1,lri-1
            lmd=lsm_inn(lrd)
            if(lmd.ne.jml) cycle
c            list=list4(lrd,lri,lrj,lra)
            iwdl=jud(lrd)
            w0ds2 =w0_ds(2)
            w1ds2 =w1_ds(2)
            w0ds3 =w0_ds(3)
            w1ds3 =w1_ds(3)
            ni=mod(norb_dz-lrj+lri-lrd,2)
            if(ni.eq.0) then
             w0ds2 =-w0ds2
              w1ds2 =-w1ds2
             w0ds3 =-w0ds3
              w1ds3 =-w1ds3
            endif
!ds(7-3) ar(23)-bl(32)-br(31)-
            iwdr=just(lri,lrj)
           vlop0=w0*w0ds3
            vlop1=w1*w1ds3
c            wl=(vlop0-vlop1)*vint_ci(list)-
c    :             2*vlop0*vint_ci(list+1)            !1.1
            wl=vlop0-vlop1
         call trans_ijkl_intpos(lra,lri,lrj,lrd,nxo)

           call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            wl=-2*vlop0
         call trans_ijkl_intpos(lra,lrj,lri,lrd,nxo)

           call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            if(jb_sys.gt.0) then
!ds(7-2) ar(23)-bl(31)-br(32)-         the symmetry problem
              iwdr=just(lrj,lri)
             vlop0=w0*w0ds2
              vlop1=w1*w1ds2
c              wl=(vlop0-vlop1)*vint_ci(list)-
c     :               2*vlop0*vint_ci(list+1)            !1.1
            wl=vlop0-vlop1
         call trans_ijkl_intpos(lra,lri,lrj,lrd,nxo)

           call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            wl=-2*vlop0
         call trans_ijkl_intpos(lra,lrj,lri,lrd,nxo)

           call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

            endif
          enddo
        enddo
      enddo
      goto 10
109   do lri=norb_frz+1,norb_dz
!d1s(9-1) ar(13)-drl(30)-
        lmi=lsm_inn(lri)
        do lrd=norb_frz+1,lri-1
          lmd=lsm_inn(lrd)
          if(lmd.eq.jml.and.jmr.eq.1) then
            iwdr=just(lri,lri)
            iwdl=jud(lrd)
            w0ds1=w0_d1s(1)
            ni=mod(norb_dz-lri+lri-lrd,2)
            if(ni.eq.0) w0ds1 =-w0ds1
            vlop0=w0*w0ds1
c            list=list3(lrd,lra,lri)
         call trans_ijkl_intpos(lra,lri,lrd,lri,nxo)

c          wl=vlop0*vint_ci(list)          !   3.2
            wl=vlop0
            call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
          endif
!d1s(9-4) drl(12)-br(31)-
          if(jml.eq.lmd.and.jmr.eq.mul_tab(lmd,lmi)) then
            iwdr=just(lrd,lri)
            iwdl=jud(lrd)
            w1ds=w1_d1s(4)
           if(mod(norb_dz-lri,2).eq.1) w1ds=-w1ds
           vlop1=w1*w1ds
c          list=list3(lri,lra,lrd)
         call trans_ijkl_intpos(lra,lrd,lri,lrd,nxo)

c          wl=-vlop1*vint_ci(list)
            wl=-vlop1
            call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
          endif
        enddo
      enddo
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
!d1s(9-3) ar(13)-bl(32)-br(31)-
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jmr) cycle
          do lrd=norb_frz+1,lri-1
            iwdr=just(lri,lrj)
            lmd=lsm_inn(lrd)
            if(lmd.ne.jml) cycle
            iwdl=jud(lrd)
            w0ds3 =w0_d1s(3)
            w1ds3 =w1_d1s(3)
            ni=mod(norb_dz-lrj+lri-lrd,2)
            if(ni.eq.0) w0ds3 =-w0ds3
            if(ni.eq.0) w1ds3 =-w1ds3
           vlop0=w0*w0ds3
            vlop1=w1*w1ds3
c          list=list4(lrd,lri,lrj,lra)
c            wl=(vlop0-vlop1)*vint_ci(list)-
c    :               2*vlop0*vint_ci(list+1)            !1.1
            wl=vlop0-vlop1
         call trans_ijkl_intpos(lra,lri,lrj,lrd,nxo)
           call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            wl=-2*vlop0
         call trans_ijkl_intpos(lra,lrj,lri,lrd,nxo)
           call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

            if(jb_sys.gt.0) then
!d1s(9-2)   ar(13)-bl(31)-br(32)-   the symmetry problem
              iwdr=just(lrj,lri)
              w0ds3 =w0_d1s(2)
              w1ds3 =w1_d1s(2)
              ni=mod(norb_dz-lrj+lri-lrd,2)
              if(ni.eq.0) w0ds3 =-w0ds3
              if(ni.eq.0) w1ds3 =-w1ds3
             vlop0=w0*w0ds3
              vlop1=w1*w1ds3
c           list=list4(lrd,lri,lrj,lra)
c              wl=(vlop0-vlop1)*vint_ci(list)-
c    :                 2*vlop0*vint_ci(list+1)            !1.1
            wl=vlop0-vlop1
         call trans_ijkl_intpos(lra,lri,lrj,lrd,nxo)
           call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            wl=-2*vlop0
         call trans_ijkl_intpos(lra,lrj,lri,lrd,nxo)
           call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

            endif
          enddo
        enddo
      enddo
      return

!dt(14) ar(23)-bl(32)-br(32)-
114   do lri=norb_frz+1,norb_dz-1
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
         lmij=mul_tab(lmi,lmj)
         if(lmij.ne.jmr) cycle
          iwdr=just(lri,lrj)        !
         do lrd=norb_frz+1,lri-1
            lmd=lsm_inn(lrd)
           if(lmd.ne.jml) cycle
            iwdl=jud(lrd)
           vlop0=w0*w0_dt
            vlop1=w1*w1_dt
            ni=mod(lri-lrd+norb_dz-lrj,2)
           if(ni.eq.0) then
             vlop0=-vlop0
              vlop1=-vlop1
            endif
c          list=list4(lrd,lri,lrj,lra)
c            wl=(vlop0-vlop1)*vint_ci(list)-2*vlop0*vint_ci(list+1) !1.1
            wl=vlop0-vlop1
         call trans_ijkl_intpos(lra,lri,lrj,lrd,nxo)
           call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            wl=-2*vlop0
         call trans_ijkl_intpos(lra,lrj,lri,lrd,nxo)
           call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

          enddo
        enddo
      enddo
      return
!d1t1(16)  ar(13)-bl(31)-br(31)-
116   do lri=norb_frz+1,norb_dz-1
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
         lmij=mul_tab(lmi,lmj)
         if(lmij.ne.jmr) cycle
          iwdr=just(lri,lrj)        !
         do lrd=norb_frz+1,lri-1
            lmd=lsm_inn(lrd)
            lmd=mul_tab(lmd,1)
           if(lmd.ne.jml) cycle
            iwdl=jud(lrd)
           vlop0=w0*w0_d1t1
            vlop1=w1*w1_d1t1
            ni=mod(lri-lrd+norb_dz-lrj,2)
           if(ni.eq.0) then
             vlop0=-vlop0
              vlop1=-vlop1
            endif
c          list=list4(lrd,lri,lrj,lra)
c            wl=(vlop0-vlop1)*vint_ci(list)-2*vlop0*vint_ci(list+1) !1.1
            wl=vlop0-vlop1
         call trans_ijkl_intpos(lra,lri,lrj,lrd,nxo)
           call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            wl=-2*vlop0
         call trans_ijkl_intpos(lra,lrj,lri,lrd,nxo)
           call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

          enddo
        enddo
      enddo
      return
!line=24:-a^r<-->sd(6),sdd(8),td(13),ttdd(15),dv(23),ddv(24)
200   goto(10,10,10,10,10,206,10,208,10,10,10,10,
     :     213,10,215,10,10,10,10,10,10,10,223,224,10,10),lpok


!sd(6-1) a&r(02)-
!sd(6-2) c(22)a&(13)-
!sd(6-3) a&r(13)c'(22)-
!sd(6-4) a&r(23)c'(12)-
!sd(6-5) a&r(23)b&r(13)b^r(32)
!sd(6-6) a&r(13)b&r(23)b^r(32)
!sd(6-7) a&r(13)b&l(32)b^l(23)
!sd(6-8) a&r(23)b&l(32)b^l(13)
!sd(6-9) d&r&r(03)b^r(32)
!sd(6-10) d&r&l(12)b^l(23)
!sd(6-11) d&r&l(22)b^l(13)
!sd(6-12) d&r&l(33)b^l(02)
!sd(6-13) (22)d&r&l(33)b^l(13)
!sd(6-14) d&r&l(33)c"(22)b^l(13)
!sd(6-15) d&r&l(33)b^l(13)c'(22)
!sd(6-16) d&r&l(33)b^l(23)c'(12)
!sd(6-1) a&r(02)-

206   call sd_head_dbl_tail_act_g(lra,lpcoe)
      return

208   call sdd_head_dbl_tail_act_g(lra,lpcoe)
      return
!td(13-1) (22)a&(23)
!td(13-1) a&(23)c'(22)
!td(13-5) (22)d&&l(33)b^l(23)
213   do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jmlr) cycle
        w0td1=w0_td(1)
        w0td4=w0_td(4)
        w0td5=w0_td(5)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1)w0td1=-w0td1
        if(ni.eq.1)w0td4=-w0td4
        if(ni.eq.1)w0td5=-w0td5

!td(13-1) a&(23)c'(22)
        do lrd=lri+1,norb_dz
          lmd=lsm_inn(lrd)
          if(lmd.ne.jmr) cycle
         iwdl=just(lri,lrd)      !
          iwdr=jud(lrd)
          vlop0=-w0*w0td1
c          list=list3(lri,lra,lri)
c          wl=voint(lri,lra)+vint_ci(list)             !310,act_coe,610,
           wl=vlop0
          call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,lri,lra)
         call trans_ijkl_intpos(lra,lri,lri,lri,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

c          list=list3(lri,lra,lrd)
c          wl=wl+vint_ci(list+1)                          !310 c'(22) co
         call trans_ijkl_intpos(lra,lri,lrd,lrd,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

         do lr=lri+1,norb_dz
            if(lr.eq.lrd) cycle
c          list =list3(lri,lra,lr)
c            wl=wl+2*vint_ci(list+1)-vint_ci(list)       !310:neoc=2,coe
            wl=2.0d0*vlop0
         call trans_ijkl_intpos(lra,lri,lr,lr,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            wl=-vlop0
         call trans_ijkl_intpos(lra,lr,lri,lr,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

          enddo
          do lrk=norb_dz+1,lra
c            list=list3(lri,lra,lrk)
            kcoe=lpcoe(lrk)
            call neoc(kcoe,nocc,tcoe)
c            wl=wl+nocc*(vint_ci(list+1)+tcoe*vint_ci(list))
            wl=vlop0*nocc
         call trans_ijkl_intpos(lra,lri,lrk,lrk,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            wl=vlop0*tcoe*nocc
         call trans_ijkl_intpos(lra,lrk,lri,lrk,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

          enddo
c          wl=wl*vlop0
!td(13-5) d&rl(33)b^l(23)c'(22)
          vlop0=-w0*w0td5
          do lrk=1,lri-1
c            list=list3(lri,lra,lrk)
c            wl=wl-vlop0*(2*vint_ci(list+1)-vint_ci(list))
            wl=-vlop0*2
         call trans_ijkl_intpos(lra,lri,lrk,lrk,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            wl=vlop0
         call trans_ijkl_intpos(lra,lrk,lri,lrk,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

          enddo
c          call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
        enddo
!-------------------------------------------------------------------
        do lrd=norb_frz+1,lri-1
          lmd=lsm_inn(lrd)
          if(lmd.ne.jmr) cycle
         iwdl=just(lrd,lri)      !
          iwdr=jud(lrd)
!td(13-1) (22)a&(23)
          vlop0=w0*w0td1
c          list=list3(lri,lra,lri)
c          wl=vlop0*(voint(lri,lra)+vint_ci(list))     !310,act_coe,610,
           wl=vlop0
          call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,lri,lra)
         call trans_ijkl_intpos(lra,lri,lri,lri,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

          do lr=lri+1,norb_dz
c            list =list3(lri,lra,lr)
c            wl=wl+vlop0*(2*vint_ci(list+1)-vint_ci(list)) !  310:neoc=2
            wl=vlop0*2
         call trans_ijkl_intpos(lra,lri,lr,lr,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            wl=-vlop0
         call trans_ijkl_intpos(lra,lr,lri,lr,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

          enddo
          do lrk=norb_dz+1,lra
c            list=list3(lri,lra,lrk)
            kcoe=lpcoe(lrk)
            call neoc(kcoe,nocc,tcoe)
c           wl=wl+vlop0*nocc*(vint_ci(list+1)+tcoe*vint_ci(list))
            wl=vlop0*nocc
         call trans_ijkl_intpos(lra,lri,lrk,lrk,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            wl=vlop0*tcoe*nocc
         call trans_ijkl_intpos(lra,lrk,lri,lrk,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

          enddo
c            wl=wl*vlop0
!td(13-4) d&r&l(22)b^l(23)
          vlop0=w0*w0td4
          vlop1=w1*w0td4
c          list=list3(lri,lra,lrd)
c          wl=wl+(vlop0-vlop1)*vint_ci(list)-2*vlop0*vint_ci(list+1)
            wl=-2*vlop0
         call trans_ijkl_intpos(lra,lri,lrd,lrd,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            wl=vlop0-vlop1
         call trans_ijkl_intpos(lra,lrd,lri,lrd,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)


!td(13-5) d&rl(33)c"(22)b^l(23)
          vlop0=w0*w0td5
          do lrk=1,lri-1
            if(lrk.eq.lrd) cycle
c          list=list3(lri,lra,lrk)
c            wl=wl+vlop0*(vint_ci(list)-2*vint_ci(list+1))      !4.3
            wl=-2*vlop0
         call trans_ijkl_intpos(lra,lri,lrk,lrk,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            wl=vlop0
         call trans_ijkl_intpos(lra,lrk,lri,lrk,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

          enddo
c          call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
        enddo
      enddo

      do lri=norb_frz+1,norb_dz-1
        lmi=lsm_inn(lri)
      do lrj=lri+1,norb_dz
        lmj=lsm_inn(lrj)
        lmij=mul_tab(lmi,lmj)
        if(lmij.ne.jml) cycle
        iwdl=just(lri,lrj)     !

!td(13-2) a&(23)b&r(23)b^r(32)
        do lrd=lrj+1,norb_dz
          lmd=lsm_inn(lrd)
          if(lmd.ne.jmr) cycle
          w0td2=w0_td(2)
          w1td2=w1_td(2)
          ni=mod(lrj-lri+norb_dz-lrd,2)
        if(ni.eq.0) w0td2=-w0td2
        if(ni.eq.0) w1td2=-w1td2

          iwdr=jud(lrd)
          vlop0=w0*w0td2
          vlop1=w1*w1td2
c          list=list4(lri,lrj,lrd,lra)
c          wl=vlop0*(vint_ci(list+2)+vint_ci(list))  !1.3
c     :       -vlop1*(vint_ci(list+2)-vint_ci(list))
         wl=vlop0-vlop1
         call trans_ijkl_intpos(lra,lri,lrj,lrd,nxo)
        call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
         wl=vlop0+vlop1
         call trans_ijkl_intpos(lra,lrj,lrd,lri,nxo)
        call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

        enddo
!td(13-3) a&(23)b&l(32)b^l(23)
        do lrd=lri+1,lrj-1
          lmd=lsm_inn(lrd)
          if(lmd.ne.jmr) cycle
          iwdr=jud(lrd)
          w0td3=w0_td(3)
          w1td3=w1_td(3)
          ni=mod(lrd-lri+norb_dz-lrj,2)
          if(ni.eq.0)   w0td3=-w0td3
          if(ni.eq.0)   w1td3=-w1td3
          vlop0=w0*w0td3                !d6-8
          vlop1=w1*w1td3
c          list=list4(lri,lrd,lrj,lra)
c       wl=vlop0*(vint_ci(list+2)-2*vint_ci(list+1))      !1.2
c     :       -vlop1*vint_ci(list+2)

          wl=vlop0-vlop1
         call trans_ijkl_intpos(lra,lri,lrd,lrj,nxo)
        call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
          wl=-vlop0*2
         call trans_ijkl_intpos(lra,lrj,lrd,lri,nxo)
        call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

        enddo
      enddo
      enddo
      goto 10
215   call ttdd_head_dbl_tail_act_g(lra,lpcoe)
      return

!dv(23-1) ar(23)-
!dv(23-2) drl(33)-bl(23)-
223   iwdr=0
      do lrd=norb_frz+1,norb_dz
        imd=lsm_inn(lrd)
        if(imd.ne.jml) cycle
        iwdl=jud(lrd)
        w0dv1=w0_dv(1)
        ni=mod(norb_dz-lrd,2)
        if(ni.eq.1) w0dv1=-w0dv1
        vlop0=w0*w0dv1                !d23-1
        vlop1=w1*w0dv1
c**********************************************************
              lr0=lrd
              lr=kk(jpel) -1
c              list=list3(lr0,lr,lr0)
c            wl=vlop0*(voint(lr0,lr)+vint_ci(list))       !310+710
              wl=vlop0

          call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,lr0,lr)
         call trans_ijkl_intpos(lr,lr0,lr0,lr0,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
          do l=lr0+1,norb_dz
c            list=list3(lr0,lr,l)
           nocc=2
            tcoe=-0.5d0
c           wl=wl+nocc*vlop0*(vint_ci(list+1)+tcoe*vint_ci(list))  !dbl_
            wl=vlop0*nocc
         call trans_ijkl_intpos(lr,lr0,l,l,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            wl=vlop0*tcoe*nocc
         call trans_ijkl_intpos(lr,l,lr0,l,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

         enddo
          do l=norb_dz+1,lr
c            list=list3(lr0,lr,l)
           kcoe=lpcoe(l)
          call neoc(kcoe,nocc,tcoe)
c         wl=wl+nocc*vlop0*(vint_ci(list+1)+tcoe*vint_ci(list))   !act_c
            wl=vlop0*nocc
         call trans_ijkl_intpos(lr,lr0,l,l,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            wl=vlop0*tcoe*nocc
         call trans_ijkl_intpos(lr,l,lr0,l,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

         enddo
c          wl_430=0.d0
          w0dv2=w0_dv(2)
          ni=mod(norb_dz-lrd,2)
          if(ni.eq.1) w0dv2=-w0dv2
         do lrk=1,lrd-1
c            list=list3(lr0,lr,lrk)
            vlop0=w0*w0dv2
c            wl_430=wl_430+vlop0*(vint_ci(list)-2*vint_ci(list+1))
            wl=vlop0
         call trans_ijkl_intpos(lr,lrk,lr0,lrk,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            wl=-vlop0*2
         call trans_ijkl_intpos(lr,lr0,lrk,lrk,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

          enddo
c          wl=wl+wl_430
c**********************************************************
c          call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
      enddo
      return
!d1v(24-1) ar(13)-
!d1v(24-2) drl(33)-bl(13)-
224   iwdr=0
      do lrd=norb_frz+1,norb_dz
        imd=lsm_inn(lrd)
        if(imd.ne.jml) cycle
        iwdl=jud(lrd)
        w0dv1=w0_d1v(1)
        ni=mod(norb_dz-lrd,2)
        if(ni.eq.1) w0dv1=-w0dv1
          vlop0=w0*w0dv1                !d24-1
        vlop1=w1*w0dv1
c**********************************************************
          lr0=lrd
          lr=kk(jpel) -1
c          list=list3(lr0,lr,lr0)
c           wl=vlop0*(voint(lr0,lr)+vint_ci(list))       !310+710
              wl=vlop0
          call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,lr0,lr)
         call trans_ijkl_intpos(lr,lr0,lr0,lr0,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

          do l=lr0+1,norb_dz
c            list=list3(lr0,lr,l)
             nocc=2
            tcoe=-0.5d0
c           wl=wl+nocc*vlop0*(vint_ci(list+1)+tcoe*vint_ci(list))  !dbl_
            wl=vlop0*nocc
         call trans_ijkl_intpos(lr,lr0,l,l,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            wl=vlop0*tcoe*nocc
         call trans_ijkl_intpos(lr,l,lr0,l,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

         enddo
          do l=norb_dz+1,lr
c            list=list3(lr0,lr,l)
           kcoe=lpcoe(l)
          call neoc(kcoe,nocc,tcoe)
c         wl=wl+nocc*vlop0*(vint_ci(list+1)+tcoe*vint_ci(list))   !act_c
            wl=vlop0*nocc
         call trans_ijkl_intpos(lr,lr0,l,l,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            wl=vlop0*tcoe*nocc
         call trans_ijkl_intpos(lr,l,lr0,l,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

         enddo
c          wl_430=0.d0
          w0dv2=w0_d1v(2)
          ni=mod(norb_dz-lrd,2)
          if(ni.eq.1) w0dv2=-w0dv2
         do lrk=1,lrd-1
c            list=list3(lr0,lr,lrk)
            vlop0=w0*w0dv2
c            wl_430=wl_430+vlop0*(vint_ci(list)-2*vint_ci(list+1))
            wl=vlop0
         call trans_ijkl_intpos(lr,lrk,lr0,lrk,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
            wl=-vlop0*2
         call trans_ijkl_intpos(lr,lr0,lrk,lrk,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

          enddo
c          wl=wl+wl_430
c**********************************************************
c          call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
      enddo
      return

!line=25:-d^r^r<-->sv(10),tv(17),ttv(18)
300   goto(10,10,10,10,10,10,10,10,10,310,10,10,
     :     10,10,10,10,317,318,10,10,10,10,10,10,10,10),lpok

310   call sv_head_dbl_tail_act_g(lra)
      return
!tv(17) ar(23)-br(23)-
317   iwdr=0
      do lri=norb_frz+1,norb_dz
        imi=lsm_inn(lri)
      do lrj=lri,norb_dz
        imj=lsm_inn(lrj)
        imij=mul_tab(imi,imj)
        if(imij.ne.jml) cycle
        iwdl=just(lri,lrj)     !
        vlop1=w1*w1_tv             !d17 vlop0=0
c       list=list3(lri,lrj,lra)
c        wl=vlop1*vint_ci(list)        !2.1                  !!!!!
        wl=vlop1
         call trans_ijkl_intpos(lrj,lra,lri,lra,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
      enddo
      enddo
      return

!t1v(18) ar(13)-br(13)-
318   iwdr=0
      do lri=norb_frz+1,norb_dz
        imi=lsm_inn(lri)
      do lrj=lri,norb_dz
        imj=lsm_inn(lrj)
        imij=mul_tab(imi,imj)
        if(imij.ne.jml) cycle
        iwdl=just(lri,lrj)     !
        vlop1=w1*w1_t1v             !d18 vlop0=0
c       list=list3(lri,lrj,lra)
c        wl=vlop1*vint_ci(list)        !2.1                  !!!!!
        wl=vlop1
         call trans_ijkl_intpos(lrj,lra,lri,lra,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

      enddo
      enddo
      return

!line=26:-d^r^l<-->ss(1),st(2),ts(3),stt(4),tts(5),tt(11),tttt(12),dd(19
400   goto(401,402,403,404,405,10,10,10,10,10,411,412,
     :     10,10,10,10,10,10,419,420,421,422,10,10,425,10),lpok

401   call ss_head_dbl_tail_act_g(lra)
      return

402   call st_head_dbl_tail_act_g(lra)
      return
c=======================================================================
403   call ts_head_dbl_tail_act_g(lra)
      return

404   call stt_head_dbl_tail_act_g(lra)
      return

405   call tts_head_dbl_tail_act_g(lra)
      return

411   call tt_head_dbl_tail_act_g(lra)
      return

412   call tttt_head_dbl_tail_act_g(lra)
      return

419   call dd_head_dbl_tail_act_g(lra)
      return

420   call dddd_head_dbl_tail_act_g(lra)
      return

421   call dd1_head_dbl_tail_act_g(lra)
      return

422   call d1d_head_dbl_tail_act_g(lra)
      return

!vv(25) drl(33)-
425   if(jwl.eq.jwr) return
      vlop0=w0*w0_vv             !d25
      wl=0.d0
      iwdl=0
      iwdr=0
      do lri=1,norb_dz
c       wl=wl+vlop0*voint(lri,lra)
        wl=vlop0
          call prodab_1(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,lri,lra)
      enddo
c      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
      return

!line=27:-b^r-a^r<-->sv(10),tv(17),ttv(18)
500   goto(10,10,10,10,10,10,10,10,10,510,10,10,
     :   10,10,10,10,517,518,10,10,10,10,10,10,10,10),lpok

510   call sv_head_dbl_tail_act_g(lra)
      return
!tv(17) ar(23)-br(23)-
517   iwdr=0
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
      do lrj=lri+1,norb_dz
        lmj=lsm_inn(lrj)
        lmij=mul_tab(lmi,lmj)
        if(lmij.ne.jml) cycle
        w1tv=w1_tv
        if(mod(lrj-lri,2).eq.0) w1tv=-w1tv
        iwdl=just(lri,lrj)    !
        vlop1=w1*w1tv             !d17
c        list=list4(lri,lrj,lrs,lra)
c        wl=vlop1*(vint_ci(list)-vint_ci(list+2)) !1.3 vlop0=0      !!!!
c        call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
          wl=vlop1
         call trans_ijkl_intpos(lra,lrj,lrs,lri,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
          wl=-vlop1
         call trans_ijkl_intpos(lra,lri,lrj,lrs,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

      enddo
      enddo
      return
!t1v(18) ar(13)-br(13)-
518   iwdr=0
      do lri=norb_frz+1,norb_dz-1
        lmi=lsm_inn(lri)
      do lrj=lri+1,norb_dz
        lmj=lsm_inn(lrj)
        lmij=mul_tab(lmi,lmj)
        if(lmij.ne.jml) cycle
        w1tv=w1_t1v
        if(mod(lrj-lri,2).eq.0) w1tv=-w1tv
        iwdl=just(lri,lrj)
        vlop1=w1*w1tv             !d17
c        list=list4(lri,lrj,lrs,lra)
c        wl=vlop1*(vint_ci(list)-vint_ci(list+2)) !1.3 vlop0=0      !!!!
c        call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
          wl=vlop1
         call trans_ijkl_intpos(lra,lrj,lrs,lri,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
          wl=-vlop1
         call trans_ijkl_intpos(lra,lri,lrj,lrs,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
      enddo
      enddo
      return


!line=28:-b^l-a^r<-->ss(1),st(2),ts(3),stt(4),tts(5),tt(11),tttt(12),dd(
600   goto(601,602,603,604,605,10,10,10,10,10,611,612,
     :     10,10,10,10,10,10,619,620,621,622,10,10,625,10),lpok

601   call ss_head_dbl_tail_act_g(lra)
      return

602   call st_head_dbl_tail_act_g(lra)
      return

603   call ts_head_dbl_tail_act_g(lra)
      return

604   call stt_head_dbl_tail_act_g(lra)
      return

605   call tts_head_dbl_tail_act_g(lra)
      return

611   call tt_head_dbl_tail_act_g(lra)
      return

612   call tttt_head_dbl_tail_act_g(lra)
      return

619   call dd_head_dbl_tail_act_g(lra)
      return

620   call dddd_head_dbl_tail_act_g(lra)
      return

621   call dd1_head_dbl_tail_act_g(lra)
      return

622   call d1d_head_dbl_tail_act_g(lra)
      return


!vv(25) drl(33)-
625   vlop0=w0*w0_vv             !d25
      wl=0.d0
      iwdl=0
      iwdr=0
        do lrk=1,norb_dz
c        list=list3(lrs,lra,lrk)
c        wl=wl+vlop0*(vint_ci(list)-2*vint_ci(list+1))   !4.3 vlop1=0
          wl=vlop0
         call trans_ijkl_intpos(lra,lrk,lrs,lrk,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
          wl=-vlop0*2
         call trans_ijkl_intpos(lra,lrs,lrk,lrk,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

        enddo
c      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
      return

!line=29:-b^r-a^l<-->ss(1),st(2),ts(3),stt(4),tts(5),tt(11),tttt(12),dd(
700   goto(701,702,703,704,705,10,10,10,10,10,711,712,
     :     10,10,10,10,10,10,719,720,721,722,10,10,725,10),lpok

701   call ss_head_dbl_tail_act_g(lra)
      return

702   call st_head_dbl_tail_act_g(lra)
      return

703   call ts_head_dbl_tail_act_g(lra)
      return

704   call stt_head_dbl_tail_act_g(lra)
      return

705   call tts_head_dbl_tail_act_g(lra)
      return

711   call tt_head_dbl_tail_act_g(lra)
      return

712   call tttt_head_dbl_tail_act_g(lra)
      return

719   call dd_head_dbl_tail_act_g(lra)
      return

720   call dddd_head_dbl_tail_act_g(lra)
      return

721   call dd1_head_dbl_tail_act_g(lra)
      return

722   call d1d_head_dbl_tail_act_g(lra)
      return

!vv(25) drl(33)-
725   if(jwl.ge.jwr) return
      vlop0=w0*w0_vv             !d25
      wl=0.d0
      iwdl=0
      iwdr=0
      do lrk=1,norb_dz
c          list=list3(lrs,lra,lrk)
c          wl=vlop0*vint_ci(list)      !4.2  vlop1=0         !!!!!
          wl=vlop0
         call trans_ijkl_intpos(lra,lrk,lrs,lrk,nxo)
         call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
      enddo
c      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
      return

!line=30:-b&r-d^r^r<-->sd(6),sdd(8),td(13),ttdd(15),dv(23),ddv(24)
800   goto(10,10,10,10,10,806,10,808,10,10,10,10,
     :     813,10,815,10,10,10,10,10,10,10,823,824,10,10),lpok

!                                    sd(6-3) a&r(13)c'(22)-
806   call dbl_sd_act_comp_g(3,lra)
      return
!sd(6-3) tmp for spin=0
!      return
808   call dbl_sdd_act_comp_g(3,lra)
      return
813   call dbl_td_act_comp_g(3,lra)
      return

815   call dbl_ttdd_act_comp_g(3,lra)
      return
!dv(23-1) ar(23)-
!dv(23-2) drl(33)-bl(23)-
823   iwdr=0
      do lrd=norb_frz+1,norb_dz
        imd=lsm_inn(lrd)
        if(imd.ne.jml) cycle
        iwdl=jud(lrd)
          w0dv1=w0_dv(1)
        ni=mod(norb_dz-lrd,2)
        if(ni.eq.1) w0dv1=-w0dv1
         vlop0=w0*w0dv1                   !d23-1
        vlop1=w1*w0dv1
c        list=list3(lrd,lrg,lra)
c        wl=(vlop0+vlop1)*vint_ci(list)        !2.1       !!!!!
c        call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
          wl= vlop0+vlop1
         call trans_ijkl_intpos(lrg,lra,lrd,lra,nxo)
         call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
      enddo
      return
!d1v(24-1)  ar(13)-
824   iwdr=0
      do lrd=norb_frz+1,norb_dz
        imd=lsm_inn(lrd)
        if(imd.ne.jml) cycle
        iwdl=jud(lrd)
        w0dv1=w0_d1v(1)
        ni=mod(norb_dz-lrd,2)
        if(ni.eq.1) w0dv1=-w0dv1
         vlop0=w0*w0dv1                   !d23-1
          vlop1=w1*w0dv1
c        list=list3(lrd,lrg,lra)
c        wl=(vlop0+vlop1)*vint_ci(list)        !2.1
c        call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
          wl= vlop0+vlop1
         call trans_ijkl_intpos(lrg,lra,lrd,lra,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
      enddo
      return

!line=31:-b&l-d^r^l<-->sd(6),sdd(8),td(13),ttdd(15),dv(23),ddv(24)
900   goto(10,10,10,10,10,906,10,908,10,10,10,10,
     :     913,10,915,10,10,10,10,10,10,10,923,924,10,10),lpok

906   call dbl_sd_act_comp_g(5,lra)
      return

908   call dbl_sdd_act_comp_g(5,lra)
      return
913   call dbl_td_act_comp_g(5,lra)
      return

915   call dbl_ttdd_act_comp_g(5,lra)
      return
!dv(23-1) ar(23)-
!dv(23-1) ar(23)-
!dv(23-2) drl(33)-bl(23)-
923   iwdr=0
      do lrd=norb_frz+1,norb_dz
        imd=lsm_inn(lrd)
        if(imd.ne.jml) cycle
        iwdl=jud(lrd)
        w0dv1=w0_dv(1)
        ni=mod(norb_dz-lrd,2)
        if(ni.eq.1) w0dv1=-w0dv1
        vlop0=w0*w0dv1                !d23-1
        vlop1=w1*w0dv1
c        list=list3(lrd,lrg,lra)
c        wl=vlop0*(vint_ci(list)-2*vint_ci(list+1))   !2.2          !!!!
c     :       -vlop1*(vint_ci(list))
c        call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
          wl= vlop0-vlop1
         call trans_ijkl_intpos(lrg,lra,lrd,lra,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
          wl= -vlop0*2
         call trans_ijkl_intpos(lrg,lrd,lra,lra,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
      enddo
      return
!d1v(24-1)  ar(13)-
924   iwdr=0
      do lrd=norb_frz+1,norb_dz
        imd=lsm_inn(lrd)
        if(imd.ne.jml) cycle
        iwdl=jud(lrd)
          w0dv1=w0_d1v(1)
        ni=mod(norb_dz-lrd,2)
        if(ni.eq.1) w0dv1=-w0dv1
        vlop0=w0*w0dv1                !d23-1
        vlop1=w1*w0dv1
c        list=list3(lrd,lrg,lra)
c        wl=vlop0*(vint_ci(list)-2*vint_ci(list+1))   !2.2          !!!!
c     :       -vlop1*(vint_ci(list))
c        call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
          wl= vlop0-vlop1
         call trans_ijkl_intpos(lrg,lra,lrd,lra,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
          wl= -vlop0*2
         call trans_ijkl_intpos(lrg,lrd,lra,lra,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

      enddo
      return


!line=32:-b&r-b^r-a^r<-->sd(6),sdd(8),td(13),ttdd(15),dv(23),ddv(24)
1000  goto(10,10,10,10,10,1006,10,1008,10,10,10,10,
     :     1013,10,1015,10,10,10,10,10,10,10,1023,1024,10,10),lpok
1006  call dbl_sd_act_comp_g(4,lra)
      return
1008  call dbl_sdd_act_comp_g(4,lra)
      return
1013  call dbl_td_act_comp_g(4,lra)
      return
1015  call dbl_ttdd_act_comp_g(4,lra)
      return
!dv(23-1) ar(23)-
!dv(23-2) drl(33)-bl(23)-
1023  iwdr=0
      do lrd=norb_frz+1,norb_dz
        imd=lsm_inn(lrd)
        if(imd.ne.jml) cycle
        iwdl=jud(lrd)
        w0dv1=w0_dv(1)
        ni=mod(norb_dz-lrd,2)
        if(ni.eq.1) w0dv1=-w0dv1
        vlop0=w0*w0dv1                !d23-1
        vlop1=w1*w0dv1
c        list=list4(lrd,lrg,lrs,lra)
c        wl=vlop0*(vint_ci(list+2)+vint_ci(list))  !1.3        !!!!!
c     :       -vlop1*(vint_ci(list+2)-vint_ci(list))
c        call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
          wl= vlop0+vlop1
         call trans_ijkl_intpos(lra,lrg,lrs,lrd,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
          wl= vlop0-vlop1
         call trans_ijkl_intpos(lra,lrd,lrg,lrs,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

      enddo
      return
!d1v(24-1)  ar(13)-
1024  iwdr=0
      do lrd=norb_frz+1,norb_dz
        imd=lsm_inn(lrd)
        if(imd.ne.jml) cycle
        iwdl=jud(lrd)
          w0dv1=w0_d1v(1)
        ni=mod(norb_dz-lrd,2)
        if(ni.eq.1) w0dv1=-w0dv1
        vlop0=w0*w0dv1                !d23-1
        vlop1=w1*w0dv1
c        list=list4(lrd,lrg,lrs,lra)
c        wl=vlop0*(vint_ci(list+2)+vint_ci(list))  !1.3        !!!!!
c     :       -vlop1*(vint_ci(list+2)-vint_ci(list))
c        call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
          wl= vlop0+vlop1
         call trans_ijkl_intpos(lra,lrg,lrs,lrd,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
          wl= vlop0-vlop1
         call trans_ijkl_intpos(lra,lrd,lrg,lrs,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
      enddo
      return

!line=33:-b&l-b^r-a^l<-->sd(6),sdd(8),td(13),ttdd(15),dv(23),ddv(24)
1100  goto(10,10,10,10,10,1106,10,1108,10,10,10,10,
     :     1113,10,1115,10,10,10,10,10,10,10,1123,1124,10,10),lpok

1106  call dbl_sd_act_comp_g(6,lra)
      return
1108  call dbl_sdd_act_comp_g(6,lra)
      return
1113  call dbl_td_act_comp_g(6,lra)
      return
1115  call dbl_ttdd_act_comp_g(6,lra)
      return

!dv(23-1) ar(23)-
!dv(23-2) drl(33)-bl(23)-
1123  iwdr=0
      do lrd=norb_frz+1,norb_dz
        imd=lsm_inn(lrd)
        if(imd.ne.jml) cycle
        iwdl=jud(lrd)
          w0dv1=w0_dv(1)
        ni=mod(norb_dz-lrd,2)
        if(ni.eq.1) w0dv1=-w0dv1
        vlop0=w0*w0dv1                !d23-1
        vlop1=w1*w0dv1
c        list=list4(lrd,lrg,lrs,lra)
c        wl=vlop0*(vint_ci(list)-2*vint_ci(list+1))  !1.1          !!!!!
c     :       -vlop1*vint_ci(list)
c          call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
          wl= vlop0-vlop1
         call trans_ijkl_intpos(lra,lrg,lrs,lrd,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
          wl= -vlop0*2
         call trans_ijkl_intpos(lra,lrs,lrg,lrd,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
      enddo
      return
!d1v(24-1)  ar(13)-
1124  iwdr=0
      do lrd=norb_frz+1,norb_dz
        imd=lsm_inn(lrd)
        if(imd.ne.jml) cycle
        iwdl=jud(lrd)
          w0dv1=w0_d1v(1)
        ni=mod(norb_dz-lrd,2)
        if(ni.eq.1) w0dv1=-w0dv1
        vlop0=w0*w0dv1                !d23-1
        vlop1=w1*w0dv1
c        list=list4(lrd,lrg,lrs,lra)
c        wl=vlop0*(vint_ci(list)-2*vint_ci(list+1))  !1.1          !!!!!
c     :       -vlop1*vint_ci(list)
c          call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
          wl= vlop0-vlop1
         call trans_ijkl_intpos(lra,lrg,lrs,lrd,nxo)
         call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
          wl= -vlop0*2
         call trans_ijkl_intpos(lra,lrs,lrg,lrd,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
      enddo
      return

!line=34:-b&l-b^l-a^r<-->sd(6),sdd(8),td(13),ttdd(15),dv(23),ddv(24)
1200  goto(10,10,10,10,10,1206,10,1208,10,10,10,10,
     :     1213,10,1215,10,10,10,10,10,10,10,1223,1224,10,10),lpok

1206  call dbl_sd_act_comp_g(7,lra)
      return
1208  call dbl_sdd_act_comp_g(7,lra)
      return
1213  call dbl_td_act_comp_g(7,lra)
      return
1215  call dbl_ttdd_act_comp_g(7,lra)
      return
!dv(23-1) ar(23)-
!dv(23-2) drl(33)-bl(23)-
1223  iwdr=0
      do lrd=norb_frz+1,norb_dz
        imd=lsm_inn(lrd)
        if(imd.ne.jml) cycle
        iwdl=jud(lrd)
        w0dv1=w0_dv(1)
        ni=mod(norb_dz-lrd,2)
        if(ni.eq.1) w0dv1=-w0dv1
        vlop0=w0*w0dv1                !d23-1
        vlop1=w1*w0dv1
c        list=list4(lrd,lrg,lrs,lra)
c        wl=vlop0*(vint_ci(list+2)-2.0d0*vint_ci(list+1)) !1.2      !!!!
c     :           -vlop1*vint_ci(list+2)
c          call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
          wl= vlop0-vlop1
         call trans_ijkl_intpos(lra,lrd,lrg,lrs,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
          wl= -vlop0*2
         call trans_ijkl_intpos(lra,lrs,lrg,lrd,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)

      enddo
      return
!d1v(24-1)  ar(13)-
1224  iwdr=0
      do lrd=norb_frz+1,norb_dz
        imd=lsm_inn(lrd)
        if(imd.ne.jml) cycle
        iwdl=jud(lrd)
          w0dv1=w0_d1v(1)
        ni=mod(norb_dz-lrd,2)
        if(ni.eq.1) w0dv1=-w0dv1
        vlop0=w0*w0dv1                !d23-1
        vlop1=w1*w0dv1
c        list=list4(lrd,lrg,lrs,lra)
c        wl=vlop0*(vint_ci(list+2)-2.0d0*vint_ci(list+1)) !1.2      !!!!
c     :           -vlop1*vint_ci(list+2)
c          call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
         call trans_ijkl_intpos(lra,lrd,lrg,lrs,nxo)
          wl= vlop0-vlop1
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
          wl= -vlop0*2
         call trans_ijkl_intpos(lra,lrs,lrg,lrd,nxo)
          call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
      enddo
      return

!line=35:-d&r^l-a^l<-->sd(6),sdd(8),td(13),ttdd(15),dv(23),ddv(24)
1300  goto(10,10,10,10,10,1306,10,1308,10,10,10,10,
     :     1313,10,1315,10,10,10,10,10,10,10,1323,1324,10,10),lpok

1306  call dbl_sd_act_comp_g(2,lra)
      return
1308  call dbl_sdd_act_comp_g(2,lra)
      return
1313  call dbl_td_act_comp_g(2,lra)
      return
1315  call dbl_ttdd_act_comp_g(2,lra)
      return
!dv(23-1) ar(23)-
!dv(23-2) drl(33)-bl(23)-
1323  iwdr=0
      do lrd=norb_frz+1,norb_dz
        imd=lsm_inn(lrd)
        if(imd.ne.jml) cycle
        iwdl=jud(lrd)
        w0dv1=w0_dv(1)
        ni=mod(norb_dz-lrd,2)
        if(ni.eq.1) w0dv1=-w0dv1
        vlop0=w0*w0dv1                !d23-1
c        list=list3(lrd,lra,lrs)
c        wl=vlop0*vint_ci(list)          !3.2         !!!!!
c        call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
         wl=vlop0
         call trans_ijkl_intpos(lra,lrs,lrd,lrs,nxo)
         call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
      enddo
      return
!d1v(24-1)  ar(13)-
1324  iwdr=0
      do lrd=norb_frz+1,norb_dz
        imd=lsm_inn(lrd)
        if(imd.ne.jml) cycle
        iwdl=jud(lrd)
          w0dv1=w0_d1v(1)
        ni=mod(norb_dz-lrd,2)
        if(ni.eq.1) w0dv1=-w0dv1
        vlop0=w0*w0dv1                !d23-1
c        list=list3(lrd,lra,lrg)
c        wl=vlop0*vint_ci(list)          !3.2         !!!!!
c        call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
         wl=vlop0
         call trans_ijkl_intpos(lra,lrg,lrd,lrg,nxo)
         call prodab_2(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper,nxo)
      enddo
      return

10    return
      end

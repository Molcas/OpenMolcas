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
c look for partial loops in active space drt and save them into disk
      subroutine guga_ploop(npl,maxplcon)
#include "drt_h.fh"
#include "files_gugaci.fh"
      common/count/mhsum,lp_count(22),mhlpmax
      common/lpdisk/idisk_lp,idisk_array(13)


      maxplcon=0
      mhlpmax=0
      idisk_lp=0
      idisk_array=0
      call idafile(luloop,1,idisk_array,13,idisk_lp)
      npl=0
      call vd_lp_search(npl)
      call dv_lp_search(npl)
      call dd_lp_search(npl)
      call dt_lp_search(npl)
      call ds_lp_search(npl)
      call tv_lp_search(npl)
      call td_lp_search(npl)
      call tt_lp_search(npl)
      call ts_lp_search(npl)
      call sv_lp_search(npl)
      call sd_lp_search(npl)
      call st_lp_search(npl)
      call ss_lp_search(npl)

      maxplcon=mhlpmax
      idisk_lp=0
      call idafile(luloop,1,idisk_array,13,idisk_lp)
      write(6,"(4x,a,i20)") "Max number of partial loops: ",mhlpmax
!      print*, "idisk_array"
!      print"(13i8)", idisk_array

      return
      end

      subroutine sv_lp_search(npl)
#include "drt_h.fh"
#include "pl_structure_h.fh"
      common/lpdisk/idisk_lp,idisk_array(13)
      common/count/mhsum,lp_count(22),mhlpmax
      mhsum=0
      lpblock=0
      lp_count(1:22)=0
      idisk_array(10)=idisk_lp

      imr=1
      do iml=1,ng_sm
        call act_lp_search(4,4,1)     !id=4
      enddo
      lpblock_sv=lpblock
      npl=npl+mhsum
      write(6,'(a15,2i10)')'sv:',lpblock,mhsum !,idisk_array(10),idisk_l
      return
      end

      subroutine sd_lp_search(npl)
#include "drt_h.fh"
#include "pl_structure_h.fh"
      common/count/mhsum,lp_count(22),mhlpmax
      common/lpdisk/idisk_lp,idisk_array(13)
      mhsum=0
      lpblock=0
      lp_count(1:22)=0
      idisk_array(11)=idisk_lp

      do iml=1,ng_sm
        do imr=1,ng_sm
          call act_lp_search(1,4,2)    !id=1
        enddo
      enddo
      lpblock_sd=lpblock
      npl=npl+mhsum
      write(6,'(a15,2i10)')'sd:',lpblock,mhsum
      return
      end

      subroutine st_lp_search(npl)
#include "drt_h.fh"
#include "pl_structure_h.fh"
      common/count/mhsum,lp_count(22),mhlpmax
      common/lpdisk/idisk_lp,idisk_array(13)
      mhsum=0
      lpblock=0
      lp_count(1:22)=0
      idisk_array(12)=idisk_lp

      do iml=1,ng_sm
        do imr=1,ng_sm
          call act_lp_search(3,4,3)   !id=3
        enddo
      enddo
      lpblock_st=lpblock
      npl=npl+mhsum
      write(6,'(a15,2i10)')'st:',lpblock,mhsum
      return
      end

      subroutine ss_lp_search(npl)
#include "drt_h.fh"
#include "pl_structure_h.fh"
      common/count/mhsum,lp_count(22),mhlpmax
      common/lpdisk/idisk_lp,idisk_array(13)
      mhsum=0
      lpblock=0
      lp_count(1:22)=0
      idisk_array(13)=idisk_lp

      do iml=1,ng_sm
        do imr=1,ng_sm
          call act_lp_search(3,4,4)    !id=3
        enddo
      enddo
      lpblock_ss=lpblock
      npl=npl+mhsum
      write(6,'(a15,2i10)')'ss:',lpblock,mhsum
      return
      end

      subroutine tv_lp_search(npl)
#include "drt_h.fh"
#include "pl_structure_h.fh"
      common/count/mhsum,lp_count(22),mhlpmax
      common/lpdisk/idisk_lp,idisk_array(13)
      mhsum=0
      lpblock=0
      lp_count(1:22)=0
      idisk_array(6)=idisk_lp


      imr=1
      do iml=1,ng_sm
        call act_lp_search(4,3,1)       !id=4
      enddo
      lpblock_tv=lpblock
      npl=npl+mhsum
      write(6,'(a15,2i10)')'tv:',lpblock,mhsum
      return
      end

      subroutine td_lp_search(npl)
#include "drt_h.fh"
#include "pl_structure_h.fh"
      common/count/mhsum,lp_count(22),mhlpmax
      common/lpdisk/idisk_lp,idisk_array(13)
      mhsum=0
      lpblock=0
      lp_count(1:22)=0
      idisk_array(7)=idisk_lp

      do iml=1,ng_sm
        do imr=1,ng_sm
          call act_lp_search(1,3,2)    !id=1
        enddo
      enddo
      lpblock_td=lpblock
      npl=npl+mhsum
      write(6,'(a15,2i10)')'td:',lpblock,mhsum
      return
      end

      subroutine tt_lp_search(npl)
#include "drt_h.fh"
#include "pl_structure_h.fh"
      common/count/mhsum,lp_count(22),mhlpmax
      common/lpdisk/idisk_lp,idisk_array(13)
      mhsum=0
      lpblock=0
      lp_count(1:22)=0
      idisk_array(8)=idisk_lp

      do iml=1,ng_sm
        do imr=1,ng_sm
          call act_lp_search(3,3,3)    !id=3
        enddo
      enddo
      lpblock_tt=lpblock
      npl=npl+mhsum
      write(6,'(a15,2i10)')'tt:',lpblock,mhsum

      return
      end

      subroutine ts_lp_search(npl)
#include "drt_h.fh"
#include "pl_structure_h.fh"
      common/count/mhsum,lp_count(22),mhlpmax
      common/lpdisk/idisk_lp,idisk_array(13)
      mhsum=0
      lpblock=0
      lp_count(1:22)=0
      idisk_array(9)=idisk_lp

      do iml=1,ng_sm
        do imr=1,ng_sm
          call act_lp_search(3,3,4)    !id=3
        enddo
      enddo
      lpblock_ts=lpblock
      npl=npl+mhsum
      write(6,'(a15,2i10)')'ts:',lpblock,mhsum
      return
      end

      subroutine dv_lp_search(npl)
#include "drt_h.fh"
#include "pl_structure_h.fh"
      common/count/mhsum,lp_count(22),mhlpmax
      common/lpdisk/idisk_lp,idisk_array(13)
      mhsum=0
      lpblock=0
      lp_count(1:22)=0
      idisk_array(2)=idisk_lp

      imr=1
      do iml=1,ng_sm
        call act_lp_search(1,2,1)      !id=1
      enddo
      lpblock_dv=lpblock
      npl=npl+mhsum
      write(6,'(a15,2i10)')'dv:',lpblock,mhsum
      return
      end

      subroutine dd_lp_search(npl)
#include "drt_h.fh"
#include "pl_structure_h.fh"
      common/count/mhsum,lp_count(22),mhlpmax
      common/lpdisk/idisk_lp,idisk_array(13)
      mhsum=0
      lpblock=0
      lp_count(1:22)=0
      idisk_array(3)=idisk_lp

      do iml=1,ng_sm
        do imr=1,ng_sm
          call act_lp_search(3,2,2)    !id=3
        enddo
      enddo
      lpblock_dd=lpblock
      npl=npl+mhsum
      write(6,'(a15,2i10)')'dd:',lpblock,mhsum

      return
      end

      subroutine dt_lp_search(npl)
#include "drt_h.fh"
#include "pl_structure_h.fh"
      common/count/mhsum,lp_count(22),mhlpmax
      common/lpdisk/idisk_lp,idisk_array(13)
      mhsum=0
      lpblock=0
      lp_count(1:22)=0
      idisk_array(4)=idisk_lp

      do iml=1,ng_sm
        do imr=1,ng_sm
          call act_lp_search(2,2,3)      !id=2
        enddo
      enddo
      lpblock_dt=lpblock
      npl=npl+mhsum
      write(6,'(a15,2i10)')'dt:',lpblock,mhsum
      return
      end

      subroutine ds_lp_search(npl)
#include "drt_h.fh"
#include "pl_structure_h.fh"
      common/count/mhsum,lp_count(22),mhlpmax
      common/lpdisk/idisk_lp,idisk_array(13)
      mhsum=0
      lpblock=0
      lp_count(1:22)=0
      idisk_array(5)=idisk_lp

      do iml=1,ng_sm
        do imr=1,ng_sm
          call act_lp_search(2,2,4)    !id=2
        enddo
      enddo
      lpblock_ds=lpblock
      npl=npl+mhsum
      write(6,'(a15,2i10)')'ds:',lpblock,mhsum
      return
      end

      subroutine vd_lp_search(npl)
#include "drt_h.fh"
#include "pl_structure_h.fh"
      common/count/mhsum,lp_count(22),mhlpmax
      common/lpdisk/idisk_lp,idisk_array(13)
      mhsum=0
      lpblock=0
      lp_count(1:22)=0
      idisk_array(1)=idisk_lp

      iml=1
      do imr=1,ng_sm
        call act_lp_search(2,1,2)       !id=2
      enddo
      lpblock_vd=lpblock
      npl=npl+mhsum
      write(6,'(a15,2i10)')'vd:',lpblock,mhsum
      return
      end

      subroutine act_lp_search(id,iptyl,iptyr)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      common /scratch/ tmpdir,len_str
      character*256 tmpdir
      call get_jp(iptyl,iml,ipael,0)
      call get_jp(iptyr,imr,ipaer,0)

      jpael=nu_ae(ipael)
      jpaer=nu_ae(ipaer)
      if(jpael.eq.0.or.jpaer.eq.0) return
      ide=0
      if(iptyl.eq.4.and.iptyr.eq.4) ide=1
      if(iptyl.eq.4.and.iptyr.eq.3) ide=2
      if(iptyl.eq.3.and.iptyr.eq.4) ide=2
      do jptyl=1,6            !1,6???
        jptyrend=6
        if(jptyl.eq.1) jptyrend=1
        jmlsta=1
        jmlend=ng_sm
        if(jptyl.eq.1) jmlsta=ns_sm
        if(jptyl.eq.1) jmlend=ns_sm
        do jml=jmlsta,jmlend
          call get_jp(jptyl,jml,jpadl,0)
          if(nu_ad(jpadl).eq.0) cycle
          jpad=jpadl
          jpae=jpael
          ipae=ipael
          call seg_drt()
          if(ndim.eq.0) cycle
          call copy_to_drtl()
          do jptyr=1,jptyrend
            jmrsta=1
            jmrend=ng_sm
            if(jptyr.eq.1) jmrsta=ns_sm
            if(jptyr.eq.1) jmrend=ns_sm
            do jmr=jmrsta,jmrend
              call get_jp(jptyr,jmr,jpadr,0)
              if(nu_ad(jpadr).eq.0) cycle
              jpad=jpadr
              jpae=jpaer
              ipae=ipaer
              call seg_drt()
              if(ndim.eq.0) cycle
              jpadlr=map_jplr(jptyl,jptyr)
              if(id.eq.1) call lp_head_in_dbl_1()
              if(id.eq.2) call lp_head_in_dbl_2()
              if(id.eq.3) call lp_head_in_dbl_3(ide)
              if(id.eq.4) call lp_head_in_dbl_4()
              if(jpadl.ne.jpadr) cycle
              if(id.eq.1) call lp_head_in_act_1()
              if(id.eq.2) call lp_head_in_act_2()
              if(id.eq.3) call lp_head_in_act_3(ide)
              if(id.eq.4) call lp_head_in_act_4()
            enddo
          enddo
        enddo
      enddo
      return
      end

      subroutine lp_head_in_act_1()      !for dv,td,sd
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      imlr=mul_tab(iml,imr)
      do lra=norb_dz+1,norb_inn
        lma=lsm_inn(lra)
        if(lma.ne.imlr) cycle
!line=1 a&r<-->a^r
        call head_ar_at_given_orb(mh,lra)
        call link_c1_to_given_orb_coe(mh,lra+1,norb_inn)
        if ( mh.eq.0 ) cycle
        call value_sort_ploop(mh,.true.,.true.,.false.)
        if(log_prod.eq.3) then
          linelp=1
          mhlp=mh
          nlg1=lra
          nlg2=0
          call ext_head_in_act()
        else
          call save_lp(1,mh,lra,0)
          lpblock=lpblock+1
        endif
      enddo

      do lrai=norb_dz+1,norb_inn
        lsmi=lsm_inn(lrai)
        do lraj=lrai+1,norb_inn
          lsmj=lsm_inn(lraj)
          lsmij=mul_tab(lsmi,lsmj)
          lmk=mul_tab(lsmij,imlr)
          do lrak=lraj+1,norb_inn
            lsmk=lsm_inn(lrak)
            if(lmk.ne.lsmk) cycle
            ijk=lrai-norb_frz+ngw2(lraj-norb_frz)+ngw3(lrak-norb_frz)
            intpos=intind_ijka(ijk)
!line=4  a&r--b&r--b^r<-->a^r
            call head_ar_at_given_orb(mh,lrai)
            call link_c1_to_given_orb(mh,lrai+1,lraj-1)
            call link_b4_at_given_orb(mh,lraj)
            logic_br(1:mh)=.true.
            call link_c2_to_given_orb(mh,lraj+1,lrak-1)
            call link_b2_at_given_orb(mh,lrak)
            call link_c1_to_given_orb(mh,lrak+1,norb_inn)
            if ( mh.eq.0 ) goto 7
            call value_sort_ploop(mh,.false.,.true.,.true.)
            if(log_prod.eq.3) then
              linelp=4
              mhlp=mh
              nlg1=intpos
              nlg2=0
              call ext_head_in_act()
            else
              call save_lp(4,mh,intpos,0)
              lpblock=lpblock+1
            endif
!line=7 a&r--b&l--b^l<-->a^r
7           call head_ar_at_given_orb(mh,lrai)
            call link_c1_to_given_orb(mh,lrai+1,lraj-1)
            call link_b3_at_given_orb(mh,lraj)
            logic_br(1:mh)=.true.
            call link_c2_to_given_orb(mh,lraj+1,lrak-1)
            call link_b1_at_given_orb(mh,lrak)
            call link_c1_to_given_orb(mh,lrak+1,norb_inn)
            if ( mh.eq.0 ) cycle
            call value_sort_ploop(mh,.false.,.true.,.true.)
            if(log_prod.eq.3) then
              linelp=7
              mhlp=mh
              nlg1=intpos
              nlg2=0
              call ext_head_in_act()
            else
              call save_lp(7,mh,intpos,0)
              lpblock=lpblock+1
            endif
          enddo
        enddo
      enddo
      do lrai=norb_dz+2,norb_inn
        lmai=lsm_inn(lrai)
        if(lmai.ne.imlr) cycle
        do lrak=norb_dz+1,lrai-1
!line=10 d&rr--b^r<-->a^r
          call head_drr_at_given_orb(mh,lrak)
          logic_br(1:mh)=.true.
          call link_c2_to_given_orb(mh,lrak+1,lrai-1)
          call link_b2_at_given_orb(mh,lrai)
          call link_c1_to_given_orb(mh,lrai+1,norb_inn)
          if ( mh.eq.0 ) goto 12
          call value_sort_ploop(mh,.false.,.true.,.true.)
          if(log_prod.eq.3) then
            linelp=10
            mhlp=mh
            nlg1=lrak
            nlg2=lrai
            call ext_head_in_act()
          else
            call save_lp(10,mh,lrak,lrai)
            lpblock=lpblock+1
          endif
!line=12 d&rl--b^l<-->a^r
12        call head_drl_at_given_orb(mh,lrak)
          call link_c2_to_given_orb(mh,lrak+1,lrai-1)
          call link_b1_at_given_orb(mh,lrai)
          call link_c1_to_given_orb(mh,lrai+1,norb_inn)
          if ( mh.eq.0 ) cycle
          call value_sort_ploop(mh,.false.,.true.,.true.)
          if(log_prod.eq.3) then
            linelp=12
            mhlp=mh
            nlg1=lrak
            nlg2=lrai
            call ext_head_in_act()
          else
            call save_lp(12,mh,lrak,lrai)
            lpblock=lpblock+1
          endif
        enddo
      enddo
      return
      end

      subroutine lp_head_in_act_2()      !for vd,dt,ds
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      imlr=mul_tab(iml,imr)

      do lra=norb_dz+1,norb_inn
        lma=lsm_inn(lra)
        if(lma.ne.imlr) cycle
        do lrd=lra+1,norb_inn
!line=2 a&r--d&l^r<-->a^l
          call head_ar_at_given_orb(mh,lra)
          call link_c1_to_given_orb(mh,lra+1,lrd-1)
          call link_d10_at_given_orb(mh,lrd)
          call link_c1_to_given_orb(mh,lrd+1,norb_inn)
          if ( mh.eq.0 ) cycle
          call value_sort_ploop(mh,.false.,.true.,.false.)
          if(log_prod.eq.3) then
            linelp=2
            mhlp=mh
            nlg1=lra
            nlg2=lrd
            !call print_lp()
            call ext_head_in_act()
          else
            linelp=2
            mhlp=mh
            nlg1=lra
            nlg2=lrd
            !call print_lp()
            call save_lp(2,mh,lra,lrd)
            lpblock=lpblock+1
          endif
        enddo
      enddo

      do lrai=norb_dz+1,norb_inn
        lsmi=lsm_inn(lrai)
        do lraj=lrai+1,norb_inn
          lsmj=lsm_inn(lraj)
          lsmij=mul_tab(lsmi,lsmj)
          lmk=mul_tab(lsmij,imlr)
          do lrak=lraj+1,norb_inn
            lsmk=lsm_inn(lrak)
            if(lmk.ne.lsmk) cycle
            ijk=lrai-norb_frz+ngw2(lraj-norb_frz)+ngw3(lrak-norb_frz)
            intpos=intind_ijka(ijk)
!line=6  a&r--b&l--b^r<-->a^l
            call head_ar_at_given_orb(mh,lrai)
            call link_c1_to_given_orb(mh,lrai+1,lraj-1)
            call link_b3_at_given_orb(mh,lraj)
            logic_br(1:mh)=.true.
            call link_c2_to_given_orb(mh,lraj+1,lrak-1)
            call link_b2_at_given_orb(mh,lrak)
            call link_c1_to_given_orb(mh,lrak+1,norb_inn)
            if ( mh.eq.0 ) cycle
            call value_sort_ploop(mh,.false.,.true.,.true.)
            if(log_prod.eq.3) then
              linelp=6
              mhlp=mh
              nlg1=intpos
              nlg2=0
              !call print_lp()
              call ext_head_in_act()
            else
              call save_lp(6,mh,intpos,0)
              lpblock=lpblock+1
            endif
          enddo
        enddo
      enddo

      do lrai=norb_dz+2,norb_inn
        lmai=lsm_inn(lrai)
        if(lmai.ne.imlr) cycle
        do lrak=norb_dz+1,lrai-1
!line=11 d&rl--b^r<-->a^l
          call head_drl_at_given_orb(mh,lrak)
          call link_c2_to_given_orb(mh,lrak+1,lrai-1)
          call link_b2_at_given_orb(mh,lrai)
          call link_c1_to_given_orb(mh,lrai+1,norb_inn)
          if ( mh.eq.0 ) cycle
          call value_sort_ploop(mh,.false.,.true.,.true.)
          if(log_prod.eq.3) then
            linelp=11
            mhlp=mh
            nlg1=lrak
            nlg2=lrai
            !call print_lp()
            call ext_head_in_act()
          else
            call save_lp(11,mh,lrak,lrai)
            lpblock=lpblock+1
          endif
        enddo
      enddo
      return
      end

      subroutine lp_head_in_act_3(ide)      !for ide=0:dd,tt,ide=1:ss,id
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      imlr=mul_tab(iml,imr)

      do lra=norb_dz+1,norb_inn
        lma=lsm_inn(lra)
!line=9 d&r&l-
          call head_drl_at_given_orb(mh,lra)
          call link_c2_to_given_orb(mh,lra+1,norb_inn)
          if ( mh.eq.0 ) cycle
          call value_sort_ploop(mh,.false.,.true.,.true.)
          call save_lp(9,mh,lra,0)
          lpblock=lpblock+1
      enddo
      do lrai=norb_dz+1,norb_inn
        lsmi=lsm_inn(lrai)
        do lraj=lrai+1,norb_inn
          lsmj=lsm_inn(lraj)
          lsmij=mul_tab(lsmi,lsmj)
          if(lsmij.ne.imlr) cycle
!line=5 a&r-b&l-
          call head_ar_at_given_orb(mh,lrai)
          call link_c1_to_given_orb(mh,lrai+1,lraj-1)
          call link_b3_at_given_orb(mh,lraj)
          logic_br(1:mh)=.true.
          call link_c2_to_given_orb(mh,lraj+1,norb_inn)
          if ( mh.eq.0 ) cycle
        if(ide.ne.0) then
          goto(22,33),ide
22          call value_sort_ploop(mh,.false.,.true.,.false.)    !ss
            goto 55
33          call value_sort_ploop(mh,.false.,.false.,.true.)    !st,ts
            goto 55
        endif
        call value_sort_ploop(mh,.false.,.true.,.true.)
55      call save_lp(5,mh,lrai,lraj)
        lpblock=lpblock+1
        enddo
      enddo
      return
      end

      subroutine lp_head_in_act_4()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      imlr=mul_tab(iml,imr)

      do lra=norb_dz+1,norb_inn
!line=8 d&rr-
        call head_drr_at_given_orb(mh,lra)
        logic_br(1:mh)=.true.
        call link_c2_to_given_orb(mh,lra+1,norb_inn)
        if ( mh.eq.0 ) cycle
        call value_sort_ploop(mh,.false.,.true.,.true.)
        if(log_prod.eq.3) then
          linelp=8
          mhlp=mh
          nlg1=lra
          nlg2=0
          call ext_head_in_act()
        else
          call save_lp(8,mh,lra,0)
          lpblock=lpblock+1
        endif
      enddo
      do lrai=norb_dz+1,norb_inn
        lsmi=lsm_inn(lrai)
        do lraj=lrai+1,norb_inn
          lsmj=lsm_inn(lraj)
          lsmij=mul_tab(lsmi,lsmj)
          if(lsmij.ne.imlr) cycle
          call head_ar_at_given_orb(mh,lrai)
          call link_c1_to_given_orb(mh,lrai+1,lraj-1)
          call link_b4_at_given_orb(mh,lraj)
          logic_br(1:mh)=.true.
          call link_c2_to_given_orb(mh,lraj+1,norb_inn)
          if ( mh.eq.0 ) cycle
          call value_sort_ploop(mh,.false.,.true.,.true.)
          if(log_prod.eq.3) then
            linelp=3
            mhlp=mh
            nlg1=lrai
            nlg2=lraj
            call ext_head_in_act()
          else
            call save_lp(3,mh,lrai,lraj)
            lpblock=lpblock+1
          endif
        enddo
      enddo
      return
      end

      subroutine lp_head_in_dbl_1()      !for dv,sd,td
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      imlr=mul_tab(iml,imr)
      jmlr=mul_tab(jml,jmr)
      lsmact=mul_tab(imlr,jmlr)
      lsta=norb_dz+1
      lend=norb_inn

      if(imlr.eq.jmlr) then       !!! imlr.eq.1

!line=13 -c'-
        call link_c1_to_given_orb_coe(mh,lsta,lend)
        if(mh.eq.0) goto 15
        call value_sort_ploop(mh,.true.,.true.,.true.)
        call save_lp(13,mh,0,0)
        lpblock=lpblock+1
      endif

!line=15 -b^l-
15    do lra=norb_dz+1,norb_inn
        lsma=lsm_inn(lra)

        if(lsma.ne.lsmact) cycle

        if(jml.ne.jmr) goto 1502
        logic_br(1)=.false.
        call link_c2_to_given_orb(mh,norb_dz+1,lra-1)
        call link_b1_at_given_orb(mh,lra)
        call link_c1_to_given_orb(mh,lra+1,norb_inn)
        if(mh.eq.0) goto 1502
        call value_sort_ploop(mh,.false.,.true.,.true.)
        call save_lp(15,mh,lra,1)
        lpblock=lpblock+1
1502    logic_br(1)=.true.
        call link_c2_to_given_orb(mh,norb_dz+1,lra-1)
        call link_b1_at_given_orb(mh,lra)        !b^l
        call link_c1_to_given_orb(mh,lra+1,norb_inn)
        if(mh.eq.0) goto 1601
        call value_sort_ploop(mh,.false.,.true.,.true.)
        call save_lp(15,mh,lra,2)
        lpblock=lpblock+1

!line=16 -b^r-
1601    if(jml.ne.jmr) goto 1602
        logic_br(1)=.false.
        call link_c2_to_given_orb(mh,norb_dz+1,lra-1)
        call link_b2_at_given_orb(mh,lra)
        call link_c1_to_given_orb(mh,lra+1,norb_inn)
        if(mh.eq.0) goto 1602
        call value_sort_ploop(mh,.false.,.true.,.true.)
        call save_lp(16,mh,lra,1)
        lpblock=lpblock+1

1602    logic_br(1)=.true.
        call link_c2_to_given_orb(mh,norb_dz+1,lra-1)
        call link_b2_at_given_orb(mh,lra)
        call link_c1_to_given_orb(mh,lra+1,norb_inn)
        if(mh.eq.0) cycle
        call value_sort_ploop(mh,.false.,.true.,.true.)
        call save_lp(16,mh,lra,2)
        lpblock=lpblock+1

      enddo

!line=20 -b&r-b^r-
      do lrai=norb_dz+1,norb_inn-1
        lsmi=lsm_inn(lrai)
        do lraj=lrai+1,norb_inn
          lsmj=lsm_inn(lraj)
          lsmij=mul_tab(lsmi,lsmj)
         if(lsmij.ne.lsmact) cycle
          jk=ngw2(lrai-norb_frz)+ngw3(lraj-norb_frz)

          call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
          call link_b4_at_given_orb(mh,lrai)        !b&r
          logic_br(1:mh)=.true.
          call link_c2_to_given_orb(mh,lrai+1,lraj-1)
          call link_b2_at_given_orb(mh,lraj)        !b^r
          call link_c1_to_given_orb(mh,lraj+1,norb_inn)
          if(mh.eq.0) goto 22
          call value_sort_ploop(mh,.false.,.true.,.true.)
          call save_lp(20,mh,jk,0)
          lpblock=lpblock+1
!line=22 -b&l-b^l-
22        call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
          call link_b3_at_given_orb(mh,lrai)        !b&l
          logic_br(1:mh)=.true.
          call link_c2_to_given_orb(mh,lrai+1,lraj-1)
          call link_b1_at_given_orb(mh,lraj)        !b^l
          call link_c1_to_given_orb(mh,lraj+1,norb_inn)
          if(mh.eq.0) cycle
          call value_sort_ploop(mh,.false.,.true.,.true.)
          call save_lp(22,mh,jk,0)
          lpblock=lpblock+1
        enddo
      enddo
      return
      end

      subroutine lp_head_in_dbl_1_mrpt2()    !for dv,sd,td       !200709

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      imlr=mul_tab(iml,imr)
      jmlr=mul_tab(jml,jmr)
      lsmact=mul_tab(imlr,jmlr)
      lsta=norb_dz+1
      lend=norb_inn

      goto(100,200,300,400),jpadlrel(jpadlr)              !200709

100   if(imlr.eq.jmlr) then               !!! imlr.eq.1
!line=13 -c'-
        call link_c1_to_given_orb_coe(mh,lsta,lend)
        if(mh.eq.0) goto 15
        call value_sort_ploop(mh,.true.,.true.,.true.)
        linelp=13
        mhlp=mh
        nlg1=0
        nlg2=0
        call ext_head_in_dbl()
       endif

!line=20 -b&r-b^r-
      do lrai=norb_dz+1,norb_inn-1
        lsmi=lsm_inn(lrai)
        do lraj=lrai+1,norb_inn
          lsmj=lsm_inn(lraj)
          lsmij=mul_tab(lsmi,lsmj)
              if(lsmij.ne.lsmact) cycle
          jk=ngw2(lrai-norb_frz)+ngw3(lraj-norb_frz)

          call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
          call link_b4_at_given_orb(mh,lrai)        !b&r
          logic_br(1:mh)=.true.
          call link_c2_to_given_orb(mh,lrai+1,lraj-1)
          call link_b2_at_given_orb(mh,lraj)        !b^r
          call link_c1_to_given_orb(mh,lraj+1,norb_inn)
          if(mh.eq.0) goto 22
          call value_sort_ploop(mh,.false.,.true.,.true.)
          linelp=20
          mhlp=mh
          nlg1=jk
          nlg2=0
          call ext_head_in_dbl()
!line=22 -b&l-b^l-
22        call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
          call link_b3_at_given_orb(mh,lrai)        !b&l
          logic_br(1:mh)=.true.
          call link_c2_to_given_orb(mh,lrai+1,lraj-1)
          call link_b1_at_given_orb(mh,lraj)        !b^l
          call link_c1_to_given_orb(mh,lraj+1,norb_inn)
          if(mh.eq.0) cycle
              call value_sort_ploop(mh,.false.,.true.,.true.)
          linelp=22
          mhlp=mh
          nlg1=jk
          nlg2=0
          call ext_head_in_dbl()
         enddo
      enddo
      return

300        continue
!line=15 -b^l-
15    do lra=norb_dz+1,norb_inn
        lsma=lsm_inn(lra)

        if(lsma.ne.lsmact) cycle

        if(jml.ne.jmr) goto 1502
        logic_br(1)=.false.
        call link_c2_to_given_orb(mh,norb_dz+1,lra-1)
        call link_b1_at_given_orb(mh,lra)
        call link_c1_to_given_orb(mh,lra+1,norb_inn)
        if(mh.eq.0) goto 1502
        call value_sort_ploop(mh,.false.,.true.,.true.)
        linelp=15
        mhlp=mh
        nlg1=lra
        nlg2=1
        call ext_head_in_dbl()
!        write(*,*) mh,linelp
1502    logic_br(1)=.true.
        call link_c2_to_given_orb(mh,norb_dz+1,lra-1)
        call link_b1_at_given_orb(mh,lra)        !b^l
        call link_c1_to_given_orb(mh,lra+1,norb_inn)
        if(mh.eq.0) goto 1601
        call value_sort_ploop(mh,.false.,.true.,.true.)
        linelp=15
        mhlp=mh
        nlg1=lra
        nlg2=2
        call ext_head_in_dbl()
!        write(*,*) mh,linelp
      enddo
      return

400       continue
!line=16 -b^r-
1601  do lra=norb_dz+1,norb_inn
        lsma=lsm_inn(lra)
        if(jml.ne.jmr) goto 1602
        logic_br(1)=.false.
        call link_c2_to_given_orb(mh,norb_dz+1,lra-1)
        call link_b2_at_given_orb(mh,lra)
        call link_c1_to_given_orb(mh,lra+1,norb_inn)
        if(mh.eq.0) goto 1602
        call value_sort_ploop(mh,.false.,.true.,.true.)
        linelp=16
        mhlp=mh
        nlg1=lra
        nlg2=1
        call ext_head_in_dbl()
!        write(*,*) mh,linelp

1602    logic_br(1)=.true.
        call link_c2_to_given_orb(mh,norb_dz+1,lra-1)
        call link_b2_at_given_orb(mh,lra)
        call link_c1_to_given_orb(mh,lra+1,norb_inn)
        if(mh.eq.0) cycle
        call value_sort_ploop(mh,.false.,.true.,.true.)
        linelp=16
        mhlp=mh
        nlg1=lra
        nlg2=2
        call ext_head_in_dbl()
!        write(*,*) mh,linelp

      enddo

200   return
      end

      subroutine lp_head_in_dbl_2()      !for vd,dt,ds
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      imlr=mul_tab(iml,imr)
      jmlr=mul_tab(jml,jmr)
      lsmact=mul_tab(imlr,jmlr)
      lsta=norb_dz+1
      lend=norb_inn

      if(imlr.eq.jmlr) then       !!! imlr.eq.1

!line=13 -c'-
        call link_c1_to_given_orb_coe(mh,lsta,lend)
        if(mh.eq.0) goto 16
        call value_sort_ploop(mh,.true.,.true.,.true.)
        if(log_prod.eq.3) then
          ! mrpt2, do not save loop
          linelp=13
          mhlp=mh
          nlg1=0
          nlg2=0
          call ext_head_in_dbl()
        else
          call save_lp(13,mh,0,0)
          lpblock=lpblock+1
        endif
      endif

!line=16 -b^r-
16    do lra=norb_dz+1,norb_inn
        lsma=lsm_inn(lra)

        if(lsma.ne.lsmact) goto 19

        if(jml.ne.jmr) goto 1602
        logic_br(1)=.false.
        call link_c2_to_given_orb(mh,norb_dz+1,lra-1)
        call link_b2_at_given_orb(mh,lra)        !b^r
        call link_c1_to_given_orb(mh,lra+1,norb_inn)
        if(mh.eq.0) goto 1602
        call value_sort_ploop(mh,.false.,.true.,.true.)
        if(log_prod.eq.3) then
          linelp=16
          mhlp=mh
          nlg1=lra
          nlg2=1
          call ext_head_in_dbl()
        else
          call save_lp(16,mh,lra,1)
          lpblock=lpblock+1
        endif

1602    logic_br(1)=.true.
        call link_c2_to_given_orb(mh,norb_dz+1,lra-1)
        call link_b2_at_given_orb(mh,lra)        !b^r
        call link_c1_to_given_orb(mh,lra+1,norb_inn)
        if(mh.eq.0) goto 19
        call value_sort_ploop(mh,.false.,.true.,.true.)
        if(log_prod.eq.3) then
          linelp=16
          mhlp=mh
          nlg1=lra
          nlg2=2
          call ext_head_in_dbl()
        else
          call save_lp(16,mh,lra,2)
          lpblock=lpblock+1
        endif

!line=19 -d&r^l-
19      call link_c1_to_given_orb(mh,norb_dz+1,lra-1)
        call link_d10_at_given_orb(mh,lra)
        call link_c1_to_given_orb(mh,lra+1,norb_inn)
        if(mh.eq.0) cycle
        call value_sort_ploop(mh,.true.,.true.,.true.)
        if(log_prod.eq.3) then
          linelp=19
          mhlp=mh
          nlg1=lra
          nlg2=0
          call ext_head_in_dbl()
        else
          call save_lp(19,mh,lra,0)
          lpblock=lpblock+1
        endif
      enddo

      do lrai=norb_dz+1,norb_inn-1
        lsmi=lsm_inn(lrai)
        do lraj=lrai+1,norb_inn
          lsmj=lsm_inn(lraj)
          lsmij=mul_tab(lsmi,lsmj)
          if(lsmij.ne.lsmact) cycle
          jk=ngw2(lrai-norb_frz)+ngw3(lraj-norb_frz)

!line=21 -b&l-b^r-
          call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
          call link_b3_at_given_orb(mh,lrai)        !b&l
          logic_br(1:mh)=.true.
          call link_c2_to_given_orb(mh,lrai+1,lraj-1)
          call link_b2_at_given_orb(mh,lraj)        !b^r
          call link_c1_to_given_orb(mh,lraj+1,norb_inn)
          if(mh.eq.0) cycle
          call value_sort_ploop(mh,.false.,.true.,.true.)
          if(log_prod.eq.3) then
            linelp=21
            mhlp=mh
            nlg1=jk
            nlg2=0
            call ext_head_in_dbl()
          else
            call save_lp(21,mh,jk,0)
            lpblock=lpblock+1
          endif
        enddo
      enddo

      return
      end

      subroutine lp_head_in_dbl_2_mrpt2()          !for vd,dt,ds
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      imlr=mul_tab(iml,imr)
      jmlr=mul_tab(jml,jmr)
      lsmact=mul_tab(imlr,jmlr)
      lsta=norb_dz+1
      lend=norb_inn

      goto(100,200,300,400),jpadlrel(jpadlr)        !200709

400   return

200   if(imlr.eq.jmlr) then               !!! imlr.eq.1

!line=13 -c'-
        call link_c1_to_given_orb_coe(mh,lsta,lend)
        if(mh.eq.0) goto 16
        call value_sort_ploop(mh,.true.,.true.,.true.)
        linelp=13
        mhlp=mh
        nlg1=0
        nlg2=0
        !call print_lp()
        call ext_head_in_dbl()

       endif
       return

300   continue
!line=16 -b^r-
16    do lra=norb_dz+1,norb_inn
        lsma=lsm_inn(lra)

        if(lsma.ne.lsmact) cycle

        if(jml.ne.jmr) goto 1602
        logic_br(1)=.false.
        call link_c2_to_given_orb(mh,norb_dz+1,lra-1)
        call link_b2_at_given_orb(mh,lra)        !b^r
        call link_c1_to_given_orb(mh,lra+1,norb_inn)
        if(mh.eq.0) goto 1602
        call value_sort_ploop(mh,.false.,.true.,.true.)
        linelp=16
        mhlp=mh
        nlg1=lra
        nlg2=1
        !call print_lp()
        call ext_head_in_dbl()

1602    logic_br(1)=.true.
        call link_c2_to_given_orb(mh,norb_dz+1,lra-1)
        call link_b2_at_given_orb(mh,lra)        !b^r
        call link_c1_to_given_orb(mh,lra+1,norb_inn)
        if(mh.eq.0) cycle
        call value_sort_ploop(mh,.false.,.true.,.true.)
        linelp=16
        mhlp=mh
        nlg1=lra
        nlg2=2
        !call print_lp()
        call ext_head_in_dbl()
      enddo
      return

100   do lra=norb_dz+1,norb_inn
        lsma=lsm_inn(lra)
!line=19 -d&l^r-
        call link_c1_to_given_orb(mh,norb_dz+1,lra-1)
        call link_d10_at_given_orb(mh,lra)
        call link_c1_to_given_orb(mh,lra+1,norb_inn)
        if(mh.eq.0) cycle
        call value_sort_ploop(mh,.true.,.true.,.true.)
        linelp=19
        mhlp=mh
        nlg1=lra
        nlg2=0
        !call print_lp()
        call ext_head_in_dbl()
      enddo

      do lrai=norb_dz+1,norb_inn-1
        lsmi=lsm_inn(lrai)
        do lraj=lrai+1,norb_inn
          lsmj=lsm_inn(lraj)
          lsmij=mul_tab(lsmi,lsmj)
          if(lsmij.ne.lsmact) cycle
          jk=ngw2(lrai-norb_frz)+ngw3(lraj-norb_frz)

!line=21 -b&l-b^r-
          call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
          call link_b3_at_given_orb(mh,lrai)        !b&l
          logic_br(1:mh)=.true.
          call link_c2_to_given_orb(mh,lrai+1,lraj-1)
          call link_b2_at_given_orb(mh,lraj)        !b^r
          call link_c1_to_given_orb(mh,lraj+1,norb_inn)
          if(mh.eq.0) cycle
          call value_sort_ploop(mh,.false.,.true.,.true.)
          linelp=21
          mhlp=mh
          nlg1=jk
          nlg2=0
          !call print_lp()
          call ext_head_in_dbl()
        enddo
      enddo
      return
      end

      subroutine lp_head_in_dbl_3(ide)      !for ide=0:dd,tt,ide=1:ss,id
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      imlr=mul_tab(iml,imr)
      jmlr=mul_tab(jml,jmr)
      lsmact=mul_tab(imlr,jmlr)
      lsta=norb_dz+1
      lend=norb_inn

      if(imlr.eq.jmlr) then       !!! imlr.eq.1

!line=14 -c"-
        if(jml.ne.jmr) goto 1402
        logic_br(1)=.false.
        call link_c2_to_given_orb(mh,lsta,lend)
        if(mh.eq.0) goto 1402
        call value_sort_ploop(mh,.true.,.true.,.true.)
        call save_lp(14,mh,0,1)
        lpblock=lpblock+1

1402    logic_br(1)=.true.
        call link_c2_to_given_orb(mh,lsta,lend)
        if(mh.eq.0) goto 15
        call value_sort_ploop(mh,.true.,.true.,.true.)
        call save_lp(14,mh,0,2)
        lpblock=lpblock+1

      endif

15    do lra=norb_dz+1,norb_inn
        lsma=lsm_inn(lra)
        if(lsma.ne.lsmact) cycle

!line=17 -b&l-
        call link_c1_to_given_orb(mh,norb_dz+1,lra-1)
        call link_b3_at_given_orb(mh,lra)        !b&l
        logic_br(1:mh)=.true.
        call link_c2_to_given_orb(mh,lra+1,norb_inn)
        if(mh.eq.0) cycle
        if(ide.ne.0) then
          goto(1,2),ide
1           call value_sort_ploop(mh,.false.,.true.,.false.)    !ss
            goto 5
2           call value_sort_ploop(mh,.false.,.false.,.true.)    !st,ts
            goto 5
        endif
        call value_sort_ploop(mh,.false.,.true.,.true.)
5       call save_lp(17,mh,lra,0)
        lpblock=lpblock+1

      enddo
      return
      end

      subroutine lp_head_in_dbl_4()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      imlr=mul_tab(iml,imr)
      jmlr=mul_tab(jml,jmr)
      lsmact=mul_tab(imlr,jmlr)
      lsta=norb_dz+1
      lend=norb_inn

      if(imlr.eq.jmlr) then       !!! imlr.eq.1

!line=14 -c"-
        if(jml.ne.jmr) goto 1402
        logic_br(1)=.false.
        call link_c2_to_given_orb(mh,lsta,lend)
        if(mh.eq.0) goto 1402
        call value_sort_ploop(mh,.true.,.true.,.true.)
        call save_lp(14,mh,0,1)
        lpblock=lpblock+1

1402    logic_br(1)=.true.
        call link_c2_to_given_orb(mh,lsta,lend)
        if(mh.eq.0) goto 15
        call value_sort_ploop(mh,.true.,.true.,.true.)
        call save_lp(14,mh,0,2)
        lpblock=lpblock+1

        endif
15      do lra=norb_dz+1,norb_inn
        lsma=lsm_inn(lra)
        if(lsma.ne.lsmact) cycle

!line=18 -b&r-
        call link_c1_to_given_orb(mh,norb_dz+1,lra-1)
        call link_b4_at_given_orb(mh,lra)        !b&r
        logic_br(1:mh)=.true.
        call link_c2_to_given_orb(mh,lra+1,norb_inn)
        if(mh.eq.0) cycle
        call value_sort_ploop(mh,.false.,.true.,.true.)
        call save_lp(18,mh,lra,0)
        lpblock=lpblock+1

      enddo

      return
      end

      subroutine lp_head_in_dbl_4_mrpt2()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      imlr=mul_tab(iml,imr)
      jmlr=mul_tab(jml,jmr)
      lsmact=mul_tab(imlr,jmlr)
      lsta=norb_dz+1
      lend=norb_inn

      goto(100,200,300,400),jpadlrel(jpadlr)         !200709

200   return
300   return
400      if(imlr.eq.jmlr) then               !!! imlr.eq.1

!line=14 -c"-
        if(jml.ne.jmr) goto 1402
        logic_br(1)=.false.
        call link_c2_to_given_orb(mh,lsta,lend)
        if(mh.eq.0) goto 1402
        call value_sort_ploop(mh,.true.,.true.,.true.)
        linelp=14
        mhlp=mh
        nlg1=0
        nlg2=1
        call ext_head_in_dbl()

1402    logic_br(1)=.true.
        call link_c2_to_given_orb(mh,lsta,lend)
        if(mh.eq.0) goto 15
        call value_sort_ploop(mh,.true.,.true.,.true.)
        linelp=14
        mhlp=mh
        nlg1=0
        nlg2=2
        call ext_head_in_dbl()

        endif
       return

100   continue

15      do lra=norb_dz+1,norb_inn
        lsma=lsm_inn(lra)
        if(lsma.ne.lsmact) cycle

!line=18 -b&r-
        call link_c1_to_given_orb(mh,norb_dz+1,lra-1)
        call link_b4_at_given_orb(mh,lra)        !b&r
        logic_br(1:mh)=.true.
        call link_c2_to_given_orb(mh,lra+1,norb_inn)
        if(mh.eq.0) cycle
        call value_sort_ploop(mh,.false.,.true.,.true.)
        linelp=18
        mhlp=mh
        nlg1=lra
        nlg2=0
        call ext_head_in_dbl()

        enddo

      return
      end

      subroutine print_lp()
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "files_gugaci.fh"
      common/count/mhsum,lp_count(22),mhlpmax
      common/lpdisk/idisk_lp,idisk_array(13)
      dimension info(10)
      line=linelp
      mh=mhlp
      lg1=nlg1
      lg2=nlg2

      mhsum=mhsum+mh
      lp_count(line)=lp_count(line)+1
      if((line.eq.14.or.line.eq.15.or.line.eq.16).and.lg2.eq.1) then
        mhsum=mhsum-mh
        lp_count(line)=lp_count(line)-1
      endif
      info=0
!===========================================================
      jpml=mul_tab(jml,ns_sm)
      jpmr=mul_tab(jmr,ns_sm)
      write(20,'(2x,10i8)')
     :   line,iml,imr,jpml,jpmr,jpadlr,mtype,mh,lg1,lg2
      write(20,"(7f12.6)")vplpnew_w0(1:mtype)
      write(20,"(7f12.6)")vplpnew_w1(1:mtype)
      write(20,"(8i4)")nstaval(1:mtype)
      write(20,"(8i4)")nvalue(1:mtype)
      write(20,"(8i6)")lpnew_lwei(1:mh)
      write(20,"(8i6)")lpnew_rwei(1:mh)
!===========================================================
      if(line.gt.12) goto 100
!upwei_record
      do m=1,mh
        jph=lpnew_head(m)
        jhyl=jphyl(jph)
        jhyr=jphy(jph)
        in=ihyl(jhyl)
        write(20,"(8(1x,i4))")jph,jhyl,jhyr,in,ihyl(jhyl+1:jhyl+in),
     *    ihy(jhyr+1:jhyr+in)
      enddo
100   if(line.ne.1.and.line.ne.13) return
!coe_record
      do m=1,mh
        write(20,*)lpnew_coe(norb_dz+1:norb_inn,m)
      enddo

      end

      subroutine save_lp(line,mh,lg1,lg2)
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "files_gugaci.fh"
      common/count/mhsum,lp_count(22),mhlpmax
      common/lpdisk/idisk_lp,idisk_array(13)
      dimension info(10)

      if(mh.gt.mhlpmax) mhlpmax=mh
      mhsum=mhsum+mh
      lp_count(line)=lp_count(line)+1
      if((line.eq.14.or.line.eq.15.or.line.eq.16).and.lg2.eq.1) then
        mhsum=mhsum-mh
        lp_count(line)=lp_count(line)-1
      endif
      info=0
!===========================================================
!      jpml=mul_tab(jml,ns_sm)
!      jpmr=mul_tab(jmr,ns_sm)
!      write(200,'(2x,10i8)')
!     :   line,iml,imr,jpml,jpmr,jpadlr,mtype,mh,lg1,lg2
!      write(200,*)vplpnew_w0(1:mtype)
!      write(200,*)vplpnew_w1(1:mtype)
!      write(200,*)nstaval(1:mtype)
!      write(200,*)nvalue(1:mtype)
!      write(200,*)lpnew_lwei(1:mh)
!      write(200,*)lpnew_rwei(1:mh)
!===========================================================
      jpml=mul_tab(jml,ns_sm)
      jpmr=mul_tab(jmr,ns_sm)
      info(1)=line
      info(2)=iml
      info(3)=imr
      info(4)=jpml
      info(5)=jpmr
      info(6)=jpadlr
      info(7)=mtype
      info(8)=mh
      info(9)=lg1
      info(10)=lg2
      call idafile(luloop,1,info,10,idisk_lp)
c      print*, "in save_lp, linelp",line,idisk_lp,mh,mtype
      call ddafile(luloop,1,vplpnew_w0,mtype,idisk_lp)
      call ddafile(luloop,1,vplpnew_w1,mtype,idisk_lp)
      call idafile(luloop,1,nstaval,mtype,idisk_lp)
      call idafile(luloop,1,nvalue,mtype,idisk_lp)
      call idafile(luloop,1,lpnew_lwei,mh,idisk_lp)
      call idafile(luloop,1,lpnew_rwei,mh,idisk_lp)
      if(line.gt.12) goto 100
      do m=1,mh
        jph=lpnew_head(m)
        jhyl=jphyl(jph)
        jhyr=jphy(jph)
        in=ihyl(jhyl)
        call idafile(luloop,1,[in],1,idisk_lp)
        call idafile(luloop,1,ihyl(jhyl+1:jhyl+in),in,idisk_lp)
        call idafile(luloop,1,ihy(jhyr+1:jhyr+in),in,idisk_lp)
c        write(20)in,ihyl(jhyl+1:jhyl+in),ihy(jhyr+1:jhyr+in)
      enddo
100   if(line.ne.1.and.line.ne.13) return
c      print*, "in read_lp, write coe",idisk_lp,norb_inn-norb_dz
!coe_record
      lenw=norb_inn-norb_dz
      do m=1,mh
        call idafile(luloop,1,lpnew_coe(norb_dz+1:norb_inn,m),lenw,
     &               idisk_lp)
      enddo
      return
      end

      subroutine read_lp()
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
#include "files_gugaci.fh"
      common/lpdisk/idisk_lp,idisk_array(13)
      dimension info(10)

      call idafile(luloop,2,info,10,idisk_lp)
      linelp=info(1)
      iml=info(2)
      imr=info(3)
      jml=info(4)
      jmr=info(5)
      jpadlr=info(6)
      mtype=info(7)
      mhlp=info(8)
      nlg1=info(9)
      nlg2=info(10)
      call ddafile(luloop,2,vplpnew_w0,mtype,idisk_lp)
      call ddafile(luloop,2,vplpnew_w1,mtype,idisk_lp)
      call idafile(luloop,2,nstaval,mtype,idisk_lp)
      call idafile(luloop,2,nvalue,mtype,idisk_lp)
      call idafile(luloop,2,lpnew_lwei,mhlp,idisk_lp)
      call idafile(luloop,2,lpnew_rwei,mhlp,idisk_lp)
!===========================================================

      if(linelp.gt.12) goto 100
!upwei_read
      ihypos=1
      do m=1,mhlp
        jphy(m)=ihypos
        call idafile(luloop,2,info,1,idisk_lp)
        ndim=info(1)
        call idafile(luloop,2,ihyl(ihypos+1:ihypos+ndim),ndim,idisk_lp)
        call idafile(luloop,2,ihy(ihypos+1:ihypos+ndim),ndim,idisk_lp)
        ihy(ihypos)=ndim
        ihypos=ihypos+ndim+1
      enddo
100   if(linelp.ne.1.and.linelp.ne.13) return
c      print*, "in read_lp, read coe",idisk_lp,norb_inn-norb_dz
!coe_read
      lenr=norb_inn-norb_dz
      do m=1,mhlp
        call idafile(luloop,2,lpnew_coe(norb_dz+1:norb_inn,m),lenr,
     &               idisk_lp)
      enddo
      return
      end

      subroutine value_sort_ploop(mh,logic_ar,logic_w0,logic_w1)
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "stdalloc.fh"
      real*8 w0,w1,s_w0,s_w1
      logical logic_ar,logic_w0,logic_w1
      allocatable ntype(:)
      call mma_allocate(ntype,nvaltype,label='nvaltype')
      zero=0.d0
      nvalue(1:nvaltype)=0
      ntype(1:mh)=0
      mtype=0
      if(.not.logic_w0) then
        do n=1,mh
          vplp_w0(n)=zero
        enddo
      endif
      if(.not.logic_w1) then
        do n=1,mh
          vplp_w1(n)=zero
        enddo
      endif
      do n=1,mh
        if(ntype(n).ne.0) cycle
        s_w0=vplp_w0(n)
        s_w1=vplp_w1(n)
        mtype=mtype+1
        if(mtype.gt.nvaltype) then
          write(6,*) " out of array boundary"
          write(6,*) " subroutine value_sort_ploop"
          write(6,*) " program stop"
#ifdef MOLPRO
#else
      call qtrace
      call abend()
#endif
#ifdef _XIANEST_
#endif
        endif
        ntype(n)=mtype
        nvalue(mtype)=nvalue(mtype)+1
        vplpnew_w0(mtype)=vplp_w0(n)
        vplpnew_w1(mtype)=vplp_w1(n)
        do m=n+1,mh
          if(ntype(m).ne.0) cycle
          w0=vplp_w0(m)
          w1=vplp_w1(m)
          if(w0.ne.s_w0.or.w1.ne.s_w1) cycle
          ntype(m)=mtype
          nvalue(mtype)=nvalue(mtype)+1
        enddo
      enddo

      nstaval(1)=0
      do mty=1,mtype-1
        nstaval(mty+1)=nstaval(mty)+nvalue(mty)
      enddo

      mnew=0
      do mty=1,mtype
        do n=1,mh
          nty=ntype(n)
          if(nty.eq.mty) then
            mnew=mnew+1
            lpnew_head(mnew)=lp_head(n)
            lpnew_lwei(mnew)=lp_lwei(n)
            lpnew_rwei(mnew)=lp_rwei(n)
            if(logic_ar) then
              do lr=norb_dz+1,norb_inn
                lpnew_coe(lr,mnew)=lp_coe(lr,n)
              enddo
            endif
          endif
        enddo
      enddo
      call mma_deallocate(ntype)
      return
      end

      subroutine ext_head_in_act()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      logic_dh=.false.

      jml=mul_tab(jml,ns_sm)
      jmr=mul_tab(jmr,ns_sm)
      goto( 10, 10, 10, 10, 10, 10, 10, 10, 10,110, 10, 10, 10,
     :      10, 10, 10,117, 10, 10, 10, 10, 10,123, 10, 10,126),ipaety
110   call sv_ext_head_in_act()
       goto 10
117   call tv_ext_head_in_act()
       goto 10
123   call dv_ext_head_in_act()
       goto 10
126   call vd_ext_head_in_act()
10    jml=mul_tab(jml,ns_sm)
      jmr=mul_tab(jmr,ns_sm)
       return
      end

      subroutine ext_head_in_dbl()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      logic_dh=.true.
      jml=mul_tab(jml,ns_sm)
      jmr=mul_tab(jmr,ns_sm)
      goto( 10, 10, 10, 10, 10, 10, 10, 10, 10,110, 10, 10, 10,
     :      10, 10, 10,117, 10, 10, 10, 10, 10,123, 10, 10,126),ipaety
110   call sv_ext_head_in_dbl()
      goto 10
117   call tv_ext_head_in_dbl()
      goto 10
123   call dv_ext_head_in_dbl()
      goto 10
126   call vd_ext_head_in_dbl()

10    jml=mul_tab(jml,ns_sm)
      jmr=mul_tab(jmr,ns_sm)
      return
      end

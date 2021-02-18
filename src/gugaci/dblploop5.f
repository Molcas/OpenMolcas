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
      subroutine ds_ardlr_act_c1(lin )
!=======================================================================
!ds(7-1) ar(23)-drl(30)-
!=======================================================================
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      do lri=norb_frz+1,norb_dz
        do lrd=norb_frz+1,lri-1
          lmd=lsm_inn(lrd)
          if(lmd.ne.jml) cycle
          iwdr=just(lri,lri)
          iwdl=jud(lrd)
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo
          w0ds1 =w0_ds(1)
          ni=mod(norb_dz-lrd,2)
          if(ni.eq.0) w0ds1 =-w0ds1
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0ds1
            vplp_w1(mpl)=0.d0
          enddo
          call ar_drl_ext_al_new(lin,lrd,lri )
        enddo
      enddo
      return
      end

      subroutine ds_arblbr_act_c1(lin )
!=======================================================================
!ds(7-3) ar(23)-bl(32)-br(31)-
!=======================================================================
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      isma=mul_tab(iml,imr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jmr) cycle
!ds(7-3) ar(23)-bl(32)-br(31)-
          do lrd=norb_frz+1,lri-1
            lmd=lsm_inn(lrd)
            if(lmd.ne.jml) cycle
            ijk=lrd-norb_frz+ngw2(lri-norb_frz)+ngw3(lrj-norb_frz)
            intpos=intind_ijka(ijk)
            w0ds =w0_ds(3)
            w1ds =w1_ds(3)
            ni=mod(norb_dz-lrj+lri-lrd,2)
            if(ni.eq.0) w0ds =-w0ds
            if(ni.eq.0) w1ds =-w1ds

            iwdr=just(lri,lrj)
            iwdl=jud(lrd)
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0ds
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1ds
            enddo
           call ar_bl_br_ext_al_new(lin,intpos,isma,1 )
          enddo
        enddo
      enddo
      return
      end

      subroutine dt_arblbr_act_c1(lin )
!=======================================================================
!dt(14) ar(23)-bl(32)-br(32)-
!=======================================================================
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      isma=mul_tab(iml,imr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jmr) cycle
!dt(14) ar(23)-bl(32)-br(32)-
          do lrd=norb_frz+1,lri-1
            lmd=lsm_inn(lrd)
            if(lmd.ne.jml) cycle
            ijk=lrd-norb_frz+ngw2(lri-norb_frz)+ngw3(lrj-norb_frz)
            intpos=intind_ijka(ijk)
            w0dt =w0_dt
            w1dt =w1_dt
            ni=mod(norb_dz-lrj+lri-lrd,2)
            if(ni.eq.0) w0dt =-w0dt
            if(ni.eq.0) w1dt =-w1dt

            iwdr=just(lri,lrj)     !
            iwdl=jud(lrd)
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0dt
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1dt
            enddo
           call ar_bl_br_ext_al_new(lin,intpos,isma,1 )
          enddo
        enddo
      enddo
      return
      end

      subroutine dv_ar_act_blbr(lin,jk )
!=======================================================================
!dv(23-1) ar(23)-
!=======================================================================
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      isma=mul_tab(iml,imr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jml) cycle
        ijk=lri-norb_frz+jk
        intpos=intind_ijka(ijk)

        w0dv1=w0_dv(1)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1) w0dv1=-w0dv1
        iwdl=jud(lri)
        iwdr=0
        do mpl=1,mtype
          vplp_w0(mpl)=vplpnew_w0(mpl)*w0dv1
          vplp_w1(mpl)=vplpnew_w1(mpl)*w0dv1
        enddo
        do mpl=1,mhlp
          iwal=lpnew_lwei(mpl)
          iwar=lpnew_rwei(mpl)
          lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
        enddo
        call ar_bl_br_ext_al_new(lin,intpos,isma,1 )
      enddo
      return
      end

      subroutine dv_ar_act_brbr(lin,lrai,lraj )
!=======================================================================
!dv(23-1) ar(23)-
!=======================================================================
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      isma=mul_tab(iml,imr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jml) cycle
        ijk=lri-norb_frz+ngw2(lrai-norb_frz)+ngw3(lraj-norb_frz)
        intpos=intind_ijka(ijk)

        w0dv1=w0_dv(1)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1) w0dv1=-w0dv1
        iwdl=jud(lri)
        iwdr=0
        do mpl=1,mtype
          vplp_w0(mpl)=vplpnew_w0(mpl)*w0dv1
          vplp_w1(mpl)=vplpnew_w1(mpl)*w0dv1
        enddo
        do mpl=1,mhlp
          iwal=lpnew_lwei(mpl)
          iwar=lpnew_rwei(mpl)
          lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
        enddo
        call ar_br_br_ext_ar_new(lin,intpos,isma )
      enddo
      return
      end

      subroutine dv_ar_act_dlr(lin,lra)
!=======================================================================
!dv(23-1) ar(23)-
!=======================================================================
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jml) cycle

        w0dv1=w0_dv(1)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1) w0dv1=-w0dv1
        iwdl=jud(lri)
        iwdr=0
        do mpl=1,mtype
          vplp_w0(mpl)=vplpnew_w0(mpl)*w0dv1
          vplp_w1(mpl)=vplpnew_w1(mpl)*w0dv1
        enddo
        do mpl=1,mhlp
          iwal=lpnew_lwei(mpl)
          iwar=lpnew_rwei(mpl)
          lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
        enddo
        call ar_drl_ext_al_new(lin,lri,lra)
      enddo
      return
      end

      subroutine tt_arbl_act_bl(lin,lra)
!=======================================================================
!tt(11-1) (22)ar(23)-bl(32)-
!tt(11-1) ar(23)-c'(22)-bl(32)-
!tt(11-1) ar(23)-bl(32)-c"(22)-
!=======================================================================
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      jmlr=mul_tab(jml,jmr)
      isma=mul_tab(iml,imr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jmlr) cycle
          ijk=lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz)
          intpos=intind_ijka(ijk)
          call tt1_ext(lri,lrj,nk,1)
          call ar_bl_bl_ext_ar_new(lin,intpos,isma,nk )
          call tt1_ext(lri,lrj,nk,-1)
          call ar_bl_bl_ext_ar_new(lin,intpos,isma,nk )
        enddo
      enddo
      return
      end

      subroutine tt_drl_act_bl(lin,lra)
!tt(11-2) (22)drl(22)-
!tt(11-2) drl(22)-c"(22)-
!tt(11-3) (22)(22)drl(33)-
!tt(11-3) (22)drl(33)-c"(22)-
!tt(11-3) drl(33)-c"(22)-c"(22)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

        do lri=norb_frz+1,norb_dz
          lmi=lsm_inn(lri)
          do lrj=lri+1,norb_dz
            lmj=lsm_inn(lrj)
            lmij=mul_tab(lmi,lmj)
            if(lmij.ne.jml) cycle
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0_tt(2)
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1_tt(2)
            enddo
!tt(11-2) (22)drl(22)-
!tt(11-2) drl(22)-c"(22)-
            iwdl=just(lri,lrj)      !
            iwdr=iwdl
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            call drl_bl_ext_ar_new(lin,lri,lra)
            call drl_bl_ext_ar_new(lin,lrj,lra)
!tt(11-3) drl(33)-c"(22)-c"(22)-
!tt(11-3) (22)drl(33)-c"(22)-
!tt(11-3) (22)(22)drl(33)-
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0_tt(3)
              vplp_w1(mpl)=0.d0
            enddo
            if(lra.gt.norb_dz) then
              call drl_bl_sum_ar_new(lin,lri,lrj,lra)
            else
              do lrk=1,norb_dz
                if(lrk.eq.lri) cycle
                if(lrk.eq.lrj) cycle
                call drl_bl_ext_ar_new(lin,lrk,lra)
              enddo
            endif
c            do lrk=1,norb_dz
c              if(lrk.eq.lri) cycle
c              if(lrk.eq.lrj) cycle
c           call drl_bl_ext_ar_new(lin,lrk,lra)
c            enddo
          enddo
        enddo
      return
      end

      subroutine tt_arbl_act_br(lin,lra)
!=======================================================================
!tt(11-1) (22)ar(23)-bl(32)-
!tt(11-1) ar(23)-c'(22)-bl(32)-
!tt(11-1) ar(23)-bl(32)-c"(22)-
!=======================================================================
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      jmlr=mul_tab(jml,jmr)
      isma=mul_tab(iml,imr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jmlr) cycle
          ijk=lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz)
          intpos=intind_ijka(ijk)
          call tt1_ext(lri,lrj,nk,1)
          call ar_bl_br_ext_al_new(lin,intpos,isma,nk)
          call tt1_ext(lri,lrj,nk,-1)
          call ar_bl_br_ext_al_new(lin,intpos,isma,nk)
        enddo
      enddo
      return
      end

      subroutine tt_drl_act_br(lin,lra)
!tt(11-2) (22)drl(22)-
!tt(11-2) drl(22)-c"(22)-
!tt(11-3) (22)(22)drl(33)-
!tt(11-3) (22)drl(33)-c"(22)-
!tt(11-3) drl(33)-c"(22)-c"(22)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

        do lri=norb_frz+1,norb_dz
          lmi=lsm_inn(lri)
          do lrj=lri+1,norb_dz
            lmj=lsm_inn(lrj)
            lmij=mul_tab(lmi,lmj)
            if(lmij.ne.jml) cycle
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0_tt(2)
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1_tt(2)
            enddo
!tt(11-2) (22)drl(22)-
!tt(11-2) drl(22)-c"(22)-
            iwdl=just(lri,lrj)       !
            iwdr=iwdl
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            call drl_br_ext_al_new(lin,lri,lra)
            call drl_br_ext_al_new(lin,lrj,lra)
!tt(11-3) drl(33)-c"(22)-c"(22)-
!tt(11-3) (22)drl(33)-c"(22)-
!tt(11-3) (22)(22)drl(33)-
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0_tt(3)
              vplp_w1(mpl)=0.d0
            enddo
            if(lra.gt.norb_dz) then
             call drl_br_sum_al_new(lin,lri,lrj,lra)
            else
             do lrk=1,norb_dz
                if(lrk.eq.lri) cycle
                if(lrk.eq.lrj) cycle
                call drl_br_ext_al_new(lin,lrk,lra)
              enddo
            endif
c            do lrk=1,norb_dz
c              if(lrk.eq.lri) cycle
c              if(lrk.eq.lrj) cycle
c           call drl_br_ext_al_new(lin,lrk,lra)
c            enddo
          enddo
        enddo
      return
      end

      subroutine st_drl_act_bl(lin,lra)
!st(2-5) (22)drl(12)-
!st(2-6) drl(22)-c"(12)-
!st(2-7) drl(12)-c"(22)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
      do lrj=lri+1,norb_dz
        lmj=lsm_inn(lrj)
        lmij=mul_tab(lmi,lmj)
        if(lmij.ne.jml) cycle
!-------------------------------------------------------------------
!st(2-5) (22)drl(12)-
          iwdl=just(lri,lrj)
          iwdr=iwdl         !
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo
          do mpl=1,mtype
            vplp_w0(mpl)=0.d0
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1_st(5)
          enddo
          call drl_bl_ext_ar_new(lin,lrj,lra)
!st(2-6) drl(22)-c"(12)-
          do mpl=1,mtype
            vplp_w0(mpl)=0.d0
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1_st(6)
          enddo
          call drl_bl_ext_ar_new(lin,lri,lra)
        enddo
      enddo
      return
      end

      subroutine st_arbl_act_bl(lin,lra)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      isma=mul_tab(iml,imr)
      do lri=norb_frz+1,norb_dz-1
        do lrj=lri+1,norb_dz
          ijk=lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz)
          intpos=intind_ijka(ijk)
          call st1_ext(lri,lrj,nk)
          if(nk.ne.0) call ar_bl_bl_ext_ar_new(lin,intpos,isma,1)
          call st2_ext(lri,lrj,nk)
          if(nk.ne.0) call ar_bl_bl_ext_ar_new(lin,intpos,isma,nk)
          call st4_ext(lri,lrj,nk,1)
          if(nk.ne.0) call ar_bl_bl_ext_ar_new(lin,intpos,isma,nk)
          call st4_ext(lri,lrj,nk,-1)
          if(nk.ne.0) call ar_bl_bl_ext_ar_new(lin,intpos,isma,nk)
        enddo
      enddo
      return
      end

      subroutine st_drl_act_br(lin,lra)
!st(2-5) (22)drl(12)-
!st(2-6) drl(22)-c"(12)-
!st(2-7) drl(12)-c"(22)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
      do lrj=lri+1,norb_dz
        lmj=lsm_inn(lrj)
        lmij=mul_tab(lmi,lmj)
        if(lmij.ne.jml) cycle
!-------------------------------------------------------------------
!st(2-5) (22)drl(12)-
          iwdl=just(lri,lrj)
          iwdr=iwdl             !
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo
          do mpl=1,mtype
            vplp_w0(mpl)=0.d0
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1_st(5)
          enddo
          call drl_br_ext_al_new(lin,lrj,lra)
!st(2-6) drl(22)-c"(12)-
          do mpl=1,mtype
            vplp_w0(mpl)=0.d0
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1_st(6)
          enddo
          call drl_br_ext_al_new(lin,lri,lra)
        enddo
      enddo
      return
      end

      subroutine st_arbl_act_br(lin,lra)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      isma=mul_tab(iml,imr)
      do lri=norb_frz+1,norb_dz-1
        do lrj=lri+1,norb_dz
          ijk=lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz)
          intpos=intind_ijka(ijk)
          call st1_ext(lri,lrj,nk)
          if(nk.ne.0) call ar_bl_br_ext_al_new(lin,intpos,isma,1)
          call st2_ext(lri,lrj,nk)
          if(nk.ne.0) call ar_bl_br_ext_al_new(lin,intpos,isma,nk)
          call st4_ext(lri,lrj,nk,1)
          if(nk.ne.0) call ar_bl_br_ext_al_new(lin,intpos,isma,nk)
          call st4_ext(lri,lrj,nk,-1)
          if(nk.ne.0) call ar_bl_br_ext_al_new(lin,intpos,isma,nk)
        enddo
      enddo
      return
      end

      subroutine ss_arbl_act_br(lin,lra)
!=======================================================================
!ss(1-3)  ar(13)-bl(20)-
!ss(1-6)  (11)-ar(23)-bl(32)-
!ss(1-7)  ar(13)-c'(21)-bl(32)-
!ss(1-8)  ar(13)-c'(22)-bl(31)-
!ss(1-9)  ar(23)-c'(11)-bl(32)-
!ss(1-11) ar(13)-bl(31)-c"(22)-
!ss(1-12) ar(13)-bl(32)-c"(21)-
!ss(1-13) ar(23)-bl(31)-c"(12)-
!=======================================================================
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      isma=mul_tab(iml,imr)
      do lri=norb_frz+1,norb_dz
        do lrj=lri+1,norb_dz
          ijk=lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz)
          intpos=intind_ijka(ijk)
!-------------------------------------------------------------------
          call ss2_ext(lri,lrj,nk)
          if(nk.ne.0) call ar_bl_br_ext_al_new(lin,intpos,isma,1)
          call ss4_ext(lri,lrj,nk)
          if(nk.ne.0) call ar_bl_br_ext_al_new(lin,intpos,isma,1)
          call ss5_ext(lri,lrj,nk)
          if(nk.ne.0) call ar_bl_br_ext_al_new(lin,intpos,isma,nk)
          call ss10_ext(lri,lrj,nk)
          if(nk.ne.0) call ar_bl_br_ext_al_new(lin,intpos,isma,nk)
          call ss14_ext(lri,lrj,nk)
          if(nk.ne.0) call ar_bl_br_ext_al_new(lin,intpos,isma,nk)
!-------------------------------------------------------------------
        enddo
      enddo
      return
      end

      subroutine ss_drl_act_br(lin,lra)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

        do lri=norb_frz+1,norb_dz
          lmi=lsm_inn(lri)
          if(jml.eq.1) then
!ss(1-20) drl(33)-c"(00)-
            iwdl=just(lri,lri)
            iwdr=iwdl
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0_ss(20)
              vplp_w1(mpl)=0.d0
            enddo
            if(lri.gt.norb_dz) then
             call drl_br_sum_al_new(lin,lri,0,lri)
            else
             do lrk=1,norb_dz
                if(lrk.eq.lri) cycle
                call drl_br_ext_al_new(lin,lrk,lra)
              enddo
            endif
         endif
          do lrj=lri+1,norb_dz
            lmj=lsm_inn(lrj)
            lmij=mul_tab(lmi,lmj)
            if(lmij.ne.jml) cycle
!ss(1-15) (22)-drl(11)-
            iwdl=just(lri,lrj)
            iwdr=iwdl
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0_ss(15)
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1_ss(15)
            enddo
            call drl_br_ext_al_new(lin,lrj,lra)
!ss(1-17) drl(22)-c"(11)-
            iwdl=just(lri,lrj)
            iwdr=iwdl
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0_ss(17)
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1_ss(17)
            enddo
            call drl_br_ext_al_new(lin,lri,lra)
!ss(1-20) drl(33)-c"(22)-c"(11)-
!ss(1-20) (22)drl(33)-c"(11)-
!ss(1-20) (22)(11)drl(33)-
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0_ss(20)
              vplp_w1(mpl)=0.d0
            enddo
            if(lra.gt.norb_dz) then
              call drl_br_sum_al_new(lin,lri,lrj,lra)
            else
              do lrk=1,norb_dz
                if(lrk.eq.lri) cycle
                if(lrk.eq.lrj) cycle
                call drl_br_ext_al_new(lin,lrk,lra)
              enddo
            endif
         enddo
        enddo
      return
      end

      subroutine ts_arbl_act_br(lin,lra)
!=======================================================================
!ts(3) a&r-b^l-  act -b&l ............................................
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      isma=mul_tab(iml,imr)
        do lri=norb_frz+1,norb_dz
          do lrj=lri+1,norb_dz
            ijk=lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz)
            intpos=intind_ijka(ijk)
            call ts1_ext(lri,lrj,nk)
            if(nk.ne.0) call ar_bl_br_ext_al_new(lin,intpos,isma,1)
            call ts2_ext(lri,lrj,nk,1)
            if(nk.ne.0) call ar_bl_br_ext_al_new(lin,intpos,isma,nk)
            call ts2_ext(lri,lrj,nk,-1)
            if(nk.ne.0) call ar_bl_br_ext_al_new(lin,intpos,isma,nk)
            call ts4_ext(lri,lrj,nk)
            if(nk.ne.0) call ar_bl_br_ext_al_new(lin,intpos,isma,nk)
          enddo
        enddo
      return
      end

      subroutine ts_arbl_act_bl(lin,lra)
!=======================================================================
!ts(3) a&r-b^l-  act -b&l ............................................
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      isma=mul_tab(iml,imr)
        do lri=norb_frz+1,norb_dz
          do lrj=lri+1,norb_dz
            ijk=lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz)
            intpos=intind_ijka(ijk)
            call ts1_ext(lri,lrj,nk)
            if(nk.ne.0) call ar_bl_bl_ext_ar_new(lin,intpos,isma,1)
            call ts2_ext(lri,lrj,nk,1)
            if(nk.ne.0) call ar_bl_bl_ext_ar_new(lin,intpos,isma,nk)
            call ts2_ext(lri,lrj,nk,-1)
            if(nk.ne.0) call ar_bl_bl_ext_ar_new(lin,intpos,isma,nk)
            call ts4_ext(lri,lrj,nk)
            if(nk.ne.0) call ar_bl_bl_ext_ar_new(lin,intpos,isma,nk)
          enddo
        enddo
      return
      end

      subroutine ss_drl_act_bl(lin,lra)
!=======================================================================
!ss(1)    act -bl-
!ss(1-16) (11)-drl(22)-
!ss(1-17) drl(22)-c"(11)-
!ss(1-18) drl(11)-c"(22)-
!ss(1-19) drl(12)-c"(21)-
!ss(1-20) drl(33)-c"(11)-c"(22)-
!ss(1-20) (11)drl(33)-c"(22)-
!ss(1-20) (11)(22)drl(33)-
!=======================================================================
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

        do lri=norb_frz+1,norb_dz
          lmi=lsm_inn(lri)
          if(jml.eq.1) then
!ss(1-20) drl(33)-c"(00)-
            iwdl=just(lri,lri)
            iwdr=iwdl
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0_ss(20)
              vplp_w1(mpl)=0.d0
            enddo
            if(lra.gt.norb_dz) then
              call drl_bl_sum_ar_new(lin,lri,0,lra)
            else
              do lrk=1,norb_dz
                if(lrk.eq.lri) cycle
                call drl_bl_ext_ar_new(lin,lrk,lra)
              enddo
            endif
c            do lrk=1,norb_dz
c              if(lrk.eq.lri) cycle
c              call drl_bl_ext_ar_new(lin,lrk,lra)
c            enddo
          endif
          do lrj=lri+1,norb_dz
            lmj=lsm_inn(lrj)
            lmij=mul_tab(lmi,lmj)
            if(lmij.ne.jml) cycle
            iwdl=just(lri,lrj)
            iwdr=iwdl
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
!ss(1-15) (22)-drl(11)-
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0_ss(15)
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1_ss(15)
            enddo
            call drl_bl_ext_ar_new(lin,lrj,lra)
!ss(1-17) drl(22)-c"(11)-
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0_ss(17)
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1_ss(17)
            enddo
            call drl_bl_ext_ar_new(lin,lri,lra)
!ss(1-20) drl(33)-c"(22)-c"(11)-
!ss(1-20) (22)drl(33)-c"(11)-
!ss(1-20) (22)(11)-drl(33)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0_ss(20)
              vplp_w1(mpl)=0.d0
            enddo
            if(lra.gt.norb_dz) then
              call drl_bl_sum_ar_new(lin,lri,lrj,lra)
            else
              do lrk=1,norb_dz
                if(lrk.eq.lri) cycle
                if(lrk.eq.lrj) cycle
                call drl_bl_ext_ar_new(lin,lrk,lra)
              enddo
            endif
c            do lrk=1,norb_dz
c              if(lrk.eq.lri) cycle
c              if(lrk.eq.lrj) cycle
c              call drl_bl_ext_ar_new(lin,lrk,lra)
c            enddo
          enddo
        enddo
      return
      end

      subroutine ss_arbl_act_bl(lin,lra)
!=======================================================================
!ss(1)    act -bl-
!ss(1-3)  ar(13)-bl(20)-
!ss(1-6)  (11)-ar(23)-bl(32)-
!ss(1-7)  ar(13)-c'(21)-bl(32)-
!ss(1-8)  ar(13)-c'(22)-bl(31)-
!ss(1-9)  ar(23)-c'(11)-bl(32)-
!ss(1-11) ar(13)-bl(31)-c"(22)-
!ss(1-12) ar(13)-bl(32)-c"(21)-
!ss(1-13) ar(23)-bl(31)-c"(12)-
!=======================================================================
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      isma=mul_tab(iml,imr)
      do lri=norb_frz+1,norb_dz
        do lrj=lri+1,norb_dz
          ijk=lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz)
          intpos=intind_ijka(ijk)
          call ss2_ext(lri,lrj,nk)
          if(nk.ne.0) call ar_bl_bl_ext_ar_new(lin,intpos,isma,1)
          call ss4_ext(lri,lrj,nk)
          if(nk.ne.0) call ar_bl_bl_ext_ar_new(lin,intpos,isma,1)
          call ss5_ext(lri,lrj,nk)
          if(nk.ne.0) call ar_bl_bl_ext_ar_new(lin,intpos,isma,nk)
          call ss10_ext(lri,lrj,nk)
          if(nk.ne.0) call ar_bl_bl_ext_ar_new(lin,intpos,isma,nk)
          call ss14_ext(lri,lrj,nk)
          if(nk.ne.0) call ar_bl_bl_ext_ar_new(lin,intpos,isma,nk)
        enddo
      enddo
      return
      end

      subroutine sd_ar_act_bl(lin,lra)
!sd(6-1) a&r(02)-
!sd(6-2) (22)a&(13)-
!sd(6-3) a&r(13)c'(22)-
!sd(6-4) a&r(23)c'(12)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      jmlr=mul_tab(jml,jmr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jmlr) cycle
        w0sd1 =w0_sd(1)
        w0sd2 =w0_sd(2)
        w0sd4 =w0_sd(4)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1) then
          w0sd1 =-w0sd1
          w0sd2 =-w0sd2
          w0sd4 =-w0sd4
        endif
        if(jml.eq.1.and.lmi.eq.jmr) then
!sd(6-1) a&r(02)-
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd1
            vplp_w1(mpl)=vplpnew_w1(mpl)*w0sd1
          enddo
          iwdl=just(lri,lri)
          iwdr=jud(lri)
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo
          if(lin.eq.1) call ar_bl_ext_ss(lri,lra,1 )
          if(lin.eq.2) call ar_bl_ext_st(lri,lra,1 )
          if(lin.eq.3) call ar_bl_ext_ts(lri,lra,1 )
          if(lin.eq.11) call ar_bl_ext_tt(lri,lra,1 )
        endif
!sd(6-2) (22)a&(13)-
        do mpl=1,mtype
          vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd2
          vplp_w1(mpl)=vplpnew_w1(mpl)*w0sd2
        enddo
        do lrk=norb_frz+1,lri-1
          lmk=lsm_inn(lrk)
          if(lmk.ne.jmr) cycle
          iwdl=just(lrk,lri)
          iwdr=jud(lrk)
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo
          if(lin.eq.1) call ar_bl_ext_ss(lri,lra,1 )
          if(lin.eq.2) call ar_bl_ext_st(lri,lra,1 )
          if(lin.eq.3) call ar_bl_ext_ts(lri,lra,1 )
          if(lin.eq.11) call ar_bl_ext_tt(lri,lra,1 )
        enddo
!sd(6-4) a&r(23)-c'(12)-
        do mpl=1,mtype
          vplp_w0(mpl)=-vplpnew_w0(mpl)*w0sd4
          vplp_w1(mpl)=-vplpnew_w1(mpl)*w0sd4
        enddo
        do lrk=lri+1,norb_dz
          lmk=lsm_inn(lrk)
          if(lmk.ne.jmr) cycle
          iwdl=just(lri,lrk)
          iwdr=jud(lrk)
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo
          if(lin.eq.1) call ar_bl_ext_ss(lri,lra,1 )
          if(lin.eq.2) call ar_bl_ext_st(lri,lra,1 )
          if(lin.eq.3) call ar_bl_ext_ts(lri,lra,1 )
          if(lin.eq.11) call ar_bl_ext_tt(lri,lra,1 )
        enddo
      enddo
      return
      end

      subroutine sd_ar_act_blbl(lin,jk)
!sd(6-1) a&r(02)-
!sd(6-2) (22)a&(13)-
!sd(6-3) a&r(13)c'(22)-
!sd(6-4) a&r(23)c'(12)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      jmlr=mul_tab(jml,jmr)
      isma=mul_tab(iml,imr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jmlr) cycle
        ijk=lri-norb_frz+jk
        intpos=intind_ijka(ijk)
        w0sd1 =w0_sd(1)
        w0sd2 =w0_sd(2)
        w0sd4 =w0_sd(4)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1) then
          w0sd1 =-w0sd1
          w0sd2 =-w0sd2
          w0sd4 =-w0sd4
        endif
!sd(6-1) a&r(02)-
        if(jml.eq.1.and.lmi.eq.jmr) then
          iwdl=just(lri,lri)
          iwdr=jud(lri)
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd1
            vplp_w1(mpl)=vplpnew_w1(mpl)*w0sd1
          enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
           call ar_bl_bl_ext_ar_new(lin,intpos,isma,1 )
          endif
!------------------------------------------------------------------
!sd(6-2) c(22)a&(13)-
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd2
              vplp_w1(mpl)=vplpnew_w1(mpl)*w0sd2
            enddo
            do lrk=norb_frz+1,lri-1
              lmk=lsm_inn(lrk)
              if(lmk.ne.jmr) cycle
              iwdl=just(lrk,lri)
              iwdr=jud(lrk)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
              call ar_bl_bl_ext_ar_new(lin,intpos,isma,1)
            enddo
!-------------------------------------------------------------------
!sd(6-4) a&r(23)c'(12)-
            do mpl=1,mtype
              vplp_w0(mpl)=-vplpnew_w0(mpl)*w0sd4
              vplp_w1(mpl)=-vplpnew_w1(mpl)*w0sd4
            enddo
            do lrk=lri+1,norb_dz
              lmk=lsm_inn(lrk)
              if(lmk.ne.jmr) cycle
              iwdl=just(lri,lrk)
              iwdr=jud(lrk)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
              call ar_bl_bl_ext_ar_new(lin,intpos,isma,1)
            enddo
          enddo
        return
      end

      subroutine sd_ar_act_brbr(lin,jk)
!sd(6-1) a&r(02)-
!sd(6-2) (22)a&(13)-
!sd(6-3) a&r(13)c'(22)-
!sd(6-4) a&r(23)c'(12)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      jmlr=mul_tab(jml,jmr)
      isma=mul_tab(iml,imr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jmlr) cycle
        ijk=lri-norb_frz+jk
        intpos=intind_ijka(ijk)
        w0sd1 =w0_sd(1)
        w0sd2 =w0_sd(2)
        w0sd4 =w0_sd(4)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1) then
          w0sd1 =-w0sd1
          w0sd2 =-w0sd2
          w0sd4 =-w0sd4
        endif
!sd(6-1) a&r(02)-
        if(jml.eq.1.and.lmi.eq.jmr) then
          iwdl=just(lri,lri)
          iwdr=jud(lri)
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd1
            vplp_w1(mpl)=vplpnew_w1(mpl)*w0sd1
          enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
         call ar_br_br_ext_ar_new(lin,intpos,isma )
          endif
!------------------------------------------------------------------
!sd(6-2) c(22)a&(13)-
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd2
              vplp_w1(mpl)=vplpnew_w1(mpl)*w0sd2
            enddo
            do lrk=norb_frz+1,lri-1
              lmk=lsm_inn(lrk)
              if(lmk.ne.jmr) cycle
              iwdl=just(lrk,lri)
              iwdr=jud(lrk)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
              call ar_br_br_ext_ar_new(lin,intpos,isma )
            enddo
!-------------------------------------------------------------------
!sd(6-4) a&r(23)c'(12)-
            do mpl=1,mtype
              vplp_w0(mpl)=-vplpnew_w0(mpl)*w0sd4
              vplp_w1(mpl)=-vplpnew_w1(mpl)*w0sd4
            enddo
            do lrk=lri+1,norb_dz
              lmk=lsm_inn(lrk)
              if(lmk.ne.jmr) cycle
              iwdl=just(lri,lrk)
              iwdr=jud(lrk)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
              call ar_br_br_ext_ar_new(lin,intpos,isma )
            enddo
          enddo
        return
      end

      subroutine sd_ar_act_blbr(lin,jk)
!sd(6-1) a&r(02)-
!sd(6-2) (22)a&(13)-
!sd(6-3) a&r(13)c'(22)-
!sd(6-4) a&r(23)c'(12)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      jmlr=mul_tab(jml,jmr)
      isma=mul_tab(iml,imr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jmlr) cycle
        ijk=lri-norb_frz+jk
        intpos=intind_ijka(ijk)
        w0sd1 =w0_sd(1)
        w0sd2 =w0_sd(2)
        w0sd4 =w0_sd(4)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1) then
          w0sd1 =-w0sd1
          w0sd2 =-w0sd2
          w0sd4 =-w0sd4
        endif
!sd(6-1) a&r(02)-
        if(jml.eq.1.and.lmi.eq.jmr) then
          iwdl=just(lri,lri)
          iwdr=jud(lri)
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd1
            vplp_w1(mpl)=vplpnew_w1(mpl)*w0sd1
          enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
           call ar_bl_br_ext_al_new(lin,intpos,isma,1 )
          endif
!------------------------------------------------------------------
!sd(6-2) c(22)a&(13)-
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd2
              vplp_w1(mpl)=vplpnew_w1(mpl)*w0sd2
            enddo
            do lrk=norb_frz+1,lri-1
              lmk=lsm_inn(lrk)
              if(lmk.ne.jmr) cycle
              iwdl=just(lrk,lri)
              iwdr=jud(lrk)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
              call ar_bl_br_ext_al_new(lin,intpos,isma,1 )
            enddo
!-------------------------------------------------------------------
!sd(6-4) a&r(23)c'(12)-
            do mpl=1,mtype
              vplp_w0(mpl)=-vplpnew_w0(mpl)*w0sd4
              vplp_w1(mpl)=-vplpnew_w1(mpl)*w0sd4
            enddo
            do lrk=lri+1,norb_dz
              lmk=lsm_inn(lrk)
              if(lmk.ne.jmr) cycle
              iwdl=just(lri,lrk)
              iwdr=jud(lrk)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
              call ar_bl_br_ext_al_new(lin,intpos,isma,1 )
            enddo
          enddo
        return
      end

      subroutine sd_ar_act_dlr(lin,lra)
!sd(6-1) a&r(02)-
!sd(6-2) (22)a&(13)-
!sd(6-3) a&r(13)c'(22)-
!sd(6-4) a&r(23)c'(12)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      jmlr=mul_tab(jml,jmr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jmlr) cycle
        w0sd1 =w0_sd(1)
        w0sd2 =w0_sd(2)
        w0sd4 =w0_sd(4)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1) then
          w0sd1 =-w0sd1
          w0sd2 =-w0sd2
          w0sd4 =-w0sd4
        endif
!sd(6-1) a&r(02)-
        if(jml.eq.1.and.lmi.eq.jmr) then
          iwdl=just(lri,lri)
          iwdr=jud(lri)
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd1
            vplp_w1(mpl)=vplpnew_w1(mpl)*w0sd1
          enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            call ar_drl_ext_al_new(lin,lri,lra)
          endif
!------------------------------------------------------------------
!sd(6-2) c(22)a&(13)-
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd2
              vplp_w1(mpl)=vplpnew_w1(mpl)*w0sd2
            enddo
            do lrk=norb_frz+1,lri-1
              lmk=lsm_inn(lrk)
              if(lmk.ne.jmr) cycle
              iwdl=just(lrk,lri)
              iwdr=jud(lrk)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
              call ar_drl_ext_al_new(lin,lri,lra)
            enddo
!-------------------------------------------------------------------
!sd(6-4) a&r(23)c'(12)-
            do mpl=1,mtype
              vplp_w0(mpl)=-vplpnew_w0(mpl)*w0sd4
              vplp_w1(mpl)=-vplpnew_w1(mpl)*w0sd4
            enddo
            do lrk=lri+1,norb_dz
              lmk=lsm_inn(lrk)
              if(lmk.ne.jmr) cycle
              iwdl=just(lri,lrk)
              iwdr=jud(lrk)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
              call ar_drl_ext_al_new(lin,lri,lra)
            enddo
          enddo
        return
      end

      subroutine td_ar_act_bl(lin,lra)
!td(13-1) (22)a(23)
!td(13-1) a(23)c'(22)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      jmlr=mul_tab(jml,jmr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jmlr) cycle
        w0td1=w0_td(1)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1) w0td1=-w0td1
!-------------------------------------------------------------------
!td(13-1) (22)a(23)
          do lrk=norb_frz+1,lri-1
            lmk=lsm_inn(lrk)
           if(lmk.ne.jmr) cycle
            iwdl=just(lrk,lri)                   !
            iwdr=jud(lrk)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0td1
              vplp_w1(mpl)=vplpnew_w1(mpl)*w0td1
            enddo
              do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
          if(lin.eq.1) call ar_bl_ext_ss(lri,lra,1 )
          if(lin.eq.2) call ar_bl_ext_st(lri,lra,1 )
          if(lin.eq.3) call ar_bl_ext_ts(lri,lra,1 )
          if(lin.eq.11) call ar_bl_ext_tt(lri,lra,1 )
          enddo
!-------------------------------------------------------------------
!td(13-1) a(23)c'(22)
          do lrk=lri+1,norb_dz
            lmk=lsm_inn(lrk)
           if(lmk.ne.jmr) cycle
            iwdl=just(lri,lrk)        !
            iwdr=jud(lrk)
            do mpl=1,mtype
              vplp_w0(mpl)=-vplpnew_w0(mpl)*w0td1
              vplp_w1(mpl)=-vplpnew_w1(mpl)*w0td1
            enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
          if(lin.eq.1) call ar_bl_ext_ss(lri,lra,1 )
          if(lin.eq.2) call ar_bl_ext_st(lri,lra,1 )
          if(lin.eq.3) call ar_bl_ext_ts(lri,lra,1 )
          if(lin.eq.11) call ar_bl_ext_tt(lri,lra,1 )
          enddo
        enddo
      return
      end

      subroutine td_ar_act_brbr(lin,jk )
!td(13-1) (22)a(23)
!td(13-1) a(23)c'(22)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      jmlr=mul_tab(jml,jmr)
      isma=mul_tab(iml,imr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jmlr) cycle
        w0td1=w0_td(1)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1) w0td1=-w0td1
        ijk=lri-norb_frz+jk
        intpos=intind_ijka(ijk)
!-------------------------------------------------------------------
!td(13-1) (22)a(23)
          do lrk=norb_frz+1,lri-1
            lmk=lsm_inn(lrk)
           if(lmk.ne.jmr) cycle
            iwdl=just(lrk,lri)     !
            iwdr=jud(lrk)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0td1
              vplp_w1(mpl)=vplpnew_w1(mpl)*w0td1
            enddo
              do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
          call ar_br_br_ext_ar_new(lin,intpos,isma )
          enddo
!-------------------------------------------------------------------
!td(13-1) a(23)c'(22)
          do lrk=lri+1,norb_dz
            lmk=lsm_inn(lrk)
           if(lmk.ne.jmr) cycle
            iwdl=just(lri,lrk)     !
            iwdr=jud(lrk)
            do mpl=1,mtype
              vplp_w0(mpl)=-vplpnew_w0(mpl)*w0td1
              vplp_w1(mpl)=-vplpnew_w1(mpl)*w0td1
            enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
          call ar_br_br_ext_ar_new(lin,intpos,isma )
          enddo
        enddo
      return
      end

      subroutine td_ar_act_blbl(lin,jk )
!td(13-1) (22)a(23)
!td(13-1) a(23)c'(22)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      jmlr=mul_tab(jml,jmr)
      isma=mul_tab(iml,imr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jmlr) cycle
        w0td1=w0_td(1)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1) w0td1=-w0td1
        ijk=lri-norb_frz+jk
        intpos=intind_ijka(ijk)
!-------------------------------------------------------------------
!td(13-1) (22)a(23)
          do lrk=norb_frz+1,lri-1
            lmk=lsm_inn(lrk)
           if(lmk.ne.jmr) cycle
            iwdl=just(lrk,lri)     !
            iwdr=jud(lrk)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0td1
              vplp_w1(mpl)=vplpnew_w1(mpl)*w0td1
            enddo
              do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
          call ar_bl_bl_ext_ar_new(lin,intpos,isma,1 )
          enddo
!-------------------------------------------------------------------
!td(13-1) a(23)c'(22)
          do lrk=lri+1,norb_dz
            lmk=lsm_inn(lrk)
           if(lmk.ne.jmr) cycle
            iwdl=just(lri,lrk)    !
            iwdr=jud(lrk)
            do mpl=1,mtype
              vplp_w0(mpl)=-vplpnew_w0(mpl)*w0td1
              vplp_w1(mpl)=-vplpnew_w1(mpl)*w0td1
            enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
          call ar_bl_bl_ext_ar_new(lin,intpos,isma,1 )
          enddo
        enddo
      return
      end

      subroutine td_ar_act_blbr(lin,jk )
!td(13-1) (22)a(23)
!td(13-1) a(23)c'(22)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      jmlr=mul_tab(jml,jmr)
      isma=mul_tab(iml,imr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jmlr) cycle
        w0td1=w0_td(1)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1) w0td1=-w0td1
        ijk=lri-norb_frz+jk
        intpos=intind_ijka(ijk)
!-------------------------------------------------------------------
!td(13-1) (22)a(23)
          do lrk=norb_frz+1,lri-1
            lmk=lsm_inn(lrk)
           if(lmk.ne.jmr) cycle
            iwdl=just(lrk,lri)    !
            iwdr=jud(lrk)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0td1
              vplp_w1(mpl)=vplpnew_w1(mpl)*w0td1
            enddo
              do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
          call ar_bl_br_ext_al_new(lin,intpos,isma,1 )
          enddo
!-------------------------------------------------------------------
!td(13-1) a(23)c'(22)
          do lrk=lri+1,norb_dz
            lmk=lsm_inn(lrk)
           if(lmk.ne.jmr) cycle
            iwdl=just(lri,lrk)    !
            iwdr=jud(lrk)
            do mpl=1,mtype
              vplp_w0(mpl)=-vplpnew_w0(mpl)*w0td1
              vplp_w1(mpl)=-vplpnew_w1(mpl)*w0td1
            enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
          call ar_bl_br_ext_al_new(lin,intpos,isma,1 )
          enddo
        enddo
      return
      end

      subroutine td_ar_act_dlr(lin,lra)
!td(13-1) (22)a(23)
!td(13-1) a(23)c'(22)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      jmlr=mul_tab(jml,jmr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jmlr) cycle
        w0td1=w0_td(1)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1) w0td1=-w0td1
!-------------------------------------------------------------------
!td(13-1) (22)a(23)
          do lrk=norb_frz+1,lri-1
            lmk=lsm_inn(lrk)
           if(lmk.ne.jmr) cycle
            iwdl=just(lrk,lri)   !
            iwdr=jud(lrk)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0td1
              vplp_w1(mpl)=vplpnew_w1(mpl)*w0td1
            enddo
              do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            call ar_drl_ext_al_new(lin,lri,lra)
          enddo
!-------------------------------------------------------------------
!td(13-1) a(23)c'(22)
          do lrk=lri+1,norb_dz
            lmk=lsm_inn(lrk)
           if(lmk.ne.jmr) cycle
            iwdl=just(lri,lrk)     !
            iwdr=jud(lrk)
            do mpl=1,mtype
              vplp_w0(mpl)=-vplpnew_w0(mpl)*w0td1
              vplp_w1(mpl)=-vplpnew_w1(mpl)*w0td1
            enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            call ar_drl_ext_al_new(lin,lri,lra)
          enddo
        enddo
      return
      end

      subroutine dd_drl_act_bl(lin,lra)
!dd(19-2) drl(22)-
!dd(19-3) drl(33)-c"(22)-
!dd(19-3) (22)drl(33)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      do lrl=norb_frz+1,norb_dz
        lml=lsm_inn(lrl)
        if(lml.ne.jml) cycle
!-------------------------------------------------------------------
         iwdl=jud(lrl)
         iwdr=iwdl
         do mpl=1,mhlp
           iwal=lpnew_lwei(mpl)
           iwar=lpnew_rwei(mpl)
           lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
           lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
         enddo
         w0=w0_dd(2)
         w1=w1_dd(2)
         do mpl=1,mtype
           vplp_w0(mpl)=vplpnew_w0(mpl)*w0
           vplp_w1(mpl)=vplpnew_w1(mpl)*w1
         enddo
         call drl_bl_ext_ar_new(lin,lrl,lra)

         w0=w0_dd(3)
         w1=0.d0
         do mpl=1,mtype
           vplp_w0(mpl)=vplpnew_w0(mpl)*w0
           vplp_w1(mpl)=vplpnew_w1(mpl)*w1
         enddo
         if(lra.gt.norb_dz) then
          call drl_bl_sum_ar_new(lin,lrl,0,lra)
         else
          do lrk=1,norb_dz
            if(lrk.eq.lrl) cycle
            call drl_bl_ext_ar_new(lin,lrk,lra)
          enddo
         endif
c        do lrk=1,norb_dz
c           if(lrk.eq.lrl) cycle
c           call drl_bl_ext_ar_new(lin,lrk,lra)
c         enddo
       enddo
      return
      end

      subroutine dd_drl_act_br(lin,lra)
!dd(19-2) drl(22)-
!dd(19-3) drl(33)-c"(22)-
!dd(19-3) (22)drl(33)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      do lrl=norb_frz+1,norb_dz
        lml=lsm_inn(lrl)
        if(lml.ne.jml) cycle
!-------------------------------------------------------------------
         iwdl=jud(lrl)
         iwdr=iwdl
         do mpl=1,mhlp
           iwal=lpnew_lwei(mpl)
           iwar=lpnew_rwei(mpl)
           lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
           lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
         enddo
         w0=w0_dd(2)
         w1=w1_dd(2)
         do mpl=1,mtype
           vplp_w0(mpl)=vplpnew_w0(mpl)*w0
           vplp_w1(mpl)=vplpnew_w1(mpl)*w1
         enddo
         call drl_br_ext_al_new(lin,lrl,lra)

         w0=w0_dd(3)
         w1=0.d0
         do mpl=1,mtype
           vplp_w0(mpl)=vplpnew_w0(mpl)*w0
           vplp_w1(mpl)=vplpnew_w1(mpl)*w1
         enddo
         do lrk=1,norb_dz
           if(lrk.eq.lrl) cycle
           call drl_br_ext_al_new(lin,lrk,lra)
         enddo
       enddo
      return
      end

      subroutine dd_arbl_act_bl(lin,lra)
!dd(19-1) a&r(23)-b&l(32)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      isma=mul_tab(iml,imr)
      do lrl=norb_frz+1,norb_dz-1
        lml=lsm_inn(lrl)
        if(lml.ne.jml) cycle
        do lrr=lrl+1,norb_dz
          lmr=lsm_inn(lrr)
          if(lmr.ne.jmr) cycle
            ijk=lrl-norb_frz+ngw2(lrr-norb_frz)+ngw3(lra-norb_frz)
            intpos=intind_ijka(ijk)
            w0dd1=w0_dd(1)
            w1dd1=w1_dd(1)
            ni=mod(lrr-lrl,2)
            if(ni.eq.0) w0dd1=-w0_dd(1)
            if(ni.eq.0) w1dd1=-w1_dd(1)
           do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0dd1
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1dd1
            enddo

            iwdl=jud(lrl)
            iwdr=jud(lrr)
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
         call ar_bl_bl_ext_ar_new(lin,intpos,isma,1 )
!-------------------------------------------------------------------
        enddo
      enddo
      return
      end

      subroutine dd_arbl_act_br(lin,lra)
!dd(19-1) a&r(23)-b&l(32)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      isma=mul_tab(iml,imr)
      do lrl=norb_frz+1,norb_dz-1
        lml=lsm_inn(lrl)
        if(lml.ne.jml) cycle
        do lrr=lrl+1,norb_dz
          lmr=lsm_inn(lrr)
          if(lmr.ne.jmr) cycle
            ijk=lrl-norb_frz+ngw2(lrr-norb_frz)+ngw3(lra-norb_frz)
            intpos=intind_ijka(ijk)
            w0dd1=w0_dd(1)
            w1dd1=w1_dd(1)
            ni=mod(lrr-lrl,2)
            if(ni.eq.0) w0dd1=-w0_dd(1)
            if(ni.eq.0) w1dd1=-w1_dd(1)
           do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0dd1
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1dd1
            enddo

            iwdl=jud(lrl)
            iwdr=jud(lrr)
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
           call ar_bl_br_ext_al_new(lin,intpos,isma,1 )
!-------------------------------------------------------------------
        enddo
      enddo
      return
      end

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
       subroutine dds_drlbr_act_c_sgt0(lin)
!=======================================================================
!dds(9-4) drl(12)-br(31)-
!=======================================================================
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jmr.or.lmi.ne.jml) cycle
          iwdr=just(lri,lrj)
          iwdl=jud(lri)
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo
          w1ds =w1_d1s(4)
          ni=mod(norb_dz-lrj,2)
          if(ni.eq.1) w1ds=-w1ds
          do mpl=1,mtype
            vplp_w0(mpl)=0.d0
            vplp_w1(mpl)=vplpnew_w0(mpl)*w1ds
          enddo
          call drl_br_ext_al_new(lin,lri,lrj)
        enddo
      enddo

      return
      end

      subroutine dds_ardlr_act_c_sgt0(lin )
!=======================================================================
!dds(9-1) ar(13)-drl(30)-
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
          w0ds1 =w0_d1s(1)
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

      subroutine dds_arblbr_act_c_sgt0(lin )
!=======================================================================
!d1s(9-2)    ar(13)-bl(31)-br(32)-
!d1s(9-3)    ar(13)-bl(32)-br(31)-
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
          do lrd=norb_frz+1,lri-1
            lmd=lsm_inn(lrd)
            if(lmd.ne.jml) cycle
            ijk=lrd-norb_frz+ngw2(lri-norb_frz)+ngw3(lrj-norb_frz)
            intpos=intind_ijka(ijk)
            w0ds2 =w0_d1s(2)
            w1ds2 =w1_d1s(2)
            w0ds3 =w0_d1s(3)
            w1ds3 =w1_d1s(3)
            ni=mod(norb_dz-lrj+lri-lrd,2)
            if(ni.eq.0) then
             w0ds2 =-w0ds2
              w1ds2 =-w1ds2
             w0ds3 =-w0ds3
              w1ds3 =-w1ds3
            endif
!d1s(9-2)    ar(13)-bl(31)-br(32)-
            iwdr=just(lrj,lri)
            iwdl=jud(lrd)
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0ds2
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1ds2
            enddo
           call ar_bl_br_ext_al_new(lin,intpos,isma,1)
!d1s(9-3)    ar(13)-bl(32)-br(31)-
            iwdr=just(lri,lrj)
            iwdl=jud(lrd)
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0ds3
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1ds3
            enddo
           call ar_bl_br_ext_al_new(lin,intpos,isma,1)
          enddo
        enddo
      enddo
      return
      end

      subroutine tttt_drl_act_br_sgt1(lin,lra)
!t1t1(12-2)  drl(11)-
!t1t1(12-2)  drl(11)-c"(11)-
!t1t1(12-3)  drl(33)-
!t1t1(12-3)  drl(33)-c"(11)-
!t1t1(12-3)  drl(33)-c"(11)-c"(11)-
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
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0_t1t1(2)
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1_t1t1(2)
            enddo
!t1t1(12-2)  drl(11)-
!t1t1(12-2)  drl(11)-c"(11)-
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
!t1t1(12-3)  drl(33)-
!t1t1(12-3)  drl(33)-c"(11)-
!t1t1(12-3)  drl(33)-c"(11)-c"(11)-
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0_t1t1(3)
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


      subroutine tttt_arbl_act_br_sgt1(lin,lra)
!=======================================================================
!t1t1(12-1)  (22)ar(13)-bl(31)-
!t1t1(12-1)  ar(13)-c'(11)-bl(31)-
!t1t1(12-1)  ar(13)-bl(31)-c"(11)-
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
            w0tt1=w0_t1t1(1)
            w1tt1=w1_t1t1(1)
            ni=mod(lrj-lri,2)
            if(ni.eq.0) then
              w0tt1=-w0tt1
              w1tt1=-w1tt1
            endif
!-------------------------------------------------------------------
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0tt1
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1tt1
          enddo
!t1t1(12-1)  (22)ar(13)-bl(31)-
          do lrk=norb_frz+1,lri-1
            lmk=lsm_inn(lrk)
            if(mul_tab(lmk,lmi).ne.jml) cycle
           iwdl=just(lrk,lri)      !
            iwdr=just(lrk,lrj)      !
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
          call ar_bl_br_ext_al_new(lin,intpos,isma,1)
          enddo
!t1t1(12-1)  ar(13)-bl(31)-c"(11)-
          do lrk=lrj+1,norb_dz
            lmk=lsm_inn(lrk)
            if(mul_tab(lmk,lmi).ne.jml) cycle
            iwdl=just(lri,lrk)    !
            iwdr=just(lrj,lrk)    !
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
          call ar_bl_br_ext_al_new(lin,intpos,isma,1)
          enddo
!t1t1(12-1)  ar(13)-c'(11)-bl(31)-
          do mpl=1,mtype
            vplp_w0(mpl)=-vplp_w0(mpl)
            vplp_w1(mpl)=-vplp_w1(mpl)
          enddo
          do lrk=lri+1,lrj-1
            lmk=lsm_inn(lrk)
            if(mul_tab(lmk,lmi).ne.jml) cycle
           iwdl=just(lri,lrk)       !
            iwdr=just(lrk,lrj)       !
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
          call ar_bl_br_ext_al_new(lin,intpos,isma,1)
          enddo
        enddo
      enddo
      return
      end

      subroutine ttdd_ar_act_dlr_sgt1(lin,lra)
!t1d1(15-1) (11)a(13)
!t1d1(15-1) a(13)c'(11)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      jmlr=mul_tab(jml,jmr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jmlr) cycle
        w0td1=w0_t1d1(1)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1) w0td1=-w0td1
        ijk=lri-norb_frz+lra
!-------------------------------------------------------------------
!t1d1(15-1) (11)a(13)
          do lrk=norb_frz+1,lri-1
            lmk=lsm_inn(lrk)
           if(lmk.ne.jmr.and.mul_tab(lmk,lmi).ne.jml) cycle
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
!t1d1(15-1) a(13)c'(11)
          do lrk=lri+1,norb_dz
            lmk=lsm_inn(lrk)
           if(lmk.ne.jmr.and.mul_tab(lmk,lmi).ne.jml) cycle
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

      subroutine ttdd_ar_act_blbr_sgt1(lin,jk )
!t1d1(15-1) (11)a(13)
!t1d1(15-1) a(13)c'(11)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      jmlr=mul_tab(jml,jmr)
      isma=mul_tab(iml,imr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jmlr) cycle
        w0td1=w0_t1d1(1)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1) w0td1=-w0td1
        ijk=lri-norb_frz+jk
        intpos=intind_ijka(ijk)
!-------------------------------------------------------------------
!t1d1(15-1) (11)a(13)
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
          call ar_bl_br_ext_al_new(lin,intpos,isma,1)
          enddo
!-------------------------------------------------------------------
!t1d1(15-1) a(13)c'(11)
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
          call ar_bl_br_ext_al_new(lin,intpos,isma,1)
          enddo
        enddo
      return
      end

      subroutine ddtt_arblbr_act_c1_sgt1(lin )
!=======================================================================
!d1t1(16) ar(13)-bl(31)-br(31)-
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
!d1t1(16) ar(13)-bl(31)-br(31)-
          do lrd=norb_frz+1,lri-1
            lmd=lsm_inn(lrd)
            if(lmd.ne.jml) cycle
            ijk=lrd-norb_frz+ngw2(lri-norb_frz)+ngw3(lrj-norb_frz)
            intpos=intind_ijka(ijk)
            w0dt =w0_d1t1
            w1dt =w1_d1t1
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
           call ar_bl_br_ext_al_new(lin,intpos,isma,1)
          enddo
        enddo
      enddo

      return
      end

      subroutine dddd_drl_act_br_sgt0(lin,lra)
!d1d1(20-2) drl(11)-
!d1d1(20-3) drl(33)-c"(11)-
!d1d1(20-3) (11)drl(33)-
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
         w0=w0_d1d1(2)
         w1=w1_d1d1(2)
         do mpl=1,mtype
           vplp_w0(mpl)=vplpnew_w0(mpl)*w0
           vplp_w1(mpl)=vplpnew_w1(mpl)*w1
         enddo
         call drl_br_ext_al_new(lin,lrl,lra)

         w0=w0_d1d1(3)
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

      subroutine dddd_arbl_act_br_sgt0(lin,lra)
!d1d1(20-1) a&r(13)-b&l(31)-
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
            w0dd1=w0_d1d1(1)
            w1dd1=w1_d1d1(1)
            ni=mod(lrr-lrl,2)
            if(ni.eq.0) w0dd1=-w0dd1
            if(ni.eq.0) w1dd1=-w1dd1
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
           call ar_bl_br_ext_al_new(lin,intpos,isma,1)
!-------------------------------------------------------------------
        enddo
      enddo
      return
      end

      subroutine ddv_ar_act_blbr_sgt0(lin,jk )
!=======================================================================
!dv(24-1) ar(13)-
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

        w0dv1=w0_d1v(1)
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
        call ar_bl_br_ext_al_new(lin,intpos,isma,1)
      enddo
      return
      end

      subroutine ddv_ar_act_dlr_sgt0(lin,lra)
!=======================================================================
!d1v(24-1) ar(13)-
!=======================================================================
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jml) cycle

        w0dv1=w0_d1v(1)
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

      subroutine dd1_arbl_act_br_sgt0(lin,lra)
!dd1(21-1) a&r(23)-b&l(31)-
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
!            w0dd1=0.d0
            w1dd1=w1_dd1
            ni=mod(lrr-lrl,2)
            if(ni.eq.0) w1dd1=-w1_dd1
           do mpl=1,mtype
              vplp_w0(mpl)=0.d0
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
           call ar_bl_br_ext_al_new(lin,intpos,isma,1)
!-------------------------------------------------------------------
        enddo
      enddo
      return
      end

      subroutine d1d_arbl_act_br_sgt0(lin,lra)
!d1d(22-1) a&r(13)-b&l(32)-
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
!            w0dd1=0.d0
            w1dd1=w1_d1d(1)
            ni=mod(lrr-lrl,2)
            if(ni.eq.0) w1dd1=-w1dd1
           do mpl=1,mtype
              vplp_w0(mpl)=0.d0
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
           call ar_bl_br_ext_al_new(lin,intpos,isma,1)
!-------------------------------------------------------------------
        enddo
      enddo
      return
      end

      subroutine d1d_drl_act_br_sgt0(lin,lra)
!d1d(22-2) drl(12)-
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
         w1=w1_d1d(2)
         do mpl=1,mtype
           vplp_w0(mpl)=0.d0
           vplp_w1(mpl)=vplpnew_w1(mpl)*w1
         enddo
         call drl_br_ext_al_new(lin,lrl,lra)
      enddo

      return
      end

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! end of vd, next dd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine ss_arbr_act_c_dd_ext_sgt0()
!ss(1-1)  ar(01)-bl(32)-        act -c"-
!ss(1-3)  ar(13)-bl(20)-        act -c"-
!ss(1-6)  (11)-ar(23)-bl(32)-   act -c"-
!ss(1-7)  ar(13)-c'(21)-bl(32)- act -c"-
!ss(1-8)  ar(13)-c'(22)-bl(31)- act -c"-
!ss(1-9)  ar(23)-c'(11)-bl(32)- act -c"-
!ss(1-11) ar(13)-bl(31)-c"(22)- act -c"-
!ss(1-12) ar(13)-bl(32)-c"(21)- act -c"-
!ss(1-13) ar(23)-bl(31)-c"(12)- act -c"-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      do lri=norb_frz+1,norb_dz-1
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          w0ss1=w0_ss(1)
          w1ss1=w1_ss(1)
          w0ss3=w0_ss(3)
          w1ss3=w1_ss(3)
          w0ss6=w0_ss(6)
          w1ss6=w1_ss(6)
          w0ss7=w0_ss(7)
          w1ss7=w1_ss(7)
          w0ss8=w0_ss(8)
          w1ss8=w1_ss(8)
          w0ss9=w0_ss(9)
          w1ss9=w1_ss(9)
          w0ss11=w0_ss(11)
          w1ss11=w1_ss(11)
          w0ss12=w0_ss(12)
          w1ss12=w1_ss(12)
          w0ss13=w0_ss(13)
          w1ss13=w1_ss(13)
          ni=mod(lrj-lri,2)
          if(ni.eq.0) then
            w0ss1=-w0ss1
            w1ss1=-w1ss1
            w0ss3=-w0ss3
            w1ss3=-w1ss3
            w0ss6=-w0ss6
            w1ss6=-w1ss6
            w0ss7=-w0ss7
            w1ss7=-w1ss7
            w0ss8=-w0ss8
            w1ss8=-w1ss8
            w0ss9=-w0ss9
            w1ss9=-w1ss9
            w0ss11=-w0ss11
            w1ss11=-w1ss11
            w0ss12=-w0ss12
            w1ss12=-w1ss12
            w0ss13=-w0ss13
            w1ss13=-w1ss13
          endif
!ss(1-1)  ar(01)-bl(32)-        act -c"-
          if(jml.eq.1.and.lmij.eq.jmr) then
            iwdl=just(lri,lri)
            iwdr=just(lrj,lri)
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
              do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0ss1
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1ss1
            enddo
          call ar_bl_dd_ext(lri,lrj,1)
          endif
!ss(1-3)  ar(13)-bl(20)-        act -c"-
          if(jmr.eq.1.and.lmij.eq.jml) then
            iwdl=just(lrj,lri)
            iwdr=just(lrj,lrj)
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
              do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0ss3
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1ss3
            enddo
          call ar_bl_dd_ext(lri,lrj,1)
          endif
!ss(1-6)  (11)-ar(23)-bl(32)-   act -c"-
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0ss6
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1ss6
          enddo
          do lrk=norb_frz+1,lri-1
            lmk=lsm_inn(lrk)
            lmki=mul_tab(lmk,lmi)
            lmkj=mul_tab(lmk,lmj)
              if(lmki.eq.jml.and.lmkj.eq.jmr) then
              iwdl=just(lri,lrk)
              iwdr=just(lrj,lrk)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
            call ar_bl_dd_ext(lri,lrj,1)
            endif
            enddo
!ss(1-7)  ar(13)-c'(21)-bl(32)- act -c"-
          do mpl=1,mtype
            vplp_w0(mpl)=-vplpnew_w0(mpl)*w0ss7
            vplp_w1(mpl)=-vplpnew_w1(mpl)*w1ss7
          enddo
          do lrk=lri+1,lrj-1
            lmk=lsm_inn(lrk)
            lmki=mul_tab(lmk,lmi)
            lmkj=mul_tab(lmk,lmj)
            if(lmki.eq.jml.and.lmkj.eq.jmr) then
              iwdl=just(lrk,lri)
              iwdr=just(lrj,lrk)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
              call ar_bl_dd_ext(lri,lrj,1)
            endif
          enddo
!ss(1-8)  ar(13)-c'(22)-bl(31)- act -c"-
          do mpl=1,mtype
            vplp_w0(mpl)=-vplpnew_w0(mpl)*w0ss8
            vplp_w1(mpl)=-vplpnew_w1(mpl)*w1ss8
          enddo
          do lrk=lri+1,lrj-1
            lmk=lsm_inn(lrk)
            lmki=mul_tab(lmk,lmi)
            lmkj=mul_tab(lmk,lmj)
            if(lmki.eq.jml.and.lmkj.eq.jmr) then
              iwdl=just(lrk,lri)
              iwdr=just(lrk,lrj)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
              call ar_bl_dd_ext(lri,lrj,1)
            endif
          enddo
!ss(1-9)  ar(23)-c'(11)-bl(32)- act -c"-
          do mpl=1,mtype
            vplp_w0(mpl)=-vplpnew_w0(mpl)*w0ss9
            vplp_w1(mpl)=-vplpnew_w1(mpl)*w1ss9
          enddo
          do lrk=lri+1,lrj-1
            lmk=lsm_inn(lrk)
            lmki=mul_tab(lmk,lmi)
            lmkj=mul_tab(lmk,lmj)
            if(lmki.eq.jml.and.lmkj.eq.jmr) then
              iwdl=just(lri,lrk)
              iwdr=just(lrj,lrk)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
              call ar_bl_dd_ext(lri,lrj,1)
            endif
          enddo
!ss(1-11) ar(13)-bl(31)-c"(22)- act -c"-
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0ss11
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1ss11
          enddo
          do lrk=lrj+1,norb_dz
            lmk=lsm_inn(lrk)
            lmki=mul_tab(lmk,lmi)
            lmkj=mul_tab(lmk,lmj)
            if(lmki.eq.jml.and.lmkj.eq.jmr) then
              iwdl=just(lrk,lri)
              iwdr=just(lrk,lrj)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
              call ar_bl_dd_ext(lri,lrj,1)
            endif
          enddo
!ss(1-12) ar(13)-bl(32)-c"(21)- act -c"-
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0ss12
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1ss12
          enddo
          do lrk=lrj+1,norb_dz
            lmk=lsm_inn(lrk)
            lmki=mul_tab(lmk,lmi)
            lmkj=mul_tab(lmk,lmj)
            if(lmki.eq.jml.and.lmkj.eq.jmr) then
              iwdl=just(lrk,lri)
              iwdr=just(lrj,lrk)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
              call ar_bl_dd_ext(lri,lrj,1)
            endif
          enddo
!ss(1-13) ar(23)-bl(31)-c"(12)- act -c"-
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0ss13
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1ss13
          enddo
          do lrk=lrj+1,norb_dz
            lmk=lsm_inn(lrk)
            lmki=mul_tab(lmk,lmi)
            lmkj=mul_tab(lmk,lmj)
            if(lmki.eq.jml.and.lmkj.eq.jmr) then
              iwdl=just(lri,lrk)
              iwdr=just(lrk,lrj)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
              call ar_bl_dd_ext(lri,lrj,1)
            endif
          enddo
        enddo
      enddo

      return
      end

      subroutine ss_arbr_act_c_dd_ext()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      do lri=norb_frz+1,norb_dz-1
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          w0ss2=w0_ss(2)
          w1ss2=w1_ss(2)
          w0ss4=w0_ss(4)
          w1ss4=w1_ss(4)
          w0ss5=w0_ss(5)
          w1ss5=w1_ss(5)
          w0ss10=-w0_ss(10)
          w1ss10=-w1_ss(10)
          w0ss14=w0_ss(14)
          w1ss14=w1_ss(14)
          ni=mod(lrj-lri,2)
          if(ni.eq.0) then
            w0ss2=-w0ss2
            w1ss2=-w1ss2
            w0ss4=-w0ss4
            w1ss4=-w1ss4
            w0ss5=-w0ss5
            w1ss5=-w1ss5
            w0ss10=-w0ss10
            w1ss10=-w1ss10
            w0ss14=-w0ss14
            w1ss14=-w1ss14
          endif
          if(jml.eq.1.and.lmij.eq.jmr) then
!ss(1-2)  ar(02)-bl(31)-        act -c"-
            iwdl=just(lri,lri)
            iwdr=just(lri,lrj)
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
              do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0ss2
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1ss2
            enddo
          call ar_bl_dd_ext(lri,lrj,1)
          endif
          if(jmr.eq.1.and.lmij.eq.jml) then
!ss(1-4)  ar(23)-bl(10)-        act -c"-                         ! iprad
            iwdl=just(lri,lrj)
            iwdr=just(lrj,lrj)
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
              do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0ss4
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1ss4
            enddo
          call ar_bl_dd_ext(lri,lrj,1)
          endif
!ss(1-5)  (22)-ar(13)-bl(31)-   act -c"-
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0ss5
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1ss5
          enddo
          do lrk=norb_frz+1,lri-1
            lmk=lsm_inn(lrk)
            lmki=mul_tab(lmk,lmi)
            lmkj=mul_tab(lmk,lmj)
              if(lmki.eq.jml.and.lmkj.eq.jmr) then
              iwdl=just(lrk,lri)
              iwdr=just(lrk,lrj)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
            call ar_bl_dd_ext(lri,lrj,1)
            endif
            enddo
!ss(1-10) ar(23)-c'(12)-bl(31)- act -c"-
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0ss10
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1ss10
          enddo
          do lrk=lri+1,lrj-1
            lmk=lsm_inn(lrk)
            lmki=mul_tab(lmk,lmi)
            lmkj=mul_tab(lmk,lmj)
              if(lmki.eq.jml.and.lmkj.eq.jmr) then
              iwdl=just(lri,lrk)
              iwdr=just(lrk,lrj)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
            call ar_bl_dd_ext(lri,lrj,1)
            endif
          enddo
!ss(1-14) ar(23)-bl(32)-c"(11)- act -c"-
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0ss14
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1ss14
          enddo
          do lrk=lrj+1,norb_dz
            lmk=lsm_inn(lrk)
            lmki=mul_tab(lmk,lmi)
            lmkj=mul_tab(lmk,lmj)
            if(lmki.eq.jml.and.lmkj.eq.jmr) then
              iwdl=just(lri,lrk)
              iwdr=just(lrj,lrk)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
            call ar_bl_dd_ext(lri,lrj,1)
            endif
          enddo
        enddo
      enddo

      return
      end


      subroutine ss_drl_act_c_dd_ext_sgt0()
!ss(1-18) drl(11)-c"(22)-       act -c"-
!ss(1-19) drl(12)-c"(21)-       act -c"-
!ss(1-20) (11)-drl(33)-c"(22)-  act -c"-
!ss(1-20) drl(33)-c"(11)-c"(22)-act -c"-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      if(jml.ne.jmr) return
      do lri=norb_frz+1,norb_dz-1
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jml) cycle
          iwdl=just(lrj,lri)
          iwdr=iwdl
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo
!ss(1-16) (11)-drl(22)-         act -c"-
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0_ss(16)
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1_ss(16)
          enddo
          call drl_dd_ext(lrj)
!ss(1-18) drl(11)-c"(22)-       act -c"-
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0_ss(18)
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1_ss(18)
          enddo
          call drl_dd_ext(lri)
!ss(1-20) (11)(22)drl(33)-  act -c"-
!ss(1-20) (11)-drl(33)-c"(22)-  act -c"-
!ss(1-20) drl(33)-c"(11)-c"(22)-act -c"-
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0_ss(20)
            vplp_w1(mpl)=0.d0
          enddo
          do lrk=1,norb_dz
            if(lrk.eq.lri) cycle
            if(lrk.eq.lrj) cycle
            call drl_dd_ext(lrk)
          enddo
        enddo
      enddo

      return
      end

      subroutine ss_s_drl_act_c_dd_ext_sgt0()
!ss(1-19) drl(12)-c"(21)-       act -c"-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      if(jml.ne.jmr) return
      do lri=norb_frz+1,norb_dz-1
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jml) cycle
!ss(1-19) drl(12)-c"(21)-       act -c"-
          iwdl=just(lrj,lri)
          iwdr=just(lri,lrj)
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo
          do mpl=1,mtype
            vplp_w0(mpl)=0.d0
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1_ss(19)
          enddo
          call drl_dd_ext(lri)
        enddo
      enddo

      return
      end

      subroutine ss_drl_act_c_dd_ext()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      w0ss15=w0_ss(15)
      w1ss15=w1_ss(15)
      w0ss17=w0_ss(17)
      w1ss17=w1_ss(17)
      w0ss20=w0_ss(20)
      if(jml.eq.1.and.jmr.eq.1) then
!ss(1-20) drl(33)-c"(00)-       act -c"-                     ! ipl(r)ad=
        do lr0=norb_frz+1,norb_dz
          iwdl=just(lr0,lr0)
          iwdr=iwdl
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo
          do lrk=1,norb_dz
            if(lrk.eq.lr0) cycle
              do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0ss20
              vplp_w1(mpl)=0.d0
            enddo
            call drl_dd_ext(lrk)
          enddo
        enddo
      endif
      do lri=norb_frz+1,norb_dz-1
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jml.or.lmij.ne.jmr) cycle
          iwdl=just(lri,lrj)
          iwdr=iwdl
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo
!ss(1-15) (22)-drl(11)-         act -c"-
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0ss15
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1ss15
          enddo
          call drl_dd_ext(lrj)
!ss(1-17) drl(22)-c"(11)-       act -c"-
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0ss17
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1ss17
          enddo
          call drl_dd_ext(lri)
!ss(1-20) (22)(11)drl(33)-      act -c"-
!ss(1-20) (22)drl(33)-c"(11)-   act -c"-
!ss(1-20) drl(33)-c"(22)-c"(11)-act -c"-
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0ss20
            vplp_w1(mpl)=0.d0
          enddo
          do lrk=1,norb_dz
            if(lrk.eq.lri) cycle
            if(lrk.eq.lrj) cycle
            call drl_dd_ext(lrk)
          enddo
        enddo
      enddo

      return
      end

      subroutine st_arbl_act_c_dd_ext_sgt0()
!st(2-3) ar(13)-c'(22)-bl(32)-   act -c"-
!st(2-3) ar(13)-bl(32)-c'(22)-   act -c"-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      do lri=norb_frz+1,norb_dz-1
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          w1st3=w1_st(3)
          ni=mod(lrj-lri,2)
          if(ni.eq.0) then
            w1st3=-w1st3
          endif
          do lrk=lri+1,lrj-1
            lmk=lsm_inn(lrk)
            lmki=mul_tab(lmk,lmi)
            lmkj=mul_tab(lmk,lmj)
            if(lmki.eq.jml.and.lmkj.eq.jmr) then
!st(2-3) ar(13)-c'(22)-bl(32)-   act -c"-
              iwdl=just(lrk,lri)
              iwdr=just(lrk,lrj)      !
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
                do mpl=1,mtype
                vplp_w0(mpl)=0.d0
                vplp_w1(mpl)=-vplpnew_w1(mpl)*w1st3
              enddo
              call ar_bl_dd_ext(lri,lrj,1)
            endif
          enddo
          do lrk=lrj+1,norb_dz
            lmk=lsm_inn(lrk)
            lmki=mul_tab(lmk,lmi)
            lmkj=mul_tab(lmk,lmj)
            if(lmki.eq.jml.and.lmkj.eq.jmr) then
!st(2-3) ar(13)-bl(32)-c'(22)-   act -c"-
              iwdl=just(lrk,lri)
              iwdr=just(lrj,lrk)      !
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
              do mpl=1,mtype
                vplp_w0(mpl)=0.d0
                vplp_w1(mpl)=vplpnew_w1(mpl)*w1st3
              enddo
              call ar_bl_dd_ext(lri,lrj,1)
            endif
          enddo
        enddo
      enddo

      return
      end

      subroutine st_arbl_act_c_dd_ext()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      do lri=norb_frz+1,norb_dz-1
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          w1st1=w1_st(1)
          w1st2=w1_st(2)
          w1st4=w1_st(4)
          ni=mod(lrj-lri,2)
          if(ni.eq.0) then
            w1st1=-w1st1
            w1st2=-w1st2
            w1st4=-w1st4
          endif
          if(jml.eq.1.and.lmij.eq.jmr) then
!st(2-1) ar(02)-bl(32)-          act -c"-
            iwdl=just(lri,lri)
            iwdr=just(lri,lrj)      !
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
              do mpl=1,mtype
              vplp_w0(mpl)=0.d0
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1st1
            enddo
            call ar_bl_dd_ext(lri,lrj,1)
          endif
          do lrk=norb_frz+1,lri-1
            lmk=lsm_inn(lrk)
            lmki=mul_tab(lmk,lmi)
            lmkj=mul_tab(lmk,lmj)
              if(lmki.eq.jml.and.lmkj.eq.jmr) then
!st(2-2) (22)ar(13)-bl(32)-      act -c"-
              iwdl=just(lrk,lri)
              iwdr=just(lrk,lrj)         !
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
                do mpl=1,mtype
                vplp_w0(mpl)=0.d0
                vplp_w1(mpl)=vplpnew_w1(mpl)*w1st2
              enddo
              call ar_bl_dd_ext(lri,lrj,1)
            endif
            enddo
            do lrk=lri+1,lrj-1
              lmk=lsm_inn(lrk)
              lmki=mul_tab(lmk,lmi)
              lmkj=mul_tab(lmk,lmj)
              if(lmki.eq.jml.and.lmkj.eq.jmr) then
!st(2-4) ar(23)-c'(12)-bl(32)-   act -c"-
              iwdl=just(lri,lrk)
              iwdr=just(lrk,lrj)      !
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
                do mpl=1,mtype
                vplp_w0(mpl)=0.d0
                vplp_w1(mpl)=vplpnew_w1(mpl)*w1st4
              enddo
              call ar_bl_dd_ext(lri,lrj,1)
            endif
          enddo
          do lrk=lrj+1,norb_dz
            lmk=lsm_inn(lrk)
            lmki=mul_tab(lmk,lmi)
            lmkj=mul_tab(lmk,lmj)
            if(lmki.eq.jml.and.lmkj.eq.jmr) then
!st(2-4) ar(23)-bl(32)-c'(12)-   act -c"-
              iwdl=just(lri,lrk)
              iwdr=just(lrj,lrk)      !
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
                do mpl=1,mtype
                vplp_w0(mpl)=0.d0
                vplp_w1(mpl)=vplpnew_w1(mpl)*w1st4
              enddo
            call ar_bl_dd_ext(lri,lrj,1)
            endif
          enddo
        enddo
      enddo

      return
      end

      subroutine st_drl_act_c_dd_ext_sgt0()
!----------------------------------------------------------
!st(2-7) drl(12)-c"(22)-         act -c"-
!----------------------------------------------------------
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      do lri=norb_frz+1,norb_dz-1
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
!st(2-7) drl(12)-c"(22)-         act -c"-
          if(lmij.ne.jml.or.lmij.ne.jmr) cycle
          iwdl=just(lrj,lri)
          iwdr=just(lri,lrj)      !
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo
          do mpl=1,mtype
            vplp_w0(mpl)=0.d0
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1_st(7)
          enddo
          call drl_dd_ext(lri)
        enddo
      enddo

      return
      end

      subroutine st_drl_act_c_dd_ext()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      do lri=norb_frz+1,norb_dz-1
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jml.or.lmij.ne.jmr) cycle
!st(2-5) (22)drl(12)-          act -c"-
          iwdl=just(lri,lrj)
          iwdr=just(lrj,lri)
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
          call drl_dd_ext(lrj)
!st(2-6) drl(22)-c"(12)-       act -c"-
          do mpl=1,mtype
            vplp_w0(mpl)=0.d0
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1_st(6)
          enddo
          call drl_dd_ext(lri)
        enddo
      enddo

      return
      end

      subroutine ts_arbl_act_c_dd_exe_sgt0()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      do lri=norb_frz+1,norb_dz-1
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          w1ts3=w1_ts(3)
          ni=mod(lrj-lri,2)
          if(ni.eq.0) then
            w1ts3=-w1ts3
          endif
!ts(3-3) ar(23)-bl(31)-c"(22)-   act -c"-
          do mpl=1,mtype
            vplp_w0(mpl)=0.d0
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1ts3
          enddo
          do lrk=lrj+1,norb_dz
            lmk=lsm_inn(lrk)
            lmki=mul_tab(lmk,lmi)
            lmkj=mul_tab(lmk,lmj)
            if(lmki.eq.jml.and.lmkj.eq.jmr) then
              iwdl=just(lri,lrk)          !
              iwdr=just(lrk,lrj)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
            call ar_bl_dd_ext(lri,lrj,1)
            endif
          enddo
        enddo
      enddo

      return
      end

      subroutine stt_arbl_act_c_dd_ext_sgt1()
!st1(4-1) ar(01)-bl(31)-
!st1(4-2) ar(23)-bl(31)-
!st1(4-3) ar(13)-c'(21)-bl(31)-
!st1(4-3) ar(13)-bl(31)-c"(21)-
!st1(4-4) ar(23)-c'(11)-bl(31)-
!st1(4-4) ar(23)-bl(31)-c"(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      do lri=norb_frz+1,norb_dz-1
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          w1ts1=w1_st1(1)
          w1ts2=w1_st1(2)
          w1ts3=w1_st1(3)
          w1ts4=w1_st1(4)
          ni=mod(lrj-lri,2)
          if(ni.eq.0) then
            w1ts1=-w1ts1
            w1ts2=-w1ts2
            w1ts3=-w1ts3
            w1ts4=-w1ts4
          endif
          if(jml.eq.1.and.lmij.eq.jmr) then
!st1(4-1) ar(01)-bl(31)-
            iwdl=just(lri,lri)   !
            iwdr=just(lri,lrj)
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            do mpl=1,mtype
              vplp_w0(mpl)=0.d0
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1ts1
            enddo
            call ar_bl_dd_ext(lri,lrj,1)
          endif
          do lrk=norb_frz+1,lri-1
            lmk=lsm_inn(lrk)
            lmki=mul_tab(lmk,lmi)
            lmkj=mul_tab(lmk,lmj)
            if(lmki.eq.jml.and.lmkj.eq.jmr) then
!st1(4-2) (11)ar(23)-bl(31)-
              iwdl=just(lri,lrk)      !
              iwdr=just(lrk,lrj)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
              do mpl=1,mtype
                vplp_w0(mpl)=0.d0
                vplp_w1(mpl)=vplpnew_w1(mpl)*w1ts2
              enddo
              call ar_bl_dd_ext(lri,lrj,1)
            endif
          enddo
          do lrk=lri+1,lrj-1
            lmk=lsm_inn(lrk)
            lmki=mul_tab(lmk,lmi)
            lmkj=mul_tab(lmk,lmj)
            if(lmki.eq.jml.and.lmkj.eq.jmr) then
!st1(4-3) ar(13)-c'(21)-bl(31)-
!st1(4-4) ar(23)-c'(11)-bl(31)-
              iwdl=just(lrk,lri)        !
              iwdr=just(lrk,lrj)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
              do mpl=1,mtype
                vplp_w0(mpl)=0.d0
                vplp_w1(mpl)=-vplpnew_w1(mpl)*w1ts3
              enddo
              call ar_bl_dd_ext(lri,lrj,1)
              iwdl=just(lri,lrk)        !
              iwdr=just(lrk,lrj)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
              do mpl=1,mtype
                vplp_w0(mpl)=0.d0
                vplp_w1(mpl)=-vplpnew_w1(mpl)*w1ts4
              enddo
              call ar_bl_dd_ext(lri,lrj,1)
            endif
          enddo
!st1(4-3) ar(13)-bl(31)-c"(21)-
!st1(4-4) ar(23)-bl(31)-c"(11)-
          do lrk=lrj+1,norb_dz
            lmk=lsm_inn(lrk)
            lmki=mul_tab(lmk,lmi)
            lmkj=mul_tab(lmk,lmj)
            if(lmki.eq.jml.and.lmkj.eq.jmr) then
              iwdl=just(lrk,lri)          !
              iwdr=just(lrj,lrk)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
              do mpl=1,mtype
               vplp_w0(mpl)=0.d0
               vplp_w1(mpl)=vplpnew_w1(mpl)*w1ts3
              enddo
              call ar_bl_dd_ext(lri,lrj,1)
              iwdl=just(lri,lrk)          !
              iwdr=just(lrj,lrk)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
              do mpl=1,mtype
               vplp_w0(mpl)=0.d0
               vplp_w1(mpl)=vplpnew_w1(mpl)*w1ts4
              enddo
              call ar_bl_dd_ext(lri,lrj,1)
            endif
          enddo
        enddo
      enddo

      return
      end

      subroutine tts_arbl_act_c_dd_ext_sgt1()
!t1s(5-1)   ar(13)-bl(10)-
!t1s(5-2)   ar(13)-bl(32)-
!t1s(5-2)   ar(13)-c'(11)-bl(32)-
!t1s(5-3)   ar(13)-bl(31)-c"(12)-
!t1s(5-4)   ar(13)-bl(32)-c"(11)-
!t1s(5-5)   drl(12)-
!t1s(5-6)   drl(12)-c"(12)-
!t1s(5-7)   drl(12)-c"(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      do lri=norb_frz+1,norb_dz-1
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          w1ts1=w1_t1s(1)
          w1ts2=w1_t1s(2)
          w1ts3=w1_t1s(3)
          w1ts4=w1_t1s(4)
          ni=mod(lrj-lri,2)
          if(ni.eq.0) then
            w1ts1=-w1ts1
            w1ts2=-w1ts2
            w1ts3=-w1ts3
            w1ts4=-w1ts4
          endif
          if(jmr.eq.1.and.lmij.eq.jml) then
!t1s(5-1)   ar(13)-bl(10)-
            iwdl=just(lri,lrj)   !
            iwdr=just(lrj,lrj)
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
              do mpl=1,mtype
              vplp_w0(mpl)=0.d0
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1ts1
            enddo
          call ar_bl_dd_ext(lri,lrj,1)
          endif
          do lrk=norb_frz+1,lri-1
            lmk=lsm_inn(lrk)
            lmki=mul_tab(lmk,lmi)
            lmkj=mul_tab(lmk,lmj)
            if(lmki.eq.jml.and.lmkj.eq.jmr) then
!t1s(5-2)   (11)ar(13)-bl(32)-
              iwdl=just(lrk,lri)      !
              iwdr=just(lrj,lrk)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
                do mpl=1,mtype
                vplp_w0(mpl)=0.d0
                vplp_w1(mpl)=vplpnew_w1(mpl)*w1ts2
              enddo
            call ar_bl_dd_ext(lri,lrj,1)
            endif
            enddo
            do lrk=lri+1,lrj-1
             lmk=lsm_inn(lrk)
              lmki=mul_tab(lmk,lmi)
              lmkj=mul_tab(lmk,lmj)
              if(lmki.eq.jml.and.lmkj.eq.jmr) then
!t1s(5-2)   ar(13)-c'(11)-bl(32)-
              iwdl=just(lri,lrk)        !
              iwdr=just(lrj,lrk)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
                do mpl=1,mtype
                vplp_w0(mpl)=0.d0
                vplp_w1(mpl)=-vplpnew_w1(mpl)*w1ts2
              enddo
            call ar_bl_dd_ext(lri,lrj,1)
            endif
          enddo
!t1s(5-3)   ar(13)-bl(31)-c"(12)-
          do mpl=1,mtype
            vplp_w0(mpl)=0.d0
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1ts3
          enddo
          do lrk=lrj+1,norb_dz
            lmk=lsm_inn(lrk)
            lmki=mul_tab(lmk,lmi)
            lmkj=mul_tab(lmk,lmj)
            if(lmki.eq.jml.and.lmkj.eq.jmr) then
              iwdl=just(lri,lrk)          !
              iwdr=just(lrk,lrj)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
            call ar_bl_dd_ext(lri,lrj,1)
            endif
          enddo
!t1s(5-4)   ar(13)-bl(32)-c"(11)-
          do mpl=1,mtype
            vplp_w0(mpl)=0.d0
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1ts4
          enddo
          do lrk=lrj+1,norb_dz
            lmk=lsm_inn(lrk)
            lmki=mul_tab(lmk,lmi)
            lmkj=mul_tab(lmk,lmj)
            if(lmki.eq.jml.and.lmkj.eq.jmr) then
              iwdl=just(lri,lrk)          !
              iwdr=just(lrj,lrk)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
            call ar_bl_dd_ext(lri,lrj,1)
            endif
          enddo
        enddo
      enddo

      return
      end

      subroutine tts_drl_act_c_dd_ext_sgt1()
!t1s(5-5)   (11)drl(12)-
!t1s(5-6)   drl(12)-c"(12)-
!t1s(5-7)   drl(12)-c"(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      if(jml.ne.jmr) return
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
      do lrj=norb_frz+1,lri-1
        lmj=lsm_inn(lrj)
        lmij=mul_tab(lmi,lmj)
        if(jml.ne.lmij) cycle
!t1s(5-5)   (11)drl(12)-
!t1s(5-6)   drl(11)-c"(12)-
         iwdl=just(lrj,lri)          !
         iwdr=just(lri,lrj)
         do mpl=1,mhlp
           iwal=lpnew_lwei(mpl)
           iwar=lpnew_rwei(mpl)
           lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
           lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
         enddo
         do mpl=1,mtype
           vplp_w0(mpl)=0.d0
           vplp_w1(mpl)=vplpnew_w1(mpl)*w1_t1s(5)
         enddo
         call drl_dd_ext(lri)
         do mpl=1,mtype
           vplp_w0(mpl)=0.d0
           vplp_w1(mpl)=vplpnew_w1(mpl)*w1_t1s(6)
         enddo
         call drl_dd_ext(lrj)
!t1s(5-7)   drl(12)-c"(11)-
         iwdl=just(lrj,lri)          !
         iwdr=iwdl
         do mpl=1,mhlp
           iwal=lpnew_lwei(mpl)
           iwar=lpnew_rwei(mpl)
           lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
           lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
         enddo
         do mpl=1,mtype
           vplp_w0(mpl)=0.d0
           vplp_w1(mpl)=vplpnew_w1(mpl)*w1_t1s(7)
         enddo
         call drl_dd_ext(lrj)
      enddo
      enddo

      return
      end

      subroutine sd_ar_act_bl_dd_ext_sgt0(lra)
!sd(6-3) a&r(13)c'(22)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      jmlr=mul_tab(jml,jmr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jmlr) cycle
        w0sd3=w0_sd(3)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1) w0sd3=-w0sd3
!sd(6-3) a&r(13)-c'(22)-
        do lrk=lri+1,norb_dz
          lmk=lsm_inn(lrk)
          lmki=mul_tab(lmk,lmi)
          if(lmki.eq.jml.and.lmk.eq.jmr) then
            iwdl=just(lrk,lri)
            iwdr=jud(lrk)
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            do mpl=1,mtype
              vplp_w0(mpl)=-vplpnew_w0(mpl)*w0sd3
              vplp_w1(mpl)=-vplpnew_w1(mpl)*w0sd3
            enddo
          call ar_bl_dd_ext(lri,lra,1)
          endif
        enddo
      enddo

      return
      end


      subroutine sdd_ar_act_bl_dd_ext_sgt0(lra)
!sd1(8-1)    ar(01)-
!sd1(8-2)    (11)ar(23)-
!sd1(8-3)    ar(13)-c'(21)-
!sd1(8-4)    ar(23)-c'(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      jmlr=mul_tab(jml,jmr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jmlr) cycle
        w0sd1=w0_sd1(1)
        w0sd2=w0_sd1(2)
        w0sd3=w0_sd1(3)
        w0sd4=w0_sd1(4)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1)w0sd1=-w0sd1
        if(ni.eq.1)w0sd2=-w0sd2
        if(ni.eq.1)w0sd3=-w0sd3
        if(ni.eq.1)w0sd4=-w0sd4
        if(jml.eq.1.and.lmi.eq.jmr) then
!sd1(8-1)    ar(01)-
          iwdl=just(lri,lri)
          iwdr=jud(lri)
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd1
            vplp_w1(mpl)=vplpnew_w1(mpl)*w0sd1
          enddo
          call ar_bl_dd_ext(lri,lra,1)
        endif
!sd1(8-2)    (11)ar(23)-
        do lrk=norb_frz+1,lri-1
          lmk=lsm_inn(lrk)
          lmki=mul_tab(lmk,lmi)
          if(lmki.eq.jml.and.lmk.eq.jmr) then
            iwdl=just(lri,lrk)
            iwdr=jud(lrk)
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd2
              vplp_w1(mpl)=vplpnew_w1(mpl)*w0sd2
            enddo
            call ar_bl_dd_ext(lri,lra,1)
          endif
        enddo
!sd1(8-3)    ar(13)-c'(21)-
        do mpl=1,mtype
          vplp_w0(mpl)=-vplpnew_w0(mpl)*w0sd3           ! - ????
          vplp_w1(mpl)=-vplpnew_w1(mpl)*w0sd3           ! - ????
        enddo
        do lrk=lri+1,norb_dz
          lmk=lsm_inn(lrk)
          lmki=mul_tab(lmk,lmi)
          if(lmki.eq.jml.and.lmk.eq.jmr) then
            iwdl=just(lrk,lri)
            iwdr=jud(lrk)
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
          call ar_bl_dd_ext(lri,lra,1)
          endif
        enddo
!sd1(8-4)    ar(23)-c'(11)-
        do mpl=1,mtype
          vplp_w0(mpl)=-vplpnew_w0(mpl)*w0sd4           ! - ????
          vplp_w1(mpl)=-vplpnew_w1(mpl)*w0sd4           ! - ????
        enddo
        do lrk=lri+1,norb_dz
          lmk=lsm_inn(lrk)
          lmki=mul_tab(lmk,lmi)
          if(lmki.eq.jml.and.lmk.eq.jmr) then
            iwdl=just(lri,lrk)
            iwdr=jud(lrk)
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
          call ar_bl_dd_ext(lri,lra,1)
          endif
        enddo
      enddo

      return
      end

      subroutine tttt_arbl_act_c_dd_ext_sgt1()
!t1t1(12-1)  ar(13)-bl(31)-
!t1t1(12-1)  ar(13)-c'(11)-bl(31)-
!t1t1(12-1)  ar(13)-bl(31)-c"(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      do lri=norb_frz+1,norb_dz-1
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          w0tt1=w0_t1t1(1)
          w1tt1=w1_t1t1(1)
          ni=mod(lrj-lri,2)
          if(ni.eq.0) then
            w0tt1=-w0tt1
            w1tt1=-w1tt1
          endif
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0tt1
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1tt1
          enddo
!t1t1(12-1)  (11)ar(13)-bl(31)-
          do lrk=norb_frz+1,lri-1
            lmk=lsm_inn(lrk)
            lmki=mul_tab(lmk,lmi)
            lmkj=mul_tab(lmk,lmj)
            if(lmki.eq.jml.and.lmkj.eq.jmr) then
              iwdl=just(lrk,lri)       !
              iwdr=just(lrk,lrj)       !
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
            call ar_bl_dd_ext(lri,lrj,1)
            endif
          enddo
!t1t1(12-1)  ar(13)-bl(31)-c"(11)-
          do lrk=lrj+1,norb_dz
            lmk=lsm_inn(lrk)
            lmki=mul_tab(lmk,lmi)
            lmkj=mul_tab(lmk,lmj)
            if(lmki.eq.jml.and.lmkj.eq.jmr) then
              iwdl=just(lri,lrk)     !
              iwdr=just(lrj,lrk)     !
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
            call ar_bl_dd_ext(lri,lrj,1)
            endif
          enddo
          do mpl=1,mtype
            vplp_w0(mpl)=-vplp_w0(mpl)
            vplp_w1(mpl)=-vplp_w1(mpl)
          enddo
!t1t1(12-1)  ar(13)-c'(11)-bl(31)-
          do lrk=lri+1,lrj-1
           lmk=lsm_inn(lrk)
            lmki=mul_tab(lmk,lmi)
            lmkj=mul_tab(lmk,lmj)
            if(lmki.eq.jml.and.lmkj.eq.jmr) then
              iwdl=just(lri,lrk)      !
              iwdr=just(lrk,lrj)      !
            do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
            call ar_bl_dd_ext(lri,lrj,1)
            endif
          enddo
        enddo
      enddo

      return
      end

      subroutine tttt_drl_act_c_dd_ext_sgt1()
!t1t1(12-2)  drl(11)-
!t1t1(12-2)  drl(11)-c"(11)-
!t1t1(12-3)  drl(33)-
!t1t1(12-3)  drl(33)-c"(11)-
!t1t1(12-3)  drl(33)-c"(11)-c"(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      do lri=norb_frz+1,norb_dz-1
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jml.or.lmij.ne.jmr) cycle
!t1t1(12-2)  drl(11)-
!t1t1(12-2)  drl(11)-c"(11)-
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0_t1t1(2)
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1_t1t1(2)
            enddo
            iwdl=just(lri,lrj)   !
            iwdr=iwdl
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            call drl_dd_ext(lri)
            call drl_dd_ext(lrj)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0_t1t1(3)
              vplp_w1(mpl)=0.d0
            enddo
            do lrk=1,norb_dz
              if(lrk.eq.lri) cycle
              if(lrk.eq.lrj) cycle
              lmk=lsm_inn(lrk)
!t1t1(12-3)  drl(33)-
!t1t1(12-3)  drl(33)-c"(11)-
!t1t1(12-3)  drl(33)-c"(11)-c"(11)-
              call drl_dd_ext(lrk)
            enddo
          enddo
        enddo

      return
      end

      subroutine ttdd_ar_act_bl_dd_ext_sgt1(lra)
!t1d1(15-1)  (11)ar(13)-
!t1d1(15-1)  ar(13)-c'(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      jmlr=mul_tab(jml,jmr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jmlr) cycle
        w0td1=w0_t1d1(1)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1)w0td1=-w0td1

        do lrk=norb_frz+1,lri-1
          lmk=lsm_inn(lrk)
          if(lmk.eq.jmr) then
!t1d1(15-1)  (11)ar(13)-
            iwdl=just(lrk,lri)    !
            iwdr=jud(lrk)
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0td1
              vplp_w1(mpl)=vplpnew_w1(mpl)*w0td1
            enddo
          call ar_bl_dd_ext(lri,lra,1)
          endif
        enddo
        do lrk=lri+1,norb_dz
          lmk=lsm_inn(lrk)
          if(lmk.eq.jmr) then
!t1d1(15-1)  ar(13)-c'(11)-
            iwdl=just(lri,lrk)    !
            iwdr=jud(lrk)
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            do mpl=1,mtype
              vplp_w0(mpl)=-vplpnew_w0(mpl)*w0td1
              vplp_w1(mpl)=-vplpnew_w1(mpl)*w0td1
            enddo
            call ar_bl_dd_ext(lri,lra,1)
          endif
        enddo
      enddo

      return
      end

      subroutine dddd_arbl_act_c_dd_ext_sgt0()
!d1d1(20-1) ar(13)-bl(31)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      jmlr=mul_tab(jml,jmr)
      do lri=norb_frz+1,norb_dz-1
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jmlr) cycle
          w0dd1=w0_d1d1(1)
          w1dd1=w1_d1d1(1)
          ni=mod(lrj-lri,2)
          if(ni.eq.0) then
            w0dd1=-w0dd1
            w1dd1=-w1dd1
          endif
          if(lmi.eq.jml.and.lmj.eq.jmr) then
!d1d1(20-1) ar(13)-bl(31)-
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0dd1
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1dd1
            enddo
            iwdl=jud(lri)
            iwdr=jud(lrj)
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            call ar_bl_dd_ext(lri,lrj,1)
          endif
        enddo
      enddo

      return
      end

      subroutine dddd_drl_act_c_dd_ext_sgt0()
!d1d1(20-2) drl(11)-
!d1d1(20-3) drl(33)-
!d1d1(20-3) drl(33)-c"(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      if(jmr.ne.jml) return
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jml) cycle
!d1d1(20-2) drl(11)-
        iwdl=jud(lri)
        iwdr=iwdl
        do mpl=1,mhlp
          iwal=lpnew_lwei(mpl)
          iwar=lpnew_rwei(mpl)
          lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
        enddo
        do mpl=1,mtype
          vplp_w0(mpl)=vplpnew_w0(mpl)*w0_d1d1(2)
          vplp_w1(mpl)=vplpnew_w1(mpl)*w1_d1d1(2)
        enddo
        call drl_dd_ext(lri)
        do mpl=1,mtype
          vplp_w0(mpl)=vplpnew_w0(mpl)*w0_d1d1(3)
          vplp_w1(mpl)=0.d0
        enddo
!d1d1(20-3) drl(33)-
!d1d1(20-3) drl(33)-c"(11)-
        do lrk=1,norb_dz
          if(lrk.eq.lri) cycle
          call drl_dd_ext(lrk)
        enddo
      enddo

      return
      end

      subroutine dd1_arbl_act_c_dd_ext_sgt0()
!dd1(21)ar(23)-bl(31)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      jmlr=mul_tab(jml,jmr)
      do lri=norb_frz+1,norb_dz-1
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jmlr) cycle
          w0dd1=w0_dd1
          w1dd1=w1_dd1
          ni=mod(lrj-lri,2)
          if(ni.eq.0) then
            w0dd1=-w0dd1
            w1dd1=-w1dd1
          endif
          if(lmi.eq.jml.and.lmj.eq.jmr) then
!dd1(21)ar(23)-bl(31)-
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0dd1
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1dd1
            enddo
            iwdl=jud(lri)
            iwdr=jud(lrj)
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            call ar_bl_dd_ext(lri,lrj,1)
          endif
        enddo
      enddo

      return
      end

      subroutine d1d_arbl_act_c_dd_ext_sgt0()
!d1d(22-1)   ar(13)-bl(32)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      jmlr=mul_tab(jml,jmr)
      do lri=norb_frz+1,norb_dz-1
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jmlr) cycle
          w0dd1=w0_d1d(1)
          w1dd1=w1_d1d(1)
          ni=mod(lrj-lri,2)
          if(ni.eq.0) then
            w0dd1=-w0dd1
            w1dd1=-w1dd1
          endif
          if(lmi.eq.jml.and.lmj.eq.jmr) then
!d1d(22-1)   ar(13)-bl(32)-
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0dd1
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1dd1
            enddo
            iwdl=jud(lri)
            iwdr=jud(lrj)
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            call ar_bl_dd_ext(lri,lrj,1)
          endif
        enddo
      enddo

      return
      end

      subroutine d1d_drl_act_c_dd_ext_sgt0()
!d1d(22-2)   drl(12)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      if(jml.ne.jmr) return
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jml) cycle
!d1d(22-2)   drl(12)-
        do mpl=1,mtype
          vplp_w0(mpl)=0.d0
          vplp_w1(mpl)=vplpnew_w1(mpl)*w1_d1d(2)
        enddo
        iwdl=jud(lri)
        iwdr=iwdl
        do mpl=1,mhlp
          iwal=lpnew_lwei(mpl)
          iwar=lpnew_rwei(mpl)
          lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
        enddo
        call drl_dd_ext(lri)
      enddo

      return
      end

      subroutine d1v_ar_act_bl_dd_ext_sgt0(lra)
!d1v(24-1)  ar(13)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      jmlr=mul_tab(jml,jmr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jmlr) cycle
        w0dv1=w0_d1v(1)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1)w0dv1=-w0dv1
!d1v(24-1)  ar(13)-
        iwdl=jud(lri)
        iwdr=0
        do mpl=1,mhlp
          iwal=lpnew_lwei(mpl)
          iwar=lpnew_rwei(mpl)
          lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
        enddo
        do mpl=1,mtype
          vplp_w0(mpl)=vplpnew_w0(mpl)*w0dv1
          vplp_w1(mpl)=vplpnew_w1(mpl)*w0dv1
        enddo
        call ar_bl_dd_ext(lri,lra,1)
      enddo

      return
      end

      subroutine sd_adb_act_c_ext_ar(lin)
!sd(6-3) a&r(13)c'(22)-
!sd(6-10) d&r&l(12)b^l(23)
!sd(6-15) d&r&l(33)b^l(13)c'(22)
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
          if(lmij.ne.jml) cycle
         if(lmi.eq.jmr) then
            iwdl=just(lrj,lri)
            iwdr=jud(lri)
            w0sd10=w0_sd(10)
            w1sd10=w1_sd(10)
            ni=mod(norb_dz-lrj,2)
            if(ni.eq.1) then
              w0sd10=-w0sd10
              w1sd10=-w1sd10
            endif
!sd(6-10) d&r&l(12)b^l(23)
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd10
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1sd10
            enddo
            call drl_bl_ext_ar_new(lin,lri,lrj)
          endif
!sd(6-3) a&r(13)c'(22)-
          w0sd3=w0_sd(3)
          w1sd3=w1_sd(3)
          w0sd15=w0_sd(15)
          w1sd15=w1_sd(15)
          ni=mod(norb_dz-lri,2)
          if(ni.eq.0) then
            w0sd3=-w0sd3
            w1sd3=-w1sd3
            w0sd15=-w0sd15
            w1sd15=-w1sd15
          endif
          if(lmj.eq.jmr) then
            iwdl=just(lrj,lri)
            iwdr=jud(lrj)
           do mpl=1,mtype
             vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd3
             vplp_w1(mpl)=vplpnew_w1(mpl)*w1sd3
           enddo
           do mpl=1,mhlp
             iwal=lpnew_lwei(mpl)
             iwar=lpnew_rwei(mpl)
             lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
             lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
           enddo
      if(lin.eq.13) then
        call ar_td_ext_ar(57,lri,lrj,isma)
        call ar_td_ext_rest(lri)
      endif
      if(lin.eq.6) then
        call ar_sd_ext_ar(57,lri,lrj,isma)
        call ar_sd_ext_rest(lri)
      endif
      if(lin.eq.23) then
        call ar_dv_ext_ar(57,isma,lri,lrj) !ar_dv
      endif
!sd(6-15) d&r&l(33)b^l(13)c'(22)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd15
              vplp_w1(mpl)=0.d0
            enddo
            do lrk=1,lri-1
              call drl_bl_ext_ar_new(lin,lrk,lri)
            enddo
          endif
!sd(6-6) a&r(13)b&r(23)b^r(32)
          do lrd=lrj+1,norb_dz
            lmd=lsm_inn(lrd)
           if(lmd.ne.jmr) cycle

            w0sd6=w0_sd(6)
            w1sd6=w1_sd(6)
            ni=mod(lrj-lri+norb_dz-lrd,2)
            if(ni.eq.0) then
              w0sd6=-w0sd6
              w1sd6=-w1sd6
            endif
            iwdl=just(lrj,lri)
            iwdr=jud(lrd)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd6
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1sd6
            enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            ijk=lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lrd-norb_frz)
            intpos=intind_ijka(ijk)
            call ar_br_br_ext_ar_new(lin,intpos,isma)
          enddo
!sd(6-7) a&r(13)b&l(32)b^l(23)
          do lrd=lri+1,lrj-1
            lmd=lsm_inn(lrd)
            if(lmd.ne.jmr) cycle
            iwdl=just(lrj,lri)
            iwdr=jud(lrd)
            w0sd7=w0_sd(7)
            w1sd7=w1_sd(7)
            ni=mod(lrj-lri+norb_dz-lrd,2)
            if(ni.eq.0) then
              w0sd7=-w0sd7
              w1sd7=-w1sd7
            endif
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd7
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1sd7
            enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            ijk=lri-norb_frz+ngw2(lrd-norb_frz)+ngw3(lrj-norb_frz)
            intpos=intind_ijka(ijk)
            call ar_bl_bl_ext_ar_new(lin,intpos,isma,1)
          enddo
        enddo
      enddo

      return
      end

      subroutine sd_ar_act_brbr_sgt0(lin,jk)
!sd(6-3) a&r(13)c'(22)-
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
        w0sd3 =w0_sd(3)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1) then
          w0sd3 =-w0sd3
        endif
!-------------------------------------------------------------------
!sd(6-3) a&r(13)c'(22)-
        if(jb_sys.gt.0) then
          do mpl=1,mtype
            vplp_w0(mpl)=-vplpnew_w0(mpl)*w0sd3
            vplp_w1(mpl)=-vplpnew_w1(mpl)*w0sd3
          enddo
          do lrk=lri+1,norb_dz
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
        endif
      enddo

      return
      end

      subroutine sd_ar_act_blbl_sgt0(lin,jk)
!sd(6-3) a&r(13)c'(22)-
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
        w0sd3 =w0_sd(3)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.0) then
          w0sd3 =-w0sd3
        endif
!sd(6-3) a&r(13)c'(22)-
         do mpl=1,mtype
           vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd3
           vplp_w1(mpl)=vplpnew_w1(mpl)*w0sd3
         enddo
         do lrk=lri+1,norb_dz
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
      enddo

      return
      end

      subroutine sv_arbr_act_br_sgt0(lin,lra)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      isma=mul_tab(iml,imr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        do lrj=lri,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jml) cycle
          w0sv1=w0_sv(1)
          w1sv1=w1_sv(1)
          ni=mod(lrj-lri,2)
          if(ni.eq.0) then
            w0sv1=-w0sv1
            w1sv1=-w1sv1
          endif
!-------------------------------------------------------------------
         if(lri.ne.lrj) then
!sv(10-1) ar(13)-br(23)-
            if(jb_sys.gt.0) then
              iwdl=just(lrj,lri)
              iwdr=0
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
              do mpl=1,mtype
                vplp_w0(mpl)=vplpnew_w0(mpl)*w0sv1
                vplp_w1(mpl)=vplpnew_w1(mpl)*w1sv1
              enddo
              ijk=lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz)
             intpos=intind_ijka(ijk)
              call ar_br_br_ext_ar_new(lin,intpos,isma)
            endif
          endif
        enddo
      enddo

      return
      end


      subroutine d1v_ar_act_c_ext(lin)
!d1v(24-1)  ar(13)-
!d1v(24-2)  drl(33)-bl(13)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      isma=mul_tab(iml,imr)
      do lrd=norb_frz+1,norb_dz
        lmd=lsm_inn(lrd)
        if(lmd.ne.jml) cycle
        ni=mod(norb_dz-lrd,2)

!d1v(24-1)  ar(13)-
        w0=w0_d1v(1)
        if(ni.eq.1)w0=-w0_d1v(1)
          iwdl=jud(lrd)
          iwdr=0
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0
          enddo
        do mpl=1,mhlp
          iwal=lpnew_lwei(mpl)
          iwar=lpnew_rwei(mpl)
          lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
        enddo
        if(lin.eq.6) then
          call ar_sd_ext_ar(51,lrd,0,isma)
          call ar_sd_ext_rest(lrd)
        endif
        if(lin.eq.13) then
          call ar_td_ext_ar(51,lrd,0,isma)
          call ar_td_ext_rest(lrd)
        endif
!d1v(24-2)  drl(33)-bl(13)-
        w0=w0_d1v(2)
        if(ni.eq.1)w0=-w0_d1v(2)
        do mpl=1,mtype
          vplp_w0(mpl)=vplpnew_w0(mpl)*w0
          vplp_w1(mpl)=0.d0
        enddo
        do lrk=1,lrd-1
          call drl_bl_ext_ar_new(lin,lrk,lrd)
        enddo
      enddo

      return
      end

      subroutine d1v_ar_act_brbr_ext(lin,jk)
!d1v(24-1)  ar(13)-
!d1v(24-2)  drl(33)-bl(13)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      isma=mul_tab(iml,imr)
      do lrd=norb_frz+1,norb_dz
        lmd=lsm_inn(lrd)
        if(lmd.ne.jml) cycle
        ni=mod(norb_dz-lrd,2)
        w0=w0_d1v(1)
        if(ni.eq.1)w0=-w0_d1v(1)
!..............................................................
!d1v(24-1) ar(13)-
        iwdl=jud(lrd)
        iwdr=0
        do mpl=1,mtype
          vplp_w0(mpl)=vplpnew_w0(mpl)*w0
          vplp_w1(mpl)=vplpnew_w1(mpl)*w0
        enddo
        do mpl=1,mhlp
          iwal=lpnew_lwei(mpl)
          iwar=lpnew_rwei(mpl)
          lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
        enddo
        ijk=lrd-norb_frz+jk
        intpos=intind_ijka(ijk)
        call ar_br_br_ext_ar_new(lin,intpos,isma)
      enddo

      return
      end

      subroutine d1v_ar_act_blbl_ext(lin,jk)
!d1v(24-1)  ar(13)-
!d1v(24-2)  drl(33)-bl(13)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      isma=mul_tab(iml,imr)
      do lrd=norb_frz+1,norb_dz
        lmd=lsm_inn(lrd)
        if(lmd.ne.jml) cycle
        ni=mod(norb_dz-lrd,2)
        w0=w0_d1v(1)
        if(ni.eq.1)w0=-w0_d1v(1)
!dv(23-1) ar(23)-
        iwdl=jud(lrd)
        iwdr=0
        do mpl=1,mtype
          vplp_w0(mpl)=vplpnew_w0(mpl)*w0
          vplp_w1(mpl)=vplpnew_w1(mpl)*w0
        enddo
        do mpl=1,mhlp
          iwal=lpnew_lwei(mpl)
          iwar=lpnew_rwei(mpl)
          lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
        enddo
        ijk=lrd-norb_frz+jk
        intpos=intind_ijka(ijk)
        call ar_bl_bl_ext_ar_new(lin,intpos,isma,1)
      enddo

      return
      end

      subroutine ts_arbl_act_c_ext_ab_sgt0(lin)
!ts(3-3) ar(23)-bl(31)-c"(22)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      do lri=norb_frz+1,norb_dz-1
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          w1ts3=w1_ts(3)
          ni=mod(lrj-lri,2)
          if(ni.eq.0) then
            w1ts3=-w1ts3
          endif
          do lrk=lrj+1,norb_dz
            lmk=lsm_inn(lrk)
            lmki=mul_tab(lmk,lmi)
            lmkj=mul_tab(lmk,lmj)
            if(lmki.eq.jml.and.lmkj.eq.jmr) then
!ts(3-3) ar(23)-bl(31)-c"(22)-
              iwdl=just(lri,lrk)
              iwdr=just(lrk,lrj)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
              do mpl=1,mtype
                vplp_w0(mpl)=0.d0
                vplp_w1(mpl)=vplpnew_w1(mpl)*w1ts3
              enddo
              call arbl_act_c_link_ext_ab(lin,lri,lrj)
            endif
          enddo
        enddo
      enddo

      return
      end

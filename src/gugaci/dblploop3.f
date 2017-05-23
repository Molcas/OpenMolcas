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
!===================================================
!  have computed vlop0 and vlop1, link to extern space
!  (ss, st, ts, tt)
!  lin type of the act-ext space node
!  lri index of ar
!  lrj index of bl
!====================================================
      subroutine arbl_act_c_link_ext_ab(lin,lri,lrj)

      if(lin.eq.1) call ar_bl_ext_ss(lri,lrj,1)
      if(lin.eq.2) call ar_bl_ext_st(lri,lrj,1)
      if(lin.eq.3) call ar_bl_ext_ts(lri,lrj,1)
      if(lin.eq.11) call ar_bl_ext_tt(lri,lrj,1)
      if(lin.eq.10) call ar_br_sv_ext_br_ar(lri,lrj)
      if(lin.eq.17) call ar_br_tv_ext_br_ar(lri,lrj)

      return
      end

!===================================================
!  have computed vlop0 and vlop1, link to extern space
!  (ss, st, ts, tt)
!  lin type of the act-ext space node
!  lri index of drl
!====================================================
      subroutine drl_act_c_link_ext_ab(lin,lri)

      if(lin.eq.1) call drl_ss_ext(lri)
      if(lin.eq.2) call drl_st_ext(lri)
      if(lin.eq.3) call drl_ts_ext(lri)
      if(lin.eq.11) call drl_tt_ext(lri)

      return
      end

      subroutine drl_act_c_link_ext_ab_sum(lin,lri,lrj)

      if(lin.eq.1) call drl_ss_sum(lri,lrj)
      if(lin.eq.11) call drl_tt_sum(lri,lrj)

      return
      end

      subroutine sd_ar_act_bl_sgt0(lin,lra)
!sd(6-3) a&r(13)c'(22)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      jmlr=mul_tab(jml,jmr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jmlr) cycle
        w0sd3 =w0_sd(3)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1) then
          w0sd3 =-w0sd3
        endif
!sd(6-3) a&r(13)c'(22)-
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
          if(lin.eq.1) call ar_bl_ext_ss(lri,lra,1)
          if(lin.eq.2) call ar_bl_ext_st(lri,lra,1)
          if(lin.eq.3) call ar_bl_ext_ts(lri,lra,1)
          if(lin.eq.11) call ar_bl_ext_tt(lri,lra,1)
        enddo
      enddo

      return
      end

!======================================================
!extern space -bl-ar or -br-al (ss,st,ts,tt)
!double occupied space ar-bl-  (ss)
!active space -c"-
!spin great than 0
!======================================================
      subroutine ss_arbl_act_c_ext_ab_sgt0(lin)
!-------------------------------------------------------------------
!ss(1-1)  ar(01)-bl(32)-        act -c"-
!ss(1-3)  ar(13)-bl(20)-        act -c"-
!ss(1-6)  (11)-ar(23)-bl(32)-   act -c"-
!ss(1-7)  ar(13)-c'(21)-bl(32)- act -c"-
!ss(1-8)  ar(13)-c'(22)-bl(31)- act -c"-
!ss(1-9)  ar(23)-c'(11)-bl(32)- act -c"-
!ss(1-11) ar(13)-bl(31)-c"(22)- act -c"-
!ss(1-12) ar(13)-bl(32)-c"(21)- act -c"-
!ss(1-13) ar(23)-bl(31)-c"(12)- act -c"-
!-------------------------------------------------------------------
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
          if(jml.eq.1.and.lmij.eq.jmr) then
!ss(1-1)  ar(01)-bl(32)-        act -c"-
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
            call arbl_act_c_link_ext_ab(lin,lri,lrj)
          endif
          if(jmr.eq.1.and.lmij.eq.jml) then
!ss(1-3)  ar(13)-bl(20)-        act -c"-
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
            call arbl_act_c_link_ext_ab(lin,lri,lrj)
          endif
!ss(1-6)  (11)-ar(23)-bl(32)-   act -c"-
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0ss6
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1ss6
          enddo
          do lrk=norb_frz+1,lri-1
            lmk=lsm_inn(lrk)
            lmki=mul_tab(lmi,lmk)
            lmkj=mul_tab(lmj,lmk)
            if(lmki.eq.jml.and.lmkj.eq.jmr) then
              iwdl=just(lri,lrk)
              iwdr=just(lrj,lrk)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
              call arbl_act_c_link_ext_ab(lin,lri,lrj)
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
            lmkj=mul_tab(lmj,lmk)
            if(lmki.eq.jml.and.lmkj.eq.jmr) then
              iwdl=just(lrk,lri)
              iwdr=just(lrj,lrk)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
              call arbl_act_c_link_ext_ab(lin,lri,lrj)
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
              call arbl_act_c_link_ext_ab(lin,lri,lrj)
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
            call arbl_act_c_link_ext_ab(lin,lri,lrj)
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
            call arbl_act_c_link_ext_ab(lin,lri,lrj)
            endif
          enddo
!ss(1-12) ar(13)-bl(32)-c"(21)- act -c"-
          do mpl=1,mtype
            vplp_w0(mpl)=0.d0
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1ss12
          enddo
          do lrk=lrj+1,norb_dz
            lmk=lsm_inn(lrk)
            lmki=mul_tab(lmk,lmi)
            lmkj=mul_tab(lmj,lmk)
            if(lmki.eq.jml.and.lmkj.eq.jmr) then
              iwdl=just(lrk,lri)
              iwdr=just(lrj,lrk)
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
              call arbl_act_c_link_ext_ab(lin,lri,lrj)
            endif
          enddo
!ss(1-13) ar(23)-bl(31)-c"(12)- act -c"-
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0ss13
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1ss13
          enddo
          do lrk=lrj+1,norb_dz
            lmk=lsm_inn(lrk)
            lmki=mul_tab(lmi,lmk)
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
              call arbl_act_c_link_ext_ab(lin,lri,lrj)
            endif
          enddo
        enddo
      enddo

      return
      end


!======================================================
!extern space -bl-ar or -br-al (ss,st,ts,tt)
!double occupied space drl-  (ss)
!active space -c"-
!spin great than 0
!======================================================
      subroutine ss_drl_act_c_ext_ab_sgt0(lin)
!ss(1-16) (11)-drl(22)-         act -c"-
!ss(1-18) drl(11)-c"(22)-       act -c"-
!ss(1-19) drl(12)-c"(21)-       act -c"-
!ss(1-20) (11)-drl(33)-c"(22)-  act -c"-
!ss(1-20) drl(33)-c"(11)-c"(22)-act -c"-
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
          call drl_act_c_link_ext_ab(lin,lrj)
!ss(1-18) drl(11)-c"(22)-       act -c"-
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0_ss(18)
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1_ss(18)
          enddo
          call drl_act_c_link_ext_ab(lin,lri)
!ss(1-20) (11)-drl(33)-c"(22)-  act -c"-
!ss(1-20) drl(33)-c"(11)-c"(22)-act -c"-
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0_ss(20)
            vplp_w1(mpl)=0.d0
          enddo
          do lrk=1,norb_dz
            if(lrk.eq.lri) cycle
            if(lrk.eq.lrj) cycle
            call drl_act_c_link_ext_ab(lin,lrk)
          enddo
        enddo
      enddo

      return
      end

!======================================================
!extern space -bl-ar or -br-al (ss,st,ts,tt)
!double occupied space drl-  (st)
!active space -c"-
!spin great than 0
!======================================================
      subroutine st_drl_act_c_ext_ab_sgt0(lin)
!st(2-7) drl(12)-c"(22)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      if(jmr.ne.jml) return
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
      do lrj=lri+1,norb_dz
        lmj=lsm_inn(lrj)
        lmij=mul_tab(lmi,lmj)
        if(jml.ne.lmij) cycle
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
        call drl_act_c_link_ext_ab(lin,lri)
      enddo
      enddo

      return
      end

!======================================================
!extern space -bl-ar or -br-al (ss,st,ts,tt)
!double occupied space ar-bl-  (st)
!active space -c"-
!spin great than 0
!======================================================
      subroutine st_arbl_act_c_ext_ab_sgt0(lin)
!st(2-3) ar(13)-c'(22)-bl(32)-
!st(2-3) ar(13)-bl(32)-c'(22)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      do lri=norb_frz+1,norb_dz-1
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
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
!st(2-3) ar(13)-c'(22)-bl(32)-
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
              call arbl_act_c_link_ext_ab(lin,lri,lrj)
            endif
          enddo

          do lrk=lrj+1,norb_dz
            lmk=lsm_inn(lrk)
            lmki=mul_tab(lmk,lmi)
            lmkj=mul_tab(lmk,lmj)
            if(lmki.eq.jml.and.lmkj.eq.jmr) then
!st(2-3) ar(13)-bl(32)-c'(22)-
              iwdl=just(lrk,lri)
              iwdr=just(lrj,lrk)
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
              call arbl_act_c_link_ext_ab(lin,lri,lrj)
            endif
          enddo
        enddo
      enddo

      return
      end

      subroutine stt_arbl_act_c_ext_ab_sgt1(lin)
!st1(4-1) ar(01)-bl(31)-
!st1(4-2) (11)ar(23)-bl(31)-
!st1(4-3) ar(13)-c'(21)-bl(31)-
!st1(4-3) ar(13)-bl(31)-c"(21)-
!st1(4-4) ar(23)-c'(11)-bl(31)-
!st1(4-4) ar(23)-bl(31)-c"(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
      do lrj=lri+1,norb_dz
        lmj=lsm_inn(lrj)
        w1st1=w1_st1(1)
        w1st2=w1_st1(2)
        w1st3=w1_st1(3)
        w1st4=w1_st1(4)
        if(mod(lrj-lri,2).eq.0) then
          w1st1=-w1st1
          w1st2=-w1st2
          w1st3=-w1st3
          w1st4=-w1st4
        endif
!st1(4-1) ar(01)-bl(31)-
        if(jml.eq.1.and.jmr.eq.mul_tab(lmi,lmj)) then
          iwdl=just(lri,lri)
          iwdr=just(lri,lrj)
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
          call arbl_act_c_link_ext_ab(lin,lri,lrj)
        endif
!st1(4-2) (11)ar(23)-bl(31)-
        do lrk=norb_frz+1,lri-1
           lmk=lsm_inn(lrk)
      if(jml.ne.mul_tab(lmk,lmi).or.jmr.ne.mul_tab(lmk,lmj)) cycle
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
             vplp_w1(mpl)=vplpnew_w1(mpl)*w1st2
           enddo
           call arbl_act_c_link_ext_ab(lin,lri,lrj)
        enddo
!st1(4-3) ar(13)-c'(21)-bl(31)-
         do lrk=lri+1,lrj-1
           lmk=lsm_inn(lrk)
      if(jml.eq.mul_tab(lmi,lmk).and.jmr.eq.mul_tab(lmk,lmj)) then
           iwdl=just(lrk,lri)
           iwdr=just(lrk,lrj)
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
           call arbl_act_c_link_ext_ab(lin,lri,lrj)
!st1(4-4) ar(23)-c'(11)-bl(31)-
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
             vplp_w1(mpl)=-vplpnew_w1(mpl)*w1st4
           enddo
           call arbl_act_c_link_ext_ab(lin,lri,lrj)
      endif
      enddo
      do lrk=lrj+1,norb_dz
         lmk=lsm_inn(lrk)
!st1(4-4) ar(23)-bl(31)-c"(11)-
      if(jml.eq.mul_tab(lmi,lmk).and.jmr.eq.mul_tab(lmj,lmk)) then
           iwdl=just(lri,lrk)
           iwdr=just(lrj,lrk)
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
           call arbl_act_c_link_ext_ab(lin,lri,lrj)
!st1(4-3) ar(13)-bl(31)-c"(21)-
           iwdl=just(lrk,lri)
           iwdr=just(lrj,lrk)
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
           call arbl_act_c_link_ext_ab(lin,lri,lrj)
      endif
         enddo
      enddo
      enddo

      return
      end

      subroutine tts_arbl_act_c_ext_ab_sgt1(lin)
!t1s(5-1)   ar(13)-bl(10)-
!t1s(5-2)   ar(13)-bl(32)-
!t1s(5-2)   ar(13)-c'(11)-bl(32)-
!t1s(5-3)   ar(13)-bl(31)-c"(12)-
!t1s(5-4)   ar(13)-bl(32)-c"(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
      do lrj=lri+1,norb_dz
        lmj=lsm_inn(lrj)
        w1ts1=w1_t1s(1)
        w1ts2=w1_t1s(2)
        w1ts3=w1_t1s(3)
        w1ts4=w1_t1s(4)
        if(mod(lrj-lri,2).eq.0) then
          w1ts1=-w1ts1
          w1ts2=-w1ts2
          w1ts3=-w1ts3
          w1ts4=-w1ts4
        endif
!t1s(5-1)   ar(13)-bl(10)-
        if(jml.eq.mul_tab(lmi,lmj).and.jmr.eq.1) then
          iwdl=just(lri,lrj)
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
          call arbl_act_c_link_ext_ab(lin,lri,lrj)
        endif
!t1s(5-2)   (11)ar(13)-bl(32)-
        do mpl=1,mtype
          vplp_w0(mpl)=0.d0
          vplp_w1(mpl)=vplpnew_w1(mpl)*w1ts2
        enddo
        do lrk=norb_frz+1,lri-1
          lmk=lsm_inn(lrk)
      if(jml.ne.mul_tab(lmk,lmi).or.jmr.ne.mul_tab(lmk,lmj)) cycle
          iwdl=just(lrk,lri)
          iwdr=just(lrj,lrk)
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo
          call arbl_act_c_link_ext_ab(lin,lri,lrj)
        enddo
!t1s(5-2)   ar(13)-c'(11)-bl(32)-
        do mpl=1,mtype
          vplp_w0(mpl)=0.d0
          vplp_w1(mpl)=-vplpnew_w1(mpl)*w1ts2
        enddo
        do lrk=lri+1,lrj-1
           lmk=lsm_inn(lrk)
      if(jml.ne.mul_tab(lmi,lmk).or.jmr.ne.mul_tab(lmk,lmj)) cycle
          iwdl=just(lri,lrk)
          iwdr=just(lrj,lrk)
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo
          call arbl_act_c_link_ext_ab(lin,lri,lrj)
        enddo
!t1s(5-3)   ar(13)-bl(31)-c"(12)-
        do mpl=1,mtype
          vplp_w0(mpl)=0.d0
          vplp_w1(mpl)=vplpnew_w1(mpl)*w1ts3
        enddo
        do lrk=lrj+1,norb_dz
          lmk=lsm_inn(lrk)
      if(jml.ne.mul_tab(lmi,lmk).or.jmr.ne.mul_tab(lmj,lmk)) cycle
          iwdl=just(lri,lrk)
          iwdr=just(lrk,lrj)
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo
          call arbl_act_c_link_ext_ab(lin,lri,lrj)
        enddo
!t1s(5-4)   ar(13)-bl(32)-c"(11)-
        do mpl=1,mtype
          vplp_w0(mpl)=0.d0
          vplp_w1(mpl)=vplpnew_w1(mpl)*w1ts4
        enddo
        do lrk=lrj+1,norb_dz
          lmk=lsm_inn(lrk)
      if(jml.ne.mul_tab(lmi,lmk).or.jmr.ne.mul_tab(lmj,lmk)) cycle
          iwdl=just(lri,lrk)
          iwdr=just(lrj,lrk)
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo
          call arbl_act_c_link_ext_ab(lin,lri,lrj)
        enddo
      enddo
      enddo

      return
      end

      subroutine tts_drl_act_c_ext_ab_sgt1(lin)
!t1s(5-5)   (11)drl(12)-
!t1s(5-6)   drl(11)-c"(12)-
!t1s(5-7)   drl(12)-c"(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      if(jml.ne.jmr) return
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
      do lrj=lri+1,norb_dz
        lmj=lsm_inn(lrj)
        lmij=mul_tab(lmi,lmj)
        if(jml.ne.lmij) cycle
!t1s(5-5)   (11)drl(12)-
        iwdl=just(lri,lrj)
        iwdr=just(lrj,lri)      !
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
        call drl_act_c_link_ext_ab(lin,lrj)
!t1s(5-6)   drl(11)-c"(12)-
        iwdl=just(lri,lrj)
        iwdr=just(lrj,lri)      !
        do mpl=1,mhlp
          iwal=lpnew_lwei(mpl)
          iwar=lpnew_rwei(mpl)
          lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
        enddo
        do mpl=1,mtype
          vplp_w0(mpl)=0.d0
          vplp_w1(mpl)=vplpnew_w1(mpl)*w1_t1s(6)
        enddo
        call drl_act_c_link_ext_ab(lin,lri)
!t1s(5-7)   drl(12)-c"(11)-
        iwdl=just(lri,lrj)
        iwdr=just(lri,lrj)      !
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
        call drl_act_c_link_ext_ab(lin,lri)
      enddo
      enddo

      return
      end

      subroutine sdd_ar_act_bl_sgt0(lin,lra)
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
          if(lin.eq.1) call ar_bl_ext_ss(lri,lra,1)
          if(lin.eq.2) call ar_bl_ext_st(lri,lra,1)
          if(lin.eq.3) call ar_bl_ext_ts(lri,lra,1)
          if(lin.eq.11) call ar_bl_ext_tt(lri,lra,1)
        endif
!sd1(8-2)    (11)ar(23)-
        do mpl=1,mtype
          vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd2
          vplp_w1(mpl)=vplpnew_w1(mpl)*w0sd2
        enddo
        do lrk=norb_frz+1,lri-1
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
          if(lin.eq.1) call ar_bl_ext_ss(lri,lra,1)
          if(lin.eq.2) call ar_bl_ext_st(lri,lra,1)
          if(lin.eq.3) call ar_bl_ext_ts(lri,lra,1)
          if(lin.eq.11) call ar_bl_ext_tt(lri,lra,1)

        enddo
!sd1(8-3)    ar(13)-c'(21)-
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
          if(lin.eq.1) call ar_bl_ext_ss(lri,lra,1)
          if(lin.eq.2) call ar_bl_ext_st(lri,lra,1)
          if(lin.eq.3) call ar_bl_ext_ts(lri,lra,1)
          if(lin.eq.11) call ar_bl_ext_tt(lri,lra,1)
        enddo
!sd1(8-4)    ar(23)-c'(11)-
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
          if(lin.eq.1) call ar_bl_ext_ss(lri,lra,1)
          if(lin.eq.2) call ar_bl_ext_st(lri,lra,1)
          if(lin.eq.3) call ar_bl_ext_ts(lri,lra,1)
          if(lin.eq.11) call ar_bl_ext_tt(lri,lra,1)
        enddo
      enddo
      return
      end

      subroutine tttt_arbl_act_c_ext_ab_sgt0(lin)
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
          lmij=mul_tab(lmi,lmj)
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
              iwdr=just(lrk,lrj)      !
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
              call arbl_act_c_link_ext_ab(lin,lri,lrj)
            endif
          enddo
!t1t1(12-1)  ar(13)-bl(31)-c"(11)-
          do lrk=lrj+1,norb_dz
            lmk=lsm_inn(lrk)
            lmki=mul_tab(lmk,lmi)
            lmkj=mul_tab(lmk,lmj)
            if(lmki.eq.jml.and.lmkj.eq.jmr) then
              iwdl=just(lri,lrk)     !
              iwdr=just(lrj,lrk)      !
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
              call arbl_act_c_link_ext_ab(lin,lri,lrj)
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
              iwdl=just(lri,lrk)     !
              iwdr=just(lrk,lrj)     !
              do mpl=1,mhlp
                iwal=lpnew_lwei(mpl)
                iwar=lpnew_rwei(mpl)
                lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
                lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
              enddo
              call arbl_act_c_link_ext_ab(lin,lri,lrj)
            endif
          enddo
        enddo
      enddo

      return
      end


      subroutine tttt_drl_act_c_ext_ab_sgt0(lin)
!t1t1(12-2)  drl(11)-
!t1t1(12-2)  drl(11)-c"(11)-
!t1t1(12-3)  drl(33)-
!t1t1(12-3)  drl(33)-c"(11)-
!t1t1(12-3)  drl(33)-c"(11)-c"(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      w0tt2=w0_t1t1(2)
      w1tt2=w1_t1t1(2)
      w0tt3=w0_t1t1(3)
      do lri=norb_frz+1,norb_dz-1
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jml.or.lmij.ne.jmr) cycle
!t1t1(12-2)  drl(11)-
!t1t1(12-2)  drl(11)-c"(11)-
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0tt2
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1tt2
            enddo
            iwdl=just(lri,lrj)   !
            iwdr=iwdl
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            call drl_act_c_link_ext_ab(lin,lri)
            call drl_act_c_link_ext_ab(lin,lrj)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0tt3
              vplp_w1(mpl)=0.d0
            enddo
!t1t1(12-3)  drl(33)-
!t1t1(12-3)  drl(33)-c"(11)-
!t1t1(12-3)  drl(33)-c"(11)-c"(11)-
c      do lrk=1,norb_dz
c        if(lrk.eq.lri) cycle
c        if(lrk.eq.lrj) cycle
c        call drl_ss_ext(lrk)
c     enddo
            call drl_act_c_link_ext_ab_sum(lin,lri,lrj)
          enddo
        enddo

      return
      end

      subroutine ttdd_ar_act_bl_sgt1(lin,lra)
!t1d1(15-1)  ar(13)-
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
        if(ni.eq.1) w0td1=-w0td1
!-------------------------------------------------------------------
!t1d1(15-1)  (11)ar(13)-
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
            if(lin.eq.1) call ar_bl_ext_ss(lri,lra,1)
            if(lin.eq.2) call ar_bl_ext_st(lri,lra,1)
            if(lin.eq.3) call ar_bl_ext_ts(lri,lra,1)
            if(lin.eq.11) call ar_bl_ext_tt(lri,lra,1)
          enddo
!-------------------------------------------------------------------
!t1d1(15-1)  ar(13)-c'(11)-
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
            if(lin.eq.1) call ar_bl_ext_ss(lri,lra,1)
            if(lin.eq.2) call ar_bl_ext_st(lri,lra,1)
            if(lin.eq.3) call ar_bl_ext_ts(lri,lra,1)
            if(lin.eq.11) call ar_bl_ext_tt(lri,lra,1)
          enddo
        enddo

      return
      end

      subroutine d1d1_arbl_act_c_ext_ab_sgt0(lin)
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
          if(jml.ne.lmi.and.jmr.ne.lmj) cycle
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
            call arbl_act_c_link_ext_ab(lin,lri,lrj)
          endif
        enddo
      enddo

      return
      end

      subroutine d1d1_drl_act_c_ext_ab_sgt0(lin)
!d1d1(20-2) drl(11)-
!d1d1(20-3) drl(33)-
!d1d1(20-3) drl(33)-c"(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      if(jml.ne.jmr) return
      w0dd2=w0_d1d1(2)
      w1dd2=w1_d1d1(2)
      w0dd3=w0_d1d1(3)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jml) cycle
!d1d1(20-2) drl(11)-
        do mpl=1,mtype
          vplp_w0(mpl)=vplpnew_w0(mpl)*w0dd2
          vplp_w1(mpl)=vplpnew_w1(mpl)*w1dd2
        enddo
        iwdl=jud(lri)
        iwdr=iwdl
        do mpl=1,mhlp
          iwal=lpnew_lwei(mpl)
          iwar=lpnew_rwei(mpl)
          lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
        enddo
        call drl_act_c_link_ext_ab(lin,lri)
        do mpl=1,mtype
          vplp_w0(mpl)=vplpnew_w0(mpl)*w0dd3
          vplp_w1(mpl)=0.d0
        enddo
!d1d1(20-3) drl(33)-
!d1d1(20-3) drl(33)-c"(11)-
      do lrk=1,norb_dz
        if(lrk.eq.lri) cycle
        call drl_act_c_link_ext_ab(lin,lrk)
      enddo
!        call drl_act_c_link_ext_ab_sum(lin,lri,0)
      enddo

      return
      end

      subroutine dd1_arbl_act_c_ext_ab_sgt0(lin)
!dd1(21) ar(23)-bl(31)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      jmlr=mul_tab(jml,jmr)
      do lri=norb_frz+1,norb_dz
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
            call arbl_act_c_link_ext_ab(lin,lri,lrj)
          endif
        enddo
      enddo

      return
      end

      subroutine d1d_arbl_act_c_ext_ab_sgt0(lin)
!d1d(22-1)   ar(13)-bl(32)-
!d1d(22-2)   drl(12)-
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
            call arbl_act_c_link_ext_ab(lin,lri,lrj)
          endif
        enddo
      enddo

      return
      end

      subroutine d1d_drl_act_c_ext_ab_sgt0(lin)
!d1d(22-2)   drl(12)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      if(jmr.ne.jml) return
      do lri=norb_frz+1,norb_dz
         lmi=lsm_inn(lri)
         if(jml.ne.lmi) cycle
         do mpl=1,mtype
           vplp_w0(mpl)=vplpnew_w0(mpl)*w0_d1d(2)
           vplp_w1(mpl)=vplpnew_w1(mpl)*w1_d1d(2)
         enddo
         iwdl=jud(lri)
         iwdr=jud(lri)
         do mpl=1,mhlp
           iwal=lpnew_lwei(mpl)
           iwar=lpnew_rwei(mpl)
           lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
           lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
         enddo
         call drl_act_c_link_ext_ab(lin,lri)
      enddo

      return
      end

      subroutine d1v_ar_act_bl_ext_ab_sgt0(lin,lra)
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
        call arbl_act_c_link_ext_ab(lin,lri,lra)
      enddo

      return
      end

      subroutine ar_br_stv_link_ext_brar(lin,lri,lrj)
        if(lin.eq.17) call ar_br_tv_ext_br_ar(lri,lrj)
        if(lin.eq.10) call ar_br_sv_ext_br_ar(lri,lrj)
      end

      subroutine sd_ar_act_br_sgt0(lin,lra)
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

!sd(6-3) a&r(13)c'(22)-
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
          call ar_br_stv_link_ext_brar(lin,lri,lra)
        enddo
      enddo

      return
      end

      subroutine sdd_ar_act_br_sgt0(lin,lra)
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
          call ar_br_stv_link_ext_brar(lin,lri,lra)
        endif
!sd1(8-2)    (11)ar(23)-
        do mpl=1,mtype
          vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd2
          vplp_w1(mpl)=vplpnew_w1(mpl)*w0sd2
        enddo
        do lrk=norb_frz+1,lri-1
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
          call ar_br_stv_link_ext_brar(lin,lri,lra)
        enddo
!sd1(8-3)    ar(13)-c'(21)-
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
          call ar_br_stv_link_ext_brar(lin,lri,lra)
        enddo
!sd1(8-4)    ar(23)-c'(11)-
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
          call ar_br_stv_link_ext_brar(lin,lri,lra)
        enddo
      enddo

      return
      end

      subroutine ttdd_ar_act_br_sgt1(lin,lra)
!t1d1(15-1)  ar(13)-
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
        if(ni.eq.1) w0td1=-w0td1
!-------------------------------------------------------------------
!t1d1(15-1)  (11)ar(13)-
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
            call ar_br_stv_link_ext_brar(lin,lri,lra)
          enddo
!-------------------------------------------------------------------
!t1d1(15-1)  ar(13)-c'(11)-
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
            call ar_br_stv_link_ext_brar(lin,lri,lra)
          enddo
        enddo

      return
      end

      subroutine ttv_arbr_act_c_stv_sgt1(lin)
!t1v(18) ar(13)-br(13)-
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
!-------------------------------------------------------------------
          iwdl=just(lri,lrj)     !
          iwdr=0
          ni=mod(lrj-lri,2)
          w1=w1_t1v
          if(ni.eq.0) w1=-w1
          do mpl=1,mtype
            vplp_w0(mpl)=0.d0
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1
          enddo
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo
          call ar_br_stv_link_ext_brar(lin,lri,lrj)
        enddo
      enddo

      return
      end

      subroutine ss_s_drl_act_c_ext_ab_sgt0(lin)
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
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0_ss(19)
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1_ss(19)
          enddo
          call drl_act_c_link_ext_ab(lin,lri)
        enddo
      enddo

      return
      end

      subroutine ds_arblbr_act_c1_sgt0(lin )
!=======================================================================
!ds(7-2) ar(23)-bl(31)-br(32)-
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
!ds(7-2) ar(23)-bl(31)-br(32)-
          do lrd=norb_frz+1,lri-1
            lmd=lsm_inn(lrd)
            if(lmd.ne.jml) cycle
            ijk=lrd-norb_frz+ngw2(lri-norb_frz)+ngw3(lrj-norb_frz)
            intpos=intind_ijka(ijk)
            w1ds =w1_ds(2)
            ni=mod(norb_dz-lrj+lri-lrd,2)
            if(ni.eq.0) w1ds =-w1ds

            iwdr=just(lrj,lri)
            iwdl=jud(lrd)
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            do mpl=1,mtype
              vplp_w0(mpl)=0.d0
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1ds
            enddo
           call ar_bl_br_ext_al_new(lin,intpos,isma,1)
          enddo
        enddo
      enddo

      return
      end

      subroutine sdd_ar_act_dlr_sgt0(lin,lra)
!sd1(8-1)    ar(01)-
!sd1(8-2)    (11)ar(23)-
!sd1(8-3)    ar(13)-c'(21)-
!sd1(8-4)    ar(23)-c'(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        w0sd1=w0_sd1(1)
        w1sd1=w1_sd1(1)
        w0sd2=w0_sd1(2)
        w1sd2=w1_sd1(2)
        w0sd3=w0_sd1(3)
        w1sd3=w1_sd1(3)
        w0sd4=w0_sd1(4)
        w1sd4=w1_sd1(4)
        if(mod(norb_dz-lri,2).eq.1) then
          w0sd1=-w0sd1
          w1sd1=-w1sd1
          w0sd2=-w0sd2
          w1sd2=-w1sd2
          w0sd3=-w0sd3
          w1sd3=-w1sd3
          w0sd4=-w0sd4
          w1sd4=-w1sd4
        endif
!sd1(8-1)    ar(01)-
        if(jml.eq.1.and.jmr.eq.lmi) then
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
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1sd1
          enddo
          call ar_drl_ext_al_new(lin,lri,lra)
        endif
!sd1(8-2)    (11)ar(23)-
        do mpl=1,mtype
          vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd2
          vplp_w1(mpl)=vplpnew_w1(mpl)*w1sd2
        enddo
        do lrj=norb_frz+1,lri-1
          lmj=lsm_inn(lrj)
          if(jml.ne.mul_tab(lmi,lmj).or.jmr.ne.lmj) cycle
          iwdl=just(lri,lrj)
          iwdr=jud(lrj)
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo
          call ar_drl_ext_al_new(lin,lri,lra)
        enddo
!sd1(8-3)    ar(13)-c'(21)-
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          if(jml.ne.mul_tab(lmi,lmj).or.jmr.ne.lmj) cycle
          iwdl=just(lrj,lri)
          iwdr=jud(lrj)
          do mpl=1,mtype
            vplp_w0(mpl)=-vplpnew_w0(mpl)*w0sd3
            vplp_w1(mpl)=-vplpnew_w1(mpl)*w1sd3
          enddo
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo
          call ar_drl_ext_al_new(lin,lri,lra)
!sd1(8-4)    ar(23)-c'(11)-
          iwdl=just(lri,lrj)
          iwdr=jud(lrj)
          do mpl=1,mtype
            vplp_w0(mpl)=-vplpnew_w0(mpl)*w0sd4
            vplp_w1(mpl)=-vplpnew_w1(mpl)*w1sd4
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

      subroutine tts_arbl_act_c_sgt1(lin,lra)
!=======================================================================
!tts(5) a&r-b^l-  act -b&l ............................................
!t1s(5-1)   ar(13)-bl(10)-
!t1s(5-2)   ar(13)-bl(32)-
!t1s(5-2)   ar(13)-c'(11)-bl(32)-
!t1s(5-3)   ar(13)-bl(31)-c"(12)-
!t1s(5-4)   ar(13)-bl(32)-c"(11)-
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
        w1ts1=w1_t1s(1)
        w1ts2=w1_t1s(2)
        w1ts3=w1_t1s(3)
        w1ts4=w1_t1s(4)
        if(mod(lrj-lri,2).eq.0) then
          w1ts1=-w1ts1
          w1ts2=-w1ts2
          w1ts3=-w1ts3
          w1ts4=-w1ts4
        endif
        ijk=lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz) !???
        intpos=intind_ijka(ijk)                         !???
!t1s(5-1)   ar(13)-bl(10)-
        if(jmr.eq.1.and.jml.eq.lmij) then
          iwdl=just(lri,lrj)
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
          call ar_bl_br_ext_al_new(lin,intpos,isma,1)
        endif
!t1s(5-2)   (11)ar(13)-bl(32)-
        do lrk=norb_frz+1,lri-1
          lmk=lsm_inn(lrk)
         if(jmr.eq.mul_tab(lmk,lmi).and.jml.eq.mul_tab(lmk,lmj)) then
            iwdl=just(lrk,lri)
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
            call ar_bl_br_ext_al_new(lin,intpos,isma,1)
          endif
        enddo
!t1s(5-2)   ar(13)-c'(11)-bl(32)-
        do lrk=lri+1,lrj-1
          lmk=lsm_inn(lrk)
          if(jmr.eq.mul_tab(lmi,lmk).and.jml.eq.mul_tab(lmk,lmj)) then
            iwdl=just(lri,lrk)
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
            call ar_bl_br_ext_al_new(lin,intpos,isma,1)
          endif
        enddo
!t1s(5-3)   ar(13)-bl(31)-c"(12)-
        do lrk=lrj+1,norb_dz
          lmk=lsm_inn(lrk)
          if(jmr.eq.mul_tab(lmi,lmk).and.jml.eq.mul_tab(lmj,lmk)) then
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
            call ar_bl_br_ext_al_new(lin,intpos,isma,1)
!t1s(5-4)   ar(13)-bl(32)-c"(11)-
            iwdl=just(lri,lrk)
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
            call ar_bl_br_ext_al_new(lin,intpos,isma,1)
          endif
        enddo
      enddo
      enddo

      return
      end

      subroutine tts_arbl_act_br_sgt1(lin,lra)
!=======================================================================
!tts(5) a&r-b^l-  act -b&l ............................................
!t1s(5-1)   ar(13)-bl(10)-
!t1s(5-2)   ar(13)-bl(32)-
!t1s(5-2)   ar(13)-c'(11)-bl(32)-
!t1s(5-3)   ar(13)-bl(31)-c"(12)-
!t1s(5-4)   ar(13)-bl(32)-c"(11)-
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
        w1ts1=w1_t1s(1)
        w1ts2=w1_t1s(2)
        w1ts3=w1_t1s(3)
        w1ts4=w1_t1s(4)
        if(mod(lrj-lri,2).eq.0) then
          w1ts1=-w1ts1
          w1ts2=-w1ts2
          w1ts3=-w1ts3
          w1ts4=-w1ts4
        endif
        ijk=lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz) !???
        intpos=intind_ijka(ijk)                         !???
!t1s(5-1)   ar(13)-bl(10)-
        if(jmr.eq.1.and.jml.eq.lmij) then
          iwdl=just(lri,lrj)
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
          call ar_bl_br_ext_al_new(lin,intpos,isma,1)
        endif
!t1s(5-2)   (11)ar(13)-bl(32)-
        do lrk=norb_frz+1,lri-1
          lmk=lsm_inn(lrk)
         if(jml.eq.mul_tab(lmk,lmi).and.jmr.eq.mul_tab(lmk,lmj)) then
            iwdl=just(lrk,lri)
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
            call ar_bl_br_ext_al_new(lin,intpos,isma,1)
          endif
        enddo
!t1s(5-2)   ar(13)-c'(11)-bl(32)-
        do lrk=lri+1,lrj-1
          lmk=lsm_inn(lrk)
          if(jml.eq.mul_tab(lmi,lmk).and.jmr.eq.mul_tab(lmk,lmj)) then
            iwdl=just(lri,lrk)
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
            call ar_bl_br_ext_al_new(lin,intpos,isma,1)
          endif
        enddo
!t1s(5-3)   ar(13)-bl(31)-c"(12)-
        do lrk=lrj+1,norb_dz
          lmk=lsm_inn(lrk)
          if(jml.eq.mul_tab(lmi,lmk).and.jmr.eq.mul_tab(lmj,lmk)) then
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
            call ar_bl_br_ext_al_new(lin,intpos,isma,1)
!t1s(5-4)   ar(13)-bl(32)-c"(11)-
            iwdl=just(lri,lrk)
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
            call ar_bl_br_ext_al_new(lin,intpos,isma,1)
          endif
        enddo
      enddo
      enddo

      return
      end

      subroutine sd_ar_act_dlr_sgt0(lin,lra)
!sd(6-3) a&r(13)c'(22)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      jmlr=mul_tab(jml,jmr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jmlr) cycle
        w0sd3 =w0_sd(3)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.0) then
          w0sd3 =-w0sd3
        endif
!------------------------------------------------------------------
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
           call ar_drl_ext_al_new(lin,lri,lra)
        enddo
      enddo

      return
      end

      subroutine sdd_ar_act_blbr_sgt0(lin,jk)
!sdd(8-1) a&r(01)-
!sdd(8-2) (11)a&(23)-
!sdd(8-3) a&r(13)c'(21)-
!sdd(8-4) a&r(23)c'(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      jmlr=mul_tab(jml,jmr)
      isma=mul_tab(iml,imr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        w0sd1=w0_sd1(1)
        w1sd1=w1_sd1(1)
        w0sd2=w0_sd1(2)
        w1sd2=w1_sd1(2)
        w0sd3=w0_sd1(3)
        w1sd3=w1_sd1(3)
        w0sd4=w0_sd1(4)
        w1sd4=w1_sd1(4)
        if(mod(norb_dz-lri,2).eq.1) then
          w0sd1=-w0sd1
          w1sd1=-w1sd1
          w0sd2=-w0sd2
          w1sd2=-w1sd2
          w0sd3=-w0sd3
          w1sd3=-w1sd3
          w0sd4=-w0sd4
          w1sd4=-w1sd4
        endif
        ijk=lri-norb_frz+jk
        intpos=intind_ijka(ijk)
!sdd(8-1) a&r(01)-
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
!sdd(8-2) c(11)a&(23)-
         do mpl=1,mtype
           vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd2
           vplp_w1(mpl)=vplpnew_w1(mpl)*w0sd2
         enddo
         do lrk=norb_frz+1,lri-1
           lmk=lsm_inn(lrk)
           if(lmk.ne.jmr.or.jml.ne.mul_tab(lmk,lmi)) cycle
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
!-------------------------------------------------------------------
!sdd(8-3) a&r(13)c'(21)-
         do lrk=lri+1,norb_dz
           lmk=lsm_inn(lrk)
           if(lmk.ne.jmr.or.jml.ne.mul_tab(lmi,lmk)) cycle
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
           call ar_bl_br_ext_al_new(lin,intpos,isma,1 )
!sdd(8-4) a&r(23)c'(11)-
           iwdl=just(lri,lrk)
           iwdr=jud(lrk)
           do mpl=1,mhlp
             iwal=lpnew_lwei(mpl)
             iwar=lpnew_rwei(mpl)
             lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
             lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
           enddo
           do mpl=1,mtype
             vplp_w0(mpl)=-vplpnew_w0(mpl)*w0sd4
             vplp_w1(mpl)=-vplpnew_w1(mpl)*w0sd4
           enddo
           call ar_bl_br_ext_al_new(lin,intpos,isma,1 )
         enddo
      enddo

      return
      end

      subroutine sd_ar_act_blbr_sgt0(lin,jk)
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
!-------------------------------------------------------------------
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
          call ar_bl_br_ext_al_new(lin,intpos,isma,1)
        enddo
      enddo

      return
      end

      subroutine sdd_ar_act_c_sd_ext_sgt0(lin)
!sd1(8-1)    ar(01)-
!sd1(8-2)    (11)ar(23)-
!sd1(8-3)    ar(13)-c'(21)-
!sd1(8-4)    ar(23)-c'(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      isma=mul_tab(iml,imr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(jml.ne.1.or.jmr.ne.lmi) cycle
        w0sd1=w0_sd1(1)
        w1sd1=w1_sd1(1)
        if(mod(norb_dz-lri,2).eq.1) then
          w0sd1=-w0sd1
          w1sd1=-w1sd1
        endif
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
           vplp_w1(mpl)=vplpnew_w1(mpl)*w1sd1
        enddo
        if(lin.eq.6) then
          call ar_sd_ext_ar(26,lri,lrj,isma)
          call ar_sd_ext_rest(lri)
        endif
        if(lin.eq.13) then
          call ar_td_ext_ar(26,lri,lrj,isma)
          call ar_td_ext_rest(lri)
        endif
        if(lin.eq.23) then
          call ar_dv_ext_ar(26,isma,lri,lrj)   !ar_dv
        endif
      enddo

       do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jml) cycle

         if(lmi.eq.jmr) then
            iwdl=just(lrj,lri)
            iwdr=jud(lri)
            w0sd2=w0_sd1(2)
!            w1sd2=w0_sd1(2)
            ni=mod(norb_dz-lrj,2)
            if(ni.eq.1) then
!              w1sd2=-w1sd2
              w0sd2=-w0sd2
            endif
!sd1(8-2)    (11)ar(23)-
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd2
              vplp_w1(mpl)=vplpnew_w1(mpl)*w0sd2
            enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            if(lin.eq.6) then
              call ar_sd_ext_ar(25,lrj,lri,isma)
              call ar_sd_ext_rest(lrj)
            endif
            if(lin.eq.13) then
              call ar_td_ext_ar(25,lrj,lri,isma)
              call ar_td_ext_rest(lrj)
            endif
           if(lin.eq.23) then
              call ar_dv_ext_ar(25,isma,lrj,lri)   !ar_dv
            endif
          endif
!sd1(8-3)    ar(13)-c'(21)-
!sd1(8-4)    ar(23)-c'(11)-
          if(lmj.ne.jmr) cycle
         w0sd3=w0_sd1(3)
          w1sd3=w1_sd1(3)
          w0sd4=w0_sd1(4)
          w1sd4=w1_sd1(4)
          ni=mod(norb_dz-lri,2)
          if(ni.eq.0) then
            w0sd3=-w0sd3
            w1sd3=-w1sd3
            w0sd4=-w0sd4
            w1sd4=-w1sd4
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
            if(lin.eq.6) then
              call ar_sd_ext_ar(29,lri,lrj,isma)
              call ar_sd_ext_rest(lri)
            endif
            if(lin.eq.13) then
              call ar_td_ext_ar(29,lri,lrj,isma)
              call ar_td_ext_rest(lri)
            endif
           if(lin.eq.23) then
              call ar_dv_ext_ar(29,isma,lri,lrj)   !ar_dv
            endif
            iwdl=just(lri,lrj)
            iwdr=jud(lrj)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd4
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1sd4
            enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            if(lin.eq.6) then
              call ar_sd_ext_ar(57,lri,lrj,isma)
              call ar_sd_ext_rest(lri)
            endif
            if(lin.eq.13) then
              call ar_td_ext_ar(57,lri,lrj,isma)
              call ar_td_ext_rest(lri)
            endif
           if(lin.eq.23) then
              call ar_dv_ext_ar(57,isma,lri,lrj)   !ar_dv
            endif
          endif
        enddo
      enddo

      return
      end

      subroutine ttdd_ar_act_c_ttdd_ext_sgt1(lin)
!t1d1(15-1)  ar(13)-
!t1d1(15-1)  ar(13)-c'(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      isma=mul_tab(iml,imr)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
      do lrj=norb_frz+1,lri-1
        lmj=lsm_inn(lrj)
        if(jml.ne.mul_tab(lmi,lmj).or.jmr.ne.lmj) cycle
        w0td1=w0_t1d1(1)
        w1td1=w1_t1d1(1)
        if(mod(norb_dz-lri,2).eq.1) then
          w0td1=-w0td1
          w1td1=-w1td1
        endif
        do mpl=1,mtype
           vplp_w0(mpl)=vplpnew_w0(mpl)*w0td1
           vplp_w1(mpl)=vplpnew_w1(mpl)*w1td1
        enddo
!t1d1(15-1)  (11)ar(13)-
        iwdl=just(lrj,lri)
        iwdr=jud(lrj)
        do mpl=1,mhlp
           iwal=lpnew_lwei(mpl)
           iwar=lpnew_rwei(mpl)
           lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
           lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
        enddo
        if(lin.eq.6) then
          call ar_sd_ext_ar(43,lri,lrj,isma)
          call ar_sd_ext_rest(lri)
        endif
        if(lin.eq.13) then
          call ar_td_ext_ar(43,lri,lrj,isma)
          call ar_td_ext_rest(lri)
        endif
        if(lin.eq.23) then
          call ar_dv_ext_ar(43,isma,lri,lrj)
        endif
      enddo
!t1d1(15-1)  ar(13)-c'(11)-
      do lrj=lri+1,norb_dz
        lmj=lsm_inn(lrj)
        if(jml.ne.mul_tab(lmi,lmj).or.jmr.ne.lmj) cycle
        w0td1=w0_t1d1(1)
        w1td1=w1_t1d1(1)
        if(mod(norb_dz-lri,2).eq.1) then
          w0td1=-w0td1
          w1td1=-w1td1
        endif
        do mpl=1,mtype
           vplp_w0(mpl)=-vplpnew_w0(mpl)*w0td1
           vplp_w1(mpl)=-vplpnew_w1(mpl)*w1td1
        enddo
        iwdl=just(lri,lrj)
        iwdr=jud(lrj)
        do mpl=1,mhlp
           iwal=lpnew_lwei(mpl)
           iwar=lpnew_rwei(mpl)
           lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
           lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
        enddo
        if(lin.eq.6) then
          call ar_sd_ext_ar(46,lri,lrj,isma)
          call ar_sd_ext_rest(lri)
        endif
        if(lin.eq.13) then
          call ar_td_ext_ar(46,lri,lrj,isma)
          call ar_td_ext_rest(lri)
        endif
        if(lin.eq.23) then
          call ar_dv_ext_ar(46,isma,lri,lrj)
        endif
      enddo
      enddo

      return
      end

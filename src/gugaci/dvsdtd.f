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
      subroutine dv_drt_ci_new()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      common/lpdisk/idisk_lp,idisk_array(13)

      call external_space_plpmode_value_dv()

      idisk_lp=idisk_array(2)
      do lpb=1,lpblock_dv
        call read_lp()
        ipael=iml+1
        ipae=1
        call get_jpty(jpadlr,jptyl,jptyr)
        call get_jp(jptyl,jml,jpadl,1)
        call get_jp(jptyr,jmr,jpad,1)
c        jmlr=mul_tab(jml,jmr)
        call gsd_determine_extarmode_paras(iml,imr,.false.)
        if(linelp.le.12)   then
          call dv_ext_head_in_act()
        else
          call dv_ext_head_in_dbl()
        endif
      enddo

      return
      end

      subroutine dv_ext_head_in_dbl()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      logic_dh=.true.
      isma=iml
      lpok=jpadlr
      jmlr=mul_tab(jml,jmr)
      goto(101,102,103,104,105,106,10,108,109,110,111,112,    ! 107???
     :     113,10,115,10,117,118,119,120,121,122,123,124,125,10),lpok
!====================================================================
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
106   if(linelp.ne.13) goto 107
      if(jb_sys.gt.0) call sd_adb_act_c_ext_ar(23)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jmlr) cycle
        w0sd1=w0_sd(1)
        w0sd9=w0_sd(9)
        w1sd9=w1_sd(9)
        w0sd12=w0_sd(12)
        w1sd12=w1_sd(12)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1) w0sd1=-w0sd1
        if(ni.eq.1) w0sd9=-w0sd9
        if(ni.eq.1) w1sd9=-w1sd9
        if(ni.eq.1) w0sd12=-w0sd12
        if(ni.eq.1) w1sd12=-w1sd12
        if(jml.eq.1.and.lmi.eq.jmr) then
          iwdl=just(lri,lri)
          iwdr=jud(lri)
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo

!sd(6-1) a&r(02)-
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd1
            vplp_w1(mpl)=vplpnew_w1(mpl)*w0sd1
          enddo
          call ar_dv_ext_ar(26,isma,lri,lrj)   !ar_dv

!sd(6-12) d&r&l(33)b^l(02)
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd12
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1sd12
          enddo
          do lrk=1,lri-1
            call drl_bl_ext_ar_new(23,lrk,lri)
          enddo

!sd(6-9) d&r&r(03)b^r(32)
          do lrk=norb_frz+1,lri-1
            iwdl=just(lrk,lrk)
            iwdr=jud(lri)
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd9
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1sd9
            enddo
            call drr_br_ext_ar(23,lrk,lri)
          enddo
        endif
      enddo

      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jml) cycle

         if(lmi.eq.jmr) then
            iwdl=just(lri,lrj)
            iwdr=jud(lri)
            w0sd2=w0_sd(2)
            w0sd11=w0_sd(11)
            w1sd11=w1_sd(11)
            w0sd14=w0_sd(14)
            ni=mod(norb_dz-lrj,2)
            if(ni.eq.1) then
              w0sd2=-w0sd2
              w0sd11=-w0sd11
              w1sd11=-w1sd11
              w0sd14=-w0sd14
            endif
!sd(6-2) c(22)-a&r(13)-
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd2
            enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            call ar_dv_ext_ar(25,isma,lrj,lri) !ar_dv
!sd(6-11) d&r&l(22)b^l(13)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd11
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1sd11
            enddo
            call drl_bl_ext_ar_new(23,lri,lrj)
!sd(6-14) (22)d&r&l(33)b^l(13)
!sd(6-14) d&r&l(33)c"(22)b^l(13)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd14
              vplp_w1(mpl)=0.d0
            enddo
            do lrk=1,lrj-1
              if(lrk.eq.lri) cycle
              call drl_bl_ext_ar_new(23,lrk,lrj)
            enddo
          endif
!sd(6-4) a&r(23)-c'(12)-
          w0sd4=w0_sd(4)
          w0sd16=w0_sd(16)
          ni=mod(norb_dz-lri,2)
          if(ni.eq.0) then
            w0sd4=-w0sd4
            w0sd16=-w0sd16
          endif
          if(lmj.eq.jmr) then
            iwdl=just(lri,lrj)
            iwdr=jud(lrj)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd4
            enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            call ar_dv_ext_ar(28,isma,lri,lrj)   !ar_dv
!sd(6-16) d&r&l(33)b^l(23)c'(12)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd16
              vplp_w1(mpl)=0.d0
            enddo
            do lrk=1,lri-1
              call drl_bl_ext_ar_new(23,lrk,lri)
            enddo
          endif
!sd(6-5) a&r(23)b&r(13)b^r(32)
          do lrd=lrj+1,norb_dz
            lmd=lsm_inn(lrd)
           if(lmd.ne.jmr) cycle

            w0sd5=w0_sd(5)
            w1sd5=w1_sd(5)
            ni=mod(lrj-lri+norb_dz-lrd,2)
            if(ni.eq.0) then
              w0sd5=-w0sd5
              w1sd5=-w1sd5
            endif
            iwdl=just(lri,lrj)
            iwdr=jud(lrd)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd5
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1sd5
            enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            ijk=lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lrd-norb_frz)
            intpos=intind_ijka(ijk)
            call ar_br_br_ext_ar_new(23,intpos,isma)
          enddo
!sd(6-8) a&r(23)b&l(32)b^l(13)
          do lrd=lri+1,lrj-1
            lmd=lsm_inn(lrd)
            if(lmd.ne.jmr) cycle
            iwdl=just(lri,lrj)
            iwdr=jud(lrd)
            w0sd8=w0_sd(8)
            w1sd8=w1_sd(8)
            ni=mod(lrj-lri+norb_dz-lrd,2)
            if(ni.eq.0) then
              w0sd8=-w0sd8
              w1sd8=-w1sd8
            endif
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd8
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1sd8
            enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            ijk=lri-norb_frz+ngw2(lrd-norb_frz)+ngw3(lrj-norb_frz)
            intpos=intind_ijka(ijk)
            call ar_bl_bl_ext_ar_new(23,intpos,isma,1)
          enddo
        enddo
      enddo
!====================================================================
!sd(6)    act -br-br-
107   if(linelp.eq.20) then
        jk=nlg1
        call sd_ar_act_brbr(23,jk)
        if(jb_sys.gt.0) call sd_ar_act_brbr_sgt0(23,jk)
      endif
!sd(6)  act -bl-bl-
      if(linelp.eq.22) then
        jk=nlg1
        call sd_ar_act_blbl(23,jk)
        if(jb_sys.gt.0) call sd_ar_act_blbl_sgt0(23,jk)
      endif
      return
!=======================================================================
!td(13)                                  act -c'-
!td(13-1) (22)a&(23)
!td(13-1) a&(23)c'(22)
!td(13-2) a&(23)b&r(23)b^r(32)
!td(13-3) a&(23)b&l(32)b^l(23)
!td(13-4) d&r&l(22)b^l(23)
!td(13-5) (22)d&&l(33)b^l(23)
!td(13-5) d&rl(33)c"(22)b^l(23)
!td(13-5) d&rl(33)b^l(23)c'(22)
113   if(linelp.ne.13) goto 132
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jml) cycle
          w0td1=w0_td(1)
          w0td4=w0_td(4)
          w1td4=w1_td(4)
          w0td5=w0_td(5)
          ni=mod(norb_dz-lrj,2)
          if(ni.eq.1) then
            w0td1=-w0td1
            w0td4=-w0td4
            w1td4=-w1td4
            w0td5=-w0td5
          endif
          iwdl=just(lri,lrj)     !
!-------------------------------------------------------------------
!td(13-1) (22)a&(23)
          if(lmi.eq.jmr) then
            iwdr=jud(lri)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0td1
            enddo
              do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            call ar_dv_ext_ar(43,isma,lrj,lri)
!td(13-4) d&r&l(22)b^l(23)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0td4
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1td4
            enddo
            call drl_bl_ext_ar_new(23,lri,lrj)
!td(13-5) (22)d&rl(33)b^l(23)
!td(13-5) d&rl(33)c"(22)b^l(23)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0td5
              vplp_w1(mpl)=0.d0
            enddo
            do lrk=1,lrj-1
              if(lrk.eq.lri) cycle
              call drl_bl_ext_ar_new(23,lrk,lrj)
            enddo
          endif
!-------------------------------------------------------------------
!td(13-1) a&(23)c'(22)
          if(lmj.eq.jmr) then
            iwdr=jud(lrj)
            w0td1=w0_td(1)
            w0td5=w0_td(5)
            ni=mod(norb_dz-lri,2)
            if(ni.eq.0) then
              w0td1=-w0td1
              w0td5=-w0td5
            endif
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0td1
            enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            call ar_dv_ext_ar(46,isma,lri,lrj) !ar_dv
!td(13-5) d&rl(33)b^l(23)c'(22)
              do mpl=1,mtype
                vplp_w0(mpl)=vplpnew_w0(mpl)*w0td5
                vplp_w1(mpl)=0.d0
              enddo
              do lrk=1,lri-1
                call drl_bl_ext_ar_new(23,lrk,lri)
              enddo
            endif
!-------------------------------------------------------------------
!td(13-2) a&(23)b&r(23)b^r(32)
         do lrd=lrj+1,norb_dz
            lmd=lsm_inn(lrd)
            if(lmd.ne.jmr) cycle
            iwdr=jud(lrd)
            ijk=lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lrd-norb_frz)
            intpos=intind_ijka(ijk)
            w0td2=w0_td(2)
            w1td2=w1_td(2)
            ni=mod(lrj-lri+norb_dz-lrd,2)
            if(ni.eq.0) then
              w0td2=-w0td2
              w1td2=-w1td2
            endif
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0td2
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1td2
            enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            call ar_br_br_ext_ar_new(23,intpos,isma)
          enddo
!-------------------------------------------------------------------
!td(13-3) a&(23)b&l(32)b^l(23)
          do lrd=lri+1,lrj-1
            lmd=lsm_inn(lrd)
            if(lmd.ne.jmr) cycle
            iwdr=jud(lrd)
            w0td3=w0_td(3)
            w1td3=w1_td(3)
            ni=mod(lrj-lri+norb_dz-lrd,2)
            if(ni.eq.0) then
              w0td3=-w0td3
              w1td3=-w1td3
            endif
            ijk=lri-norb_frz+ngw2(lrd-norb_frz)+ngw3(lrj-norb_frz)
            intpos=intind_ijka(ijk)
              do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0td3
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1td3
            enddo
              do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            call ar_bl_bl_ext_ar_new(23,intpos,isma,1)
          enddo
!-------------------------------------------------------------------
        enddo
      enddo
!=======================================================================
!td(13) act -br-br-
132   if(linelp.eq.20) then
        jk=nlg1
        call td_ar_act_brbr(23,jk)
      endif
!td(13)  act -bl-bl-
      if(linelp.eq.22) then
        jk=nlg1
        call td_ar_act_blbl(23,jk)
      endif
      return
!=======================================================================
!dv(23) act -c'-..................................................
!dv(23-1) ar(23)-
!dv(23-2) drl(33)-bl(23)-
123   if(linelp.ne.13) goto 232
      imap_1=0
      do lrd=norb_frz+1,norb_dz
        lmd=lsm_inn(lrd)
        if(lmd.ne.jml) cycle
        ni=mod(norb_dz-lrd,2)

!dv(23-1) ar(23)-
        w0=w0_dv(1)
        if(ni.eq.1)w0=-w0_dv(1)
          iwdl=jud(lrd)
          iwdr=0
          imap_1=imap_1+1
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0
          enddo
        do mpl=1,mhlp
          iwal=lpnew_lwei(mpl)
          iwar=lpnew_rwei(mpl)
          lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
        enddo
        call ar_dv_ext_ar(51,isma,lrd,0)   !ar_dv
!dv(23-2) drl(33)-bl(23)-
        w0=w0_dv(2)
        if(ni.eq.1)w0=-w0_dv(2)
        do mpl=1,mtype
          vplp_w0(mpl)=vplpnew_w0(mpl)*w0
          vplp_w1(mpl)=0.d0
        enddo
        do lrk=1,lrd-1
          call drl_bl_ext_ar_new(23,lrk,lrd)
        enddo
      enddo
!=======================================================================
!dv(23) act -br-br-..................................................
232   if(linelp.ne.20) goto 233
      jk=nlg1
      do lrd=norb_frz+1,norb_dz
        lmd=lsm_inn(lrd)
        if(lmd.ne.jml) cycle
        ni=mod(norb_dz-lrd,2)
        w0=w0_dv(1)
        if(ni.eq.1)w0=-w0_dv(1)
!..............................................................
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
        call ar_br_br_ext_ar_new(23,intpos,isma)
      enddo
!dv(23) act -bl-bl-.....................................................
233   if(linelp.ne.22) return
      jk=nlg1
      do lrd=norb_frz+1,norb_dz
        lmd=lsm_inn(lrd)
        if(lmd.ne.jml) cycle
        ni=mod(norb_dz-lrd,2)
        w0=w0_dv(1)
        if(ni.eq.1)w0=-w0_dv(1)
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
        call ar_bl_bl_ext_ar_new(23,intpos,isma,1)
      enddo
      return
!=======================================================================
!ss(1)     act -b^l-
101   if(linelp.ne.15) return
      lra=nlg1
      if(nlg2.eq.1) then
        call ss_drl_act_bl(23,lra)
        if(jb_sys.gt.0) then
           call ss_drl_act_bl_sgt0(23,lra)
        endif
      else
        call ss_arbl_act_bl(23,lra)
        if(jb_sys.gt.0) then
          call ss_arbl_act_bl_sgt0(23,lra)
          call ss_s_drl_act_bl_sgt0(23,lra)
        endif
      endif
      return
!=======================================================================
!st(2)      act -b^l
102   if(linelp.ne.15.or.nlg2.ne.2) return
      lra=nlg1
      call st_arbl_act_bl(23,lra)
      if(jb_sys.gt.0) call st_arbl_act_bl_sgt0(23,lra)
      if(jml.ne.jmr) return
      call st_drl_act_bl(23,lra)
      if(jb_sys.gt.0) call st_drl_act_bl_sgt0(23,lra)
      return
!=======================================================================
!ts(3)   act -b^l ............................................
103   if(linelp.ne.15) return
      lra=nlg1
      if(nlg2.ne.2) return
      call ts_arbl_act_bl(23,lra)
      if(jb_sys.gt.0) call ts_arbl_act_bl_sgt0(23,lra)
      return
!=======================================================================
!=======================================================================
!stt(4) act-b^l ..............................................
104   if(linelp.ne.15) return
      lra=nlg1
      if(nlg2.eq.1) return
      call stt_arbl_act_bl_sgt1(23,lra)
      return
!=======================================================================
!tts(5) act-b^l ..............................................
105   if(linelp.ne.15.or.nlg2.ne.2) return
      lra=nlg1
      call tts_arbl_act_bl_sgt1(23,lra)
      call tts_drl_act_bl_sgt1(23,lra)
      return
!=======================================================================
!sdd(8) act-b^l -c'- -c"-..............................................

108   if(linelp.eq.13) then
        lra=nlg1
         call sdd_drlbl_act_c_sgt0(23)
         call sdd_drrbr_act_c_sgt0(23)
         call sdd_abb_act_c_sgt0(23)
         call sdd_ar_act_c_sd_ext_sgt0(23)
      else
        lra=nlg1
        if(linelp.eq.20) call sdd_ar_act_brbr_sgt0(23,lra)
        if(linelp.eq.22) call sdd_ar_act_blbl_sgt0(23,lra)
      endif
      return
!=======================================================================
!dds(9) act -c'- ..............................................
109   if(linelp.eq.13) then
        lra=nlg1
        call dds_abb_act_c_sgt0(23)
      endif
      return
!=======================================================================
!tt(11)   act -b^l
111   if(linelp.ne.15) return
      lra=nlg1
      if(nlg2.eq.1) then
        call tt_drl_act_bl(23,lra)
      else
        call tt_arbl_act_bl(23,lra)
      endif
      return
!=======================================================================
!tttt(12)   act -b^l
112   if(linelp.ne.15) return
      lra=nlg1
      if(nlg2.eq.1) then
        call tttt_drl_act_bl_sgt1(23,lra)
      else
        call tttt_arbl_act_bl_sgt1(23,lra)
      endif
      return
!=======================================================================
!ttv(19) act -b^r- ....................................................
118   if(linelp.ne.16.or.nlg2.ne.2) return
      lra=nlg1
      call ttv_arbr_act_c_sgt1(23,lra)
      return
!=======================================================================
!dd(19) act -b^l- ....................................................
119   if(linelp.ne.15) return
      lra=nlg1
      if(nlg2.eq.1) then
        call dd_drl_act_bl(23,lra)
      else
        call dd_arbl_act_bl(23,lra)
      endif
      return
!=======================================================================
!dddd(19) act -b^l- -br- ...............................................
120   lra=nlg1
      if(linelp.eq.15) then
        if(nlg2.eq.2) then
          call  dddd_arbl_act_bl_sgt0(23,lra)
        else
          call dddd_drl_act_bl_sgt0(23,lra)
        endif
      endif
      return
!=======================================================================
!dd1(21) act -b^l-  ....................................................
121   if(linelp.ne.15.or.nlg2.ne.2) return
      lra=nlg1
      call dd1_arbl_act_bl_sgt0(23,lra)
      return
!=======================================================================
!d1d(21) act -b^l-  ....................................................
122   lra=nlg1
      if(linelp.eq.15.and.nlg2.eq.2) then
        call  d1d_arbl_act_bl_sgt0(23,lra)
        call d1d_drl_act_bl_sgt0(23,lra)
      endif
      return
!=======================================================================
!d1v(24) act -bl-bl -br-br- -c'- .......................................
124   lra=nlg1
      if(linelp.eq.13) then
        call d1v_drl_bl_act_c_sgt0(23)    ! drl-bl- ar-
      endif
      if(linelp.eq.22) then
        call d1v_ar_act_blbl_sgt0(23,lra)
      endif
      if(linelp.eq.20) then
        call d1v_ar_act_brbr_sgt0(23,lra)
      endif
      return
!=======================================================================
!vv(25) act -b^l- ....................................................
125   if(linelp.ne.15) return
      lra=nlg1
      if(nlg2.eq.1) then
        iwdl=0
        iwdr=0
        do mpl=1,mtype
          vplp_w0(mpl)=vplpnew_w0(mpl)*w0_vv
          vplp_w1(mpl)=0.d0
        enddo
        do mpl=1,mhlp
          iwal=lpnew_lwei(mpl)
          iwar=lpnew_rwei(mpl)
          lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
        enddo
        call drl_bl_sum_ar_new(23,0,0,lra)
c        do lrk=1,norb_dz
c          call drl_bl_ext_ar_new(23,lrk,lra)
c        enddo
       endif
      return
!=======================================================================
!sv(10)     act -b^r-
!sv(10-1) ar(13)-br(23)-
!sv(10-2) ar(23)-br(13)-
!sv(10-3) drr(03)-
110   if(linelp.ne.16.or.nlg2.ne.2) return
      lra=nlg1
      if(jb_sys.gt.0) then
        call sv_arbr_act_br_sgt0(23,lra)
      endif
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        do lrj=lri,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jml) cycle
          w0sv2=w0_sv(2)
          w1sv2=w1_sv(2)
          ni=mod(lrj-lri,2)
          if(ni.eq.0) then
            w0sv2=-w0sv2
            w1sv2=-w1sv2
          endif
!-------------------------------------------------------------------
          iwdl=just(lri,lrj)
          iwdr=0
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo
         if(lri.ne.lrj) then
!sv(10-2) ar(23)-br(13)-
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sv2
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1sv2
            enddo
            ijk=lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz)
           intpos=intind_ijka(ijk)
            call ar_br_br_ext_ar_new(23,intpos,isma)
          else
!sv(10-3) drr(03)-
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0_sv(3)
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1_sv(3)
            enddo
            call drr_br_ext_ar(23,lri,lra)
          endif
        enddo
      enddo
      return
!=======================================================================
!ttdd(15)    act -c'-   -br-br- -bl-bl
115   lra=nlg1
      if(linelp.eq.13) then
          call  ttdd_drlbl_act_c_sgt1(23)
          call ttdd_abb_act_c_sgt1(23)             ! ????
          call ttdd_ar_act_c_ttdd_ext_sgt1(23)
      endif
      if(linelp.eq.20) then
        call ttdd_ar_act_brbr_sgt0(23,lra)
      endif
      if(linelp.eq.22) then
        call ttdd_ar_act_blbl_sgt0(23,lra)
      endif
      return
!=======================================================================
!tv(17)    act -b^r-
117   if(linelp.ne.16.or.nlg2.ne.2) return
      lra=nlg1
      do lri=norb_frz+1,norb_dz-1
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jml) cycle

          w1tv=w1_tv
          ni=mod(lrj-lri,2)
          if(ni.eq.0)  w1tv=-w1tv

!-------------------------------------------------------------------
!tv(17) ar(23)-br(23)-
          iwdl=just(lri,lrj)     !
          iwdr=0
          do mpl=1,mtype
            vplp_w0(mpl)=0.d0
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1tv
          enddo
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo
          ijk=lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz)
          intpos=intind_ijka(ijk)
          call ar_br_br_ext_ar_new(23,intpos,isma)
        enddo
      enddo
      return
!=======================================================================
10    return
      end

      subroutine dv_ext_head_in_act()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      logic_dh=.false.
      lri=nlg1
      lrj=nlg2
      intpos=nlg1
      isma=iml

      goto (1,2,3,4,5,6,7,8,9,10,11,12),linelp
!line=1 a&r<-->a^r
1     call ar_dv_ext_ar(51,isma,lri,lrj)      ! ar_dv_ext_ar(100_???
      goto 100
!line=4 a&r--b&r--b^r<-->a^r
4     call ar_br_br_ext_ar_new(23,intpos,isma)
      goto 100
!line=7 a&r--b&l--b^l<-->a^r
7     call ar_bl_bl_ext_ar_new(23,intpos,isma,1)
      goto 100
!line=10 d&rr--b^r<-->a^r
10    call drr_br_ext_ar(23,lri,lrj)
      goto 100
!line=12 d&rl--b^l<-->a^r
12    call drl_bl_ext_ar_new(23,lri,lrj)
      goto 100
2     goto 100
3     goto 100
5     goto 100
6     goto 100
8     goto 100
9     goto 100
11    goto 100
100   return
      end

      subroutine sd_drt_ci_new_den()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      common/lpdisk/idisk_lp,idisk_array(13)
      data dsq2/1.414213562373d0/vsq2/0.7071067811865d0/
c      data dsq3/ 1.732050807569d0/
c      write(6,*)'  sd_wyb'

!      logic_sd=.true.
      call external_space_plpmode_value_sd()

      idisk_lp=idisk_array(11)

      do lpb=1,lpblock_sd
        call read_lp()
        ipael=iml+17
        ipae =imr+1
        call get_jpty(jpadlr,jptyl,jptyr)
        call get_jp(jptyl,jml,jpadl,1)
        call get_jp(jptyr,jmr,jpad,1)
c        jmlr=mul_tab(jml,jmr)
        call gsd_determine_extarmode_paras(iml,imr,.true.)
        if(linelp.le.12)   then
          call sd_ext_head_in_act()
        else
          call sd_ext_head_in_dbl()
        endif
      enddo
      return
      end


      subroutine sd_drt_ci_new()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      common/lpdisk/idisk_lp,idisk_array(13)
      data dsq2/1.414213562373d0/vsq2/0.7071067811865d0/
c      data dsq3/ 1.732050807569d0/
c      write(6,*)'  sd_wyb'

!      logic_sd=.true.
      call external_space_plpmode_value_sd()

      idisk_lp=idisk_array(11)

      do lpb=1,lpblock_sd
        call read_lp()
        ipael=iml+17
        ipae =imr+1
        call get_jpty(jpadlr,jptyl,jptyr)
        call get_jp(jptyl,jml,jpadl,1)
        call get_jp(jptyr,jmr,jpad,1)
c        jmlr=mul_tab(jml,jmr)
        call gsd_determine_extarmode_paras(iml,imr,.true.)
        if(linelp.le.12)   then
          call sd_ext_head_in_act()
        else
          call sd_ext_head_in_dbl()
        endif
      enddo
      return
      end

      subroutine sd_ext_head_in_dbl()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      logic_dh=.true.
      isma=mul_tab(iml,imr)
      jmlr=mul_tab(jml,jmr)
      lpok=jpadlr
      goto(101,102,103,104,105,106,10,108,10,110,111,112, ! 109???
     :     113,10,115,10,117,118,119,120,121,122,123,124,125,10),lpok
!====================================================================
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
106   if(linelp.ne.13) goto 107
      if(jb_sys.gt.0) then
        call sd_adb_act_c_ext_ar(6)
      endif
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jmlr) cycle
        w0sd1=w0_sd(1)
        w0sd9=w0_sd(9)
        w1sd9=w1_sd(9)
        w0sd12=w0_sd(12)
        w1sd12=w1_sd(12)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1) w0sd1=-w0sd1
        if(ni.eq.1) w0sd9=-w0sd9
        if(ni.eq.1) w1sd9=-w1sd9
        if(ni.eq.1) w0sd12=-w0sd12
        if(ni.eq.1) w1sd12=-w1sd12
        if(jml.eq.1.and.lmi.eq.jmr) then
          iwdl=just(lri,lri)
          iwdr=jud(lri)
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo

!sd(6-1) a&r(02)-
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd1
            vplp_w1(mpl)=vplpnew_w1(mpl)*w0sd1
          enddo
          call ar_sd_ext_ar(26,lri,lrj,isma)
          call ar_sd_ext_rest(lri)

!sd(6-12) d&r&l(33)b^l(02)
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd12
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1sd12
          enddo
          do lrk=1,lri-1
            call drl_bl_ext_ar_new(6,lrk,lri)
          enddo
!sd(6-9) d&r&r(03)b^r(32)
          do lrk=norb_frz+1,lri-1
            iwdl=just(lrk,lrk)
            iwdr=jud(lri)
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd9
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1sd9
            enddo
            call drr_br_ext_ar(6,lrk,lri)
          enddo
        endif
      enddo

      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jml) cycle

         if(lmi.eq.jmr) then
            iwdl=just(lri,lrj)
            iwdr=jud(lri)
            w0sd2=w0_sd(2)
            w0sd11=w0_sd(11)
            w1sd11=w1_sd(11)
            w0sd14=w0_sd(14)
            ni=mod(norb_dz-lrj,2)
            if(ni.eq.1) then
              w0sd2=-w0sd2
              w0sd11=-w0sd11
              w1sd11=-w1sd11
              w0sd14=-w0sd14
            endif
!sd(6-2) c(22)-a&r(13)-
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd2
            enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            call ar_sd_ext_ar(25,lrj,lri,isma)
            call ar_sd_ext_rest(lrj)
!sd(6-11) d&r&l(22)b^l(13)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd11
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1sd11
            enddo
            call drl_bl_ext_ar_new(6,lri,lrj)
!sd(6-14) (22)d&r&l(33)b^l(13)
!sd(6-14) d&r&l(33)c"(22)b^l(13)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd14
              vplp_w1(mpl)=0.d0
            enddo
            do lrk=1,lrj-1
              if(lrk.eq.lri) cycle
              call drl_bl_ext_ar_new(6,lrk,lrj)
            enddo
          endif
!sd(6-4) a&r(23)-c'(12)-
          w0sd4=w0_sd(4)
          w0sd16=w0_sd(16)
          ni=mod(norb_dz-lri,2)
          if(ni.eq.0) then
            w0sd4=-w0sd4
            w0sd16=-w0sd16
          endif
          if(lmj.eq.jmr) then
            iwdl=just(lri,lrj)
            iwdr=jud(lrj)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd4
            enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            call ar_sd_ext_ar(28,lri,lrj,isma)
            call ar_sd_ext_rest(lri)
!sd(6-16) d&r&l(33)b^l(23)c'(12)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd16
              vplp_w1(mpl)=0.d0
            enddo
            do lrk=1,lri-1
              call drl_bl_ext_ar_new(6,lrk,lri)
            enddo
          endif
!sd(6-5) a&r(23)b&r(13)b^r(32)
          do lrd=lrj+1,norb_dz
            lmd=lsm_inn(lrd)
           if(lmd.ne.jmr) cycle

            w0sd5=w0_sd(5)
            w1sd5=w1_sd(5)
            ni=mod(lrj-lri+norb_dz-lrd,2)
            if(ni.eq.0) then
              w0sd5=-w0sd5
              w1sd5=-w1sd5
            endif
            iwdl=just(lri,lrj)
            iwdr=jud(lrd)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd5
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1sd5
            enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            ijk=lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lrd-norb_frz)
            intpos=intind_ijka(ijk)
            call ar_br_br_ext_ar_new(6,intpos,isma)
          enddo
!sd(6-8) a&r(23)b&l(32)b^l(13)
          do lrd=lri+1,lrj-1
            lmd=lsm_inn(lrd)
            if(lmd.ne.jmr) cycle
            iwdl=just(lri,lrj)
            iwdr=jud(lrd)
            w0sd8=w0_sd(8)
            w1sd8=w1_sd(8)
            ni=mod(lrj-lri+norb_dz-lrd,2)
            if(ni.eq.0) then
              w0sd8=-w0sd8
              w1sd8=-w1sd8
            endif
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd8
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1sd8
            enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            ijk=lri-norb_frz+ngw2(lrd-norb_frz)+ngw3(lrj-norb_frz)
            intpos=intind_ijka(ijk)
            call ar_bl_bl_ext_ar_new(6,intpos,isma,1)
          enddo
        enddo
      enddo
!====================================================================
!sd(6)    act -br-br-
107   if(linelp.eq.20) then
        jk=nlg1
        call sd_ar_act_brbr(6,jk)
        if(jb_sys.gt.0) call sd_ar_act_brbr_sgt0(6,jk)
      endif
!sd(6)  act -bl-bl-
      if(linelp.eq.22) then
        jk=nlg1
        call sd_ar_act_blbl(6,jk)
        if(jb_sys.gt.0) call sd_ar_act_blbl_sgt0(6,jk)
      endif
      return
!=======================================================================
!td(13)                                  act -c'-
!td(13-1) (22)a&(23)
!td(13-1) a&(23)c'(22)
!td(13-2) a&(23)b&r(23)b^r(32)
!td(13-3) a&(23)b&l(32)b^l(23)
!td(13-4) d&r&l(22)b^l(23)
!td(13-5) (22)d&&l(33)b^l(23)
!td(13-5) d&rl(33)c"(22)b^l(23)
!td(13-5) d&rl(33)b^l(23)c'(22)
113   if(linelp.ne.13) goto 132
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jml) cycle
          w0td1=w0_td(1)
          w0td4=w0_td(4)
          w1td4=w1_td(4)
          w0td5=w0_td(5)
          ni=mod(norb_dz-lrj,2)
          if(ni.eq.1) then
            w0td1=-w0td1
            w0td4=-w0td4
            w1td4=-w1td4
            w0td5=-w0td5
          endif
          iwdl=just(lri,lrj)       !
!-------------------------------------------------------------------
!td(13-1) (22)a&(23)
          if(lmi.eq.jmr) then
            iwdr=jud(lri)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0td1
            enddo
              do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            call ar_sd_ext_ar(43,lrj,lri,isma)
            call ar_sd_ext_rest(lrj)
!td(13-4) d&r&l(22)b^l(23)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0td4
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1td4
            enddo
            call drl_bl_ext_ar_new(6,lri,lrj)
!td(13-5) (22)d&rl(33)b^l(23)
!td(13-5) d&rl(33)c"(22)b^l(23)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0td5
              vplp_w1(mpl)=0.d0
            enddo
            do lrk=1,lrj-1
              if(lrk.eq.lri) cycle
              call drl_bl_ext_ar_new(6,lrk,lrj)
            enddo
          endif
!-------------------------------------------------------------------
!td(13-1) a&(23)c'(22)
          if(lmj.eq.jmr) then
            iwdr=jud(lrj)
            w0td1=w0_td(1)
            w0td5=w0_td(5)
            ni=mod(norb_dz-lri,2)
            if(ni.eq.0) then
              w0td1=-w0td1
              w0td5=-w0td5
            endif
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0td1
            enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            call ar_sd_ext_ar(46,lri,lrj,isma)
            call ar_sd_ext_rest(lri)
!td(13-5) d&rl(33)b^l(23)c'(22)
              do mpl=1,mtype
                vplp_w0(mpl)=vplpnew_w0(mpl)*w0td5
                vplp_w1(mpl)=0.d0
              enddo
              do lrk=1,lri-1
                call drl_bl_ext_ar_new(6,lrk,lri)
              enddo
            endif
!-------------------------------------------------------------------
!td(13-2) a&(23)b&r(23)b^r(32)
         do lrd=lrj+1,norb_dz
            lmd=lsm_inn(lrd)
            if(lmd.ne.jmr) cycle
            iwdr=jud(lrd)
            ijk=lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lrd-norb_frz)
            intpos=intind_ijka(ijk)
            w0td2=w0_td(2)
            w1td2=w1_td(2)
            ni=mod(lrj-lri+norb_dz-lrd,2)
            if(ni.eq.0) then
              w0td2=-w0td2
              w1td2=-w1td2
            endif
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0td2
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1td2
            enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            call ar_br_br_ext_ar_new(6,intpos,isma)
          enddo
!-------------------------------------------------------------------
!td(13-3) a&(23)b&l(32)b^l(23)
          do lrd=lri+1,lrj-1
            lmd=lsm_inn(lrd)
            if(lmd.ne.jmr) cycle
            iwdr=jud(lrd)
            w0td3=w0_td(3)
            w1td3=w1_td(3)
            ni=mod(lrj-lri+norb_dz-lrd,2)
            if(ni.eq.0) then
              w0td3=-w0td3
              w1td3=-w1td3
            endif
            ijk=lri-norb_frz+ngw2(lrd-norb_frz)+ngw3(lrj-norb_frz)
            intpos=intind_ijka(ijk)
              do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0td3
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1td3
            enddo
              do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            call ar_bl_bl_ext_ar_new(6,intpos,isma,1)
          enddo
!-------------------------------------------------------------------
        enddo
      enddo
!=======================================================================
!td(13) act -br-br-
132   if(linelp.eq.20) then
        jk=nlg1
        call td_ar_act_brbr(6,jk)
      endif
!td(13)  act -bl-bl-
      if(linelp.eq.22) then
        jk=nlg1
        call td_ar_act_blbl(6,jk)
      endif
      return
!=======================================================================
!dv(23) act -c'-..................................................
!dv(23-1) ar(23)-
!dv(23-2) drl(33)-bl(23)-
123   if(linelp.ne.13) goto 232
      imap_1=0
      do lrd=norb_frz+1,norb_dz
        lmd=lsm_inn(lrd)
        if(lmd.ne.jml) cycle
        ni=mod(norb_dz-lrd,2)

!dv(23-1) ar(23)-
        w0=w0_dv(1)
        if(ni.eq.1)w0=-w0_dv(1)
          iwdl=jud(lrd)
          iwdr=0
          imap_1=imap_1+1
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0
          enddo
        do mpl=1,mhlp
          iwal=lpnew_lwei(mpl)
          iwar=lpnew_rwei(mpl)
          lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
        enddo
        call ar_sd_ext_ar(51,lrd,0,isma)
        call ar_sd_ext_rest(lrd)
!dv(23-2) drl(33)-bl(23)-
        w0=w0_dv(2)
        if(ni.eq.1)w0=-w0_dv(2)
        do mpl=1,mtype
          vplp_w0(mpl)=vplpnew_w0(mpl)*w0
          vplp_w1(mpl)=0.d0
        enddo
        do lrk=1,lrd-1
          call drl_bl_ext_ar_new(6,lrk,lrd)
        enddo
      enddo
!=======================================================================
!dv(23) act -br-br-..................................................
232   if(linelp.ne.20) goto 233
      jk=nlg1
      do lrd=norb_frz+1,norb_dz
        lmd=lsm_inn(lrd)
        if(lmd.ne.jml) cycle
        ni=mod(norb_dz-lrd,2)
        w0=w0_dv(1)
        if(ni.eq.1)w0=-w0_dv(1)
!..............................................................
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
        call ar_br_br_ext_ar_new(6,intpos,isma)
      enddo
!dv(23) act -bl-bl-.....................................................
233   if(linelp.ne.22) return
      jk=nlg1
      do lrd=norb_frz+1,norb_dz
        lmd=lsm_inn(lrd)
        if(lmd.ne.jml) cycle
        ni=mod(norb_dz-lrd,2)
        w0=w0_dv(1)
        if(ni.eq.1)w0=-w0_dv(1)
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
        call ar_bl_bl_ext_ar_new(6,intpos,isma,1)
      enddo
      return
!=======================================================================
!d1v(24) act -c'-..................................................
!d1v(24-1) ar(13)-
!d1v(24-2) drl(33)-bl(13)-
124   if(linelp.eq.13) then
        call d1v_ar_act_c_ext(6)
      endif
      if(linelp.eq.20) then
        jk=nlg1
        call d1v_ar_act_brbr_ext(6,jk)
      endif
      if(linelp.eq.22) then
        jk=nlg1
        call d1v_ar_act_blbl_ext(6,jk)
      endif
      return
!=======================================================================
!ss(1)     act -b^l-
101   if(linelp.ne.15) return
      lra=nlg1
      if(nlg2.eq.1) then
        call ss_drl_act_bl(6,lra)
        if(jb_sys.gt.0) call ss_drl_act_bl_sgt0(6,lra)
      else
        call ss_arbl_act_bl(6,lra)
        if(jb_sys.gt.0) then
          call ss_arbl_act_bl_sgt0(6,lra)
          call ss_s_drl_act_bl_sgt0(6,lra)
        endif
      endif
      return
!=======================================================================
!st(2)      act -b^l
102   if(linelp.ne.15.or.nlg2.ne.2) return
      lra=nlg1
      call st_arbl_act_bl(6,lra)
      if(jml.eq.jmr) then
        call st_drl_act_bl(6,lra)
      endif
      if(jb_sys.gt.0) then
        call st_arbl_act_bl_sgt0(6,lra)
        if(jml.eq.jmr) then
          call st_drl_act_bl_sgt0(6,lra)
        endif
      endif
      return
!=======================================================================
!ts(3)   act -b^l ............................................
103   if(linelp.ne.15) return
      lra=nlg1
      if(nlg2.eq.2) then
        call ts_arbl_act_bl(6,lra)
      endif
      return
!=======================================================================
!stt(4) act-b^l ..............................................
104   if(linelp.ne.15) return
      lra=nlg1
      if(nlg2.eq.1) return
      call stt_arbl_act_bl_sgt1(6,lra)
      return
!=======================================================================
!tts(5) act-b^l ..............................................
105   if(linelp.ne.15.or.nlg2.ne.2) return
      lra=nlg1
      call tts_arbl_act_bl_sgt1(6,lra)
      call tts_drl_act_bl_sgt1(6,lra)
      return
!=======================================================================
!sdd(8) act-b^l -c'- -c"-..............................................
108   lra=nlg1
      if(linelp.eq.13) then
        call sdd_drlbl_act_c_sgt0(6)
        call sdd_drrbr_act_c_sgt0(6)
        call sdd_abb_act_c_sgt0(6)
        call sdd_ar_act_c_sd_ext_sgt0(6)
      else
        lra=nlg1
        if(linelp.eq.20) call sdd_ar_act_brbr_sgt0(6,lra)
        if(linelp.eq.22) call sdd_ar_act_blbl_sgt0(6,lra)
      endif
      return
!=======================================================================
!dds(9) act -c'- ..............................................
c109   if(linelp.eq.13) then
c        lra=nlg1
c        call dds_abb_act_c_sgt0(6)
c      endif
c      return
!=======================================================================
!tt(11)   act -b^l
111   if(linelp.ne.15) return
      lra=nlg1
      if(nlg2.eq.1) then
        call tt_drl_act_bl(6,lra)
      else
        call tt_arbl_act_bl(6,lra)
      endif
      return
!=======================================================================
!tttt(12)   act -b^l
112   if(linelp.ne.15) return
      lra=nlg1
      if(nlg2.eq.1) then
        call tttt_drl_act_bl_sgt1(6,lra)
      else
        call tttt_arbl_act_bl_sgt1(6,lra)
      endif
      return
!=======================================================================
!ttdd(15)    act -c'-   -br-br- -bl-bl
115   lra=nlg1
      if(linelp.eq.13) then                  ! no complete ar-
!        if(nlg2.eq.1) then
          call  ttdd_drlbl_act_c_sgt1(6)
          call ttdd_abb_act_c_sgt1(6)             ! ????
          call ttdd_ar_act_c_ttdd_ext_sgt1(6)
!        endif
      endif
      if(linelp.eq.20) then
        call ttdd_ar_act_brbr_sgt0(6,lra)
      endif
      if(linelp.eq.22) then
        call ttdd_ar_act_blbl_sgt0(6,lra)
      endif
      return
!=======================================================================
!ttv(19) act -b^r- ....................................................
118   if(linelp.ne.16.or.nlg2.ne.2) return
      lra=nlg1
      call ttv_arbr_act_c_sgt1(6,lra)
      return
!=======================================================================
!dd(19) act -b^l- ....................................................
119   if(linelp.ne.15) return
      lra=nlg1
      if(nlg2.eq.1) then
        call dd_drl_act_bl(6,lra)
      else
        call dd_arbl_act_bl(6,lra)
      endif
      return
!=======================================================================
!dddd(19) act -b^l- -br- ...............................................
120   lra=nlg1
      if(linelp.eq.15) then
        if(nlg2.eq.2) then
          call  dddd_arbl_act_bl_sgt0(6,lra)
        else
          call dddd_drl_act_bl_sgt0(6,lra)
        endif
      endif
      return
!=======================================================================
!dd1(21) act -b^l-  ....................................................
121   if(linelp.ne.15.or.nlg2.ne.2) return
      lra=nlg1
      call dd1_arbl_act_bl_sgt0(6,lra)
      return
!=======================================================================
!d1d(22) act -b^l-  ....................................................
122   lra=nlg1
      if(linelp.eq.15.and.nlg2.eq.2) then
        call d1d_arbl_act_bl_sgt0(6,lra)
        call d1d_drl_act_bl_sgt0(6,lra)
      endif
      return
!=======================================================================
!vv(25) act -b^l- ....................................................
125   if(linelp.ne.15) return
      lra=nlg1
      if(nlg2.eq.1) then
        iwdl=0
        iwdr=0
        do mpl=1,mtype
          vplp_w0(mpl)=vplpnew_w0(mpl)*w0_vv
          vplp_w1(mpl)=0.d0
        enddo
        do mpl=1,mhlp
          iwal=lpnew_lwei(mpl)
          iwar=lpnew_rwei(mpl)
          lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
        enddo
        call drl_bl_sum_ar_new(6,0,0,lra)
c        do lrk=1,norb_dz
c          call drl_bl_ext_ar_new(6,lrk,lra)
c        enddo
      endif
      return
!=======================================================================
!sv(10)     act -b^r-
!sv(10-1) ar(13)-br(23)-
!sv(10-2) ar(23)-br(13)-
!sv(10-3) drr(03)-
110   if(linelp.ne.16.or.nlg2.ne.2) return
      lra=nlg1
      if(jb_sys.gt.0) then
        call sv_arbr_act_br_sgt0(6,lra)
      endif
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        do lrj=lri,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jml) cycle
          w0sv2=w0_sv(2)
          w1sv2=w1_sv(2)
          ni=mod(lrj-lri,2)
          if(ni.eq.0) then
            w0sv2=-w0sv2
            w1sv2=-w1sv2
          endif
!-------------------------------------------------------------------
          iwdl=just(lri,lrj)
          iwdr=0
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo
         if(lri.ne.lrj) then
!sv(10-2) ar(23)-br(13)-
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sv2
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1sv2
            enddo
            ijk=lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz)
           intpos=intind_ijka(ijk)
            call ar_br_br_ext_ar_new(6,intpos,isma)
          else
!sv(10-3) drr(03)-
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0_sv(3)
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1_sv(3)
            enddo
            call drr_br_ext_ar(6,lri,lra)
          endif
        enddo
      enddo
      return
!=======================================================================
!tv(17)    act -b^r-
117   if(linelp.ne.16.or.nlg2.ne.2) return
      lra=nlg1
      do lri=norb_frz+1,norb_dz-1
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jml) cycle

          w1tv=w1_tv
          ni=mod(lrj-lri,2)
          if(ni.eq.0)  w1tv=-w1tv

!-------------------------------------------------------------------
!tv(17) ar(23)-br(23)-
          iwdl=just(lri,lrj)     !
          iwdr=0
          do mpl=1,mtype
            vplp_w0(mpl)=0.d0
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1tv
          enddo
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo
          ijk=lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz)
          intpos=intind_ijka(ijk)
          call ar_br_br_ext_ar_new(6,intpos,isma)
        enddo
      enddo
      return
!=======================================================================
10    return
      end


      subroutine sd_ext_head_in_act()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      logic_dh=.false.
      lri=nlg1
      lrj=nlg2
      intpos=nlg1
      isma=mul_tab(iml,imr)

      goto (1,2,3,4,5,6,7,8,9,10,11,12),linelp
!line=1 a&r<-->a^r
1     call ar_sd_ext_ar(100,lri,lrj,isma)
      call ar_sd_ext_rest(lri)
      goto 100
!line=4 a&r--b&r--b^r<-->a^r
4     call ar_br_br_ext_ar_new(6,intpos,isma)
      goto 100
!line=7 a&r--b&l--b^l<-->a^r
7     call ar_bl_bl_ext_ar_new(6,intpos,isma,1)
      goto 100
!line=10 d&rr--b^r<-->a^r
10    call drr_br_ext_ar(6,lri,lrj)
      goto 100
!line=12 d&rl--b^l<-->a^r
12    call drl_bl_ext_ar_new(6,lri,lrj)
      goto 100
2     goto 100
3     goto 100
5     goto 100
6     goto 100
8     goto 100
9     goto 100
11    goto 100
100   return
      end

      subroutine sd_ext_space_w01plp_value()
#include "drt_h.fh"
#include "lpextmode_h.fh"
      w0plp25=w0_sdplp*w0g25a
      w0plp26=w0_sdplp*w0g26a
      w0plp27=w0_sdplp*w0g27
      w1plp27=w0_sdplp*w1g27
      w0plp28=w0_sdplp*w0g28a
      w0plp29=w0_sdplp*w0g29
      w0plp30=w0_sdplp*w0g30
      w0plp31=w0_sdplp*w0g31
      w1plp31=w0_sdplp*w1g31
      w0plp32=w0_sdplp*w0g32
      w1plp32=w0_sdplp*w1g32

      end

      subroutine gsd_samesym_aaa(lri,isma)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei

      ia0=(lri-1)*norb_ext
      iabc0=(lri-1)*nabc

      ic=m_jd
        lrc=norb_number(ic)
      jcoffset=lrc*2-2

      iasta=ibsm_ext(isma)
      ibend=iesm_ext(isma)
      ibsta=iasta+1

      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)
      do ib=ibsta,ic-1
        do ia=iasta,ib-1
!     g31   type_12 arbrb^ra^r
          iabc=iabc0+ia+ngw2(ib)+ngw3(ic)
          iposint=intind_iabc(iabc)
          value_lpext(ilwei)=vint_ci(iposint+1)*w0plp31
     *               +vint_ci(iposint+2)*w1plp31
         ilwei=ilwei+1
      enddo
      enddo

      ib=ic
      ilwei=icnt_base+iwt_orb_ext(iasta,ib)
      do ia=iasta,ib-1
!        g28     cw-ar                    330
         iaqq=ia0+ia
         intpos=intind_iaqq(iaqq)
         iposint=intpos+jcoffset
          value_lpext(ilwei)=(vint_ci(iposint)/w0g28a
     :             +vint_ci(iposint+1))*w0plp28
         ilwei=ilwei+1
      enddo

      ia=ic
      do ib=ic+1,ibend
!        g25,g27: bl(20)-drl(11)         220
         iaqq=ia0+ib
         intpos=intind_iaqq(iaqq)
         iposint=intpos+jcoffset
         ilwei=icnt_base+iwt_orb_ext(ic,ib)
         value_lpext(ilwei)=vint_ci(iposint)*w0plp27
     :                    -vint_ci(iposint+1)*w1plp27
      enddo

      do ib=ic+1,ibend
        ilwei=icnt_base+iwt_orb_ext(iasta,ib)
        do ia=iasta,ic-1
!     g32a   type g12
          iabc=iabc0+ia+ngw2(ic)+ngw3(ib)
          iposint=intind_iabc(iabc)
          value_lpext(ilwei)=vint_ci(iposint+2)*w0plp32   !severe_new_er
     *                  -vint_ci(iposint)*w1plp32
          ilwei=ilwei+1
        enddo
      enddo

      do ib=ic+2,ibend
        ilwei=icnt_base+iwt_orb_ext(ic+1,ib)
        do ia=ic+1,ib-1
          iabc=iabc0+ic+ngw2(ia)+ngw3(ib)
          iposint=intind_iabc(iabc)
!     g32b  type 11
         value_lpext(ilwei)=vint_ci(iposint+1)*w0plp32      !severe_new_
     *                  -vint_ci(iposint)*w1plp32
         ilwei=ilwei+1
        enddo
      enddo
      end

      subroutine gsd_diffsamesym_abb(lri,isma,ismb)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei

      ia0=(lri-1)*norb_ext
      iabc0=(lri-1)*nabc

      ic=m_jd
        lrc=norb_number(ic)
      jcoffset=lrc*2-2

      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)
      ibsta=ibsm_ext(ismb)
      ibend=iesm_ext(ismb)

      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)
      do ib=ibsta,ic-1
        do ia=iasta,iaend
         !g31
!     g31   type_12 arbrb^ra^r
          iabc=iabc0+ia+ngw2(ib)+ngw3(ic)
          iposint=intind_iabc(iabc)
         value_lpext(ilwei)=vint_ci(iposint+1)*w0plp31      !severe_new_
     *               +vint_ci(iposint+2)*w1plp31
         ilwei=ilwei+1
        enddo
      enddo

      ilwei=icnt_base+iwt_orb_ext(iasta,ic+1)
      jb=m_jc
      do ib=ic+1,ibend
         jb=jb+1
        do ia=iasta,iaend
         !g32a
          iabc=iabc0+ia+ngw2(ic)+ngw3(ib)
          iposint=intind_iabc(iabc)
         value_lpext(ilwei)=vint_ci(iposint+2)*w0plp32   !severe_new_err
     *               -vint_ci(iposint)*w1plp32
         ilwei=ilwei+1
        enddo
      enddo

      ib=ic
      ilwei=icnt_base+iwt_orb_ext(iasta,ib)
      do ia=iasta,iaend      !ib-1      !severe_error_1020
!        g28
         iaqq=ia0+ia
         intpos=intind_iaqq(iaqq)
         iposint=intpos+jcoffset
c        value_lpext(ilwei)=(vint_ci(iposint2)+
         value_lpext(ilwei)=w0plp28*(vint_ci(iposint)/w0g28a+
     *                              vint_ci(iposint+1))
         ilwei=ilwei+1
      enddo
      end

      subroutine gsd_diffsamesym_aab(lri,isma,ismb)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei

      ia0=(lri-1)*norb_ext
      iabc0=(lri-1)*nabc

      ic=m_jd
        lrc=norb_number(ic)
      jcoffset=lrc*2-2

      ibsta=ibsm_ext(ismb)
      ibend=iesm_ext(ismb)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)

      !   if ( ic-1.gt.iasta) then
      do ib=ibsta,ibend
         ilwei=icnt_base+iwt_orb_ext(iasta,ib)
        do ia=iasta,ic-1
         !g32a
          iabc=iabc0+ia+ngw2(ic)+ngw3(ib)
          iposint=intind_iabc(iabc)
         value_lpext(ilwei)=vint_ci(iposint+2)*w0plp32      !severe_new_
     *                  -vint_ci(iposint)*w1plp32
         ilwei=ilwei+1
        enddo
      enddo

      do ib=ibsta,ibend
         ilwei=icnt_base+iwt_orb_ext(ic+1,ib)
        do ia=ic+1,iaend      !ib-1      !severe_error_1020
          iabc=iabc0+ic+ngw2(ia)+ngw3(ib)
          iposint=intind_iabc(iabc)
         !g32b
         value_lpext(ilwei)=vint_ci(iposint+1)*w0plp32   !severe_new_err
     *                  -vint_ci(iposint)*w1plp32
         ilwei=ilwei+1
        enddo
      enddo

      ia=ic
      do ib=ibsta,ibend
!        g25
         iaqq=ia0+ib
         intpos=intind_iaqq(iaqq)
         iposint=intpos+jcoffset
         ilwei=icnt_base+iwt_orb_ext(ic,ib)
         value_lpext(ilwei)=vint_ci(iposint)*w0plp27
     :                    -vint_ci(iposint+1)*w1plp27
      enddo
      end

      subroutine gsd_arlp_s1(lri)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei

      ia0=(lri-1)*norb_ext
      ic=m_jd

      iaqq=ia0+ic
      intpos=intind_iaqq(iaqq)

      ilwei=icnt_base+isegdownwei-norb_ext+1
      do is1orb=1,ic-1
!g30 -b^rd^rr
         lrk=norb_number(is1orb)
         intoffset=(lrk-1)*2
         iposint=intpos+intoffset
         value_lpext(ilwei)=vint_ci(iposint)*w0plp30
         ilwei=ilwei+1
      enddo

!g26 -a^r     610
      lrk=norb_number(ic)
      intoffset=(lrk-1)*2
      iposint=intpos+intoffset
      value_lpext(ilwei)=vint_ci(iposint)*w0plp26
      ilwei=ilwei+1

!g29 -dl^ra^l
      do is1orb=ic+1,norb_ext
      lrk=norb_number(is1orb)
      intoffset=(lrk-1)*2
      iposint=intpos+intoffset
        value_lpext(ilwei)=vint_ci(iposint)*w0plp29
        ilwei=ilwei+1
      enddo
      end

      subroutine td_drt_ci_new()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      common/lpdisk/idisk_lp,idisk_array(13)
      data dzero/0.d0/dsq2/1.414213562373d0/
      data dsq3/ 1.732050807569d0/
      data dsq3vsq2/1.224744871392d0/
c      write(6,*)'  td_wyb'

!      logic_sd=.true.
      call external_space_plpmode_value_td()
      w0_sdplp=1.0d0
      call sd_ext_space_w01plp_value()

      idisk_lp=idisk_array(7)
      do lpb=1,lpblock_td
        call read_lp()
        ipael=iml+9
        ipae =imr+1
        call get_jpty(jpadlr,jptyl,jptyr)
        call get_jp(jptyl,jml,jpadl,1)
        call get_jp(jptyr,jmr,jpad,1)
c        jmlr=mul_tab(jml,jmr)
        call gsd_determine_extarmode_paras(iml,imr,.false.)
        if(linelp.le.12)   then
          call td_ext_head_in_act()
        else
          call td_ext_head_in_dbl()
        endif
      enddo
      return
      end

      subroutine td_ext_head_in_dbl()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      logic_dh=.true.
      isma=mul_tab(iml,imr)
      jmlr=mul_tab(jml,jmr)
      lpok=jpadlr
!      goto(101,102,103,104,105,106,10,108,109,110,111,112,
!     :     113,10,115,10,117,118,119,120,121,122,123,124,125,10),lpok
      goto(101,102,103,104,105,106,10,108,109,110,111,112,
     :     113,10,115,10,117,118,119,120,121,122,123,124,125,10),lpok
!====================================================================
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
106   if(linelp.ne.13) goto 107
      if(jb_sys.gt.0) then
        call sd_adb_act_c_ext_ar(13)
      endif
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jmlr) cycle
        w0sd1=w0_sd(1)
        w0sd9=w0_sd(9)
        w1sd9=w1_sd(9)
        w0sd12=w0_sd(12)
        w1sd12=w1_sd(12)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1) w0sd1=-w0sd1
        if(ni.eq.1) w0sd9=-w0sd9
        if(ni.eq.1) w1sd9=-w1sd9
        if(ni.eq.1) w0sd12=-w0sd12
        if(ni.eq.1) w1sd12=-w1sd12
        if(jml.eq.1.and.lmi.eq.jmr) then
          iwdl=just(lri,lri)
          iwdr=jud(lri)
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo

!sd(6-1) a&r(02)-
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd1
            vplp_w1(mpl)=vplpnew_w1(mpl)*w0sd1
          enddo
          call ar_td_ext_ar(26,lri,lrj,isma)
          call ar_td_ext_rest(lri)

!sd(6-12) d&r&l(33)b^l(02)
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd12
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1sd12
          enddo
          do lrk=1,lri-1
            call drl_bl_ext_ar_new(13,lrk,lri)
          enddo

!sd(6-9) d&r&r(03)b^r(32)
          do lrk=norb_frz+1,lri-1
            iwdl=just(lrk,lrk)
            iwdr=jud(lri)
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd9
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1sd9
            enddo
            call drr_br_ext_ar(13,lrk,lri)
          enddo
        endif
      enddo

      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jml) cycle

         if(lmi.eq.jmr) then
            iwdl=just(lri,lrj)
            iwdr=jud(lri)
            w0sd2=w0_sd(2)
            w0sd11=w0_sd(11)
            w1sd11=w1_sd(11)
            w0sd14=w0_sd(14)
            ni=mod(norb_dz-lrj,2)
            if(ni.eq.1) then
              w0sd2=-w0sd2
              w0sd11=-w0sd11
              w1sd11=-w1sd11
              w0sd14=-w0sd14
            endif
!sd(6-2) c(22)-a&r(13)-
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd2
            enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            call ar_td_ext_ar(25,lrj,lri,isma)
            call ar_td_ext_rest(lrj)
!sd(6-11) d&r&l(22)b^l(13)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd11
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1sd11
            enddo
            call drl_bl_ext_ar_new(13,lri,lrj)
!sd(6-14) (22)d&r&l(33)b^l(13)
!sd(6-14) d&r&l(33)c"(22)b^l(13)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd14
              vplp_w1(mpl)=0.d0
            enddo
            do lrk=1,lrj-1
              if(lrk.eq.lri) cycle
              call drl_bl_ext_ar_new(13,lrk,lrj)
            enddo
          endif
!sd(6-4) a&r(23)-c'(12)-
          w0sd4=w0_sd(4)
          w0sd16=w0_sd(16)
          ni=mod(norb_dz-lri,2)
          if(ni.eq.0) then
            w0sd4=-w0sd4
            w0sd16=-w0sd16
          endif
          if(lmj.eq.jmr) then
            iwdl=just(lri,lrj)
            iwdr=jud(lrj)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd4
 !             vplp_w1(mpl)=vplpnew_w1(mpl)*w0sd4
            enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            call ar_td_ext_ar(28,lri,lrj,isma)
            call ar_td_ext_rest(lri)
!sd(6-16) d&r&l(33)b^l(23)c'(12)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd16
              vplp_w1(mpl)=0.d0
            enddo
            do lrk=1,lri-1
              call drl_bl_ext_ar_new(13,lrk,lri)
            enddo
          endif
!sd(6-5) a&r(23)b&r(13)b^r(32)
          do lrd=lrj+1,norb_dz
            lmd=lsm_inn(lrd)
           if(lmd.ne.jmr) cycle

            w0sd5=w0_sd(5)
            w1sd5=w1_sd(5)
            ni=mod(lrj-lri+norb_dz-lrd,2)
            if(ni.eq.0) then
              w0sd5=-w0sd5
              w1sd5=-w1sd5
            endif
            iwdl=just(lri,lrj)
            iwdr=jud(lrd)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd5
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1sd5
            enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            ijk=lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lrd-norb_frz)
            intpos=intind_ijka(ijk)
            call ar_br_br_ext_ar_new(13,intpos,isma)
          enddo
!sd(6-8) a&r(23)b&l(32)b^l(13)
          do lrd=lri+1,lrj-1
            lmd=lsm_inn(lrd)
            if(lmd.ne.jmr) cycle
            iwdl=just(lri,lrj)
            iwdr=jud(lrd)
            w0sd8=w0_sd(8)
            w1sd8=w1_sd(8)
            ni=mod(lrj-lri+norb_dz-lrd,2)
            if(ni.eq.0) then
              w0sd8=-w0sd8
              w1sd8=-w1sd8
            endif
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sd8
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1sd8
            enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            ijk=lri-norb_frz+ngw2(lrd-norb_frz)+ngw3(lrj-norb_frz)
            intpos=intind_ijka(ijk)
            call ar_bl_bl_ext_ar_new(13,intpos,isma,1)
          enddo
        enddo
      enddo
!====================================================================
!sd(6)    act -br-br-
107   if(linelp.eq.20) then
        jk=nlg1
        call sd_ar_act_brbr(13,jk)
        if(jb_sys.gt.0) call sd_ar_act_brbr_sgt0(13,jk)
      endif
!sd(6)  act -bl-bl-
      if(linelp.eq.22) then
        jk=nlg1
        call sd_ar_act_blbl(13,jk)
        if(jb_sys.gt.0) call sd_ar_act_blbl_sgt0(13,jk)
      endif
      return
!=======================================================================
!td(13)                                  act -c'-
!td(13-1) (22)a&(23)
!td(13-1) a&(23)c'(22)
!td(13-2) a&(23)b&r(23)b^r(32)
!td(13-3) a&(23)b&l(32)b^l(23)
!td(13-4) d&r&l(22)b^l(23)
!td(13-5) (22)d&&l(33)b^l(23)
!td(13-5) d&rl(33)c"(22)b^l(23)
!td(13-5) d&rl(33)b^l(23)c'(22)
113   if(linelp.ne.13) goto 132
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jml) cycle
          w0td1=w0_td(1)
          w0td4=w0_td(4)
          w1td4=w1_td(4)
          w0td5=w0_td(5)
          ni=mod(norb_dz-lrj,2)
          if(ni.eq.1) then
            w0td1=-w0td1
            w0td4=-w0td4
            w1td4=-w1td4
            w0td5=-w0td5
          endif
          iwdl=just(lri,lrj)      !
!-------------------------------------------------------------------
!td(13-1) (22)a&(23)
          if(lmi.eq.jmr) then
            iwdr=jud(lri)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0td1
            enddo
              do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            call ar_td_ext_ar(43,lrj,lri,isma)
            call ar_td_ext_rest(lrj)
!td(13-4) d&r&l(22)b^l(23)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0td4
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1td4
            enddo
            call drl_bl_ext_ar_new(13,lri,lrj)
!td(13-5) (22)d&rl(33)b^l(23)
!td(13-5) d&rl(33)c"(22)b^l(23)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0td5
              vplp_w1(mpl)=0.d0
            enddo
            do lrk=1,lrj-1
              if(lrk.eq.lri) cycle
              call drl_bl_ext_ar_new(13,lrk,lrj)
            enddo
          endif
!-------------------------------------------------------------------
!td(13-1) a&(23)c'(22)
          if(lmj.eq.jmr) then
            iwdr=jud(lrj)
            w0td1=w0_td(1)
            w0td5=w0_td(5)
            ni=mod(norb_dz-lri,2)
            if(ni.eq.0) then
              w0td1=-w0td1
              w0td5=-w0td5
            endif
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0td1
            enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            call ar_td_ext_ar(46,lri,lrj,isma)
            call ar_td_ext_rest(lri)
!td(13-5) d&rl(33)b^l(23)c'(22)
              do mpl=1,mtype
                vplp_w0(mpl)=vplpnew_w0(mpl)*w0td5
                vplp_w1(mpl)=0.d0
              enddo
              do lrk=1,lri-1
                call drl_bl_ext_ar_new(13,lrk,lri)
              enddo
            endif
!-------------------------------------------------------------------
!td(13-2) a&(23)b&r(23)b^r(32)
         do lrd=lrj+1,norb_dz
            lmd=lsm_inn(lrd)
            if(lmd.ne.jmr) cycle
            iwdr=jud(lrd)
            ijk=lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lrd-norb_frz)
            intpos=intind_ijka(ijk)
            w0td2=w0_td(2)
            w1td2=w1_td(2)
            ni=mod(lrj-lri+norb_dz-lrd,2)
            if(ni.eq.0) then
              w0td2=-w0td2
              w1td2=-w1td2
            endif
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0td2
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1td2
            enddo
            do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            call ar_br_br_ext_ar_new(13,intpos,isma)
          enddo
!-------------------------------------------------------------------
!td(13-3) a&(23)b&l(32)b^l(23)
          do lrd=lri+1,lrj-1
            lmd=lsm_inn(lrd)
            if(lmd.ne.jmr) cycle
            iwdr=jud(lrd)
            w0td3=w0_td(3)
            w1td3=w1_td(3)
            ni=mod(lrj-lri+norb_dz-lrd,2)
            if(ni.eq.0) then
              w0td3=-w0td3
              w1td3=-w1td3
            endif
            ijk=lri-norb_frz+ngw2(lrd-norb_frz)+ngw3(lrj-norb_frz)
            intpos=intind_ijka(ijk)
              do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0td3
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1td3
            enddo
              do mpl=1,mhlp
              iwal=lpnew_lwei(mpl)
              iwar=lpnew_rwei(mpl)
              lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
              lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
            enddo
            call ar_bl_bl_ext_ar_new(13,intpos,isma,1)
          enddo
!-------------------------------------------------------------------
        enddo
      enddo
!=======================================================================
!td(13) act -br-br-
132   if(linelp.eq.20) then
        jk=nlg1
        call td_ar_act_brbr(13,jk)
      endif
!td(13)  act -bl-bl-
      if(linelp.eq.22) then
        jk=nlg1
        call td_ar_act_blbl(13,jk)
      endif
      return
!=======================================================================
!dv(23) act -c'-..................................................
!dv(23-1) ar(23)-
!dv(23-2) drl(33)-bl(23)-
123   if(linelp.ne.13) goto 232
      imap_1=0
      do lrd=norb_frz+1,norb_dz
        lmd=lsm_inn(lrd)
        if(lmd.ne.jml) cycle
        ni=mod(norb_dz-lrd,2)

!dv(23-1) ar(23)-
        w0=w0_dv(1)
        if(ni.eq.1)w0=-w0_dv(1)
          iwdl=jud(lrd)
          iwdr=0
          imap_1=imap_1+1
          do mpl=1,mtype
            vplp_w0(mpl)=vplpnew_w0(mpl)*w0
          enddo
        do mpl=1,mhlp
          iwal=lpnew_lwei(mpl)
          iwar=lpnew_rwei(mpl)
          lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
        enddo
        call ar_td_ext_ar(51,lrd,0,isma)
        call ar_td_ext_rest(lrd)
!dv(23-2) drl(33)-bl(23)-
        w0=w0_dv(2)
        if(ni.eq.1)w0=-w0_dv(2)
        do mpl=1,mtype
          vplp_w0(mpl)=vplpnew_w0(mpl)*w0
          vplp_w1(mpl)=0.d0
        enddo
        do lrk=1,lrd-1
          call drl_bl_ext_ar_new(13,lrk,lrd)
        enddo
      enddo
!=======================================================================
!dv(23) act -br-br-..................................................
232   if(linelp.ne.20) goto 233
      jk=nlg1
      do lrd=norb_frz+1,norb_dz
        lmd=lsm_inn(lrd)
        if(lmd.ne.jml) cycle
        ni=mod(norb_dz-lrd,2)
        w0=w0_dv(1)
        if(ni.eq.1)w0=-w0_dv(1)
!..............................................................
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
        call ar_br_br_ext_ar_new(13,intpos,isma)
      enddo
!dv(23) act -bl-bl-.....................................................
233   if(linelp.ne.22) return
      jk=nlg1
      do lrd=norb_frz+1,norb_dz
        lmd=lsm_inn(lrd)
        if(lmd.ne.jml) cycle
        ni=mod(norb_dz-lrd,2)
        w0=w0_dv(1)
        if(ni.eq.1)w0=-w0_dv(1)
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
        call ar_bl_bl_ext_ar_new(13,intpos,isma,1)
      enddo
      return
!=======================================================================
!ss(1)     act -b^l-
101   if(linelp.ne.15) return
      lra=nlg1
      if(nlg2.eq.1) then
        call ss_drl_act_bl(13,lra)
        if(jb_sys.gt.0) call ss_drl_act_bl_sgt0(13,lra)
      else
        call ss_arbl_act_bl(13,lra)
        if(jb_sys.gt.0) then
          call ss_arbl_act_bl_sgt0(13,lra)
          call ss_s_drl_act_bl_sgt0(13,lra)
        endif
      endif
      return
!=======================================================================
!st(2)      act -b^l
102   if(linelp.ne.15.or.nlg2.ne.2) return
      lra=nlg1
      call st_arbl_act_bl(13,lra)
      if(jb_sys.gt.0) call st_arbl_act_bl_sgt0(13,lra)
      if(jml.eq.jmr) then
        call st_drl_act_bl(13,lra)
        if(jb_sys.gt.0) call st_drl_act_bl_sgt0(13,lra)
      endif
      return
!=======================================================================
!ts(3)   act -b^l ............................................
103   if(linelp.ne.15) return
      lra=nlg1
      if(nlg2.eq.2) then
        call ts_arbl_act_bl(13,lra)
        if(jb_sys.gt.0) call ts_arbl_act_bl_sgt0(13,lra)
      endif
      return
!=======================================================================
!stt(4) act-b^l ..............................................
104   if(linelp.ne.15) return
      lra=nlg1
      if(nlg2.eq.1) return
      call stt_arbl_act_bl_sgt1(13,lra)
      return
!=======================================================================
!tts(5) act-b^l ..............................................
105   if(linelp.ne.15.or.nlg2.ne.2) return
      lra=nlg1
      call tts_arbl_act_bl_sgt1(13,lra)
      call tts_drl_act_bl_sgt1(13,lra)
      return
!=======================================================================
!sdd(8) act-b^l -c'- -c"-..............................................
108   if(linelp.eq.13) then
        lra=nlg1
        call sdd_drlbl_act_c_sgt0(13)
        call sdd_drrbr_act_c_sgt0(13)
        call sdd_abb_act_c_sgt0(13)
        call sdd_ar_act_c_sd_ext_sgt0(13)
      else
        lra=nlg1
        if(linelp.eq.20) call sdd_ar_act_brbr_sgt0(13,lra)
        if(linelp.eq.22) call sdd_ar_act_blbl_sgt0(13,lra)
      endif
      return
!=======================================================================
!dds(9) act -c'- ..............................................
109   if(linelp.eq.13) then
        lra=nlg1
        call dds_abb_act_c_sgt0(13)
      endif
      return
!=======================================================================
!tt(11)   act -b^l
111   if(linelp.ne.15) return
      lra=nlg1
      if(nlg2.eq.1) then
        call tt_drl_act_bl(13,lra)
      else
        call tt_arbl_act_bl(13,lra)
      endif
      return
!=======================================================================
!tttt(12)   act -b^l
112   if(linelp.ne.15) return
      lra=nlg1
      if(nlg2.eq.1) then
        call tttt_drl_act_bl_sgt1(13,lra)
      else
        call tttt_arbl_act_bl_sgt1(13,lra)
      endif
      return
!=======================================================================
!ttdd(15)    act -c'-   -br-br- -bl-bl
115   lra=nlg1
      if(linelp.eq.13) then                  ! no complete ar-
!        if(nlg2.eq.1) then
          call ttdd_drlbl_act_c_sgt1(13)
          call ttdd_abb_act_c_sgt1(13)             ! ????
          call ttdd_ar_act_c_ttdd_ext_sgt1(13)
!        endif
      endif
      if(linelp.eq.20) then
        call ttdd_ar_act_brbr_sgt0(13,lra)
      endif
      if(linelp.eq.22) then
        call ttdd_ar_act_blbl_sgt0(13,lra)
      endif
      return
!======================================================
!t1v(18) arbr act -c"- ext -br-ar -ar  ! no complete
118   if(nlg2.ne.2.or.linelp.ne.16) return
      lra=nlg1
      call ttv_arbr_act_c_sgt1(13,lra)
      return
!=======================================================================
!dd(19) act -b^l- ....................................................
119   if(linelp.ne.15) return
      lra=nlg1
      if(nlg2.eq.1) then
        call dd_drl_act_bl(13,lra)
      else
        call dd_arbl_act_bl(13,lra)
      endif
      return
!=======================================================================
!dddd(20) act -b^l-  ...................................................
120   if(linelp.ne.15) return
      lra=nlg1
      if(nlg2.eq.1) then
        call dddd_drl_act_bl_sgt0(13,lra)
      else
        call  dddd_arbl_act_bl_sgt0(13,lra)
      endif
      return
!=======================================================================
!dd1(21) act -b^l-  ....................................................
121   if(linelp.ne.15.or.nlg2.ne.2) return
      lra=nlg1
      call dd1_arbl_act_bl_sgt0(13,lra)
      return
!=======================================================================
!d1d(22) act -b^l-  ....................................................
122   lra=nlg1
      if(linelp.eq.15.and.nlg2.eq.2) then
        call d1d_arbl_act_bl_sgt0(13,lra)
        call d1d_drl_act_bl_sgt0(13,lra)
      endif
      return
!=======================================================================
!d1v(24) act -c'-..................................................
!d1v(24-1) ar(13)-
!d1v(24-2) drl(33)-bl(13)-
124   if(linelp.eq.13) then
        call d1v_ar_act_c_ext(13)
      endif
      if(linelp.eq.20) then
        jk=nlg1
        call d1v_ar_act_brbr_ext(13,jk)
      endif
      if(linelp.eq.22) then
        jk=nlg1
        call d1v_ar_act_blbl_ext(13,jk)
      endif
      return
!=======================================================================
!vv(25) act -b^l- ....................................................
125   if(linelp.ne.15) return
      lra=nlg1
      if(nlg2.eq.1) then
        iwdl=0
        iwdr=0
        do mpl=1,mtype
          vplp_w0(mpl)=vplpnew_w0(mpl)*w0_vv
          vplp_w1(mpl)=0.d0
        enddo
        do mpl=1,mhlp
          iwal=lpnew_lwei(mpl)
          iwar=lpnew_rwei(mpl)
          lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
          lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
        enddo
        call drl_bl_sum_ar_new(13,0,0,lra)
c        do lrk=1,norb_dz
c          call drl_bl_ext_ar_new(13,lrk,lra)
c        enddo
      endif
      return
!=======================================================================
!sv(10)     act -b^r-
!sv(10-1) ar(13)-br(23)-
!sv(10-2) ar(23)-br(13)-
!sv(10-3) drr(03)-
110   if(linelp.ne.16.or.nlg2.ne.2) return
      lra=nlg1
      if(jb_sys.gt.0) then
        call sv_arbr_act_br_sgt0(13,lra)
      endif
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        do lrj=lri,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jml) cycle
          w0sv2=w0_sv(2)
          w1sv2=w1_sv(2)
          ni=mod(lrj-lri,2)
          if(ni.eq.0) then
            w0sv2=-w0sv2
            w1sv2=-w1sv2
          endif
!-------------------------------------------------------------------
          iwdl=just(lri,lrj)
          iwdr=0
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo
         if(lri.ne.lrj) then
!sv(10-2) ar(23)-br(13)-
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sv2
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1sv2
            enddo
            ijk=lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz)
           intpos=intind_ijka(ijk)
            call ar_br_br_ext_ar_new(13,intpos,isma)
          else
!sv(10-3) drr(03)-
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0_sv(3)
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1_sv(3)
            enddo
            call drr_br_ext_ar(13,lri,lra)
          endif
        enddo
      enddo
      return
!=======================================================================
!tv(17)    act -b^r-
117   if(linelp.ne.16.or.nlg2.ne.2) return
      lra=nlg1
      do lri=norb_frz+1,norb_dz-1
        lmi=lsm_inn(lri)
        do lrj=lri+1,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jml) cycle

          w1tv=w1_tv
          ni=mod(lrj-lri,2)
          if(ni.eq.0)  w1tv=-w1tv

!-------------------------------------------------------------------
!tv(17) ar(23)-br(23)-
          iwdl=just(lri,lrj)      !
          iwdr=0
          do mpl=1,mtype
            vplp_w0(mpl)=0.d0
            vplp_w1(mpl)=vplpnew_w1(mpl)*w1tv
          enddo
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo
          ijk=lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lra-norb_frz)
          intpos=intind_ijka(ijk)
          call ar_br_br_ext_ar_new(13,intpos,isma)
        enddo
      enddo
      return
!=======================================================================
10    return
      end


      subroutine td_ext_head_in_act()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      logic_dh=.false.
      lri=nlg1
      lrj=nlg2
      intpos=nlg1
      isma=mul_tab(iml,imr)

      goto (1,2,3,4,5,6,7,8,9,10,11,12),linelp
!line=1 a&r<-->a^r
1     call ar_td_ext_ar(100,lri,lrj,isma)
      call ar_td_ext_rest(lri)
      goto 100
!line=4 a&r--b&r--b^r<-->a^r
4     call ar_br_br_ext_ar_new(13,intpos,isma)
      goto 100
!line=7 a&r--b&l--b^l<-->a^r
7     call ar_bl_bl_ext_ar_new(13,intpos,isma,1)
      goto 100
!line=10 d&rr--b^r<-->a^r
10    call drr_br_ext_ar(13,lri,lrj)
      goto 100
!line=12 d&rl--b^l<-->a^r
12    call drl_bl_ext_ar_new(13,lri,lrj)
      goto 100
2     goto 100
3     goto 100
5     goto 100
6     goto 100
8     goto 100
9     goto 100
11    goto 100
100   return
      end


      subroutine gsd_ext_sequence(iltype,ilsm,irsm,lri)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *      m_jc,m_jd, isegsta,isegupwei,isegdownwei
      ismnodes=ilsm
      ismnoded=irsm
      indl=0 !?
      if(iltype.eq.2)indl= 1+ismnodes
      if(iltype.eq.3)indl= 9+ismnodes
      if(iltype.eq.4)indl=17+ismnodes
      ilnodedownwei=iseg_downwei(indl)
      isegdownwei  =ilnodedownwei
      icano_nnsta=1
      icnt_base=0
      icsta=ibsm_ext(ismnoded)
      icend=iesm_ext(ismnoded)
      m_jc=0

      do ic=icsta,icend
            m_jd=ic
            m_jc=ic-icsta+1
            icano_nn=m_jc
            icano_nnend=icano_nn
            do ismb=1,ismnoded-1
                  isma=mul_tab(ismnodes,ismb)
                  if ( isma .gt. ismb ) cycle
                  call g31_diffsym(lri,isma,ismb)
            enddo

            ismb=ismnoded
            isma=mul_tab(ismnodes,ismb)
            if ( isma .eq. ismb ) then
               call gsd_samesym_aaa(lri,isma)
            elseif ( isma .lt. ismb ) then
               call gsd_diffsamesym_abb(lri,isma,ismb)
            endif

            do ismb=ismnoded+1,ng_sm
                  isma=mul_tab(ismnodes,ismb)
                  if ( isma .gt. ismb ) cycle
                  if ( ismnoded .gt. isma ) then
                     call g32a_diffsym(lri,isma,ismb)
                  elseif ( ismnoded .eq. isma ) then
                    call gsd_diffsamesym_aab(lri,isma,ismb)
                  else
                     call g32b_diffsym(lri,isma,ismb)
                  endif
            enddo

            if ( ismnodes.eq.1 .and. iltype.eq.4 ) then
                  call gsd_arlp_s1(lri)
            endif
            icnt_base=icnt_base+ilnodedownwei

      enddo
      end

      subroutine ar_sd_ext_rest(lri)
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *      m_jc,m_jd, isegsta,isegupwei,isegdownwei

      iwuplwei=jpad_upwei(jpadl)
      ilsegdownwei=iseg_downwei(ipael)
      irsegdownwei=iseg_downwei(ipae)

      if(logic_grad) then
        call gsd_ext_sequence_g(4,iml,imr,lri)

        w0_old=1.0d0
        do iw0=1,mtype
          w0_sdplp=vplpnew_w0(iw0)
          if(logic_dh) w0_sdplp=vplp_w0(iw0)
          w0multi=w0_sdplp/w0_old
          w0_old=w0_sdplp
          do iiext=1,icnt_base
            value_lpext(iiext)=value_lpext(iiext)*w0multi
            value_lpext1(iiext)=value_lpext1(iiext)*w0multi
          enddo
          ilpsta=nstaval(iw0)+1
          ilpend=nstaval(iw0)+nvalue(iw0)
          do iplp=ilpsta,ilpend
            if(logic_dh) then              !lp_head is in dbl_space
              ilw=lp_lwei(iplp)
              irw=lp_rwei(iplp)
              call complete_sd_ar_ext_loop_g(ilw,irw,ilsegdownwei)
            else
              ihypos=jphy(iplp)
              ndim =ihy(ihypos)
              iwal0=lpnew_lwei(iplp)
              iwar0=lpnew_rwei(iplp)
              do in=1,ndim
                iwal=iwal0+ihyl(ihypos+in)
                iwar=iwar0+ihy(ihypos+in)
                do iwd=0,iwuplwei-1
                  ilw=iwalk_ad(jpadl,ipael,iwal,iwd)
                  irw=iwalk_ad(jpad,ipae,iwar,iwd)
                  call complete_sd_ar_ext_loop_g(ilw,irw,ilsegdownwei)
                enddo
              enddo
            endif
          enddo
        enddo

      else
        call gsd_ext_sequence(4,iml,imr,lri)

        w0_old=1.0d0
        do iw0=1,mtype
          w0_sdplp=vplpnew_w0(iw0)
          if(logic_dh) w0_sdplp=vplp_w0(iw0)
          w0multi=w0_sdplp/w0_old
          w0_old=w0_sdplp
          do iiext=1,icnt_base
            value_lpext(iiext)=value_lpext(iiext)*w0multi
          enddo
          ilpsta=nstaval(iw0)+1
          ilpend=nstaval(iw0)+nvalue(iw0)
          do iplp=ilpsta,ilpend
            if(logic_dh) then              !lp_head is in dbl_space
              ilw=lp_lwei(iplp)
              irw=lp_rwei(iplp)
              call complete_sd_ar_ext_loop(ilw,irw,ilsegdownwei)
            else                                    !lp_head is in act_s
              ihypos=jphy(iplp)
              ndim =ihy(ihypos)
              iwal0=lpnew_lwei(iplp)
              iwar0=lpnew_rwei(iplp)
              do in=1,ndim
                iwal=iwal0+ihyl(ihypos+in)
                iwar=iwar0+ihy(ihypos+in)
                do iwd=0,iwuplwei-1
                  ilw=iwalk_ad(jpadl,ipael,iwal,iwd)
                  irw=iwalk_ad(jpad,ipae,iwar,iwd)
                  call complete_sd_ar_ext_loop(ilw,irw,ilsegdownwei)
                enddo
              enddo
            endif
          enddo
        enddo
      endif
      return
      end

      subroutine ar_td_ext_rest(lri)
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *      m_jc,m_jd, isegsta,isegupwei,isegdownwei

      iwuplwei=jpad_upwei(jpadl)
      ilsegdownwei=iseg_downwei(ipael)
      irsegdownwei=iseg_downwei(ipae)

      if(logic_grad) then
      call gsd_ext_sequence_g(3,iml,imr,lri)

      w0_old=1.0d0
      do iw0=1,mtype
        w0_sdplp=vplpnew_w0(iw0)
        if(logic_dh) w0_sdplp=vplp_w0(iw0)
        w0multi=w0_sdplp/w0_old
        w0_old=w0_sdplp
        do iiext=1,icnt_base
          value_lpext(iiext)=value_lpext(iiext)*w0multi
          value_lpext1(iiext)=value_lpext1(iiext)*w0multi
        enddo
        ilpsta=nstaval(iw0)+1
        ilpend=nstaval(iw0)+nvalue(iw0)
        do iplp=ilpsta,ilpend
          if(logic_dh) then              !lp_head is in dbl_space
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call complete_sd_ar_ext_loop_g(ilw,irw,ilsegdownwei)
          else                                    !lp_head is in act_spa
            ihypos=jphy(iplp)
            ndim =ihy(ihypos)
            iwal0=lpnew_lwei(iplp)
            iwar0=lpnew_rwei(iplp)
            do in=1,ndim
              iwal=iwal0+ihyl(ihypos+in)
              iwar=iwar0+ihy(ihypos+in)
              do iwd=0,iwuplwei-1
                ilw=iwalk_ad(jpadl,ipael,iwal,iwd)
                irw=iwalk_ad(jpad,ipae,iwar,iwd)
            call complete_sd_ar_ext_loop_g(ilw,irw,ilsegdownwei)
              enddo
            enddo
          endif
        enddo
      enddo

      else

      call gsd_ext_sequence(3,iml,imr,lri)

      w0_old=1.0d0
      do iw0=1,mtype
        w0_sdplp=vplpnew_w0(iw0)
        if(logic_dh) w0_sdplp=vplp_w0(iw0)
        w0multi=w0_sdplp/w0_old
        w0_old=w0_sdplp
        do iiext=1,icnt_base
          value_lpext(iiext)=value_lpext(iiext)*w0multi
        enddo
        ilpsta=nstaval(iw0)+1
        ilpend=nstaval(iw0)+nvalue(iw0)
        do iplp=ilpsta,ilpend
          if(logic_dh) then              !lp_head is in dbl_space
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call complete_sd_ar_ext_loop(ilw,irw,ilsegdownwei)
          else                                    !lp_head is in act_spa
            ihypos=jphy(iplp)
            ndim =ihy(ihypos)
            iwal0=lpnew_lwei(iplp)
            iwar0=lpnew_rwei(iplp)
            do in=1,ndim
              iwal=iwal0+ihyl(ihypos+in)
              iwar=iwar0+ihy(ihypos+in)
              do iwd=0,iwuplwei-1
                ilw=iwalk_ad(jpadl,ipael,iwal,iwd)
                irw=iwalk_ad(jpad,ipae,iwar,iwd)
                call complete_sd_ar_ext_loop(ilw,irw,ilsegdownwei)
              enddo
            enddo
          endif
        enddo
      enddo
      endif
      return
      end

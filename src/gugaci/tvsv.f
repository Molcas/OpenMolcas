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
      subroutine tv_drt_ci_new()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      common/lpdisk/idisk_lp,idisk_array(13)
      dsq3=1.732050807568877d0
c      iltype =3
c      irtype =1
c      w0_d48=0.0d0
c      w1_d48=dsq3

      idisk_lp=idisk_array(6)

      do lpblock=1,lpblock_tv
        call read_lp()
        ipael=iml+9
        ipae=1
        call get_jpty(jpadlr,jptyl,jptyr)
        call get_jp(jptyl,jml,jpadl,1)
        call get_jp(jptyr,jmr,jpad,1)
c        jmlr=mul_tab(jml,jmr)
        if(linelp.le.12)   then
          call tv_ext_head_in_act()
        else
          call tv_ext_head_in_dbl()
        endif
      enddo
      return
      end

      subroutine tv_ext_head_in_dbl()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      logic_dh=.true.
      lpok=jpadlr
      jmlr=mul_tab(jml,jmr)
      goto(10,10,10,10,10,106,10,108,10,110,10,10,
     :     113,10,115,10,117,118,10,10,10,10,123,124,10,10),lpok
!=======================================================
!sd1(8) ar- act -br-
108   if(linelp.ne.18) return
      lra=nlg1
      call sdd_ar_act_br_sgt0(17,lra)
      return
!=======================================================
!t1d1(15) ar- act -br-
115   if(linelp.ne.18) return
      lra=nlg1
      call ttdd_ar_act_br_sgt1(17,lra)
      return
!=======================================================
!t1v(18) ar-br- act -c"-
118   if(linelp.ne.14.or.nlg2.ne.2) return
      call ttv_arbr_act_c_stv_sgt1(17)
      return
!=================================================
!d1v(24-1) ar- act -br-
124   if(linelp.ne.18) return
      lra=nlg1
      call d1v_ar_act_bl_ext_ab_sgt0(17,lra)
      return
!=======================================================
!sd(6-1) ar(02)-    act: -b&r-  ext tv: -b^r-a^r
!sd(6-2) c(22)-ar(13)-
!sd(6-3) ar(13)-c'(22)-
!sd(6-4) ar(23)-c'(12)-
106   if(linelp.ne.18) return
      lra=nlg1
      if(jb_sys.gt.0) call sd_ar_act_br_sgt0(17,lra)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jmlr) cycle
        w0sd1=w0_sd(1)
        w0sd2=w0_sd(2)
        w0sd4=w0_sd(4)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1) w0sd1=-w0sd1
        if(ni.eq.1) w0sd2=-w0sd2
        if(ni.eq.1) w0sd4=-w0sd4
!sd(6-1) ar(02)-         act: -br-  ext tv: -br-ar
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
          call ar_br_tv_ext_br_ar(lri,lra)
        endif

!sd(6-2) c(22)-ar(13)-
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
          call ar_br_tv_ext_br_ar(lri,lra)
        enddo
!sd(6-4)   ar(23)-c'(12)-
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
          call ar_br_tv_ext_br_ar(lri,lra)
        enddo
      enddo
      goto 10
!td(13-1) (22)-ar(23)-         act: -b&r-  ext tv: -br-ar
!td(13-2) ar(23)-c'(22)-
113   if(linelp.ne.18) return
      lra=nlg1
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jmlr) cycle
        w0td1=w0_td(1)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1) w0td1=-w0td1
!td(13-1) (22)-ar(23)-        act: -br-  ext tv: -br-ar
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
          call ar_br_tv_ext_br_ar(lri,lra)
        enddo
!td(13-1)   ar(23)-c'(22)-
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
          call ar_br_tv_ext_br_ar(lri,lra)
        enddo
      enddo
      goto 10
!=======================================================================
!dv(23-1)  ar(23)-    act: -b&r-  ext tv: -br-ar
123   if(linelp.ne.18) return
      lra=nlg1
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jml) cycle
        iwdl=jud(lri)
        iwdr=0
        ni=mod(norb_dz-lri,2)
        w0=w0_dv(1)
        if(ni.eq.1) w0=-w0
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
        call ar_br_tv_ext_br_ar(lri,lra)
      enddo
      goto 10
!=======================================================================
!sv(10-1) ar(13)br(23)  act -c"-  tv_ext -br-ar
!sv(10-2) ar(23)br(13)  act -c"-  tv_ext -br-ar
!sv(10-3) dr(03)        act -c"-  tv_ext -br-ar
110   if(linelp.ne.14.or.nlg2.ne.2) return
      if(jb_sys.gt.0) then
        call sv_arbr_act_c_ext_stv_sgt0(17)
      endif
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        do lrj=lri,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
          if(lmij.ne.jmlr) cycle
!-------------------------------------------------------------------
          w0sv2=w0_sv(2)
          w1sv2=w1_sv(2)
          ni=mod(lrj-lri,2)
          if(ni.eq.0) then
            w0sv2=-w0sv2
            w1sv2=-w1sv2
         endif
          iwdl=just(lri,lrj)
          iwdr=0
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo
          if(lri.ne.lrj) then
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sv2
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1sv2
            enddo
            call ar_br_tv_ext_br_ar(lri,lrj)
          else
!sv(10-3) dr(03)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0_sv(3)
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1_sv(3)
            enddo
            call drr_tv_ext_br_ar(lri)
          endif
        enddo
      enddo
      goto 10
!=======================================================================
!tv(17) ar(23)br(23) act -c"- tv_ext -br-ar
117   if(linelp.ne.14.or.nlg2.ne.2) return
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
          w1=w1_tv
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
          call ar_br_tv_ext_br_ar(lri,lrj)
        enddo
      enddo
      goto 10
!=======================================================================
10    return
      end

      subroutine tv_ext_head_in_act()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      logic_dh=.false.

      if(linelp.eq.3) then
!line=3 a&r--b&r<-->b^r-a^r
        lrai=nlg1
        lraj=nlg2
        call ar_br_tv_ext_br_ar(lrai,lraj)
      endif
      if(linelp.eq.8) then
!line=8 d&rr<-->b^r-a^r
        lra=nlg1
        call drr_tv_ext_br_ar(lra)
      endif
       return
      end

      subroutine sv_drt_ci_new()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      common/lpdisk/idisk_lp,idisk_array(13)
      dsq2=1.414213562373095d0
      w0g36a=-1.d0
      w1g36a=0.d0
      w0g13a=-dsq2

c      iltype =4
c      irtype =1
      logic_g36a=.true.
      logic_g13=.false.


      idisk_lp=idisk_array(10)

      do lpblock=1,lpblock_sv
        call read_lp()
        ipael=iml+17
        ipae=1
        if ( iml.eq.1)logic_g13=.true.
        call get_jpty(jpadlr,jptyl,jptyr)
        call get_jp(jptyl,jml,jpadl,1)
        call get_jp(jptyr,jmr,jpad,1)
c        jmlr=mul_tab(jml,jmr)
        if(linelp.le.12)   then
          call sv_ext_head_in_act()
        else
          call sv_ext_head_in_dbl()
        endif
      enddo
      return
      end

      subroutine sv_ext_head_in_dbl()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      logic_dh=.true.
      lpok=jpadlr
      jmlr=mul_tab(jml,jmr)
!      goto(10,10,10,10,10,106,10,10,10,110,10,10,
!     :     113,10,10,10,117,10,10,10,10,10,123,10,10,10),lpok
      goto(10,10,10,10,10,106,10,108,10,110,10,10,
     :     113,10,115,10,117,118,10,10,10,10,123,124,10,10),lpok
!=======================================================
!sd1(8) ar- act -br-
108   if(linelp.ne.18) return
      lra=nlg1
      call sdd_ar_act_br_sgt0(10,lra)
      return
!=======================================================
!t1d1(15) ar- act -br-
115   if(linelp.ne.18) return
      lra=nlg1
      call ttdd_ar_act_br_sgt1(10,lra)
      return
!=======================================================
!t1v(18) ar-br- act -c"-
118   if(linelp.ne.14) return
      call ttv_arbr_act_c_stv_sgt1(10)
      return
!=================================================
!d1v(24-1) ar- act -br-
124   if(linelp.ne.18) return
      lra=nlg1
      call d1v_ar_act_bl_ext_ab_sgt0(10,lra)
      return
!sd(6-1) ar(02)-    act: -b&r-  ext tv: -b^r-a^r
!sd(6-2) c(22)-ar(13)-
!sd(6-3) ar(13)-c'(22)-
!sd(6-4) ar(23)-c'(12)-
106   if(linelp.ne.18) return
      lra=nlg1
      if(jb_sys.gt.0) call sd_ar_act_br_sgt0(10,lra)
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jmlr) cycle
        w0sd1=w0_sd(1)
        w0sd2=w0_sd(2)
        w0sd4=w0_sd(4)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1) w0sd1=-w0sd1
        if(ni.eq.1) w0sd2=-w0sd2
        if(ni.eq.1) w0sd4=-w0sd4
!sd(6-1) ar(02)-         act: -br-  ext tv: -br-ar
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
          call ar_br_sv_ext_br_ar(lri,lra)
        endif

!sd(6-2) c(22)-ar(13)-
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
          call ar_br_sv_ext_br_ar(lri,lra)
        enddo
!sd(6-4)   ar(23)-c'(12)-
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
          call ar_br_sv_ext_br_ar(lri,lra)
        enddo
      enddo
      goto 10
!td(13-1) (22)-ar(23)-         act: -b&r-  ext tv: -br-ar
!td(13-2) ar(23)-c'(22)-
113   if(linelp.ne.18) return
      lra=nlg1
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jmlr) cycle
        w0td1=w0_td(1)
        ni=mod(norb_dz-lri,2)
        if(ni.eq.1) w0td1=-w0td1
!td(13-1) (22)-ar(23)-        act: -br-  ext tv: -br-ar
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
          call ar_br_sv_ext_br_ar(lri,lra)
        enddo
!td(13-1)   ar(23)-c'(22)-
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
          call ar_br_sv_ext_br_ar(lri,lra)
        enddo
      enddo
      goto 10
!=======================================================================
!dv(23-1)  ar(23)-    act: -b&r-  ext tv: -br-ar
123   if(linelp.ne.18) return
      lra=nlg1
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        if(lmi.ne.jml) cycle
        iwdl=jud(lri)
        iwdr=0
        ni=mod(norb_dz-lri,2)
        w0=w0_dv(1)
        if(ni.eq.1) w0=-w0
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
        call ar_br_sv_ext_br_ar(lri,lra)
      enddo
      goto 10
!=======================================================================
!sv(10-1) ar(13)br(23)  act -c"-  tv_ext -br-ar
!sv(10-2) ar(23)br(13)  act -c"-  tv_ext -br-ar
!sv(10-1) dr(03)        act -c"-  tv_ext -br-ar
110   if(linelp.ne.14.or.nlg2.ne.2) return
      if(jb_sys.gt.0) then
        call sv_arbr_act_c_ext_stv_sgt0(10)
      endif
      do lri=norb_frz+1,norb_dz
        lmi=lsm_inn(lri)
        do lrj=lri,norb_dz
          lmj=lsm_inn(lrj)
          lmij=mul_tab(lmi,lmj)
         if(lmij.ne.jmlr) cycle
!-------------------------------------------------------------------
          w0sv2=w0_sv(2)
          w1sv2=w1_sv(2)
          ni=mod(lrj-lri,2)
          if(ni.eq.0) then
            w0sv2=-w0sv2
            w1sv2=-w1sv2
         endif
          iwdl=just(lri,lrj)
          iwdr=0
          do mpl=1,mhlp
            iwal=lpnew_lwei(mpl)
            iwar=lpnew_rwei(mpl)
            lp_lwei(mpl)=iwalk_ad(jpadl,ipael,iwal,iwdl)
            lp_rwei(mpl)=iwalk_ad(jpad,ipae,iwar,iwdr)
          enddo
          if(lri.ne.lrj) then
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0sv2
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1sv2
            enddo
            call ar_br_sv_ext_br_ar(lri,lrj)
          else
!sv(10-1) dr(03)
            do mpl=1,mtype
              vplp_w0(mpl)=vplpnew_w0(mpl)*w0_sv(3)
              vplp_w1(mpl)=vplpnew_w1(mpl)*w1_sv(3)
            enddo
            call drr_sv_ext_br_ar(lri)
          endif
        enddo
      enddo
      goto 10
!=======================================================================
!tv(17) ar(23)br(23) act -c"- tv_ext -br-ar
117   if(linelp.ne.14.or.nlg2.ne.2) return
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
          w1=w1_tv
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
          call ar_br_sv_ext_br_ar(lri,lrj)
        enddo
      enddo
      goto 10
!=======================================================================
10    return
      end

      subroutine sv_ext_head_in_act()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      logic_dh=.false.

      if(linelp.eq.3) then
!line=3 a&r--b&r<-->b^r-a^r
        lrai=nlg1
        lraj=nlg2
        call ar_br_sv_ext_br_ar(lrai,lraj)
      endif
      if(linelp.eq.8) then
!line=8 d&rr<-->b^r-a^r
        lra=nlg1
        call drr_sv_ext_br_ar(lra)
      endif
       return
      end

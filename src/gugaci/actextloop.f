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
!---------------------------------------bbs act-ext
      subroutine ar_bl_dd_ext(lri,lrj,nk)
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
c      write(nf2,*) 'ar_bl_dd_ext'
      logic_g49b=.true.
      iwuplwei=jpad_upwei(jpadl)
      ilsegdownwei=iseg_downwei(ipael)
      irsegdownwei=iseg_downwei(ipae)
      ildownwei_segdd=iseg_downwei(ipael)
      irdownwei_segdd=iseg_downwei(ipae)
      iml0= iml
      imr0= imr
      do iw0=1,mtype
        w0_plp=vplpnew_w0(iw0)
        w1_plp=vplpnew_w1(iw0)
        if(logic_dh) w0_plp=vplp_w0(iw0)
        if(logic_dh) w1_plp=vplp_w1(iw0)

      if(logic_grad) then
        call lp_arbl_ext_dd_calcuvalue_g
     :                      (lri,lrj,iml0,imr0,nlp_value)
        ilpsta=nstaval(iw0)*nk+1
        ilpend=(nstaval(iw0)+nvalue(iw0))*nk
        do iplp=ilpsta,ilpend
          if(logic_dh) then                    !lp_head is in dbl_space
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
           call inn_ext_dd_loop_unpack_g(ilw,irw)
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
                call inn_ext_dd_loop_unpack_g(ilw,irw)
              enddo
            enddo
          endif
        enddo
      else
        call lp_arbl_ext_dd_calcuvalue(lri,lrj,iml0,imr0,nlp_value)
        ilpsta=nstaval(iw0)*nk+1
        ilpend=(nstaval(iw0)+nvalue(iw0))*nk
        do iplp=ilpsta,ilpend
          if(logic_dh) then                    !lp_head is in dbl_space
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
           call inn_ext_dd_loop_unpack(ilw,irw)
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
                call inn_ext_dd_loop_unpack(ilw,irw)
              enddo
            enddo
          endif
        enddo
      endif
      enddo
      return
      end

      subroutine drl_dd_ext(lri)
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
      iwuplwei=jpad_upwei(jpadl)
      ildownwei_segdd=iseg_downwei(ipael)
      irdownwei_segdd=iseg_downwei(ipae)

      do iw0=1,mtype
        w0_plp=vplpnew_w0(iw0)
        w1_plp=vplpnew_w1(iw0)
        if(logic_dh) w0_plp=vplp_w0(iw0)
        if(logic_dh) w1_plp=vplp_w1(iw0)
        ilpsta=nstaval(iw0)+1
        ilpend=nstaval(iw0)+nvalue(iw0)
        do iplp=ilpsta,ilpend
          ilw=lp_lwei(iplp)
          irw=lp_rwei(iplp)
          logic_g49b=.false.
          if(ilw.ne.irw) logic_g49b=.true.

      if(logic_grad) then
!          call lp_drl_ext_dd_calcuvalue_g(lri,iml,nlp_value)
          if(logic_dh) then                       !lp_head is in dbl_spa
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            logic_g49b=.false.
            if(ilw.ne.irw) logic_g49b=.true.
            call lp_drl_ext_dd_calcuvalue_g(lri,iml,nlp_value)
            call inn_ext_dd_loop_unpack_g(ilw,irw)
          else                                    !lp_head is in act_spa
            ihypos=jphy(iplp)
            ndim =ihy(ihypos)
            iwal0=lpnew_lwei(iplp)
            iwar0=lpnew_rwei(iplp)
            logic_g49b=.false.
            if(iwal0.ne.iwar0) logic_g49b=.true.
            call lp_drl_ext_dd_calcuvalue_g(lri,iml,nlp_value)
            do in=1,ndim
              iwal=iwal0+ihyl(ihypos+in)
              iwar=iwar0+ihy(ihypos+in)
              do iwd=0,iwuplwei-1
                ilw=iwalk_ad(jpadl,ipael,iwal,iwd)
                irw=iwalk_ad(jpad,ipae,iwar,iwd)
                call inn_ext_dd_loop_unpack_g(ilw,irw)
            enddo
            enddo
          endif
      else
!          call lp_drl_ext_dd_calcuvalue_wyb(lri,iml,nlp_value)
          if(logic_dh) then                       !lp_head is in dbl_spa
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            logic_g49b=.false.
            if(ilw.ne.irw) logic_g49b=.true.
            call lp_drl_ext_dd_calcuvalue_wyb(lri,iml,nlp_value)
            call inn_ext_dd_loop_unpack(ilw,irw)
          else                                    !lp_head is in act_spa
            ihypos=jphy(iplp)
            ndim =ihy(ihypos)
            iwal0=lpnew_lwei(iplp)
            iwar0=lpnew_rwei(iplp)
            logic_g49b=.false.
            if(iwal0.ne.iwar0) logic_g49b=.true.
            call lp_drl_ext_dd_calcuvalue_wyb(lri,iml,nlp_value)
            do in=1,ndim
              iwal=iwal0+ihyl(ihypos+in)
              iwar=iwar0+ihy(ihypos+in)
              do iwd=0,iwuplwei-1
                ilw=iwalk_ad(jpadl,ipael,iwal,iwd)
                irw=iwalk_ad(jpad,ipae,iwar,iwd)
                call inn_ext_dd_loop_unpack(ilw,irw)
            enddo
            enddo
          endif
       endif
        enddo
      enddo
      return
      end

      subroutine drl_ss_ext(lri)
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      data crl/1.0e-8/

      logic_g1415=.false.
      logic_g2g4b=.false.
      logic_g36b=.false.
      logic_g35b=.false.
      logic_g34b=.false.
      iwuplwei=jpad_upwei(jpadl)
      ilsegdownwei=iseg_downwei(ipael)
      irsegdownwei=iseg_downwei(ipae)
      ildownwei_segdd=iseg_downwei(ipael)
      irdownwei_segdd=iseg_downwei(ipae)

      w0_plp=vplpnew_w0(1)
      if(logic_dh) w0_plp=vplp_w0(1)

      if(logic_grad) then
      call lp_drl_ext_ss_calcuvalue_g(lri,nlp_value)
      w0_old=w0_plp
      do iw0=1,mtype
        if(iw0.eq.1) goto 201
        w0_plp=vplpnew_w0(iw0)
        if(logic_dh) w0_plp=vplp_w0(iw0)
        if(abs(w0_plp).lt.crl) cycle
        w0multi=w0_plp/w0_old
        w0_old=w0_plp
        do iiext=1,nlp_value
         value_lpext(iiext)=value_lpext(iiext)*w0multi
         value_lpext1(iiext)=value_lpext1(iiext)*w0multi
        enddo

201     ilpsta=nstaval(iw0)+1
        ilpend=nstaval(iw0)+nvalue(iw0)
        do iplp=ilpsta,ilpend
          if(logic_dh) then                    !lp_head is in dbl_space
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call inn_ext_ss_drl_loop_unpack_g(ilw,irw)
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
                call inn_ext_ss_drl_loop_unpack_g(ilw,irw)
              enddo
            enddo
          endif
        enddo
      enddo

      else
      call lp_drl_ext_ss_calcuvalue(lri,nlp_value)
      w0_old=w0_plp
      do iw0=1,mtype
        if(iw0.eq.1) goto 202
        w0_plp=vplpnew_w0(iw0)
        if(logic_dh) w0_plp=vplp_w0(iw0)
        if(abs(w0_plp).lt.crl) cycle
        w0multi=w0_plp/w0_old
        w0_old=w0_plp
        do iiext=1,nlp_value
         value_lpext(iiext)=value_lpext(iiext)*w0multi
        enddo

202     ilpsta=nstaval(iw0)+1
        ilpend=nstaval(iw0)+nvalue(iw0)
        do iplp=ilpsta,ilpend
          if(logic_dh) then                    !lp_head is in dbl_space
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call inn_ext_ss_drl_loop_unpack(ilw,irw)
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
                call inn_ext_ss_drl_loop_unpack(ilw,irw)
              enddo
            enddo
          endif
        enddo
      enddo
      endif
      return
      end


      subroutine drl_ss_sum(lri,lrj)
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
                  logic_g1415=.false.
                  logic_g2g4b=.false.
                  logic_g36b=.false.
                  logic_g35b=.false.
                  logic_g34b=.false.
      data crl/1.0e-8/
      if(logic_grad) then
      do lrk=1,norb_dz
         if(lri.eq.lrk)   cycle
         if(lrj.eq.lrk)   cycle
      ilsegdownwei=iseg_downwei(ipael)
      irsegdownwei=iseg_downwei(ipae)
      ildownwei_segdd=iseg_downwei(ipael)
      irdownwei_segdd=iseg_downwei(ipae)
      w0_plp=vplp_w0(1)
      call lp_drl_ext_ss_calcuvalue_g(lrk,nlp_value)
      w0_old=w0_plp
      do iw0=1,mtype
        if(iw0.eq.1) goto 201
        w0_plp=vplp_w0(iw0)
        if(abs(w0_plp).lt.crl) cycle
        w0multi=w0_plp/w0_old
        w0_old=w0_plp
        do iiext=1,nlp_value
         value_lpext(iiext)=value_lpext(iiext)*w0multi
         value_lpext1(iiext)=value_lpext1(iiext)*w0multi
        enddo
201     ilpsta=nstaval(iw0)+1
        ilpend=nstaval(iw0)+nvalue(iw0)
        do iplp=ilpsta,ilpend
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call inn_ext_ss_drl_loop_unpack_g(ilw,irw)
        enddo
      enddo
      enddo

      else
      ilsegdownwei=iseg_downwei(ipael)
      irsegdownwei=iseg_downwei(ipae)
      ildownwei_segdd=iseg_downwei(ipael)
      irdownwei_segdd=iseg_downwei(ipae)
      w0_plp=vplp_w0(1)
      call lp_drl_sum_ss_calcuvalue(lri,lrj,nlp_value)
      w0_old=w0_plp
      do iw0=1,mtype
        if(iw0.eq.1) goto 202
        w0_plp=vplp_w0(iw0)
        if(abs(w0_plp).lt.crl) cycle
        w0multi=w0_plp/w0_old
        w0_old=w0_plp
        do iiext=1,nlp_value
         value_lpext(iiext)=value_lpext(iiext)*w0multi
        enddo
202     ilpsta=nstaval(iw0)+1
        ilpend=nstaval(iw0)+nvalue(iw0)
        do iplp=ilpsta,ilpend
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call inn_ext_ss_drl_loop_unpack(ilw,irw)
        enddo
      enddo
      endif

      return
      end

      subroutine drl_st_ext(lri)
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      data crl/1.0e-8/

      iwuplwei=jpad_upwei(jpadl)
      ilsegdownwei=iseg_downwei(ipael)
      irsegdownwei=iseg_downwei(ipae)
      ildownwei_segdd=iseg_downwei(ipael)
      irdownwei_segdd=iseg_downwei(ipae)
c      logic_g1415=.false.
c      logic_g34b=.false.
c      logic_g35b=.false.
c      logic_g36b=.false.

      w1_plp=vplpnew_w1(1)
      if(logic_dh) w1_plp=vplp_w1(1)

      if(logic_grad) then
      call lp_drl_ext_st_calcuvalue_g(lri,nlp_value)
      w1_old=w1_plp
      do iw0=1,mtype
        if(iw0.eq.1) goto 201
        w1_plp=vplpnew_w1(iw0)
        if(logic_dh) w1_plp=vplp_w1(iw0)
        if(abs(w1_plp).lt.crl) cycle
        w1multi=w1_plp/w1_old
        w1_old=w1_plp

        do iiext=1,nlp_value
         value_lpext(iiext)=value_lpext(iiext)*w1multi
         value_lpext1(iiext)=value_lpext1(iiext)*w1multi
        enddo
201     ilpsta=nstaval(iw0)+1
        ilpend=nstaval(iw0)+nvalue(iw0)
        do iplp=ilpsta,ilpend

         if(logic_dh) then                    !lp_head is in dbl_space
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call inn_ext_st_drl_loop_unpack_g(ilw,irw)
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
                call inn_ext_st_drl_loop_unpack_g(ilw,irw)
              enddo
            enddo
          endif
        enddo
      enddo

      else
      call lp_drl_ext_st_calcuvalue(lri,nlp_value)
      w1_old=w1_plp
      do iw0=1,mtype
        if(iw0.eq.1) goto 202
        w1_plp=vplpnew_w1(iw0)
        if(logic_dh) w1_plp=vplp_w1(iw0)
        if(abs(w1_plp).lt.crl) cycle
        w1multi=w1_plp/w1_old
        w1_old=w1_plp
c        call calcu_drl2_value_wyb(1,lri,lrj)
c        call lp_drl_ext_st_calcuvalue(lri,nlp_value)
        do iiext=1,nlp_value
         value_lpext(iiext)=value_lpext(iiext)*w1multi
        enddo
202     ilpsta=nstaval(iw0)+1
        ilpend=nstaval(iw0)+nvalue(iw0)
        do iplp=ilpsta,ilpend
c      write(6,*) '      calcuvalue_w', intentry

         if(logic_dh) then                    !lp_head is in dbl_space
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call inn_ext_st_drl_loop_unpack(ilw,irw)
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
                call inn_ext_st_drl_loop_unpack(ilw,irw)
              enddo
            enddo
          endif
        enddo
      enddo

      endif
      return
      end

      subroutine drl_tt_ext(lri)
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

c                  logic_g1415=.false.
c                  logic_g2g4b=.false.
c                  logic_g36b=.false.
c                  logic_g35b=.false.
c                  logic_g34b=.false.
      iwuplwei=jpad_upwei(jpadl)
      ilsegdownwei=iseg_downwei(ipael)
      irsegdownwei=iseg_downwei(ipae)
      ildownwei_segdd=iseg_downwei(ipael)
      irdownwei_segdd=iseg_downwei(ipae)

      do iw0=1,mtype
        w0_plp=vplpnew_w0(iw0)
        w1_plp=vplpnew_w1(iw0)
        if(logic_dh) w0_plp=vplp_w0(iw0)
        if(logic_dh) w1_plp=vplp_w1(iw0)

       if(logic_grad) then
        call lp_drl_ext_tt_calcuvalue_g(lri,n1415,nlp_value)
        ilpsta=nstaval(iw0)+1
        ilpend=nstaval(iw0)+nvalue(iw0)
        do iplp=ilpsta,ilpend

          if(logic_dh) then                    !lp_head is in dbl_space
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call inn_ext_tt_drl_loop_unpack_g(ilw,irw,n1415)
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
               call inn_ext_tt_drl_loop_unpack_g(ilw,irw,n1415)
              enddo
            enddo
          endif
        enddo
       else
        call lp_drl_ext_tt_calcuvalue(lri,n1415,nlp_value)
        ilpsta=nstaval(iw0)+1
        ilpend=nstaval(iw0)+nvalue(iw0)
        do iplp=ilpsta,ilpend

          if(logic_dh) then                    !lp_head is in dbl_space
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call inn_ext_tt_drl_loop_unpack(ilw,irw,n1415)
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
               call inn_ext_tt_drl_loop_unpack(ilw,irw,n1415)
              enddo
            enddo
          endif
        enddo
       endif
      enddo
      return
      end

      subroutine drl_tt_sum(lri,lrj)
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      data crl/1.0e-8/
c                  logic_g1415=.false.
c                  logic_g2g4b=.false.
c                  logic_g36b=.false.
c                  logic_g35b=.false.
c                  logic_g34b=.false.
c      write(nf2,*) 'drl_tt_sum'
      if(logic_grad) then
      do lrk=1,norb_dz
         if(lri.eq.lrk)   cycle
         if(lrj.eq.lrk)   cycle
      ilsegdownwei=iseg_downwei(ipael)
      irsegdownwei=iseg_downwei(ipae)
      ildownwei_segdd=iseg_downwei(ipael)
      irdownwei_segdd=iseg_downwei(ipae)
      w0_plp=vplp_w0(1)
      call lp_drl_sum_tt_calcuvalue_g(lrk,n1415,nlp_value)
      w0_old=w0_plp
      do iw0=1,mtype
        if(iw0.eq.1) goto 201
        w0_plp=vplp_w0(iw0)
        if(abs(w0_plp).lt.crl) cycle
        w0multi=w0_plp/w0_old
        w0_old=w0_plp
        do iiext=1,nlp_value
         value_lpext(iiext)=value_lpext(iiext)*w0multi
         value_lpext1(iiext)=value_lpext1(iiext)*w0multi
        enddo

201     ilpsta=nstaval(iw0)+1
        ilpend=nstaval(iw0)+nvalue(iw0)
        do iplp=ilpsta,ilpend
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call inn_ext_tt_drl_loop_unpack_g(ilw,irw,n1415)
        enddo
      enddo
      enddo

      else
      ilsegdownwei=iseg_downwei(ipael)
      irsegdownwei=iseg_downwei(ipae)
      ildownwei_segdd=iseg_downwei(ipael)
      irdownwei_segdd=iseg_downwei(ipae)
      w0_plp=vplp_w0(1)
      call lp_drl_sum_tt_calcuvalue(lri,lrj,n1415,nlp_value)
      w0_old=w0_plp
      do iw0=1,mtype
        if(iw0.eq.1) goto 202
        w0_plp=vplp_w0(iw0)
        if(abs(w0_plp).lt.crl) cycle
        w0multi=w0_plp/w0_old
        w0_old=w0_plp
        do iiext=1,nlp_value
         value_lpext(iiext)=value_lpext(iiext)*w0multi
        enddo
c      call lp_drl_sum_tt_calcuvalue(lri,lrj,n1415,nlp_value)
202     ilpsta=nstaval(iw0)+1
        ilpend=nstaval(iw0)+nvalue(iw0)
        do iplp=ilpsta,ilpend
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call inn_ext_tt_drl_loop_unpack(ilw,irw,n1415)
        enddo
      enddo
      endif
      return
      end

      subroutine drl_ts_ext(lri)
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      data crl/1.0e-8/

c                  logic_g1415=.false.
c                  logic_g2g4b=.false.
c                  logic_g36b=.false.
c                  logic_g35b=.false.
c                  logic_g34b=.false.
      ildownwei_segdd=iseg_downwei(ipael)
      irdownwei_segdd=iseg_downwei(ipae)
      iwuplwei=jpad_upwei(jpadl)
      ilsegdownwei=iseg_downwei(ipael)
      irsegdownwei=iseg_downwei(ipae)
      w1_plp=vplpnew_w1(1)
      if(logic_dh) w1_plp=vplp_w1(1)

      if(logic_grad) then
      call lp_drl_ext_ts_calcuvalue_g(lri,nlp_value)
      w1_old=w1_plp
      do iw0=1,mtype
        if(iw0.eq.1) goto 201
        w1_plp=vplpnew_w1(iw0)
        if(logic_dh) w1_plp=vplp_w1(iw0)
        if(abs(w1_plp).lt.crl) cycle
        w1multi=w1_plp/w1_old
        w1_old=w1_plp
        do iiext=1,nlp_value
         value_lpext(iiext)=value_lpext(iiext)*w1multi
         value_lpext1(iiext)=value_lpext1(iiext)*w1multi
        enddo
201     ilpsta=nstaval(iw0)+1
        ilpend=nstaval(iw0)+nvalue(iw0)
        do iplp=ilpsta,ilpend

          if(logic_dh) then                    !lp_head is in dbl_space
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call inn_ext_ts_drl_loop_unpack_g(ilw,irw)
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
                call inn_ext_ts_drl_loop_unpack_g(ilw,irw)
              enddo
            enddo
          endif
        enddo
      enddo

      else
      call lp_drl_ext_ts_calcuvalue(lri,nlp_value)
      w1_old=w1_plp
      do iw0=1,mtype
        if(iw0.eq.1) goto 202
        w1_plp=vplpnew_w1(iw0)
        if(logic_dh) w1_plp=vplp_w1(iw0)
        if(abs(w1_plp).lt.crl) cycle
        w1multi=w1_plp/w1_old
        w1_old=w1_plp
        do iiext=1,nlp_value
         value_lpext(iiext)=value_lpext(iiext)*w1multi
        enddo
202     ilpsta=nstaval(iw0)+1
        ilpend=nstaval(iw0)+nvalue(iw0)
        do iplp=ilpsta,ilpend

          if(logic_dh) then                    !lp_head is in dbl_space
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call inn_ext_ts_drl_loop_unpack(ilw,irw)
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
                call inn_ext_ts_drl_loop_unpack(ilw,irw)
              enddo
            enddo
          endif
        enddo
      enddo

      endif
      return
      end

      subroutine ar_br_tv_ext_br_ar(lri,lrj)
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"

      ilsegdownwei=iseg_downwei(ipael)
      irsegdownwei=iseg_downwei(ipae)
      iwuplwei=jpad_upwei(jpadl)
      w0g36a=0.d0
      w1g36a=-1.d0
      do iw0=1,mtype
        w0_plp=vplpnew_w0(iw0)
        w1_plp=vplpnew_w1(iw0)
        if(logic_dh) w0_plp=vplp_w0(iw0)
        if(logic_dh) w1_plp=vplp_w1(iw0)

      if(logic_grad) then
       call lp_arbr_ext_svtv_calcuvalue_g(lri,lrj,nlp_value)
        ilpsta=nstaval(iw0)+1
        ilpend=nstaval(iw0)+nvalue(iw0)
        do iplp=ilpsta,ilpend
          if(logic_dh) then
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call inn_ext_sv_loop_unpack_g(ilw,irw)
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
                call inn_ext_sv_loop_unpack_g(ilw,irw)
              enddo
            enddo
          endif
        enddo
      else
        call lp_arbr_ext_svtv_calcuvalue_wyb(lri,lrj,nlp_value)
        ilpsta=nstaval(iw0)+1
        ilpend=nstaval(iw0)+nvalue(iw0)
        do iplp=ilpsta,ilpend
          if(logic_dh) then
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call inn_ext_sv_loop_unpack(ilw,irw)
          else
            if(log_prod.eq.3) then
              lphead=lpnew_head(iplp)
              ihyposl=jphyl(lphead)
              ihyposr=jphy(lphead)
              ndim =ihyl(ihyposl)
            else
              ihyposl=jphy(iplp)
              ihyposr=ihyposl
              ndim =ihy(ihyposl)
            endif
            iwal0=lpnew_lwei(iplp)
            iwar0=lpnew_rwei(iplp)
            do in=1,ndim
              iwal=iwal0+ihyl(ihyposl+in)
              iwar=iwar0+ihy(ihyposr+in)
              do iwd=0,iwuplwei-1
                ilw=iwalk_ad(jpadl,ipael,iwal,iwd)
                irw=iwalk_ad(jpad,ipae,iwar,iwd)
                call inn_ext_sv_loop_unpack(ilw,irw)
              enddo
            enddo
          endif
        enddo
      endif
      enddo
      return
      end

      subroutine drr_sv_ext_br_ar(lri)
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"

      ndorb=norb_dz-lri
      ilsegdownwei=iseg_downwei(ipael)
      irsegdownwei=iseg_downwei(ipae)
      iwuplwei=jpad_upwei(jpadl)
      do iw0=1,mtype
        w0_plp=vplpnew_w0(iw0)
        w1_plp=vplpnew_w1(iw0)
        if(ndorb.ge.0) w0_plp=vplp_w0(iw0)
        if(ndorb.ge.0) w1_plp=vplp_w1(iw0)

      if(logic_grad) then
        call lp_drr_ext_svtv_calcuvalue_g(lri,nlp_value)
        ilpsta=nstaval(iw0)+1
        ilpend=nstaval(iw0)+nvalue(iw0)
        do iplp=ilpsta,ilpend
          if(logic_dh) then
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
           call inn_ext_sv_loop_unpack_g(ilw,irw)
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
                call inn_ext_sv_loop_unpack_g(ilw,irw)
              enddo
            enddo
          endif
        enddo
      else
        call lp_drr_ext_svtv_calcuvalue_wyb(lri,nlp_value)
        ilpsta=nstaval(iw0)+1
        ilpend=nstaval(iw0)+nvalue(iw0)
        do iplp=ilpsta,ilpend
          if(logic_dh) then
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
           call inn_ext_sv_loop_unpack(ilw,irw)
          else
            if(log_prod.eq.3) then
              lphead=lpnew_head(iplp)
              ihyposl=jphyl(lphead)
              ihyposr=jphy(lphead)
              ndim =ihyl(ihyposl)
            else
              ihyposl=jphy(iplp)
              ihyposr=ihyposl
              ndim =ihy(ihyposl)
            endif

            iwal0=lpnew_lwei(iplp)
            iwar0=lpnew_rwei(iplp)
            do in=1,ndim
              iwal=iwal0+ihyl(ihyposl+in)
              iwar=iwar0+ihy(ihyposr+in)
              do iwd=0,iwuplwei-1
                ilw=iwalk_ad(jpadl,ipael,iwal,iwd)
                irw=iwalk_ad(jpad,ipae,iwar,iwd)
                call inn_ext_sv_loop_unpack(ilw,irw)
              enddo
            enddo
          endif
        enddo
      endif
      enddo
      return
      end
      subroutine drr_tv_ext_br_ar(lri)
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"

      ndorb=norb_dz-lri
      ilsegdownwei=iseg_downwei(ipael)
      irsegdownwei=iseg_downwei(ipae)
      iwuplwei=jpad_upwei(jpadl)
      do iw0=1,mtype
        w0_plp=vplpnew_w0(iw0)
        w1_plp=vplpnew_w1(iw0)
        if(ndorb.ge.0) w0_plp=vplp_w0(iw0)
        if(ndorb.ge.0) w1_plp=vplp_w1(iw0)

      if(logic_grad) then
        call lp_drr_ext_svtv_calcuvalue_g(lri,nlp_value)
        ilpsta=nstaval(iw0)+1
        ilpend=nstaval(iw0)+nvalue(iw0)
        do iplp=ilpsta,ilpend
          if(logic_dh) then
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
           call inn_ext_sv_loop_unpack_g(ilw,irw)
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
                call inn_ext_sv_loop_unpack_g(ilw,irw)
              enddo
            enddo
          endif
        enddo
      else
        call lp_drr_ext_svtv_calcuvalue_wyb(lri,nlp_value)
        ilpsta=nstaval(iw0)+1
        ilpend=nstaval(iw0)+nvalue(iw0)
        do iplp=ilpsta,ilpend
          if(logic_dh) then
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
           call inn_ext_sv_loop_unpack(ilw,irw)
          else
            if(log_prod.eq.3) then
              lphead=lpnew_head(iplp)
              ihyposl=jphyl(lphead)
              ihyposr=jphy(lphead)
              ndim =ihyl(ihyposl)
            else
              ihyposl=jphy(iplp)
              ihyposr=ihyposl
              ndim =ihy(ihyposl)
            endif

            iwal0=lpnew_lwei(iplp)
            iwar0=lpnew_rwei(iplp)
            do in=1,ndim
              iwal=iwal0+ihyl(ihyposl+in)
              iwar=iwar0+ihy(ihyposr+in)
              do iwd=0,iwuplwei-1
                ilw=iwalk_ad(jpadl,ipael,iwal,iwd)
                irw=iwalk_ad(jpad,ipae,iwar,iwd)
                call inn_ext_sv_loop_unpack(ilw,irw)
              enddo
            enddo
          endif
        enddo
      endif
      enddo
      return
      end

      subroutine ar_dv_ext_ar(idtu,isma,lri,lrj)
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
      dimension lpcoe(norb_dz+1:norb_inn)
      iwuplwei=jpad_upwei(jpadl)
      ilsegdownwei=iseg_downwei(ipael)
      irsegdownwei=iseg_downwei(ipae)
      do iw0=1,mtype
        w0_sdplp=vplpnew_w0(iw0)
        if(logic_dh) w0_sdplp=vplp_w0(iw0)
        ilpsta=nstaval(iw0)+1
        ilpend=nstaval(iw0)+nvalue(iw0)
        do iplp=ilpsta,ilpend
c                        ipcoe=lp_arpos(iplp)
          do norb=norb_dz+1,norb_inn
            lpcoe(norb)=lpnew_coe(norb,iplp)       !01.12.25
          enddo

      if(logic_grad) then
         call lp_ar_coe_calcuvalue_g
     *            (idtu,isma,lri,lrj,nlp_value,lpcoe,nvalue1)
          if(logic_dh) then              !lp_head is in dbl_space
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call gdv_sequence_extspace1_g(ilw,irw,nvalue1)
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
                call gdv_sequence_extspace1_g(ilw,irw,nvalue1)
              enddo
            enddo
         endif
       else
         call lp_ar_coe_calcuvalue_wyb
     *            (idtu,isma,lri,lrj,nlp_value,lpcoe)
         if(logic_dh) then              !lp_head is in dbl_space
           ilw=lp_lwei(iplp)
           irw=lp_rwei(iplp)
           call gdv_sequence_extspace(ilw,irw)
         else           !lp_head is in act_space
           if(log_prod.eq.3) then
             lphead=lpnew_head(iplp)
             ihyposl=jphyl(lphead)
             ihyposr=jphy(lphead)
             ndim =ihyl(ihyposl)
           else
             ihyposl=jphy(iplp)
             ihyposr=ihyposl
             ndim =ihy(ihyposl)
           endif

           iwal0=lpnew_lwei(iplp)
           iwar0=lpnew_rwei(iplp)
           do in=1,ndim
             iwal=iwal0+ihyl(ihyposl+in)
             iwar=iwar0+ihy(ihyposr+in)
             do iwd=0,iwuplwei-1
               ilw=iwalk_ad(jpadl,ipael,iwal,iwd)
               irw=iwalk_ad(jpad,ipae,iwar,iwd)
               call gdv_sequence_extspace(ilw,irw)
             enddo
           enddo
         endif
      endif
        enddo
      enddo
      return
      end

      subroutine ar_sd_ext_ar(idtu,lri,lrj,isma)
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
      dimension lpcoe(norb_dz+1:norb_inn)

      iwuplwei=jpad_upwei(jpadl)
      ilsegdownwei=iseg_downwei(ipael)
      irsegdownwei=iseg_downwei(ipae)
      do iw0=1,mtype
        w0_sdplp=vplpnew_w0(iw0)
        if(logic_dh) w0_sdplp=vplp_w0(iw0)
        ilpsta=nstaval(iw0)+1
        ilpend=nstaval(iw0)+nvalue(iw0)
        do iplp=ilpsta,ilpend

          do norb=norb_dz+1,norb_inn
            lpcoe(norb)=lpnew_coe(norb,iplp)
          enddo

      if(logic_grad) then
          call lp_ar_coe_calcuvalue_g
     *            (idtu,isma,lri,lrj,nlp_value,lpcoe,nvalue1)
          if(logic_dh) then              !lp_head is in dbl_space
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call gsd_sequence_extspace1_g(ilw,irw,nvalue1)
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
                call gsd_sequence_extspace1_g(ilw,irw,nvalue1)
              enddo
            enddo
         endif
      else
          call lp_ar_coe_calcuvalue_wyb
     *            (idtu,isma,lri,lrj,nlp_value,lpcoe)
          if(logic_dh) then              !lp_head is in dbl_space
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call gsd_sequence_extspace(ilw,irw)
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
                call gsd_sequence_extspace(ilw,irw)
              enddo
            enddo
         endif
      endif
        enddo
      enddo
      return
      end

      subroutine ar_td_ext_ar(idtu,lri,lrj,isma)
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
      dimension lpcoe(norb_dz+1:norb_inn)

      iwuplwei=jpad_upwei(jpadl)
      ilsegdownwei=iseg_downwei(ipael)
      irsegdownwei=iseg_downwei(ipae)
      do iw0=1,mtype
        w0_sdplp=vplpnew_w0(iw0)
        if(logic_dh) w0_sdplp=vplp_w0(iw0)
        ilpsta=nstaval(iw0)+1
        ilpend=nstaval(iw0)+nvalue(iw0)
        do iplp=ilpsta,ilpend
          do norb=norb_dz+1,norb_inn
            lpcoe(norb)=lpnew_coe(norb,iplp)
          enddo

      if(logic_grad) then
         call lp_ar_coe_calcuvalue_g
     *            (idtu,isma,lri,lrj,nlp_value,lpcoe,nvalue1)
          if(logic_dh) then              !lp_head is in dbl_space
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)

            call gtd_sequence_extspace1_g(ilw,irw,nvalue1)
          else                           !lp_head is in act_space
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
                call gtd_sequence_extspace1_g(ilw,irw,nvalue1)
              enddo
            enddo
          endif
      else
         call lp_ar_coe_calcuvalue_wyb
     *            (idtu,isma,lri,lrj,nlp_value,lpcoe)
          if(logic_dh) then              !lp_head is in dbl_space
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call gtd_sequence_extspace(ilw,irw)
          else                           !lp_head is in act_space
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
                call gtd_sequence_extspace(ilw,irw)
              enddo
            enddo
          endif
      endif
        enddo
      enddo
      return
      end

      subroutine ar_br_sv_ext_br_ar(lri,lrj)
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"

      ilsegdownwei=iseg_downwei(ipael)
      irsegdownwei=iseg_downwei(ipae)
      iwuplwei=jpad_upwei(jpadl)
      do iw0=1,mtype
        w0_plp=vplpnew_w0(iw0)
        w1_plp=vplpnew_w1(iw0)
        if(logic_dh) w0_plp=vplp_w0(iw0)
        if(logic_dh) w1_plp=vplp_w1(iw0)

      if(logic_grad) then
        call lp_arbr_ext_svtv_calcuvalue_g(lri,lrj,nlp_value)
        ilpsta=nstaval(iw0)+1
        ilpend=nstaval(iw0)+nvalue(iw0)
        do iplp=ilpsta,ilpend
          if(logic_dh) then
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call inn_ext_sv_loop_unpack_g(ilw,irw)
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
                call inn_ext_sv_loop_unpack_g(ilw,irw)
              enddo
            enddo
          endif
        enddo
      else
        call lp_arbr_ext_svtv_calcuvalue_wyb(lri,lrj,nlp_value)
        ilpsta=nstaval(iw0)+1
        ilpend=nstaval(iw0)+nvalue(iw0)
        do iplp=ilpsta,ilpend
          if(logic_dh) then
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call inn_ext_sv_loop_unpack(ilw,irw)
          else
            if(log_prod.eq.3) then
              lphead=lpnew_head(iplp)
              ihyposl=jphyl(lphead)
              ihyposr=jphy(lphead)
              ndim =ihyl(ihyposl)
            else
              ihyposl=jphy(iplp)
              ihyposr=ihyposl
              ndim =ihy(ihyposl)
            endif
            iwal0=lpnew_lwei(iplp)
            iwar0=lpnew_rwei(iplp)
            do in=1,ndim
              iwal=iwal0+ihyl(ihyposl+in)
              iwar=iwar0+ihy(ihyposr+in)
              do iwd=0,iwuplwei-1
                ilw=iwalk_ad(jpadl,ipael,iwal,iwd)
                irw=iwalk_ad(jpad,ipae,iwar,iwd)
                call inn_ext_sv_loop_unpack(ilw,irw)
              enddo
            enddo
          endif
        enddo
      endif
      enddo
      return
      end

      subroutine logicg_dd(ilnodesm,irnodesm)
#include "drt_h.fh"
#include "lpextmode_h.fh"
      logic_g50 =.false.
      logic_g49a=.false.
      logic_g49b=.false.
      if (ilnodesm.lt.irnodesm) then
            logic_g49a=.true.
      elseif ( ilnodesm.eq.irnodesm) then
            logic_g49a=.true.
            logic_g49b=.true.
            logic_g50 =.true.
      else
            logic_g49b=.true.
      endif
      end

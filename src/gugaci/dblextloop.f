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
      subroutine drl_bl_ext_ar_new(lin,lrk,lri)
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"

      isma=lsm_inn(lri)
      iwuplwei=jpad_upwei(jpadl)
      ilsegdownwei=iseg_downwei(ipael)
      irsegdownwei=iseg_downwei(ipae)
      do iw0=1,mtype
        w0_sdplp=vplpnew_w0(iw0)
        w1_sdplp=vplpnew_w1(iw0)
        if(logic_dh) w0_sdplp=vplp_w0(iw0)
        if(logic_dh) w1_sdplp=vplp_w1(iw0)
        if(logic_grad)   then
          call lp9_drlbl_ext_calcuvalue_g(lri,lrk,isma)
          ilpsta=nstaval(iw0)+1
          ilpend=nstaval(iw0)+nvalue(iw0)
          do iplp=ilpsta,ilpend
            if(logic_dh) then
              ilw=lp_lwei(iplp)
              irw=lp_rwei(iplp)
              if(lin.eq.6)  call gsd_sequence_extspace_g(ilw,irw)
              if(lin.eq.13) call gtd_sequence_extspace_g(ilw,irw)
              if(lin.eq.23) call gdv_sequence_extspace_g(ilw,irw)
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
              if(lin.eq.6)  call gsd_sequence_extspace_g(ilw,irw)
              if(lin.eq.13) call gtd_sequence_extspace_g(ilw,irw)
              if(lin.eq.23) call gdv_sequence_extspace_g(ilw,irw)
                enddo
              enddo
            endif
          enddo
        else
          call lp9_drlbl_ext_calcuvalue_wyb(lri,lrk,isma)
          ilpsta=nstaval(iw0)+1
          ilpend=nstaval(iw0)+nvalue(iw0)
          do iplp=ilpsta,ilpend
            if(logic_dh) then
              ilw=lp_lwei(iplp)
              irw=lp_rwei(iplp)
              if(lin.eq.6)  call gsd_sequence_extspace(ilw,irw)
              if(lin.eq.13) call gtd_sequence_extspace(ilw,irw)
              if(lin.eq.23) call gdv_sequence_extspace(ilw,irw)
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
                  if(lin.eq.6)  call gsd_sequence_extspace(ilw,irw)
                  if(lin.eq.13) call gtd_sequence_extspace(ilw,irw)
                  if(lin.eq.23) call gdv_sequence_extspace(ilw,irw)
                enddo
              enddo
            endif
          enddo

        endif
      enddo
      return
      end

      subroutine drl_br_ext_al_new(lin,lrk,lri )   !to be reviced
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
      data crl/1.0e-8/

      iwuplwei=jpad_upwei(jpad)
      ilsegdownwei=iseg_downwei(ipae)
      irsegdownwei=iseg_downwei(ipael)
      isma=lsm_inn(lri)

      do iw0=1,mtype
        w0_sdplp=vplpnew_w0(iw0)-vplpnew_w1(iw0)
        if(logic_dh) w0_sdplp=vplp_w0(iw0)-vplp_w1(iw0)
        if(abs(w0_sdplp).lt.crl) cycle
        if(logic_grad) then
          call lp678_ext_calcuvalue_g(lri,lrk,isma,nlp_value)
          ilpsta=nstaval(iw0)+1
          ilpend=nstaval(iw0)+nvalue(iw0)
          do iplp=ilpsta,ilpend
            if(logic_dh) then
              ilw=lp_lwei(iplp)
              irw=lp_rwei(iplp)
              if(lin.eq.7)  call gsd_sequence_extspace_g(irw,ilw)
              if(lin.eq.14) call gtd_sequence_extspace_g(irw,ilw)
              if(lin.eq.26) call gdv_sequence_extspace_g(irw,ilw)
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
              if(lin.eq.7)  call gsd_sequence_extspace_g(irw,ilw)
              if(lin.eq.14) call gtd_sequence_extspace_g(irw,ilw)
              if(lin.eq.26) call gdv_sequence_extspace_g(irw,ilw)
                enddo
              enddo
            endif
          enddo

        else
          call lp678_ext_wyb_calcuvalue(lri,lrk,isma,nlp_value)
          ilpsta=nstaval(iw0)+1
          ilpend=nstaval(iw0)+nvalue(iw0)
          do iplp=ilpsta,ilpend
            if(logic_dh) then
              ilw=lp_lwei(iplp)
              irw=lp_rwei(iplp)
              if(lin.eq.7)  call gsd_sequence_extspace(irw,ilw)
              if(lin.eq.14) call gtd_sequence_extspace(irw,ilw)
              if(lin.eq.26) call gdv_sequence_extspace(irw,ilw)
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
                  if(lin.eq.7)  call gsd_sequence_extspace(irw,ilw)
                  if(lin.eq.14) call gtd_sequence_extspace(irw,ilw)
                  if(lin.eq.26) call gdv_sequence_extspace(irw,ilw)
                enddo
              enddo
            endif
          enddo
        endif
      enddo
      return
      end

      subroutine ar_bl_br_ext_al_new(lin,intentry,isma,nk)
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
c      if(lin.eq.14.and.intentry.eq.15195.and.isma
c     :    .eq.1.and.jpad.eq.20.and.jpad.eq.jpadl) then
c      write(6,*) "bbs_tmp"
c      endif

      iwuplwei=jpad_upwei(jpad)
      ilsegdownwei=iseg_downwei(ipae)
      irsegdownwei=iseg_downwei(ipael)
      do iw0=1,mtype
        w0_sdplp=vplpnew_w0(iw0)
        w1_sdplp=vplpnew_w1(iw0)
        if(logic_dh) w0_sdplp=vplp_w0(iw0)
        if(logic_dh) w1_sdplp=vplp_w1(iw0)

        if(logic_grad) then
          call lp11_arblbr_ext_calcuvalue_g(intentry,isma,nlp_value)
          ilpsta=nstaval(iw0)*nk+1
          ilpend=(nstaval(iw0)+nvalue(iw0))*nk
          do iplp=ilpsta,ilpend
            if(logic_dh) then                     !lp_head is in dbl_spa
              ilw=lp_lwei(iplp)
              irw=lp_rwei(iplp)
              if(lin.eq.7)  call gsd_sequence_extspace_g(irw,ilw)
              if(lin.eq.14) call gtd_sequence_extspace_g(irw,ilw)
              if(lin.eq.26) call gdv_sequence_extspace_g(irw,ilw)
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
              if(lin.eq.7)  call gsd_sequence_extspace_g(irw,ilw)
              if(lin.eq.14) call gtd_sequence_extspace_g(irw,ilw)
              if(lin.eq.26) call gdv_sequence_extspace_g(irw,ilw)
                enddo
              enddo
            endif
          enddo
        else
          call lp11_arblbr_ext_calcuvalue(intentry,isma,nlp_value)
          ilpsta=nstaval(iw0)*nk+1
          ilpend=(nstaval(iw0)+nvalue(iw0))*nk
          do iplp=ilpsta,ilpend
            if(logic_dh) then                     !lp_head is in dbl_spa
              ilw=lp_lwei(iplp)
              irw=lp_rwei(iplp)
              if(lin.eq.7)  call gsd_sequence_extspace(irw,ilw)
              if(lin.eq.14) call gtd_sequence_extspace(irw,ilw)
              if(lin.eq.26) call gdv_sequence_extspace(irw,ilw)
            else                                    !lp_head is in act_s
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
                  if(lin.eq.7)  call gsd_sequence_extspace(irw,ilw)
                  if(lin.eq.14) call gtd_sequence_extspace(irw,ilw)
                  if(lin.eq.26) call gdv_sequence_extspace(irw,ilw)
                enddo
              enddo
            endif
          enddo
        endif
      enddo
      return
      end


      subroutine ar_drl_ext_al_new(lin,lri,lrk )
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
      data crl/1.0e-8/

      isma=lsm_inn(lri)
      iwuplwei=jpad_upwei(jpad)
      ilsegdownwei=iseg_downwei(ipae)
      irsegdownwei=iseg_downwei(ipael)
      w0_sdplp=vplpnew_w0(1)
      if(logic_dh) w0_sdplp=vplp_w0(1)

      if(logic_grad) then
      call lp678_ext_calcuvalue_g(lri,lrk,isma,nlp_value)
      w0_sdold=w0_sdplp
      do iw0=1,mtype
        if(iw0.eq.1) goto 201
        w0_sdplp=vplpnew_w0(iw0)
        if(logic_dh) w0_sdplp=vplp_w0(iw0)
        if(abs(w0_sdplp).lt.crl) cycle
        w0multi=w0_sdplp/w0_sdold
        w0_sdold=w0_sdplp
        do iiext=1,nlp_value
         value_lpext(iiext)=value_lpext(iiext)*w0multi
         value_lpext1(iiext)=value_lpext1(iiext)*w0multi
        enddo
201     ilpsta=nstaval(iw0)+1
        ilpend=nstaval(iw0)+nvalue(iw0)
        do iplp=ilpsta,ilpend
          if(logic_dh) then
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            if(lin.eq.7)  call gsd_sequence_extspace_g(irw,ilw)
            if(lin.eq.14) call gtd_sequence_extspace_g(irw,ilw)
            if(lin.eq.26) call gdv_sequence_extspace_g(irw,ilw)
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
            if(lin.eq.7)  call gsd_sequence_extspace_g(irw,ilw)
            if(lin.eq.14) call gtd_sequence_extspace_g(irw,ilw)
            if(lin.eq.26) call gdv_sequence_extspace_g(irw,ilw)
              enddo
            enddo
          endif
        enddo
      enddo

      else
      call lp678_ext_wyb_calcuvalue(lri,lrk,isma,nlp_value)
      w0_sdold=w0_sdplp
      do iw0=1,mtype
        if(iw0.eq.1) goto 202
        w0_sdplp=vplpnew_w0(iw0)
        if(logic_dh) w0_sdplp=vplp_w0(iw0)
        if(abs(w0_sdplp).lt.crl) cycle
        w0multi=w0_sdplp/w0_sdold
        w0_sdold=w0_sdplp
        do iiext=1,nlp_value
         value_lpext(iiext)=value_lpext(iiext)*w0multi
        enddo
202     ilpsta=nstaval(iw0)+1
        ilpend=nstaval(iw0)+nvalue(iw0)
        do iplp=ilpsta,ilpend
          if(logic_dh) then
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            if(lin.eq.7)  call gsd_sequence_extspace(irw,ilw)
            if(lin.eq.14) call gtd_sequence_extspace(irw,ilw)
            if(lin.eq.26) call gdv_sequence_extspace(irw,ilw)
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
                if(lin.eq.7)  call gsd_sequence_extspace(irw,ilw)
                if(lin.eq.14) call gtd_sequence_extspace(irw,ilw)
                if(lin.eq.26) call gdv_sequence_extspace(irw,ilw)
              enddo
            enddo
          endif
        enddo
      enddo
      endif
      return
      end

      subroutine drr_br_ext_ar(lin,lrk,lri)
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
      data crl/1.0e-8/

      isma=lsm_inn(lri)
      iwuplwei=jpad_upwei(jpadl)
      ilsegdownwei=iseg_downwei(ipael)
      irsegdownwei=iseg_downwei(ipae)

      do iw0=1,mtype
        w0_sdplp=vplpnew_w0(iw0)
        if(logic_dh) w0_sdplp=vplp_w0(iw0)
        if(abs(w0_sdplp).lt.crl) cycle
        if(logic_grad) then
          call lp678_ext_calcuvalue_g(lri,lrk,isma,nlp_value)
          ilpsta=nstaval(iw0)+1
          ilpend=nstaval(iw0)+nvalue(iw0)
          do iplp=ilpsta,ilpend
            if(logic_dh) then
              ilw=lp_lwei(iplp)
              irw=lp_rwei(iplp)
              if(lin.eq.6)  call gsd_sequence_extspace_g(ilw,irw)
              if(lin.eq.13) call gtd_sequence_extspace_g(ilw,irw)
              if(lin.eq.23) call gdv_sequence_extspace_g(ilw,irw)
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
              if(lin.eq.6)  call gsd_sequence_extspace_g(ilw,irw)
              if(lin.eq.13) call gtd_sequence_extspace_g(ilw,irw)
              if(lin.eq.23) call gdv_sequence_extspace_g(ilw,irw)
                enddo
              enddo
            endif
          enddo

        else
          call lp678_ext_wyb_calcuvalue(lri,lrk,isma,nlp_value)
          ilpsta=nstaval(iw0)+1
          ilpend=nstaval(iw0)+nvalue(iw0)
          do iplp=ilpsta,ilpend
            if(logic_dh) then
              ilw=lp_lwei(iplp)
              irw=lp_rwei(iplp)
              if(lin.eq.6)  call gsd_sequence_extspace(ilw,irw)
              if(lin.eq.13) call gtd_sequence_extspace(ilw,irw)
              if(lin.eq.23) call gdv_sequence_extspace(ilw,irw)
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
                  if(lin.eq.6)  call gsd_sequence_extspace(ilw,irw)
                  if(lin.eq.13) call gtd_sequence_extspace(ilw,irw)
                  if(lin.eq.23) call gdv_sequence_extspace(ilw,irw)
                enddo
              enddo
            endif
          enddo
        endif
      enddo
      return
      end

      subroutine ar_br_br_ext_ar_new(lin,intentry,isma)
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"

      iwuplwei=jpad_upwei(jpadl)
      ilsegdownwei=iseg_downwei(ipael)
      irsegdownwei=iseg_downwei(ipae)
      do iw0=1,mtype
        w0_sdplp=vplpnew_w0(iw0)
        w1_sdplp=vplpnew_w1(iw0)
        if(logic_dh) w0_sdplp=vplp_w0(iw0)
        if(logic_dh) w1_sdplp=vplp_w1(iw0)
      if(logic_grad) then
          call lp10_arbrbr_ext_calcuvalue_g(intentry,isma,nlp_vlue)
          ilpsta=nstaval(iw0)+1
          ilpend=nstaval(iw0)+nvalue(iw0)
          do iplp=ilpsta,ilpend
            if(logic_dh) then
              ilw=lp_lwei(iplp)
              irw=lp_rwei(iplp)
              if(lin.eq.6)  call gsd_sequence_extspace_g(ilw,irw)
              if(lin.eq.13) call gtd_sequence_extspace_g(ilw,irw)
              if(lin.eq.23) call gdv_sequence_extspace_g(ilw,irw)
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
              if(lin.eq.6)  call gsd_sequence_extspace_g(ilw,irw)
              if(lin.eq.13) call gtd_sequence_extspace_g(ilw,irw)
              if(lin.eq.23) call gdv_sequence_extspace_g(ilw,irw)
                enddo
              enddo
            endif
          enddo

        else
          call lp10_arbrbr_ext_calcuvalue(intentry,isma,nlp_vlue)
          ilpsta=nstaval(iw0)+1
          ilpend=nstaval(iw0)+nvalue(iw0)
          do iplp=ilpsta,ilpend
            if(logic_dh) then
              ilw=lp_lwei(iplp)
              irw=lp_rwei(iplp)
              if(lin.eq.6)  call gsd_sequence_extspace(ilw,irw)
              if(lin.eq.13) call gtd_sequence_extspace(ilw,irw)
              if(lin.eq.23) call gdv_sequence_extspace(ilw,irw)
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
                  if(lin.eq.6)  call gsd_sequence_extspace(ilw,irw)
                  if(lin.eq.13) call gtd_sequence_extspace(ilw,irw)
                  if(lin.eq.23) call gdv_sequence_extspace(ilw,irw)
                enddo
              enddo
            endif
          enddo
        endif
      enddo
      return
      end

      subroutine ar_bl_bl_ext_ar_new(lin,intentry,isma,nk)
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"

      iwuplwei=jpad_upwei(jpadl)
      ilsegdownwei=iseg_downwei(ipael)
      irsegdownwei=iseg_downwei(ipae)
      do iw0=1,mtype
        w0_sdplp=vplpnew_w0(iw0)
        w1_sdplp=vplpnew_w1(iw0)
        if(logic_dh) w0_sdplp=vplp_w0(iw0)
        if(logic_dh) w1_sdplp=vplp_w1(iw0)
        if(logic_grad) then
          call lp12_arblbl_ext_calcuvalue_g(intentry,isma,nlp_value)
          ilpsta=nstaval(iw0)*nk+1
          ilpend=(nstaval(iw0)+nvalue(iw0))*nk
          do iplp=ilpsta,ilpend
            if(logic_dh) then                     !lp_head is in dbl_spa
              ilw=lp_lwei(iplp)
              irw=lp_rwei(iplp)
              if(lin.eq.6)  call gsd_sequence_extspace_g(ilw,irw)
              if(lin.eq.13) call gtd_sequence_extspace_g(ilw,irw)
              if(lin.eq.23) call gdv_sequence_extspace_g(ilw,irw)
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
              if(lin.eq.6)  call gsd_sequence_extspace_g(ilw,irw)
              if(lin.eq.13) call gtd_sequence_extspace_g(ilw,irw)
              if(lin.eq.23) call gdv_sequence_extspace_g(ilw,irw)
                enddo
              enddo
            endif
          enddo

        else
          call lp12_arblbl_ext_calcuvalue(intentry,isma,nlp_value)
          ilpsta=nstaval(iw0)*nk+1
          ilpend=(nstaval(iw0)+nvalue(iw0))*nk
          do iplp=ilpsta,ilpend
            if(logic_dh) then                     !lp_head is in dbl_spa
              ilw=lp_lwei(iplp)
              irw=lp_rwei(iplp)
              if(lin.eq.6)  call gsd_sequence_extspace(ilw,irw)
              if(lin.eq.13) call gtd_sequence_extspace(ilw,irw)
              if(lin.eq.23) call gdv_sequence_extspace(ilw,irw)
            else                                    !lp_head is in act_s
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
                  if(lin.eq.6)  call gsd_sequence_extspace(ilw,irw)
                  if(lin.eq.13) call gtd_sequence_extspace(ilw,irw)
                  if(lin.eq.23) call gdv_sequence_extspace(ilw,irw)
                enddo
              enddo
            endif
          enddo
        endif
      enddo
      return
      end

      subroutine drl_br_sum_al_new(lin,lrp,lrq,lri)   !to be reviced
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
      data crl/1.0e-8/
      if(logic_grad) then
      do lrk=1,norb_dz
        if(lrp.eq.lrk) cycle
        if(lrq.eq.lrk) cycle
        ilsegdownwei=iseg_downwei(ipae)
        irsegdownwei=iseg_downwei(ipael)
        isma=lsm_inn(lri)
        w0_sdplp=vplp_w0(1)
        call lp8_drlbr_sum_calcuvalue_g(lri,lrk,isma,nlp_value)
        w0_sdold=w0_sdplp
        do iw0=1,mtype
          if(iw0.eq.1) goto 201
          w0_sdplp=vplp_w0(iw0)
          if(abs(w0_sdplp).lt.crl) cycle
          w0multi=w0_sdplp/w0_sdold
          w0_sdold=w0_sdplp
          do iiext=1,nlp_value
           value_lpext(iiext)=value_lpext(iiext)*w0multi
           value_lpext1(iiext)=value_lpext1(iiext)*w0multi
          enddo
201       ilpsta=nstaval(iw0)+1
          ilpend=nstaval(iw0)+nvalue(iw0)
          do iplp=ilpsta,ilpend
              ilw=lp_lwei(iplp)
              irw=lp_rwei(iplp)
              if(lin.eq.7)  call gsd_sequence_extspace_g(irw,ilw)
              if(lin.eq.14) call gtd_sequence_extspace_g(irw,ilw)
              if(lin.eq.26) call gdv_sequence_extspace_g(irw,ilw)
          enddo
        enddo
      enddo

      else
        ilsegdownwei=iseg_downwei(ipae)
        irsegdownwei=iseg_downwei(ipael)
        isma=lsm_inn(lri)
        w0_sdplp=vplp_w0(1)
        call lp8_drlbr_sum_calcuvalue_wyb(lri,lrp,lrq,isma,nlp_value)
        w0_sdold=w0_sdplp
        do iw0=1,mtype
          if(iw0.eq.1) goto 202
          w0_sdplp=vplp_w0(iw0)
          if(abs(w0_sdplp).lt.crl) cycle
          w0multi=w0_sdplp/w0_sdold
          w0_sdold=w0_sdplp
          do iiext=1,nlp_value
           value_lpext(iiext)=value_lpext(iiext)*w0multi
          enddo
202     ilpsta=nstaval(iw0)+1
        ilpend=nstaval(iw0)+nvalue(iw0)
        do iplp=ilpsta,ilpend
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            if(lin.eq.7)  call gsd_sequence_extspace(irw,ilw)
            if(lin.eq.14) call gtd_sequence_extspace(irw,ilw)
            if(lin.eq.26) call gdv_sequence_extspace(irw,ilw)
        enddo
      enddo
      endif

      return
      end


      subroutine drl_bl_sum_ar_new(lin,lrp,lrq,lri)
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
      data crl/1.0e-8/
      if(logic_grad) then
        do lrk=1,norb_dz
          if(lrp.eq.lrk) cycle
          if(lrq.eq.lrk) cycle
          isma=lsm_inn(lri)
          ilsegdownwei=iseg_downwei(ipael)
          irsegdownwei=iseg_downwei(ipae)
          w0_sdplp=vplp_w0(1)
          call lp9_drlbl_sum_calcuvalue_g(lri,lrk,isma,nlp_value)
          w0_sdold=w0_sdplp
          do iw0=1,mtype
            if(iw0.eq.1) goto 201
            w0_sdplp=vplp_w0(iw0)
            if(abs(w0_sdplp).lt.crl) cycle
            w0multi=w0_sdplp/w0_sdold
            w0_sdold=w0_sdplp
            do iiext=1,nlp_value
             value_lpext(iiext)=value_lpext(iiext)*w0multi
             value_lpext1(iiext)=value_lpext1(iiext)*w0multi
            enddo
201         ilpsta=nstaval(iw0)+1
            ilpend=nstaval(iw0)+nvalue(iw0)
            do iplp=ilpsta,ilpend
                ilw=lp_lwei(iplp)
                irw=lp_rwei(iplp)
                if(lin.eq.6)  call gsd_sequence_extspace_g(ilw,irw)
                if(lin.eq.13) call gtd_sequence_extspace_g(ilw,irw)
                if(lin.eq.23) call gdv_sequence_extspace_g(ilw,irw)
            enddo
          enddo
        enddo

      else
        isma=lsm_inn(lri)
        ilsegdownwei=iseg_downwei(ipael)
        irsegdownwei=iseg_downwei(ipae)
        w0_sdplp=vplp_w0(1)

        call lp9_drlbl_sum_calcuvalue_wyb(lri,lrp,lrq,isma,nlp_value)
        w0_sdold=w0_sdplp
        do iw0=1,mtype
          if(iw0.eq.1) goto 202
          w0_sdplp=vplp_w0(iw0)
          if(abs(w0_sdplp).lt.crl) cycle
          w0multi=w0_sdplp/w0_sdold
          w0_sdold=w0_sdplp
          do iiext=1,nlp_value
           value_lpext(iiext)=value_lpext(iiext)*w0multi
          enddo
202       ilpsta=nstaval(iw0)+1
          ilpend=nstaval(iw0)+nvalue(iw0)
          do iplp=ilpsta,ilpend
              ilw=lp_lwei(iplp)
              irw=lp_rwei(iplp)
              if(lin.eq.6)  call gsd_sequence_extspace(ilw,irw)
              if(lin.eq.13) call gtd_sequence_extspace(ilw,irw)
              if(lin.eq.23) call gdv_sequence_extspace(ilw,irw)
          enddo
        enddo
      endif
      return
      end

      subroutine ar_bl_ext_ss(lri,lrj,nk)
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
      data crl/1.0e-8/

      iwuplwei=jpad_upwei(jpadl)
      ilsegdownwei=iseg_downwei(ipael)
      irsegdownwei=iseg_downwei(ipae)
      imlr=mul_tab(iml,imr)
      if (imlr.eq.1) then
        logic_g1415=.true.
        if(iml.eq.1) logic_g13=.true.
        if(iml.eq.1) logic_g2g4b=.true.
        logic_g36b=.true.
        logic_g35b=.true.
        logic_g34b=.true.
      endif

      w0_plp=vplpnew_w0(1)
      w1_plp=0.d0
      if(logic_dh) w0_plp=vplp_w0(1)

      if(logic_grad) then
      call lp_arbl_ext_st_calcuvalue_g(lri,lrj,nlp_value)
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

201     ilpsta=nstaval(iw0)*nk+1
        ilpend=(nstaval(iw0)+nvalue(iw0))*nk
        do iplp=ilpsta,ilpend
          if(logic_dh) then                    !lp_head is in dbl_space
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call inn_ext_ss_loop_unpack_g(ilw,irw)
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
                call inn_ext_ss_loop_unpack_g(ilw,irw)
              enddo
            enddo
          endif
        enddo
      enddo
      else
      call lp_arbl_ext_st_calcuvalue(lri,lrj,nlp_value)
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

202     ilpsta=nstaval(iw0)*nk+1
        ilpend=(nstaval(iw0)+nvalue(iw0))*nk
        do iplp=ilpsta,ilpend
          if(logic_dh) then                    !lp_head is in dbl_space
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call inn_ext_ss_loop_unpack(ilw,irw)
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
                call inn_ext_ss_loop_unpack(ilw,irw)
              enddo
            enddo
          endif
        enddo
      enddo
      endif

      return
      end

      subroutine ar_bl_ext_st(lri,lrj,nk)     !w0=0
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
      data crl/1.0e-8/
      iwuplwei=jpad_upwei(jpadl)
      ilsegdownwei=iseg_downwei(ipael)
      irsegdownwei=iseg_downwei(ipae)
      w0_plp=0.d0
      w1_plp=vplpnew_w1(1)
      if(logic_dh) w1_plp=vplp_w1(1)

      if(logic_grad) then
      call lp_arbl_ext_st_calcuvalue_g(lri,lrj,nlp_value)
      w1_old=w1_plp
      do iw0=1,mtype
        if(iw0.eq.1) goto 201
        w1_plp=vplpnew_w1(iw0)
        if(logic_dh) w1_plp=vplp_w1(iw0)
        if(abs(w1_plp).lt.crl) cycle
        w0multi=w1_plp/w1_old
        w1_old=w1_plp
        do iiext=1,nlp_value
         value_lpext(iiext)=value_lpext(iiext)*w0multi
         value_lpext1(iiext)=value_lpext1(iiext)*w0multi
        enddo

201     ilpsta=nstaval(iw0)*nk+1
        ilpend=(nstaval(iw0)+nvalue(iw0))*nk
        do iplp=ilpsta,ilpend
          if(logic_dh) then                    !lp_head is in dbl_space
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call inn_ext_st_loop_unpack_g(ilw,irw)
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
                call inn_ext_st_loop_unpack_g(ilw,irw)
              enddo
            enddo
          endif
        enddo
      enddo
      else
      call lp_arbl_ext_st_calcuvalue(lri,lrj,nlp_value)
      w1_old=w1_plp
      do iw0=1,mtype
        if(iw0.eq.1) goto 202
        w1_plp=vplpnew_w1(iw0)
        if(logic_dh) w1_plp=vplp_w1(iw0)
        if(abs(w1_plp).lt.crl) cycle
        w0multi=w1_plp/w1_old
        w1_old=w1_plp
        do iiext=1,nlp_value
         value_lpext(iiext)=value_lpext(iiext)*w0multi
        enddo

202     ilpsta=nstaval(iw0)*nk+1
        ilpend=(nstaval(iw0)+nvalue(iw0))*nk
        do iplp=ilpsta,ilpend
          if(logic_dh) then                    !lp_head is in dbl_space
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call inn_ext_st_loop_unpack(ilw,irw)
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
                call inn_ext_st_loop_unpack(ilw,irw)
              enddo
            enddo
          endif
        enddo
      enddo
      endif

      return
      end

      subroutine ar_bl_ext_ts(lri,lrj,nk)
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
      data crl/1.0e-8/
      iwuplwei=jpad_upwei(jpadl)
      ilsegdownwei=iseg_downwei(ipael)
      irsegdownwei=iseg_downwei(ipae)
      w0_plp=0.d0
      w1_plp=vplpnew_w1(1)
      if(logic_dh) w1_plp=vplp_w1(1)

      if(logic_grad) then
      call lp_arbl_ext_st_calcuvalue_g(lri,lrj,nlp_value)
      w1_old=w1_plp
      do iw0=1,mtype
        if(iw0.eq.1) goto 201
        w1_plp=vplpnew_w1(iw0)
        if(logic_dh) w1_plp=vplp_w1(iw0)
        if(abs(w1_plp).lt.crl) cycle
        w0multi=w1_plp/w1_old
        w1_old=w1_plp
        do iiext=1,nlp_value
         value_lpext(iiext)=value_lpext(iiext)*w0multi
         value_lpext1(iiext)=value_lpext1(iiext)*w0multi
        enddo

201     ilpsta=nstaval(iw0)*nk+1
        ilpend=(nstaval(iw0)+nvalue(iw0))*nk
        do iplp=ilpsta,ilpend
          if(logic_dh) then                    !lp_head is in dbl_space
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call inn_ext_ts_loop_unpack_g(ilw,irw)
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
                call inn_ext_ts_loop_unpack_g(ilw,irw)
              enddo
            enddo
          endif
        enddo
      enddo
      else
      call lp_arbl_ext_st_calcuvalue(lri,lrj,nlp_value)
      w1_old=w1_plp
      do iw0=1,mtype
        if(iw0.eq.1) goto 202
        w1_plp=vplpnew_w1(iw0)
        if(logic_dh) w1_plp=vplp_w1(iw0)
        if(abs(w1_plp).lt.crl) cycle
        w0multi=w1_plp/w1_old
        w1_old=w1_plp
        do iiext=1,nlp_value
         value_lpext(iiext)=value_lpext(iiext)*w0multi
        enddo

202     ilpsta=nstaval(iw0)*nk+1
        ilpend=(nstaval(iw0)+nvalue(iw0))*nk
        do iplp=ilpsta,ilpend
          if(logic_dh) then                    !lp_head is in dbl_space
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call inn_ext_ts_loop_unpack(ilw,irw)
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
                call inn_ext_ts_loop_unpack(ilw,irw)
              enddo
            enddo
          endif
        enddo
      enddo
      endif

      return
      end

      subroutine ar_bl_ext_tt(lri,lrj,nk)
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
      iwuplwei=jpad_upwei(jpadl)
      ilsegdownwei=iseg_downwei(ipael)
      irsegdownwei=iseg_downwei(ipae)
      imlr=mul_tab(iml,imr)
      if (imlr.eq.1) then
        logic_g1415=.true.
        logic_g36b=.true.
        logic_g35b=.true.
        logic_g34b=.true.
      endif
      do iw0=1,mtype
        w0_plp=vplpnew_w0(iw0)
        w1_plp=vplpnew_w1(iw0)
        if(logic_dh) w0_plp=vplp_w0(iw0)
        if(logic_dh) w1_plp=vplp_w1(iw0)

      if(logic_grad) then
        call lp_arbl_ext_st_calcuvalue_g(lri,lrj,nlp_value)
        ilpsta=nstaval(iw0)*nk+1
        ilpend=(nstaval(iw0)+nvalue(iw0))*nk
        do iplp=ilpsta,ilpend
          if(logic_dh) then                    !lp_head is in dbl_space
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call inn_ext_tt_loop_unpack_g(ilw,irw)
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
                call inn_ext_tt_loop_unpack_g(ilw,irw)
              enddo
            enddo
          endif
        enddo
      else
        call lp_arbl_ext_st_calcuvalue(lri,lrj,nlp_value)
        ilpsta=nstaval(iw0)*nk+1
        ilpend=(nstaval(iw0)+nvalue(iw0))*nk
        do iplp=ilpsta,ilpend
          if(logic_dh) then                    !lp_head is in dbl_space
            ilw=lp_lwei(iplp)
            irw=lp_rwei(iplp)
            call inn_ext_tt_loop_unpack(ilw,irw)
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
                call inn_ext_tt_loop_unpack(ilw,irw)
              enddo
            enddo
          endif
        enddo

       endif
      enddo
      return
      end

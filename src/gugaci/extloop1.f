************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2009, Bingbing Suo                                     *
************************************************************************
C Jul. 3, 2009 -BSUO- External space loops
      subroutine dd_ext_plpmode(ilnodesm,irnodesm)
#include "drt_h.fh"
#include "lpextmode_h.fh"
      logic_g36a=.false.
      logic_g36b=.false.
      logic_g1415=.false.
      logic_g13=.false.
      if (ilnodesm.lt.irnodesm) then
         logic_g36a=.true.
      elseif ( ilnodesm.eq.irnodesm) then
         logic_g36a=.true.
         logic_g36b=.true.
         logic_g1415=.true.
         logic_g13=.true.
      else
         logic_g36b=.true.
      endif
      end

      subroutine external_space_plpmode_value_dv()
#include "drt_h.fh"
#include "paraconstants_h.fh"
#include "lpextmode_h.fh"
!     sd   lpmode_value
      w0g25=-1.d0
      w0g25a=-1.d0
      w1g25a=-1.d0
      end

      subroutine external_space_plpmode_value_vd()
#include "drt_h.fh"
#include "paraconstants_h.fh"
#include "lpextmode_h.fh"
!     sd   lpmode_value
      w0g25=-v_sqtwo
      w0g25a=-v_sqtwo
      w1g25a=-v_sqtwo
      end

      subroutine g12_t_diffsym(isma,ismb,ismc)
#include "drt_h.fh"
#include "intsort_h.fh"

      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      if ( isma .lt. ismb ) then
      ip3smabc=isma+jp2(ismb)+jp3(ismc)
      ip2cd=m_jc+(m_jd-1)*nlsm_ext(ismc)   !severe_new_error_1206
      num_smab=nlsm_ext(isma)*nlsm_ext(ismb)  !need checking
      ipos_intbasetmp=ip4_abcd_ext_base(ip3smabc)+
     *            (ip2cd-1)*num_smab*3
      ibsta=ibsm_ext(ismb)
      ibend=iesm_ext(ismb)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)
      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)    !need checking
      iposint=ipos_intbasetmp
      do ib=ibsta,ibend
      do ia=iasta,iaend
         value_lpext(ilwei)=vint_ci(iposint+1)-vint_ci(iposint+2)
         ilwei=ilwei+1
         iposint=iposint+3
      enddo
      enddo
      else
      ip3smabc=isma+jp2(ismb)+jp3(ismc)
      ip2cd=m_jc+NGW2(m_jd)
      nsma=nlsm_ext(isma)
      num_smab=nsma*(nsma-1)/2   !need checking
      ipos_intbasetmp=ip4_abcd_ext_base(ip3smabc)+
     *            (ip2cd-1)*num_smab*3
      ibsta=ibsm_ext(ismb)+1
      ibend=iesm_ext(ismb)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)
      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)    !need checking
      iposint=ipos_intbasetmp
      do ib=ibsta,ibend
      do ia=iasta,ib-1
         value_lpext(ilwei)=vint_ci(iposint+1)-vint_ci(iposint+2)
         ilwei=ilwei+1
         iposint=iposint+3
      enddo
      enddo
      endif
      end

      subroutine g11a_t_diffsym(isma,ismb,ismc)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      ipsmabc=isma+jp2(ismc)+jp3(ismb)
      ipos_intbasetmp=ip4_abcd_ext_base(ipsmabc)
      jdoffset=(m_jd-1)*nlsm_ext(ismb)
      jcoffset=(m_jc-1)*nlsm_ext(isma)
      num_smab=nlsm_ext(isma)*nlsm_ext(ismc)*3
      ibsta=ibsm_ext(ismb)
      ibend=iesm_ext(ismb)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)
      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)    !need checking
      ipos_intbasetmp=ipos_intbasetmp+jdoffset*num_smab+jcoffset*3
      do ib=ibsta,ibend
         iposint=ipos_intbasetmp
      do ia=iasta,iaend
         value_lpext(ilwei)=vint_ci(iposint)-vint_ci(iposint+2)
         iposint=iposint+3
         ilwei=ilwei+1
      enddo
         ipos_intbasetmp=ipos_intbasetmp+num_smab
      enddo
      end

      subroutine g11b_t_diffsym(isma,ismb,ismc)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      ipsmabc=ismc+jp2(isma)+jp3(ismb)
      ipos_intbasetmp=ip4_abcd_ext_base(ipsmabc)
      jdoffset=(m_jd-1)*nlsm_ext(ismb)
      jcoffset=m_jc-1
      num_smab=nlsm_ext(isma)*nlsm_ext(ismc)*3
      numint_jc=nlsm_ext(ismc)*3
      ibsta=ibsm_ext(ismb)
      ibend=iesm_ext(ismb)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)
      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)    !need checking
      ipos_intbasetmp=ipos_intbasetmp+jdoffset*num_smab+jcoffset*3
      do ib=ibsta,ibend
         iposint=ipos_intbasetmp      !+jcoffset*3
      do ia=iasta,iaend
         value_lpext(ilwei)=vint_ci(iposint)-vint_ci(iposint+1)
         ilwei=ilwei+1
         iposint=iposint+numint_jc
      enddo
         ipos_intbasetmp=ipos_intbasetmp+num_smab
      enddo
      end

      subroutine g1112_t_symaaaa(isma,ic,id)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      ipsmabc=isma+jp2(isma)+jp3(isma)
      ipos_intbasetmp=ip4_abcd_ext_base(ipsmabc)
!
      icdpos_12=NGW4(m_jd)+NGW3(m_jc)
      iasta=ibsm_ext(isma)
      ibsta=iasta+1
      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)    !need checking
      iposint=ipos_intbasetmp+icdpos_12*3         !-3 unnecessary ?
      do   ib=ibsta,ic-1
         do ia=iasta,ib-1
            value_lpext(ilwei)=vint_ci(iposint+1)-vint_ci(iposint+2)
            ilwei=ilwei+1
            iposint=iposint+3
         enddo
      enddo

      icdpos_11a=NGW4(m_jd)+NGW2(m_jc)
      jb=m_jc
      do ib=ic+1,id-1
         jb=jb+1
         ibcdpos_11a=icdpos_11a+NGW3(jb)
         iposint=ipos_intbasetmp+ibcdpos_11a*3
         ilwei=icnt_base+iwt_orb_ext(iasta,ib)
         do ia=iasta,ic-1
            value_lpext(ilwei)=vint_ci(iposint)-vint_ci(iposint+2)
            ilwei=ilwei+1
            iposint=iposint+3
         enddo
      enddo

      icdpos_11b=NGW4(m_jd)+m_jc
      jb=m_jc+1
      do ib=ic+2,id-1
         jb=jb+1
         ibcdpos_11b=icdpos_11b+NGW3(jb)
         ja=m_jc
         ilwei=icnt_base+iwt_orb_ext(ic+1,ib)
         do ia=ic+1,ib-1
            ja=ja+1
            iabcdpos_11b=ibcdpos_11b+NGW2(ja)
            iposint=ipos_intbasetmp+iabcdpos_11b*3-3
            value_lpext(ilwei)=vint_ci(iposint)-vint_ci(iposint+1)
            ilwei=ilwei+1
         enddo
      enddo
      end

      subroutine g11a11b_t_symaacc(isma,ismc,ic,id)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      ipsmabc=isma+jp2(isma)+jp3(ismc)
      ipos_intbasetmp=ip4_abcd_ext_base(ipsmabc)
      num_smab=nlsm_ext(isma)
      num_smab=(num_smab*(num_smab-1))/2*3
      ibsta=ibsm_ext(ismc)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)

      ibdpos=NGW2(m_jd)
      ipos_intbasetmp=ipos_intbasetmp+
     *            ibdpos*num_smab+NGW2(m_jc)*3
      do ib=ibsta,id-1
         iposint=ipos_intbasetmp
         ilwei=icnt_base+iwt_orb_ext(iasta,ib)
         do ia=iasta,ic-1
            value_lpext(ilwei)=vint_ci(iposint)-vint_ci(iposint+2)
            ilwei=ilwei+1
            iposint=iposint+3
         enddo
         ipos_intbasetmp=ipos_intbasetmp+num_smab
      enddo

      ipos_intbasetmp=ip4_abcd_ext_base(ipsmabc)
      ibdpos=NGW2(m_jd)
      ipos_intbasetmp=ipos_intbasetmp+
     *            ibdpos*num_smab+m_jc*3
      do ib=ibsta,id-1
         iposint=ipos_intbasetmp
         ja=m_jc
         ilwei=icnt_base+iwt_orb_ext(ic+1,ib)
         do ia=ic+1,iaend
            ja=ja+1
            iposint=ipos_intbasetmp+NGW2(ja)*3-3
            value_lpext(ilwei)=vint_ci(iposint)-vint_ci(iposint+1)
            ilwei=ilwei+1
         enddo
         ipos_intbasetmp=ipos_intbasetmp+num_smab
      enddo
      end

      subroutine g36_t_ext(ismc,ic,id)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      ipos_intbasetmp=ip3_abd_ext_base+(id-1)*np3_abd_ext   !severe_new_
      iasta=ibsm_ext(ismc)
      ilwei=icnt_base+iwt_orb_ext(iasta,id)
      do ia=iasta,ic-1
         ipos_g36=iwt_orb_ext(ia,ic)-1
         iposint=ipos_intbasetmp+ipos_g36*2         !severe_new_error
         value_lpext(ilwei)=vint_ci(iposint+1)-vint_ci(iposint)
     *            +vint_ci(ip2_aa_ext_base+ipos_g36)
         ilwei=ilwei+1
      enddo
      end

      subroutine g5_t_ext(ismd,ic,id)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      ipos_intbasetmp=ip3_abd_ext_base+(ic-1)*np3_abd_ext
      ibsta=ibsm_ext(ismd)
      do ib=max(ic+1,ibsta),id-1
         ipos_g5=iwt_orb_ext(ib,id)-1            !severe_new_error
         iposint=ipos_intbasetmp+ipos_g5*2         !severe_new_error
         ilwei=icnt_base+iwt_orb_ext(ic,ib)
         value_lpext(ilwei)=vint_ci(iposint+1)-vint_ci(iposint)
     *            +vint_ci(ip2_aa_ext_base+ipos_g5)
      enddo
      end

      subroutine   g9_t_ext(ismd,ic,id)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      ipos_intbasetmp=ip3_abd_ext_base+(ic-1)*np3_abd_ext
      iasta=ibsm_ext(ismd)
      ilwei=icnt_base+iwt_orb_ext(iasta,ic)
      do ia=iasta,ic-1                     !severe_new_error
         ipos_g9=iwt_orb_ext(ia,id)-1
         iposint=ipos_intbasetmp+ipos_g9*2      !severe_new_error
         value_lpext(ilwei)=vint_ci(iposint)-vint_ci(iposint+1)
     *            -vint_ci(ip2_aa_ext_base+ipos_g9)         !severe_new_
         ilwei=ilwei+1
      enddo
      end

      subroutine gsd_determine_extarmode_paras(
     *                           ismnodes,ismnoded,logic_sd)
#include "drt_h.fh"
#include "lpextmode_h.fh"
      logical logic_sd
      ismnodesd=mul_tab(ismnodes,ismnoded)
      numsmd=nlsm_ext(ismnoded)
      numsmsd=nlsm_ext(ismnodesd)
      iorbid=ibsm_ext(ismnoded)
      iorbisd=ibsm_ext(ismnodesd)

      logic_g25a=.false.
      logic_g25b=.false.
      logic_g28a=.false.
      logic_g26=.false.
      if ( ismnoded .gt.ismnodesd ) then
         logic_g28a=.true.
         iweista_g28=iwt_orb_ext(iorbisd,iorbid)
         nwei_g28=numsmd
         nint_g28=numsmsd
      elseif ( ismnoded .eq. ismnodesd ) then
         logic_g25b=.true.
         iweista_g25=iwt_orb_ext(iorbisd,iorbisd+1)
         nwei_g25=numsmd
         nint_g25=numsmsd
         iweista_g28=iwt_orb_ext(iorbisd,iorbisd+1)
         nwei_g28=numsmd
         nint_g28=numsmsd
      else
         logic_g25a=.true.
         iweista_g25=iwt_orb_ext(iorbid,iorbisd)
         nwei_g25=numsmd
         nint_g25=numsmsd
      endif
      if ( ismnodes .eq. 1 .and. logic_sd ) then
         logic_g26=.true.
         iweista_g26=iwt_sm_s_ext+iorbid
         nwei_g26=numsmd
         ivaluesta_g26=0
      endif
      end

      subroutine g_dd_ext_sequence(ism)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      icano_nnsta=2
      icnt_base=0
      iasta=ibsm_ext(ism)
      ibsta=iasta+1
      ibend=iesm_ext(ism)
      ilwei=0
      do ib=ibsta,ibend
         lrb=norb_number(ib)
         do ia=iasta,ib-1
           lra=norb_number(ia)
           ilwei=ilwei+1
          value_lpext(ilwei)=voint(lrb,lra)
         enddo
      enddo
      icano_nnend=ibend-iasta+1
      call complete_ext_loop()
      end


      subroutine g_ss_ext_sequence(ism,itype)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      icano_nnsta=2
      icnt_base=0
      do ismd=1,ng_sm
         ismc=mul_tab(ism,ismd)
         if ( ismc .gt. ismd ) cycle
         id_sta=ibsm_ext(ismd)
         idsta=id_sta
         idend=iesm_ext(ismd)
         ic_sta=ibsm_ext(ismc)
         icend=iesm_ext(ismc)
         if ( ismd.eq.ismc ) idsta=idsta+1
         do id=idsta,idend
            m_jd=id-id_sta+1
         do ic=ic_sta,min(icend,id-1)
            m_jc=ic-ic_sta+1
            icano_nn=iwt_orb_ext(ic,id)
            if ( icnt_base+icano_nn-1.gt. max_tmpvalue ) then
               call complete_ext_loop()
               icnt_base=0
               icano_nnsta=icano_nn
            endif
            icano_nnend=icano_nn
            do ismb=1,ismd-1
               isma=mul_tab(ism,ismb)
               if ( isma .gt. ismb ) cycle
               if ( ismc .gt. ismb ) then
                  call g12_diffsym(isma,ismb,ismc)
               elseif ( ismc .gt. isma ) then
                  call g11a_diffsym(isma,ismb,ismc)
               else
                  call g11b_diffsym(isma,ismb,ismc)
               endif
            enddo
            if ( ism.eq.1 )   then
               isma=ismd
               call g1112_symaaaa(isma,ic,id)
            else
               isma=mul_tab(ism,ismd)
               call g11a11b_symaacc(isma,ismd,ic,id)
            endif
            call g10_ext(ismc,ic,id)
            call g5_ext(ismd,ic,id)
            if ( ism.eq.1 ) call g9_ext(ismd,ic,id)
            icnt_base=icnt_base+icano_nn-1
         enddo
         enddo
      enddo
      if ( ism.eq.1 .and. itype.eq.4 ) then
      do id=1,norb_ext
         icano_nn=id+iwt_sm_s_ext
            if ( icnt_base+icano_nn-1.gt. max_tmpvalue ) then
               call complete_ext_loop()
               icnt_base=0
               icano_nnsta=icano_nn
            endif
            icano_nnend=icano_nn
            call ext_lp_ab_s1(id)
            icnt_base=icnt_base+icano_nn-1
      enddo
      endif
      call complete_ext_loop()
      end

      subroutine g12_diffsym(isma,ismb,ismc)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      if ( isma .lt. ismb ) then
      ip3smabc=isma+jp2(ismb)+jp3(ismc)
      ip2cd=m_jc+(m_jd-1)*nlsm_ext(ismc)      !severe_new_error_1206
      num_smab=nlsm_ext(isma)*nlsm_ext(ismb)  !need checking
      ipos_intbasetmp=ip4_abcd_ext_base(ip3smabc)+
     *            (ip2cd-1)*num_smab*3
      ibsta=ibsm_ext(ismb)
      ibend=iesm_ext(ismb)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)
      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)    !need checking
      iposint=ipos_intbasetmp
      do ib=ibsta,ibend
      do ia=iasta,iaend
         value_lpext(ilwei)=vint_ci(iposint+1)+vint_ci(iposint+2)
         ilwei=ilwei+1
         iposint=iposint+3
      enddo
      enddo
      else
      ip3smabc=isma+jp2(ismb)+jp3(ismc)
      ip2cd=m_jc+NGW2(m_jd)
      nsma=nlsm_ext(isma)
      num_smab=nsma*(nsma-1)/2   !need checking
      ipos_intbasetmp=ip4_abcd_ext_base(ip3smabc)+
     *            (ip2cd-1)*num_smab*3
      ibsta=ibsm_ext(ismb)+1
      ibend=iesm_ext(ismb)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)
      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)    !need checking
      iposint=ipos_intbasetmp
      do ib=ibsta,ibend
      do ia=iasta,ib-1
         value_lpext(ilwei)=vint_ci(iposint+1)+vint_ci(iposint+2)
         ilwei=ilwei+1
         iposint=iposint+3
      enddo
      enddo
      endif
      end

      subroutine g11a_diffsym(isma,ismb,ismc)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      ipsmabc=isma+jp2(ismc)+jp3(ismb)
      ipos_intbasetmp=ip4_abcd_ext_base(ipsmabc)
      jdoffset=(m_jd-1)*nlsm_ext(ismb)
      jcoffset=(m_jc-1)*nlsm_ext(isma)
      num_smab=nlsm_ext(isma)*nlsm_ext(ismc)*3
      ibsta=ibsm_ext(ismb)
      ibend=iesm_ext(ismb)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)
      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)    !need checking
      ipos_intbasetmp=ipos_intbasetmp+jdoffset*num_smab+jcoffset*3
      do ib=ibsta,ibend
         iposint=ipos_intbasetmp
      do ia=iasta,iaend
         value_lpext(ilwei)=vint_ci(iposint)+vint_ci(iposint+2)
         iposint=iposint+3
         ilwei=ilwei+1
      enddo
         ipos_intbasetmp=ipos_intbasetmp+num_smab
      enddo
      end

      subroutine g11b_diffsym(isma,ismb,ismc)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      ipsmabc=ismc+jp2(isma)+jp3(ismb)
      ipos_intbasetmp=ip4_abcd_ext_base(ipsmabc)
      jdoffset=(m_jd-1)*nlsm_ext(ismb)
      jcoffset=m_jc-1
      num_smab=nlsm_ext(isma)*nlsm_ext(ismc)*3
      numint_jc=nlsm_ext(ismc)*3
      ibsta=ibsm_ext(ismb)
      ibend=iesm_ext(ismb)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)
      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)    !need checking
      ipos_intbasetmp=ipos_intbasetmp+jdoffset*num_smab+jcoffset*3
      do ib=ibsta,ibend
         iposint=ipos_intbasetmp      !+jcoffset*3
      do ia=iasta,iaend
         value_lpext(ilwei)=vint_ci(iposint)+vint_ci(iposint+1)
         ilwei=ilwei+1
         iposint=iposint+numint_jc
      enddo
         ipos_intbasetmp=ipos_intbasetmp+num_smab
      enddo
      end

      subroutine g1112_symaaaa(isma,ic,id)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      ipsmabc=isma+jp2(isma)+jp3(isma)
      ipos_intbasetmp=ip4_abcd_ext_base(ipsmabc)
!
      icdpos_12=NGW4(m_jd)+NGW3(m_jc)
      iasta=ibsm_ext(isma)
      ibsta=iasta+1
      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)    !need checking
      iposint=ipos_intbasetmp+icdpos_12*3         !-3 unnecessary ?
      do   ib=ibsta,ic-1
!        jb=jb+1
         do ia=iasta,ib-1
            value_lpext(ilwei)=vint_ci(iposint+1)+vint_ci(iposint+2)
            ilwei=ilwei+1
            iposint=iposint+3
         enddo
      enddo

      icdpos_11a=NGW4(m_jd)+NGW2(m_jc)
      jb=m_jc
      do ib=ic+1,id-1
         jb=jb+1
         ibcdpos_11a=icdpos_11a+NGW3(jb)
         iposint=ipos_intbasetmp+ibcdpos_11a*3
         ilwei=icnt_base+iwt_orb_ext(iasta,ib)
         do ia=iasta,ic-1
            value_lpext(ilwei)=vint_ci(iposint)+vint_ci(iposint+2)
            ilwei=ilwei+1
            iposint=iposint+3
         enddo
      enddo

      icdpos_11b=NGW4(m_jd)+m_jc
      jb=m_jc+1
      do ib=ic+2,id-1
         jb=jb+1
         ibcdpos_11b=icdpos_11b+NGW3(jb)
         ja=m_jc
         ilwei=icnt_base+iwt_orb_ext(ic+1,ib)
         do ia=ic+1,ib-1
            ja=ja+1
            iabcdpos_11b=ibcdpos_11b+NGW2(ja)
            iposint=ipos_intbasetmp+iabcdpos_11b*3-3
            value_lpext(ilwei)=vint_ci(iposint)+vint_ci(iposint+1)
            ilwei=ilwei+1
         enddo
      enddo
      end

      subroutine g11a11b_symaacc(isma,ismc,ic,id)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      ipsmabc=isma+jp2(isma)+jp3(ismc)
      ipos_intbasetmp=ip4_abcd_ext_base(ipsmabc)
      num_smab=nlsm_ext(isma)
      num_smab=(num_smab*(num_smab-1))/2*3
      ibsta=ibsm_ext(ismc)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)

      ibdpos=NGW2(m_jd)
      ipos_intbasetmp=ipos_intbasetmp+
     *            ibdpos*num_smab+NGW2(m_jc)*3
      do ib=ibsta,id-1
         iposint=ipos_intbasetmp
         ilwei=icnt_base+iwt_orb_ext(iasta,ib)
         do ia=iasta,ic-1
            value_lpext(ilwei)=vint_ci(iposint)+vint_ci(iposint+2)
            ilwei=ilwei+1
            iposint=iposint+3
         enddo
         ipos_intbasetmp=ipos_intbasetmp+num_smab
      enddo

      ipos_intbasetmp=ip4_abcd_ext_base(ipsmabc)
      ibdpos=NGW2(m_jd)
      ipos_intbasetmp=ipos_intbasetmp+
     *            ibdpos*num_smab+m_jc*3
      do ib=ibsta,id-1
         iposint=ipos_intbasetmp
         ja=m_jc
         ilwei=icnt_base+iwt_orb_ext(ic+1,ib)
         do ia=ic+1,iaend
            ja=ja+1
            iposint=ipos_intbasetmp+NGW2(ja)*3-3
            value_lpext(ilwei)=vint_ci(iposint)+vint_ci(iposint+1)
            ilwei=ilwei+1
         enddo
         ipos_intbasetmp=ipos_intbasetmp+num_smab
      enddo
      end

      subroutine g10_ext(ismc,ic,id)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      ipos_intbasetmp=ip3_abd_ext_base+(id-1)*np3_abd_ext   !severe_new_
      iasta=ibsm_ext(ismc)
      ilwei=icnt_base+iwt_orb_ext(iasta,id)
      do ia=iasta,ic-1
         ipos_g10=iwt_orb_ext(ia,ic)-1
         iposint=ipos_intbasetmp+ipos_g10*2
         value_lpext(ilwei)=vint_ci(iposint+1)+vint_ci(iposint)
     *            +vint_ci(ip2_aa_ext_base+ipos_g10)
         ilwei=ilwei+1
      enddo
      end

      subroutine g5_ext(ismd,ic,id)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      ipos_intbasetmp=ip3_abd_ext_base+(ic-1)*np3_abd_ext
      ibsta=ibsm_ext(ismd)
      do ib=max(ic+1,ibsta),id-1
         ipos_g5=iwt_orb_ext(ib,id)-1            !severe_new_error
         iposint=ipos_intbasetmp+ipos_g5*2         !severe_new_error
         ilwei=icnt_base+iwt_orb_ext(ic,ib)
         value_lpext(ilwei)=vint_ci(iposint+1)+vint_ci(iposint)
     *            +vint_ci(ip2_aa_ext_base+ipos_g5)
      enddo
      end

      subroutine   g9_ext(ismd,ic,id)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      ipos_intbasetmp=ip3_abd_ext_base+(ic-1)*np3_abd_ext
      iasta=ibsm_ext(ismd)
      ilwei=icnt_base+iwt_orb_ext(iasta,ic)
      do ia=iasta,ic-1                     !severe_new_error
         ipos_g9=iwt_orb_ext(ia,id)-1
         iposint=ipos_intbasetmp+ipos_g9*2      !severe_new_error
         value_lpext(ilwei)=vint_ci(iposint+1)+vint_ci(iposint)
     *            +vint_ci(ip2_aa_ext_base+ipos_g9)
         ilwei=ilwei+1
      enddo
      end

      subroutine   ext_lp_ab_s1(id)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      parameter (   v_sqtwo=1.414213562373095d0 )
      ipos_intbasetmp=ip3_abd_ext_base+(id-1)*np3_abd_ext
      ilwei=icnt_base
      iposintaa=ip2_aa_ext_base-1
      iposint=ipos_intbasetmp
      do ismb=1,ng_sm
         iasta=ibsm_ext(ismb)
         ibend=iesm_ext(ismb)
         ibsta=iasta+1
         do ib=ibsta,ibend
         do ia=iasta,ib-1
            iposintaa=iposintaa+1
            if ( ib .eq. id .or. ia .eq. id ) then
!              g2   arar+ar(head)ar
               ilwei=ilwei+1
               value_lpext(ilwei)=
     *         v_sqtwo*(vint_ci(iposintaa)+vint_ci(iposint))
            else
!              g6g7g8
               ilwei=ilwei+1
               value_lpext(ilwei)=v_sqtwo*vint_ci(iposint)
            endif
            iposint=iposint+2
         enddo
         enddo
      enddo
!           g1
!     lw=lw0
      iposint=ip2_dd_ext_base+NGW2(id)
      do ic=1,id-1
         ilwei=ilwei+1
         value_lpext(ilwei)=vint_ci(iposint)
         iposint=iposint+1
      enddo
      end

      subroutine g31_diffsym(lri,isma,ismb)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei

      iabc0=(lri-1)*nabc
      ic=m_jd
      ibsta=ibsm_ext(ismb)
      ibend=iesm_ext(ismb)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)
      if ( isma .eq. ismb ) ibsta=ibsta+1
      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)
      do ib=ibsta,ibend
        if ( isma .eq. ismb ) iaend=ib-1
        do ia=iasta,iaend
!     g31   type_12 arbrb^ra^r
          iabc=iabc0+ia+ngw2(ib)+ngw3(ic)
          iposint=intind_iabc(iabc)
         value_lpext(ilwei)=vint_ci(iposint+1)*w0plp31
     *               +vint_ci(iposint+2)*w1plp31
         ilwei=ilwei+1
        enddo
      enddo
      end

      subroutine g32a_diffsym(lri,isma,ismb)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei

      iabc0=(lri-1)*nabc

      ic=m_jd

      ibsta=ibsm_ext(ismb)
      ibend=iesm_ext(ismb)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)
      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)
      do ib=ibsta,ibend
        do ia=iasta,iaend
          iabc=iabc0+ia+ngw2(ic)+ngw3(ib)
          iposint=intind_iabc(iabc)
         value_lpext(ilwei)=vint_ci(iposint+2)*w0plp32
     *                  -vint_ci(iposint)*w1plp32
         ilwei=ilwei+1
        enddo
      enddo
      end

      subroutine g32b_diffsym(lri,isma,ismb)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei

      iabc0=(lri-1)*nabc

      ic=m_jd

      ibsta=ibsm_ext(ismb)
      ibend=iesm_ext(ismb)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)
      if ( ismb.eq.isma ) ibsta=ibsta+1         !severe_new_error
      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)    !need checking
      do ib=ibsta,ibend
        do ia=iasta,min(iaend,ib-1)
      !g32b
          iabc=iabc0+ic+ngw2(ia)+ngw3(ib)
          iposint=intind_iabc(iabc)
         value_lpext(ilwei)=vint_ci(iposint+1)*w0plp32   !severe_new_err
     *                  -vint_ci(iposint)*w1plp32
         ilwei=ilwei+1
        enddo
      enddo
      end


      subroutine assign_segmode_paras(indl,indr,ilrdivnodesm)
#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "intsort_h.fh"
      ilsegdownwei=iseg_downwei(indl)
      irsegdownwei=iseg_downwei(indr)
      ilsegstawei=iseg_sta(indl)
      irsegstawei=iseg_sta(indr)
      ip2_intbase=ip2_ab_inn_base(ilrdivnodesm)
      ip2_intspace=int2ind_space_extsmab(ilrdivnodesm)
      ip2_intsymspace=int2ind_numijkl_extsmab(ilrdivnodesm)
      end

      subroutine external_space_plpmode_value_st()   !tt
#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "paraconstants_h.fh"
      w0g36a=0.d0
      w0g36b=0.d0
      w0g34a=0.d0
      w0g34b=0.d0
      w0g35a=0.d0
      w0g35b=0.d0
      w0g2a=0.d0
      w0g2b=0.d0
      w0g4a=0.d0
      w0g4b=0.d0
      w0g14a=0.d0
      w0g15a=0.d0
!     w0g13a=0.d0

      !st
      w1g36a=v_sqthreevsqtwo
      w1g36b=v_sqthreevsqtwo
      w1g34a=-w1g36a
      w1g34b=-w1g36a
      w1g35a=w1g36a
      w1g35b=-w1g36a
      w1g2a=0.d0      !g2a   ab_d
      w1g2b=-v_sqthree   !severe_error
      w1g4a=0.d0
      w1g4b=v_sqthree
      w1g14a=w1g36a
      w1g15a=-w1g36a
      end

      subroutine external_space_plpmode_value_ts()   !tt
#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "paraconstants_h.fh"
      w0g36a=0.d0
      w0g36b=0.d0
      w0g34a=0.d0
      w0g34b=0.d0
      w0g35a=0.d0
      w0g35b=0.d0
      w0g2a=0.d0
      w0g2b=0.d0
      w0g4a=0.d0
      w0g4b=0.d0
      w0g14a=0.d0
      w0g15a=0.d0
      !ts
      w1g36a=-v_onevsqtwo
      w1g36b=-v_onevsqtwo
      w1g34a=-w1g36a
      w1g34b=-w1g36a
      w1g35a=-w1g36a   !severe_error
      w1g35b=w1g36a
      w1g2a=1.d0      !g2a   ab_d
      w1g2b=0.d0
      w1g4a=-1.d0
      w1g4b=0.d0
      w1g14a=w1g36a
      w1g15a=-w1g36a
      end

      subroutine external_space_plpmode_value_ss()   !ss
#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "paraconstants_h.fh"
      w0g36a=-v_onevsqtwo
      w0g36b=-v_onevsqtwo
      w0g34a=w0g36a
      w0g34b=w0g36a
      w0g35a=w0g36a
      w0g35b=w0g36a
      w0g2a=-1.d0
      w0g2b=-1.d0
      w0g4a=-1.d0
      w0g4b=-1.d0
      w0g14a=w0g36a
      w0g15a=w0g36a
      w0g13a=-v_sqtwo
      w1g36a=0.d0
      w1g36b=0.d0
      w1g34a=0.d0
      w1g34b=0.d0
      w1g35a=0.d0
      w1g35b=0.d0
      w1g2a=0.d0
      w1g2b=0.d0
      w1g4a=0.d0
      w1g4b=0.d0
      w1g14a=0.d0
      w1g15a=0.d0
      end

      subroutine external_space_plpmode_value_tt()   !tt
#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "paraconstants_h.fh"
      w0g36a=-v_onevsqtwo
      w0g36b=-v_onevsqtwo
      w0g34a=w0g36a
      w0g34b=w0g36a
      w0g35a=-w0g36a
      w0g35b=-w0g36a
      w0g2a=0.d0
      w0g2b=0.d0
      w0g4a=0.d0
      w0g4b=0.d0
      w0g14a=w0g36a
      w0g15a=w0g36a
!     w0g13a=-v_sqtwo

      w1g36a=1.d0
      w1g36b=1.d0
      w1g34a=w1g36a
      w1g34b=w1g36a
      w1g35a=-w1g36a
      w1g35b=-w1g36a
      w1g2a=0.d0      !g2a   ab_d
      w1g2b=0.d0
      w1g4a=0.d0
      w1g4b=0.d0
      w1g14a=w1g36a
      w1g15a=w1g36a
      end

      subroutine external_space_plpmode_value_sd()
#include "drt_h.fh"
#include "paraconstants_h.fh"
#include "lpextmode_h.fh"
!     sd   lpmode_value
      w0g29=v_sqtwo
      w0g30=v_sqtwo
      w0g26a=v_sqtwo
      w1g26a=v_sqtwo
      w0g25=1.d0
      w0g25a=1.d0
      w1g25a=1.d0
      w0g27=1.d0
      w1g27=-1.d0
      w0g28a=1.d0
      w1g28a=1.d0
      w0g31=1.d0
      w1g31=1.d0
      w0g32=1.d0
      w1g32=-1.d0

      w0plp25=w0g25a
      w0plp26=w0g26a
      w0plp27=w0g27
      w1plp27=w1g27
      w0plp28=w0g28a
      w0plp29=w0g29
      w0plp30=w0g30
      w0plp31=w0g31
      w1plp31=w1g31
      w0plp32=w0g32
      w1plp32=w1g32

      end

      subroutine external_space_plpmode_value_td()
#include "drt_h.fh"
#include "paraconstants_h.fh"
#include "lpextmode_h.fh"
!     td   lpmode_value
      w0g29=0.d0
      w0g30=0.d0
      w0g26a=0.d0
      w1g26a=0.d0
      w0g25=1.d0
      w0g25a=1.d0
      w1g25a=1.d0
      w0g27=-1.d0
      w1g27=-1.d0
      w0g28a=-1.d0
      w1g28a=-1.d0
      w0g31=1.d0
      w1g31=-1.d0
      w0g32=-1.d0
      w1g32=-1.d0

      w0plp25=w0g25a
      w0plp26=w0g26a
      w0plp27=w0g27
      w1plp27=w1g27
      w0plp28=w0g28a
      w0plp31=w0g31
      w1plp31=w1g31
      w0plp32=w0g32
      w1plp32=w1g32

      end

      subroutine external_space_plpmode_value_ds()
#include "drt_h.fh"
#include "paraconstants_h.fh"
#include "lpextmode_h.fh"
!     ds   lpmode_value
      w0g25=-v_onevsqtwo
      w0g25b=-v_onevsqtwo
      w1g25b=-v_onevsqtwo
      w0g28b=-v_onevsqtwo
      w1g28b=-v_onevsqtwo
      w0g26b=-1.d0
      w1g26b=-1.d0

      w0plp25=w0g25
      w0plp28=w0g28b
      w0plp26=w0g26b

      end

      subroutine external_space_plpmode_value_dt()
#include "drt_h.fh"
#include "paraconstants_h.fh"
#include "lpextmode_h.fh"
!     dt   lpmode_value
      w0g25=v_sqthreevsqtwo   !v_onevsqtwo      !severe_new_error
      w0g25b=v_sqthreevsqtwo
      w1g25b=v_sqthreevsqtwo
      w0g28b=-v_sqthreevsqtwo
      w1g28b=-v_sqthreevsqtwo
      w0g26b=0.d0
      w1g26b=0.d0
      end

      subroutine ss_ext_plpmode(ilnodesm,irnodesm,iltype,irtype,lptype)
#include "drt_h.fh"
#include "lpextmode_h.fh"
      ilrsm=mul_tab(ilnodesm,irnodesm)
      iii=1   !index to determine lwei rwei iposint and nlinkorb
!G2G4a G2G4b G1415 G13
      logic_g36a=.false.
      logic_g36b=.false.
      logic_g35a=.false.
      logic_g35b=.false.
      logic_g34a=.false.
      logic_g34b=.false.
!     if ( ilnodesm.eq.1 .and. irnodesm.eq.3 ) then
!        write(6,*)
!     endif

!G36b
      lpsta36a=iii
      call do_g36mode(ilrsm,ilnodesm,iii)
      lpend36a=iii-4
      if ( lpend36a .ge. lpsta36a ) logic_g36a=.true.   !severe_new_erro
      lpsta35a=iii
      call do_g35mode(ilrsm,ilnodesm,iii)
      lpend35a=iii-4
      if ( lpend35a .ge. lpsta35a ) logic_g35a=.true.
      lpsta34a=iii
      call do_g34mode(ilrsm,ilnodesm,iii)
      lpend34a=iii-4
      if ( lpend34a .ge. lpsta34a ) logic_g34a=.true.
      if ( ilrsm .ne. 1 ) then
         lpsta36b=iii
         call do_g36mode(ilrsm,irnodesm,iii)
         lpend36b=iii-4
         if ( lpend36b .ge. lpsta36b ) logic_g36b=.true.
         lpsta35b=iii
         call do_g35mode(ilrsm,irnodesm,iii)
         lpend35b=iii-4
         if ( lpend35b .ge. lpsta35b ) logic_g35b=.true.
         lpsta34b=iii
         call do_g34mode(ilrsm,irnodesm,iii)
         lpend34b=iii-4
         if ( lpend34b .ge. lpsta34b ) logic_g34b=.true.
      else
         logic_g36b=logic_g36a
         lpsta36b=lpsta36a
         lpend36b=lpend36a
         logic_g35b=logic_g35a
         lpsta35b=lpsta35a
         lpend35b=lpend35a
         logic_g34b=logic_g34a
         lpsta34b=lpsta34a
         lpend34b=lpend34a
      endif

!     G2G4a G2G4b G1415 G13
      logic_g2g4a=.false.
      logic_g2g4b=.false.
      logic_g1415=.false.
      logic_g13=.false.

      if ( irnodesm .eq. 1 .and. irtype .eq. 4 ) then
         logic_g2g4a=.true.
         ism_g2g4=ilnodesm
      endif
      if ( ilnodesm .eq. 1 .and. iltype .eq. 4 ) then
         logic_g2g4b=.true.
         ism_g2g4=irnodesm
      endif

      if ( lptype .ne. 2 .and. ilrsm .eq. 1 ) then
         logic_g1415=.true.
         ism_g1415=ilnodesm
         if ( ilnodesm .eq. 1 .and. iltype .eq. 4
     *      .and.irtype .eq. 4 ) then
            logic_g13=.true.
         endif
         if ( iltype.eq.4 ) then          !severe_error
            idownwei_g131415=irsegdownwei
         else
            idownwei_g131415=ilsegdownwei
         endif
      endif

      nvalue_space_ss=ip2_intsymspace/3

      end

      subroutine do_g36mode(ilrsm,ilnodesm,iii)
#include "drt_h.fh"
      do ismb=1,ng_sm
         isma=mul_tab(ismb,ilrsm)
         if (isma .gt. ismb ) cycle
         ismlink=mul_tab(isma,ilnodesm)
         if ( ismlink .gt. isma ) cycle
         call g36_form(isma,ismb,ismlink,iii)
      enddo
      end
      subroutine do_g34mode(ilrsm,ilnodesm,iii)
#include "drt_h.fh"
      do ismb=1,ng_sm
         isma=mul_tab(ismb,ilrsm)
         if (isma .gt. ismb ) cycle
         ismlink=mul_tab(isma,ilnodesm)
         if ( ismlink .lt. ismb ) cycle
         call g34_form(isma,ismb,ismlink,iii)
      enddo
      end
      subroutine do_g35mode(ilrsm,ilnodesm,iii)
#include "drt_h.fh"
      do ismb=1,ng_sm
         isma=mul_tab(ismb,ilrsm)
         if (isma .gt. ismb ) cycle
         ismlink=mul_tab(isma,ilnodesm)
         if ( ismlink .gt. ismb .or. ismlink .lt. isma ) cycle
         call g35_form(isma,ismb,ismlink,iii)
      enddo
      end

      subroutine g36_form(isma,ismb,ismlink,iii)
#include "drt_h.fh"
#include "lpextmode_h.fh"
      ibsta=ibsm_ext(ismb)
      ibend=iesm_ext(ismb)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)
      ilinksta=ibsm_ext(ismlink)
      ilinkend=iesm_ext(ismlink)
      if ( isma .eq. ismlink ) iasta=iasta+1
      if ( ismb .eq. isma )    ibsta=ibsta+1
      if ( ismb .eq. ismlink ) ibsta=ibsta+1

      do ib=ibsta,ibend
         do ia=iasta,min(iaend,ib-1)
         icsta=ilinksta
         nlinkorb=min(ilinkend,ia-1)-ilinksta+1
         if ( nlinkorb .le. 0 ) cycle
         lpext_wei(iii)=iwt_orb_ext(icsta,ia)
         lpext_wei(iii+1)=iwt_orb_ext(icsta,ib)
         lpext_wei(iii+2)=iwt_orb_ext(ia,ib)
         lpext_wei(iii+3)=nlinkorb
         iii=iii+4
         enddo
      enddo
      end

      subroutine g34_form(isma,ismb,ismlink,iii)
#include "drt_h.fh"
#include "lpextmode_h.fh"
      ibsta=ibsm_ext(ismb)
      ibend=iesm_ext(ismb)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)
      ilinksta=ibsm_ext(ismlink)
      ilinkend=iesm_ext(ismlink)
      if ( ismlink .eq. ismb ) ilinksta=ilinksta+1
      if ( ismlink .eq. isma ) ilinksta=ilinksta+1
      if ( ismb .eq. isma ) ibsta=ibsta+1

      do ilink=ilinksta,ilinkend
         do ib=ibsta,min(ibend,ilink-1)
         naorb=min(iaend,ib-1)-iasta+1
         if ( naorb .le. 0 ) cycle
         lpext_wei(iii)=iwt_orb_ext(iasta,ilink)
         lpext_wei(iii+1)=iwt_orb_ext(ib,ilink)
         lpext_wei(iii+2)=iwt_orb_ext(iasta,ib)
         lpext_wei(iii+3)=naorb
         iii=iii+4
         enddo
      enddo
      end

      subroutine g35_form(isma,ismb,ismlink,iii)
#include "drt_h.fh"
#include "lpextmode_h.fh"
      ibsta=ibsm_ext(ismb)
      ibend=iesm_ext(ismb)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)
      ilinksta=ibsm_ext(ismlink)
      ilinkend=iesm_ext(ismlink)
      if ( ismb .eq. ismlink ) ibsta=ibsta+1
      if ( ismb .eq. isma ) ibsta=ibsta+1
      if ( ismlink .eq. isma ) ilinksta=ilinksta+1

      do ib=ibsta,ibend
         do ilink=ilinksta,min(ilinkend,ib-1)
         naorb=min(iaend,ilink-1)-iasta+1
         if ( naorb .le. 0 ) cycle
         lpext_wei(iii)=iwt_orb_ext(iasta,ilink)
         lpext_wei(iii+1)=iwt_orb_ext(ilink,ib)
         lpext_wei(iii+2)=iwt_orb_ext(iasta,ib)
         lpext_wei(iii+3)=naorb
         iii=iii+4
         enddo
      enddo
      end
      function ibfunction(ib,istep)
      SELECT CASE (istep)
         CASE (2)
            ibfunction=ib-1
         CASE (3)
            ibfunction=ib+1
         CASE DEFAULT
            ibfunction=ib
      END SELECT
      if ( ibfunction .lt. 0 ) then
         write(6,*) "error"
      endif
      end

      subroutine determine_para_array_for_int1ind()
#include "drt_h.fh"
#include "intsort_h.fh"
      do ismabc=1,ng_sm
        nintcount=0
        do ismc=1,ng_sm
          ismab=mul_tab(ismabc,ismc)
          numc=nlsm_ext(ismc)
          do ismb=1,ismc
            isma=mul_tab(ismab,ismb)
            if ( isma .gt. ismb ) cycle
            numb=nlsm_ext(ismb)
            numa=nlsm_ext(isma)
            numint=0
            if ( isma .eq. ismc ) then         !aaa
C              numint=numa-2+ip2(numb-1)+ip3(numc)
              if(numb.gt.1) then
                numint=numa-2+ngw2(numb-1)+ngw3(numc)
              endif
            elseif ( isma .eq. ismb ) then      !aac
C              numint=(numa-1+ip2(numb))*numc
              if(numb.gt.0) then
                if((numa-1+ngw2(numb)).gt.0) then
                  numint=(numa-1+ngw2(numb))*numc
                endif
              endif
            elseif ( ismb .eq. ismc ) then      !acc
C              numint=(numb-1+ip2(numc))*numa
              if(numc.gt.0) then
                if((numb-1+ngw2(numc)).gt.0) then
                  numint=(numb-1+ngw2(numc))*numa
                endif
              endif
            else                        !abc
              numint=numa*numb*numc
            endif
            if ( numint .le. 0 ) cycle            !severe_new_error
            numint=numint*3
            jpsmabc=isma+jp2(ismb)+jp3(ismc)
            intoset_of_ext_jpsmabc(jpsmabc)=nintcount
            nintcount=nintcount+numint
          enddo
        enddo
        intoset_of_aad_isma(ismabc)=nintcount
        numint=nlsm_ext(ismabc)*norb_ext*2
        nintcount=nintcount+numint
        numint=nlsm_ext(ismabc)
        intnum_of_innext_isma(ismabc)=nintcount+numint
      enddo
      end

      subroutine g_dd_ext_sequence_G(ism)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      COMMON /IAIB/ ICAN_A(MAX_ORB),ICAN_B(MTMP+MAX_ORB)
      icano_nnsta=2
      icnt_base=0
      iasta=ibsm_ext(ism)
      ibsta=iasta+1
      ibend=iesm_ext(ism)
      ilwei=0
      do ib=ibsta,ibend
         lrb=norb_number(ib)
         do ia=iasta,ib-1
            lra=norb_number(ia)
            ilwei=ilwei+1
            index_lpext(ilwei)=0
            index_lpext1(ilwei)=0

            NAC=ICAN_A(LRA)+LRB
            index_lpext2(ilwei)=NAC
            value_lpext2(ilwei)=1.0D+00
         enddo
      enddo
      icano_nnend=ibend-iasta+1
      call complete_ext_loop_G()
      end

      subroutine g_tt_ext_sequence_G(ism)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *      m_jc,m_jd, isegsta,isegupwei,isegdownwei
      icano_nnsta=2
      icnt_base=0
      do ismd=1,ng_sm
            ismc=mul_tab(ism,ismd)
            if ( ismc .gt. ismd ) cycle
            id_sta=ibsm_ext(ismd)
            idsta=id_sta
            idend=iesm_ext(ismd)
            ic_sta=ibsm_ext(ismc)
            icend=iesm_ext(ismc)
            if ( ismd.eq.ismc ) idsta=idsta+1
            do id=idsta,idend
                  m_jd=id-id_sta+1
            do ic=ic_sta,min(icend,id-1)
                  m_jc=ic-ic_sta+1
                  icano_nn=iwt_orb_ext(ic,id)
                  if ( icnt_base+icano_nn-1.gt. max_tmpvalue ) then
                        call complete_ext_loop_G()
                        icnt_base=0
                        icano_nnsta=icano_nn
                  endif
                  icano_nnend=icano_nn
                  do ismb=1,ismd-1
                        isma=mul_tab(ism,ismb)
                        if ( isma .gt. ismb ) cycle
                        if ( ismc .gt. ismb ) then
                          call g12_t_diffsym_G(isma,ismb,ismc,ic,id)
                        elseif ( ismc .gt. isma ) then
                          call g11a_t_diffsym_G(isma,ismb,ismc,ic,id)
                        else
                          call g11b_t_diffsym_G(isma,ismb,ismc,ic,id)
                        endif
                  enddo
                  if ( ism.eq.1 ) then
                        isma=ismd
                        call g1112_t_symaaaa_G(isma,ic,id)
                  else
                        isma=mul_tab(ism,ismd)
                        call g11a11b_t_symaacc_G(isma,ismd,ic,id)
                  endif
                  call g36_t_ext_G(ismc,ic,id)
                  call g5_t_ext_G(ismd,ic,id)
                  if ( ism.eq.1 ) call g9_t_ext_G(ismd,ic,id)
                  icnt_base=icnt_base+icano_nn-1
            enddo
            enddo
      enddo
      call complete_ext_loop_G()
      end

      subroutine g_ss_ext_sequence_G(ism,itype)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      icano_nnsta=2
      icnt_base=0
      do ismd=1,ng_sm
         ismc=mul_tab(ism,ismd)
         if ( ismc .gt. ismd ) cycle
         id_sta=ibsm_ext(ismd)
         idsta=id_sta
         idend=iesm_ext(ismd)
         ic_sta=ibsm_ext(ismc)
         icend=iesm_ext(ismc)
         if ( ismd.eq.ismc ) idsta=idsta+1
         do id=idsta,idend
            m_jd=id-id_sta+1
         do ic=ic_sta,min(icend,id-1)
            m_jc=ic-ic_sta+1
            icano_nn=iwt_orb_ext(ic,id)
            if ( icnt_base+icano_nn-1.gt. max_tmpvalue ) then
               call complete_ext_loop_G()
               icnt_base=0
               icano_nnsta=icano_nn
            endif
            icano_nnend=icano_nn
            do ismb=1,ismd-1
               isma=mul_tab(ism,ismb)
               if ( isma .gt. ismb ) cycle
               if ( ismc .gt. ismb ) then
                  call g12_diffsym_G(isma,ismb,ismc,ic,id)
               elseif ( ismc .gt. isma ) then
                  call g11a_diffsym_G(isma,ismb,ismc,ic,id)
               else
                  call g11b_diffsym_G(isma,ismb,ismc,ic,id)
               endif
            enddo
            if ( ism.eq.1 ) then
               isma=ismd
               call g1112_symaaaa_G(isma,ic,id)
            else
               isma=mul_tab(ism,ismd)
               call g11a11b_symaacc_G(isma,ismd,ic,id)
            endif
            call g10_ext_G(ismc,ic,id)
            call g5_ext_G(ismd,ic,id)
            if ( ism.eq.1 ) then
               call g9_ext_G(ismd,ic,id)
            endif
            icnt_base=icnt_base+icano_nn-1
         enddo
         enddo
      enddo
      if ( ism.eq.1 .and. itype.eq.4 ) then
      do id=1,norb_ext
         icano_nn=id+iwt_sm_s_ext
            if ( icnt_base+icano_nn-1.gt. max_tmpvalue ) then
               call complete_ext_loop_G()
               icnt_base=0
               icano_nnsta=icano_nn
            endif
            icano_nnend=icano_nn
            call ext_lp_ab_s1_G(id)
            icnt_base=icnt_base+icano_nn-1
      enddo
      endif
      call complete_ext_loop_G()
      end

      subroutine g12_diffsym_G(isma,ismb,ismc,ic,id)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
!       write(nf2,*) 'g12_diff'
      lrc=norb_number(ic)
      lrd=norb_number(id)

      if ( isma .lt. ismb ) then
      ibsta=ibsm_ext(ismb)
      ibend=iesm_ext(ismb)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)
      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)

      do ib=ibsta,ibend
         lrb=norb_number(ib)
         do ia=iasta,iaend
            lra=norb_number(ia)
            CALL TRANS_IJKL_INTPOS(lra,lrc,lrb,lrd,NXO)
            index_lpext(ilwei)=NXO
            value_lpext(ilwei)=1.0d+00
            CALL TRANS_IJKL_INTPOS(lra,lrd,lrc,lrb,NXO)
            index_lpext1(ilwei)=NXO
            value_lpext1(ilwei)=1.0d+00
            index_lpext2(ilwei)=0
            ilwei=ilwei+1

         enddo
      enddo
      else

      ibsta=ibsm_ext(ismb)+1
      ibend=iesm_ext(ismb)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)
      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)

      do ib=ibsta,ibend
         lrb=norb_number(ib)
         do ia=iasta,ib-1
            lra=norb_number(ia)
            CALL TRANS_IJKL_INTPOS(lra,lrc,lrb,lrd,NXO)
            index_lpext(ilwei)=NXO
            value_lpext(ilwei)=1.0d+00
            CALL TRANS_IJKL_INTPOS(lra,lrd,lrc,lrb,NXO)
            index_lpext1(ilwei)=NXO
            value_lpext1(ilwei)=1.0d+00
            index_lpext2(ilwei)=0
            ilwei=ilwei+1
         enddo
      enddo
      endif

c Avoid unused argument warnings
      if (.false.) call Unused_integer(ismc)
      end

      subroutine g11a_diffsym_G(isma,ismb,ismc,ic,id)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei

!       write(nf2,*) 'g11a_diff'
      lrc=norb_number(ic)
      lrd=norb_number(id)

      ibsta=ibsm_ext(ismb)
      ibend=iesm_ext(ismb)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)
      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)

      do ib=ibsta,ibend
         lrb=norb_number(ib)
         do ia=iasta,iaend
            lra=norb_number(ia)
            CALL TRANS_IJKL_INTPOS(lra,lrc,lrb,lrd,NXO)
            index_lpext(ilwei)=NXO
            value_lpext(ilwei)=1.0d+00
            CALL TRANS_IJKL_INTPOS(lra,lrd,lrb,lrc,NXO)
            index_lpext1(ilwei)=NXO
            value_lpext1(ilwei)=1.0d+00
            index_lpext2(ilwei)=0
            ilwei=ilwei+1
         enddo
      enddo
c Avoid unused argument warnings
      if (.false.) call Unused_integer(ismc)
      end

      subroutine g11b_diffsym_G(isma,ismb,ismc,ic,id)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
!       write(nf2,*) 'g11b_diff'
      lrc=norb_number(ic)
      lrd=norb_number(id)

      ibsta=ibsm_ext(ismb)
      ibend=iesm_ext(ismb)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)
      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)    !need checking
      do ib=ibsta,ibend
         lrb=norb_number(ib)
         do ia=iasta,iaend
            lra=norb_number(ia)
            CALL TRANS_IJKL_INTPOS(lrc,lra,lrb,lrd,NXO)
            index_lpext(ilwei)=NXO
            value_lpext(ilwei)=1.0d+00
            CALL TRANS_IJKL_INTPOS(lrc,lrb,lra,lrd,NXO)
            index_lpext1(ilwei)=NXO
            value_lpext1(ilwei)=1.0d+00
            index_lpext2(ilwei)=0
            ilwei=ilwei+1
         enddo
      enddo
c Avoid unused argument warnings
      if (.false.) call Unused_integer(ismc)
      end

      subroutine g1112_symaaaa_G(isma,ic,id)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
!       write(nf2,*) 'g1112_symaaaa'
      lrc=norb_number(ic)
      lrd=norb_number(id)

      iasta=ibsm_ext(isma)
      ibsta=iasta+1
      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)

      do ib=ibsta,ic-1
         lrb=norb_number(ib)
         do ia=iasta,ib-1
            lra=norb_number(ia)
            CALL TRANS_IJKL_INTPOS(lra,lrc,lrb,lrd,NXO)
            index_lpext(ilwei)=NXO
            value_lpext(ilwei)=1.0d+00
            CALL TRANS_IJKL_INTPOS(lra,lrd,lrc,lrb,NXO)
            index_lpext1(ilwei)=NXO
            value_lpext1(ilwei)=1.0d+00
            index_lpext2(ilwei)=0
            ilwei=ilwei+1
         enddo
      enddo


      do ib=ic+1,id-1
         lrb=norb_number(ib)
         ilwei=icnt_base+iwt_orb_ext(iasta,ib)
         do ia=iasta,ic-1
            lra=norb_number(ia)
            CALL TRANS_IJKL_INTPOS(lra,lrc,lrb,lrd,NXO)
            index_lpext(ilwei)=NXO
            value_lpext(ilwei)=1.0d+00
            CALL TRANS_IJKL_INTPOS(lra,lrd,lrb,lrc,NXO)
            index_lpext1(ilwei)=NXO
            value_lpext1(ilwei)=1.0d+00
            index_lpext2(ilwei)=0
            ilwei=ilwei+1
         enddo
      enddo

      do ib=ic+2,id-1
         lrb=norb_number(ib)
         ilwei=icnt_base+iwt_orb_ext(ic+1,ib)
         do ia=ic+1,ib-1
            lra=norb_number(ia)
            CALL TRANS_IJKL_INTPOS(lrc,lra,lrb,lrd,NXO)
            index_lpext(ilwei)=NXO
            value_lpext(ilwei)=1.0d+00
            CALL TRANS_IJKL_INTPOS(lrc,lrb,lra,lrd,NXO)
            index_lpext1(ilwei)=NXO
            value_lpext1(ilwei)=1.0d+00
            index_lpext2(ilwei)=0
            ilwei=ilwei+1
         enddo
      enddo
      end

      subroutine g11a11b_symaacc_G(isma,ismc,ic,id)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei

!       write(nf2,*) 'g11a11b_symaacc'
      lrc=norb_number(ic)
      lrd=norb_number(id)

      ibsta=ibsm_ext(ismc)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)

      do ib=ibsta,id-1
         lrb=norb_number(ib)
         ilwei=icnt_base+iwt_orb_ext(iasta,ib)
         do ia=iasta,ic-1
            lra=norb_number(ia)
            CALL TRANS_IJKL_INTPOS(lra,lrc,lrb,lrd,NXO)
            index_lpext(ilwei)=NXO
            value_lpext(ilwei)=1.0d+00
            CALL TRANS_IJKL_INTPOS(lra,lrd,lrb,lrc,NXO)
            index_lpext1(ilwei)=NXO
            value_lpext1(ilwei)=1.0d+00
            index_lpext2(ilwei)=0
            ilwei=ilwei+1
         enddo
      enddo

      do ib=ibsta,id-1
         lrb=norb_number(ib)
         ilwei=icnt_base+iwt_orb_ext(ic+1,ib)
         do ia=ic+1,iaend
            lra=norb_number(ia)
            CALL TRANS_IJKL_INTPOS(lrc,lra,lrb,lrd,NXO)
            index_lpext(ilwei)=NXO
            value_lpext(ilwei)=1.0d+00
            CALL TRANS_IJKL_INTPOS(lrc,lrb,lra,lrd,NXO)
            index_lpext1(ilwei)=NXO
            value_lpext1(ilwei)=1.0d+00
            index_lpext2(ilwei)=0
            ilwei=ilwei+1
         enddo
      enddo
      end

      subroutine g10_ext_G(ismc,ic,id)
#include "drt_h.fh"
#include "intsort_h.fh"
      COMMON /IAIB/ ICAN_A(MAX_ORB),ICAN_B(MTMP+MAX_ORB)
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
!        write(nf2,*) 'g10_ext'
      lrc=norb_number(ic)
      lrd=norb_number(id)
      iasta=ibsm_ext(ismc)
      ilwei=icnt_base+iwt_orb_ext(iasta,id)
      do ia=iasta,ic-1
         lra=norb_number(ia)
         CALL TRANS_IJKL_INTPOS(lra,lrc,lrd,lrd,NXO)
         index_lpext(ilwei)=NXO
         value_lpext(ilwei)=1.0d+00
         CALL TRANS_IJKL_INTPOS(lra,lrd,lrc,lrd,NXO)
         index_lpext1(ilwei)=NXO
         value_lpext1(ilwei)=1.0d+00
         NAC=ICAN_A(LRA)+LRC
         index_lpext2(ilwei)=NAC
         value_lpext2(ilwei)=1.0d+00
         ilwei=ilwei+1
      enddo
      end

      subroutine g5_ext_G(ismd,ic,id)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      COMMON /IAIB/ ICAN_A(MAX_ORB),ICAN_B(MTMP+MAX_ORB)
!        write(nf2,*) 'g5_ext'
      lrc=norb_number(ic)
      lrd=norb_number(id)
      ibsta=ibsm_ext(ismd)
      do ib=max(ic+1,ibsta),id-1
         lrb=norb_number(ib)
         ilwei=icnt_base+iwt_orb_ext(ic,ib)
         CALL TRANS_IJKL_INTPOS(lrb,lrd,lrc,lrc,NXO)
         index_lpext(ilwei)=NXO
         value_lpext(ilwei)=1.0d+00
         CALL TRANS_IJKL_INTPOS(lrb,lrc,lrd,lrc,NXO)
         index_lpext1(ilwei)=NXO
         value_lpext1(ilwei)=1.0d+00
         NAC=ICAN_A(LRB)+LRD
         index_lpext2(ilwei)=NAC
         value_lpext2(ilwei)=1.0d+00
      enddo
      end

      subroutine g9_ext_G(ismd,ic,id)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      COMMON /IAIB/ ICAN_A(MAX_ORB),ICAN_B(MTMP+MAX_ORB)
!         write(nf2,*) 'g9_ext'
      lrc=norb_number(ic)
      lrd=norb_number(id)
      iasta=ibsm_ext(ismd)
      ilwei=icnt_base+iwt_orb_ext(iasta,ic)
      do ia=iasta,ic-1
         lra=norb_number(ia)
         CALL TRANS_IJKL_INTPOS(lra,lrd,lrc,lrc,NXO)
         index_lpext(ilwei)=NXO
         value_lpext(ilwei)=1.0d+00
         CALL TRANS_IJKL_INTPOS(lra,lrc,lrd,lrc,NXO)
         index_lpext1(ilwei)=NXO
         value_lpext1(ilwei)=1.0d+00
         NAC=ICAN_A(LRA)+LRD
         index_lpext2(ilwei)=NAC
         value_lpext2(ilwei)=1.0d+00
         ilwei=ilwei+1
      enddo
      end

      subroutine ext_lp_ab_s1_G(id)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      COMMON /IAIB/ ICAN_A(MAX_ORB),ICAN_B(MTMP+MAX_ORB)
      parameter (   v_sqtwo=1.414213562373095d0 )
!         write(nf2,*) 'ext_lp_ab_s1'
      lrd=norb_number(id)
      ilwei=icnt_base

      do ismb=1,ng_sm
         iasta=ibsm_ext(ismb)
         ibend=iesm_ext(ismb)
         ibsta=iasta+1
         do ib=ibsta,ibend
            lrb=norb_number(ib)
         do ia=iasta,ib-1
            lra=norb_number(ia)
            if ( ib .eq. id .or. ia .eq. id ) then
!              g2   arar+ar(head)ar
               ilwei=ilwei+1
               CALL TRANS_IJKL_INTPOS(lra,lrd,lrb,lrd,NXO)
               index_lpext(ilwei)=NXO
               value_lpext(ilwei)=v_sqtwo
               index_lpext1(ilwei)=0
               NAC=ICAN_A(LRA)+LRB
               index_lpext2(ilwei)=NAC
               value_lpext2(ilwei)=v_sqtwo
            else
!              g6g7g8
               ilwei=ilwei+1
               CALL TRANS_IJKL_INTPOS(lra,lrd,lrb,lrd,NXO)
               index_lpext(ilwei)=NXO
               value_lpext(ilwei)=v_sqtwo
               index_lpext1(ilwei)=0
               index_lpext2(ilwei)=0
            endif
         enddo
         enddo
      enddo
!           g1
C===================================
C  Drr-DRR
C  WL=2.0D0 but not 1.0D0 is based on that the non-diagonal just uses th
C  non-triangle <Ci|H|Cj> which designates that I > J.

      do ic=1,id-1
         lrc=norb_number(ic)
         ilwei=ilwei+1
         CALL TRANS_IJKL_INTPOS(lrc,lrd,lrc,lrd,NXO)
         index_lpext(ilwei)=NXO
         value_lpext(ilwei)=2.0D+00
         index_lpext1(ilwei)=0
         index_lpext2(ilwei)=0
      enddo
      end

      subroutine g12_t_diffsym_G(isma,ismb,ismc,ic,id)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      lrc=norb_number(ic)
      lrd=norb_number(id)

      if ( isma .lt. ismb ) then
      ibsta=ibsm_ext(ismb)
      ibend=iesm_ext(ismb)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)
      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)

      do ib=ibsta,ibend
         lrb=norb_number(ib)
         do ia=iasta,iaend
            lra=norb_number(ia)
            CALL TRANS_IJKL_INTPOS(lra,lrc,lrb,lrd,NXO)
            index_lpext(ilwei)=NXO
            value_lpext(ilwei)=1.0d+00
            CALL TRANS_IJKL_INTPOS(lra,lrd,lrc,lrb,NXO)
            index_lpext1(ilwei)=NXO
            value_lpext1(ilwei)=-1.0d+00
            index_lpext2(ilwei)=0
            ilwei=ilwei+1

         enddo
      enddo
      else

      ibsta=ibsm_ext(ismb)+1
      ibend=iesm_ext(ismb)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)
      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)

      do ib=ibsta,ibend
         lrb=norb_number(ib)
         do ia=iasta,ib-1
            lra=norb_number(ia)
            CALL TRANS_IJKL_INTPOS(lra,lrc,lrb,lrd,NXO)
            index_lpext(ilwei)=NXO
            value_lpext(ilwei)=1.0d+00
            CALL TRANS_IJKL_INTPOS(lra,lrd,lrc,lrb,NXO)
            index_lpext1(ilwei)=NXO
            value_lpext1(ilwei)=-1.0d+00
            index_lpext2(ilwei)=0
            ilwei=ilwei+1
         enddo
      enddo
      endif
c Avoid unused argument warnings
      if (.false.) call Unused_integer(ismc)
      end

      subroutine g11a_t_diffsym_G(isma,ismb,ismc,ic,id)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      lrc=norb_number(ic)
      lrd=norb_number(id)

      ibsta=ibsm_ext(ismb)
      ibend=iesm_ext(ismb)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)
      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)

      do ib=ibsta,ibend
         lrb=norb_number(ib)
         do ia=iasta,iaend
            lra=norb_number(ia)
            CALL TRANS_IJKL_INTPOS(lra,lrc,lrb,lrd,NXO)
            index_lpext(ilwei)=NXO
            value_lpext(ilwei)=1.0d+00
            CALL TRANS_IJKL_INTPOS(lra,lrd,lrb,lrc,NXO)
            index_lpext1(ilwei)=NXO
            value_lpext1(ilwei)=-1.0d+00
            index_lpext2(ilwei)=0
            ilwei=ilwei+1
         enddo
      enddo
c Avoid unused argument warnings
      if (.false.) call Unused_integer(ismc)
      end
      subroutine g11b_t_diffsym_G(isma,ismb,ismc,ic,id)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      lrc=norb_number(ic)
      lrd=norb_number(id)

      ibsta=ibsm_ext(ismb)
      ibend=iesm_ext(ismb)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)
      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)    !need checking
      do ib=ibsta,ibend
         lrb=norb_number(ib)
         do ia=iasta,iaend
            lra=norb_number(ia)
            CALL TRANS_IJKL_INTPOS(lrc,lra,lrb,lrd,NXO)
            index_lpext(ilwei)=NXO
            value_lpext(ilwei)=1.0d+00
            CALL TRANS_IJKL_INTPOS(lrc,lrb,lra,lrd,NXO)
            index_lpext1(ilwei)=NXO
            value_lpext1(ilwei)=-1.0d+00
            index_lpext2(ilwei)=0
            ilwei=ilwei+1
         enddo
      enddo
c Avoid unused argument warnings
      if (.false.) call Unused_integer(ismc)
      end


      subroutine g1112_t_symaaaa_G(isma,ic,id)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      lrc=norb_number(ic)
      lrd=norb_number(id)

      iasta=ibsm_ext(isma)
      ibsta=iasta+1
      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)

      do ib=ibsta,ic-1
         lrb=norb_number(ib)
         do ia=iasta,ib-1
            lra=norb_number(ia)
            CALL TRANS_IJKL_INTPOS(lra,lrc,lrb,lrd,NXO)
            index_lpext(ilwei)=NXO
            value_lpext(ilwei)=1.0d+00
            CALL TRANS_IJKL_INTPOS(lra,lrd,lrc,lrb,NXO)
            index_lpext1(ilwei)=NXO
            value_lpext1(ilwei)=-1.0d+00
            index_lpext2(ilwei)=0
            ilwei=ilwei+1
         enddo
      enddo


      do ib=ic+1,id-1
         lrb=norb_number(ib)
         ilwei=icnt_base+iwt_orb_ext(iasta,ib)
         do ia=iasta,ic-1
            lra=norb_number(ia)
            CALL TRANS_IJKL_INTPOS(lra,lrc,lrb,lrd,NXO)
            index_lpext(ilwei)=NXO
            value_lpext(ilwei)=1.0d+00
            CALL TRANS_IJKL_INTPOS(lra,lrd,lrb,lrc,NXO)
            index_lpext1(ilwei)=NXO
            value_lpext1(ilwei)=-1.0d+00
            index_lpext2(ilwei)=0
            ilwei=ilwei+1
         enddo
      enddo

      do ib=ic+2,id-1
         lrb=norb_number(ib)
         ilwei=icnt_base+iwt_orb_ext(ic+1,ib)
         do ia=ic+1,ib-1
            lra=norb_number(ia)
            CALL TRANS_IJKL_INTPOS(lrc,lra,lrb,lrd,NXO)
            index_lpext(ilwei)=NXO
            value_lpext(ilwei)=1.0d+00
            CALL TRANS_IJKL_INTPOS(lrc,lrb,lra,lrd,NXO)
            index_lpext1(ilwei)=NXO
            value_lpext1(ilwei)=-1.0d+00
            index_lpext2(ilwei)=0
            ilwei=ilwei+1
         enddo
      enddo
      end

      subroutine g11a11b_t_symaacc_G(isma,ismc,ic,id)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      lrc=norb_number(ic)
      lrd=norb_number(id)

      ibsta=ibsm_ext(ismc)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)

      do ib=ibsta,id-1
         lrb=norb_number(ib)
         ilwei=icnt_base+iwt_orb_ext(iasta,ib)
         do ia=iasta,ic-1
            lra=norb_number(ia)
            CALL TRANS_IJKL_INTPOS(lra,lrc,lrb,lrd,NXO)
            index_lpext(ilwei)=NXO
            value_lpext(ilwei)=1.0d+00
            CALL TRANS_IJKL_INTPOS(lra,lrd,lrb,lrc,NXO)
            index_lpext1(ilwei)=NXO
            value_lpext1(ilwei)=-1.0d+00
            index_lpext2(ilwei)=0
            ilwei=ilwei+1
         enddo
      enddo

      do ib=ibsta,id-1
         lrb=norb_number(ib)
         ilwei=icnt_base+iwt_orb_ext(ic+1,ib)
         do ia=ic+1,iaend
            lra=norb_number(ia)
            CALL TRANS_IJKL_INTPOS(lrc,lra,lrb,lrd,NXO)
            index_lpext(ilwei)=NXO
            value_lpext(ilwei)=1.0d+00
            CALL TRANS_IJKL_INTPOS(lrc,lrb,lra,lrd,NXO)
            index_lpext1(ilwei)=NXO
            value_lpext1(ilwei)=-1.0d+00
            index_lpext2(ilwei)=0
            ilwei=ilwei+1
         enddo
      enddo
      end

      subroutine g36_t_ext_G(ismc,ic,id)
#include "drt_h.fh"
#include "intsort_h.fh"
      COMMON /IAIB/ ICAN_A(MAX_ORB),ICAN_B(MTMP+MAX_ORB)
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      lrc=norb_number(ic)
      lrd=norb_number(id)
      iasta=ibsm_ext(ismc)
      ilwei=icnt_base+iwt_orb_ext(iasta,id)
      do ia=iasta,ic-1
         lra=norb_number(ia)
         CALL TRANS_IJKL_INTPOS(lra,lrc,lrd,lrd,NXO)
         index_lpext(ilwei)=NXO
         value_lpext(ilwei)=1.0d+00
         CALL TRANS_IJKL_INTPOS(lra,lrd,lrc,lrd,NXO)
         index_lpext1(ilwei)=NXO
         value_lpext1(ilwei)=-1.0d+00
         NAC=ICAN_A(LRA)+LRC
         index_lpext2(ilwei)=NAC
         value_lpext2(ilwei)=1.0d+00
         ilwei=ilwei+1
      enddo
      end

      subroutine g5_t_ext_G(ismd,ic,id)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      COMMON /IAIB/ ICAN_A(MAX_ORB),ICAN_B(MTMP+MAX_ORB)
      lrc=norb_number(ic)
      lrd=norb_number(id)
      ibsta=ibsm_ext(ismd)
      do ib=max(ic+1,ibsta),id-1
         lrb=norb_number(ib)
         ilwei=icnt_base+iwt_orb_ext(ic,ib)
         CALL TRANS_IJKL_INTPOS(lrb,lrd,lrc,lrc,NXO)
         index_lpext(ilwei)=NXO
         value_lpext(ilwei)=1.0d+00
         CALL TRANS_IJKL_INTPOS(lrb,lrc,lrd,lrc,NXO)
         index_lpext1(ilwei)=NXO
         value_lpext1(ilwei)=-1.0d+00
         NAC=ICAN_A(LRB)+LRD
         index_lpext2(ilwei)=NAC
         value_lpext2(ilwei)=1.0d+00
      enddo
      end

      subroutine g9_t_ext_G(ismd,ic,id)
#include "drt_h.fh"
#include "intsort_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      COMMON /IAIB/ ICAN_A(MAX_ORB),ICAN_B(MTMP+MAX_ORB)
      lrc=norb_number(ic)
      lrd=norb_number(id)
      iasta=ibsm_ext(ismd)
      ilwei=icnt_base+iwt_orb_ext(iasta,ic)
      do ia=iasta,ic-1
         lra=norb_number(ia)
         CALL TRANS_IJKL_INTPOS(lra,lrd,lrc,lrc,NXO)
         index_lpext(ilwei)=NXO
         value_lpext(ilwei)=-1.0d+00
         CALL TRANS_IJKL_INTPOS(lra,lrc,lrd,lrc,NXO)
         index_lpext1(ilwei)=NXO
         value_lpext1(ilwei)=1.0d+00
         NAC=ICAN_A(LRA)+LRD
         index_lpext2(ilwei)=NAC
         value_lpext2(ilwei)=-1.0d+00
         ilwei=ilwei+1
      enddo
      end


      subroutine g31_diffsym_G(lri,isma,ismb)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei

      LRC=norb_number(m_jd)

      ibsta=ibsm_ext(ismb)
      ibend=iesm_ext(ismb)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)
      if ( isma .eq. ismb ) ibsta=ibsta+1
      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)
      do ib=ibsta,ibend
         LRB=norb_number(IB)
         if ( isma .eq. ismb ) iaend=ib-1
         do ia=iasta,iaend
            LRA=norb_number(IA)
!     g31   type_12 arbrb^ra^r

         CALL TRANS_IJKL_INTPOS(LRA,LRC,LRB,LRI,NXO)
         index_lpext(ilwei)=NXO
         value_lpext(ilwei)=w0plp31
         CALL TRANS_IJKL_INTPOS(LRA,LRI,LRC,LRB,NXO)
         index_lpext1(ilwei)=NXO
         value_lpext1(ilwei)=w1plp31

         ilwei=ilwei+1
        enddo
      enddo
      end

      subroutine g32a_diffsym_G(lri,isma,ismb)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei

      LRC=norb_number(m_jd)
      ibsta=ibsm_ext(ismb)
      ibend=iesm_ext(ismb)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)
      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)
      do ib=ibsta,ibend
         LRB=norb_number(IB)
        do ia=iasta,iaend
           LRA=norb_number(IA)
           CALL TRANS_IJKL_INTPOS(LRA,LRI,LRB,LRC,NXO)
           index_lpext(ilwei)=NXO
           value_lpext(ilwei)=w0plp32
           CALL TRANS_IJKL_INTPOS(LRA,LRC,LRB,LRI,NXO)
           index_lpext1(ilwei)=NXO
           value_lpext1(ilwei)=-w1plp32

         ilwei=ilwei+1
        enddo
      enddo
      end

      subroutine g32b_diffsym_G(lri,isma,ismb)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei

      LRC=norb_number(m_jd)
      ibsta=ibsm_ext(ismb)
      ibend=iesm_ext(ismb)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)
      if ( ismb.eq.isma ) ibsta=ibsta+1
      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)
      do ib=ibsta,ibend
         LRB=norb_number(IB)
        do ia=iasta,min(iaend,ib-1)
           LRA=norb_number(IA)
      !g32b
           CALL TRANS_IJKL_INTPOS(LRC,LRB,LRA,LRI,NXO)
           index_lpext(ilwei)=NXO
           value_lpext(ilwei)=w0plp32
           CALL TRANS_IJKL_INTPOS(LRC,LRA,LRB,LRI,NXO)
           index_lpext1(ilwei)=NXO
           value_lpext1(ilwei)=-w1plp32

         ilwei=ilwei+1
        enddo
      enddo
      end

      subroutine gsd_samesym_aaa_G(lri,isma)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei

      ic=m_jd
      lrc=norb_number(ic)

      iasta=ibsm_ext(isma)
      ibend=iesm_ext(isma)
      ibsta=iasta+1

      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)
      do ib=ibsta,ic-1
         LRB=norb_number(IB)
        do ia=iasta,ib-1
           LRA=norb_number(IA)
!     g31   type_12 arbrb^ra^r
         CALL TRANS_IJKL_INTPOS(LRA,LRC,LRB,LRI,NXO)
         index_lpext(ilwei)=NXO
         value_lpext(ilwei)=w0plp31
         CALL TRANS_IJKL_INTPOS(LRA,LRI,LRC,LRB,NXO)
         index_lpext1(ilwei)=NXO
         value_lpext1(ilwei)=w1plp31

         ilwei=ilwei+1
      enddo
      enddo

      ib=ic
      ilwei=icnt_base+iwt_orb_ext(iasta,ib)
      do ia=iasta,ib-1
         LRA=norb_number(IA)
!        g28     Cw-Ar                    330
         CALL TRANS_IJKL_INTPOS(LRA,LRC,LRI,LRC,NXO)
         index_lpext(ilwei)=NXO
         value_lpext(ilwei)=w0plp28/w0g28a
         CALL TRANS_IJKL_INTPOS(LRA,LRI,LRC,LRC,NXO)
         index_lpext1(ilwei)=NXO
         value_lpext1(ilwei)=w0plp28

         ilwei=ilwei+1
      enddo

      ia=ic
      do ib=ic+1,ibend
         LRB=norb_number(IB)
!        g25,g27: Bl(20)-Drl(11)         220

         ilwei=icnt_base+iwt_orb_ext(ic,ib)

         CALL TRANS_IJKL_INTPOS(LRB,LRC,LRI,LRC,NXO)
         index_lpext(ilwei)=NXO
         value_lpext(ilwei)=w0plp27
         CALL TRANS_IJKL_INTPOS(LRB,LRI,LRC,LRC,NXO)
         index_lpext1(ilwei)=NXO
         value_lpext1(ilwei)=-w1plp27

      enddo

      do ib=ic+1,ibend
         LRB=norb_number(IB)

        ilwei=icnt_base+iwt_orb_ext(iasta,ib)
        do ia=iasta,ic-1
           LRA=norb_number(IA)

!     g32a   type g12
           CALL TRANS_IJKL_INTPOS(LRA,LRI,LRB,LRC,NXO)
           index_lpext(ilwei)=NXO
           value_lpext(ilwei)=w0plp32
           CALL TRANS_IJKL_INTPOS(LRA,LRC,LRB,LRI,NXO)
           index_lpext1(ilwei)=NXO
           value_lpext1(ilwei)=-w1plp32

          ilwei=ilwei+1
        enddo
      enddo

      do ib=ic+2,ibend
         LRB=norb_number(IB)

        ilwei=icnt_base+iwt_orb_ext(ic+1,ib)
        do ia=ic+1,ib-1
           LRA=norb_number(IA)

!     g32b  type 11
           CALL TRANS_IJKL_INTPOS(LRC,LRB,LRA,LRI,NXO)
           index_lpext(ilwei)=NXO
           value_lpext(ilwei)=w0plp32
           CALL TRANS_IJKL_INTPOS(LRC,LRA,LRB,LRI,NXO)
           index_lpext1(ilwei)=NXO
           value_lpext1(ilwei)=-w1plp32

         ilwei=ilwei+1
        enddo
      enddo
      end

      subroutine gsd_diffsamesym_abb_G(lri,isma,ismb)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei


      ic=m_jd
      lrc=norb_number(ic)

      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)
      ibsta=ibsm_ext(ismb)
      ibend=iesm_ext(ismb)

      ilwei=icnt_base+iwt_orb_ext(iasta,ibsta)
      do ib=ibsta,ic-1
         LRB=norb_number(IB)
        do ia=iasta,iaend
           LRA=norb_number(IA)

!     g31   type_12 arbrb^ra^r
         CALL TRANS_IJKL_INTPOS(LRA,LRC,LRB,LRI,NXO)
         index_lpext(ilwei)=NXO
         value_lpext(ilwei)=w0plp31
         CALL TRANS_IJKL_INTPOS(LRA,LRI,LRC,LRB,NXO)
         index_lpext1(ilwei)=NXO
         value_lpext1(ilwei)=w1plp31

         ilwei=ilwei+1
        enddo
      enddo

      ilwei=icnt_base+iwt_orb_ext(iasta,ic+1)
      do ib=ic+1,ibend
         LRB=norb_number(IB)
        do ia=iasta,iaend
           LRA=norb_number(IA)
         !g32a
           CALL TRANS_IJKL_INTPOS(LRA,LRI,LRB,LRC,NXO)
           index_lpext(ilwei)=NXO
           value_lpext(ilwei)=w0plp32
           CALL TRANS_IJKL_INTPOS(LRA,LRC,LRB,LRI,NXO)
           index_lpext1(ilwei)=NXO
           value_lpext1(ilwei)=-w1plp32

         ilwei=ilwei+1
        enddo
      enddo

      ib=ic
      ilwei=icnt_base+iwt_orb_ext(iasta,ib)
      do ia=iasta,iaend      !ib-1      !severe_error_1020
         LRA=norb_number(IA)

!        g28
         CALL TRANS_IJKL_INTPOS(LRA,LRC,LRI,LRC,NXO)
         index_lpext(ilwei)=NXO
         value_lpext(ilwei)=w0plp28/w0g28a
         CALL TRANS_IJKL_INTPOS(LRA,LRI,LRC,LRC,NXO)
         index_lpext1(ilwei)=NXO
         value_lpext1(ilwei)=w0plp28

         ilwei=ilwei+1
      enddo
      end

      subroutine gsd_diffsamesym_aab_G(lri,isma,ismb)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei

      ic=m_jd
      lrc=norb_number(ic)

      ibsta=ibsm_ext(ismb)
      ibend=iesm_ext(ismb)
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)

      do ib=ibsta,ibend
         LRB=norb_number(IB)
         ilwei=icnt_base+iwt_orb_ext(iasta,ib)
        do ia=iasta,ic-1
         LRA=norb_number(IA)
         !g32a
           CALL TRANS_IJKL_INTPOS(LRA,LRI,LRB,LRC,NXO)
           index_lpext(ilwei)=NXO
           value_lpext(ilwei)=w0plp32
           CALL TRANS_IJKL_INTPOS(LRA,LRC,LRB,LRI,NXO)
           index_lpext1(ilwei)=NXO
           value_lpext1(ilwei)=-w1plp32

         ilwei=ilwei+1
        enddo
      enddo

      do ib=ibsta,ibend
         LRB=norb_number(IB)
         ilwei=icnt_base+iwt_orb_ext(ic+1,ib)
        do ia=ic+1,iaend
           LRA=norb_number(IA)
         !g32b
           CALL TRANS_IJKL_INTPOS(LRC,LRB,LRA,LRI,NXO)
           index_lpext(ilwei)=NXO
           value_lpext(ilwei)=w0plp32
           CALL TRANS_IJKL_INTPOS(LRC,LRA,LRB,LRI,NXO)
           index_lpext1(ilwei)=NXO
           value_lpext1(ilwei)=-w1plp32

         ilwei=ilwei+1
        enddo
      enddo

      ia=ic
      do ib=ibsta,ibend
         LRB=norb_number(IB)
!        g25
         ilwei=icnt_base+iwt_orb_ext(ic,ib)

         CALL TRANS_IJKL_INTPOS(LRB,LRC,LRI,LRC,NXO)
         index_lpext(ilwei)=NXO
         value_lpext(ilwei)=w0plp27
         CALL TRANS_IJKL_INTPOS(LRB,LRI,LRC,LRC,NXO)
         index_lpext1(ilwei)=NXO
         value_lpext1(ilwei)=-w1plp27

      enddo
      end

      subroutine gsd_arlp_s1_G(lri)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei

      ic=m_jd
      lrc=norb_number(ic)

      ilwei=icnt_base+isegdownwei-norb_ext+1
      do is1orb=1,ic-1
!g30 -B^rD^rr
         lrk=norb_number(is1orb)
         CALL TRANS_IJKL_INTPOS(LRC,LRK,LRI,LRK,NXO)
         index_lpext(ilwei)=NXO
         value_lpext(ilwei)=w0plp30
         index_lpext1(ilwei)=0

         ilwei=ilwei+1
      enddo

!g26 -A^r     610
         lrk=norb_number(ic)
         CALL TRANS_IJKL_INTPOS(LRC,LRK,LRI,LRK,NXO)
         index_lpext(ilwei)=NXO
         value_lpext(ilwei)=w0plp26
         index_lpext1(ilwei)=0

      ilwei=ilwei+1

!g29 -Dl^rA^l
      do is1orb=ic+1,norb_ext
         lrk=norb_number(is1orb)
         CALL TRANS_IJKL_INTPOS(LRC,LRK,LRI,LRK,NXO)
         index_lpext(ilwei)=NXO
         value_lpext(ilwei)=w0plp29
         index_lpext1(ilwei)=0

        ilwei=ilwei+1
      enddo
      end

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
! Module_calcu  Completed LOOP multiply MOs and partly complete DM1 and

! SUBROUTINE lp10_arbrbr_ext_calcuvalue_G
! SUBROUTINE lp11_arblbr_ext_calcuvalue_G
! SUBROUTINE lp12_arblbl_ext_calcuvalue_G
! SUBROUTINE lp9_drlbl_ext_calcuvalue_wyb_G
! SUBROUTINE lp_ar_coe_sd_calcuvalue_wyb_G
! SUBROUTINE lp_ar_coe_td_calcuvalue_wyb_G
! SUBROUTINE lp_ar_coe_dv_calcuvalue_wyb_G
! SUBROUTINE lp_arbr_ext_svtv_calcuvalue_wyb_G
! SUBROUTINE lp_drr_ext_svtv_calcuvalue_wyb_G
! SUBROUTINE lp_arbl_ext_ss_calcuvalue_G
! SUBROUTINE lp_arbl_ext_st_calcuvalue_G
! SUBROUTINE lp_arbl_ext_ts_calcuvalue_G
! SUBROUTINE lp_arbl_ext_tt_calcuvalue_G
! SUBROUTINE lp_arbl_ext_dd_calcuvalue_G
! SUBROUTINE lp_drl_ext_SS_calcuvalue_G
! SUBROUTINE lp_drl_ext_ST_calcuvalue_G
! SUBROUTINE lp_drl_ext_TS_calcuvalue_G
! SUBROUTINE lp_drl_ext_TT_calcuvalue_G
! SUBROUTINE lp_drl_ext_dd_calcuvalue_wyb_G
! SUBROUTINE lp678_ext_wyb_calcuvalue_G_1
! SUBROUTINE lp678_ext_wyb_calcuvalue_G_2
! SUBROUTINE lp_drl_sum_SS_calcuvalue_G
! SUBROUTINE lp_drl_sum_TT_calcuvalue_G
! SUBROUTINE lp8_drlbr_sum_calcuvalue_wyb_G
! SUBROUTINE lp9_drlbl_sum_calcuvalue_wyb_G
! SUBROUTINE gsd_ext_sequence_G


      subroutine lp_drl_ext_SS_calcuvalue_G(lri,nlp_value)
#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "intsort_h.fh"

      ivalue=0

!G2G4b
      if ( logic_g2g4a ) then
        w0lp=w0_plp*w0g2a
        w1lp=w1_plp*w1g2a


      do lmb=1,ng_sm
         ibsta=ibsm_ext(lmb)
         ibend=iesm_ext(lmb)
         do irb=ibsta,ibend
            lrb=norb_number(irb)
            do ira=ibsta,irb-1
               lra=norb_number(ira)

         ivalue=ivalue+2
!        Drl -- B^lA^r =4_3
         CALL TRANS_IJKL_INTPOS(lra,lrb,lri,lri,NXO)
         index_lpext(ivalue)=NXO
         value_lpext(ivalue)=-2.0D0*w0lp

         CALL TRANS_IJKL_INTPOS(lra,lri,lrb,lri,NXO)
         index_lpext1(ivalue)=NXO
         value_lpext1(ivalue)=(w0lp-w1lp)

               enddo
            enddo
        enddo
      endif


      w0lp=w0_plp*w0g36a
      w1lp=w1_plp*w1g36a

      do lmb=1,ng_sm
         ibsta=ibsm_ext(lmb)
         ibend=iesm_ext(lmb)
         do irb=ibsta,ibend
            lrb=norb_number(irb)
            do ira=ibsta,irb-1
               lra=norb_number(ira)
               ivalue=ivalue+1

               CALL TRANS_IJKL_INTPOS(lra,lrb,lri,lri,NXO)
               index_lpext(ivalue)=NXO
               value_lpext(ivalue)=-2.0D0*w0lp

               CALL TRANS_IJKL_INTPOS(lra,lri,lrb,lri,NXO)
               index_lpext1(ivalue)=NXO
               value_lpext1(ivalue)=(w0lp-w1lp)

            enddo
         enddo
      enddo
      nlp_value=ivalue
      end


      subroutine lp_drl_ext_ST_calcuvalue_G(lri,nlp_value)
#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "intsort_h.fh"

      ivalue=0
!G1415
      if ( logic_g1415 ) then

      w1lp=w1_plp*w1g14a
C=========================lyb====================
C THE W1LP SHOULD BE MULTIPLED BY TWO,BECAUSE WE JUST USE HALF OF THEM
C IS H*C CALCULATIONS, ???

      w1lp=w1lp*2.0D+00

      do ismb=1,ng_sm
         isma=mul_tab(ismb,ism_g1415)
         if ( isma .gt. ismb ) cycle
         ibsta=ibsm_ext(ismb)
         ibend=iesm_ext(ismb)
         iasta=ibsm_ext(isma)
         iaend=iesm_ext(isma)
         if ( ismb .eq. isma )    ibsta=ibsta+1
         do ib=ibsta,ibend
            lrb=norb_number(ib)
            do ia=iasta,min(iaend,ib-1)
              lra=norb_number(ia)

            ivalue=ivalue+1
            CALL TRANS_IJKL_INTPOS(lra,LRI,lra,LRI,NXO)
            index_lpext(ivalue)=NXO
            value_lpext(ivalue)=w1lp

            CALL TRANS_IJKL_INTPOS(lrb,LRI,lrb,LRI,NXO)

           index_lpext1(ivalue)=NXO
           value_lpext1(ivalue)=-w1lp

            enddo
         enddo
      enddo
      endif

!G2G4b
      if ( logic_g2g4b ) then
        w1lp=w1_plp*w1g4b
      do lmb=1,ng_sm
         ibsta=ibsm_ext(lmb)
         ibend=iesm_ext(lmb)
         do irb=ibsta,ibend
            lrb=norb_number(irb)
            do ira=ibsta,irb-1
               lra=norb_number(ira)

           ivalue=ivalue+1
           CALL TRANS_IJKL_INTPOS(lra,lri,lrb,lri,NXO)
          index_lpext(ivalue)=NXO
          value_lpext(ivalue)=-w1lp

           ivalue=ivalue+1
          index_lpext(ivalue)=NXO
          value_lpext(ivalue)=w1lp

            enddo
         enddo
      enddo
      endif

      w1lp=w1_plp*w1g36a
      do lmb=1,ng_sm
         ibsta=ibsm_ext(lmb)
         ibend=iesm_ext(lmb)
         do irb=ibsta,ibend
            lrb=norb_number(irb)
            do ira=ibsta,irb-1
               lra=norb_number(ira)
              ivalue=ivalue+1
              CALL TRANS_IJKL_INTPOS(lra,lri,lrb,lri,NXO)
              index_lpext(ivalue)=NXO
              value_lpext(ivalue)=-w1lp
            enddo
         enddo
      enddo
      nlp_value=ivalue

      end

      subroutine lp_drl_ext_TT_calcuvalue_G
     :                           (lri,n1415_value,nlp_value)
#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "intsort_h.fh"

      ivalue=0
!G1415
      if ( logic_g1415 ) then

      w014=w0_plp*w0g14a
      w114=w1_plp*w1g14a
      w015=w0_plp*w0g15a
      w115=w1_plp*w1g15a

      w14lp=w014-w114
      w15lp=w015-w115
C=========================lyb====================
C THE W1LP SHOULD BE MULTIPLED BY TWO,BECAUSE WE JUST USE HALF OF THEM
C IS H*C CALCULATIONS, ???
C

      w14lp=w14lp*2.0D+00
      w15lp=w15lp*2.0D+00

      do ismb=1,ng_sm
         isma=mul_tab(ismb,ism_g1415)
         if ( isma .gt. ismb ) cycle
         ibsta=ibsm_ext(ismb)
         ibend=iesm_ext(ismb)
         iasta=ibsm_ext(isma)
         iaend=iesm_ext(isma)
         if ( ismb .eq. isma )    ibsta=ibsta+1
         do ib=ibsta,ibend
            lrb=norb_number(ib)
            do ia=iasta,min(iaend,ib-1)
              lra=norb_number(ia)
              ivalue=ivalue+1
              CALL TRANS_IJKL_INTPOS(lra,LRI,lra,LRI,NXO)
              index_lpext(ivalue)=NXO
              value_lpext(ivalue)=w15lp
              CALL TRANS_IJKL_INTPOS(lrb,LRI,lrb,LRI,NXO)
              index_lpext1(ivalue)=NXO
              value_lpext1(ivalue)=w14lp
            enddo
         enddo
      enddo
      endif

      n1415_value=ivalue

      w0lp=w0_plp*w0g36a
      w1lp=w1_plp*w1g36a

      do lmb=1,ng_sm
         ibsta=ibsm_ext(lmb)
         ibend=iesm_ext(lmb)
         do irb=ibsta,ibend
            lrb=norb_number(irb)
            do ira=ibsta,irb-1
               lra=norb_number(ira)
               ivalue=ivalue+1

               CALL TRANS_IJKL_INTPOS(lra,lrb,lri,lri,NXO)
               index_lpext(ivalue)=NXO
               value_lpext(ivalue)=-2.0D0*w0lp

               CALL TRANS_IJKL_INTPOS(lra,lri,lrb,lri,NXO)
               index_lpext1(ivalue)=NXO
               value_lpext1(ivalue)=(w0lp-w1lp)

            enddo
         enddo
      enddo
      nlp_value=ivalue
      end

      subroutine lp_drl_SUM_TT_calcuvalue_G
     :                           (lri,n1415_value,nlp_value)
#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "intsort_h.fh"

      ivalue=0
!G1415
      if ( logic_g1415 ) then

      w014=w0_plp*w0g14a
      w114=w1_plp*w1g14a
      w015=w0_plp*w0g15a
      w115=w1_plp*w1g15a

      w14lp=w014-w114
      w15lp=w015-w115
C=========================lyb====================
C THE W1LP SHOULD BE MULTIPLED BY TWO,BECAUSE WE JUST USE HALF OF THEM
C IS H*C CALCULATIONS, ???
C
      w14lp=w14lp*2.0D+00
      w15lp=w15lp*2.0D+00
      do ismb=1,ng_sm
         isma=mul_tab(ismb,ism_g1415)
         if ( isma .gt. ismb ) cycle
         ibsta=ibsm_ext(ismb)
         ibend=iesm_ext(ismb)
         iasta=ibsm_ext(isma)
         iaend=iesm_ext(isma)
         if ( ismb .eq. isma )    ibsta=ibsta+1
         do ib=ibsta,ibend
            lrb=norb_number(ib)
            do ia=iasta,min(iaend,ib-1)
              lra=norb_number(ia)
            ivalue=ivalue+1
            CALL TRANS_IJKL_INTPOS(lra,LRI,lra,LRI,NXO)
            index_lpext(ivalue)=NXO
            value_lpext(ivalue)=w15lp
            CALL TRANS_IJKL_INTPOS(lrb,LRI,lrb,LRI,NXO)
           index_lpext1(ivalue)=NXO
           value_lpext1(ivalue)=w14lp
         enddo
         enddo
      enddo
      endif

      n1415_value=ivalue

      w0lp=w0_plp*w0g36a

      do lmb=1,ng_sm
         ibsta=ibsm_ext(lmb)
         ibend=iesm_ext(lmb)
         do irb=ibsta,ibend
            lrb=norb_number(irb)
            do ira=ibsta,irb-1
               lra=norb_number(ira)
               ivalue=ivalue+1
               CALL TRANS_IJKL_INTPOS(lra,lrb,lri,lri,NXO)
               index_lpext(ivalue)=NXO
               value_lpext(ivalue)=-2.0D0*w0lp
               CALL TRANS_IJKL_INTPOS(lra,lri,lrb,lri,NXO)
               index_lpext1(ivalue)=NXO
               value_lpext1(ivalue)=w0lp
            enddo
         enddo
      enddo
      nlp_value=ivalue
      end
      subroutine lp_drl_ext_TS_calcuvalue_G(lri,nlp_value)
#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "intsort_h.fh"

      ivalue=0
!G1415
      if ( logic_g1415 ) then

      w1lp=w1_plp*w1g14a
C=========================lyb====================
C THE W1LP SHOULD BE MULTIPLED BY TWO,BECAUSE WE JUST USE HALF OF THEM
C IS H*C CALCULATIONS, ???

      w1lp=w1lp*2.0D+00
      do ismb=1,ng_sm
         isma=mul_tab(ismb,ism_g1415)
         if ( isma .gt. ismb ) cycle
         ibsta=ibsm_ext(ismb)
         ibend=iesm_ext(ismb)
         iasta=ibsm_ext(isma)
         iaend=iesm_ext(isma)
         if ( ismb .eq. isma )    ibsta=ibsta+1
         do ib=ibsta,ibend
            lrb=norb_number(ib)
            do ia=iasta,min(iaend,ib-1)
              lra=norb_number(ia)
              ivalue=ivalue+1
              CALL TRANS_IJKL_INTPOS(lra,LRI,lra,LRI,NXO)
              index_lpext(ivalue)=NXO
              value_lpext(ivalue)=w1lp
              CALL TRANS_IJKL_INTPOS(lrb,LRI,lrb,LRI,NXO)
              index_lpext1(ivalue)=NXO
              value_lpext1(ivalue)=-w1lp
           enddo
         enddo
      enddo
      endif

!G2G4b
      if ( logic_g2g4a ) then
c       w0lp=w0_plp*w0g2a
        w1lp=w1_plp*w1g2a

C       valuelptmp1=w0lp
C       w0lp=w0lp-w1lp
C       w1lp=-valuelptmp1*2.d0
C       valuelptmp1=ww0lp
C       ww0lp=ww0lp-ww1lp
C       ww1lp=-valuelptmp1*2.d0


      do lmb=1,ng_sm
         ibsta=ibsm_ext(lmb)
         ibend=iesm_ext(lmb)
         do irb=ibsta,ibend
            lrb=norb_number(irb)
            do ira=ibsta,irb-1
               lra=norb_number(ira)

           ivalue=ivalue+1
           CALL TRANS_IJKL_INTPOS(lra,lri,lrb,lri,NXO)
          index_lpext(ivalue)=NXO
          value_lpext(ivalue)=-w1lp

           ivalue=ivalue+1
          index_lpext(ivalue)=NXO
          value_lpext(ivalue)=w1lp
             enddo
           enddo
        enddo
      endif


      w1lp=w1_plp*w1g36a

      do lmb=1,ng_sm
         ibsta=ibsm_ext(lmb)
         ibend=iesm_ext(lmb)
         do irb=ibsta,ibend
            lrb=norb_number(irb)
            do ira=ibsta,irb-1
               lra=norb_number(ira)
               ivalue=ivalue+1
               CALL TRANS_IJKL_INTPOS(lra,lri,lrb,lri,NXO)
               index_lpext(ivalue)=NXO
               value_lpext(ivalue)=-w1lp
            enddo
         enddo
      enddo
      nlp_value=ivalue
      end

      subroutine lp_arbl_ext_st_calcuvalue_G(lri,lrj,nlp_value)
#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "intsort_h.fh"

      ivalue=0
!G1415
      if ( logic_g1415 ) then

      w0lp=w0_plp*w0g14a
      w1lp=w1_plp*w1g14a
      ww0lp=w0_plp*w0g15a
      ww1lp=w1_plp*w1g15a

      valuelptmp1=w0lp
      w0lp=w0lp-w1lp
      w1lp=-valuelptmp1*2.d0
      valuelptmp1=ww0lp
      ww0lp=ww0lp-ww1lp
      ww1lp=-valuelptmp1*2.d0

      do ismb=1,ng_sm
         isma=mul_tab(ismb,ism_g1415)
         if ( isma .gt. ismb ) cycle
         ibsta=ibsm_ext(ismb)
         ibend=iesm_ext(ismb)
         iasta=ibsm_ext(isma)
         iaend=iesm_ext(isma)
         if ( ismb .eq. isma )    ibsta=ibsta+1
         do ib=ibsta,ibend
            lrb=norb_number(ib)
            do ia=iasta,min(iaend,ib-1)
               lra=norb_number(ia)

               ivalue=ivalue+1
               CALL TRANS_IJKL_INTPOS(LRJ,lrb,LRI,lrb,NXO)
               index_lpext(ivalue)=NXO
               value_lpext(ivalue)=w0lp
               CALL TRANS_IJKL_INTPOS(LRJ,LRI,lrb,lrb,NXO)
               index_lpext1(ivalue)=NXO
               value_lpext1(ivalue)=w1lp

               ivalue=ivalue+1
               CALL TRANS_IJKL_INTPOS(LRJ,lra,LRI,lra,NXO)
               index_lpext(ivalue)=NXO
               value_lpext(ivalue)=ww0lp
               CALL TRANS_IJKL_INTPOS(LRJ,LRI,lra,lra,NXO)
               index_lpext1(ivalue)=NXO
               value_lpext1(ivalue)=ww1lp
            enddo
         enddo
      enddo
      endif
!G13
      if ( logic_g13 ) then
         w0lp=w0g13a*w0_plp
         do ia=1,norb_ext
            lra=norb_number(ia)

            ivalue=ivalue+1
            CALL TRANS_IJKL_INTPOS(LRJ,lra,LRI,lra,NXO)
            index_lpext(ivalue)=NXO
           value_lpext(ivalue)=w0lp
            CALL TRANS_IJKL_INTPOS(LRJ,LRI,lra,lra,NXO)
            index_lpext1(ivalue)=NXO
           value_lpext1(ivalue)=-w0lp*2.d0

            ivalue=ivalue+1
            index_lpext(ivalue)=0
            index_lpext1(ivalue)=0

         enddo
      endif

      lsmi=lsm_inn(lri)
      lsmj=lsm_inn(lrj)
      lsmij=mul_tab(lsmi,lsmj)

!G2G4a
      if ( logic_g2g4a ) then

        w0lp=w0_plp*w0g2a
        w1lp=w1_plp*w1g2a
        ww0lp=w0_plp*w0g4a
        ww1lp=w1_plp*w1g4a

        valuelptmp1=w0lp
        w0lp=w0lp-w1lp
        w1lp=-valuelptmp1*2.d0
        valuelptmp1=ww0lp
        ww0lp=ww0lp-ww1lp
        ww1lp=-valuelptmp1*2.d0

      do lsmb=1,ng_sm
        lsma=mul_tab(lsmij,lsmb)
        if ( lsma .gt. lsmb ) cycle
        ibsta=ibsm_ext(lsmb)
        ibend=iesm_ext(lsmb)
        iasta=ibsm_ext(lsma)
        iaend=iesm_ext(lsma)
        if ( lsmb .eq. lsma ) ibsta=ibsta+1
        do ib=ibsta,ibend
          LRB=norb_number(ib)
          do ia=iasta,min(iaend,ib-1)
            LRA=norb_number(ia)
!        ArBl -- B^lA^r =12
            ivalue=ivalue+1
            CALL TRANS_IJKL_INTPOS(lra,lri,lrj,lrb,NXO)
            index_lpext(ivalue)=NXO
            value_lpext(ivalue)=w0lp
            CALL TRANS_IJKL_INTPOS(lra,lrb,lrj,lri,NXO)
            index_lpext1(ivalue)=NXO
            value_lpext1(ivalue)=w1lp
!        ArBl -- B^rA^l =11
            ivalue=ivalue+1
            CALL TRANS_IJKL_INTPOS(lra,lrj,lrb,lri,NXO)
            index_lpext(ivalue)=NXO
            value_lpext(ivalue)=ww0lp
            CALL TRANS_IJKL_INTPOS(lra,lrb,lrj,lri,NXO)
            index_lpext1(ivalue)=NXO
            value_lpext1(ivalue)=ww1lp
                  enddo
              enddo
         enddo
      else

!G2G4b
      if ( logic_g2g4b ) then

        w0lp=w0_plp*w0g2b
        w1lp=w1_plp*w1g2b
        ww0lp=w0_plp*w0g4b
        ww1lp=w1_plp*w1g4b

        valuelptmp1=w0lp
        w0lp=w0lp-w1lp
        w1lp=-valuelptmp1*2.d0
        valuelptmp1=ww0lp
        ww0lp=ww0lp-ww1lp
        ww1lp=-valuelptmp1*2.d0

      do lsmb=1,ng_sm
        lsma=mul_tab(lsmij,lsmb)
        if ( lsma .gt. lsmb ) cycle
        ibsta=ibsm_ext(lsmb)
        ibend=iesm_ext(lsmb)
        iasta=ibsm_ext(lsma)
        iaend=iesm_ext(lsma)
        if ( lsmb .eq. lsma ) ibsta=ibsta+1
        do ib=ibsta,ibend
          LRB=norb_number(ib)
          do ia=iasta,min(iaend,ib-1)
            LRA=norb_number(ia)
!        ArBl -- B^lA^r =12
            ivalue=ivalue+1
            CALL TRANS_IJKL_INTPOS(lra,lri,lrj,lrb,NXO)
            index_lpext(ivalue)=NXO
            value_lpext(ivalue)=ww0lp
            CALL TRANS_IJKL_INTPOS(lra,lrb,lrj,lri,NXO)
            index_lpext1(ivalue)=NXO
            value_lpext1(ivalue)=ww1lp
!        ArBl -- B^rA^l =11
            ivalue=ivalue+1
            CALL TRANS_IJKL_INTPOS(lra,lrj,lrb,lri,NXO)
            index_lpext(ivalue)=NXO
            value_lpext(ivalue)=w0lp
            CALL TRANS_IJKL_INTPOS(lra,lrb,lrj,lri,NXO)
            index_lpext1(ivalue)=NXO
            value_lpext1(ivalue)=w1lp
                  enddo
              enddo
         enddo

      endif

      endif
!G36a

      w0lp=w0_plp*w0g36a
      w1lp=w1_plp*w1g36a
      valuelptmp1=w0lp
      w0lp=w0lp-w1lp
      w1lp=-valuelptmp1*2.d0
      do lsmb=1,ng_sm
        lsma=mul_tab(lsmij,lsmb)
        if ( lsma .gt. lsmb ) cycle
        ibsta=ibsm_ext(lsmb)
        ibend=iesm_ext(lsmb)
        iasta=ibsm_ext(lsma)
        iaend=iesm_ext(lsma)
        if ( lsmb .eq. lsma ) ibsta=ibsta+1
        do ib=ibsta,ibend
          LRB=norb_number(ib)
          do ia=iasta,min(iaend,ib-1)
            LRA=norb_number(ia)
!        ArBl -- B^lA^r =12
            ivalue=ivalue+1
            CALL TRANS_IJKL_INTPOS(lra,lri,lrj,lrb,NXO)
            index_lpext(ivalue)=NXO
            value_lpext(ivalue)=w0lp
            CALL TRANS_IJKL_INTPOS(lra,lrb,lrj,lri,NXO)
            index_lpext1(ivalue)=NXO
            value_lpext1(ivalue)=w1lp
          enddo
        enddo
      enddo
!G36b

      w0lp=w0_plp*w0g36b
      w1lp=w1_plp*w1g36b
      valuelptmp1=w0lp
      w0lp=w0lp-w1lp
      w1lp=-valuelptmp1*2.d0
      do lsmb=1,ng_sm
        lsma=mul_tab(lsmij,lsmb)
        if ( lsma .gt. lsmb ) cycle
        ibsta=ibsm_ext(lsmb)
        ibend=iesm_ext(lsmb)
        iasta=ibsm_ext(lsma)
        iaend=iesm_ext(lsma)
        if ( lsmb .eq. lsma ) ibsta=ibsta+1
        do ib=ibsta,ibend
          LRB=norb_number(ib)
          do ia=iasta,min(iaend,ib-1)
            LRA=norb_number(ia)
            ivalue=ivalue+1
!        ArBl -- B^rA^l =11
            CALL TRANS_IJKL_INTPOS(lra,lrj,lrb,lri,NXO)
            index_lpext(ivalue)=NXO
            value_lpext(ivalue)=w0lp
            CALL TRANS_IJKL_INTPOS(lra,lrb,lrj,lri,NXO)
            index_lpext1(ivalue)=NXO
            value_lpext1(ivalue)=w1lp
          enddo
        enddo
      enddo
      nlp_value=ivalue
      end

      subroutine lp10_arbrbr_ext_calcuvalue_G(intentry,
     :                                      isma,nlp_value)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"

      w0_sdplp25=(w0_sdplp-w1_sdplp)*w0g25
      w1_sdplp25=(w0_sdplp+w1_sdplp)*w0g25

      DO LRITMP=NORB_FRZ+1,norb_inn-2
         lsmi=lsm_inn(LRITMP)
         DO LRJTMP=LRITMP+1,norb_inn-1
            lsmj=lsm_inn(LRJTMP)
            lsmij=mul_tab(lsmi,lsmj)
            DO LRKTMP=LRJTMP+1,norb_inn
              lsmk=lsm_inn(LRKTMP)
               IF(mul_tab(lsmij,lsmk).NE.isma) CYCLE
               IJK=LRITMP-NORB_FRZ+NGW2(LRJTMP-NORB_FRZ)
     :                          +NGW3(LRKTMP-NORB_FRZ)
               IF(INTIND_IJKA(IJK).EQ.intentry) THEN
                  LRI=LRITMP
                  LRJ=LRJTMP
                  LRK=LRKTMP
                  GOTO 100
             ENDIF
           ENDDO
         ENDDO
      ENDDO

100   next_sta=ibsm_ext(isma)-1

      ivalue=0
      do m_ia=1,nlsm_ext(isma)
         ira=m_ia+next_sta
         LRA=norb_number(ira)
         ivalue=ivalue+1
         CALL TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRJ,NXO)
         index_lpext(ivalue)=NXO
         value_lpext(ivalue)=w0_sdplp25
         CALL TRANS_IJKL_INTPOS(LRA,LRJ,LRK,LRI,NXO)
         index_lpext1(ivalue)=NXO
         value_lpext1(ivalue)=w1_sdplp25
      enddo
      nlp_value=ivalue
      end

      subroutine lp11_arblbr_ext_calcuvalue_G(intentry,isma,
     :                                                nlp_value)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"


      w0_sdplp25=(w0_sdplp-w1_sdplp)*w0g25
      w1_sdplp25=-2.d0*w0_sdplp*w0g25

      DO LRITMP=NORB_FRZ+1,norb_inn-2
         lsmi=lsm_inn(LRITMP)
         DO LRJTMP=LRITMP+1,norb_inn-1
            lsmj=lsm_inn(LRJTMP)
            lsmij=mul_tab(lsmi,lsmj)
            DO LRKTMP=LRJTMP+1,norb_inn
              lsmk=lsm_inn(LRKTMP)
               IF(mul_tab(lsmij,lsmk).NE.isma) CYCLE
               IJK=LRITMP-NORB_FRZ+NGW2(LRJTMP-NORB_FRZ)
     :                          +NGW3(LRKTMP-NORB_FRZ)
               IF(INTIND_IJKA(IJK).EQ.intentry) THEN
                  LRI=LRITMP
                  LRJ=LRJTMP
                  LRK=LRKTMP
                  GOTO 100
             ENDIF
           ENDDO
         ENDDO
      ENDDO
100   next_sta=ibsm_ext(isma)-1

      ivalue=0
      do m_ia=1,nlsm_ext(isma)
         ira=m_ia+next_sta
         LRA=norb_number(ira)
         ivalue=ivalue+1
         CALL TRANS_IJKL_INTPOS(LRA,LRJ,LRK,LRI,NXO)
         index_lpext(ivalue)=NXO
         value_lpext(ivalue)=w0_sdplp25
         CALL TRANS_IJKL_INTPOS(LRA,LRK,LRJ,LRI,NXO)
         index_lpext1(ivalue)=NXO
         value_lpext1(ivalue)=w1_sdplp25
      enddo
      nlp_value=ivalue
      end

      subroutine lp12_arblbl_ext_calcuvalue_G(intentry,isma,
     :                                                nlp_value)
#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"

      w0_sdplp25=(w0_sdplp-w1_sdplp)*w0g25
      w1_sdplp25=-2.d0*w0_sdplp*w0g25

      DO LRITMP=NORB_FRZ+1,norb_inn-2
         lsmi=lsm_inn(LRITMP)
         DO LRJTMP=LRITMP+1,norb_inn-1
            lsmj=lsm_inn(LRJTMP)
            lsmij=mul_tab(lsmi,lsmj)
            DO LRKTMP=LRJTMP+1,norb_inn
              lsmk=lsm_inn(LRKTMP)
               IF(mul_tab(lsmij,lsmk).NE.isma) CYCLE
               IJK=LRITMP-NORB_FRZ+NGW2(LRJTMP-NORB_FRZ)
     :                          +NGW3(LRKTMP-NORB_FRZ)
               IF(INTIND_IJKA(IJK).EQ.intentry) THEN
                  LRI=LRITMP
                  LRJ=LRJTMP
                  LRK=LRKTMP
                  GOTO 100
             ENDIF
           ENDDO
         ENDDO
      ENDDO

100   next_sta=ibsm_ext(isma)-1

      ivalue=0
      do m_ia=1,nlsm_ext(isma)
         ira=m_ia+next_sta
         LRA=norb_number(ira)
         ivalue=ivalue+1
         CALL TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRJ,NXO)
         index_lpext(ivalue)=NXO
         value_lpext(ivalue)=w0_sdplp25
         CALL TRANS_IJKL_INTPOS(LRA,LRK,LRJ,LRI,NXO)
         index_lpext1(ivalue)=NXO
         value_lpext1(ivalue)=w1_sdplp25
      enddo
      nlp_value=ivalue
      end

      subroutine lp9_drlbl_ext_calcuvalue_G(lri,lrk,isma)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"

      next_sta=ibsm_ext(isma)-1
      w0_sdplp25=(w0_sdplp-w1_sdplp)*w0g25
      w1_sdplp25=-2.d0*w0_sdplp*w0g25

      ivalue=0

      do m_ia=1,nlsm_ext(isma)
          ira=m_ia+next_sta
          LRA=norb_number(ira)
          ivalue=ivalue+1
          CALL TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
          index_lpext(ivalue)=NXO
          value_lpext(ivalue)=w0_sdplp25
          CALL TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
          index_lpext1(ivalue)=NXO
          value_lpext1(ivalue)=w1_sdplp25
      enddo
      end

      subroutine lp8_drlbr_sum_calcuvalue_G(lri,LRK,isma,nv)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"

      next_sta=ibsm_ext(isma)-1

      w0_sdplp25=w0_sdplp*w0g25
      w1_sdplp25=2*w0_sdplp*w0g25

      ivalue=0
      do m_ia=1,nlsm_ext(isma)
          ira=m_ia+next_sta
          LRA=norb_number(ira)
          ivalue=ivalue+1
          CALL TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
          index_lpext(ivalue)=NXO
          value_lpext(ivalue)=w0_sdplp25
          CALL TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
          index_lpext1(ivalue)=NXO
          value_lpext1(ivalue)=-w1_sdplp25
      enddo
      nv=ivalue
      end

      subroutine lp9_drlbl_sum_calcuvalue_G(lri,LRK,isma,nv)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"

      next_sta=ibsm_ext(isma)-1
      w0_sdplp25=w0_sdplp*w0g25
      w1_sdplp25=-2.d0*w0_sdplp*w0g25

      ivalue=0
      do m_ia=1,nlsm_ext(isma)
          ira=m_ia+next_sta
          LRA=norb_number(ira)
          ivalue=ivalue+1
          CALL TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
          index_lpext(ivalue)=NXO
          value_lpext(ivalue)=w0_sdplp25
          CALL TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
          index_lpext1(ivalue)=NXO
          value_lpext1(ivalue)=w1_sdplp25
      enddo
      nv=ivalue
      end


      subroutine gsd_ext_sequence_G(iltype,ilsm,irsm,lri)
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
                  call g31_diffsym_G(lri,isma,ismb)
            enddo

            ismb=ismnoded
            isma=mul_tab(ismnodes,ismb)
            if ( isma .eq. ismb ) then
               call gsd_samesym_aaa_G(lri,isma)
            elseif ( isma .lt. ismb ) then
               call gsd_diffsamesym_abb_G(lri,isma,ismb)
            endif

            do ismb=ismnoded+1,ng_sm
                  isma=mul_tab(ismnodes,ismb)
                  if ( isma .gt. ismb ) cycle
                  if ( ismnoded .gt. isma ) then
                     call g32a_diffsym_G(lri,isma,ismb)
                  elseif ( ismnoded .eq. isma ) then
                    call gsd_diffsamesym_aab_G(lri,isma,ismb)
                  else
                     call g32b_diffsym_G(lri,isma,ismb)
                  endif
            enddo

            if ( ismnodes.eq.1 .and. iltype.eq.4 ) then
                  call gsd_arlp_s1_G(lri)
            endif
            icnt_base=icnt_base+ilnodedownwei
      enddo
      end


      subroutine lp678_ext_calcuvalue_G(lri,lrk,isma,nlp_value)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"

      w0_sdplp25=w0_sdplp*w0g25
      ilpvalue=0
      iasta=ibsm_ext(isma)
      iaend=iesm_ext(isma)
      do ia=iasta,iaend
         LRA=norb_number(ia)
         ilpvalue=ilpvalue+1
         CALL TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
         index_lpext(ilpvalue)=NXO
         value_lpext(ilpvalue)=w0_sdplp25
         index_lpext1(ilpvalue)=0
      enddo
      nlp_value=nlsm_ext(isma)
      end


      subroutine lp_ar_coe_calcuvalue_G
     *          (idtu,isma,lri,lrj,nlp_value,lpcoe,nlp_value1)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      COMMON /IAIB/ ICAN_A(MAX_ORB),ICAN_B(MTMP+MAX_ORB)
      dimension lpcoe(norb_dz+1:norb_inn)

      next_sta=ibsm_ext(isma)-1
      w0_sdplp25=w0_sdplp*w0g25
      ilpvalue=0

      if(idtu.eq.100) then
        lsta=lri
        lend=norb_inn
        ilpvalue1=0
        do m_ia=1,nlsm_ext(isma)
          ira=m_ia+next_sta
          lra=norb_number(ira)
          NIA=ICAN_A(LRA)+LRI
          ilpvalue=ilpvalue+1
          index_lpext5(ilpvalue)=NIA
          value_lpext5(ilpvalue)=w0_sdplp25

          ilpvalue1=0
          do iorb=lsta,lend
            KCOE=lpcoe(iorb)
            CALL NEOC(KCOE,NOCC,TCOE)

            ilpvalue1=ilpvalue1+1
            CALL TRANS_IJKL_INTPOS(LRA,iorb,LRI,iorb,NXO)
            index_lpext3(ilpvalue,ilpvalue1)=NXO
            value_lpext3(ilpvalue,ilpvalue1)=w0_sdplp25*nocc*tcoe
            CALL TRANS_IJKL_INTPOS(LRA,LRI,iorb,iorb,NXO)
            index_lpext4(ilpvalue,ilpvalue1)=NXO
            value_lpext4(ilpvalue,ilpvalue1)=w0_sdplp25*nocc

          enddo
        enddo
        nlp_value=nlsm_ext(isma)
        nlp_value1=ilpvalue1
        return
      endif

      if(idtu.eq.51) then
         lsta=lri
         lend=norb_DZ
         ilpvalue1=0
      do m_ia=1,nlsm_ext(isma)
          ira=m_ia+next_sta
            lra=norb_number(ira)
            NIA=ICAN_A(LRA)+LRI
            ilpvalue=ilpvalue+1
            index_lpext5(ilpvalue)=NIA
            value_lpext5(ilpvalue)=w0_sdplp25

           ilpvalue1=0
            if(LOGIC_DH) then
            ilpvalue1=ilpvalue1+1
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRI,LRI,NXO)
            index_lpext3(ilpvalue,ilpvalue1)=NXO
            value_lpext3(ilpvalue,ilpvalue1)=w0_sdplp25
            index_lpext4(ilpvalue,ilpvalue1)=0

              do idorb=lsta+1,lend
                 ilpvalue1=ilpvalue1+1
                 CALL TRANS_IJKL_INTPOS(LRA,idorb,LRI,idorb,NXO)
                 index_lpext3(ilpvalue,ilpvalue1)=NXO
                 value_lpext3(ilpvalue,ilpvalue1)=-w0_sdplp25
                 CALL TRANS_IJKL_INTPOS(LRA,LRI,idorb,idorb,NXO)
                 index_lpext4(ilpvalue,ilpvalue1)=NXO
                 value_lpext4(ilpvalue,ilpvalue1)=w0_sdplp25*2.D0

              enddo
            endif
            iorbs=max(norb_dz+1,lri)
            do iorb=iorbs,norb_inn
              Kcoe=lpcoe(iorb)
              CALL NEOC(KCOE,NOCC,TCOE)
              ilpvalue1=ilpvalue1+1
              CALL TRANS_IJKL_INTPOS(LRA,iorb,LRI,iorb,NXO)
              index_lpext3(ilpvalue,ilpvalue1)=NXO
              value_lpext3(ilpvalue,ilpvalue1)=w0_sdplp25*nocc*tcoe
              CALL TRANS_IJKL_INTPOS(LRA,LRI,iorb,iorb,NXO)
              index_lpext4(ilpvalue,ilpvalue1)=NXO
              value_lpext4(ilpvalue,ilpvalue1)=w0_sdplp25*nocc
            enddo
        enddo
        nlp_value=nlsm_ext(isma)
        nlp_value1=ilpvalue1
        return
      endif

      if(idtu.eq.25.or.idtu.eq.43) then
      ilpvalue1=0
      do m_ia=1,nlsm_ext(isma)
         ira=m_ia+next_sta
         lra=norb_number(ira)
         NIA=ICAN_A(LRA)+LRI
         ilpvalue=ilpvalue+1
         index_lpext5(ilpvalue)=NIA
         value_lpext5(ilpvalue)=w0_sdplp25

         ilpvalue1=0

         ilpvalue1=ilpvalue1+1
         CALL TRANS_IJKL_INTPOS(LRA,LRI,LRI,LRI,NXO)
         index_lpext3(ilpvalue,ilpvalue1)=NXO
         value_lpext3(ilpvalue,ilpvalue1)=w0_sdplp25
         index_lpext4(ilpvalue,ilpvalue1)=0

        do idorb=lri+1,norb_dz
           ilpvalue1=ilpvalue1+1
           CALL TRANS_IJKL_INTPOS(LRA,idorb,LRI,idorb,NXO)
           index_lpext3(ilpvalue,ilpvalue1)=NXO
           value_lpext3(ilpvalue,ilpvalue1)=-w0_sdplp25
           CALL TRANS_IJKL_INTPOS(LRA,LRI,idorb,idorb,NXO)
           index_lpext4(ilpvalue,ilpvalue1)=NXO
           value_lpext4(ilpvalue,ilpvalue1)=w0_sdplp25*2.D0
        enddo
        iorbs=max(norb_dz+1,lri)
        do iorb=iorbs,norb_inn
          Kcoe=lpcoe(iorb)
          CALL NEOC(KCOE,NOCC,TCOE)
          ilpvalue1=ilpvalue1+1
          CALL TRANS_IJKL_INTPOS(LRA,iorb,LRI,iorb,NXO)
          index_lpext3(ilpvalue,ilpvalue1)=NXO
          value_lpext3(ilpvalue,ilpvalue1)=w0_sdplp25*nocc*tcoe
          CALL TRANS_IJKL_INTPOS(LRA,LRI,iorb,iorb,NXO)
          index_lpext4(ilpvalue,ilpvalue1)=NXO
          value_lpext4(ilpvalue,ilpvalue1)=w0_sdplp25*nocc
        enddo
      enddo
      nlp_value=nlsm_ext(isma)
      nlp_value1=ilpvalue1
      return
      endif

      if(idtu.eq.26) then
         lsta=lri
         lend=norb_DZ
         ilpvalue1=0
        do m_ia=1,nlsm_ext(isma)
          ira=m_ia+next_sta
          lra=norb_number(ira)
          NIA=ICAN_A(LRA)+LRI
          ilpvalue=ilpvalue+1
          index_lpext5(ilpvalue)=NIA
          value_lpext5(ilpvalue)=w0_sdplp25

         ilpvalue1=0

          if(LOGIC_DH) then
            do idorb=lsta+1,lend
               ilpvalue1=ilpvalue1+1
               CALL TRANS_IJKL_INTPOS(LRA,idorb,LRI,idorb,NXO)
               index_lpext3(ilpvalue,ilpvalue1)=NXO
               value_lpext3(ilpvalue,ilpvalue1)=-w0_sdplp25
               CALL TRANS_IJKL_INTPOS(LRA,LRI,idorb,idorb,NXO)
               index_lpext4(ilpvalue,ilpvalue1)=NXO
               value_lpext4(ilpvalue,ilpvalue1)=w0_sdplp25*2.D0
            enddo
          endif
          do iorb=norb_dz+1,norb_inn
              Kcoe=lpcoe(iorb)
              CALL NEOC(KCOE,NOCC,TCOE)
              ilpvalue1=ilpvalue1+1
              CALL TRANS_IJKL_INTPOS(LRA,iorb,LRI,iorb,NXO)
              index_lpext3(ilpvalue,ilpvalue1)=NXO
              value_lpext3(ilpvalue,ilpvalue1)=w0_sdplp25*nocc*tcoe
              CALL TRANS_IJKL_INTPOS(LRA,LRI,iorb,iorb,NXO)
              index_lpext4(ilpvalue,ilpvalue1)=NXO
              value_lpext4(ilpvalue,ilpvalue1)=w0_sdplp25*nocc
          enddo
        enddo
        nlp_value=nlsm_ext(isma)
        nlp_value1=ilpvalue1
      return
      endif
      if(idtu.eq.28.OR.idtu.eq.46) then
      nsorb=1
      ndorb=norb_DZ-lri-1
        icoe=0
        if(ndorb.gt.0) nsorb=1
        if(idtu.eq.28) icoe=-(JB_SYS+2)
        if(idtu.eq.46) icoe=0
        ilpvalue1=0
        do m_ia=1,nlsm_ext(isma)
            ira=m_ia+next_sta
            lra=norb_number(ira)
            NIA=ICAN_A(LRA)+LRI
            ilpvalue=ilpvalue+1
            index_lpext5(ilpvalue)=NIA
            value_lpext5(ilpvalue)=w0_sdplp25

              ilpvalue1=0

            ilpvalue1=ilpvalue1+1
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRI,LRI,NXO)
            index_lpext3(ilpvalue,ilpvalue1)=NXO
            value_lpext3(ilpvalue,ilpvalue1)=w0_sdplp25
            index_lpext4(ilpvalue,ilpvalue1)=0

            if(nsorb.gt.0) then
              do iorb=lri+1,norb_dz
                 if(iorb.eq.lrj) THEN
                    ilpvalue1=ilpvalue1+1
                    CALL TRANS_IJKL_INTPOS(LRA,iorb,LRI,iorb,NXO)
                    index_lpext3(ilpvalue,ilpvalue1)=NXO
                    value_lpext3(ilpvalue,ilpvalue1)=w0_sdplp25*icoe
                    CALL TRANS_IJKL_INTPOS(LRA,LRI,iorb,iorb,NXO)
                    index_lpext4(ilpvalue,ilpvalue1)=NXO
                    value_lpext4(ilpvalue,ilpvalue1)=w0_sdplp25
                ENDIF
                 if(iorb.ne.lrj) THEN
                    ilpvalue1=ilpvalue1+1
                    CALL TRANS_IJKL_INTPOS(LRA,iorb,LRI,iorb,NXO)
                    index_lpext3(ilpvalue,ilpvalue1)=NXO
                    value_lpext3(ilpvalue,ilpvalue1)=-w0_sdplp25
                    CALL TRANS_IJKL_INTPOS(LRA,LRI,iorb,iorb,NXO)
                    index_lpext4(ilpvalue,ilpvalue1)=NXO
                    value_lpext4(ilpvalue,ilpvalue1)=w0_sdplp25*2.D0
                 ENDIF
              enddo
            endif
            do iorb=norb_dz+1,norb_inn
              Kcoe=lpcoe(iorb)
              CALL NEOC(KCOE,NOCC,TCOE)
              ilpvalue1=ilpvalue1+1
              CALL TRANS_IJKL_INTPOS(LRA,iorb,LRI,iorb,NXO)
              index_lpext3(ilpvalue,ilpvalue1)=NXO
              value_lpext3(ilpvalue,ilpvalue1)=w0_sdplp25*nocc*tcoe
              CALL TRANS_IJKL_INTPOS(LRA,LRI,iorb,iorb,NXO)
              index_lpext4(ilpvalue,ilpvalue1)=NXO
              value_lpext4(ilpvalue,ilpvalue1)=w0_sdplp25*nocc
            enddo
        enddo
        nlp_value=nlsm_ext(isma)
        nlp_value1=ilpvalue1
      return
      endif
      if(idtu.eq.57.or.idtu.eq.29) then
      nsorb=1
      ndorb=norb_DZ-lri-1
        if(ndorb.gt.0) nsorb=1
        if(idtu.eq.29) then
          icoe=JB_SYS
        else ! if(idtu.eq.57)
          icoe=-1
        endif
        ilpvalue1=0
        do m_ia=1,nlsm_ext(isma)
            ira=m_ia+next_sta
            lra=norb_number(ira)
            NIA=ICAN_A(LRA)+LRI
            ilpvalue=ilpvalue+1
            index_lpext5(ilpvalue)=NIA
            value_lpext5(ilpvalue)=w0_sdplp25

              ilpvalue1=0

            ilpvalue1=ilpvalue1+1
            CALL TRANS_IJKL_INTPOS(LRA,LRI,LRI,LRI,NXO)
            index_lpext3(ilpvalue,ilpvalue1)=NXO
            value_lpext3(ilpvalue,ilpvalue1)=w0_sdplp25
            index_lpext4(ilpvalue,ilpvalue1)=0

            if(nsorb.gt.0) then
              do iorb=lri+1,norb_dz
                 if(iorb.eq.lrj) THEN
                    ilpvalue1=ilpvalue1+1
                    CALL TRANS_IJKL_INTPOS(LRA,iorb,LRI,iorb,NXO)
                    index_lpext3(ilpvalue,ilpvalue1)=NXO
                    value_lpext3(ilpvalue,ilpvalue1)=w0_sdplp25*icoe
                    CALL TRANS_IJKL_INTPOS(LRA,LRI,iorb,iorb,NXO)
                    index_lpext4(ilpvalue,ilpvalue1)=NXO
                    value_lpext4(ilpvalue,ilpvalue1)=w0_sdplp25
                ENDIF
                 if(iorb.ne.lrj) THEN
                    ilpvalue1=ilpvalue1+1
                    CALL TRANS_IJKL_INTPOS(LRA,iorb,LRI,iorb,NXO)
                    index_lpext3(ilpvalue,ilpvalue1)=NXO
                    value_lpext3(ilpvalue,ilpvalue1)=-w0_sdplp25
                    CALL TRANS_IJKL_INTPOS(LRA,LRI,iorb,iorb,NXO)
                    index_lpext4(ilpvalue,ilpvalue1)=NXO
                    value_lpext4(ilpvalue,ilpvalue1)=w0_sdplp25*2.D0
                 ENDIF
              enddo
            endif

            do iorb=norb_dz+1,norb_inn
              Kcoe=lpcoe(iorb)
              CALL NEOC(KCOE,NOCC,TCOE)
              ilpvalue1=ilpvalue1+1
              CALL TRANS_IJKL_INTPOS(LRA,iorb,LRI,iorb,NXO)
              index_lpext3(ilpvalue,ilpvalue1)=NXO
              value_lpext3(ilpvalue,ilpvalue1)=w0_sdplp25*nocc*tcoe
              CALL TRANS_IJKL_INTPOS(LRA,LRI,iorb,iorb,NXO)
              index_lpext4(ilpvalue,ilpvalue1)=NXO
              value_lpext4(ilpvalue,ilpvalue1)=w0_sdplp25*nocc
            enddo
        enddo
        nlp_value=nlsm_ext(isma)
        nlp_value1=ilpvalue1
      return
      endif
!=====================
      if(idtu.eq.55.or.idtu.eq.73) then    !(11)(23)=55 (11)(13)=75
      ilpvalue1=0
      do m_ia=1,nlsm_ext(isma)
         ira=m_ia+next_sta
         lra=norb_number(ira)
         NIA=ICAN_A(LRA)+LRI
         ilpvalue=ilpvalue+1
         index_lpext5(ilpvalue)=NIA
         value_lpext5(ilpvalue)=w0_sdplp25

         ilpvalue1=0

         ilpvalue1=ilpvalue1+1
         CALL TRANS_IJKL_INTPOS(LRA,LRI,LRI,LRI,NXO)
         index_lpext3(ilpvalue,ilpvalue1)=NXO
         value_lpext3(ilpvalue,ilpvalue1)=w0_sdplp25
         index_lpext4(ilpvalue,ilpvalue1)=0
        do idorb=lri+1,norb_dz
           ilpvalue1=ilpvalue1+1
           CALL TRANS_IJKL_INTPOS(LRA,idorb,LRI,idorb,NXO)
           index_lpext3(ilpvalue,ilpvalue1)=NXO
           value_lpext3(ilpvalue,ilpvalue1)=-w0_sdplp25
           CALL TRANS_IJKL_INTPOS(LRA,LRI,idorb,idorb,NXO)
           index_lpext4(ilpvalue,ilpvalue1)=NXO
           value_lpext4(ilpvalue,ilpvalue1)=w0_sdplp25*2.D0
        enddo
        iorbs=max(norb_dz+1,lri)
        do iorb=iorbs,norb_inn
          Kcoe=lpcoe(iorb)
          CALL NEOC(KCOE,NOCC,TCOE)
          ilpvalue1=ilpvalue1+1
          CALL TRANS_IJKL_INTPOS(LRA,iorb,LRI,iorb,NXO)
          index_lpext3(ilpvalue,ilpvalue1)=NXO
          value_lpext3(ilpvalue,ilpvalue1)=w0_sdplp25*nocc*tcoe
          CALL TRANS_IJKL_INTPOS(LRA,LRI,iorb,iorb,NXO)
          index_lpext4(ilpvalue,ilpvalue1)=NXO
          value_lpext4(ilpvalue,ilpvalue1)=w0_sdplp25*nocc
        enddo
      enddo
      nlp_value=nlsm_ext(isma)
      nlp_value1=ilpvalue1
      return
      endif

      end


      subroutine lp_arbl_ext_dd_calcuvalue_G(lri,lrj,
     :                                 iml,imr,nlp_value)
#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "intsort_h.fh"

      nlbf=ibsm_ext(iml)
      nlef=iesm_ext(iml)
      nrbf=ibsm_ext(imr)
      nref=iesm_ext(imr)

      w0lp=w0_plp*w0gdd
      w1lp=w1_plp*w1gdd
      valuetmp1=w0lp
      w0lp=w0lp-w1lp
      w1lp=-valuetmp1*2.d0
      ivalue=0
!G50
      if ( logic_g50 ) then
       if ( logic_g49b) then

         do ira=nlbf,nlef
           lra=norb_number(ira)
         ivalue=ivalue+1
         CALL TRANS_IJKL_INTPOS(lrj,lra,lri,lra,NXO)
         index_lpext(ivalue)=NXO
         value_lpext(ivalue)=w0lp
         CALL TRANS_IJKL_INTPOS(lrj,lri,lra,lra,NXO)
         index_lpext1(ivalue)=NXO
         value_lpext1(ivalue)=w1lp
        enddo
       endif


      ivalue=ivalue+int_dd_drl

       do irb=nlbf,nlef
          lrb=norb_number(irb)
C            lsmb=lsm(irb)
          do ira=nlbf,irb-1
             lra=norb_number(ira)
C             lsma=lsm(ira)
C             lsmba=mul_tab(lsmb,lsma)
             ivalue=ivalue+1
             CALL TRANS_IJKL_INTPOS(lra,lri,lrj,lrb,NXO)
             index_lpext(ivalue)=NXO
             value_lpext(ivalue)=w0lp
             CALL TRANS_IJKL_INTPOS(lra,lrb,lrj,lri,NXO)
             index_lpext1(ivalue)=NXO
             value_lpext1(ivalue)=w1lp
          enddo
       enddo
       do irb=nlbf,nlef
          lrb=norb_number(irb)
C          lsmb=lsm(irb)
          do ira=nlbf,irb-1
             lra=norb_number(ira)
C             lsma=lsm(ira)
C             lsmba=mul_tab(lsmb,lsma)
             ivalue=ivalue+1
             CALL TRANS_IJKL_INTPOS(lra,lrj,lrb,lri,NXO)
             index_lpext(ivalue)=NXO
             value_lpext(ivalue)=w0lp
             CALL TRANS_IJKL_INTPOS(lra,lrb,lrj,lri,NXO)
             index_lpext1(ivalue)=NXO
             value_lpext1(ivalue)=w1lp
          enddo
       enddo

      else
       ivalue=ivalue+int_dd_drl
       if (logic_g49a ) then
!G49a:Bl_Ar  line=12
       do irb=nrbf,nref
          lrb=norb_number(irb)
          do ira=nlbf,nlef
             lra=norb_number(ira)

             ivalue=ivalue+1
             CALL TRANS_IJKL_INTPOS(lra,lri,lrj,lrb,NXO)
             index_lpext(ivalue)=NXO
             value_lpext(ivalue)=w0lp
             CALL TRANS_IJKL_INTPOS(lra,lrb,lrj,lri,NXO)
             index_lpext1(ivalue)=NXO
             value_lpext1(ivalue)=w1lp
          enddo
       enddo
       else
!G49b:Br_Al  line=11
       do irb=nlbf,nlef
          lrb=norb_number(irb)
          do ira=nrbf,nref
             lra=norb_number(ira)

             ivalue=ivalue+1
             CALL TRANS_IJKL_INTPOS(lra,lrj,lrb,lri,NXO)
             index_lpext(ivalue)=NXO
             value_lpext(ivalue)=w0lp
             CALL TRANS_IJKL_INTPOS(lra,lrb,lrj,lri,NXO)
             index_lpext1(ivalue)=NXO
             value_lpext1(ivalue)=w1lp
          enddo
       enddo
       endif
      endif
      nlp_value=ivalue
      end


      subroutine lp_drl_ext_dd_calcuvalue_G(lri,iml,nlp_value)
#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "intsort_h.fh"

      w0lp=w0_plp*w0gdd
      w1lp=w1_plp*w1gdd
      nliml=nlsm_ext(iml)
      iasta=ibsm_ext(iml)
      iaend=iesm_ext(iml)

      ivalue=0
      if(.not.logic_g49b) goto 100
      DO ira=iasta,iaend
         lra=norb_number(ira)
         ivalue=ivalue+1
         CALL TRANS_IJKL_INTPOS(lra,LRI,lra,LRI,NXO)
         index_lpext(ivalue)=NXO
         value_lpext(ivalue)=-w1lp*2.0D0
         index_lpext1(ivalue)=0
      ENDDO
100   mloop=nliml*(nliml-1)/2

      ivalue=ivalue+int_dd_drl
      jvalue=0
      do irb=iasta,iaend
         lrb=norb_number(irb)
         do ira=iasta,irb-1
            lra=norb_number(ira)
            ivalue=ivalue+1
            jvalue=ivalue+mloop
            CALL TRANS_IJKL_INTPOS(lra,LRI,lrb,LRI,NXO)
            index_lpext(ivalue)=NXO
            value_lpext(ivalue)=w0lp-w1lp
            index_lpext(jvalue)=index_lpext(ivalue)
            value_lpext(jvalue)=value_lpext(ivalue)

            CALL TRANS_IJKL_INTPOS(lra,lrb,LRI,LRI,NXO)
            index_lpext1(ivalue)=NXO
            value_lpext1(ivalue)=-2.0D0*w0lp
            index_lpext1(jvalue)=index_lpext1(ivalue)
            value_lpext1(jvalue)=value_lpext1(ivalue)
         enddo
      enddo
      nlp_value=jvalue
      end

      subroutine lp_arbr_ext_svtv_calcuvalue_G(LRI,LRJ,nlp_value)
#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "intsort_h.fh"

      ivalue=0
      lsmi=lsm_inn(lri)
      lsmj=lsm_inn(lrj)
      lsmij=mul_tab(lsmi,lsmj)
!G36a
      w0lp=w0_plp*w0g36a
      w1lp=w1_plp*w1g36a
      valuelptmp1=w0lp
      w0lp=w0lp-w1lp
      w1lp=valuelptmp1+w1lp
!        ArBr -- B^rA^r =10
      do lsmc=1,ng_sm
        lsmd=mul_tab(lsmij,lsmc)
        if ( lsmd .gt. lsmc ) cycle
        icsta=ibsm_ext(lsmc)
        icend=iesm_ext(lsmc)
        idsta=ibsm_ext(lsmd)
        idend=iesm_ext(lsmd)
        if ( lsmc .eq. lsmd ) icsta=icsta+1
        do ic=icsta,icend
          jc=norb_number(ic)
          do id=idsta,min(idend,ic-1)
             jd=norb_number(id)
             ivalue=ivalue+1
             CALL TRANS_IJKL_INTPOS(jd,lri,lrj,jc,NXO)
             index_lpext(ivalue)=NXO
             value_lpext(ivalue)=w0lp
             CALL TRANS_IJKL_INTPOS(jd,lrj,jc,lri,NXO)
             index_lpext1(ivalue)=NXO
             value_lpext1(ivalue)=w1lp
          enddo
        enddo
      enddo
!G36b
!G1415
      if ( logic_g13 ) then

      w0lp=w0g13a*(w0_plp+w1_plp)
      do ic=1,norb_ext
         lrc=norb_number(ic)
         ivalue=ivalue+1
         CALL TRANS_IJKL_INTPOS(lrj,lrc,lri,lrc,NXO)
         index_lpext(ivalue)=NXO
         value_lpext(ivalue)=w0lp
         index_lpext1(ivalue)=0
      enddo
      endif
      nlp_value=ivalue
      end

      subroutine lp_drr_ext_svtv_calcuvalue_G(lri,nlp_value)
#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "intsort_h.fh"

      ivalue=0
!G36a
      w0lp=(w0_plp+w1_plp)*w0g36a
      do lmb=1,ng_sm
         ibsta=ibsm_ext(lmb)
         ibend=iesm_ext(lmb)
         do irb=ibsta,ibend
            lrb=norb_number(irb)
            do ira=ibsta,irb-1
              lra=norb_number(ira)
              ivalue=ivalue+1
              CALL TRANS_IJKL_INTPOS(lra,lri,lrb,lri,NXO)
              index_lpext(ivalue)=NXO
              value_lpext(ivalue)=w0lp
              index_lpext1(ivalue)=0
            enddo
         enddo
      enddo

      if ( logic_g13 ) then
C===================================
C  Drr-DRR
C  w0lp=2.0D0*w0lp but not 1.0D0*w0lp is based on that the non-diagonal
C  just uses the non-triangle <Ci|H|Cj> which designates that I > J.

        w0lp=0.5d0*w0_plp*w0g13a
        w0lp=2.0d0*w0lp
        do LRA=NORB_ALL,NORB_INN+1,-1
           ivalue=ivalue+1
           CALL TRANS_IJKL_INTPOS(LRA,lri,LRA,lri,NXO)
           index_lpext(ivalue)=NXO
           value_lpext(ivalue)=w0lp
           index_lpext1(ivalue)=0
        enddo
      endif
      nlp_value=ivalue
      end

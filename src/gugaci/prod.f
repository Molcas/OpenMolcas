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
* Copyright (C) 2007, Bingbing Suo                                     *
************************************************************************
c 26 feb 2007 -bsuo- revised by suo bing for multi-root calculation
c***********************************************************
      subroutine gdv_sequence_extspace(ilw,irw)
c***********************************************************
c     26 feb 2007 - revised
c
#include "drt_h.fh"
#include "lpextmode_h.fh"
c      write(6,*) '  dv_test ','  vd_test '

      ! mrpt2
      if(log_prod.eq.3) then
        call gdv_sequence_extspace_pt(ilw,irw)
        return
      endif

      do irot=1,mcroot
        irtidx=indx(irot)
        mm=ilw+irtidx
        nn=irw+1+irtidx
        valuelp=vector2(nn)
        valuelptmp1=vector1(nn)
        do iij=1,ilsegdownwei
          vlptmp=value_lpext(iij)
          mm=mm+1
          vector2(mm)=vector2(mm)+valuelptmp1*vlptmp
          valuelp=valuelp+vector1(mm)*vlptmp
          !if(mm.eq.86.and.nn.eq.5) then
          !  write(6,*) "dv_test 1 ",vlptmp,vector2(mm)
          !endif
          !if(mm.eq.5.and.nn.eq.86) then
          !  write(6,*) "dv_test 2 ",valuelp,vector2(nn)
          !endif
        enddo
        vector2(nn)=valuelp
      enddo
      end

      subroutine gdv_sequence_extspace_pt(ilw,irw)  !log_prod=
#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "pl_structure_h.fh"
      mm=ilw
      nn=irw+1
      do iij=1,ilsegdownwei
        mm=mm+1
        wl=value_lpext(iij)
        vector2(mm)=vector2(mm)+vcm(nn)*wl
!        if(jpad.eq.1.and.ipae.eq.2)then
!          write(6,'(a,4i6,3f18.8)') "dv_test_pt 1 ",jpadl,ipael,mm,nn,
!     *      wl,vcm(nn),
!     *      vector2(mm)
!        endif
!        if(jpadl.eq.1.and.ipael.eq.2)then
!          write(6,'(a,4i6,3f18.8)') "dv_test_pt 2 ",jpad,ipae,mm,nn,
!     *      wl,vcm(nn),
!     *      vector2(mm)
!        endif
      enddo
      end

c********************************************************************
      subroutine gtd_sequence_extspace(iplplweiin,iplprweiin)
c********************************************************************
c     26 feb 2007 - revised
c
#include "drt_h.fh"
#include "lpextmode_h.fh"
      common /cont_tmp/icount_ext
      parameter (   v_sqtwo=1.414213562373095d0 )
c     write(6,*) ' td_test _1/2',' dt_test '

      do irot=1,mcroot
        irtidx=indx(irot)
        iplplwei=iplplweiin+irtidx
        iplprwei=iplprweiin+irtidx

        ilpvalue=0
        if ( logic_g25a ) then
          mm=iplplwei+iweista_g25-1
          nn0=iplprwei
          do itmp=1,nint_g25
            ilpvalue=ilpvalue+1
            vlptmp=value_lpext(ilpvalue)
            nn=nn0
            do i=1,nwei_g25
              mm=mm+1
              nn=nn+1
              vector2(mm)=vector2(mm)+vector1(nn)*vlptmp
              vector2(nn)=vector2(nn)+vector1(mm)*vlptmp
            enddo
          enddo
        elseif ( logic_g25b ) then
          mm=iplplwei+iweista_g25-1
          nn0=iplprwei
          ilpvalue=ilpvalue+1
          do itmp=2,nint_g25
            ilpvalue=ilpvalue+1
            vlptmp=value_lpext(ilpvalue)
            nn=nn0
            do i=1,itmp-1
              mm=mm+1
              nn=nn+1
              vector2(mm)=vector2(mm)+vector1(nn)*vlptmp
              vector2(nn)=vector2(nn)+vector1(mm)*vlptmp
            enddo
          enddo
          mm=iplplwei+iweista_g28-1
          nn=iplprwei
          nn=nn+1
          do itmp=2,nwei_g28
            nn=nn+1
            vlptmp=vector2(nn)
            vlptmp1=vector1(nn)
            ilpvalue=0
            do i=1,itmp-1
              ilpvalue=ilpvalue+1
              vlptmp2=-value_lpext(ilpvalue)
              mm=mm+1
              vector2(mm)=vector2(mm)+vlptmp1*vlptmp2
              vlptmp=vlptmp+vector1(mm)*vlptmp2
            enddo
            vector2(nn)=vlptmp
          enddo
        elseif ( logic_g28a ) then
          mm=iplplwei+iweista_g28-1
          nn0=iplprwei
          do nn=nn0+1,nn0+nwei_g28
            vlptmp=vector2(nn)
            vlptmp1=vector1(nn)
            ilpvalue=0
            do i=1,nint_g28
              ilpvalue=ilpvalue+1
              vlptmp2=-value_lpext(ilpvalue)
              mm=mm+1
              vector2(mm)=vector2(mm)+vlptmp1*vlptmp2
              vlptmp=vlptmp+vector1(mm)*vlptmp2
            enddo
            vector2(nn)=vlptmp
          enddo
        endif
      enddo

c...end of gtd_sequence_extspace
      end

c**************************************************************
      subroutine inn_ext_ss_loop_unpack(iplplweiin,iplprweiin)
c**************************************************************
c     26 feb 2007 - revised
c
#include "drt_h.fh"
#include "lpextmode_h.fh"

c     write(6,*) ' ss_test 1/2'
      do irot=1,mcroot
        irtidx=indx(irot)
        iplplwei=iplplweiin+irtidx
        iplprwei=iplprweiin+irtidx

        ii=1
        if ( logic_g1415 ) then
          mm=iplplwei
          nn=iplprwei
          do i=1,idownwei_g131415
            valuelp=value_lpext(ii)
            mm=mm+1
            nn=nn+1
            vector2(mm)=vector2(mm)+vector1(nn)*valuelp
            vector2(nn)=vector2(nn)+vector1(mm)*valuelp
            ii=ii+1
          enddo
        endif

        ii0=ii
        if ( logic_g2g4a ) then
          ii=ii0
          mm0=iplplwei
          nn0=iplprwei+iwt_sm_s_ext
          mm=mm0
          do ismb=1,ng_sm
            isma=mul_tab(ismb,ism_g2g4)
            if ( isma .gt. ismb ) cycle
            ibsta=ibsm_ext(ismb)
            ibend=iesm_ext(ismb)
            iasta=ibsm_ext(isma)
            iaend=iesm_ext(isma)
            if ( ismb .eq. isma ) ibsta=ibsta+1
            nna=nn0+ibsta-1
            do ib=ibsta,ibend
              nna=nna+1
              nnb=nn0+iasta-1
              valuelptmp1=vector1(nna)
              valuelptmp2=vector2(nna)
              do ia=iasta,min(iaend,ib-1)
                valuelp1=value_lpext(ii)
                valuelp2=value_lpext(ii+1)
                mm=mm+1
                nnb=nnb+1
                vector2(mm)=vector2(mm)+valuelptmp1*valuelp1
     *                     +vector1(nnb)*valuelp2
                vector2(nnb)=vector2(nnb)+vector1(mm)*valuelp2
                valuelptmp2=valuelptmp2+vector1(mm)*valuelp1
                ii=ii+2
              enddo
              vector2(nna)=valuelptmp2
            enddo
          enddo
        endif
        if ( logic_g2g4b ) then
          ii=ii0
          mm0=iplprwei
          nn0=iplplwei+iwt_sm_s_ext
          mm=mm0
          do ismb=1,ng_sm
            isma=mul_tab(ismb,ism_g2g4)
            if ( isma .gt. ismb ) cycle
            ibsta=ibsm_ext(ismb)
            ibend=iesm_ext(ismb)
            iasta=ibsm_ext(isma)
            iaend=iesm_ext(isma)
            if ( ismb .eq. isma ) ibsta=ibsta+1
            nna=nn0+ibsta-1
            do ib=ibsta,ibend
              nna=nna+1
              nnb=nn0+iasta-1
              valuelptmp1=vector1(nna)
              valuelptmp2=vector2(nna)
              do ia=iasta,min(iaend,ib-1)
                valuelp2=value_lpext(ii)
                valuelp1=value_lpext(ii+1)
                mm=mm+1
                nnb=nnb+1
                vector2(mm)=vector2(mm)+valuelptmp1*valuelp1
     *                     +vector1(nnb)*valuelp2
                vector2(nnb)=vector2(nnb)+vector1(mm)*valuelp2
                valuelptmp2=valuelptmp2+vector1(mm)*valuelp1
                ii=ii+2
              enddo
              vector2(nna)=valuelptmp2
            enddo
          enddo
        endif
!     iaddii=(ii-ii0)/2
        ii0=ii-1

        do icle=1,2
          if ( icle.eq.1 .and. logic_g36a ) then
            mm0=iplplwei
            nn0=iplprwei
            lpsta36=lpsta36a
            lpend36=lpend36a
          elseif (  icle.eq.2 .and. logic_g36b ) then
            mm0=iplprwei
            nn0=iplplwei
            lpsta36=lpsta36b
            lpend36=lpend36b
          else
            goto 936
          endif
          do iii=lpsta36,lpend36,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            ii=ii0+lpext_wei(iii+2)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            do i=1,lpext_wei(iii+3)
              valuetmp=value_lpext(ii)
              vector2(mm)=vector2(mm)+vector1(nn)*valuetmp
              vector2(nn)=vector2(nn)+vector1(mm)*valuetmp
              mm=mm+1
              nn=nn+1
            enddo
          enddo

936       continue

          if ( icle.eq.1 .and. logic_g35a ) then
            mm0=iplplwei
            nn0=iplprwei
            lpsta35=lpsta35a
            lpend35=lpend35a
          elseif (  icle.eq.2 .and. logic_g35b ) then
            mm0=iplprwei
            nn0=iplplwei
            lpsta35=lpsta35b
            lpend35=lpend35b
          else
            goto 935
          endif

          do iii=lpsta35,lpend35,4
            mm=mm0+lpext_wei(iii)
            nn=nn0+lpext_wei(iii+1)
            iij=ii0+lpext_wei(iii+2)
            valuelp=vector2(nn)
            valuelptmp1=vector1(nn)
            do i=1,lpext_wei(iii+3)
              valuetmp=value_lpext(iij)
              vector2(mm)=vector2(mm)+valuelptmp1*valuetmp
              valuelp=valuelp+vector1(mm)*valuetmp
              iij=iij+1
              mm=mm+1
            enddo
            vector2(nn)=valuelp
          enddo

935       continue
!     cycle

          if ( icle.eq.1 .and. logic_g34a ) then
            mm0=iplplwei
            nn0=iplprwei
            lpsta34=lpsta34a
            lpend34=lpend34a
          elseif (  icle.eq.2 .and. logic_g34b ) then
            mm0=iplprwei
            nn0=iplplwei
            lpsta34=lpsta34b
            lpend34=lpend34b
          else
            goto 934
          endif
          do iii=lpsta34,lpend34,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            iij=ii0+lpext_wei(iii+2)
            valuelp=vector2(nn)
            valuelptmp1=vector1(nn)
            do i=1,lpext_wei(iii+3)
              valuetmp=value_lpext(iij)
              vector2(mm)=vector2(mm)+valuelptmp1*valuetmp
              valuelp=valuelp+vector1(mm)*valuetmp
              iij=iij+1
              mm=mm+1
            enddo
            vector2(nn)=valuelp
          enddo

934       continue

          ii0=ii0+nvalue_space_ss
        enddo

      enddo

c...end of inn_ext_ss_loop_unpack
      end

c********************************************************************
      subroutine inn_ext_ss_drl_loop_unpack(iplplweiin,iplprweiin)
c********************************************************************
c     26 feb 2007 - revised
c
#include "drt_h.fh"
#include "lpextmode_h.fh"

c     write(6,*) ' ss_test 2/2'
      do irot=1,mcroot
        irtidx=indx(irot)
        iplplwei=iplplweiin+irtidx
        iplprwei=iplprweiin+irtidx

        ii=1
        if ( logic_g1415 ) then
          mm=iplplwei
          nn=iplprwei
          do i=1,idownwei_g131415
            valuelp=value_lpext(ii)
            mm=mm+1
            nn=nn+1
            vector2(mm)=vector2(mm)+vector1(nn)*valuelp
            vector2(nn)=vector2(nn)+vector1(mm)*valuelp
            ii=ii+1
          enddo
        endif

        ii0=ii
        if ( logic_g2g4a ) then
          ii=ii0
          mm0=iplplwei
          nn0=iplprwei+iwt_sm_s_ext
          mm=mm0
          do ismb=1,ng_sm
            isma=mul_tab(ismb,ism_g2g4)
            if ( isma .gt. ismb ) cycle
            ibsta=ibsm_ext(ismb)
            ibend=iesm_ext(ismb)
            iasta=ibsm_ext(isma)
            iaend=iesm_ext(isma)
            if ( ismb .eq. isma ) ibsta=ibsta+1
            nna=nn0+ibsta-1
            do ib=ibsta,ibend
              nna=nna+1
              nnb=nn0+iasta-1
              valuelptmp1=vector1(nna)
              valuelptmp2=vector2(nna)
              do ia=iasta,min(iaend,ib-1)
                valuelp1=value_lpext(ii)
                valuelp2=value_lpext(ii+1)
                mm=mm+1
                nnb=nnb+1
                vector2(mm)=vector2(mm)+valuelptmp1*valuelp2
     *                     +vector1(nnb)*valuelp2
                vector2(nnb)=vector2(nnb)+vector1(mm)*valuelp2
                valuelptmp2=valuelptmp2+vector1(mm)*valuelp2
                ii=ii+2
              enddo
              vector2(nna)=valuelptmp2
            enddo
          enddo
        endif

        if ( logic_g2g4b ) then
          ii=ii0
          mm0=iplprwei
          nn0=iplplwei+iwt_sm_s_ext
          mm=mm0
          do ismb=1,ng_sm
            isma=mul_tab(ismb,ism_g2g4)
            if ( isma .gt. ismb ) cycle
            ibsta=ibsm_ext(ismb)
            ibend=iesm_ext(ismb)
            iasta=ibsm_ext(isma)
            iaend=iesm_ext(isma)
            if ( ismb .eq. isma ) ibsta=ibsta+1
            nna=nn0+ibsta-1
            do ib=ibsta,ibend
              nna=nna+1
              nnb=nn0+iasta-1
              valuelptmp1=vector1(nna)
              valuelptmp2=vector2(nna)
              do ia=iasta,min(iaend,ib-1)
                valuelp2=value_lpext(ii)
                valuelp1=value_lpext(ii+1)
                mm=mm+1
                nnb=nnb+1
                vector2(mm)=vector2(mm)+valuelptmp1*valuelp1
     *                     +vector1(nnb)*valuelp1
                vector2(nnb)=vector2(nnb)+vector1(mm)*valuelp1
                valuelptmp2=valuelptmp2+vector1(mm)*valuelp1
                ii=ii+2
              enddo
              vector2(nna)=valuelptmp2
            enddo
          enddo
        endif
!     iaddii=(ii-ii0)/2
        ii0=ii-1

        do icle=1,2
          if ( icle.eq.1 .and. logic_g36a ) then
            mm0=iplplwei
            nn0=iplprwei
            lpsta36=lpsta36a
            lpend36=lpend36a
          elseif (  icle.eq.2 .and. logic_g36b ) then
            mm0=iplprwei
            nn0=iplplwei
            lpsta36=lpsta36b
            lpend36=lpend36b
          else
            goto 936
          endif
          do iii=lpsta36,lpend36,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            ii=ii0+lpext_wei(iii+2)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            do i=1,lpext_wei(iii+3)
              valuetmp=value_lpext(ii)
              vector2(mm)=vector2(mm)+vector1(nn)*valuetmp
              vector2(nn)=vector2(nn)+vector1(mm)*valuetmp
              mm=mm+1
              nn=nn+1
            enddo
          enddo

936       continue

          if ( icle.eq.1 .and. logic_g35a ) then
            mm0=iplplwei
            nn0=iplprwei
            lpsta35=lpsta35a
            lpend35=lpend35a
          elseif (  icle.eq.2 .and. logic_g35b ) then
            mm0=iplprwei
            nn0=iplplwei
            lpsta35=lpsta35b
            lpend35=lpend35b
          else
            goto 935
          endif

          do iii=lpsta35,lpend35,4
            mm=mm0+lpext_wei(iii)
            nn=nn0+lpext_wei(iii+1)
            iij=ii0+lpext_wei(iii+2)
            valuelp=vector2(nn)
            valuelptmp1=vector1(nn)
            do i=1,lpext_wei(iii+3)
              valuetmp=value_lpext(iij)
              vector2(mm)=vector2(mm)+valuelptmp1*valuetmp
              valuelp=valuelp+vector1(mm)*valuetmp
              iij=iij+1
              mm=mm+1
            enddo
            vector2(nn)=valuelp
          enddo

935       continue
!     cycle

          if ( icle.eq.1 .and. logic_g34a ) then
            mm0=iplplwei
            nn0=iplprwei
            lpsta34=lpsta34a
            lpend34=lpend34a
          elseif (  icle.eq.2 .and. logic_g34b ) then
            mm0=iplprwei
            nn0=iplplwei
            lpsta34=lpsta34b
            lpend34=lpend34b
          else
            goto 934
          endif
          do iii=lpsta34,lpend34,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            iij=ii0+lpext_wei(iii+2)
            valuelp=vector2(nn)
            valuelptmp1=vector1(nn)
            do i=1,lpext_wei(iii+3)
              valuetmp=value_lpext(iij)
              vector2(mm)=vector2(mm)+valuelptmp1*valuetmp
              valuelp=valuelp+vector1(mm)*valuetmp
              iij=iij+1
              mm=mm+1
            enddo
            vector2(nn)=valuelp
          enddo

934       continue
        enddo

      enddo

c...end of inn_ext_ss_drl_loop_unpack
      end

c********************************************************************
      subroutine inn_ext_sv_loop_unpack(ilw,irw)   !,ilsegdownwei)
c********************************************************************
c     26 feb 2007 - revised
c
#include "drt_h.fh"
#include "lpextmode_h.fh"
c     write(6,*) 'sv_test,   tv_test '
      if(log_prod.eq.3) then
        call inn_ext_svloop_unpack_pt(ilw,irw)
        return
      endif

      do irot=1,mcroot
        irtidx=indx(irot)
        mm=ilw+irtidx
        nn=irw+1+irtidx
        valuelp=vector2(nn)
        valuelptmp1=vector1(nn)
        do iij=1,ilsegdownwei
           mm=mm+1
           vector2(mm)=vector2(mm)+valuelptmp1*value_lpext(iij)
           valuelp=valuelp+vector1(mm)*value_lpext(iij)
        enddo
        vector2(nn)=valuelp
      enddo

c...end of inn_ext_sv_loop_unpack
      end

      subroutine inn_ext_svloop_unpack_pt(ilw,irw)
#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "pl_structure_h.fh"
      !character*16 loop_type
      !loop_type=' sv_test_pt, tv_test_pt'
      mm=ilw
      nn=irw+1
      do iij=1,ilsegdownwei
        mm=mm+1
        wl=value_lpext(iij)
        vector2(mm)=vector2(mm)+vcm(nn)*wl
!            vector2(mm)=vector2(mm)+vcm0(nn)*wl
!        if(mm.eq.mtest)then
!           write(nf2,'(2i5,3f18.10)')
!     :         1, nn,vcm0(nn),wl,vector2(mm-ndim_h0)
!         endif
      enddo
      end

c*******************************************************************
      subroutine inn_ext_tt_loop_unpack(iplplweiin,iplprweiin)
c*******************************************************************
c     26 feb 2007 - revised
c
#include "drt_h.fh"
#include "lpextmode_h.fh"
c     write(6,*) '  tt_test 1/2'


      do irot=1,mcroot
        irtidx=indx(irot)
        iplplwei=iplplweiin+irtidx
        iplprwei=iplprweiin+irtidx

        ii=1
        if ( logic_g1415 ) then
          mm=iplplwei
          nn=iplprwei
          do i=1,idownwei_g131415
            valuelp=value_lpext(ii)
            mm=mm+1
            nn=nn+1
            vector2(mm)=vector2(mm)+vector1(nn)*valuelp
            vector2(nn)=vector2(nn)+vector1(mm)*valuelp
            ii=ii+1
          enddo
        endif

        ii0=ii-1         !severe_error

        do icle=1,2
          if ( icle.eq.1 .and. logic_g36a ) then
            mm0=iplplwei
            nn0=iplprwei
            lpsta36=lpsta36a
            lpend36=lpend36a
          elseif (  icle.eq.2 .and. logic_g36b ) then
            mm0=iplprwei
            nn0=iplplwei
            lpsta36=lpsta36b
            lpend36=lpend36b
          else
            goto 936
          endif
          do iii=lpsta36,lpend36,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            ii=ii0+lpext_wei(iii+2)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            do i=1,lpext_wei(iii+3)
              valuelp=value_lpext(ii)
              vector2(mm)=vector2(mm)+vector1(nn)*valuelp
              vector2(nn)=vector2(nn)+vector1(mm)*valuelp
              mm=mm+1
              nn=nn+1
            enddo
          enddo
936       continue

          if ( icle.eq.1 .and. logic_g35a ) then
            mm0=iplplwei
            nn0=iplprwei
            lpsta35=lpsta35a
            lpend35=lpend35a
          elseif (  icle.eq.2 .and. logic_g35b ) then
            mm0=iplprwei
            nn0=iplplwei
            lpsta35=lpsta35b
            lpend35=lpend35b
          else
            goto 935
          endif

          do iii=lpsta35,lpend35,4
            mm=mm0+lpext_wei(iii)
            nn=nn0+lpext_wei(iii+1)
            iij=ii0+lpext_wei(iii+2)
            valuelp=vector2(nn)
            valuelptmp1=vector1(nn)
            do i=1,lpext_wei(iii+3)
              valuetmp=-value_lpext(iij)
              vector2(mm)=vector2(mm)+valuelptmp1*valuetmp
              valuelp=valuelp+vector1(mm)*valuetmp
              iij=iij+1
              mm=mm+1
            enddo
            vector2(nn)=valuelp      !severe_error_1202
          enddo
935       continue

          if ( icle.eq.1 .and. logic_g34a ) then
            mm0=iplplwei
            nn0=iplprwei
            lpsta34=lpsta34a
            lpend34=lpend34a
          elseif (  icle.eq.2 .and. logic_g34b ) then
            mm0=iplprwei
            nn0=iplplwei
            lpsta34=lpsta34b
            lpend34=lpend34b
          else
            goto 934
          endif
          do iii=lpsta34,lpend34,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            iij=ii0+lpext_wei(iii+2)
            valuelp=vector2(nn)
            valuelptmp1=vector1(nn)
            do i=1,lpext_wei(iii+3)
              vector2(mm)=vector2(mm)+valuelptmp1*value_lpext(iij)
              valuelp=valuelp+vector1(mm)*value_lpext(iij)
              iij=iij+1
              mm=mm+1
            enddo
            vector2(nn)=valuelp
          enddo
934       continue
          ii0=ii0+nvalue_space_ss
        enddo

      enddo

c...end of inn_ext_tt_loop_unpack
      end

c*******************************************************************
      subroutine inn_ext_ts_loop_unpack(iplplweiin,iplprweiin)
c*******************************************************************
c     26 feb 2007 - revised
c
#include "drt_h.fh"
#include "lpextmode_h.fh"
c     write(6,*) '  ts_test 1/2 '

      do irot=1,mcroot
        irtidx=indx(irot)
        iplplwei=iplplweiin+irtidx
        iplprwei=iplprweiin+irtidx

        ii=1
        if ( logic_g1415 ) then
          mm=iplplwei
          nn=iplprwei
          do i=1,idownwei_g131415
            valuelp=value_lpext(ii)
            mm=mm+1
            nn=nn+1
            vector2(mm)=vector2(mm)+vector1(nn)*valuelp
            vector2(nn)=vector2(nn)+vector1(mm)*valuelp
            ii=ii+1
          enddo
        endif

        ii0=ii
        if ( logic_g2g4a ) then
          ii=ii0
          mm0=iplplwei
          nn0=iplprwei+iwt_sm_s_ext
          mm=mm0            !severe_error
          do ismb=1,ng_sm
            isma=mul_tab(ismb,ism_g2g4)
            if ( isma .gt. ismb ) cycle
            ibsta=ibsm_ext(ismb)
            ibend=iesm_ext(ismb)
            iasta=ibsm_ext(isma)
            iaend=iesm_ext(isma)
            if ( ismb .eq. isma ) ibsta=ibsta+1
            nna=nn0+ibsta-1
            do ib=ibsta,ibend
              nna=nna+1
              nnb=nn0+iasta-1
              valuelptmp1=vector1(nna)
              valuelptmp2=vector2(nna)
              do ia=iasta,min(iaend,ib-1)
                valuelp1=value_lpext(ii)
                valuelp2=value_lpext(ii+1)
                mm=mm+1
                nnb=nnb+1
                vector2(mm)=vector2(mm)+valuelptmp1*valuelp1
     *                     +vector1(nnb)*valuelp2
                vector2(nnb)=vector2(nnb)+vector1(mm)*valuelp2
                valuelptmp2=valuelptmp2+vector1(mm)*valuelp1
                ii=ii+2

              enddo
              vector2(nna)=valuelptmp2
            enddo
          enddo
        endif
        if ( logic_g2g4b ) then
          ii=ii0
          mm0=iplprwei
          nn0=iplplwei+iwt_sm_s_ext
          mm=mm0             !severe_error
          do ismb=1,ng_sm
            isma=mul_tab(ismb,ism_g2g4)
            if ( isma .gt. ismb ) cycle
            ibsta=ibsm_ext(ismb)
            ibend=iesm_ext(ismb)
            iasta=ibsm_ext(isma)
            iaend=iesm_ext(isma)
            if ( ismb .eq. isma ) ibsta=ibsta+1
            nna=nn0+ibsta-1
            do ib=ibsta,ibend
              nna=nna+1
              nnb=nn0+iasta-1
              valuelptmp1=vector1(nna)
              valuelptmp2=vector2(nna)
              do ia=iasta,min(iaend,ib-1)
                valuelp2=value_lpext(ii)
                valuelp1=value_lpext(ii+1)
                mm=mm+1
                nnb=nnb+1
                vector2(mm)=vector2(mm)+valuelptmp1*valuelp1
     *                     +vector1(nnb)*valuelp2
                vector2(nnb)=vector2(nnb)+vector1(mm)*valuelp2
                valuelptmp2=valuelptmp2+vector1(mm)*valuelp1
                ii=ii+2
              enddo
              vector2(nna)=valuelptmp2
            enddo
          enddo
        endif
        ii0=ii-1         !severe_error

        do icle=1,2
          if ( icle.eq.1 .and. logic_g36a ) then
            mm0=iplplwei
            nn0=iplprwei
            lpsta36=lpsta36a
            lpend36=lpend36a
          elseif (  icle.eq.2 .and. logic_g36b ) then
            mm0=iplprwei
            nn0=iplplwei
            lpsta36=lpsta36b
            lpend36=lpend36b
          else
            goto 936
          endif
          do iii=lpsta36,lpend36,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            ii=ii0+lpext_wei(iii+2)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            do i=1,lpext_wei(iii+3)
              vector2(mm)=vector2(mm)+vector1(nn)*value_lpext(ii)
              vector2(nn)=vector2(nn)+vector1(mm)*value_lpext(ii)
              mm=mm+1
              nn=nn+1
            enddo
          enddo
936       continue

          if ( icle.eq.1 .and. logic_g35a ) then
            mm0=iplplwei
            nn0=iplprwei
            lpsta35=lpsta35a
            lpend35=lpend35a
            do iii=lpsta35,lpend35,4
              mm=mm0+lpext_wei(iii)
              nn=nn0+lpext_wei(iii+1)
              iij=ii0+lpext_wei(iii+2)
              valuelp=vector2(nn)
              valuelptmp1=vector1(nn)
              do i=1,lpext_wei(iii+3)
                valuetmp=-value_lpext(iij)
                vector2(mm)=vector2(mm)+valuelptmp1*valuetmp
                valuelp=valuelp+vector1(mm)*valuetmp
                iij=iij+1
                mm=mm+1
              enddo
              vector2(nn)=valuelp
            enddo
          elseif (  icle.eq.2 .and. logic_g35b ) then
            mm0=iplprwei
            nn0=iplplwei
            lpsta35=lpsta35b
            lpend35=lpend35b
            do iii=lpsta35,lpend35,4
              mm=mm0+lpext_wei(iii)
              nn=nn0+lpext_wei(iii+1)
              iij=ii0+lpext_wei(iii+2)
              valuelp=vector2(nn)         !severe_error
              valuelptmp1=vector1(nn)
              do i=1,lpext_wei(iii+3)
                vector2(mm)=vector2(mm)+valuelptmp1*value_lpext(iij)
                valuelp=valuelp+vector1(mm)*value_lpext(iij)
                iij=iij+1
                mm=mm+1
              enddo
              vector2(nn)=valuelp      !severe_error
            enddo
          endif

          if ( icle.eq.1 .and. logic_g34a ) then
            mm0=iplplwei
            nn0=iplprwei
            lpsta34=lpsta34a
            lpend34=lpend34a
          elseif (  icle.eq.2 .and. logic_g34b ) then
            mm0=iplprwei
            nn0=iplplwei
            lpsta34=lpsta34b
            lpend34=lpend34b
          else
            goto 934
          endif
          do iii=lpsta34,lpend34,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            iij=ii0+lpext_wei(iii+2)
            valuelp=vector2(nn)
            valuelptmp1=vector1(nn)
            do i=1,lpext_wei(iii+3)
              valuetmp=-value_lpext(iij)
              vector2(mm)=vector2(mm)+valuelptmp1*valuetmp
              valuelp=valuelp+vector1(mm)*valuetmp
              iij=iij+1
              mm=mm+1
            enddo
            vector2(nn)=valuelp      !severe_error
          enddo
934       continue
          ii0=ii0+nvalue_space_ss
        enddo
      enddo

c...end of inn_ext_ts_loop_unpack
      end

c*********************************************************************
      subroutine inn_ext_ts_drl_loop_unpack(iplplweiin,iplprweiin)
c*********************************************************************
c     26 feb 2007 - revised
c
#include "drt_h.fh"
#include "lpextmode_h.fh"
c     write(6,*) '  ts_test 2/2'

      do irot=1,mcroot
        irtidx=indx(irot)
        iplplwei=iplplweiin+irtidx
        iplprwei=iplprweiin+irtidx

        ii=1
        if ( logic_g1415 ) then
          mm=iplplwei
          nn=iplprwei
          do i=1,idownwei_g131415
            valuelp=value_lpext(ii)
            mm=mm+1
            nn=nn+1
            vector2(mm)=vector2(mm)+vector1(nn)*valuelp
            vector2(nn)=vector2(nn)+vector1(mm)*valuelp
            ii=ii+1
          enddo
        endif

        ii0=ii
        if ( logic_g2g4a ) then
           ii=ii0
           mm0=iplplwei
           nn0=iplprwei+iwt_sm_s_ext
           mm=mm0            !severe_error
           do ismb=1,ng_sm
             isma=mul_tab(ismb,ism_g2g4)
             if ( isma .gt. ismb ) cycle
               ibsta=ibsm_ext(ismb)
               ibend=iesm_ext(ismb)
               iasta=ibsm_ext(isma)
               iaend=iesm_ext(isma)
               if ( ismb .eq. isma ) ibsta=ibsta+1
               nna=nn0+ibsta-1
               do ib=ibsta,ibend
                 nna=nna+1
                 nnb=nn0+iasta-1
                 valuelptmp1=vector1(nna)
                 valuelptmp2=vector2(nna)
                 do ia=iasta,min(iaend,ib-1)
                   valuelp1=value_lpext(ii)
                   valuelp2=value_lpext(ii+1)
                   mm=mm+1
                   nnb=nnb+1
                   vector2(mm)=vector2(mm)+valuelptmp1*valuelp1
     *                        +vector1(nnb)*valuelp2
                   vector2(nnb)=vector2(nnb)+vector1(mm)*valuelp2
                   valuelptmp2=valuelptmp2+vector1(mm)*valuelp1
                   ii=ii+2
                enddo
                vector2(nna)=valuelptmp2
              enddo
            enddo
          endif

          if ( logic_g2g4b ) then
            ii=ii0
            mm0=iplprwei
            nn0=iplplwei+iwt_sm_s_ext
            mm=mm0             !severe_error
            do ismb=1,ng_sm
              isma=mul_tab(ismb,ism_g2g4)
              if ( isma .gt. ismb ) cycle
              ibsta=ibsm_ext(ismb)
              ibend=iesm_ext(ismb)
              iasta=ibsm_ext(isma)
              iaend=iesm_ext(isma)
              if ( ismb .eq. isma ) ibsta=ibsta+1
              nna=nn0+ibsta-1
              do ib=ibsta,ibend
                nna=nna+1
                nnb=nn0+iasta-1
                valuelptmp1=vector1(nna)
                valuelptmp2=vector2(nna)
                do ia=iasta,min(iaend,ib-1)
                  valuelp2=value_lpext(ii)
                  valuelp1=value_lpext(ii+1)
                  mm=mm+1
                  nnb=nnb+1
                  vector2(mm)=vector2(mm)+valuelptmp1*valuelp1
     *                       +vector1(nnb)*valuelp2
                  vector2(nnb)=vector2(nnb)+vector1(mm)*valuelp2
                  valuelptmp2=valuelptmp2+vector1(mm)*valuelp1
                  ii=ii+2
                enddo
                vector2(nna)=valuelptmp2
              enddo
            enddo
          endif

          ii0=ii-1         !severe_error

          do icle=1,2
            if ( icle.eq.1 .and. logic_g36a ) then
            mm0=iplplwei
            nn0=iplprwei
            lpsta36=lpsta36a
            lpend36=lpend36a
          elseif (  icle.eq.2 .and. logic_g36b ) then
            mm0=iplprwei
            nn0=iplplwei
            lpsta36=lpsta36b
            lpend36=lpend36b
          else
            goto 936
          endif
          do iii=lpsta36,lpend36,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            ii=ii0+lpext_wei(iii+2)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            do i=1,lpext_wei(iii+3)
              vector2(mm)=vector2(mm)+vector1(nn)*value_lpext(ii)
              vector2(nn)=vector2(nn)+vector1(mm)*value_lpext(ii)
              mm=mm+1
              nn=nn+1
            enddo
          enddo
936       continue

          if ( icle.eq.1 .and. logic_g35a ) then
            mm0=iplplwei
            nn0=iplprwei
            lpsta35=lpsta35a
            lpend35=lpend35a
            do iii=lpsta35,lpend35,4
              mm=mm0+lpext_wei(iii)
              nn=nn0+lpext_wei(iii+1)
              iij=ii0+lpext_wei(iii+2)
              valuelp=vector2(nn)
              valuelptmp1=vector1(nn)
              do i=1,lpext_wei(iii+3)
                valuetmp=-value_lpext(iij)
                vector2(mm)=vector2(mm)+valuelptmp1*valuetmp
                valuelp=valuelp+vector1(mm)*valuetmp
                iij=iij+1
                mm=mm+1
              enddo
              vector2(nn)=valuelp
            enddo
          elseif (  icle.eq.2 .and. logic_g35b ) then
            mm0=iplprwei
            nn0=iplplwei
            lpsta35=lpsta35b
            lpend35=lpend35b
            do iii=lpsta35,lpend35,4
              mm=mm0+lpext_wei(iii)
              nn=nn0+lpext_wei(iii+1)
              iij=ii0+lpext_wei(iii+2)
              valuelp=vector2(nn)         !severe_error
              valuelptmp1=vector1(nn)
              do i=1,lpext_wei(iii+3)
                vector2(mm)=vector2(mm)+valuelptmp1*value_lpext(iij)
                valuelp=valuelp+vector1(mm)*value_lpext(iij)
                iij=iij+1
                mm=mm+1
              enddo
              vector2(nn)=valuelp      !severe_error
            enddo
          endif

          if ( icle.eq.1 .and. logic_g34a ) then
            mm0=iplplwei
            nn0=iplprwei
            lpsta34=lpsta34a
            lpend34=lpend34a
          elseif (  icle.eq.2 .and. logic_g34b ) then
            mm0=iplprwei
            nn0=iplplwei
            lpsta34=lpsta34b
            lpend34=lpend34b
          else
            goto 934
          endif
          do iii=lpsta34,lpend34,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            iij=ii0+lpext_wei(iii+2)
            valuelp=vector2(nn)
            valuelptmp1=vector1(nn)
            do i=1,lpext_wei(iii+3)
              valuetmp=-value_lpext(iij)
              vector2(mm)=vector2(mm)+valuelptmp1*valuetmp
              valuelp=valuelp+vector1(mm)*valuetmp
              iij=iij+1
              mm=mm+1
            enddo
            vector2(nn)=valuelp      !severe_error
          enddo
934       continue
        enddo

      enddo

c...end of inn_ext_ts_drl_loop_unpack
      end

c*****************************************************************
      subroutine inn_ext_st_loop_unpack(iplplweiin,iplprweiin)
c*****************************************************************
c     26 feb 2007 - revised
c
#include "drt_h.fh"
#include "lpextmode_h.fh"
c     write(6,*) ' st_test 1/2'

      do irot=1,mcroot
        irtidx=indx(irot)
        iplplwei=iplplweiin+irtidx
        iplprwei=iplprweiin+irtidx

        ii=1
        if ( logic_g1415 ) then
          mm=iplplwei
          nn=iplprwei
          do i=1,idownwei_g131415
            valuelp=value_lpext(ii)
            mm=mm+1
            nn=nn+1
            vector2(mm)=vector2(mm)+vector1(nn)*valuelp
            vector2(nn)=vector2(nn)+vector1(mm)*valuelp
            ii=ii+1
          enddo
        endif

c      write(6,*)'st_g2g4a',iplplwei,iplprwei,vector2(137)
        ii0=ii
        if ( logic_g2g4a ) then
          ii=ii0
          mm0=iplplwei
          nn0=iplprwei+iwt_sm_s_ext
          mm=mm0            !severe_error
          do ismb=1,ng_sm
            isma=mul_tab(ismb,ism_g2g4)
            if ( isma .gt. ismb ) cycle
            ibsta=ibsm_ext(ismb)
            ibend=iesm_ext(ismb)
            iasta=ibsm_ext(isma)
            iaend=iesm_ext(isma)
            if ( ismb .eq. isma ) ibsta=ibsta+1
            nna=nn0+ibsta-1
            do ib=ibsta,ibend
              nna=nna+1
              nnb=nn0+iasta-1
              valuelptmp1=vector1(nna)
              valuelptmp2=vector2(nna)
              do ia=iasta,min(iaend,ib-1)
                valuelp1=value_lpext(ii)
                valuelp2=value_lpext(ii+1)
                mm=mm+1
                nnb=nnb+1
                vector2(mm)=vector2(mm)+valuelptmp1*valuelp1
     *                     +vector1(nnb)*valuelp2
                vector2(nnb)=vector2(nnb)+vector1(mm)*valuelp2
                valuelptmp2=valuelptmp2+vector1(mm)*valuelp1
                ii=ii+2
              enddo
              vector2(nna)=valuelptmp2
            enddo
          enddo
        endif
c      write(6,*)'st_g2g4b',iplplwei,iplprwei,vector2(137)
        if ( logic_g2g4b ) then
          ii=ii0
          mm0=iplprwei
          nn0=iplplwei+iwt_sm_s_ext
          mm=mm0             !severe_error
          do ismb=1,ng_sm
            isma=mul_tab(ismb,ism_g2g4)
            if ( isma .gt. ismb ) cycle
            ibsta=ibsm_ext(ismb)
            ibend=iesm_ext(ismb)
            iasta=ibsm_ext(isma)
            iaend=iesm_ext(isma)
            if ( ismb .eq. isma ) ibsta=ibsta+1
            nna=nn0+ibsta-1
            do ib=ibsta,ibend
              nna=nna+1
              nnb=nn0+iasta-1
              valuelptmp1=vector1(nna)
              valuelptmp2=vector2(nna)
              do ia=iasta,min(iaend,ib-1)
                valuelp2=value_lpext(ii)
                valuelp1=value_lpext(ii+1)
                mm=mm+1
                nnb=nnb+1
                vector2(mm)=vector2(mm)+valuelptmp1*valuelp1
     *                     +vector1(nnb)*valuelp2
                vector2(nnb)=vector2(nnb)+vector1(mm)*valuelp2
                valuelptmp2=valuelptmp2+vector1(mm)*valuelp1
                ii=ii+2
              enddo
              vector2(nna)=valuelptmp2
            enddo
          enddo
        endif
        ii0=ii-1

c      write(6,*)'st_g36a',iplplwei,iplprwei,vector2(137)
        do icle=1,2
          if ( icle.eq.1 .and. logic_g36a ) then
            mm0=iplplwei
            nn0=iplprwei
            lpsta36=lpsta36a
            lpend36=lpend36a
          elseif (  icle.eq.2 .and. logic_g36b ) then
            mm0=iplprwei
            nn0=iplplwei
            lpsta36=lpsta36b
            lpend36=lpend36b
          else
            goto 936
          endif
          do iii=lpsta36,lpend36,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            ii=ii0+lpext_wei(iii+2)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            do i=1,lpext_wei(iii+3)
              vector2(mm)=vector2(mm)+vector1(nn)*value_lpext(ii)
              vector2(nn)=vector2(nn)+vector1(mm)*value_lpext(ii)
              mm=mm+1
              nn=nn+1
            enddo
          enddo
936       continue

c      write(6,*)'st_g35a',iplplwei,iplprwei,vector2(137)
          if ( icle.eq.1 .and. logic_g35a ) then
            mm0=iplplwei
            nn0=iplprwei
            lpsta35=lpsta35a
            lpend35=lpend35a
            do iii=lpsta35,lpend35,4
              mm=mm0+lpext_wei(iii)
              nn=nn0+lpext_wei(iii+1)
              iij=ii0+lpext_wei(iii+2)
              valuelp=vector2(nn)
              valuelptmp1=vector1(nn)
              do i=1,lpext_wei(iii+3)
                valtmp=value_lpext(iij)
                vector2(mm)=vector2(mm)+valuelptmp1*valtmp
                valuelp=valuelp+vector1(mm)*valtmp
                iij=iij+1
                mm=mm+1
              enddo
              vector2(nn)=valuelp
            enddo
          elseif (  icle.eq.2 .and. logic_g35b ) then
            mm0=iplprwei
            nn0=iplplwei
            lpsta35=lpsta35b
            lpend35=lpend35b
            do iii=lpsta35,lpend35,4
              mm=mm0+lpext_wei(iii)
              nn=nn0+lpext_wei(iii+1)
              iij=ii0+lpext_wei(iii+2)
              valuelp=vector2(nn)
              valuelptmp1=vector1(nn)
              do i=1,lpext_wei(iii+3)
                valtmp=-value_lpext(iij)
                vector2(mm)=vector2(mm)+valuelptmp1*valtmp
                valuelp=valuelp+vector1(mm)*valtmp
                iij=iij+1
                mm=mm+1
              enddo
              vector2(nn)=valuelp      !severe_error
            enddo
          endif

c      write(6,*)'st_g34',iplplwei,iplprwei,vector2(137)
          if ( icle.eq.1 .and. logic_g34a ) then
            mm0=iplplwei
            nn0=iplprwei
            lpsta34=lpsta34a
            lpend34=lpend34a
          elseif (  icle.eq.2 .and. logic_g34b ) then
            mm0=iplprwei
            nn0=iplplwei
            lpsta34=lpsta34b
            lpend34=lpend34b
          else
            goto 934
          endif
          do iii=lpsta34,lpend34,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            iij=ii0+lpext_wei(iii+2)
            valuelp=vector2(nn)
            valuelptmp1=vector1(nn)
            do i=1,lpext_wei(iii+3)
              valuetmp=-value_lpext(iij)
              vector2(mm)=vector2(mm)+valuelptmp1*valuetmp
              valuelp=valuelp+vector1(mm)*valuetmp
              iij=iij+1
              mm=mm+1
            enddo
            vector2(nn)=valuelp
          enddo
934       continue
          ii0=ii0+nvalue_space_ss
        enddo

      enddo

c...end of inn_ext_st_loop_unpack
      end

c*******************************************************************
      subroutine inn_ext_st_drl_loop_unpack(iplplweiin,iplprweiin)
c*******************************************************************
c     26 feb 2007 - revised
c
#include "drt_h.fh"
#include "lpextmode_h.fh"
c     write(6,*) ' st_test 2/2'

      do irot=1,mcroot
        irtidx=indx(irot)
        iplplwei=iplplweiin+irtidx
        iplprwei=iplprweiin+irtidx

        ii=1
        if ( logic_g1415 ) then
          mm=iplplwei
          nn=iplprwei
          do i=1,idownwei_g131415
            valuelp=value_lpext(ii)
            mm=mm+1
            nn=nn+1
            vector2(mm)=vector2(mm)+vector1(nn)*valuelp
            vector2(nn)=vector2(nn)+vector1(mm)*valuelp
            ii=ii+1
          enddo
        endif

        ii0=ii
        if ( logic_g2g4a ) then
          ii=ii0
          mm0=iplplwei
          nn0=iplprwei+iwt_sm_s_ext
          mm=mm0            !severe_error
          do ismb=1,ng_sm
            isma=mul_tab(ismb,ism_g2g4)
            if ( isma .gt. ismb ) cycle
            ibsta=ibsm_ext(ismb)
            ibend=iesm_ext(ismb)
            iasta=ibsm_ext(isma)
            iaend=iesm_ext(isma)
            if ( ismb .eq. isma ) ibsta=ibsta+1
            nna=nn0+ibsta-1
            do ib=ibsta,ibend
              nna=nna+1
              nnb=nn0+iasta-1
              valuelptmp1=vector1(nna)
              valuelptmp2=vector2(nna)
              do ia=iasta,min(iaend,ib-1)
                valuelp1=value_lpext(ii)
                valuelp2=value_lpext(ii+1)
                mm=mm+1
                nnb=nnb+1
                vector2(mm)=vector2(mm)+valuelptmp1*valuelp1
     *                     +vector1(nnb)*valuelp2
                vector2(nnb)=vector2(nnb)+vector1(mm)*valuelp2
                valuelptmp2=valuelptmp2+vector1(mm)*valuelp1
                ii=ii+2
              enddo
              vector2(nna)=valuelptmp2
            enddo
          enddo
        endif
        if ( logic_g2g4b ) then
          ii=ii0
          mm0=iplprwei
          nn0=iplplwei+iwt_sm_s_ext
          mm=mm0             !severe_error
          do ismb=1,ng_sm
            isma=mul_tab(ismb,ism_g2g4)
            if ( isma .gt. ismb ) cycle
            ibsta=ibsm_ext(ismb)
            ibend=iesm_ext(ismb)
            iasta=ibsm_ext(isma)
            iaend=iesm_ext(isma)
            if ( ismb .eq. isma ) ibsta=ibsta+1
            nna=nn0+ibsta-1
            do ib=ibsta,ibend
              nna=nna+1
              nnb=nn0+iasta-1
              valuelptmp1=vector1(nna)
              valuelptmp2=vector2(nna)
              do ia=iasta,min(iaend,ib-1)
                valuelp2=value_lpext(ii)
                valuelp1=value_lpext(ii+1)
                mm=mm+1
                nnb=nnb+1
                vector2(mm)=vector2(mm)+valuelptmp1*valuelp1
     *                     +vector1(nnb)*valuelp2
                vector2(nnb)=vector2(nnb)+vector1(mm)*valuelp2
                valuelptmp2=valuelptmp2+vector1(mm)*valuelp1
                ii=ii+2
              enddo
              vector2(nna)=valuelptmp2
            enddo
          enddo
        endif
        ii0=ii-1

        do icle=1,2
          if ( icle.eq.1 .and. logic_g36a ) then
            mm0=iplplwei
            nn0=iplprwei
            lpsta36=lpsta36a
            lpend36=lpend36a
          elseif (  icle.eq.2 .and. logic_g36b ) then
            mm0=iplprwei
            nn0=iplplwei
            lpsta36=lpsta36b
            lpend36=lpend36b
          else
            goto 936
          endif
          do iii=lpsta36,lpend36,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            ii=ii0+lpext_wei(iii+2)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            do i=1,lpext_wei(iii+3)
              valuetmp=value_lpext(ii)
              vector2(mm)=vector2(mm)+vector1(nn)*valuetmp
              vector2(nn)=vector2(nn)+vector1(mm)*valuetmp
              mm=mm+1
              nn=nn+1
            enddo
          enddo
936       continue
          if ( icle.eq.1 .and. logic_g35a ) then
            mm0=iplplwei
            nn0=iplprwei
            lpsta35=lpsta35a
            lpend35=lpend35a
            do iii=lpsta35,lpend35,4
              mm=mm0+lpext_wei(iii)
              nn=nn0+lpext_wei(iii+1)
              iij=ii0+lpext_wei(iii+2)
              valuelp=vector2(nn)
              valuelptmp1=vector1(nn)
              do i=1,lpext_wei(iii+3)
                valuetmp=value_lpext(iij)
                vector2(mm)=vector2(mm)+valuelptmp1*valuetmp
                valuelp=valuelp+vector1(mm)*valuetmp
                iij=iij+1
                mm=mm+1
              enddo
              vector2(nn)=valuelp
            enddo
          elseif (  icle.eq.2 .and. logic_g35b ) then
            mm0=iplprwei
            nn0=iplplwei
            lpsta35=lpsta35b
            lpend35=lpend35b
            do iii=lpsta35,lpend35,4
              mm=mm0+lpext_wei(iii)
              nn=nn0+lpext_wei(iii+1)
              iij=ii0+lpext_wei(iii+2)
              valuelp=vector2(nn)         !severe_error
              valuelptmp1=vector1(nn)
              do i=1,lpext_wei(iii+3)
                valuetmp=-value_lpext(iij)
                vector2(mm)=vector2(mm)+valuelptmp1*valuetmp
                valuelp=valuelp+vector1(mm)*valuetmp
                iij=iij+1
                mm=mm+1
              enddo
              vector2(nn)=valuelp      !severe_error
            enddo
          endif

          if ( icle.eq.1 .and. logic_g34a ) then
            mm0=iplplwei
            nn0=iplprwei
            lpsta34=lpsta34a
            lpend34=lpend34a
          elseif (  icle.eq.2 .and. logic_g34b ) then
            mm0=iplprwei
            nn0=iplplwei
            lpsta34=lpsta34b
            lpend34=lpend34b
          else
            goto 934
          endif
          do iii=lpsta34,lpend34,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            iij=ii0+lpext_wei(iii+2)
            valuelp=vector2(nn)
            valuelptmp1=vector1(nn)
            do i=1,lpext_wei(iii+3)
              valuetmp=-value_lpext(iij)
              vector2(mm)=vector2(mm)+valuelptmp1*valuetmp
              valuelp=valuelp+vector1(mm)*valuetmp
              iij=iij+1
              mm=mm+1
            enddo
            vector2(nn)=valuelp
          enddo
934       continue
        enddo
      enddo

c...end of inn_ext_st_drl_loop_unpack
      end


c**************************************************************
      subroutine gsd_sequence_extspace(iplplweiin,iplprweiin)
c**************************************************************
c     26 feb 2007 - revised
c
#include "drt_h.fh"
#include "lpextmode_h.fh"
      parameter (v_sqtwo=1.414213562373095d0 )
c      write(6,*) '  sd_test 1/2','  ds_test 011'

      do irot=1,mcroot
        irtidx=indx(irot)
        iplplwei=iplplweiin+irtidx
        iplprwei=iplprweiin+irtidx
        ilpvalue=0
        if ( logic_g25a ) then
          mm=iplplwei+iweista_g25-1
          nn0=iplprwei
          do itmp=1,nint_g25
            ilpvalue=ilpvalue+1
            vlptmp=value_lpext(ilpvalue)
            nn=nn0
            do i=1,nwei_g25
              mm=mm+1
              nn=nn+1
              vector2(mm)=vector2(mm)+vector1(nn)*vlptmp
              vector2(nn)=vector2(nn)+vector1(mm)*vlptmp
            enddo
          enddo
        elseif ( logic_g25b ) then
          mm=iplplwei+iweista_g25-1
          nn0=iplprwei
          ilpvalue=ilpvalue+1
          do itmp=2,nint_g25
            ilpvalue=ilpvalue+1
            vlptmp=value_lpext(ilpvalue)
            nn=nn0
            do i=1,itmp-1
              mm=mm+1
              nn=nn+1
              vector2(mm)=vector2(mm)+vector1(nn)*vlptmp
              vector2(nn)=vector2(nn)+vector1(mm)*vlptmp
            enddo
          enddo

          mm=iplplwei+iweista_g28-1
          nn=iplprwei
          nn=nn+1
          do itmp=2,nwei_g28
            nn=nn+1
            vlptmp=vector2(nn)
            vlptmp1=vector1(nn)
            ilpvalue=0
            do i=1,itmp-1
              ilpvalue=ilpvalue+1
              vtmp=value_lpext(ilpvalue)
              mm=mm+1
              vector2(mm)=vector2(mm)+vlptmp1*vtmp
              vlptmp=vlptmp+vector1(mm)*vtmp
            enddo
            vector2(nn)=vlptmp
          enddo
        elseif ( logic_g28a ) then
          mm=iplplwei+iweista_g28-1
          nn0=iplprwei
          do nn=nn0+1,nn0+nwei_g28
            vlptmp=vector2(nn)
            vlptmp1=vector1(nn)
            ilpvalue=0
            do i=1,nint_g28
              ilpvalue=ilpvalue+1
              mm=mm+1
              vector2(mm)=vector2(mm)+
     *                    vlptmp1*value_lpext(ilpvalue)
              vlptmp=vlptmp+
     *               vector1(mm)*value_lpext(ilpvalue)
            enddo
            vector2(nn)=vlptmp
          enddo
        endif

        if ( logic_g26 ) then
          ilpvalue=ivaluesta_g26
          mm=iplplwei+iweista_g26
          nn0=iplprwei
          do nn=nn0+1,nn0+nwei_g26
            ilpvalue=ilpvalue+1
            vlptmp=value_lpext(ilpvalue)*v_sqtwo   !for g26
            vector2(mm)=vector2(mm)+vector1(nn)*vlptmp
            vector2(nn)=vector2(nn)+vector1(mm)*vlptmp
            mm=mm+1
          enddo
        endif
      enddo

c...end of gsd_sequence_extspace
      end

c***********************************************************************
      subroutine complete_sd_ar_ext_loop(ilweiin,irweiin,isdownwei)
c***********************************************************************
c     26 feb 2007 - revised
c
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
c      write(6,*) 'sd_test 2/2','  td_test_2/2 012'


      do irot=1,mcroot
        irtidx=indx(irot)
        ilwei=ilweiin+irtidx
        irwei=irweiin+irtidx

        ilpvalue=0
        mm0=ilwei
        nn=irwei+icano_nnsta-1
        do nntmp=icano_nnsta,icano_nnend
          nn=nn+1
          mm=mm0
          vlptmp1=vector1(nn)
          vlptmp=vector2(nn)
          do mmtmp=1,isdownwei
            ilpvalue=ilpvalue+1
            mm=mm+1
            vector2(mm)=vector2(mm)+vlptmp1*value_lpext(ilpvalue)
            vlptmp=vlptmp+vector1(mm)*value_lpext(ilpvalue)
          enddo
          vector2(nn)=vlptmp
        enddo
      enddo
      return

c...end of complete_sd_ar_ext_loop
      end

c******************************************************************
      subroutine inn_ext_dd_loop_unpack(iplplweiin,iplprweiin)
c******************************************************************
c     26 feb 2007 - revised
c
#include "drt_h.fh"
#include "lpextmode_h.fh"
c        'dd_test '

      do irot=1,mcroot
        irtidx=indx(irot)
        iplplwei=iplplweiin+irtidx
        iplprwei=iplprweiin+irtidx

        ii=1
        if ( logic_g50 ) then         !arbl=.true.
          if ( logic_g49b) then
            mm=iplplwei
            nn=iplprwei
            do i=1,ildownwei_segdd
              valuelp=value_lpext(ii)
              mm=mm+1
              nn=nn+1
              vector2(mm)=vector2(mm)+vector1(nn)*valuelp
              vector2(nn)=vector2(nn)+vector1(mm)*valuelp
              ii=ii+1
            enddo
          endif

          ii=ii+int_dd_drl
          mm0=iplplwei
          nn=iplprwei+1
          do icle=1,2
            do j=2,ildownwei_segdd
              nn=nn+1
              valuelp=vector2(nn)
              valuelptmp1=vector1(nn)
              mm=mm0
              do i=1,j-1
                mm=mm+1
                valuetmp=value_lpext(ii)
                vector2(mm)=vector2(mm)+valuelptmp1*valuetmp
                valuelp=valuelp+vector1(mm)*valuetmp
                ii=ii+1
              enddo
              vector2(nn)=valuelp
            enddo
            if ( .not. logic_g49b ) exit
            mm0=iplprwei
            nn=iplplwei+1
          enddo

        else               !drl=.true.


          ii=ii+int_dd_drl
          if ( logic_g49a ) then
            mm0=iplplwei
            nn=iplprwei
            ildownwei=ildownwei_segdd
            irdownwei=irdownwei_segdd
          else
            mm0=iplprwei
            nn=iplplwei
            ildownwei=irdownwei_segdd
            irdownwei=ildownwei_segdd
          endif
          do j=1,irdownwei
            nn=nn+1
            valuelp=vector2(nn)
            valuelptmp1=vector1(nn)
            mm=mm0
            do i=1,ildownwei
              mm=mm+1
              valuetmp=value_lpext(ii)
              vector2(mm)=vector2(mm)+valuelptmp1*valuetmp
              valuelp=valuelp+vector1(mm)*valuetmp
              ii=ii+1
            enddo
            vector2(nn)=valuelp
          enddo
        endif
      enddo

c...end of inn_ext_dd_loop_unpack
      end

c***********************************************************************
      subroutine inn_ext_tt_drl_loop_unpack(iplplweiin,iplprweiin,n1415)
c***********************************************************************
c     26 feb 2007 - revised
c
#include "drt_h.fh"
#include "lpextmode_h.fh"
      logical logic_g14150,logic_g34b0,logic_g35b0,logic_g36b0
c     write(6,*) '  tt_test 2/2'

      logic_g14150=logic_g1415
      logic_g36b0=logic_g36b
      logic_g35b0=logic_g35b
      logic_g34b0=logic_g34b
      if(iplplweiin.eq.iplprweiin) then
        logic_g14150=.false.
        logic_g36b0=.false.
        logic_g35b0=.false.
        logic_g34b0=.false.
      endif

      do irot=1,mcroot
        irtidx=indx(irot)
        iplplwei=iplplweiin+irtidx
        iplprwei=iplprweiin+irtidx

        ii=1
        if ( logic_g14150 ) then
          mm=iplplwei
          nn=iplprwei
          do i=1,idownwei_g131415
            valuelp=value_lpext(ii)
            mm=mm+1
            nn=nn+1
            vector2(mm)=vector2(mm)+vector1(nn)*valuelp
            vector2(nn)=vector2(nn)+vector1(mm)*valuelp
            ii=ii+1
          enddo
        endif

        ii0=n1415         !severe_error

        do icle=1,2
          if ( icle.eq.1 .and. logic_g36a ) then
            mm0=iplplwei
            nn0=iplprwei
            lpsta36=lpsta36a
            lpend36=lpend36a
          elseif (  icle.eq.2 .and. logic_g36b0 ) then
            mm0=iplprwei
            nn0=iplplwei
            lpsta36=lpsta36b
            lpend36=lpend36b
          else
            goto 936
          endif
          do iii=lpsta36,lpend36,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            ii=ii0+lpext_wei(iii+2)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            do i=1,lpext_wei(iii+3)
              valuelp=value_lpext(ii)
              vector2(mm)=vector2(mm)+vector1(nn)*valuelp
              vector2(nn)=vector2(nn)+vector1(mm)*valuelp
              mm=mm+1
              nn=nn+1
            enddo
          enddo
936       continue

          if ( icle.eq.1 .and. logic_g35a ) then
            mm0=iplplwei
            nn0=iplprwei
            lpsta35=lpsta35a
            lpend35=lpend35a
          elseif (  icle.eq.2 .and. logic_g35b0 ) then
            mm0=iplprwei
            nn0=iplplwei
            lpsta35=lpsta35b
            lpend35=lpend35b
          else
            goto 935
          endif

          do iii=lpsta35,lpend35,4
            mm=mm0+lpext_wei(iii)
            nn=nn0+lpext_wei(iii+1)
            iij=ii0+lpext_wei(iii+2)
            valuelp=vector2(nn)
            valuelptmp1=vector1(nn)
            do i=1,lpext_wei(iii+3)
              valuetmp=-value_lpext(iij)
              vector2(mm)=vector2(mm)+valuelptmp1*valuetmp
              valuelp=valuelp+vector1(mm)*valuetmp
              iij=iij+1
              mm=mm+1
            enddo
            vector2(nn)=valuelp      !severe_error_1202
          enddo
935       continue

          if ( icle.eq.1 .and. logic_g34a ) then
            mm0=iplplwei
            nn0=iplprwei
            lpsta34=lpsta34a
            lpend34=lpend34a
          elseif (  icle.eq.2 .and. logic_g34b0 ) then
            mm0=iplprwei
            nn0=iplplwei
            lpsta34=lpsta34b
            lpend34=lpend34b
          else
            goto 934
          endif
          do iii=lpsta34,lpend34,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            iij=ii0+lpext_wei(iii+2)
            valuelp=vector2(nn)
            valuelptmp1=vector1(nn)
            do i=1,lpext_wei(iii+3)
              vector2(mm)=vector2(mm)+valuelptmp1*value_lpext(iij)
              valuelp=valuelp+vector1(mm)*value_lpext(iij)
              iij=iij+1
              mm=mm+1
            enddo
            vector2(nn)=valuelp
          enddo
934       continue
        enddo
      enddo

c...end of inn_ext_tt_drl_loop_unpack
      end

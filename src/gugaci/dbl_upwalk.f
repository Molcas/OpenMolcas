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
      subroutine dbl_upwalk()
#include "drt_h.fh"
#include "intsort_h.fh"

      if(norb_dbl.eq.1) then
c     v(1),d(2-9),s(18-25)           for s=0
c     v(1),d(2-9),s(18-25),d'(26-33)   for s<>0

        mxnode=17+ng_sm
        lri=norb_frz+1
        lsmi=lsm_inn(lri)
        lsmid=mul_tab(lsmi,ns_sm)
!for node_v
        nu_ad(1)=1
        jpad_upwei(1)=1
!for node_d
        nu_ad(1+lsmid)=1+lsmid
        jpad_upwei(1+lsmid)=1
!for node_s
        nu_ad(17+ns_sm)=17+ns_sm
        jpad_upwei(17+ns_sm)=1

        if(jroute_sys.eq.1) return
        mxnode=25+ng_sm
!for node_d'
        nu_ad(25+lsmid)=25+lsmid
        jpad_upwei(25+lsmid)=1
        return
      endif
      nu_ad=0
      jpad_upwei=0

      nu_ad(1)=1
      jpad_upwei(1)=1
      if(norb_dbl.eq.0) then
        mxnode=1
        return
      endif
      do lri=norb_frz+1,norb_dz
        lsmi=lsm_inn(lri)
        lsmid=mul_tab(lsmi,ns_sm)
        no_d=lsmid+1
        jpad_upwei(no_d)=jpad_upwei(no_d)+1
        do lrj=lri+1,norb_dz
          lsmj=lsm_inn(lrj)
          lsmij=mul_tab(lsmi,lsmj)
          lsmit=mul_tab(lsmij,ns_sm)
          no_t =lsmit+9
          jpad_upwei(no_t)=jpad_upwei(no_t)+1
        enddo
      enddo
c     v(1),d(2-9),t(10-17),s(18-25),d'(26-33),t'(34-41)
      goto(100,200,300) jroute_sys
100     mxnode=25                     !v,d,t,s
        jpad_upwei(18:25)=jpad_upwei(10:17)
        jpad_upwei(17+ns_sm)=jpad_upwei(17+ns_sm)+norb_dbl
        goto 500
200     mxnode=25+8
        jpad_upwei(18:25)=jpad_upwei(10:17)+jpad_upwei(10:17)
        jpad_upwei(17+ns_sm)=jpad_upwei(17+ns_sm)+norb_dbl
        jpad_upwei(26:33)=jpad_upwei(2:9)
        goto 500
300     mxnode=25+8+8
        jpad_upwei(18:25)=jpad_upwei(10:17)+jpad_upwei(10:17)
        jpad_upwei(17+ns_sm)=jpad_upwei(17+ns_sm)+norb_dbl
        jpad_upwei(26:33)=jpad_upwei(2:9)
        jpad_upwei(34:41)=jpad_upwei(10:17)

500   do node=2,mxnode
        iw=jpad_upwei(node)
        if(iw.eq.0) cycle
        nu_ad(node)=node
      enddo
      return
      end

      subroutine ext_downwalk()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
!      common/casrst/ja(max_node),jb(max_node),jm(0:max_node)
!     :    ,jj(4,0:max_node),kk(0:max_node),no(0:max_innorb)
!     :    ,jv,jd(8),jt(8),js(8)
      dimension iwmij(8)
      nu_ae(1)=1
      do im=1,ng_sm
        nu_ae(1+im)=1+im
        nu_ae(9+im)=9+im
        nu_ae(17+im)=17+im
      enddo

      iwmij=0
      iseg_downwei(nu_ae(1))=1
      do imi=1,ng_sm
        iseg_downwei(nu_ae(1+imi))=nlsm_ext(imi)
        do imj=imi,ng_sm
          imij=mul_tab(imi,imj)
          if(imij.ne.1) then
           iwmij(imij)=iwmij(imij)+nlsm_ext(imi)*nlsm_ext(imj)
             cycle
          endif
         iwmij(1)=iwmij(1)+nlsm_ext(imi)*(nlsm_ext(imi)-1)/2
        enddo
      enddo
      do im=1,ng_sm
        iseg_downwei(nu_ae(9+im))=iwmij(im)
        iseg_downwei(nu_ae(17+im))=iwmij(im)
      enddo
      iseg_downwei(nu_ae(18))=iseg_downwei(nu_ae(18))+norb_ext
      return
      end

      subroutine readdrt(ludrt)
#include "drt_h.fh"
      !include "files_gugaci.fh"
#include "pl_structure_h.fh"
!      common/casrst/ja(max_node),jb(max_node),jm(0:max_node)
!     :    ,jj(4,0:max_node),kk(0:max_node),no(0:max_innorb)
!     :    ,jv,jd(8),jt(8),js(8)
      dimension idx(2)

      idisk=0
      !  number of nodes
      call idafile(ludrt,2,idx,2,idisk)
      idisk=idx(2)
      call idafile(ludrt,2,idx,1,idisk)
      id=idx(1)
      call idafile(ludrt,2,ja,id,idisk)
      call idafile(ludrt,2,jb,id,idisk)
      call idafile(ludrt,2,jm,id,idisk)
      call idafile(ludrt,2,jj,4*(id+1),idisk)
      call idafile(ludrt,2,kk,1+id,idisk)
      call idafile(ludrt,2,no(0),norb_inn+2,idisk)
      call idafile(ludrt,2,idx,1,idisk)
      jv=idx(1)
      call idafile(ludrt,2,jd,8,idisk)
      call idafile(ludrt,2,jt,8,idisk)
      call idafile(ludrt,2,js,8,idisk)

      return
      end

c     juv,just(nost,nost),jud(nost)
c     |  \  1         |
c     | d,dd,s(i=i)   |
c     |    \ s,t,tt(i<j)|
c     |     \       1 2 |     deal with inner of dbl_space
c     |ss(i>j)\       |
c     |  2 1  \       |
      subroutine dbl_downwalk()
#include "drt_h.fh"
#include "intsort_h.fh"
c     integer lsml(10,10)       !to del
      if(norb_dbl.ne.0) goto 200
c----------- norb_dbl=0 ------------------------------------------------
      do im=1,ng_sm
        nnd=iseg_sta(1+im)
        nnt=iseg_sta(9+im)
        nns=iseg_sta(17+im)
        do lri=norb_dz,norb_frz+1,-1
          ismi=lsm_inn(lri)
         if(ismi.ne.im) cycle
          jud(lri)=nnd
          nnd=nnd+iseg_downwei(1+im)
        enddo
        do lrj=norb_dz,norb_frz+1,-1
         ismj=lsm_inn(lrj)
          do lri=lrj,1,-1
           ismi=lsm_inn(lri)
           ismij=mul_tab(ismi,ismj)
            if(ismij.ne.im) cycle
            just(lri,lrj)=nns
            nns=nns+iseg_downwei(17+im)
            if(lri.eq.lrj) cycle
           just(lrj,lri)=nnt
            nnt=nnt+iseg_downwei(9+im)
           enddo
        enddo
      enddo
c----------- norb_dbl=0 ------------------------------------------------
c----------- norb_dbl<>0 -----------------------------------------------
200   continue
      do im=1,ng_sm
        nnd=0
        nns=0
        do lri=norb_frz+1,norb_dz
          ismi=mul_tab(lsm_inn(lri),ns_sm)
          if(ismi.ne.im) cycle
          jud(lri)=nnd
          nnd=nnd+1
        enddo
        do lri=norb_frz+1,norb_dz-1
          ismi=mul_tab(lsm_inn(lri),ns_sm)
          do lrj=lri+1,norb_dz      !tmp
            ismj=lsm_inn(lrj)
            ismij=mul_tab(ismi,ismj)
            if(ismij.ne.im) cycle
            just(lri,lrj)=nns
            nns=nns+1
          enddo
        enddo
        if(im.eq.ns_sm) then
          do lr0=norb_frz+1,norb_dz
            just(lr0,lr0)=nns
            nns=nns+1
          enddo
        endif
        do lri=norb_frz+1,norb_dz-1
          ismi=mul_tab(lsm_inn(lri),ns_sm)
          do lrj=lri+1,norb_dz      !tmp
           ismj=lsm_inn(lrj)
           ismij=mul_tab(ismi,ismj)
            if(ismij.ne.im) cycle
            just(lrj,lri)=nns
            nns=nns+1
         enddo
        enddo
      enddo
!      print*, "ns_sm",ns_sm
!      print*, just(1:4,1:4)
c     do i=norb_frz+1,norb_dz           !to del
c       lmi=lsm_inn(i)               !to del
c        do j=norb_frz+1,norb_dz         !to del
c       lmj=lsm_inn(j)                !to del
c       lsml(i,j)=mul_tab(lmi,lmj)       !to del
c        enddo                         !to del
c     enddo                         !to del
c     write(nf2,*)'   jud ...'             !to del
c     write(nf2,'(2x,12i8)')(jud(lr),lr=norb_frz+1,norb_dz)    !to del
c     write(nf2,*)'   just ...'           !to del
c     do i=norb_frz+1,norb_dz               !to del
c      write(nf2,'(2x,8i3)')             !to del
c    :         (lsml(i,lr),lr=norb_frz+1,norb_dz)  !to del
c      write(nf2,'(2x,8i3)')             !to del
c    :         (just(i,lr),lr=norb_frz+1,norb_dz)  !to del
c      write(nf2,*)
c     enddo                            !to del
c----------- norb_dbl<>0 -----------------------------------------------
      return
      end

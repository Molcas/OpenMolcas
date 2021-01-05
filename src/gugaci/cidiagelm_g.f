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
      subroutine diagonal_loop_wyb_g()  !  for norb_act<>0
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
!      common/casrst/ja(max_node),jb(max_node),jm(0:max_node)
!     :    ,jj(4,0:max_node),kk(0:max_node),no(0:max_innorb)
!     :    ,jv,jd(8),jt(8),js(8)
!      common/sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),
!     :     jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)
c      do lr0=2,norb_all
c        do lr=1,lr0-1
c          vdint(lr,lr0)=voint(lr0,lr)
c     :      -vdint(lr0,lr)-vdint(lr0,lr)   ! 520
c      write(*,'(2i4,3f14.8)')
c     :   lr,lr0,voint(lr0,lr),vdint(lr0,lr),vdint(lr,lr0)
c        enddo
c      enddo
c     write(*,*)'               ***** start h-diaelm *****'
c**************lyb***************
c     vector1(1:nci_dim)=vpotnuc
c********************************
!      wl8 = hnil*(hnil-1)*vmd(lr,lr)*0.5d0+hnil*vo(lr,lr)   800
c      write(*,*)'         jpad,     jpae,     ndim,      nohy'

      ndimsum=0
      jpae=jv
      ipae=1
      jaedownwei=iseg_downwei(ipae)
      do jpad_=1,mxnode
        jpad=jpad_ ! jpad is in common block, is this necessary?
        iw_sta(jpad,ipae)=ndimsum
        if(nu_ad(jpad).eq.0) cycle
        call seg_drt()
        iwupwei=jpad_upwei(jpad)
        iw_downwei(jpad,ipae)=ndim
        ndimsum=ndimsum+ndim*jaedownwei*iwupwei
        if(ndim .eq. 0) cycle
        call diagonal_act_d_g()
        call diagonal_act_c_g()
      enddo
      do im=1,ng_sm
        jpae=jd(im)
        ipae=1+im
        if(nu_ae(ipae).eq.0) cycle
        jaedownwei=iseg_downwei(ipae)
        do jpad_=1,mxnode
          jpad=jpad_ ! jpad is in common block, is this necessary?
          iw_sta(jpad,ipae)=ndimsum
          if(nu_ad(jpad).eq.0) cycle
          call seg_drt()
          iwupwei=jpad_upwei(jpad)
          iw_downwei(jpad,ipae)=ndim
c       if(jpad.ge.26) then
c     write(*,*)
c     endif
          ndimsum=ndimsum+ndim*jaedownwei*iwupwei
          if(ndim .eq. 0) cycle
          call diagonal_act_d_g()
          call diagonal_act_c_g()
        enddo
      enddo

      do im=1,ng_sm
        jpae=jt(im)
        ipae=9+im
        if(nu_ae(ipae).eq.0) cycle
        jaedownwei=iseg_downwei(ipae)
        do jpad_=1,mxnode
          jpad=jpad_ ! jpad is in common block, is this necessary?
          iw_sta(jpad,ipae)=ndimsum
          if(nu_ad(jpad).eq.0) cycle
          call seg_drt()
          iwupwei=jpad_upwei(jpad)
         iw_downwei(jpad,ipae)=ndim
          ndimsum=ndimsum+ndim*jaedownwei*iwupwei
          if(ndim .eq. 0) cycle
          call diagonal_act_d_g()
          call diagonal_act_c_g()
        enddo
      enddo
      do im=1,ng_sm
        jpae=js(im)
        ipae=17+im
        jaedownwei=iseg_downwei(ipae)
        if(nu_ae(ipae).eq.0) cycle
        do jpad_=1,mxnode
          jpad=jpad_ ! jpad is in common block, is this necessary?
          iw_sta(jpad,ipae)=ndimsum
          if(nu_ad(jpad).eq.0) cycle
          call seg_drt()
          iwupwei=jpad_upwei(jpad)
          iw_downwei(jpad,ipae)=ndim
          ndimsum=ndimsum+ndim*jaedownwei*iwupwei
          if(ndim .eq. 0) cycle
          call diagonal_act_d_g()
          call diagonal_act_c_g()
        enddo
      enddo
      call diagonal_dbl_g()
      call diagonal_ext_g()
      return
      end

      subroutine diagonal_act_c_g()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
!     common/casrst/ja(max_node),jb(max_node),jm(0:max_node)
!     :    ,jj(4,0:max_node),kk(0:max_node),no(0:max_innorb)
!     :    ,jv,jd(8),jt(8),js(8)
!      common/sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),
!     :     jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)
      dimension ndr(max_innorb)
      integer, pointer :: jph(:),jeh(:),jwh(:)
      real*8, pointer :: th(:),thh(:)
      common/ptlph/jph,jeh,jwh
      common/ptlphv/th,thh
      real*8, allocatable :: te(:), tee(:)
      integer, allocatable :: jpe(:), jee(:), jwe(:)
c     write(*,*)'               ***** start h-diaelm *****'
c      write(*,*)  jpad,jpae
      allocate(te(maxpl),tee(maxpl),jpe(maxpl),jee(maxpl),jwe(maxpl))
      ndr=0
      if(norb_act.eq.0) then
        mh=1
        th(1)=1.d0
        thh(1)=1.d0
        call diagonal_link_dae_g(mh)
        return
      endif
      mp=0
c 520
          mh=0
          me=0
          jp=jpad
          jpb=jb(jp)
          do idl=1,4
            if(jj_sub(idl,jp).eq.0) cycle
           ind1=idl
c     link c"
            call smidc2(isq,w,ww,mw,ind1,jpb)
            mh=mh+1
            jeh(mh)=jj_sub(idl,jp)
            th(mh)=w
            thh(mh)=ww
            jph(mh)=0
            jwh(mh)=0
            if(idl.ne.1) jwh(mh)=iy(idl,jp)
c     complete a loop 'v'
            if(ind1.eq.1) cycle
            call stml(isq,w,ww,mw,ind1-1,jpb)
            vlop0=w
            vlop1=ww
            if(vlop0.eq.0.0d0.and.vlop1.eq.0.0d0) cycle
            mpe=jj_sub(idl,jp)
            iwa=iy(idl,jp)
            call diagonal_link_ad_g(mpe,iwa,vlop0,vlop1)
          enddo
c************************************************************
c      write(*,*)ad(i)
c************************************************************
        lr=norb_dz+1
40      if(lr.eq.norb_inn) then
          if(mh.ne.0) call diagonal_link_dae_g(mh)
          return
        end if
        lr=lr+1
        me=0
        do 160 m=1,mh
          je=jeh(m)
          jeb=jb(je)
c          jp=jph(m)
          do idl=1,4
            if(jj_sub(idl,je).eq.0) cycle
           ind1=idl
c            if(lr.eq.1) goto 20
c     link loop
            call smidc2(isq,w,ww,mw,ind1,jeb)
            me=me+1
            jwe(me)=jwh(m)
            if(idl.ne.1) jwe(me)=jwe(me)+iy(idl,je)
            jee(me)=jj_sub(idl,je)
            te(me)=th(m)*w
            tee(me)=thh(m)*ww
            jpe(me)=jph(m)
c     complete a loop 'v'
c20          continue
            if(ind1.eq.1) cycle
            call stml(isq,w,ww,mw,ind1-1,jeb)
            vlop0=th(m)*w
            vlop1=thh(m)*ww
            if(vlop0.eq.0.0d0.and.vlop1.eq.0.0d0) cycle
            mp=mp+1
            mpe=jj_sub(idl,je)
            iwa=jwh(m)
           if(idl.ne.1) iwa=iwa+iy(idl,je)
            call diagonal_link_ad_g(mpe,iwa,vlop0,vlop1)
c*****   520  ******************************************************
          enddo
160     continue
c***********************************************************
c      write(*,*) ad(i)
c************************************************************
        do m=1,me
          th(m)=te(m)
          te(m)=0.0d0
          thh(m)=tee(m)
          tee(m)=0.0d0
          jwh(m)=jwe(m)
          jwe(m)=0
          jeh(m)=jee(m)
          jee(m)=0
          jph(m)=jpe(m)
          jpe(m)=0
        enddo
        mh=me
        if(ndr(lr).lt.mh) ndr(lr)=mh
        goto 40
        do m=1,mh
          th(m)=0.0d0
          thh(m)=0.0d0
          jwh(m)=0
          jph(m)=0
          jeh(m)=0
        enddo
      return
      end

      subroutine diagonal_act_d_g()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
!     common/casrst/ja(max_node),jb(max_node),jm(0:max_node)
!     :    ,jj(4,0:max_node),kk(0:max_node),no(0:max_innorb)
!     :    ,jv,jd(8),jt(8),js(8)
!      common/sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),
!     :     jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)
      dimension ndr(max_innorb)
      integer, pointer :: jph(:),jeh(:),jwh(:)
      real*8, pointer :: th(:),thh(:)
      common/ptlph/jph,jeh,jwh
      common/ptlphv/th,thh
      real*8, allocatable :: te(:), tee(:)
      integer, allocatable :: jpe(:), jee(:), jwe(:)
c     write(*,*)'               ***** start h-diaelm *****'
c      write(*,*)   '   diagonal_act_d:',jpad,ipae
      allocate(te(maxpl),tee(maxpl),jpe(maxpl),jee(maxpl),jwe(maxpl))
      ndr=0
      do lr=norb_dz+1,norb_inn
        jp0=no(lr-1)+1
        jp1=no(lr)
        lr0=lr
        do jp=jp0,jp1
          if(iy(1,jp).eq.0) cycle
c  800
          do idl=2,3
            mpe=jj_sub(idl,jp)
            if(mpe.eq.0) cycle
c            wt = voint(lr0,lr0)    ! hnil=1
             wt=1.0d0
            jw=iy(idl,jp)
            call prodel_1(3,wt,jp,mpe,jw,lr0,lr0)
          enddo
          mpe=jj_sub(4,jp)
          if(mpe.ne.0) then
c            wt = vdint(lr0,lr0)+2.d0*voint(lr0,lr0)     !idl=4 hnil=2
            jw=iy(4,jp)
            wt=2.0d0
            call prodel_1(3,wt,jp,mpe,jw,lr0,lr0)
            wt=1.0d0
            call trans_ijkl_intpos(lr0,lr0,lr0,lr0,nxo)
            call prodel_2(3,wt,jp,mpe,jw,nxo)
          endif
        enddo
      enddo
c 520
      do 5 lr0=norb_dz+1,norb_inn
        mh=0
        me=0
        jp0=no(lr0-1)+1
        jp1=no(lr0)
        if(jp0.gt.jp1) jp0=jp1
        do 1 jp=jp0,jp1
          if(iy(1,jp).eq.0) cycle
c     start '^'
          do idl=1,4
            if(jj_sub(idl,jp).eq.0) cycle
            jbl=jb(jp)
            ind1=idl-1
            if(ind1.eq.0) cycle
            call stmh(isq,w,ww,mw,ind1,jbl)
            mh=mh+1
            jeh(mh)=jj_sub(idl,jp)
            th(mh)=w
            thh(mh)=ww
            jph(mh)=jp
            jwh(mh)=0
            if(idl.ne.1) jwh(mh)=iy(idl,jp)
          enddo
1       continue
c************************************************************
c      write(*,*)ad(i)
c************************************************************
        lr=lr0
        if(ndr(lr).lt.mh) ndr(lr)=mh
40      if(lr.eq.norb_inn) then
        call diagonal_link_ae_g(mh)
        goto 5
        end if
        lr=lr+1
        me=0

        mp=0
        do 160 m=1,mh
          je=jeh(m)
          jeb=jb(je)
          jp=jph(m)
          do 17 idl=1,4
            if(jj_sub(idl,je).eq.0) goto 17
            ind1=idl
            if(lr.eq.1) goto 20
c     link loop
            call smidc2(isq,w,ww,mw,ind1,jeb)
            me=me+1
            jwe(me)=jwh(m)
            if(idl.ne.1) jwe(me)=jwe(me)+iy(idl,je)
            jee(me)=jj_sub(idl,je)
            te(me)=th(m)*w
            tee(me)=thh(m)*ww
            jpe(me)=jph(m)
c     complete a loop 'v'
20          if(ind1.eq.1) goto 17
            call stml(isq,w,ww,mw,ind1-1,jeb)
            vlop0=th(m)*w
            vlop1=thh(m)*ww
            if(vlop0.eq.0.0d0.and.vlop1.eq.0.0d0) goto 17
            mp=mp+1
            mpe=jj_sub(idl,je)
            iwa=jwh(m)
           if(idl.ne.1) iwa=iy(idl,je)+iwa
c            wt=(vlop0-vlop1)*voint(lr,lr0)-2.d0*vlop0*vdint(lr,lr0)
            wt=vlop0-vlop1
            call trans_ijkl_intpos(lr,lr0,lr,lr0,nxo)
c            write(nf2,'(4i4,i8)') lr,lr0,lr,lr0,nxo
            call prodel_2(3,wt,jp,mpe,iwa,nxo)
            wt=-2.d0*vlop0
            call trans_ijkl_intpos(lr,lr,lr0,lr0,nxo)
            call prodel_2(3,wt,jp,mpe,iwa,nxo)

c*****   520  ******************************************************
17        continue
160     continue
c***********************************************************
c      write(*,*) ad(i)
c************************************************************
      do  m=1,me
          th(m)=te(m)
          te(m)=0.0d0
          thh(m)=tee(m)
          tee(m)=0.0d0
          jwh(m)=jwe(m)
          jwe(m)=0
          jeh(m)=jee(m)
          jee(m)=0
          jph(m)=jpe(m)
          jpe(m)=0
        enddo
        mh=me
        if(ndr(lr).lt.mh) ndr(lr)=mh
        goto 40
        do m=1,mh
          th(m)=0.0d0
          thh(m)=0.0d0
          jwh(m)=0
          jph(m)=0
          jeh(m)=0
        enddo
5     continue
      return
      end

      subroutine diagonal_link_ae_g(mh)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
!     common/casrst/ja(max_node),jb(max_node),jm(0:max_node)
!     :    ,jj(4,0:max_node),kk(0:max_node),no(0:max_innorb)
!     :    ,jv,jd(8),jt(8),js(8)
!      common/sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),
!     :     jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)
      integer, pointer :: jph(:),jeh(:),jwh(:)
      real*8, pointer :: th(:),thh(:)
      common/ptlph/jph,jeh,jwh
      common/ptlphv/th,thh
      data dsq2,vsq2/1.414213562373095d0,0.7071067811865d0/
!     dsq2=sqrt(2.d0)
!     vsq2=1/sqrt(2.d0)
!     dsq3vsq2=sqrt(3.d0)/sqrt(2.d0)
      if(ipae.eq.1) return   !could not  link
      ityae=(ipae-1)/8+1
      imae=mod(ipae-1,8)
      if(imae.eq.0) then
        ityae=ityae-1
        imae=8
      endif
      do 108 ip=1,mh
        iwe=0
        jp=jph(ip)
        lr0=kk(jp)
c      write(*,*)'ip,jpe,ind0',ip,jpe,ind0
        iwa=jwh(ip)
        vlop0=th(ip)
        vlop1=thh(ip)

c      wl5=(vlop0-vlop1)*vo(lr0,lr)-2.0d0*vlop0*vmd(lr0,lr)
c      wl8=vlop0*(vo(lr0,lr0)+(vlop0-1)*0.5*vmd(lr0,lr0))
c             two-index,one-loop 520
c   520=<a,j,k,a>:13,14(ss=3),38(tt=2),50(dd=1)
c       write(nf2,*) 'ityae',ityae
      goto(100,200,300),ityae

!link arc_d
!100     zz='  g50  '
100     wg50 = vlop0*vsq2
        wwg50=-vlop1*sqrt(3.d0)/sqrt(2.d0)
        do la=1,norb_ext
         ma=lsm(la)
          if(ma.ne.imae) cycle
          lra=norb_all-la+1
          iwe=iwe+1
c          wld =(wg50-wwg50)*voint(lra,lr0)-2.d0*wg50*vdint(lra,lr0)
          wld=wg50-wwg50
          call trans_ijkl_intpos(lra,lr0,lra,lr0,nxo)
          call prodel_2(5,wld,jp,iwa,iwe,nxo)
          wld=-2.d0*wg50
          call trans_ijkl_intpos(lra,lra,lr0,lr0,nxo)
          call prodel_2(5,wld,jp,iwa,iwe,nxo)

c       write(*,'(a11,2i3,i6,1x,5f10.4)')
c     :   zz,lr0,la,jwl,vo(lr0,la),vmd(lr0,la),wg50,wwg50,wl
        enddo
        goto 108
!200     zz='  g38,39  '
200     wg38 =-vlop0*vsq2
        wwg38= vlop1
        do ima=1,ng_sm
          imb=mul_tab(ima,imae)
          if(imb.gt.ima) cycle
          do la=ibsm_ext(ima),iesm_ext(ima)
            lra=norb_all-la+1
            lbsta=ibsm_ext(imb)
            lbend=iesm_ext(imb)
            if(ima.eq.imb) lbend=la-1
            do lb=lbsta,lbend
              lrb=norb_all-lb+1
              iwe=iwe+1
c              wlt =(wg38-wwg38)*(voint(lra,lr0)+voint(lrb,lr0))
c     :          -2.d0*wg38*(vdint(lra,lr0)+vdint(lrb,lr0))
              wlt=wg38-wwg38
              call trans_ijkl_intpos(lra,lr0,lra,lr0,nxo)
              call prodel_2(5,wlt,jp,iwa,iwe,nxo)
              wlt=wg38-wwg38
              call trans_ijkl_intpos(lrb,lr0,lrb,lr0,nxo)
              call prodel_2(5,wlt,jp,iwa,iwe,nxo)
              wlt=-2.d0*wg38
              call trans_ijkl_intpos(lra,lra,lr0,lr0,nxo)
              call prodel_2(5,wlt,jp,iwa,iwe,nxo)
              wlt=-2.d0*wg38
              call trans_ijkl_intpos(lrb,lrb,lr0,lr0,nxo)
              call prodel_2(5,wlt,jp,iwa,iwe,nxo)
            enddo
         enddo
        enddo
        goto 108
!300     zz='  g14,15  '
300     wg14 =-vlop0*vsq2
        do ima=1,ng_sm
          imb=mul_tab(ima,imae)
          if(imb.gt.ima) cycle
          if(nlsm_ext(ima).eq.0) cycle
          if(nlsm_ext(imb).eq.0) cycle
         do la=ibsm_ext(ima),iesm_ext(ima)
            lra=norb_all-la+1
            lbsta=ibsm_ext(imb)
            lbend=iesm_ext(imb)
            if(ima.eq.imb) lbend=la-1
            do lb=lbsta,lbend
              lrb=norb_all-lb+1
              iwe=iwe+1
c              wls =wg14*(voint(lra,lr0)+voint(lrb,lr0))
c     :           -2.d0*wg14*(vdint(lra,lr0)+vdint(lrb,lr0))
              wls=wg14
              call trans_ijkl_intpos(lra,lr0,lra,lr0,nxo)
              call prodel_2(5,wls,jp,iwa,iwe,nxo)
              wls=wg14
              call trans_ijkl_intpos(lrb,lr0,lrb,lr0,nxo)
              call prodel_2(5,wls,jp,iwa,iwe,nxo)
              wls=-2.d0*wg14
              call trans_ijkl_intpos(lra,lra,lr0,lr0,nxo)
              call prodel_2(5,wls,jp,iwa,iwe,nxo)
              wls=-2.d0*wg14
              call trans_ijkl_intpos(lrb,lrb,lr0,lr0,nxo)
              call prodel_2(5,wls,jp,iwa,iwe,nxo)
            enddo
         enddo
        enddo
        if(ipae.ne.18) goto 108
!       zz='  g13     '
        wg13 =-vlop0*dsq2
        do la=1,norb_ext
          lra=norb_all-la+1
          iwe=iwe+1
c          wls =wg13*(voint(lra,lr0)-2.d0*vdint(lra,lr0))
          wls=wg13
          call trans_ijkl_intpos(lra,lr0,lra,lr0,nxo)
          call prodel_2(5,wls,jp,iwa,iwe,nxo)
          wls=-2.d0*wg13
          call trans_ijkl_intpos(lra,lra,lr0,lr0,nxo)
          call prodel_2(5,wls,jp,iwa,iwe,nxo)

c        write(*,*)' g13 ',
c     :   vo(lr0,la),vo(lr0,lb),vmd(lr0,la),vmd(lr0,lb)
        enddo
108   continue
      mh=0
      return
      end

      subroutine diagonal_link_ad_g(mpe,iwa,vlop0,vlop1)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
!     common/casrst/ja(max_node),jb(max_node),jm(0:max_node)
!     :    ,jj(4,0:max_node),kk(0:max_node),no(0:max_innorb)
!     :    ,jv,jd(8),jt(8),js(8)
!      common/sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),
!     :     jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)
      integer, pointer :: jph(:),jeh(:),jwh(:)
      real*8, pointer :: th(:),thh(:)
      common/ptlph/jph,jeh,jwh
      common/ptlphv/th,thh
c      common/delm_value/iq,pd,pdd,pt,ptt,ps1,ps2,ps3,ps4
      data dsq2,vsq2/1.414213562373095d0, 0.7071067811865d0/
!     dsq2=sqrt(2.d0)
!     vsq2=1/sqrt(2.d0)
!     dsq3vsq2=sqrt(3.d0)/sqrt(2.d0)
      if(norb_dz.eq.0) return
      ityad=1
      if(jpad.ne.1) ityad=(jpad-1)/8+2
      imad=mod(jpad-1,8)
      if(imad.eq.0) then
        ityad=ityad-1
        imad=8
      endif
      lra=kk(mpe)-1
      goto(100,200,300,400,500,600),ityad
! v: d&r&l(3)
100   fqi=-fg
          vl0=fqi*dsq2*vlop0
          wlv=0.d0
           do lr=1,norb_dz
c           wlv=wlv-vl0*vdint(lr,lra)
            wlv=-vl0
            call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
c            write(nf2,'(4i4,i8)') lra,lr,lra,lr,nxo
            call prodel_2(4,wlv,mpe,0,iwa,nxo)
            wlv=2.0d0*vl0
            call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
            call prodel_2(4,wlv,mpe,0,iwa,nxo)
           enddo
c          call prodel(4,wlv,mpe,0,iwa)
      return
!  jpad=jd(im)
200   fqi=-fg
      do lri=norb_frz+1,norb_dz
        imd=mul_tab(lsm_inn(lri),ns_sm)
        if(imd.ne.imad) cycle
        iwd=jud(lri)

! d: d&r&l(2)
        vl0=fqi*vsq2*vlop0
        vl1=pd*vlop1
c       wld=-2.0d0*vl0*vdint(lra,lri)+(vl0-vl1)*voint(lra,lri)
        wld=-2.0d0*vl0
        call trans_ijkl_intpos(lra,lra,lri,lri,nxo)
        call prodel_2(4,wld,mpe,iwd,iwa,nxo)
        wld=vl0-vl1
        call trans_ijkl_intpos(lra,lri,lra,lri,nxo)
        call prodel_2(4,wld,mpe,iwd,iwa,nxo)

! d: d&r&l(3)+c"(2)
        vl0=fqi*dsq2*vlop0
        do lr=1,lri-1
c          wld=wld+vl0*(voint(lra,lr)-2.d0*vdint(lra,lr))
        wld=-2.0d0*vl0
        call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
        call prodel_2(4,wld,mpe,iwd,iwa,nxo)
        wld=vl0
        call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
        call prodel_2(4,wld,mpe,iwd,iwa,nxo)
        enddo
! d: d&r&l(3)
        do lr=lri+1,norb_dz
c          wld=wld+vl0*(voint(lra,lr)-2.d0*vdint(lra,lr))
        wld=-2.0d0*vl0
        call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
        call prodel_2(4,wld,mpe,iwd,iwa,nxo)
        wld=vl0
        call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
        call prodel_2(4,wld,mpe,iwd,iwa,nxo)
        enddo

c       call prodel(4,wld,mpe,iwd,iwa)
      enddo
      return
!  jpad=jt(im)
300   fqi=fg
      iwt=0
      do lri=norb_frz+1,norb_dz
        imi=mul_tab(lsm_inn(lri),ns_sm)
        do lrj=lri+1,norb_dz
          imj=lsm_inn(lrj)
          imij=mul_tab(imi,imj)
          if(imij.ne.imad) cycle
          iwt=just(lri,lrj)    !

! t: d&r&l(2)
! t: d&r&l(2)+c"(2)
          vl0=fqi*vsq2*vlop0
          vl1=pt*vlop1
c          wlt=   -2.0d0*vl0*vdint(lra,lri)+(vl0-vl1)*voint(lra,lri)
          wlt=-2.0d0*vl0
          call trans_ijkl_intpos(lra,lra,lri,lri,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
          wlt=vl0-vl1
          call trans_ijkl_intpos(lra,lri,lra,lri,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)


c          wlt=wlt-2.0d0*vl0*vdint(lra,lrj)+(vl0-vl1)*voint(lra,lrj)
          wlt=-2.0d0*vl0
          call trans_ijkl_intpos(lra,lra,lrj,lrj,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
          wlt=vl0-vl1
          call trans_ijkl_intpos(lra,lrj,lra,lrj,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)


          vl0=fqi*dsq2*vlop0
          do lr=1,lri-1
! t: d&r&l(3)+c"(2)+c"(2)
c            wlt=wlt+vl0*vdint(lr,lra)
          wlt=vl0
          call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
          wlt=-2.0d0*vl0
          call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
          enddo
! t: d&r&l(3)+c"(2)
          do lr=lri+1,lrj-1
c            wlt=wlt+vl0*vdint(lr,lra)
          wlt=vl0
          call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
          wlt=-2.0d0*vl0
          call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
          enddo
! t: d&r&l(3)
          do lr=lrj+1,norb_dz
c            wlt=wlt+vl0*vdint(lr,lra)
          wlt=vl0
          call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
          wlt=-2.0d0*vl0
          call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
         enddo
c         call prodel(4,wlt,mpe,iwt,iwa)
        enddo
      enddo
      return
400   fqi=fg
      iws=0
      do lri=norb_frz+1,norb_dz
        if(imad.eq.ns_sm)  then
          lrj=lri
          iws=just(lri,lri)
! s: d&r&l(3)+c"(0)
! s: d&r&l(3)
        vl0=fqi*dsq2*vlop0
c         wls=0.d0
          do lr=1,norb_dz
            if(lr.eq.lri) cycle
c           wls=wls+vl0*(voint(lra,lr)-2.d0*vdint(lra,lr))
          wls=vl0
          call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
          wls=-2.0d0*vl0
          call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
          enddo
c          call prodel(4,wls,mpe,iws,iwa)
        endif

        imi=mul_tab(lsm_inn(lri),ns_sm)
        lrjsta=lri+1
        do lrj=lrjsta,norb_dz
          imj=lsm_inn(lrj)
          imij=mul_tab(imi,imj)
          if(imij.ne.imad) cycle
          iws=just(lri,lrj)
! s1: d&r&l(1)
          vl0=fqi*vsq2*vlop0
          vl1=ps1*vlop1
c          wls=-2.0d0*vl0*vdint(lra,lrj)+(vl0-vl1)*voint(lra,lrj)
          wls=-2.0d0*vl0
          call trans_ijkl_intpos(lra,lra,lrj,lrj,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
          wls=vl0-vl1
          call trans_ijkl_intpos(lra,lrj,lra,lrj,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)

! s4: d&r&l(2)+c"(1)
          vl0=fqi*vsq2*vlop0
          vl1=ps4*vlop1
c          wls=wls-2.0d0*vl0*vdint(lra,lri)+(vl0-vl1)*voint(lra,lri)
          wls=-2.0d0*vl0
          call trans_ijkl_intpos(lra,lra,lri,lri,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
          wls=vl0-vl1
          call trans_ijkl_intpos(lra,lri,lra,lri,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)

          vl0=fqi*dsq2*vlop0
! s: d&r&l(3)+c"(2)+c"(1)
          do lr=1,lri-1
c            wls=wls+vl0*vdint(lr,lra)
          wls=vl0
          call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
          wls=-2.0d0*vl0
          call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
          enddo
! s: d&r&l(3)+c"(1)
          do lr=lri+1,lrj-1
c            wls=wls+vl0*vdint(lr,lra)
          wls=vl0
          call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
          wls=-2.0d0*vl0
          call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
          enddo
! s: d&r&l(3)
          do lr=lrj+1,norb_dz
c            wls=wls+vl0*vdint(lr,lra)
          wls=vl0
          call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
          wls=-2.0d0*vl0
          call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
         enddo
c          call prodel(4,wls,mpe,iws,iwa)
        enddo
      enddo

      if(jb_sys.eq.0) return      !any diffrence when jb_sys=1 and jb_sy
      fqi=fg
      do lri=norb_frz+1,norb_dz
        imi=mul_tab(lsm_inn(lri),ns_sm)
        do lrj=lri+1,norb_dz
          imj=lsm_inn(lrj)
          imij=mul_tab(imi,imj)
          if(imij.ne.imad) cycle
       !   iws=iws+1
          iws=just(lrj,lri)
! s1: d&r&l(1)-c"(2)
          vl0=fqi*vsq2*vlop0
          vl1=ps3*vlop1
c          wls=-2.0d0*vl0*vdint(lra,lri)+(vl0-vl1)*voint(lra,lri)
          wls=-2.0d0*vl0
          call trans_ijkl_intpos(lra,lra,lri,lri,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
          wls=vl0-vl1
          call trans_ijkl_intpos(lra,lri,lra,lri,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)

! s3: (11)d&r&l(2)
          vl0=fqi*vsq2*vlop0
          vl1=ps2*vlop1
c          wls=wls-2.0d0*vl0*vdint(lra,lrj)+(vl0-vl1)*voint(lra,lrj)
          wls=-2.0d0*vl0
          call trans_ijkl_intpos(lra,lra,lrj,lrj,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
          wls=vl0-vl1
          call trans_ijkl_intpos(lra,lrj,lra,lrj,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)

          vl0=fqi*dsq2*vlop0
! s: d&r&l(3)+c"(1)+c"(2)
          do lr=1,lri-1
c            wls=wls+vl0*vdint(lr,lra)
          wls=vl0
          call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
          wls=-2.0d0*vl0
          call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
          enddo
! s: d&r&l(3)+c"(2)
          do lr=lri+1,lrj-1
c            wls=wls+vl0*vdint(lr,lra)
          wls=vl0
          call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
          wls=-2.0d0*vl0
          call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
          enddo
! s: d&r&l(3)
          do lr=lrj+1,norb_dz
c            wls=wls+vl0*vdint(lr,lra)
          wls=vl0
          call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
          wls=-2.0d0*vl0
          call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
           enddo
c          call prodel(4,wls,mpe,iws,iwa)
        enddo
      enddo
      return

500   fqi=-fg
      iwd=0
      do lri=norb_frz+1,norb_dz
        imd=mul_tab(lsm_inn(lri),ns_sm)
        if(imd.ne.imad) cycle
        iwd=jud(lri)

! dd1: d&r&l(1)
        vl0=fqi*vsq2*vlop0
        vl1=pdd*vlop1
c        wld=-2.0d0*vl0*vdint(lra,lri)+(vl0-vl1)*voint(lra,lri)
          wld=vl0-vl1
          call trans_ijkl_intpos(lra,lri,lra,lri,nxo)
          call prodel_2(4,wld,mpe,iwd,iwa,nxo)
          wld=-2.0d0*vl0
          call trans_ijkl_intpos(lra,lra,lri,lri,nxo)
          call prodel_2(4,wld,mpe,iwd,iwa,nxo)

        vl0=fqi*dsq2*vlop0
! d: d&r&l(3)+c"(1)
        do lr=1,lri-1
c          wld=wld+vl0*vdint(lr,lra)
          wld=vl0
          call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
          call prodel_2(4,wld,mpe,iwd,iwa,nxo)
          wld=-2.0d0*vl0
          call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
          call prodel_2(4,wld,mpe,iwd,iwa,nxo)
        enddo
! d: d&r&l(3)
        do lr=lri+1,norb_dz
c          wld=wld+vl0*vdint(lr,lra)
          wld=vl0
          call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
          call prodel_2(4,wld,mpe,iwd,iwa,nxo)
          wld=-2.0d0*vl0
          call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
          call prodel_2(4,wld,mpe,iwd,iwa,nxo)
        enddo

c       call prodel(4,wld,mpe,iwd,iwa)
      enddo
      return

600   fqi=fg
      iwt=0
      do lri=norb_frz+1,norb_dz
        imi=mul_tab(lsm_inn(lri),ns_sm)
        do lrj=lri+1,norb_dz
          imj=lsm_inn(lrj)
          imij=mul_tab(imi,imj)
          if(imij.ne.imad) cycle
          iwt=just(lri,lrj)

! tt: d&r&l(1)
! tt: d&r&l(1)+c"(1)
          vl0=fqi*vsq2*vlop0
          vl1=ptt*vlop1
c          wlt=   -2.0d0*vl0*vdint(lra,lri)+(vl0-vl1)*voint(lra,lri)
          wlt=-2.0d0*vl0
          call trans_ijkl_intpos(lra,lra,lri,lri,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
          wlt=vl0-vl1
          call trans_ijkl_intpos(lra,lri,lra,lri,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)

c          wlt=wlt-2.0d0*vl0*vdint(lra,lrj)+(vl0-vl1)*voint(lra,lrj)
          wlt=-2.0d0*vl0
          call trans_ijkl_intpos(lra,lra,lrj,lrj,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
          wlt=vl0-vl1
          call trans_ijkl_intpos(lra,lrj,lra,lrj,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)

          vl0=fqi*dsq2*vlop0
          do lr=1,lri-1
! t: d&r&l(3)+c"(1)+c"(1)
c            wlt=wlt+vl0*vdint(lr,lra)
          wlt=vl0
          call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
          wlt=-2.0d0*vl0
          call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
          enddo
! t: d&r&l(3)+c"(1)
          do lr=lri+1,lrj-1
c            wlt=wlt+vl0*vdint(lr,lra)
          wlt=vl0
          call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
          wlt=-2.0d0*vl0
          call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
          enddo
! t: d&r&l(3)
          do lr=lrj+1,norb_dz
c            wlt=wlt+vl0*vdint(lr,lra)
          wlt=vl0
          call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
          wlt=-2.0d0*vl0
          call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
         enddo

c         call prodel(4,wlt,mpe,iwt,iwa)
        enddo
      enddo
      return
      end


      subroutine diagonal_link_dae_g(mh)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
!     common/casrst/ja(max_node),jb(max_node),jm(0:max_node)
!     :    ,jj(4,0:max_node),kk(0:max_node),no(0:max_innorb)
!     :    ,jv,jd(8),jt(8),js(8)
!      common/sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),
!     :     jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)
      integer, pointer :: jph(:),jeh(:),jwh(:)
      real*8, pointer :: th(:),thh(:)
      common/ptlph/jph,jeh,jwh
      common/ptlphv/th,thh
      data dsq2,vsq2/1.414213562373095d0, 0.7071067811865d0/
!     dsq2=sqrt(2.d0)
!     vsq2=1/sqrt(2.d0)
!     dsq3vsq2=sqrt(3.d0)/sqrt(2.d0)
      ityad=1
      if(jpad.ne.1) ityad= (jpad-1)/8+2
      imad=mod(jpad-1,8)
      if(imad.eq.0) then
        ityad=ityad-1
        imad=8
      endif
      do 108 ip=1,mh
        iwa=jwh(ip)
        vlop0=th(ip)
        vlop1=thh(ip)

      goto(100,200,300,400,500,600),ityad

!  jpad=1
100   if(abs(vlop0).lt.1e-30) goto 108
        fqi=fg
! v: d&r&l(3)
        lri=0
        lrj=0
        iwd=0
        vij0=0.d0
        vij1=0.d0
        vij2=0.d0
        vl0=fqi*dsq2*vlop0
        call diagonal_call_dae_g(lri,lrj,iwd,iwa,vij0,vij1,vij2,vl0)
        goto 108
!  jpad=jd(im)
200   fqi=-fg
      lrj=0
      do lri=norb_frz+1,norb_dz
        imd=mul_tab(lsm_inn(lri),ns_sm)
        if(imd.ne.imad) cycle
        iwd=jud(lri)

! d: d&r&l(2)
        vij0=fqi*vsq2*vlop0
        vij1=pd*vlop1
        vij2=0.d0
! d: d&r&l(3)+c"(2)
! d: d&r&l(3)
        vl0=fqi*dsq2*vlop0

        call diagonal_call_dae_g(lri,lrj,iwd,iwa,vij0,vij1,vij2,vl0)

      enddo
      goto 108
!  jpad=jt(im)
300   fqi=fg
      iwt=0
      do lri=norb_frz+1,norb_dz
        imi=mul_tab(lsm_inn(lri),ns_sm)
        do lrj=lri+1,norb_dz
          imj=lsm_inn(lrj)
          imij=mul_tab(imi,imj)
          if(imij.ne.imad) cycle
          iwt=just(lri,lrj)      !

! t: d&r&l(2)
! t: d&r&l(2)+c"(2)
          vij0=fqi*vsq2*vlop0
          vij1=pt*vlop1
          vij2=vij1
! t: d&r&l(3)+c"(2)+c"(2)
! t: d&r&l(3)+c"(2)
! t: d&r&l(3)
          vl0=fqi*dsq2*vlop0

          call diagonal_call_dae_g(lri,lrj,iwt,iwa,vij0,vij1,vij2,vl0)
        enddo
      enddo
      goto 108
400   fqi=fg
      iws=0
      do lri=norb_frz+1,norb_dz           !cc
        if(imad.eq.ns_sm) then
          lrj=lri
         iws=just(lri,lri)
         vij0=0.d0
          vij1=0.d0
          vij2=0.d0
! s: d&r&l(3)+c"(0)
! s: d&r&l(3)
          vl0=fqi*dsq2*vlop0
          call diagonal_call_dae_g(lri,lrj,iws,iwa,vij0,vij1,vij2,vl0)
        endif

        imi=mul_tab(lsm_inn(lri),ns_sm)
        lrjsta=lri+1
        do lrj=lrjsta,norb_dz
          imj=lsm_inn(lrj)
          imij=mul_tab(imi,imj)
          if(imij.ne.imad) cycle
          iws=just(lri,lrj)
! s2: d&r&l(2)
! s4: d&r&l(2)+c"(1)
          vij0=fqi*vsq2*vlop0
          vij1=ps1*vlop1
          vij2=ps4*vlop1
! s: d&r&l(3)+c"(2)+c"(1)
! s: d&r&l(3)+c"(1)
! s: d&r&l(3)
          vl0=fqi*dsq2*vlop0
          call diagonal_call_dae_g(lri,lrj,iws,iwa,vij0,vij1,vij2,vl0)
        enddo
      enddo
      if(jb_sys.eq.0) goto 108
      fqi=fg
      do lri=norb_frz+1,norb_dz
        imi=mul_tab(lsm_inn(lri),ns_sm)
        if(imad.eq.ns_sm) lrjsta=lri
        do lrj=lri+1,norb_dz
          imj=lsm_inn(lrj)
          imij=mul_tab(imi,imj)
          if(imij.ne.imad) cycle
          iws=just(lrj,lri)
! s1: d&r&l(1)
! s3: d&r&l(1)+c"(2)
          vij0=fqi*vsq2*vlop0
          vij1=ps2*vlop1
          vij2=ps3*vlop1
! s: d&r&l(3)+c"(1)+c"(2)
! s: d&r&l(3)+c"(2)
! s: d&r&l(3)
          vl0=fqi*dsq2*vlop0

          call diagonal_call_dae_g(lri,lrj,iws,iwa,vij0,vij1,vij2,vl0)

        enddo
      enddo
      goto 108
500   fqi=-fg
      lrj=0
      do lri=norb_frz+1,norb_dz
        imd=mul_tab(lsm_inn(lri),ns_sm)
        if(imd.ne.imad) cycle
        iwd=jud(lri)

! dd1: d&r&l(1)
        vij0=fqi*vsq2*vlop0
        vij1=pdd*vlop1
        vij2=0.d0
! d: d&r&l(3)+c"(1)
! d: d&r&l(3)
        vl0=fqi*dsq2*vlop0

        call diagonal_call_dae_g(lri,lrj,iwd,iwa,vij0,vij1,vij2,vl0)

      enddo
      goto 108

600   fqi=fg               !aa
      iwt=0
      do lri=norb_frz+1,norb_dz
        imi=mul_tab(lsm_inn(lri),ns_sm)
        do lrj=lri+1,norb_dz
          imj=lsm_inn(lrj)
          imij=mul_tab(imi,imj)
          if(imij.ne.imad) cycle
          iwt=just(lri,lrj)

          vij0=fqi*vsq2*vlop0
! tt: d&r&l(1)
          vij1=ptt*vlop1
! tt: d&r&l(1)+c"(1)
          vij2=ptt*vlop1
! t: d&r&l(3)+c"(1)+c"(1)
! t: d&r&l(3)+c"(1)
! t: d&r&l(3)
          vl0=fqi*dsq2*vlop0
          call diagonal_call_dae_g(lri,lrj,iwt,iwa,vij0,vij1,vij2,vl0)
        enddo
      enddo

108   continue
      return
      end

      subroutine diagonal_call_dae_g(lri,lrj,iwd,iwa,vij0,vij1,vij2,vl0)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
!     common/casrst/ja(max_node),jb(max_node),jm(0:max_node)
!     :    ,jj(4,0:max_node),kk(0:max_node),no(0:max_innorb)
!     :    ,jv,jd(8),jt(8),js(8)
!      common/sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),
!     :     jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)
      integer, pointer :: jph(:),jeh(:),jwh(:)
      real*8, pointer :: th(:),thh(:)
      common/ptlph/jph,jeh,jwh
      common/ptlphv/th,thh
      data dsq2,vsq2/1.414213562373095d0, 0.7071067811865d0 /
      data dsq3vsq2/1.224744871392d0/
!     dsq2=sqrt(2.d0)
!     vsq2=1/sqrt(2.d0)
!     dsq3vsq2=sqrt(3.d0)/sqrt(2.d0)

      if(norb_dz.eq.0) return
      if(ipae.eq.1) return   !could not  link
      ityae=(ipae-1)/8+1
      imae=mod(ipae-1,8)
      if(imae.eq.0) then
        ityae=ityae-1
        imae=8
      endif
      iwe=0
c   520=<a,j,k,a>:13,14(ss=3),38(tt=2),50(dd=1)
      goto(100,200,300),ityae
!link arc_d
!100     zz='  g50  '

100   do la=ibsm_ext(imae),iesm_ext(imae)
          lra=norb_all-la+1
         iwe=iwe+1
         vlop0=vl0*vsq2
c        wl=0.d0
         do lr=1,norb_dz
            if(lr.eq.lri) cycle
            if(lr.eq.lrj) cycle
c            wl =wl+vlop0*vdint(lr,lra)    !db space drl(33)- ext space
            wl=vlop0
            call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            wl=-2.0d0*vlop0
            call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
         enddo
      if(lri.ge.lrj) then
         vlop0= vij0*vsq2
          vlop1=-vij1*dsq3vsq2
c         wl=wl+(vlop0-vlop1)*voint(lra,lri)-2.d0*vlop0*vdint(lra,lri)
          if(lri.ne.0) then
            wl=vlop0-vlop1
            call trans_ijkl_intpos(lra,lri,lra,lri,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            wl=-2.0d0*vlop0
            call trans_ijkl_intpos(lra,lra,lri,lri,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
          endif
          vlop1=-vij2*dsq3vsq2
c         wl=wl+(vlop0-vlop1)*voint(lra,lrj)-2.d0*vlop0*vdint(lra,lrj)
          if(lrj.ne.0) then
            wl=vlop0-vlop1
            call trans_ijkl_intpos(lra,lrj,lra,lrj,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            wl=-2.0d0*vlop0
            call trans_ijkl_intpos(lra,lra,lrj,lrj,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
          endif
c       write(*,'(a11,2i3,i6,1x,5f10.4)')
c     :   zz,lr0,la,jwl,vo(lr0,la),vmd(lr0,la),wg50,wwg50,wl
      else
          vlop0= vij0*vsq2
          vlop1=-vij1*dsq3vsq2       !db space (22)drl(11)- ext space -g
c          wl=wl+(vlop0-vlop1)*voint(lra,lrj)-2.d0*vlop0*vdint(lra,lrj)
          if(lrj.ne.0) then
            wl=vlop0-vlop1
            call trans_ijkl_intpos(lra,lrj,lra,lrj,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            wl=-2.0d0*vlop0
            call trans_ijkl_intpos(lra,lra,lrj,lrj,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
          endif
          vlop1=-vij2*dsq3vsq2        !db space drl(22)c"(11)- ext space
c          wl=wl+(vlop0-vlop1)*voint(lra,lri)-2.d0*vlop0*vdint(lra,lri)
          if(lri.ne.0) then
            wl=vlop0-vlop1
            call trans_ijkl_intpos(lra,lri,lra,lri,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            wl=-2.0d0*vlop0
            call trans_ijkl_intpos(lra,lra,lri,lri,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
          endif
      endif
c         call prodel(6,wl,iwd,iwa,iwe)
        enddo
        goto 108
!200     zz='  g38,39  '
200   do ima=1,ng_sm
        imb=mul_tab(ima,imae)
        if(imb.gt.ima) cycle
        do la=ibsm_ext(ima),iesm_ext(ima)
          lra=norb_all-la+1
          lbsta=ibsm_ext(imb)
          lbend=iesm_ext(imb)
          if(ima.eq.imb) lbend=la-1
          do lb=lbsta,lbend
            lrb=norb_all-lb+1
            iwe=iwe+1
c            volalb=0.d0
c            vd2lalb=0.d0
            vlop0=-vl0*vsq2

           do lr=1,norb_dz
              if(lr.eq.lri) cycle
              if(lr.eq.lrj) cycle
c            volalb=volalb+(voint(lra,lr)+voint(lrb,lr))
c            vd2lalb=vd2lalb-2.d0*(vdint(lra,lr)+vdint(lrb,lr))
            wl=vlop0
            call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            call trans_ijkl_intpos(lrb,lr,lrb,lr,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            wl=-2.0d0*vlop0
            call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            call trans_ijkl_intpos(lrb,lrb,lr,lr,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)

            enddo

c           wl=vlop0*(volalb+vd2lalb)
      if(lri.ge.lrj) then
            vlop0=-vij0*vsq2
            vlop1=vij1
c           wl=wl+(vlop0-vlop1)*(voint(lra,lri)+voint(lrb,lri))
c     :          -2.d0*vlop0*(vdint(lra,lri)+vdint(lrb,lri))
        if(lri.ne.0) then
            wl=vlop0-vlop1
            call trans_ijkl_intpos(lra,lri,lra,lri,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            call trans_ijkl_intpos(lrb,lri,lrb,lri,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            wl=-2.0d0*vlop0
            call trans_ijkl_intpos(lra,lra,lri,lri,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            call trans_ijkl_intpos(lrb,lrb,lri,lri,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
        endif
            vlop1=vij2
c           wl=wl+(vlop0-vlop1)*(voint(lra,lrj)+voint(lrb,lrj))
c     :          -2.d0*vlop0*(vdint(lra,lrj)+vdint(lrb,lrj))
        if(lrj.ne.0) then
            wl=vlop0-vlop1
            call trans_ijkl_intpos(lra,lrj,lra,lrj,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            call trans_ijkl_intpos(lrb,lrj,lrb,lrj,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            wl=-2.0d0*vlop0
            call trans_ijkl_intpos(lra,lra,lrj,lrj,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            call trans_ijkl_intpos(lrb,lrb,lrj,lrj,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
        endif
      else
            vlop0=-vij0*vsq2
            vlop1=vij1
c           wl=wl+(vlop0-vlop1)*(voint(lra,lrj)+voint(lrb,lrj))
c     :          -2.d0*vlop0*(vdint(lra,lrj)+vdint(lrb,lrj))
        if(lrj.ne.0) then
            wl=vlop0-vlop1
            call trans_ijkl_intpos(lra,lrj,lra,lrj,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            call trans_ijkl_intpos(lrb,lrj,lrb,lrj,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            wl=-2.0d0*vlop0
            call trans_ijkl_intpos(lra,lra,lrj,lrj,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            call trans_ijkl_intpos(lrb,lrb,lrj,lrj,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
        endif
            vlop1=vij2
c           wl=wl+(vlop0-vlop1)*(voint(lra,lri)+voint(lrb,lri))
c     :          -2.d0*vlop0*(vdint(lra,lri)+vdint(lrb,lri))
        if(lri.ne.0) then
            wl=vlop0-vlop1
            call trans_ijkl_intpos(lra,lri,lra,lri,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            call trans_ijkl_intpos(lrb,lri,lrb,lri,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            wl=-2.0d0*vlop0
            call trans_ijkl_intpos(lra,lra,lri,lri,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            call trans_ijkl_intpos(lrb,lrb,lri,lri,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
        endif
      endif
c              call prodel(6,wl,iwd,iwa,iwe)
            enddo
         enddo
        enddo
        goto 108
!300     zz='  g14,15  '
c===========================lyb=======
c300   do ima=1,8
c      the 8 should be changed to ng_sm

300   do ima=1,ng_sm
        imb=mul_tab(ima,imae)
        if(imb.gt.ima) cycle
        do la=ibsm_ext(ima),iesm_ext(ima)
          lra=norb_all-la+1
          lbsta=ibsm_ext(imb)
          lbend=iesm_ext(imb)
          if(ima.eq.imb) lbend=la-1
          do lb=lbsta,lbend
            lrb=norb_all-lb+1
            iwe=iwe+1
c            volalb=0.d0
c            vd2lalb=0.d0

            wg14 =-vl0*vsq2

           do lr=1,norb_dz
              if(lr.eq.lri) cycle
              if(lr.eq.lrj) cycle
c            volalb=volalb+(voint(lra,lr)+voint(lrb,lr))
c            vd2lalb=vd2lalb-2.d0*(vdint(lra,lr)+vdint(lrb,lr))
            wl=wg14
            call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            call trans_ijkl_intpos(lrb,lr,lrb,lr,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            wl=-2.0d0*wg14
            call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            call trans_ijkl_intpos(lrb,lrb,lr,lr,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            enddo

c           wl=wg14*(volalb+vd2lalb)

            if(jpad.eq.18.and.lri.eq.lrj) goto 301
            wg14 =-vij0*vsq2
c           wl=wl+wg14*(voint(lra,lri)+voint(lrb,lri))
c     :          -2.d0*wg14*(vdint(lra,lri)+vdint(lrb,lri))
        if(lri.ne.0) then
            wl=wg14
            call trans_ijkl_intpos(lra,lri,lra,lri,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            call trans_ijkl_intpos(lrb,lri,lrb,lri,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            wl=-2.0d0*wg14
            call trans_ijkl_intpos(lra,lra,lri,lri,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            call trans_ijkl_intpos(lrb,lrb,lri,lri,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
        endif
c           wl=wl+wg14*(voint(lra,lrj)+voint(lrb,lrj))
c     :          -2.d0*wg14*(vdint(lra,lrj)+vdint(lrb,lrj))
        if(lrj.ne.0) then
            wl=wg14
            call trans_ijkl_intpos(lra,lrj,lra,lrj,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            call trans_ijkl_intpos(lrb,lrj,lrb,lrj,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            wl=-2.0d0*wg14
            call trans_ijkl_intpos(lra,lra,lrj,lrj,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            call trans_ijkl_intpos(lrb,lrb,lrj,lrj,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
        endif
c301
c         call prodel(6,wl,iwd,iwa,iwe)
301        continue
        enddo
         enddo
        enddo

        if(ipae.ne.18) return
!       zz='  g13     '
        do la=1,norb_ext
           lra=norb_all-la+1

           iwe=iwe+1
c          vovdla=0.d0

           wg13 =-vl0*dsq2

          do lr=1,norb_dz
            if(lr.eq.lri) cycle
            if(lr.eq.lrj) cycle
c           vovdla=vovdla+vdint(lr,lra)
            wl=wg13
            call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            wl=-2.0d0*wg13
            call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
          enddo

c           wl=wg13*vovdla

            wg13 =-vij0*dsq2
c           wl=wl+wg13*(vdint(lri,lra)+vdint(lrj,lra))
        if(lri.ne.0) then
            wl=wg13
            call trans_ijkl_intpos(lra,lri,lra,lri,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            wl=-2.0d0*wg13
            call trans_ijkl_intpos(lra,lra,lri,lri,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
        endif
        if(lrj.ne.0) then
            wl=wg13
            call trans_ijkl_intpos(lra,lrj,lra,lrj,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            wl=-2.0d0*wg13
            call trans_ijkl_intpos(lra,lra,lrj,lrj,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
         endif
c            call prodel(6,wl,iwd,iwa,iwe)
        enddo
108   continue
      return
      end

      subroutine diagonal_dbl_g()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
!     common/casrst/ja(max_node),jb(max_node),jm(0:max_node)
!     :    ,jj(4,0:max_node),kk(0:max_node),no(0:max_innorb)
!     :    ,jv,jd(8),jt(8),js(8)
!      common/sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),
!     :     jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)
!     data dzero/0.d0/
      logical logic_lij

      if(norb_dz.eq.0) return
c        wt0=dzero
      do lr=1,norb_dz                   ! ........800...
c        wt0=wt0+voint(lr,lr)+voint(lr,lr)+vdint(lr,lr)
         wt0_1=2.0d0
         wt0_2=1.0d0
         call trans_ijkl_intpos(lr,lr,lr,lr,nxo)
        do ipae_=1,25
          ipae=ipae_ ! ipae is in common block, is this necessary?
          if(nu_ae(ipae).eq.0) cycle
          iwdownv=iw_downwei(1,ipae)
          do iwa=0,iwdownv-1
c       zz=' doub_800_v'
             iwad=iwalk_ad(1,ipae,iwa,0)
             call prodel_1(1,wt0_1,0,ipae,iwad,lr,lr)
             call prodel_2(1,wt0_2,0,ipae,iwad,nxo)
c         zz=' doub_800_s'
          enddo
        enddo
      enddo

c       jps=js(1)

      mr0=1
      do 100 lr0=norb_frz+1,norb_dz
      do lrd=1,norb_dz
         mr0=mul_tab(lsm_inn(lr0),ns_sm)
         iwd=jud(lr0)
         jpad =1+mr0
         jpad1=jpad+24

         call trans_ijkl_intpos(lrd,lrd,lrd,lrd,nxo)

        if(lr0.eq.lrd) then
! ......d_800.
c         wld=wt0-voint(lr0,lr0)-vdint(lr0,lr0)
c          wls=wld-voint(lr0,lr0)
           wld_1=1.0d0
          do ipae_=1,25
            ipae=ipae_ ! ipae is in common block, is this necessary?
            if(nu_ae(ipae).eq.0) cycle
            iwdownv=iw_downwei(jpad,ipae)
            do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpad,ipae,iwa,iwd)
             call prodel_1(1,wld_1,0,ipae,iwad,lrd,lrd)
c             call prodel(1,wld,0,ipae,iwad)
            enddo
          enddo

          if(jb_sys.gt.0) then
            do ipae_=1,25
            ipae=ipae_ ! ipae is in common block, is this necessary?
            if(nu_ae(ipae).eq.0) cycle
            iwdownv=iw_downwei(jpad1,ipae)
            do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpad1,ipae,iwa,iwd)
             call prodel_1(1,wld_1,0,ipae,iwad,lrd,lrd)
c             call prodel(1,wld,0,ipae,iwad)
            enddo
            enddo
          endif

        else
           wld_1=2.0d0
           wld_2=1.0d0
           wls_1=2.0d0
           wls_2=1.0d0

          do ipae_=1,25
            ipae=ipae_ ! ipae is in common block, is this necessary?
            if(nu_ae(ipae).eq.0) cycle
            iwdownv=iw_downwei(jpad,ipae)
            do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpad,ipae,iwa,iwd)
             call prodel_1(1,wld_1,0,ipae,iwad,lrd,lrd)
             call prodel_2(1,wld_2,0,ipae,iwad,nxo)
c             call prodel(1,wld,0,ipae,iwad)
            enddo
          enddo

          if(jb_sys.gt.0) then
            do ipae_=1,25
            ipae=ipae_ ! ipae is in common block, is this necessary?
            if(nu_ae(ipae).eq.0) cycle
            iwdownv=iw_downwei(jpad1,ipae)
            do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpad1,ipae,iwa,iwd)
             call prodel_1(1,wld_1,0,ipae,iwad,lrd,lrd)
             call prodel_2(1,wld_2,0,ipae,iwad,nxo)
c             call prodel(1,wld,0,ipae,iwad)
            enddo
            enddo
          endif
! 0.....s_800.
          jpad=17+ns_sm
          iwd=just(lr0,lr0)
          do ipae_=1,25
          ipae=ipae_ ! ipae is in common block, is this necessary?
          if(nu_ae(ipae).eq.0) cycle
          iwdownv=iw_downwei(jpad,ipae)
          do iwa=0,iwdownv-1
            iwad=iwalk_ad(jpad,ipae,iwa,iwd)
             call prodel_1(1,wls_1,0,ipae,iwad,lrd,lrd)
             call prodel_2(1,wls_2,0,ipae,iwad,nxo)
c         call prodel(1,wls,0,ipae,iwad)
          enddo
          enddo
        endif
      enddo
      if(lr0.eq.norb_dz) goto 100
c        wld0=wld
! ........800...
      do 200 lr=lr0+1,norb_dz
        mr=mul_tab(mr0,lsm_inn(lr))
        jpat=9+mr
        jpas=17+mr
        jpat1=jpat+24
        iws=just(lr0,lr)
        iwt=iws            !
c          wld=wld0-voint(lr,lr)-vdint(lr,lr)
c=============
c      lr0 and lr

         wld_1=1.0d0
         wld_2=0.0d0
         call trans_ijkl_intpos(lr0,lr0,lr0,lr0,nxo)
         nxo_1=nxo
         call trans_ijkl_intpos(lr,lr,lr,lr,nxo)
         nxo_2=nxo

         do ipae_=1,25
            ipae=ipae_ ! ipae is in common block, is this necessary?
            if(nu_ae(ipae).eq.0) cycle
            iwdownv=iw_downwei(jpat,ipae)
            do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpat,ipae,iwa,iwt)
             call prodel_1(1,wld_1,0,ipae,iwad,lr0,lr0)
             call prodel_2(1,wld_2,0,ipae,iwad,nxo_1)
             call prodel_1(1,wld_1,0,ipae,iwad,lr,lr)
             call prodel_2(1,wld_2,0,ipae,iwad,nxo_2)
c           call prodel(1,wld,0,ipae,iwad)
            enddo
            if(jb_sys.gt.1) then
              iwdownv=iw_downwei(jpat1,ipae)
            do iwa=0,iwdownv-1
                iwad=iwalk_ad(jpat1,ipae,iwa,iwt)         !t1
             call prodel_1(1,wld_1,0,ipae,iwad,lr0,lr0)
             call prodel_2(1,wld_2,0,ipae,iwad,nxo_1)
             call prodel_1(1,wld_1,0,ipae,iwad,lr,lr)
             call prodel_2(1,wld_2,0,ipae,iwad,nxo_2)
c             call prodel(1,wld,0,ipae,iwad)
            enddo
            endif
            iwdownv=iw_downwei(jpas,ipae)
            do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpas,ipae,iwa,iws)
             call prodel_1(1,wld_1,0,ipae,iwad,lr0,lr0)
             call prodel_2(1,wld_2,0,ipae,iwad,nxo_1)
             call prodel_1(1,wld_1,0,ipae,iwad,lr,lr)
             call prodel_2(1,wld_2,0,ipae,iwad,nxo_2)
c           call prodel(1,wld,0,ipae,iwad)
            enddo
            if(jb_sys.gt.0) then
            iws1=just(lr,lr0)
            do iwa=0,iwdownv-1
                iwad=iwalk_ad(jpas,ipae,iwa,iws1)
             call prodel_1(1,wld_1,0,ipae,iwad,lr0,lr0)
             call prodel_2(1,wld_2,0,ipae,iwad,nxo_1)
             call prodel_1(1,wld_1,0,ipae,iwad,lr,lr)
             call prodel_2(1,wld_2,0,ipae,iwad,nxo_2)
c             call prodel(1,wld,0,ipae,iwad)
            enddo
            endif
         enddo

c==============
c     1 <= lrd < norb_dz and lrd.ne.lr0.and.lrd.ne.lr
      do lrd=1,norb_dz

       if(lrd.ne.lr0.and.lrd.ne.lr) then
         wld_1=2.0d0
         wld_2=1.0d0
         call trans_ijkl_intpos(lrd,lrd,lrd,lrd,nxo)
          do ipae_=1,25
            ipae=ipae_ ! ipae is in common block, is this necessary?
            if(nu_ae(ipae).eq.0) cycle
            iwdownv=iw_downwei(jpat,ipae)
            do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpat,ipae,iwa,iwt)
             call prodel_1(1,wld_1,0,ipae,iwad,lrd,lrd)
             call prodel_2(1,wld_2,0,ipae,iwad,nxo)
c           call prodel(1,wld,0,ipae,iwad)
          enddo
            if(jb_sys.gt.1) then
              iwdownv=iw_downwei(jpat1,ipae)
              do iwa=0,iwdownv-1
                iwad=iwalk_ad(jpat1,ipae,iwa,iwt)         !t1
             call prodel_1(1,wld_1,0,ipae,iwad,lrd,lrd)
             call prodel_2(1,wld_2,0,ipae,iwad,nxo)
c             call prodel(1,wld,0,ipae,iwad)
            enddo
          endif
            iwdownv=iw_downwei(jpas,ipae)
            do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpas,ipae,iwa,iws)
             call prodel_1(1,wld_1,0,ipae,iwad,lrd,lrd)
             call prodel_2(1,wld_2,0,ipae,iwad,nxo)
c           call prodel(1,wld,0,ipae,iwad)
          enddo
            if(jb_sys.gt.0) then
            iws1=just(lr,lr0)
              do iwa=0,iwdownv-1
                iwad=iwalk_ad(jpas,ipae,iwa,iws1)
             call prodel_1(1,wld_1,0,ipae,iwad,lrd,lrd)
             call prodel_2(1,wld_2,0,ipae,iwad,nxo)
c             call prodel(1,wld,0,ipae,iwad)
            enddo
          endif
        enddo
       endif
      enddo
200     continue
100   continue

!     wl8 = 1/2*hnil*(hnil-1)*vmd(lr,lr)+hnil*voint(lr,lr)   800
!     wl5=(vlop0-vlop1)*vo(lr0,lr)-2*vlop0*vmd(lr0,lr)
!     wl5= vlop0*vmd(lr,lr0)   -vlop1*vo(lr0,lr)   2000.11.26
c     wt0=dzero
      do lr0=2,norb_dz
        do lr=1,lr0-1
                        ! .........520...   vlop0=-2 vlop1=0
c         wt0=wt0-2.d0*vdint(lr,lr0)
           wt0_1=-2.0d0
           wt0_2=4.0d0
           call trans_ijkl_intpos(lr0,lr,lr0,lr,nxo)
           nxo1_0=nxo
           call trans_ijkl_intpos(lr0,lr0,lr,lr,nxo)
           nxo2_0=nxo
           jpad=1
           iwd=0
           do ipae_=1,25
              ipae=ipae_ ! ipae is in common block, is this necessary?
              if(nu_ae(ipae).eq.0) cycle
              iwdownv=iw_downwei(jpad,ipae)
              do iwa=0,iwdownv-1
                 iwad=iwalk_ad(jpad,ipae,iwa,iwd)
                 call prodel_2(1,wt0_1,0,ipae,iwad,nxo1_0)
                 call prodel_2(1,wt0_2,0,ipae,iwad,nxo2_0)
c                 call prodel(1,wt0,0,ipae,iwad)
              enddo
           enddo
        enddo
      enddo
      do 300 lrm=norb_frz+1,norb_dz
          mrm=mul_tab(lsm_inn(lrm),ns_sm)
          iws=just(lrm,lrm)
          iwd=jud(lrm)
          jpad=1+mrm
          jpad1=jpad+24
          jpas=17+ns_sm
          do lrd=2,norb_dz
             do lrds=1,lrd-1

             logic_lij=.false.
             call trans_ijkl_intpos(lrd,lrds,lrd,lrds,nxo)
             nxod_1=nxo
             call trans_ijkl_intpos(lrd,lrd,lrds,lrds,nxo)
             nxod_2=nxo

             do lr=1,lrm-1
               if(lr.eq.lrds.and.lrm.eq.lrd)   then
                   logic_lij=.true.
c              wld=wld+vdint(lr,lrm)
               wld_1=-1.0d0
               wld_2=2.0d0
               do ipae_=1,25
                  ! ipae is in common block, is this necessary?
                  ipae=ipae_
                  if(nu_ae(ipae).eq.0) cycle
                     iwdownv=iw_downwei(jpad,ipae)
                  do iwa=0,iwdownv-1
                     iwad=iwalk_ad(jpad,ipae,iwa,iwd)
                     call prodel_2(1,wld_1,0,ipae,iwad,nxod_1)
                     call prodel_2(1,wld_2,0,ipae,iwad,nxod_2)

c                     call prodel(1,wld,0,ipae,iwad)
                  enddo
                  if(jb_sys.gt.0) then
                     iwdownv=iw_downwei(jpad1,ipae)
                     do iwa=0,iwdownv-1
                        iwad=iwalk_ad(jpad1,ipae,iwa,iwd)
                     call prodel_2(1,wld_1,0,ipae,iwad,nxod_1)
                     call prodel_2(1,wld_2,0,ipae,iwad,nxod_2)

c                        call prodel(1,wld,0,ipae,iwad)
                     enddo
                  endif
               enddo
               exit
               endif
             enddo

             do lr0=lrm+1,norb_dz
               if(lrm.eq.lrds.and.lr0.eq.lrd)   then
                   logic_lij=.true.
c              wld=wld+vdint(lrm,lr0)
               wld_1=-1.0d0
               wld_2=2.0d0
               do ipae_=1,25
                  ! ipae is in common block, is this necessary?
                  ipae=ipae_
                  if(nu_ae(ipae).eq.0) cycle
                     iwdownv=iw_downwei(jpad,ipae)
                  do iwa=0,iwdownv-1
                     iwad=iwalk_ad(jpad,ipae,iwa,iwd)
                     call prodel_2(1,wld_1,0,ipae,iwad,nxod_1)
                     call prodel_2(1,wld_2,0,ipae,iwad,nxod_2)

c                     call prodel(1,wld,0,ipae,iwad)
                  enddo
                  if(jb_sys.gt.0) then
                     iwdownv=iw_downwei(jpad1,ipae)
                     do iwa=0,iwdownv-1
                        iwad=iwalk_ad(jpad1,ipae,iwa,iwd)
                     call prodel_2(1,wld_1,0,ipae,iwad,nxod_1)
                     call prodel_2(1,wld_2,0,ipae,iwad,nxod_2)

c                        call prodel(1,wld,0,ipae,iwad)
                     enddo
                  endif
               enddo
               exit
               endif
             enddo

c      wls=wld
c      do lr=lrm+1,norb_dz
c       wls=wls+vdint(lrm,lr)
c      enddo
c      do lr0=1,lrm-1
c       wls=wls+vdint(lr0,lrm)
c      enddo

      if(.not.logic_lij) then
           wt0_1=-2.0d0
           wt0_2=4.0d0
      do ipae_=1,25
        ipae=ipae_ ! ipae is in common block, is this necessary?
        if(nu_ae(ipae).eq.0) cycle
        iwdownv=iw_downwei(jpad,ipae)
        do iwa=0,iwdownv-1
          iwad=iwalk_ad(jpad,ipae,iwa,iwd)
          call prodel_2(1,wt0_1,0,ipae,iwad,nxod_1)
          call prodel_2(1,wt0_2,0,ipae,iwad,nxod_2)

c          call prodel(1,wld,0,ipae,iwad)
        enddo
        if(jb_sys.gt.0) then
        iwdownv=iw_downwei(jpad1,ipae)
        do iwa=0,iwdownv-1
          iwad=iwalk_ad(jpad1,ipae,iwa,iwd)
          call prodel_2(1,wt0_1,0,ipae,iwad,nxod_1)
          call prodel_2(1,wt0_2,0,ipae,iwad,nxod_2)

c            call prodel(1,wld,0,ipae,iwad)
        enddo
        endif
! ....s_520.
        iwdownv=iw_downwei(jpas,ipae)
        do iwa=0,iwdownv-1
          iwad=iwalk_ad(jpas,ipae,iwa,iws)
          call prodel_2(1,wt0_1,0,ipae,iwad,nxod_1)
          call prodel_2(1,wt0_2,0,ipae,iwad,nxod_2)

c          call prodel(1,wls,0,ipae,iwad)
        enddo
      enddo
      endif
      enddo
      enddo

300   continue

                                             ! ........520.
      do 400 lr0=norb_frz+1,norb_dz-1
        mr0=mul_tab(lsm_inn(lr0),ns_sm)
        do 401 lr=lr0+1,norb_dz
         mr=mul_tab(mr0,lsm_inn(lr))

         jpat=9+mr
         jpas=17+mr
          jpat1=jpat+24
         iws=just(lr0,lr)
         iwt=iws            !
            iws1=just(lr,lr0)

          do lrd=2,norb_dz
             do lrds=1,lrd-1

             logic_lij=.false.
             call trans_ijkl_intpos(lrd,lrds,lrd,lrds,nxo)
             nxos_1=nxo
             call trans_ijkl_intpos(lrd,lrd,lrds,lrds,nxo)
             nxos_2=nxo
             if(lr0.eq.lrds.and.lr.eq.lrd) then
                logic_lij=.true.
               if(jb_sys.eq.0) then
c           wls=wt0+3.d0*(voint(lr,lr0)-vdint(lr,lr0))
c          wlt=wt0+voint(lr,lr0)-3.d0*vdint(lr,lr0)
                   wls_1=1.0d0
                   wls_2=1.0d0
                   wlt_1=-1.0d0
                   wlt_2=1.0d0
                endif
           !  ..3.3  2*vo(lr0,lr)-4*vmd(lr0,lr)  w0=-2  w1=0
           !  ..2.1    vo(lr0,lr)+  vmd(lr0,lr)  w0=-1/2  w1=-3/2
           !  ..3.3  2*vo(lr0,lr)-4*vmd(lr0,lr)  w0=-2  w1=0
           !  ..2.2  - vo(lr0,lr)+  vmd(lr0,lr)  w0=-1/2  w1=1/2
               if(jb_sys.gt.0) then
                  db=jb_sys
                 w1=-(db+3)/(2.d0*db+2)
c           wls=wt0+(1.50d0-w1)*voint(lr,lr0)-3.d0*vdint(lr,lr0)
                  wls_1=-0.5d0-w1
                  wls_2=1.0d0
c          wlt=wt0+voint(lr,lr0)-3.d0*vdint(lr,lr0)
                  wlt_1=-1.0d0
                  wlt_2=1.0d0
                 w1=-(db-1)/(2.d0*db+2)
c          wls1=wt0+(1.50d0-w1)*voint(lr,lr0)-3.d0*vdint(lr,lr0)
                  wls1_1=-0.5d0-w1
                  wls1_2=1.0d0
            !    ..3.3  2*vo(lr0,lr)-4*vmd(lr0,lr)  w0=-2  w1=0
            !    ..1.2  (-1/2-w1)*vo(lr0,lr)+vmd(lr0,lr)
           !  w0=-1/2    w1=-(db-1)/(2*db+2)
               endif
          do ipae_=1,25
            ipae=ipae_ ! ipae is in common block, is this necessary?
            if(nu_ae(ipae).eq.0) cycle
            iwdownv=iw_downwei(jpat,ipae)
            do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpat,ipae,iwa,iwt)
              call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
              call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
c             call prodel(1,wlt,0,ipae,iwad)
            enddo
            iwdownv=iw_downwei(jpas,ipae)
            do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpas,ipae,iwa,iws)
              call prodel_2(1,wls_1,0,ipae,iwad,nxos_1)
              call prodel_2(1,wls_2,0,ipae,iwad,nxos_2)
c             call prodel(1,wls,0,ipae,iwad)
           enddo
            if(jb_sys.gt.0) then
              iwdownv=iw_downwei(jpas,ipae)
              do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpas,ipae,iwa,iws1)     ! s 1-2
                call prodel_2(1,wls1_1,0,ipae,iwad,nxos_1)
                call prodel_2(1,wls1_2,0,ipae,iwad,nxos_2)
c             call prodel(1,wls1,0,ipae,iwad)
            enddo
           endif
           if(jb_sys.gt.1) then
              iwdownv=iw_downwei(jpat1,ipae)
               do iwa=0,iwdownv-1
                 iwad=iwalk_ad(jpat1,ipae,iwa,iwt)     !t 1-1
                 call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
                 call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
c                 call prodel(1,wlt,0,ipae,iwad)
               enddo
           endif
          enddo
             cycle
            endif

        do lrg=1,lr0-1
           if(lrg.eq.lrds.and.lr0.eq.lrd) then
               logic_lij=.true.
               wls_1=-1.0d0
               wls_2=2.0d0
               wlt_1=-1.0d0
               wlt_2=2.0d0
               if(jb_sys.gt.0) then
                  wls1_1=-1.0d0
                  wls1_2=2.0d0
               endif
          do ipae_=1,25
            ipae=ipae_ ! ipae is in common block, is this necessary?
            if(nu_ae(ipae).eq.0) cycle
            iwdownv=iw_downwei(jpat,ipae)
            do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpat,ipae,iwa,iwt)
              call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
              call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
c             call prodel(1,wlt,0,ipae,iwad)
           enddo
            iwdownv=iw_downwei(jpas,ipae)
            do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpas,ipae,iwa,iws)
              call prodel_2(1,wls_1,0,ipae,iwad,nxos_1)
              call prodel_2(1,wls_2,0,ipae,iwad,nxos_2)
c             call prodel(1,wls,0,ipae,iwad)
           enddo
            if(jb_sys.gt.0) then
              iwdownv=iw_downwei(jpas,ipae)
              do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpas,ipae,iwa,iws1)     ! s 1-2
                call prodel_2(1,wls1_1,0,ipae,iwad,nxos_1)
                call prodel_2(1,wls1_2,0,ipae,iwad,nxos_2)
c             call prodel(1,wls1,0,ipae,iwad)
            enddo
           endif
           if(jb_sys.gt.1) then
              iwdownv=iw_downwei(jpat1,ipae)
               do iwa=0,iwdownv-1
                 iwad=iwalk_ad(jpat1,ipae,iwa,iwt)     !t 1-1
                 call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
                 call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
c                 call prodel(1,wlt,0,ipae,iwad)
               enddo
           endif
          enddo
          exit
         endif

         if(lrg.eq.lrds.and.lr.eq.lrd) then
               logic_lij=.true.
               wls_1=-1.0d0
               wls_2=2.0d0
               wlt_1=-1.0d0
               wlt_2=2.0d0
               if(jb_sys.gt.0) then
                  wls1_1=-1.0d0
                  wls1_2=2.0d0
               endif
          do ipae_=1,25
            ipae=ipae_ ! ipae is in common block, is this necessary?
            if(nu_ae(ipae).eq.0) cycle
            iwdownv=iw_downwei(jpat,ipae)
            do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpat,ipae,iwa,iwt)
              call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
              call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
c             call prodel(1,wlt,0,ipae,iwad)
           enddo
            iwdownv=iw_downwei(jpas,ipae)
            do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpas,ipae,iwa,iws)
              call prodel_2(1,wls_1,0,ipae,iwad,nxos_1)
              call prodel_2(1,wls_2,0,ipae,iwad,nxos_2)
c             call prodel(1,wls,0,ipae,iwad)
           enddo
            if(jb_sys.gt.0) then
              iwdownv=iw_downwei(jpas,ipae)
              do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpas,ipae,iwa,iws1)     ! s 1-2
                call prodel_2(1,wls1_1,0,ipae,iwad,nxos_1)
                call prodel_2(1,wls1_2,0,ipae,iwad,nxos_2)
c             call prodel(1,wls1,0,ipae,iwad)
            enddo
           endif
           if(jb_sys.gt.1) then
              iwdownv=iw_downwei(jpat1,ipae)
               do iwa=0,iwdownv-1
                 iwad=iwalk_ad(jpat1,ipae,iwa,iwt)     !t 1-1
                 call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
                 call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
c                 call prodel(1,wlt,0,ipae,iwad)
               enddo
           endif
          enddo
          exit
         endif
        enddo

        do lrg=lr0+1,lr-1
           if(lr0.eq.lrds.and.lrg.eq.lrd) then
               logic_lij=.true.
               wls_1=-1.0d0
               wls_2=2.0d0
               wlt_1=-1.0d0
               wlt_2=2.0d0
               if(jb_sys.gt.0) then
                  wls1_1=-1.0d0
                  wls1_2=2.0d0
               endif
          do ipae_=1,25
            ipae=ipae_ ! ipae is in common block, is this necessary?
            if(nu_ae(ipae).eq.0) cycle
            iwdownv=iw_downwei(jpat,ipae)
            do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpat,ipae,iwa,iwt)
              call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
              call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
c             call prodel(1,wlt,0,ipae,iwad)
           enddo
            iwdownv=iw_downwei(jpas,ipae)
            do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpas,ipae,iwa,iws)
              call prodel_2(1,wls_1,0,ipae,iwad,nxos_1)
              call prodel_2(1,wls_2,0,ipae,iwad,nxos_2)
c             call prodel(1,wls,0,ipae,iwad)
           enddo
            if(jb_sys.gt.0) then
              iwdownv=iw_downwei(jpas,ipae)
              do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpas,ipae,iwa,iws1)     ! s 1-2
                call prodel_2(1,wls1_1,0,ipae,iwad,nxos_1)
                call prodel_2(1,wls1_2,0,ipae,iwad,nxos_2)
c             call prodel(1,wls1,0,ipae,iwad)
            enddo
           endif
           if(jb_sys.gt.1) then
              iwdownv=iw_downwei(jpat1,ipae)
               do iwa=0,iwdownv-1
                 iwad=iwalk_ad(jpat1,ipae,iwa,iwt)     !t 1-1
                 call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
                 call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
c                 call prodel(1,wlt,0,ipae,iwad)
               enddo
           endif
          enddo
          exit
         endif

         if(lrg.eq.lrds.and.lr.eq.lrd) then
               logic_lij=.true.
               wls_1=-1.0d0
               wls_2=2.0d0
               wlt_1=-1.0d0
               wlt_2=2.0d0
               if(jb_sys.gt.0) then
                  wls1_1=-1.0d0
                  wls1_2=2.0d0
               endif
          do ipae_=1,25
            ipae=ipae_ ! ipae is in common block, is this necessary?
            if(nu_ae(ipae).eq.0) cycle
            iwdownv=iw_downwei(jpat,ipae)
            do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpat,ipae,iwa,iwt)
              call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
              call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
c             call prodel(1,wlt,0,ipae,iwad)
           enddo
            iwdownv=iw_downwei(jpas,ipae)
            do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpas,ipae,iwa,iws)
              call prodel_2(1,wls_1,0,ipae,iwad,nxos_1)
              call prodel_2(1,wls_2,0,ipae,iwad,nxos_2)
c             call prodel(1,wls,0,ipae,iwad)
           enddo
            if(jb_sys.gt.0) then
              iwdownv=iw_downwei(jpas,ipae)
              do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpas,ipae,iwa,iws1)     ! s 1-2
                call prodel_2(1,wls1_1,0,ipae,iwad,nxos_1)
                call prodel_2(1,wls1_2,0,ipae,iwad,nxos_2)
c             call prodel(1,wls1,0,ipae,iwad)
            enddo
           endif
           if(jb_sys.gt.1) then
              iwdownv=iw_downwei(jpat1,ipae)
               do iwa=0,iwdownv-1
                 iwad=iwalk_ad(jpat1,ipae,iwa,iwt)     !t 1-1
                 call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
                 call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
c                 call prodel(1,wlt,0,ipae,iwad)
               enddo
           endif
          enddo
          exit
         endif
        enddo


        do lrg=lr+1,norb_dz
           if(lr.eq.lrds.and.lrg.eq.lrd) then
               logic_lij=.true.
               wls_1=-1.0d0
               wls_2=2.0d0
               wlt_1=-1.0d0
               wlt_2=2.0d0
               if(jb_sys.gt.0) then
                  wls1_1=-1.0d0
                  wls1_2=2.0d0
               endif
          do ipae_=1,25
            ipae=ipae_ ! ipae is in common block, is this necessary?
            if(nu_ae(ipae).eq.0) cycle
            iwdownv=iw_downwei(jpat,ipae)
            do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpat,ipae,iwa,iwt)
              call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
              call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
c             call prodel(1,wlt,0,ipae,iwad)
           enddo
            iwdownv=iw_downwei(jpas,ipae)
            do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpas,ipae,iwa,iws)
              call prodel_2(1,wls_1,0,ipae,iwad,nxos_1)
              call prodel_2(1,wls_2,0,ipae,iwad,nxos_2)
c             call prodel(1,wls,0,ipae,iwad)
           enddo
            if(jb_sys.gt.0) then
              iwdownv=iw_downwei(jpas,ipae)
              do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpas,ipae,iwa,iws1)     ! s 1-2
                call prodel_2(1,wls1_1,0,ipae,iwad,nxos_1)
                call prodel_2(1,wls1_2,0,ipae,iwad,nxos_2)
c             call prodel(1,wls1,0,ipae,iwad)
            enddo
           endif
           if(jb_sys.gt.1) then
              iwdownv=iw_downwei(jpat1,ipae)
               do iwa=0,iwdownv-1
                 iwad=iwalk_ad(jpat1,ipae,iwa,iwt)     !t 1-1
                 call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
                 call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
c                 call prodel(1,wlt,0,ipae,iwad)
               enddo
           endif
          enddo
          exit
         endif

         if(lr0.eq.lrds.and.lrg.eq.lrd) then
               logic_lij=.true.
               wls_1=-1.0d0
               wls_2=2.0d0
               wlt_1=-1.0d0
               wlt_2=2.0d0
               if(jb_sys.gt.0) then
                  wls1_1=-1.0d0
                  wls1_2=2.0d0
               endif
          do ipae_=1,25
            ipae=ipae_ ! ipae is in common block, is this necessary?
            if(nu_ae(ipae).eq.0) cycle
            iwdownv=iw_downwei(jpat,ipae)
            do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpat,ipae,iwa,iwt)
              call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
              call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
c             call prodel(1,wlt,0,ipae,iwad)
           enddo
            iwdownv=iw_downwei(jpas,ipae)
            do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpas,ipae,iwa,iws)
              call prodel_2(1,wls_1,0,ipae,iwad,nxos_1)
              call prodel_2(1,wls_2,0,ipae,iwad,nxos_2)
c             call prodel(1,wls,0,ipae,iwad)
           enddo
            if(jb_sys.gt.0) then
              iwdownv=iw_downwei(jpas,ipae)
              do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpas,ipae,iwa,iws1)     ! s 1-2
                call prodel_2(1,wls1_1,0,ipae,iwad,nxos_1)
                call prodel_2(1,wls1_2,0,ipae,iwad,nxos_2)
c             call prodel(1,wls1,0,ipae,iwad)
            enddo
           endif
           if(jb_sys.gt.1) then
              iwdownv=iw_downwei(jpat1,ipae)
               do iwa=0,iwdownv-1
                 iwad=iwalk_ad(jpat1,ipae,iwa,iwt)     !t 1-1
                 call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
                 call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
c                 call prodel(1,wlt,0,ipae,iwad)
               enddo
           endif
          enddo
          exit
         endif
        enddo

        if(.not.logic_lij) then
           if(jb_sys.eq.0) then
              wls_1=-2.0d0
              wls_2=4.0d0
              wlt_1=-2.0d0
              wlt_2=4.0d0
           endif

           if(jb_sys.gt.0) then
              wls_1=-2.0d0
              wls_2=4.0d0
              wlt_1=-2.0d0
              wlt_2=4.0d0
              wls1_1=-2.0d0
              wls1_2=4.0d0
           endif

          do ipae_=1,25
            ipae=ipae_ ! ipae is in common block, is this necessary?
            if(nu_ae(ipae).eq.0) cycle
            iwdownv=iw_downwei(jpat,ipae)
            do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpat,ipae,iwa,iwt)
              call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
              call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
c             call prodel(1,wlt,0,ipae,iwad)
           enddo
            iwdownv=iw_downwei(jpas,ipae)
            do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpas,ipae,iwa,iws)
              call prodel_2(1,wls_1,0,ipae,iwad,nxos_1)
              call prodel_2(1,wls_2,0,ipae,iwad,nxos_2)
c             call prodel(1,wls,0,ipae,iwad)
           enddo
            if(jb_sys.gt.0) then
              iwdownv=iw_downwei(jpas,ipae)
              do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpas,ipae,iwa,iws1)     ! s 1-2
                call prodel_2(1,wls1_1,0,ipae,iwad,nxos_1)
                call prodel_2(1,wls1_2,0,ipae,iwad,nxos_2)
c             call prodel(1,wls1,0,ipae,iwad)
            enddo
           endif
           if(jb_sys.gt.1) then
              iwdownv=iw_downwei(jpat1,ipae)
               do iwa=0,iwdownv-1
                 iwad=iwalk_ad(jpat1,ipae,iwa,iwt)     !t 1-1
                 call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
                 call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
c                 call prodel(1,wlt,0,ipae,iwad)
               enddo
           endif
          enddo

        endif
        enddo
      enddo
401     continue
400   continue
c ------------- end of .......h_delm --------------
      return
      end



      subroutine diagonal_ext_g()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
!     common/casrst/ja(max_node),jb(max_node),jm(0:max_node)
!     :    ,jj(4,0:max_node),kk(0:max_node),no(0:max_innorb)
!     :    ,jv,jd(8),jt(8),js(8)
!      common/sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),
!     :     jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)
        jws0=0
c      do mra=1,8
       do mra=1,ng_sm
        ipae=1+mra
        lasta=ibsm_ext(mra)
        laend=iesm_ext(mra)
        lrzz=laend-lasta+1
        lrzz=lrzz*(lrzz-1)/2
        jws0=jws0+lrzz
        jw=0
        do la=lasta,laend
c        jpd=jd(mra)
          jw=jw+1
          lra=norb_all-la+1
c        wld=voint(lra,lra)
          wld=1.0d0
          call prodel_1(2,wld,0,ipae,jw,lra,lra)
c         call prodel(2,wld,0,ipae,jw)
        enddo
      enddo

!      zz=' out_800_s'
c       jps=js(1)

      jweis=jws0
      do la=1,norb_ext
c        jpd=jd(mra)
        lra=norb_all-la+1
        jweis=jweis+1
c        wls=2.d0*voint(lra,lra)+vdint(lra,lra)
c       call prodel(2,wls,0,18,jweis)
        wls=2.0d0
        call prodel_1(2,wls,0,18,jweis,lra,lra)
        wls=1.0d0
        call trans_ijkl_intpos(lra,lra,lra,lra,nxo)
        call prodel_2(2,wls,0,18,jweis,nxo)
      enddo
c      do im=1,8
       do im=1,ng_sm
        jws=0
        jwt=0
        ipat=9+im
        ipas=17+im
        do la=2,norb_ext
         lra=norb_all-la+1
         ima=lsm(la)
         do 600 lb=1,la-1
           lrb=norb_all-lb+1
           imb=lsm(lb)
           mr=mul_tab(ima,imb)
           if(mr.ne.im) goto 600
c          jps=js(mr)
c          jpt=jt(mr)
           jws=jws+1
           jwt=jwt+1
c            wls=voint(lra,lra)+voint(lrb,lrb)
            wls=1.0d0
            wlt=wls
           call prodel_1(2,wls,0,ipas,jws,lra,lra)
            call prodel_1(2,wlt,0,ipat,jwt,lra,lra)
           call prodel_1(2,wls,0,ipas,jws,lrb,lrb)
            call prodel_1(2,wlt,0,ipat,jwt,lrb,lrb)

c           wls=wls+voint(lrb,lra)+vdint(lrb,lra)   ! w0=-1/2  w1=-3/2
           wls=1.0d0
           call trans_ijkl_intpos(lrb,lra,lrb,lra,nxo)
           call prodel_2(2,wls,0,ipas,jws,nxo)
           call trans_ijkl_intpos(lrb,lrb,lra,lra,nxo)
           call prodel_2(2,wls,0,ipas,jws,nxo)
c           wlt=wlt-voint(lrb,lra)+vdint(lrb,lra)   ! w0=-1/2  w1=1/2
           wlt=-1.0d0
           call trans_ijkl_intpos(lrb,lra,lrb,lra,nxo)
           call prodel_2(2,wlt,0,ipat,jwt,nxo)
           wlt=1.0d0
           call trans_ijkl_intpos(lrb,lrb,lra,lra,nxo)
           call prodel_2(2,wlt,0,ipat,jwt,nxo)

c          call prodel(2,wls,0,ipas,jws)
c           call prodel(2,wlt,0,ipat,jwt)
600       continue
          enddo
      enddo
      return
      end

!  idb=1  in dbl_space      ity_up=0-5             jd_type,jd_im,iwd
!  idb=2  in ext_space      ity_down=0-3          je_type,je_im,iwe
!  idb=3  in act_space      ity_up=0-5,itdown=0,3      jp ,mpe,iwa
!  idb=4  betwin dbl and act   ity_up=0-5,itdown=0,3      mpe,iwd,iwa
!  idb=5  betwin act and ext  ity_down=0-3            jp ,iwa,iwe
!  idb=6  betwin dbl and ext   ity_down=0-3            iwd,iwa,iwe

!  this subroutine prodel_1 does the dm1 part, which corresponds to voin

      subroutine prodel_1(idb,wl,mg1,mg2,mg3,mg6,mg7)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "grad_h.fh"
!      common/sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),
!     :     jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)
      common /iaib/ ican_a(max_orb),ican_b(mtmp+max_orb)

!      ndr=88
      goto(100,200,300,400,500,600),idb
! in dbl_space
100   ipae=mg2
      iwad=mg3
        isegdownwei=iseg_downwei(ipae)
        do mm=iwad+1,iwad+isegdownwei
c        vector1(mm)=vector1(mm)+wl
            mg67=ican_a(mg7)+mg6
            dm1tmp(mg67)=dm1tmp(mg67)+vector1(mm)*wl*vector1(mm)
c         if(mg7.eq.4)
c     :      write(nf2,'(3i4,f18.10)')mg7,mg6,mm,wl

c         write(nf2,'(3i4,3f18.10)')mg7,mg6,mm,wl,vector1(mm),vector1(mm
!     if(mm.eq.ndr) then
!       write(nf2,'(a8,3i6,2f20.14)')
!     :     ' in dbl _',mg1,mg2,mg3,wl,vector1(mm)
!       write(*,'(a8,3i6,2f20.14)')
!     :     ' in dbl _',mg1,mg2,mg3,wl,vector1(mm)
!      endif
        enddo
      goto 1000
! in ext_space
 200  ipae=mg2
      iwe=mg3
      do jdbl=1,mxnode
        if(nu_ad(jdbl).eq.0) cycle
        iw=iw_downwei(jdbl,ipae)
        iwupwei=jpad_upwei(jdbl)
        do iwa=0,iw-1
          do iwd=0,iwupwei-1
            mm=iwalk_ad(jdbl,ipae,iwa,iwd)+iwe
c          vector1(mm)=vector1(mm)+wl
            mg67=ican_a(mg7)+mg6
            dm1tmp(mg67)=dm1tmp(mg67)+vector1(mm)*wl*vector1(mm)

!     if(mm.eq.ndr) then
!       write(nf2,'(a8,3i6,2f20.14)')
!     :     ' in ext _',mg1,mg2,mg3,wl,vector1(mm)
!       write(*,'(a8,3i6,2f20.14)')
!     :     ' in ext _',mg1,mg2,mg3,wl,vector1(mm)
!      endif
          enddo
        enddo
      enddo
      goto 1000
! in act_space
300   jp=mg1
      mpe=mg2
      jw=mg3
      iwupwei=jpad_upwei(jpad)
      isegdownwei=iseg_downwei(ipae)
      jph=jphy(jp)
      in=ihy(jph)
      lwnu=iy(1,mpe)
      do jwu=jph+1,jph+in
        iwa=jw+ihy(jwu)-1
        do jwd=1,lwnu
          iwa=iwa+1
          do iwd=0,iwupwei-1
            iwad=iwalk_ad(jpad,ipae,iwa,iwd)
            do iwe=1,isegdownwei
              mm=iwe+iwad
c            vector1(mm)=vector1(mm)+wl
            mg67=ican_a(mg7)+mg6
            dm1tmp(mg67)=dm1tmp(mg67)+vector1(mm)*wl*vector1(mm)

!     if(mm.eq.ndr) then
!       write(nf2,'(a8,3i6,2f20.14)')
!     :     ' in act _',mg1,mg2,mg3,wl,vector1(mm)
!       write(*,'(a8,3i6,2f20.14)')
!     :     ' in act _',mg1,mg2,mg3,wl,vector1(mm)
!      endif
            enddo
          enddo
        enddo
      enddo
      goto 1000
! betwin dbl and act
400   mpe=mg1
      iwd=mg2
      iwa=mg3-1
      isegdownwei=iseg_downwei(ipae)   ! betwin dbl and act
      jwnu=iy(1,mpe)
      do ii=1,jwnu
        iwa=iwa+1
        iwad=iwalk_ad(jpad,ipae,iwa,iwd)
        do iwe=1,isegdownwei
          mm=iwe+iwad                  ! iwl=iwalk_ad
c             vector1(mm)=vector1(mm)+wl
            mg67=ican_a(mg7)+mg6
            dm1tmp(mg67)=dm1tmp(mg67)+vector1(mm)*wl*vector1(mm)

!     if(mm.eq.ndr) then
!       write(nf2,'(a8,3i6,2f20.14)')
!     :     ' dbl_act ',mg1,mg2,mg3,wl,vector1(mm)
!       write(*,'(a8,3i6,2f20.14)')
!     :     ' dbl_act ',mg1,mg2,mg3,wl,vector1(mm)
!      endif
        enddo
      enddo
      goto 1000
! betwin act and ext
500   jp=mg1
      iwa0=mg2
      iwe=mg3
      iwupwei=jpad_upwei(jpad)
      jph=jphy(jp)
      in=ihy(jph)
      do jwu=jph+1,jph+in
        iwa=iwa0+ihy(jwu)
        do iwd=0,iwupwei-1
          iwad=iwalk_ad(jpad,ipae,iwa,iwd)
          mm=iwe+iwad
c          vector1(mm)=vector1(mm)+wl
            mg67=ican_a(mg7)+mg6
            dm1tmp(mg67)=dm1tmp(mg67)+vector1(mm)*wl*vector1(mm)

!     if(mm.eq.ndr) then
!       write(nf2,'(a8,3i6,2f20.14)')
!     :     ' act_ext ',mg1,mg2,mg3,wl,vector1(mm)
!       write(*,'(a8,3i6,2f20.14)')
!     :     ' act_ext ',mg1,mg2,mg3,wl,vector1(mm)
!      endif
        enddo
      enddo
      goto 1000
! betwin dbl and ext
600   iwd=mg1
      iwa=mg2
      iwe=mg3
      iwad=iwalk_ad(jpad,ipae,iwa,iwd)   ! betwin dbl,act and ext
      mm=iwe+iwad
c          vector1(mm)=vector1(mm)+wl
            mg67=ican_a(mg7)+mg6
            dm1tmp(mg67)=dm1tmp(mg67)+vector1(mm)*wl*vector1(mm)

!     if(mm.eq.ndr) then
!       write(nf2,'(a8,3i6,2f20.14)')
!     :     ' dbl_ext ',mg1,mg2,mg3,wl,vector1(mm)
!       write(*,'(a8,3i6,2f20.14)')
!     :     ' dbl_ext ',mg1,mg2,mg3,wl,vector1(mm)
!      endif
1000  return
      end

!this subroutine prodel_2 does the dm2 part, which corresponds to vint_c

      subroutine prodel_2(idb,wl,mg1,mg2,mg3,mg6)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "grad_h.fh"
!      common/sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),
!     :     jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)
c      ndr=88
c      ndr=6
      goto(100,200,300,400,500,600),idb
! in dbl_space
100   ipae=mg2
      iwad=mg3
        isegdownwei=iseg_downwei(ipae)
        do mm=iwad+1,iwad+isegdownwei
c        vector1(mm)=vector1(mm)+wl
           vector2(mg6)=vector2(mg6)+vector1(mm)*wl*vector1(mm)
c      if(mg6.eq.45)
c     :     write(nf2,'(i8,2i4,4f18.10)')mg6,mm,mm,vector2(mg6),
c     :                          vector1(mm),wl,vector1(mm)

c      if(mg6.eq.ndr) then
c       write(nf2,'(a9,3i6,3f20.14)')
c     :     ' in dbl_ ',mg6,mm,mm,wl,vector1(mm),vector2(mg6)
c        write(*,'(a8,3i6,2f20.14)')
c     :     ' in dbl _',mg1,mg2,mg3,wl,vector1(mm)
c      endif
        enddo
      goto 1000
! in ext_space
 200  ipae=mg2
      iwe=mg3
      do jdbl=1,mxnode
        if(nu_ad(jdbl).eq.0) cycle
        iw=iw_downwei(jdbl,ipae)
        iwupwei=jpad_upwei(jdbl)
        do iwa=0,iw-1
          do iwd=0,iwupwei-1
            mm=iwalk_ad(jdbl,ipae,iwa,iwd)+iwe
c          vector1(mm)=vector1(mm)+wl
           vector2(mg6)=vector2(mg6)+vector1(mm)*wl*vector1(mm)
c          write(nf2,'(i8,2i4,4f18.10)')mg6,mm,mm,vector2(mg6),
c     :                          vector1(mm),wl,vector1(mm)


c      if(mg6.eq.2926) then
c       write(nf2,'(a9,3i6,3f20.14)')
c     :     ' in ext_ ',mg6,mm,mm,wl,vector1(mm),vector2(mg6)
!       write(*,'(a8,3i6,2f20.14)')
!     :     ' in ext _',mg1,mg2,mg3,wl,vector1(mm)

c      endif
          enddo
        enddo
      enddo
      goto 1000
! in act_space
300   jp=mg1
      mpe=mg2
      jw=mg3
      iwupwei=jpad_upwei(jpad)
      isegdownwei=iseg_downwei(ipae)
      jph=jphy(jp)
      in=ihy(jph)
      lwnu=iy(1,mpe)
      do jwu=jph+1,jph+in
        iwa=jw+ihy(jwu)-1
        do jwd=1,lwnu
          iwa=iwa+1
          do iwd=0,iwupwei-1
            iwad=iwalk_ad(jpad,ipae,iwa,iwd)
            do iwe=1,isegdownwei
              mm=iwe+iwad
c            vector1(mm)=vector1(mm)+wl
           vector2(mg6)=vector2(mg6)+vector1(mm)*wl*vector1(mm)
c      if(mg6.eq.231)  write(nf2,'(i4,3f18.10)') mm,wl,vector1(mm),
c     :     vector2(mg6)
c           write(nf2,'(a3,i4)') 'act',mm

c          write(nf2,'(i8,2i4,4f18.10)')mg6,mm,mm,vector2(mg6),
c     :                          vector1(mm),wl,vector1(mm)

c      if(mg6.eq.ndr) then
c       write(nf2,'(a9,3i6,3f20.14)')
c     :     ' in act_ ',mg6,mm,mm,wl,vector1(mm),vector2(mg6)
!       write(*,'(a8,3i6,2f20.14)')
!     :     ' in act _',mg1,mg2,mg3,wl,vector1(mm)
c      endif
            enddo
          enddo
        enddo
      enddo
      goto 1000
! betwin dbl and act
400   mpe=mg1
      iwd=mg2
      iwa=mg3-1
      isegdownwei=iseg_downwei(ipae)   ! betwin dbl and act
      jwnu=iy(1,mpe)
      do ii=1,jwnu
        iwa=iwa+1
        iwad=iwalk_ad(jpad,ipae,iwa,iwd)
        do iwe=1,isegdownwei
          mm=iwe+iwad                  ! iwl=iwalk_ad
c             vector1(mm)=vector1(mm)+wl
           vector2(mg6)=vector2(mg6)+vector1(mm)*wl*vector1(mm)
c           if(mg6.eq.91)
c     :      write(nf2,'(i8,i4,3f18.10)')mg6,mm,vector2(mg6),
c     :                          wl,vector1(mm)
c        write(nf2,'(a7,i4)') 'dbl_act',mm

c      if(mg6.eq.3)  write(nf2,'(i4,2f18.10)') mm,wl,vector1(mm)
c      if(mg6.eq.ndr) then
c       write(nf2,'(a9,3i6,3f20.14)')
c     :     ' dbl_act ',mg6,mm,mm,wl,vector1(mm),vector2(mg6)
c!       write(*,'(a8,3i6,2f20.14)')
!     :     ' dbl_act ',mg1,mg2,mg3,wl,vector1(mm)
c      endif
        enddo
      enddo
      goto 1000
! betwin act and ext
500   jp=mg1
      iwa0=mg2
      iwe=mg3
      iwupwei=jpad_upwei(jpad)
      jph=jphy(jp)
      in=ihy(jph)
      do jwu=jph+1,jph+in
        iwa=iwa0+ihy(jwu)
        do iwd=0,iwupwei-1
          iwad=iwalk_ad(jpad,ipae,iwa,iwd)
          mm=iwe+iwad
c          vector1(mm)=vector1(mm)+wl
           vector2(mg6)=vector2(mg6)+vector1(mm)*wl*vector1(mm)
c           write(nf2,*) 'ae',mm,wl,vector1(mm)
c          write(nf2,'(i8,2i4,4f18.10)')mg6,mm,mm,vector2(mg6),
c     :                          vector1(mm),wl,vector1(mm)


c      if(mg6.eq.ndr) then
c       write(nf2,'(a9,3i6,3f20.14)')
c     :     ' act_ext ',mg6,mm,mm,wl,vector1(mm),vector2(mg6)
!       write(*,'(a8,3i6,2f20.14)')
!     :     ' act_ext ',mg1,mg2,mg3,wl,vector1(mm)
c      endif
        enddo
      enddo
      goto 1000
! betwin dbl and ext
600   iwd=mg1
      iwa=mg2
      iwe=mg3
      iwad=iwalk_ad(jpad,ipae,iwa,iwd)   ! betwin dbl,act and ext
      mm=iwe+iwad
c          vector1(mm)=vector1(mm)+wl
           vector2(mg6)=vector2(mg6)+vector1(mm)*wl*vector1(mm)
c          write(nf2,'(i8,2i4,4f18.10)')mg6,mm,mm,vector2(mg6),
c     :                          vector1(mm),wl,vector1(mm)


c      if(mg6.eq.ndr) then
c       write(nf2,'(a9,3i6,3f20.14)')
c     :     ' dbl_ext ',mg6,mm,mm,wl,vector1(mm),vector2(mg6)
!       write(*,'(a8,3i6,2f20.14)')
!     :     ' dbl_ext ',mg1,mg2,mg3,wl,vector1(mm)
c      endif
1000  return
      end

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
c generate and print csfs
      subroutine found_a_config(ndl,de,npr)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      common /mcorb/ lsmorb(max_orb),noidx(8)
      common/config/ndr,nwalk(0:max_orb)
!      common/casrst/ja(max_node),jb(max_node),jm(0:max_node)
!     :    ,jj(4,0:max_node),kk(0:max_node),no(0:max_innorb)
!     :    ,jv,jd(8),jt(8),js(8)
!      common/sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),
!     :     jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)
      dimension  nwalk_gamess(max_orb),norbindex(max_orb),
     :     norbsymm(max_orb),norbwalk(max_orb)
      character*18 form1

      ndr=ndl
      nst=norb_all
      ndy=norb_ext
      do l=1,norb_ext+norb_act
        nwalk(l)=0
      enddo
      no_dz=nst-norb_dz+1
      do l=no_dz,nst
        nwalk(l)=3
      enddo
      ndimsum=0
      if(npr.eq.1.or.npr.eq.0) then
        jpae=jv
        ipae=1
        jaedownwei=iseg_downwei(ipae)
        jpad=1
        iw_sta(jpad,ipae)=ndimsum
        call seg_drt()                     !   jpad_upwei(*)        = jp
        iwupwei=jpad_upwei(jpad)               !   iw_sta(jpad,jpae)
        iw_downwei(jpad,ipae)=ndim             !   iw_downwei(jpad,jpae)
        ndimsum=ndimsum+ndim*jaedownwei*iwupwei  !   iseg_dim(jpae)    =
        call config_act()                       !   jpae_upwei(jpae)  =
        goto 110
      endif

      jpae=jv
      ipae=1
      jaedownwei=iseg_downwei(ipae)
      do jpad=1,mxnode
        iw_sta(jpad,ipae)=ndimsum
        if(nu_ad(jpad).eq.0) cycle
        call seg_drt()                     !   jpad_upwei(*)        = jp
        iwupwei=jpad_upwei(jpad)               !   iw_sta(jpad,jpae)
        iw_downwei(jpad,ipae)=ndim             !   iw_downwei(jpad,jpae)
        ndimsum=ndimsum+ndim*jaedownwei*iwupwei  !   iseg_dim(jpae)    =
        if(ndim .eq. 0) cycle                !   jpae_sta(jpae)    = jpa
        call config_act()                       !   jpae_upwei(jpae)  =
      enddo
      do im=1,8
        jpae=jd(im)
        ipae=1+im
        iw_sta(1:mxnode,ipae)=ndimsum
        if(nu_ae(ipae).eq.0) cycle
        jaedownwei=iseg_downwei(ipae)
        do jpad=1,mxnode
          iw_sta(jpad,ipae)=ndimsum
          if(nu_ad(jpad).eq.0) cycle
          call seg_drt()
          iwupwei=jpad_upwei(jpad)
          iw_downwei(jpad,ipae)=ndim
          ndimsum=ndimsum+ndim*jaedownwei*iwupwei
          if(ndim .eq. 0) cycle
          call config_act()
          ista=iw_sta(jpad,ipae)+1
        enddo
      enddo
      do im=1,8
        jpae=jt(im)
        ipae=9+im
        iw_sta(1:mxnode,ipae)=ndimsum
        if(nu_ae(ipae).eq.0) cycle
        jaedownwei=iseg_downwei(ipae)
        do jpad=1,mxnode
          iw_sta(jpad,ipae)=ndimsum
          if(nu_ad(jpad).eq.0) cycle
          call seg_drt()
          iwupwei=jpad_upwei(jpad)
          iw_downwei(jpad,ipae)=ndim
          ndimsum=ndimsum+ndim*jaedownwei*iwupwei
          if(ndim .eq. 0) cycle
          call config_act()
        enddo
      enddo
      do im=1,8
        jpae=js(im)
        ipae=17+im
        if(nu_ae(ipae).eq.0) cycle
        do jpad=1,mxnode
          if(nu_ad(jpad).eq.0) cycle
          call seg_drt()
          if(ndim .eq. 0) cycle
          call config_act()
        enddo
      enddo
      call config_dbl()
      call config_ext()

110   continue
      if(npr.eq.0) return

      nwalk_gamess(1:norb_all)=0
      do i=1,nst
        nwalk_gamess(i)=nwalk(nst-i+1)
      enddo
      if(npr.eq.1) then
        write(6,1000) ndr,de
      else
        sqde=de*de
        write(6,1001) ndr,de,sqde
      endif

      if(intgen.eq.1) then
        ! gamess integral
        i=0
        j=0
        do i=norb_frz+1,nst
          if(nwalk_gamess(map_orb_order(i)).gt.0) then
            j=j+1
            norbindex(j)=i
            norbsymm(j)=lsmorb(i)
            norbwalk(j)=nwalk_gamess(map_orb_order(i))
          endif
        enddo

        !write(form1,2010) "(4x,i4,1x,",ng_sm,"(1x","i"
        write(form1,1010) "(4x,a4,",j,"(1x,i3))"
        write(6,form1) "norb",(norbindex(i),i=1,j)
        write(6,form1) "sym ",(norbsymm(i),i=1,j)
        write(6,form1) "walk",(norbwalk(i),i=1,j)
      else
        i=0
        j=0
        do im=1,ng_sm
          ns=noidx(im)+nlsm_frz(im)+1
          if(im.eq.ng_sm) then
            ne=norb_all
          else
            ne=noidx(im+1)
          endif
          do i=ns,ne
            if(nwalk_gamess(map_orb_order(i)).gt.0) then
              j=j+1
              noi=i-noidx(im)
              norbindex(j)=noi
              norbsymm(j)=im
              norbwalk(j)=nwalk_gamess(map_orb_order(i))
            endif
          enddo
        enddo

        write(form1,1010) "(4x,a4,",j,"(1x,i3))"
        write(6,form1) "norb",(norbindex(i),i=1,j)
        write(6,form1) "sym ",(norbsymm(i),i=1,j)
        write(6,form1) "walk",(norbwalk(i),i=1,j)
      endif

      return
1000  format(/4x,"csf",i8,6x,"scf energy",f18.8,/)
1001  format(/4x,"csf",i8,2x,"coef",f10.6,2x,"weight",1x,f8.6/)
1010  format(a7,i3,a8)
      end

      subroutine config_act()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      common/config/ndr,nwalk(0:max_orb)
!      common/casrst/ja(max_node),jb(max_node),jm(0:max_node)
!     :    ,jj(4,0:max_node),kk(0:max_node),no(0:max_innorb)
!     :    ,jv,jd(8),jt(8),js(8)
!      common/sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),
!     :     jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)
!      dimension ndr(max_innorb)
      REAL*8, pointer :: jph(:),jeh(:),jwh(:)
c      common/ptlph/jph,jeh,jwh,th,thh
      common/ptlph/jph,jeh,jwh
c     write(6,*)'               ***** start h-diaelm *****'
c      write(6,*)   '   diagonal_act_d:',jpad,ipae
      nsum=0
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
            jw=iy(idl,jp)
c           call prodel(3,wt,jp,mpe,jw)
            call prodel_conf(3,jp,mpe,jw,lr0,0,idl-1)
          enddo
          mpe=jj_sub(4,jp)
          if(mpe.ne.0) then
c            wt = vdint(lr0,lr0)+2.d0*voint(lr0,lr0)     !idl=4 hnil=2
            jw=iy(4,jp)
c           call prodel(3,wt,jp,mpe,jw)
            call prodel_conf(3,jp,mpe,jw,lr0,0,3)
          endif
        enddo
      enddo
      return
      end


      subroutine config_dbl()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      common/config/ndr,nwalk(0:max_orb)
!      common/casrst/ja(max_node),jb(max_node),jm(0:max_node)
!     :    ,jj(4,0:max_node),kk(0:max_node),no(0:max_innorb)
!     :    ,jv,jd(8),jt(8),js(8)
!      common/sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),
!     :     jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)
      data dzero/0.d0/
      if(norb_dbl.eq.0) return
      do ipae=1,25
        if(nu_ae(ipae).eq.0) cycle
        iwdownv=iw_downwei(1,ipae)
        do iwa=0,iwdownv-1
c       zz=' doub_800_v'
          iwad=iwalk_ad(1,ipae,iwa,0)
c          call prodel(1,wt0,0,ipae,iwad)
          call prodel_conf(1,0,ipae,iwad,0,0,1)
        enddo
      enddo
c       jps=js(1)
      do 100 lr0=norb_frz+1,norb_dz
        mr0=mul_tab(lsm_inn(lr0),ns_sm)
        iwd=jud(lr0)
          jpad =1+mr0
          jpad1=jpad+24
c         wld=wt0-voint(lr0,lr0)-vdint(lr0,lr0)
c          wls=wld-voint(lr0,lr0)
          do ipae=1,25
            if(nu_ae(ipae).eq.0) cycle
            iwdownv=iw_downwei(jpad,ipae)
            do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpad,ipae,iwa,iwd)
c             call prodel(1,wld,0,ipae,iwad)
              call prodel_conf(1,0,ipae,iwad,lr0,0,2)     !d
            enddo
          enddo

          if(jb_sys.gt.0) then
            do ipae=1,25
            if(nu_ae(ipae).eq.0) cycle
            iwdownv=iw_downwei(jpad1,ipae)
            do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpad1,ipae,iwa,iwd)
c             call prodel(1,wld,0,ipae,iwad)
              call prodel_conf(1,0,ipae,iwad,lr0,0,5)     !dd
            enddo
          enddo
        endif
         jpad=17+ns_sm
          iwd=just(lr0,lr0)
          do ipae=1,25
          if(nu_ae(ipae).eq.0) cycle
          iwdownv=iw_downwei(jpad,ipae)
          do iwa=0,iwdownv-1
            iwad=iwalk_ad(jpad,ipae,iwa,iwd)
c         call prodel(1,wls,0,ipae,iwad)
            call prodel_conf(1,0,ipae,iwad,lr0,0,4)     !s
        enddo
        enddo

      if(lr0.eq.norb_dz) goto 100
c           wld0=wld

      do 200 lr=lr0+1,norb_dz
        mr=mul_tab(mr0,lsm_inn(lr))
        jpat=9+mr
        jpas=17+mr
        jpat1=jpat+24
        iws=just(lr0,lr)
        iwt=iws            !
c          wld=wld0-voint(lr,lr)-vdint(lr,lr)
          do ipae=1,25
            if(nu_ae(ipae).eq.0) cycle
            iwdownv=iw_downwei(jpat,ipae)
            do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpat,ipae,iwa,iwt)
c           call prodel(1,wld,0,ipae,iwad)
              call prodel_conf(1,0,ipae,iwad,lr0,lr,3)   !t
          enddo
            if(jb_sys.gt.1) then
              iwdownv=iw_downwei(jpat1,ipae)
              do iwa=0,iwdownv-1
                iwad=iwalk_ad(jpat1,ipae,iwa,iwt)         !tt
c             call prodel(1,wld,0,ipae,iwad)
              call prodel_conf(1,0,ipae,iwad,lr0,lr,6)
            enddo
          endif
            iwdownv=iw_downwei(jpas,ipae)
            do iwa=0,iwdownv-1
              iwad=iwalk_ad(jpas,ipae,iwa,iws)
c           call prodel(1,wld,0,ipae,iwad)
              call prodel_conf(1,0,ipae,iwad,lr0,lr,4)   !s
          enddo
            if(jb_sys.gt.0) then
            iws1=just(lr,lr0)
              do iwa=0,iwdownv-1
                iwad=iwalk_ad(jpas,ipae,iwa,iws1)
c             call prodel(1,wld,0,ipae,iwad)
              call prodel_conf(1,0,ipae,iwad,lr0,lr,7)      !ss
            enddo
          endif
        enddo
200     continue
100   continue
      return
      end

      subroutine config_ext()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      common/config/ndr,nwalk(0:max_orb)
!      common/casrst/ja(max_node),jb(max_node),jm(0:max_node)
!     :    ,jj(4,0:max_node),kk(0:max_node),no(0:max_innorb)
!     :    ,jv,jd(8),jt(8),js(8)
!      common/sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),
!     :     jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)
      character*10 cc
        cc=' out_800_d'
        jws0=0
      do mra=1,8
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
c         call prodel(2,wld,0,ipae,jw)
          call prodel_conf(2,0,ipae,jw,la,0,2)      !d
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
        call prodel_conf(2,0,18,jweis,la,0,4)      !s
      enddo
      do im=1,8
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
c            wlt=wls
c           wls=wls+voint(lrb,lra)+vdint(lrb,lra)   ! w0=-1/2  w1=-3/2
c            wlt=wlt-voint(lrb,lra)+vdint(lrb,lra)   ! w0=-1/2  w1=1/2

c          call prodel(2,wls,0,ipas,jws)
c           call prodel(2,wlt,0,ipat,jwt)
            call prodel_conf(2,0,ipas,jws,la,lb,4)      !s
            call prodel_conf(2,0,ipat,jwt,la,lb,3)      !t
600       continue
          enddo
      enddo
      return
      end

!  idb=1  in dbl_space      ity_up=0-5             jd_type,jd_im,iwd
!  idb=2  in ext_space      ity_down=0-3          je_type,je_im,iwe
!  idb=3  in act_space      ity_up=0-5,itdown=0,3      jp ,mpe,iwa
      subroutine prodel_conf(idb,mg1,mg2,mg3,lr01,lr02,jpty)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      common/config/ndr,nwalk(0:max_orb)
!      common/sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),
!     :     jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)

      goto(100,200,300),idb
! in dbl_space
! jpty=  1,  2,  3,  4,  5,  6,  7
!        v   d   t   s  dd  tt  ss
100   ipae=mg2
      iwad=mg3
      lr1=norb_all-lr01+1
      lr2=norb_all-lr02+1
        isegdownwei=iseg_downwei(ipae)
        do mm=iwad+1,iwad+isegdownwei
          if(mm.eq.ndr) then
         goto(1,2,3,4,5,6,7),jpty
1         goto 1000
2         nwalk(lr1)=2
          goto 1000
3         nwalk(lr1)=2
          nwalk(lr2)=2
          goto 1000
4         nwalk(lr1)=2
          nwalk(lr2)=1
          if(lr02.eq.0) then
            nwalk(lr1)=0
          endif
          if(lr01.eq.0) then
            nwalk(lr2)=0
          endif
          goto 1000
5         nwalk(lr1)=1
          goto 1000
6         nwalk(lr1)=1
          nwalk(lr2)=1
          goto 1000
7         nwalk(lr1)=1
          nwalk(lr2)=2
          goto 1000
          endif
        enddo
      goto 1000

! in ext_space
! jpty=  1,  2,  3,  4
!        v   d   t   s
 200  ipae=mg2
      iwe=mg3
      lr1=lr01
      lr2=lr02
      do jdbl=1,mxnode
        if(nu_ad(jdbl).eq.0) cycle
        iw=iw_downwei(jdbl,ipae)
        iwupwei=jpad_upwei(jdbl)
        do iwa=0,iw-1
          do iwd=0,iwupwei-1
            mm=iwalk_ad(jdbl,ipae,iwa,iwd)+iwe
            if(mm.eq.ndr) then
              goto(10,20,30,40),jpty
10            goto 1000
20            nwalk(lr1)=1
              goto 1000
30            nwalk(lr1)=1
              nwalk(lr2)=1
              goto 1000
40            nwalk(lr1)=2
              nwalk(lr2)=1
              if(lr02.eq.0) then
                nwalk(lr1)=3
              endif
              goto 1000
            endif
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
      lr1=norb_all-lr01+1
      do jwu=jph+1,jph+in
        iwa=jw+ihy(jwu)-1
        do jwd=1,lwnu
          iwa=iwa+1
          do iwd=0,iwupwei-1
            iwad=iwalk_ad(jpad,ipae,iwa,iwd)
            do iwe=1,isegdownwei
              mm=iwe+iwad
              if(mm.eq.ndr) then
              nwalk(lr1)=jpty
                goto 1000
              endif
            enddo
          enddo
        enddo
      enddo
1000  return
      end

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
c complete double occpied space loops
      subroutine dbl_space_loop()
#include "drt_h.fh"
      if(norb_dbl.eq.0) return
      call dbl_space_loop_ijkk_sgezero()
      call dbl_space_loop_ijkl_sgezero()

      return
      end

      subroutine dbl_space_loop_ijkk_sgezero()
#include "drt_h.fh"
#include "intsort_h.fh"
      data dzero/0.d0/
c =============================  g1,2,4,6,7,8 ==========================
c       zz=' doub_800_v'
      dzero=0.0d0
      wls1=dzero
      db=jb_sys
      im=ns_sm
      jpds0=im+17
      do lri=norb_frz+1,norb_dz-1
      !  mi=lsm_inn(lri)
        iwdl=just(lri,lri)         !drr(03)-drr(30)
        do lrj=lri+1,norb_dz
          wl=voint(lrj,lri)          ! vdint(itailorb,iheadorb)
          iwdr=just(lrj,lrj)
          call prodab(1,0,jpds0,iwdl,iwdr,0,wl,0)
        enddo
      enddo
! bbs debug 20090722 jb_sys >0
      if(jb_sys.gt.0) then
!     drl(12)-drl(21)
        vlp_1=-sqrt(db*(db+2))/(db+1)
        do lri=norb_frz+1,norb_dz
          lmi=lsm_inn(lri)
          do lrj=lri+1,norb_dz
            lmj=lsm_inn(lrj)
            lmij=mul_tab(lmi,lmj)
            jpds=17+mul_tab(lmij,ns_sm)
            iwdl=just(lrj,lri)
            iwdr=just(lri,lrj)
            wl=-vlp_1*voint(lrj,lri)
!            print*, "lri,lrj",lri,lrj,jpds,iwdl,iwdr
            call prodab(1,0,jpds,iwdl,iwdr,0,wl,0)
          enddo
        enddo
      endif

      do 200 lri=norb_frz+1,norb_dz-1           !frz
        imi=lsm_inn(lri)
c       n2=ngw2(lri-2)
        do 201 lrj=lri+1,norb_dz
          mij=mul_tab(imi,lsm_inn(lrj))
          if(mij.ne.1) goto 201
          ni=mod(lrj-lri,2)
c=============== down comm for 2 4 =====================================
          vl_0=sqrt((db+2)/(db+1))
          vl_1=sqrt(db/(db+1))
          if(ni.eq.0)     vl_0=-vl_0
          if(ni.eq.0)     vl_1=-vl_1
          vl0_2=1.0d0
          vls0_2=1/(db+1)
          vls10_2=sqrt(db*(db+2))/(db+1)
      !    if(ni.eq.1)vl0_2=-vl0_2
          vls0=-0.5d0
          vls1=(db+3)/(2*db+2)
          vls1_a=(db-1)/(2*db+2)
          vls_c=db/(db+1)
          vlt0=-0.5d0
          vlt1=-0.5d0
          vls10_2b=sqrt(db*(db+2))/(db+1)
          if(ni.eq.1) then
           vl0_2=-vl0_2
            vls0_2=-vls0_2
            vls10_2=-vls10_2
            vls_c=-vls_c
          endif
          if(ni.eq.0) then
           vls0=-vls0
            vls1=-vls1
            vls1_a=-vls1_a
            vls10_2b=-vls10_2b
            vlt0=-vlt0
            vlt1=-vlt1
          endif
!        wl0=dzero
!         wl10=dzero
! for 2,4 the common
! ar(23)-c"(33)-ar(10) ar(13)-c"(33)-ar(20)
! ar(23)-drl(33)-ar(10)  ar(13)-drl(33)-ar(20)    310
          wltmp=dzero
         do lr=lri+1,lrj-1
            list =list3(lri,lrj,lr)
            wltmp=wltmp+2*vint_ci(list+1)-vint_ci(list)
!          wl0 =wl0 +vlt0*(2*vint_ci(list+1)-vint_ci(list))
!          wl10=wl10+vlt1*(2*vint_ci(list+1)-vint_ci(list))
!           wl0=wl0+vl0_1*(2*vint_ci(list+1)-vint_ci(list)) !310:neoc=2,
          enddo
          wl0 =vl_0*wltmp
          wl10=vl_1*wltmp
!     neoc(k)*wg43*(vint(list+2)+coe(k)*vint(list+1))
!    drl(33)-bl(23)-ar(10)
!    drl(33)-bl(13)-ar(20)
          wltmp=dzero
          do lr=1,lri-1
            list =list3(lri,lrj,lr)
            wltmp=wltmp+2*vint_ci(list+1)-vint_ci(list) ! 430 w0=-vl0 w1
          enddo
         wl0 =wl0+vl_0*wltmp
         wl10=wl10+vl_1*wltmp
!    ar(23)-bl(10)-drl(33)
!    ar(13)-bl(20)-drl(33)
          wltmp=dzero
          do lr=lrj+1,norb_dz
            list =list3(lri,lrj,lr)
            wltmp=wltmp+2*vint_ci(list+1)-vint_ci(list)  ! 220 w0=-vl0 w
          enddo
         wl0 =wl0+vl_0*wltmp
         wl10=wl10+vl_1*wltmp
c=============== start comm for 2 4 ====================================
          do 300 lrm=norb_frz+1,norb_dz        !ic=1,norb_act   !frz
           imm=lsm_inn(lrm)
            im=mul_tab(imm,imi)
            im=mul_tab(im,ns_sm)
           kij=0
            if(lrm.eq.lrj) kij=2
            if(lrm.eq.lri) kij=4
            if(lrm.lt.lri) kij=7
            if(lrm.gt.lrj) kij=6
            if(lrm.gt.lri.and.lrm.lt.lrj) kij=8
            goto(300,2,300,4,300,6,7,8),kij
!ss(1-c1)  ar(23)-ar(10)-    ar(13)-ar(20)
2           iwdl=just(lri,lrj)
            iwdr=just(lrj,lrj)
            list =list3(lri,lrj,lri)
            wl=wl0+vl_0*(voint(lri,lrj)+vint_ci(list))
            call prodab(1,0,jpds0,iwdl,iwdr,0,wl,0)
          if(jb_sys.gt.0) then
           wl=wl10+vl_1*(voint(lri,lrj)+vint_ci(list))
            iwdl=just(lrj,lri)
            call prodab(1,0,jpds0,iwdl,iwdr,0,wl,0)
          endif
           goto 300
!ss(1-c2)  ar(02)-ar(31)    ar(01)-ar(32)
4           iwdl=just(lri,lri)
            iwdr=just(lri,lrj)
            list =list3(lri,lrj,lrj)
            wl=wl0+vl_0*(voint(lri,lrj)+vint_ci(list))
            call prodab(1,0,jpds0,iwdl,iwdr,0,wl,0)
          if(jb_sys.gt.0) then
           iwdr=just(lrj,lri)
            wl=wl10+vl_1*(voint(lri,lrj)+vint_ci(list))
            call prodab(1,0,jpds0,iwdl,iwdr,0,wl,0)
          endif
           goto 300
!        ar(23)-br(13)-drr(30)  w0=-sqrt((db+2)/(db+1)) w1=0
!        ar(13)-br(23)-drr(30)  w0=-sqrt(db/(db+1)) w1=0
6            iwdl=just(lri,lrj)
             iwdr=just(lrm,lrm)
             list =list3(lri,lrj,lrm)
             wl=-vl_0*vint_ci(list)
             call prodab(1,0,jpds0,iwdl,iwdr,0,wl,0)
           if(jb_sys.gt.0) then
            iwdl=just(lrj,lri)
            wl=-vl_1*vint_ci(list)
            call prodab(1,0,jpds0,iwdl,iwdr,0,wl,0)
           endif
            goto 300
!           drr(03)br(32)br(31)   w0=-sqrt((db+2)/(db+1)) w1=0
!          drr(03)br(31)br(32)   w0=-sqrt(db/(db+1)) w1=0
7            iwdl=just(lrm,lrm)
             iwdr=just(lri,lrj)
             list =list3(lri,lrj,lrm)
             wl =-vl_0*vint_ci(list)
             call prodab(1,0,jpds0,iwdl,iwdr,0,wl,0)
           if(jb_sys.gt.0) then
            iwdr=just(lrj,lri)
             iwdl=just(lrm,lrm)
            wl =-vl_1*vint_ci(list)
             call prodab(1,0,jpds0,iwdl,iwdr,0,wl,0)
           endif
            goto 300
!  ar(23)-drl(30)-al(13)   w0=sqrt((db+2)/(db+1)) w1=0
!  ar(13)-drl(30)-al(23)   w0=sqrt(db/(db+1)) w1=0
8            iwdl=just(lri,lrj)
             iwdr=just(lrm,lrm)
             list =list3(lri,lrj,lrm)
             wl =-vl_0*vint_ci(list)
             call prodab(1,0,jpds0,iwdl,iwdr,0,wl,0)
           if(jb_sys.gt.0) then
!  ar(13)-drl(30)-al(23)   w0=sqrt(db/(db+1)) w1=0
            iwdl=just(lrj,lri)
            wl =-vl_1*vint_ci(list)
            call prodab(1,0,jpds0,iwdl,iwdr,0,wl,0)
           endif
c=============== start  d9(ss) d35(tt) =================================
             jpds=17+im
             jpdt=9+im
             jpdt1=jpdt+24
            iwls=just(lri,lrm)
             iwrs=just(lrm,lrj)
             iwlt=iwls     !
             iwrt=iwrs    !
!t    ar(23)-c'(22)-cw(33)-ar(32)   w0=1 w1=0
!t    ar(13)-c'(11)-cw(33)-ar(31)   w0=1 w1=0
!s    ar(23)-c'(12)-cw(33)-ar(31)   w0=-1/(db+1) w1=0
!s    ar(23)-c'(11)-cw(33)-ar(32)   w0=-sqrt(db*(db+2))/(db+1) w1=0
!s    ar(13)-c'(21)-cw(33)-ar(32)   w0=1/(db+1) w1=0
!s    ar(13)-c'(22)-cw(33)-ar(31)   w0=-sqrt(db*(db+2))/(db+1) w1=0
             wlt0=dzero
             wls0=dzero
             wls0_1=dzero                   ! bbs_tmp
             wls0_2=dzero
             wltmp=dzero
            do lr=lrm+1,lrj-1
               list=list3(lri,lrj,lr)
               wltmp=wltmp+2*vint_ci(list+1)-vint_ci(list)
!         wlt0=wlt0+vl0_2*(2*vint_ci(list+1)-vint_ci(list))
             enddo
             wlt0=wlt0+vl0_2*wltmp
             wls0=wls0-vls0_2*wltmp
             wls0_1=wls0_1-vls10_2*wltmp
             wls0_2=wls0_2+vls0_2*wltmp
!       wls2=wls0
!t    ar(23)-cw(33)-c'(22)-ar(32)      w0=1 w1=0
!t    ar(13)-cw(33)-c'(11)-ar(31)      w0=1 w1=0
!s    ar(23)-cw(33)-c'(12)-ar(31)      w0=-1/(db+1) w1=0
!s    ar(23)-cw(33)-c'(11)-ar(32)      w0=-sqrt(db*(db+2))/(db+1) w1=0
!s    ar(13)-cw(33)-c'(21)-ar(32)      w0=1/(db+1) w1=0
!s    ar(13)-cw(33)-c'(22)-ar(31)      w0=-sqrt(db*(db+2))/(db+1) w1=0
             wltmp=dzero
             do lr=lri+1,lrm-1
               list=list3(lri,lrj,lr)
               wltmp=wltmp+2*vint_ci(list+1)-vint_ci(list)
!         wlt0=wlt0+vl0_2*(2*vint_ci(list+1)-vint_ci(list))
             enddo
             wlt0=wlt0+vl0_2*wltmp
             wls0=wls0-vls0_2*wltmp
             wls0_1=wls0_1-vls10_2*wltmp
             wls0_2=wls0_2+vls0_2*wltmp
!t    drl(33)-bl(23)-c'(22)-ar(32) w0=-1 w1=0
!t    drl(33)-bl(13)-c'(11)-ar(31) w0=-1 w1=0
!s    drl(33)-bl(23)-c'(12)-ar(31) w0=1/(db+1)  w1=0
!s    drl(33)-bl(23)-c'(11)-ar(32) w0=sqrt(db*(db+2))/(db+1) w1=0
!s    drl(33)-bl(13)-c'(21)-ar(32) w0=-1/(db+1)  w1=0
!s    drl(33)-bl(13)-c'(22)-ar(31) w0=sqrt(db*(db+2))/(db+1)   w1=0
             wltmp=dzero
             do lr=1,lri-1
               list=list3(lri,lrj,lr)
               wltmp=wltmp+(vint_ci(list)-2*vint_ci(list+1))
             enddo
             wlt0=wlt0-vl0_2*wltmp                  !bbs_tmp
             wls0=wls0+vls0_2*wltmp
             wls0_1=wls0_1+vls10_2*wltmp
             wls0_2=wls0_2-vls0_2*wltmp
!t    ar(23)-c'(22)-bl(32)-drl(33) w0=-1 w1=0
!t    ar(13)-c'(11)-bl(31)-drl(33) w0=-1 w1=0
!s    ar(23)-c'(12)-br(31)-drl(33) w0=1/(db+1)  w1=0
!s    ar(23)-c'(11)-br(32)-drl(33) w0=sqrt(db*(db+2))/(db+1) w1=0
!s    ar(13)-c'(21)-br(32)-drl(33) w0=-1/(db+1)  w1=0
!s    ar(13)-c'(22)-br(31)-drl(33) w0=sqrt(db*(db+2))/(db+1) w1=0
             wltmp=dzero
             do lr=lrj+1,norb_dz
               list=list3(lri,lrj,lr)
               wltmp=wltmp+(vint_ci(list)-2*vint_ci(list+1))
             enddo
!            wlt0=wlt0+vl0_2*wltmp
!            wls0=wls0-vls0_2*wltmp
!            wls0_1=wls0_1+vls10_2*wltmp
!            wls0_2=wls0_2+vls0_2*wltmp
             wlt0=wlt0-vl0_2*wltmp
             wls0=wls0+vls0_2*wltmp
             wls0_1=wls0_1+vls10_2*wltmp
             wls0_2=wls0_2-vls0_2*wltmp
      !       wls0=-wlt0
c=============== up comm for 9 35 ======================================
!t   ar(23)-c'(22)-arw(32) w0=1 w1=0
!t   ar(13)-c'(11)-arw(31) w0=1 w1=0
!s   ar(23)-c'(12)-arw(31) w0=-1/(db+1)  w1=0
!s   ar(23)-c'(11)-arw(31) w0=-sqrt(db*(db+2))/(db+1) w1=0
!s   ar(13)-c'(21)-arw(32) w0=1/(db+1)  w1=0
!s   ar(13)-c'(22)-arw(31) w0=-sqrt(db*(db+2))/(db+1) w1=0
!             list=list3(lri,lrj,lrj)
!             wlt=wlt0+vl0_2*vint_ci(list)
!             wls=wls0-vl0_2*vint_ci(list)
             list=list3(lri,lrj,lrj)
             wlt=wlt0+vl0_2*vint_ci(list)
             wls=wls0-vls0_2*vint_ci(list)
             wls1=wls0_1-vls10_2*vint_ci(list)
             wls2=wls0_2+vls0_2*vint_ci(list)
!     wl_bbs=-vls0_2*vint_ci(list)
!t   ar(23)-cw(22)-ar(32)
!t   ar(13)-cw(11)-ar(31)
!s   ar(23)-cw(12)-ar(31)    w0=-1/(db+1)    ar(12)-drr(22)-ar(32)
!s   ar(23)-cw(11)-ar(32)    w0=-sqrt(db*(db+1))/(db+1)
!s   ar(13)-cw(21)-ar(32)    w0=1/(db+1)
!s   ar(13)-cw(22)-ar(31)  w0=-sqrt(db*(db+1))/(db+1)
             list=list3(lri,lrj,lrm)
             wlt=wlt+vl0_2*vint_ci(list+1)
             wls=wls-vls0_2*(vint_ci(list+1)-(db+2)*vint_ci(list)) !coe=
             wls1=wls1-vls10_2*(vint_ci(list+1)-vint_ci(list))  !coe=-1
             wls2=wls2+vls0_2*(vint_ci(list+1)+db*vint_ci(list))  !coe=d
!             wls2=wls2+vls0_2*(vint_ci(list+1)-2*vint_ci(list))
!     wl_bbs=wl_bbs-vls0_2*(vint_ci(list+1)-2*vint_ci(list))
!t   arw(23)-c'(22)-ar(32)
!t   arw(13)-c'(11)-ar(31)
!s   arw(23)-c'(12)-ar(31)
!s   arw(23)-c'(11)-ar(32)
!s   arw(13)-c'(21)-ar(32)
!s   arw(13)-c'(22)-ar(31)
!t   ar(23)-c'(22)-ar(32)    w0=1
!t   ar(13)-c'(11)-ar(31)    w0=1
!s   ar(23)-c'(12)-ar(31)   w0=-1/(db+1)  w1=0
!s   ar(23)-c'(11)-ar(32)   w0=-sqrt(db*(db+2))/(db+1) w1=0
!s   ar(13)-c'(21)-ar(32)   w0=1/(db+1)  w1=0
!s   ar(13)-c'(22)-ar(31)   w0=-sqrt(db*(db+2))/(db+1) w1=0
             list=list3(lri,lrj,lri)
             wlt=wlt+vl0_2*(voint(lri,lrj)+ vint_ci(list))
             wls=wls-vls0_2*(voint(lri,lrj)+ vint_ci(list))
             wls1=wls1-vls10_2*(voint(lri,lrj)+ vint_ci(list))
             wls2=wls2+vls0_2*(voint(lri,lrj)+ vint_ci(list))
!     wl_bbs=wl_bbs-vls0_2*(voint(lri,lrj)+ vint_ci(list))

!ar(23)-drr(12)-ar(31)  w0=(db+2)/(db+1)   ???? complete
             call prodab(1,0,jpds,iwls,iwrs,0,wls,0)
             call prodab(1,0,jpdt,iwlt,iwrt,0,wlt,0)
             if(jb_sys.gt.0) then
               iwls=just(lri,lrm)
               iwrs=just(lrj,lrm)
               call prodab(1,0,jpds,iwls,iwrs,0,wls1,0)
               iwls=just(lrm,lri)
               iwrs=just(lrj,lrm)
               call prodab(1,0,jpds,iwls,iwrs,0,wls2,0)
               iwls=just(lrm,lri)
               iwrs=just(lrm,lrj)
               call prodab(1,0,jpds,iwls,iwrs,0,wls1,0)
             endif
            if(jb_sys.gt.1) then
               call prodab(1,0,jpdt1,iwlt,iwrt,0,wlt,0)    !bbs_tmp
             endif
300        continue
c=============== start  d10_ss(dd,ss,tt)  ==============================
!   ar(23)-drr(33)-ar(32)  ar(23)-cw(33)-ar(23)
!   ar(13)-drr(33)-ar(31)  ar(13)-cw(33)-ar(13) w0=1 w1=0
!      if(lri.eq.2.and.lrj.eq.4) then
!       write(6,*) "bbs_tmp err"
!     endif
           wl0=dzero
           wltmp=dzero
           do lr=lri+1,lrj-1
             list=list3(lri,lrj,lr)
             wltmp=wltmp+(2*vint_ci(list+1)-vint_ci(list))
           enddo
           wl0 =wl0-vl0_2*wltmp
!   drl(33)-bl(23)-ar(32)
!   drl(33)-bl(13)-ar(31)
           wltmp=dzero
           do lr=1,lri-1
             list=list3(lri,lrj,lr)
             wltmp=wltmp+(vint_ci(list)-2*vint_ci(list+1))
           enddo
           wl0 =wl0+vl0_2*wltmp
!   ar(23)-bl(32)-drl(33)
!   ar(13)-bl(31)-drl(33)
           wltmp=dzero
           do lr=lrj+1,norb_dz
             list=list3(lri,lrj,lr)
             wltmp=wltmp+(vint_ci(list)-2*vint_ci(list+1))
           enddo
           wl0 =wl0+vl0_2*wltmp
! ar-arw    arw-al   ar-ar
           list=list3(lri,lrj,lrj)
           wl0=wl0-vl0_2*vint_ci(list)
           list=list3(lri,lrj,lri)
           wl0=wl0-vl0_2*(voint(lri,lrj)+vint_ci(list))

           im=mul_tab(lsm_inn(lri),ns_sm)    !!!
           jpdd=1+im
          iwld=jud(lri)
           iwrd=jud(lrj)
           call prodab(1,0,jpdd,iwld,iwrd,0,wl0,0)
           if(jb_sys.gt.0) then
            jpdd=jpdd+24
             call prodab(1,0,jpdd,iwld,iwrd,0,wl0,0)
          endif
           do lr=lrj+1,norb_dz
             list =list3(lri,lrj,lr)
             wls=wl0+(vls0-vl0_2)*(vint_ci(list)-2*vint_ci(list+1))
     :         -vls1*vint_ci(list)
      !ar(23)-bl(32)-drl(22) ar(13)-bl(31)-drl(11)
             wlt=wl0+(vlt0-vl0_2)*(vint_ci(list)-2*vint_ci(list+1))
     :         -vlt1*vint_ci(list)
             im=mul_tab(lsm_inn(lri),lsm_inn(lr))
             im=mul_tab(im,ns_sm)
             jpds=17+im
             jpdt=9+im
             iwls=just(lri,lr)
             iwrs=just(lrj,lr)
             iwlt=iwls             !
             iwrt=iwrs             !
             call prodab(1,0,jpds,iwls,iwrs,0,wls,0)
             call prodab(1,0,jpdt,iwlt,iwrt,0,wlt,0)
             if(jb_sys.gt.0) then
               wls_a=wl0+(vls0-vl0_2)*(vint_ci(list)-2*vint_ci(list+1))
     :           -vls1_a*vint_ci(list)
               wls_b=-vls10_2b*vint_ci(list)
               iwls=just(lr,lri)
               iwrs=just(lr,lrj)
               call prodab(1,0,jpds,iwls,iwrs,0,wls_a,0)
               iwls=just(lri,lr)
               iwrs=just(lr,lrj)
               call prodab(1,0,jpds,iwls,iwrs,0,wls_b,0)
               iwls=just(lr,lri)
               iwrs=just(lrj,lr)
               call prodab(1,0,jpds,iwls,iwrs,0,wls_b,0)
            endif
             if(jb_sys.gt.1) then
               jpdt=jpdt+24
               call prodab(1,0,jpdt,iwlt,iwrt,0,wlt,0)
             endif
           enddo
c======= start d5(ss),d40(tt) =================================
           do lr=norb_frz+1,lri-1
             list =list3(lri,lrj,lr)
      ! drl(22)-bl(13)-ar(31)
             wls=wl0+(vls0-vl0_2)*(vint_ci(list)-2*vint_ci(list+1))
     :             -vls1*vint_ci(list)
      ! drl(22)-bl(23)-ar(32)  drl(11)-bl(13)-ar(31)
             wlt=wl0+(vlt0-vl0_2)*(vint_ci(list)-2*vint_ci(list+1))
     :             -vlt1*vint_ci(list)
             im=mul_tab(lsm_inn(lri),lsm_inn(lr))
             im=mul_tab(im,ns_sm)
             jpds=17+im                     !bbs_tmp
             jpdt=9+im
             iwls=just(lr,lri)
             iwrs=just(lr,lrj)
             iwlt=iwls             !
             iwrt=iwrs             !
             call prodab(1,0,jpds,iwls,iwrs,0,wls,0)
             call prodab(1,0,jpdt,iwlt,iwrt,0,wlt,0)
             if(jb_sys.gt.0) then
!drl(11)-bl(23)-ar(32)
               wls=wl0+(vls0-vl0_2)*(vint_ci(list)-2*vint_ci(list+1))
     :             -vls1_a*vint_ci(list)
               iwls=just(lri,lr)
               iwrs=just(lrj,lr)
               call prodab(1,0,jpds,iwls,iwrs,0,wls,0)
!drl(12)-bl(23)-ar(31) drl(21)-bl(13)-ar(32)
               wls=-vls10_2b*vint_ci(list)
               iwls=just(lri,lr)
               iwrs=just(lr,lrj)
               call prodab(1,0,jpds,iwls,iwrs,0,wls,0)
               iwls=just(lr,lri)
               iwrs=just(lrj,lr)
               call prodab(1,0,jpds,iwls,iwrs,0,wls,0)
             endif
             if(jb_sys.gt.1) then
               jpat=jpdt+24
               call prodab(1,0,jpat,iwlt,iwrt,0,wlt,0)
             endif
           enddo
c======= end g5,40 =================================
201      continue
200    continue
      continue
      return
      end

      subroutine dbl_space_loop_ijkl_sgezero()
#include "drt_h.fh"
#include "intsort_h.fh"
c =============================  g11,12  == (v-s)=======================
c =============================  g41,42  == (v-t)=======================
      wls1=0.d0
      wls2=0.d0
      db=jb_sys
      do 10 lrl=norb_frz+1,norb_dz-3
         iml=lsm_inn(lrl)
         do 20 lrk=lrl+1,norb_dz-2
           imk=lsm_inn(lrk)
           do 30 lrj=lrk+1,norb_dz-1
            imj=lsm_inn(lrj)
             do 40 lri=norb_dz,lrj+1,-1
               imi=lsm_inn(lri)
c               list=list4(lri,lrj,lrk,lrl)
               list=list4(lrl,lrk,lrj,lri)
               ni =mod(lrk-lrl+lri-lrj,2)
               imik=mul_tab(imi,imk)
             if(imik.eq.mul_tab(imj,iml)) then
               im=mul_tab(imik,ns_sm)
               jpds=17+im
               jpdt=9+im
               jpdt1=jpdt+24
               iwls=just(lrl,lrj)
               iwrs=just(lrk,lri)
               iwlt=iwls             !
               iwrt=iwrs             !
              if(jb_sys.eq.0) then
c                   w0g11=-1/2 w1g11=-3/2                 === g11 ===  1
                 wls =vint_ci(list+1)+vint_ci(list+2)
c                   w0g42=-1/2 w1g42=1/2                  === g42 ===  1
                 wlt =vint_ci(list+1)-vint_ci(list+2)
               endif
              if(jb_sys.gt.0) then
!                  wog11=-1/2 w1g11=-(db+3)/(2*db+2)
!      ar(23)bl(32)bl(13)ar(31)
                 w1=-(db+3)/(2*db+2)
                 wls =vint_ci(list+1)+(-0.5d0-w1)*vint_ci(list+2)
!                  wog42=-1/2 w1g42=1/2     ar(23)bl(32)bl(23)ar(32)
               wlt =vint_ci(list+1)-vint_ci(list+2)
!                  w0=0 w1=-sqrt(db*(db+2)/(db+1))
!      ar(13)bl(32)bl(23)ar(31)     ar(23)bl(31)bl(13)ar(32)
                 w1=-sqrt(db*(db+2))/(db+1)
               wls1=-w1*vint_ci(list+2)
!      ar(13)bl(31)bl(23)ar(32)      w0=-1/2
                 w1=-(db-1)/(2*db+2)
               wls2=vint_ci(list+1)+(-0.5d0-w1)*vint_ci(list+2)
!      ar(13)bl(31)bl(13)ar(31)
!                 w0=-1/2 w1=1/2
!               wlt1=wlt
              endif
!            if(ni.eq.1)       wls=-wls
               if(ni.eq.1)  then
                 wls=-wls
                 wlt=-wlt
                 wls1=-wls1
                 wls2=-wls2
               endif
             call prodab(1,0,jpds,iwls,iwrs,0,wls,0)
               call prodab(1,0,jpdt,iwlt,iwrt,0,wlt,0)
               if(jb_sys.gt.0) then
                 iwrs1=just(lrk,lri)
                 iwls1=just(lrj,lrl)
                 call prodab(1,0,jpds,iwls1,iwrs1,0,wls1,0)
                 iwrs1=just(lri,lrk)
                 iwls1=just(lrl,lrj)
                 call prodab(1,0,jpds,iwls1,iwrs1,0,wls1,0)
                 iwrs1=just(lri,lrk)
                 iwls1=just(lrj,lrl)
                 call prodab(1,0,jpds,iwls1,iwrs1,0,wls2,0)
               endif
               if(jb_sys.gt.1) then
                  call prodab(1,0,jpdt1,iwlt,iwrt,0,wlt,0)
             endif
              endif
            imil=mul_tab(imi,iml)
            if(imil.eq.mul_tab(imj,imk)) then
              im=mul_tab(imil,ns_sm)
              jpds=17+im
              jpdt=9+im
              jpdt1=jpdt+24
              iwls=just(lrl,lri)
            iwrs=just(lrk,lrj)
              iwlt=iwls       !
              iwrt=iwrs       !
              if(jb_sys.eq.0) then
c                   w0g11=-1/2 w1g11=-3/2                  === g11 ===
                wls =vint_ci(list)+vint_ci(list+1)
c                   w0g42=-1/2 w1g42=1/2                   === g42 ===
                wlt =vint_ci(list+1)-vint_ci(list)
              endif
              if(jb_sys.gt.0) then
!                  wog11=-1/2 w1g11=-(db+3)/(2*db+2)
!      ar(23)bl(32)br(31)al(13)
                 w1=-(db+3)/(2*db+2)
                 wls =(-0.5d0-w1)*vint_ci(list)+vint_ci(list+1)
!                  wog42=-1/2 w1g42=1/2
!      ar(23)bl(32)br(23)al(32)
               wlt =vint_ci(list+1)-vint_ci(list)
!                  w0=0 w1=-sqrt(db*(db+2)/(db+1))
!      ar(23)bl(31)br(32)al(13)
                 w1=-sqrt(db*(db+2))/(db+1)
               wls1=-w1*vint_ci(list)
!      ar(13)bl(31)br(32)al(23)
!                w0=-1/2 w1=-(db-1)/(2*db+2)
                 w1=-(db-1)/(2*db+2)
               wls2=(-0.5d0-w1)*vint_ci(list)+vint_ci(list+1)
              endif
               if(ni.eq.1) then
               wls=-wls
                 wlt=-wlt
                 wls1=-wls1
                 wls2=-wls2
             endif
             call prodab(1,0,jpds,iwls,iwrs,0,wls,0)
               call prodab(1,0,jpdt,iwlt,iwrt,0,wlt,0)
               if(jb_sys.gt.0) then
                 iwls=just(lri,lrl)
                 iwrs=just(lrk,lrj)
                 call prodab(1,0,jpds,iwls,iwrs,0,wls1,0)
                 iwls=just(lrl,lri)
                 iwrs=just(lrj,lrk)
               call prodab(1,0,jpds,iwls,iwrs,0,wls1,0)
               iwls=just(lri,lrl)
               call prodab(1,0,jpds,iwls,iwrs,0,wls2,0)
               endif
               if(jb_sys.gt.1) then
                 call prodab(1,0,jpdt1,iwlt,iwrt,0,wlt,0)
               endif
              endif
              imij=mul_tab(imi,imj)
             if(imij.eq.mul_tab(imk,iml)) then
                im=mul_tab(imij,ns_sm)
                jpds=17+im
                jpdt=9+im
                jpdt1=24+jpdt
                iwls=just(lrl,lrk)
                iwrs=just(lrj,lri)
                iwlt=iwls             !
                iwrt=iwrs             !
              if(jb_sys.eq.0) then
c           w0g12=1, w1g12=0                                === g12 ===
                  wls =vint_ci(list)+vint_ci(list+2)
c           w0g41=0, w1g41=1                                === g41 ===
                  wlt =vint_ci(list)-vint_ci(list+2)
              endif
                if(jb_sys.gt.0) then
                  w0=(db+2)/(2*db+2)
                  w1=db/(2*db+2)
                  wls =(w0+w1)*vint_ci(list)+(w0-w1)*vint_ci(list+2)
                  wlt =vint_ci(list)-vint_ci(list+2)
               w0=sqrt(db*(db+2))/(2*db+2)
                  w1=-w0
                  wls1=(w0+w1)*vint_ci(list)+(w0-w1)*vint_ci(list+2)
                  w0=db/(2*db+2)
                  w1=(db+2)/(2*db+2)
                  wls2=(w0+w1)*vint_ci(list)+(w0-w1)*vint_ci(list+2)
              endif
                if(ni.eq.1) then
                wls=-wls
                  wlt=-wlt
                  wls1=-wls1
                  wls2=-wls2
              endif
              call prodab(1,0,jpds,iwls,iwrs,0,wls,0)
                call prodab(1,0,jpdt,iwlt,iwrt,0,wlt,0)
                if(jb_sys.gt.0) then
                  iwls=just(lrl,lrk)
                  iwrs=just(lri,lrj)
                  call prodab(1,0,jpds,iwls,iwrs,0,wls1,0)
                  iwls=just(lrk,lrl)
                  iwrs=just(lrj,lri)
                  call prodab(1,0,jpds,iwls,iwrs,0,wls1,0)
      !         iwls=
               iwrs=just(lri,lrj)
                  call prodab(1,0,jpds,iwls,iwrs,0,wls2,0)
                endif
                if(jb_sys.gt.1) then
                  call prodab(1,0,jpdt1,iwlt,iwrt,0,wlt,0)
              endif
              endif
40          continue
30        continue
20      continue
10    continue
      return
      end

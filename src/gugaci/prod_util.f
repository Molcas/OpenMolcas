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
      subroutine jl_ne_jr(mp,jl,jr,jwl,jwr,lopu)
#include "drt_h.fh"
#include "pl_structure_h.fh"
!-------------------------------------------------
! on entry:
!     jl  - left DRT node
!     jr  - right DRT node
!     jwl - left weight
!     jwr - right weight
! output
!     mp   - number of partial loop tails
!     lopu(1,i) - left weight of parital loop i
!     lopu(2,i) - right weight of partial loop i
!     lopu(3,i) - loop tail node
!-------------------------------------------------
      dimension lopu(4,loputmp),lopi(4,loputmp),lopj(4,loputmp)
      !lopu(1,*)=jwl lopu(2,*)=jwr,lopu(3,*)=jl,lopu(4,*)=jr
      !write(6,*) "in subroutine jl_ne_jr",jl,jr
      mp=0
      lpi=1
      lopi(1,1)=jwl
      lopi(2,1)=jwr
      lopi(3,1)=jl
      lopi(4,1)=jr
103   lpj=0
      do 114 nlp=1,lpi
        if(lopi(3,nlp).eq.lopi(4,nlp)) then
          mp=mp+1
          lopu(1:4,mp)=lopi(1:4,nlp)
          cycle
        endif
        jlp=lopi(3,nlp)
        jrp=lopi(4,nlp)
        do 115 idlr=1,4
          ml=jjl_sub(idlr,jlp)
          mr=jj_sub(idlr,jrp)
          !write(6,*) "ml,mr",ml,mr
          if(ml.eq.0.or.mr.eq.0) goto 115
          jwlp=lopi(1,nlp)
          jwrp=lopi(2,nlp)
          if(idlr.ne.1)jwlp=jwlp+iyl(idlr,jlp)
          if(idlr.ne.1)jwrp=jwrp+iy(idlr,jrp)
          lpj=lpj+1
          lopj(1,lpj)=jwlp
          lopj(2,lpj)=jwrp
          lopj(3,lpj)=ml
          lopj(4,lpj)=mr
115     continue
114   continue
      if(lpj.eq.0) goto 200
      lpi=lpj
      do i=1,lpj
        lopi(1:4,i)=lopj(1:4,i)
      enddo
      goto 103
200   return
      end

!  idb=1  in dbl_space      ity_up=0-5               0 ,jpad,iwdl,iwdr,
!  idb=2  in act_space      ity_up=0-5,itdown=0,3      jph, jpe,iwal,iwa
!  idb=3  betwin dbl and act   ity_up=0-5,itdown=0,3      jpe,iwdl,iwdr,
c**************************************************************
      subroutine prodab_h0(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr)
c**************************************************************
c     calculate h matrix elements which were contributed
c     only by whole inner space loops and store them into
c     vector2
c
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      dimension lopu(4,loputmp)
c     write(6,*) 'prodab_02 '

      goto(100,200,300),idb
! in dbl_space
c100   jpad=mg2
100   iwdl=mg3
      iwdr=mg4
      ipae=1
      jpad=1
      mntmp=0
      iwdown=iw_downwei(jpad,ipae)
      lwnu=iseg_downwei(ipae)
      do iwa=0,iwdown-1
        iwadl=iwalk_ad(jpad,ipae,iwa,iwdl)
        iwadr=iwalk_ad(jpad,ipae,iwa,iwdr)
        mm=iwadl
        nn=iwadr
        do m=1,lwnu
          mm=mm+1
          nn=nn+1
          if(mm.gt.nn) mntmp=mm*(mm-1)/2+nn
          if(nn.gt.mm) mntmp=nn*(nn-1)/2+mm
          vector2(mntmp)=vector2(mntmp)+wl
!            if(mntmp.eq.2) then
!                write(6,*)'  102',vector2(mntmp),wl
!           endif
        enddo
      enddo
c      enddo
      goto 1000
! in act_space
200   if(jpad.ne.jpadl) return
      jph=mg1
      jpl=mg2
c     iwal=mg3
c     iwar=mg4
      iwupwei=jpad_upwei(jpad)
      isegdownwei=iseg_downwei(ipae)
      jpy=jphy(jph)
      in=ihy(jpy)

      call jl_ne_jr(mp,jpl,jpr,mg3,mg4,lopu)
      do lp=1,mp
        iwl=lopu(1,lp)-1
        iwr=lopu(2,lp)-1
        jpe=lopu(3,lp)
        lwnu=iy(1,jpe)
        do jwu=jpy+1,jpy+in
          iwal=iwl+ihy(jwu)
          iwar=iwr+ihy(jwu)
          do jwd=1,lwnu
            iwal=iwal+1
            iwar=iwar+1
            do iwd=0,iwupwei-1
              iwadl=iwalk_ad(jpadl,ipael,iwal,iwd)
              iwadr=iwalk_ad(jpad,ipae,iwar,iwd)
              do iwe=1,isegdownwei
                mm=iwadl+iwe
                nn=iwadr+iwe
                if(mm.gt.nn) then
                  mntmp=mm*(mm-1)/2+nn
                else
                  mntmp=nn*(nn-1)/2+mm
                endif
                vector2(mntmp)=vector2(mntmp)+wl
                if(mntmp.eq.7) then
                  write(6,*)'  202',vector2(mntmp),wl
                endif
              enddo
            enddo
          enddo
        enddo
      enddo
      goto 1000
! betwin act and dbl
300   jpl =mg1
      iwdl=mg2
      iwdr=mg3
      isegdownwei=iseg_downwei(ipae)

      call jl_ne_jr(mp,jpl,jpr,mg4,mg5,lopu)
      do lp=1,mp
        iwal=lopu(1,lp)-1
        iwar=lopu(2,lp)-1
        jpe= lopu(3,lp)
        jwnu=iy(1,jpe)
        do ii=1,jwnu
          iwal=iwal+1
          iwar=iwar+1
          mm=iwalk_ad(jpadl,ipael,iwal,iwdl)
          nn=iwalk_ad(jpad,ipae,iwar,iwdr)
          do iwe=1,isegdownwei
            mm=mm+1                  ! iwl=iwalk_ad
            nn=nn+1                  ! iwl=iwalk_ad
              if(mm.gt.nn) then
                mntmp=mm*(mm-1)/2+nn
              else
                mntmp=nn*(nn-1)/2+mm
              endif
              vector2(mntmp)=vector2(mntmp)+wl
!              if(mntmp.eq.2) then
!                write(6,*)'  302',vector2(mntmp),wl
!           endif
          enddo
        enddo
      enddo
      goto 1000
1000  return
      end

********************************************************************
      subroutine prodab_h(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr)
c*******************************************************************
c     whole inner space loops - h*c
c     26 feb 2007 - revised by suo bing for multi-root calculation
c
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      dimension lopu(4,loputmp)
!     write(6,*) 'prodab_02 ',mcroot,indx(1),iw_downwei(jpad,ipae)
      goto(100,200,300),idb
! in dbl_space
100   jpad=mg2
      iwdl=mg3
      iwdr=mg4
      ipaeend=25
      do ipae=1,ipaeend
        if(nu_ae(ipae).eq.0) cycle
        iwdown=iw_downwei(jpad,ipae)
        if(iwdown.eq.0) cycle
        lwnu=iseg_downwei(ipae)
        do iwa=0,iwdown-1
          iwadl=iwalk_ad(jpad,ipae,iwa,iwdl)
          iwadr=iwalk_ad(jpad,ipae,iwa,iwdr)
          do irot=1,mcroot
            irtidx=indx(irot)
            mm=iwadl+irtidx
            nn=iwadr+irtidx
            do m=1,lwnu
              mm=mm+1
              nn=nn+1
!              if(mm.gt.nci_dim.or.nn.gt.nci_dim) then
!                print*,
!     *            jpad,ipae,iw_downwei(jpad,ipae),iseg_downwei(ipae),
!     *            jpad_upwei(jpad),iw_sta(jpad,ipae),iwadl,iwadr
!              endif
              vector2(mm)=vector2(mm)+vector1(nn)*wl
              vector2(nn)=vector2(nn)+vector1(mm)*wl
            enddo
          enddo
        enddo
      enddo
      goto 1000

! in act_space
200   if(jpad.ne.jpadl) return
      jph=mg1
      jpl=mg2
c     iwal=mg3
c     iwar=mg4
      iwupwei=jpad_upwei(jpad)
      isegdownwei=iseg_downwei(ipae)
      jpy=jphy(jph)
      in=ihy(jpy)

      call jl_ne_jr(mp,jpl,jpr,mg3,mg4,lopu)
      do lp=1,mp
        iwl=lopu(1,lp)-1
        iwr=lopu(2,lp)-1
        jpe=lopu(3,lp)
        lwnu=iy(1,jpe)
        do jwu=jpy+1,jpy+in
          iwal=iwl+ihy(jwu)
          iwar=iwr+ihy(jwu)
          do jwd=1,lwnu
            iwal=iwal+1
            iwar=iwar+1
            do iwd=0,iwupwei-1
              iwadl=iwalk_ad(jpadl,ipael,iwal,iwd)
              iwadr=iwalk_ad(jpad,ipae,iwar,iwd)
              do irot=1,mcroot
                irtidx=indx(irot)
                mm=iwadl+irtidx
                nn=iwadr+irtidx
                do iwe=1,isegdownwei
                  mm=mm+1
                  nn=nn+1
                  vector2(mm)=vector2(mm)+vector1(nn)*wl
                  vector2(nn)=vector2(nn)+vector1(mm)*wl
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      goto 1000

! betwin act and dbl
300   jpl =mg1
      iwdl=mg2
      iwdr=mg3
      isegdownwei=iseg_downwei(ipae)

      call jl_ne_jr(mp,jpl,jpr,mg4,mg5,lopu)

      do lp=1,mp
        iwal=lopu(1,lp)-1
        iwar=lopu(2,lp)-1
        jpe= lopu(3,lp)
        jwnu=iy(1,jpe)
        do ii=1,jwnu
          iwal=iwal+1
          iwar=iwar+1
          iwadl=iwalk_ad(jpadl,ipael,iwal,iwdl)
          iwadr=iwalk_ad(jpad,ipae,iwar,iwdr)
          do irot=1,mcroot
            irtidx=indx(irot)
            mm=iwadl+irtidx
            nn=iwadr+irtidx
            do iwe=1,isegdownwei
              mm=mm+1                  ! iwl=iwalk_ad
              nn=nn+1                  ! iwl=iwalk_ad
              vector2(mm)=vector2(mm)+vector1(nn)*wl
              vector2(nn)=vector2(nn)+vector1(mm)*wl
            enddo
          enddo
        enddo
      enddo
      goto 1000
1000  return
      end

      subroutine prodab(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr)
#include "drt_h.fh"

      select case (log_prod)
        case (1)
          call prodab_h(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr)
        case (2)  ! save H0 CI matrix, use for mrci
          call prodab_h0(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr)
        case (3)
          call prodab_h0_t(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr)
        case (4)  ! H0C, use for mrpt2
          call prodab_h0_d(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr)
        case default ! HC, use for mrci
      end select

      return
      end

      subroutine prodab_h0_d(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      dimension lopu(4,loputmp)
!log_prod=2:directly no_formh0
      character*16 loop_type
      loop_type=' prod_h0 '
c       write(*,*) 'prodab_h0 '

      goto(100,200,300),idb
! in dbl_space
100   continue
      iwdl=mg3
      iwdr=mg4
      iwdown=iw_downwei(jpad,ipae)
      lwnu=jpae_downwei(ipae)
      do iwa=0,iwdown-1
        iwadl=iwalk_ad(jpad,ipae,iwa,iwdl)
        iwadr=iwalk_ad(jpad,ipae,iwa,iwdr)
        do irot=1,mcroot
          irtidx=indx(irot)
          mm=iwadl+irtidx
          nn=iwadr+irtidx
          do m=1,lwnu
            mm=mm+1
            nn=nn+1
!            if(mm.gt.nci_dim.or.nn.gt.nci_dim) then
!              print*,
!     *          jpad,ipae,iw_downwei(jpad,ipae),iseg_downwei(ipae),
!     *          jpad_upwei(jpad),iw_sta(jpad,ipae),iwadl,iwadr
!            endif
            vector2(mm)=vector2(mm)+vector1(nn)*wl
            vector2(nn)=vector2(nn)+vector1(mm)*wl
          enddo
        enddo
        !do m=1,lwnu
        !  mm=mm+1
        !  nn=nn+1
        !  vector2(mm)=vector2(mm)+vector1(nn)*wl
        !  vector2(nn)=vector2(nn)+vector1(mm)*wl
        !enddo
      enddo
      goto 1000
! in act_space
200   if(jpad.ne.jpadl) return
      jph=mg1
      jpl=mg2
c     iwal=mg3
c     iwar=mg4
      iwupwei=jpad_upwei(jpad)
      jpaedownwei=jpae_downwei(ipae)
      jpy=jphy(jph)
      in=ihy(jpy)

      call jl_ne_jr(mp,jpl,jpr,mg3,mg4,lopu)
      do lp=1,mp
        iwl=lopu(1,lp)-1
        iwr=lopu(2,lp)-1
        jpe=lopu(3,lp)
        lwnu=iy(1,jpe)
        do jwu=jpy+1,jpy+in
          iwal=iwl+ihy(jwu)
          iwar=iwr+ihy(jwu)
          do jwd=1,lwnu
            iwal=iwal+1
            iwar=iwar+1
            do iwd=0,iwupwei-1
              iwadl=iwalk_ad(jpadl,ipael,iwal,iwd)
              iwadr=iwalk_ad(jpad,ipae,iwar,iwd)
              do irot=1,mcroot
                irtidx=indx(irot)
                mm=iwadl+irtidx
                nn=iwadr+irtidx
                do iwe=1, jpaedownwei
                  mm=mm+1
                  nn=nn+1
                  vector2(mm)=vector2(mm)+vector1(nn)*wl
                  vector2(nn)=vector2(nn)+vector1(mm)*wl
                enddo
              enddo
              !do iwe=1,jpaedownwei
              !  mm=iwadl+iwe
              !  nn=iwadr+iwe
              !  vector2(mm)=vector2(mm)+vector1(nn)*wl
              !  vector2(nn)=vector2(nn)+vector1(mm)*wl
              !enddo
            enddo
          enddo
        enddo
      enddo
      goto 1000
! between act and dbl
300   jpl =mg1
      iwdl=mg2
      iwdr=mg3
      jpaedownwei=jpae_downwei(ipae)

      call jl_ne_jr(mp,jpl,jpr,mg4,mg5,lopu)
      do lp=1,mp
        iwal=lopu(1,lp)-1
        iwar=lopu(2,lp)-1
        jpe= lopu(3,lp)
        jwnu=iy(1,jpe)
        do ii=1,jwnu
          iwal=iwal+1
          iwar=iwar+1
          iwadl=iwalk_ad(jpadl,ipael,iwal,iwdl)
          iwadr=iwalk_ad(jpad,ipae,iwar,iwdr)
          do irot=1,mcroot
            irtidx=indx(irot)
            mm=iwadl+irtidx
            nn=iwadr+irtidx
            do iwe=1,jpaedownwei
              mm=mm+1                  ! iwl=iwalk_ad
              nn=nn+1                  ! iwl=iwalk_ad
              vector2(mm)=vector2(mm)+vector1(nn)*wl
              vector2(nn)=vector2(nn)+vector1(mm)*wl
            enddo
          enddo
        enddo
      enddo
1000  return
      end

      subroutine prodab_h0_t(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      dimension lopu(4,loputmp)
!log_prod=1:truanditional formh0
      character*16 loop_type
      loop_type=' prod_h0 '

c       write(*,*) 'prodab_h0 '

      goto(100,200,300),idb
! in dbl_space
100   continue
      iwdl=mg3
      iwdr=mg4
      iwdown=iw_downwei(jpad,ipae)
      lwnu=jpae_downwei(ipae)
      do iwa=0,iwdown-1
          mm=iwalk_ad(jpad,ipae,iwa,iwdl)
          nn=iwalk_ad(jpad,ipae,iwa,iwdr)
          do m=1,lwnu
            mm=mm+1
            nn=nn+1
            mnh0=mm*(mm-1)/2+nn
            if(nn.gt.mm) mnh0=nn*(nn-1)/2+mm
             vector2(mnh0)=vector2(mnh0)+wl
!            if((mm.eq.1.and.nn.eq.9).or.(mm.eq.9.and.nn.eq.1)) then
!           write(nf2,'(a10,2i5,2f18.8)')' in dbl ',mm,nn,wl,vector2(mnh
!             endif
         enddo
       enddo
      goto 1000
! in act_space
200   if(jpad.ne.jpadl) return
      jph=mg1
      jpl=mg2
c     iwal=mg3
c     iwar=mg4
      iwupwei=jpad_upwei(jpad)
      jpaedownwei=jpae_downwei(ipae)
      jpy=jphy(jph)
      in=ihy(jpy)

       call jl_ne_jr(mp,jpl,jpr,mg3,mg4,lopu)
      do lp=1,mp
         iwl=lopu(1,lp)-1
         iwr=lopu(2,lp)-1
         jpe=lopu(3,lp)
        lwnu=iy(1,jpe)
        do jwu=jpy+1,jpy+in
          iwal=iwl+ihy(jwu)
           iwar=iwr+ihy(jwu)
          do jwd=1,lwnu
             iwal=iwal+1
             iwar=iwar+1
             do iwd=0,iwupwei-1
               iwadl=iwalk_ad(jpadl,ipael,iwal,iwd)
               iwadr=iwalk_ad(jpad,ipae,iwar,iwd)
               do iwe=1,jpaedownwei
                 mm=iwadl+iwe
                 nn=iwadr+iwe
                mnh0=mm*(mm-1)/2+nn
                if(nn.gt.mm) mnh0=nn*(nn-1)/2+mm
                 vector2(mnh0)=vector2(mnh0)+wl
!             if((mm.eq.1.and.nn.eq.9).or.(mm.eq.9.and.nn.eq.1)) then
!           write(nf2,'(a10,2i5,2f18.8)')' in act ',mm,nn,wl,vector2(mnh
!          endif
                     enddo
             enddo
           enddo
         enddo
       enddo
      goto 1000
! between act and dbl
300       jpl =mg1
      iwdl=mg2
       iwdr=mg3
      jpaedownwei=jpae_downwei(ipae)

      call jl_ne_jr(mp,jpl,jpr,mg4,mg5,lopu)
      do lp=1,mp
         iwal=lopu(1,lp)-1
         iwar=lopu(2,lp)-1
         jpe= lopu(3,lp)
         jwnu=iy(1,jpe)
         do ii=1,jwnu
           iwal=iwal+1
           iwar=iwar+1
           mm=iwalk_ad(jpadl,ipael,iwal,iwdl)
          nn=iwalk_ad(jpad,ipae,iwar,iwdr)
           do iwe=1,jpaedownwei
             mm=mm+1                  ! iwl=iwalk_ad
             nn=nn+1                  ! iwl=iwalk_ad
             mnh0=mm*(mm-1)/2+nn
             if(nn.gt.mm) mnh0=nn*(nn-1)/2+mm
             vector2(mnh0)=vector2(mnh0)+wl
!             if((mm.eq.1.and.nn.eq.9).or.(mm.eq.9.and.nn.eq.1)) then
!       write(nf2,'(a10,2i5,2f18.8)')' act-dbl ',mm,nn,wl,vector2(mnh0)
!                 endif
           enddo
         enddo
      enddo
1000  return
      end

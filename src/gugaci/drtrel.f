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
c drt related subroutines
c
      subroutine active_drt()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "files_gugaci.fh"
!      common/casrst/ja(max_node),jb(max_node),jm(0:max_node)
!     :    ,jj(4,0:max_node),kk(0:max_node),no(0:max_innorb)
!     :    ,jv,jd(8),jt(8),js(8)
      dimension iin(0:max_node)
      nci_dim=0
      if(norb_act.ne.0) goto 100
!====================  norb_act=0 ========================
      iseg_sta(1)=0
      iseg_dim(1)=1
      do im=1,ng_sm
        jdim=nu_ae(1+im)
        jtim=nu_ae(9+im)
        jsim=nu_ae(17+im)
        jd(im)=jdim
        jt(im)=jtim
        js(im)=jsim
        iseg_dim(jdim)=iseg_downwei(jdim)*jpad_upwei(jdim)
        iseg_dim(jtim)=iseg_downwei(jtim)*jpad_upwei(jtim)
        iseg_dim(jsim)=iseg_downwei(jsim)*jpad_upwei(jsim)
        if(iseg_dim(jdim).eq.0) then
          jd(im)=0
          nu_ad(jdim)=0
          nu_ae(jdim)=0
        endif
        if(iseg_dim(jtim).eq.0) then
          jt(im)=0
          nu_ad(jtim)=0
          nu_ae(jtim)=0
        endif
        if(iseg_dim(jsim).eq.0) then
          js(im)=0
          nu_ad(jsim)=0
          nu_ae(jsim)=0
        endif
      enddo
      do jp=2,25
        iseg_sta(jp)=iseg_sta(jp-1)+iseg_dim(jp-1)
        enddo
      nci_dim=iseg_sta(25)+iseg_dim(25)
      iseg_sta(26)=nci_dim
      do jp=1,25
        if(iseg_downwei(jp).eq.0) cycle
        iseg_upwei(jp)=iseg_dim(jp)/iseg_downwei(jp)
      enddo
      goto 200
!====================  norb_act<>0 ========================
!100   if(logic_mr)      call rst(ndd,indd)    !npp=2
!      if(logic_mrelcas) call rcas(ndd,indd)   !npp=3

100   write(6,*)' '
      write(6,*) ' now reading distinct row tableau'
      call readdrt(ludrt)
c
      nu_ae(1)=jv
      do im=1,ng_sm
        nu_ae(1+im)=jd(im)
        nu_ae(9+im)=jt(im)
        nu_ae(17+im)=js(im)
      enddo
      jds=1
      jde=mxnode
      jps=no(norb_inn)+1
      jpe=no(norb_inn+1)

      jp=jv
c     iin(1:jpe)=0
c     iin(0)=0
      iin(:)=0
      iin(jp)=1
      iseg_sta(1)=0       ! for node_ext
      iseg_dim(1)=0
      do jpn=jpe,1,-1
        do i=1,4
          jji=jj(i,jpn)
          if(iin(jji).eq.0) cycle
          iin(jpn)=iin(jpn)+iin(jji)
        enddo
      enddo
      do jdn=jds,jde
        if(nu_ad(jdn).eq.0) cycle
        ndi=iin(jdn)*iseg_downwei(1)*jpad_upwei(jdn)
        iseg_dim(1)=iseg_dim(1)+ndi
      enddo
      do im=1,ng_sm
        jp=jd(im)
        iseg_sta(1+im)=nci_dim       ! for node_d_ext
        iseg_dim(1+im)=0
        if(jp.eq.0) cycle
        iin(1:jpe)=0
        iin(0)=0
        iin(jp)=1
        do jpn=jpe,1,-1
          do i=1,4
            jji=jj(i,jpn)
            if(iin(jji).eq.0) cycle
            iin(jpn)=iin(jpn)+iin(jji)
          enddo
        enddo
        do jdn=jds,jde
          if(nu_ad(jdn).eq.0) cycle
          ndi=iin(jdn)*iseg_downwei(1+im)*jpad_upwei(jdn)
          iseg_dim(1+im)=iseg_dim(1+im)+ndi
        enddo
      enddo

      do im=1,ng_sm
        jp=jt(im)
        iseg_sta(9+im)=nci_dim       ! for node_t_ext
        iseg_dim(9+im)=0
        if(jp.eq.0) cycle
        iin(1:jpe)=0
        iin(0)=0
        iin(jp)=1
        do jpn=jpe,1,-1
          do i=1,4
            jji=jj(i,jpn)
            if(iin(jji).eq.0) cycle
            iin(jpn)=iin(jpn)+iin(jji)
          enddo
        enddo
        do jdn=jds,jde
          if(nu_ad(jdn).eq.0) cycle
          ndi=iin(jdn)*iseg_downwei(9+im)*jpad_upwei(jdn)
          iseg_dim(9+im)=iseg_dim(9+im)+ndi
        enddo
      enddo
      do im=1,ng_sm
        jp=js(im)
        iseg_sta(17+im)=nci_dim       ! for node_s_ext
        iseg_dim(17+im)=0
        if(jp.eq.0) cycle
        iin(1:jpe)=0
        iin(0)=0
        iin(jp)=1
        do jpn=jpe,1,-1
          do i=1,4
            jji=jj(i,jpn)
            if(iin(jji).eq.0) cycle
            iin(jpn)=iin(jpn)+iin(jji)
          enddo
        enddo
        do jdn=jds,jde
          if(nu_ad(jdn).eq.0) cycle
          ndi=iin(jdn)*iseg_downwei(17+im)*jpad_upwei(jdn)
          iseg_dim(17+im)=iseg_dim(17+im)+ndi
        enddo
      enddo
      do jp=2,25
        iseg_sta(jp)=iseg_sta(jp-1)+iseg_dim(jp-1)
c        write(nf2,'(3i10)')jp-1,iseg_sta(jp-1),iseg_dim(jp-1)      !to
      enddo
c        write(nf2,'(3i10)')25,iseg_sta(25),iseg_dim(25)          !to de
      nci_dim=iseg_sta(25)+iseg_dim(25)
      iseg_sta(26)=nci_dim
      do jp=1,25
        if(iseg_downwei(jp).eq.0) cycle
        iseg_upwei(jp)=iseg_dim(jp)/iseg_downwei(jp)
      enddo
c      call get_ivy()

! to the end of dbl,act,ext parts
200   call dbl_downwalk()
      write(6,*)'  end of drt,nci_dim= ',nci_dim
      !write(6,*)'-----------------------------------------------'
      !write(6,*)'         ci drt-information'
      !write(6,*)'       num.of orbitals:       ',nst
      !write(6,*)'       num.of froz-orbitals:  ',nzst
      !write(6,*)'       num.of active-orbitals:',ndh
      !write(6,*)'       num.of electrons:      ',nel
      !write(6,*)'       multiplicity:          ',mult
      !write(6,*)'       num.of configurations: ',ndime
      !write(6,*)'       symmetry:              ',zm
      !write(6,*)'-----------------------------------------------'
      return
      end

      subroutine rst(id,indd)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "files_gugaci.fh"
!      common/casrst/ja(max_node),jb(max_node),jm(0:max_node)
!     :    ,jj(4,0:max_node),kk(0:max_node),no(0:max_innorb)
!     :    ,jv,jd(8),jt(8),js(8)
      write(6,*)' '
      write(6,*) ' now reading distinct row tableau'
      call readdrt(ludrt)
c      open(21,file="fort.drt",form="unformatted")
c      read(21) id
c      write(6,*) " id=",id
c      read(21) ja(1:id),jb(1:id),jm(1:id)
c      read(21) jj(1:4,0:id)
c      read(21) kk(0:id)
c      read(21) no(0:norb_inn+1)
c      read(21) jv,jd(1:8),jt(1:8),js(1:8)
c      close(21)

      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(id)
        call Unused_integer(indd)
      end if
      end

      subroutine ref_gfs(nel,ndj,locu,nm)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
!      common/casrst/ja(max_node),jb(max_node),jm(0:max_node)
!     :    ,jj(4,0:max_node),kk(0:max_node),no(0:max_innorb)
!     :    ,jv,jd(8),jt(8),js(8)
      dimension lhsm(8),locu(8,max_ref),lscu(0:8,max_ref)
      nl_act=norb_act
      ne_act=nel-2*norb_dz
      ne_s=nint(spin*2)
      lhs=nstart_act
      lhe=norb_inn
      lhsm(1:8)=0
      do lh=lhs,lhe
          lm=lsm_inn(lh)
        lhsm(lm)=lhsm(lm)+1
      enddo
      mdj=0
      do 500 nes=ne_s,ne_act,2
        do 100 l1=0,lhsm(1)
          do 101 l2=0,lhsm(2)
            do 102 l3=0,lhsm(3)
              do 103 l4=0,lhsm(4)
                do 104 l5=0,lhsm(5)
                  do 105 l6=0,lhsm(6)
                    do 106 l7=0,lhsm(7)
                      do 107 l8=0,lhsm(8)
                        lpsum=l1+l2+l3+l4+l5+l6+l7+l8
                        if(lpsum.ne.nes) goto 107
                        mys=1
                        if(mod(l1,2).eq.1)mys=mul_tab(mys,1)
                        if(mod(l2,2).eq.1)mys=mul_tab(mys,2)
                        if(mod(l3,2).eq.1)mys=mul_tab(mys,3)
                        if(mod(l4,2).eq.1)mys=mul_tab(mys,4)
                        if(mod(l5,2).eq.1)mys=mul_tab(mys,5)
                        if(mod(l6,2).eq.1)mys=mul_tab(mys,6)
                        if(mod(l7,2).eq.1)mys=mul_tab(mys,7)
                        if(mod(l8,2).eq.1)mys=mul_tab(mys,8)
                        if(mys.ne.nm) goto 107
                        mdj=mdj+1
                        lscu(0,mdj)=lpsum
                        lscu(1,mdj)=l1
                        lscu(2,mdj)=l2
                        lscu(3,mdj)=l3
                        lscu(4,mdj)=l4
                        lscu(5,mdj)=l5
                        lscu(6,mdj)=l6
                        lscu(7,mdj)=l7
                        lscu(8,mdj)=l8
107                   continue
106                 continue
105               continue
104             continue
103           continue
102         continue
101       continue
100     continue
500   continue
      ndj=0
      do 300 m=1,mdj
        npair=(ne_act-lscu(0,m))/2
        do 200 l1=0,lhsm(1)-lscu(1,m)
          do 201 l2=0,lhsm(2)-lscu(2,m)
            do 202 l3=0,lhsm(3)-lscu(3,m)
              do 203 l4=0,lhsm(4)-lscu(4,m)
                do 204 l5=0,lhsm(5)-lscu(5,m)
                  do 205 l6=0,lhsm(6)-lscu(6,m)
                    do 206 l7=0,lhsm(7)-lscu(7,m)
                      do 207 l8=0,lhsm(8)-lscu(8,m)
                        lpsum=l1+l2+l3+l4+l5+l6+l7+l8
                        if(lpsum.eq.npair) then
                          m1=l1*2+lscu(1,m)
                          m2=l2*2+lscu(2,m)
                          m3=l3*2+lscu(3,m)
                          m4=l4*2+lscu(4,m)
                          m5=l5*2+lscu(5,m)
                          m6=l6*2+lscu(6,m)
                          m7=l7*2+lscu(7,m)
                          m8=l8*2+lscu(8,m)
                          do 600 ldj=1,ndj
                            if(m1.ne.locu(1,ldj)) goto 600
                            if(m2.ne.locu(2,ldj)) goto 600
                            if(m3.ne.locu(3,ldj)) goto 600
                            if(m4.ne.locu(4,ldj)) goto 600
                            if(m5.ne.locu(5,ldj)) goto 600
                            if(m6.ne.locu(6,ldj)) goto 600
                            if(m7.ne.locu(7,ldj)) goto 600
                            if(m8.ne.locu(8,ldj)) goto 600
                            goto 207
600                       continue
                          ndj=ndj+1
                          locu(1,ndj)=m1
                          locu(2,ndj)=m2
                          locu(3,ndj)=m3
                          locu(4,ndj)=m4
                          locu(5,ndj)=m5
                          locu(6,ndj)=m6
                          locu(7,ndj)=m7
                          locu(8,ndj)=m8
                        endif
207                   continue
206                 continue
205               continue
204             continue
203           continue
202         continue
201       continue
200     continue
300   continue

      do nre=1,ndj
      write(6,'(5x,i6,8i3)')nre,(locu(i,nre),i=1,8)
      enddo
      return
      end

      subroutine rcas(id,indd)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "files_gugaci.fh"
!      common/casrst/ja(max_node),jb(max_node),jm(0:max_node)
!     :    ,jj(4,0:max_node),kk(0:max_node),no(0:max_innorb)
!     :    ,jv,jd(8),jt(8),js(8)
      write(6,*)' '
      write(6,*) ' now reading distinct row tableau'
      call readdrt(ludrt)

!      print*, "bbs debug rcas,kk(27)",kk(27)

c      open(21,file="fort.drt",form="unformatted")
c      read(21) id
c      read(21) ja(1:id),jb(1:id),jm(1:id)
c      read(21) jj(1:4,0:id)
c      read(21) kk(0:id)
c      read(21) no(0:norb_inn+1)
c      read(21) jv,jd(1:8),jt(1:8),js(1:8)
c      close(21)
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(id)
        call Unused_integer(indd)
      end if

      end

      subroutine check_rcas3(jk,ind,inb,ndj,locu)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
!     common/casrst/ja(max_node),jb(max_node),jm(0:max_node)
!     :    ,jj(4,0:max_node),kk(0:max_node),no(0:max_innorb)
!     :    ,jv,jd(8),jt(8),js(8)
      dimension ind(8,max_node),lsym(8),iexcit(ndj),locu(8,ndj)
      inb=0
      nsumel=0
      do i=1,8
        lsym(i)=ind(i,jk)
        nsumel=nsumel+lsym(i)
      enddo
      do i=1,ndj
        iexcit(i)=0
        do m=1,8
          iex=lsym(m)-locu(m,i)
         if(iex.gt.0) then
           iexcit(i)=iexcit(i)+iex
          endif
        enddo
      enddo
      inb=iexcit(1)
      do i=2,ndj
         inb=min(inb,iexcit(i))
      enddo
      inb=inb+ja(jk)*2+jb(jk)
      return
      end

      subroutine irfrst(iselcsf_occ)
#include "drt_h.fh"
#include "pl_structure_h.fh"
c ifrno(j)=i
c irfno(i)=j no. i ref is no. j cfs in h0
      common/config/ndr,nwalk(0:max_orb)
      dimension iselcsf_occ(max_innorb,max_ref)
      dimension iwalktmp(max_orb)
      logical log_exist
      icsfwlk=0
      ndimh0=nci_h0 !iw_sta(2,1)
      icount=0
      do i=1,n_ref
        do j=1,ndimh0
          call found_a_config(j,1.0,0)
          do im=1,norb_all
            iwalktmp(im)=nwalk(norb_all-im+1)
          enddo

          log_exist=.false.
          ij=norb_dz
          do ii=1,norb_act
            ij=ij+1
            icsfwlk=iwalktmp(ij)
            if(icsfwlk.eq.3) then
              icsfocc=2
            elseif(icsfwlk.eq.2) then
              icsfocc=1
            elseif(icsfwlk.eq.1) then
              icsfocc=1
            elseif(icsfwlk.eq.0) then
              icsfocc=0
            else
              icsfocc=-1
            endif
            if(icsfocc.ne.iref_occ(ij,i)) then
              goto 30
            endif
          enddo
          log_exist=.true.
          icount=icount+1
          irfno(icount)=j
          ifrno(j)=icount
c          write(6,2000) iwalktmp(norb_dz+1:norb_dz+norb_act)
c          write(6,*) "icount",icount,j
30        continue
        enddo
      enddo

      irf=icount
      write(6,3000) irf

      return
c1000  format(1x,"warnning!the selected csf is not in references states")
c2000  format(1x,"the selected csf is :",2x,32(i1))
3000  format(1x,"number of gelfand states in referance space:",1x,i4)
c Avoid unused argument warnings
      if (.false.) call Unused_integer_array(iselcsf_occ)
c...end of irfrst
      end

      subroutine irfrst_bak(iselcsf_occ)
#include "drt_h.fh"
#include "pl_structure_h.fh"
      common/config/ndr,nwalk(0:max_orb)
      dimension iselcsf_occ(max_innorb,max_ref)
      dimension iwalktmp(max_orb)
      logical log_exist
      nocc=0
      do i=1,mroot
        log_exist=.false.
        do ire=1,n_ref
          ij=norb_dz
          do io=1,norb_act
            ij=ij+1
            if(iselcsf_occ(io,i).eq.3) nocc=2
            if(iselcsf_occ(io,i).eq.2) nocc=1
            if(iselcsf_occ(io,i).eq.1) nocc=1
            if(iselcsf_occ(io,i).eq.0) nocc=0
            if(nocc.ne.iref_occ(io+norb_dz,ire)) then
               goto 10
            endif
          enddo
          log_exist=.true.
          goto 20
10        continue
        enddo
20      continue
        if(.not.log_exist) then
          write(6,1000)
          write(6,2000) iselcsf_occ(1:norb_act,i)
          write(6,*) " please select this state as reference state"
c          stop 777
        endif
      enddo
      icsfocc=0
      ndimh0=nci_h0 !iw_sta(2,1)
      icount=0
      do i=1,n_ref
        do j=1,ndimh0
          call found_a_config(j,1.0,0)
          do im=1,norb_all
            iwalktmp(im)=nwalk(norb_all-im+1)
          enddo

          log_exist=.false.
          ij=norb_dz
          do ii=1,norb_act
            ij=ij+1
            icsfwlk=iwalktmp(ij)
            if(icsfwlk.eq.3) icsfocc=2
            if(icsfwlk.eq.2) icsfocc=1
            if(icsfwlk.eq.1) icsfocc=1
            if(icsfwlk.eq.0) icsfocc=0
            if(icsfocc.ne.iref_occ(ij,i)) then
              goto 30
            endif
          enddo
          log_exist=.true.
          icount=icount+1
          irfno(icount)=j
          ifrno(j)=icount
c          write(6,2000) iwalktmp(norb_dz+1:norb_dz+norb_act)
c          write(6,*) "icount",icount,j
30        continue
        enddo
      enddo

      irf=icount
      do i=1,2*mroot
        iocsf=mjn(i)
        log_exist=.false.
        do j=1,icount
          if(iocsf.eq.irfno(j)) then
            log_exist=.true.
            goto 40
          endif
        enddo
40      continue
        if(.not.log_exist) then
          irf=irf+1
          irfno(irf)=iocsf
          ifrno(iocsf)=irf
        endif
      enddo

      write(6,3000) irf

      return
1000  format(1x,"warnning!the selected csf is not in references states")
2000  format(1x,"the selected csf is :",2x,32(i1))
3000  format(1x,"number of gelfand states in referance space:",1x,i4)
c...end of irfrst
      end

      function min_itexcit(indjk)
      common/ref/ndj,ndjgrop,ndjmod
      dimension indjk(4),itexcit(ndj)
c      integer*4 indjk  =  00 00 00 00 00 00 00 00 00 00  00 00 00 00 00
!    indexcit=  ir1 ir2 ir3 ir4 ir5 ir6 ir7 ir8 ......... ir15
!--------------------------------------------------------
      min_itexcit=3
      nj=0
      do ngrop=1,ndjgrop-1         !1-(ndjgrop-1) grop
        indexcit=indjk(ngrop)
        do lref=1,15
          ixcit=ishft(indexcit,-2*(lref-1))
          ixcit=mod(ixcit,4)
!          if(ixcit.ne.0) ixcit=ixcit-1
!          if(ixcit.eq.0) ixcit=3
            min_itexcit=min(min_itexcit,ixcit)
          if(min_itexcit.eq.0) return
          itexcit(nj+lref)=ixcit
        enddo
        nj=nj+15
      enddo
        indexcit=indjk(ndjgrop)       !last grop
        do lref=1,ndjmod
          ixcit=ishft(indexcit,-2*(lref-1))
          ixcit=mod(ixcit,4)
!          if(ixcit.ne.0) ixcit=ixcit-1
!          if(ixcit.eq.0) ixcit=3
          min_itexcit=min(min_itexcit,ixcit)
          if(min_itexcit.eq.0) return
          itexcit(nj+lref)=ixcit
        enddo
!--------------------------------------------------------
      return
      end

      subroutine njexcit(idcc,indjk,locuk0,n_ref)
      common/ref/ndj,ndjgrop,ndjmod
      dimension indjk(4),locuk0(n_ref),itexcit(n_ref)
      nj=0
      do ngrop=1,ndjgrop-1         !1-(ndjgrop-1) grop
        indexcit=indjk(ngrop)
        indjk(ngrop)=0
        do lref=1,15
          ixcit=ishft(indexcit,-2*(lref-1))
          ixcit=mod(ixcit,4)
         if(idcc-locuk0(nj+lref).eq.1) ixcit=ixcit+1
         if(idcc-locuk0(nj+lref).eq.2) ixcit=ixcit+2
         if(ixcit.ge.3) ixcit=3
         itexcit(nj+lref)=ixcit
          indjk(ngrop)=ishft(ixcit,2*(lref-1))+indjk(ngrop)
        enddo
        nj=nj+15
      enddo
        indexcit=indjk(ndjgrop)       !last grop
        indjk(ndjgrop)=0
        do lref=1,ndjmod
          ixcit=ishft(indexcit,-2*(lref-1))
          ixcit=mod(ixcit,4)
         if(idcc-locuk0(nj+lref).eq.1) ixcit=ixcit+1
         if(idcc-locuk0(nj+lref).eq.2) ixcit=ixcit+2
         if(ixcit.ge.3) ixcit=3
         itexcit(nj+lref)=ixcit
         indjk(ngrop)=ishft(ixcit,2*(lref-1))+indjk(ngrop)
       enddo

      return
      end

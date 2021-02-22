!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      subroutine gugadrt_rcas(id,indd)
#include "gendrt.fh"
#include "Sysdrt.fh"
#include "stdalloc.fh"
#include "casrst_drt.fh"
      dimension locu(8,max_ref),jc(max_node)
      dimension noh(max_innorb),itm(0:max_node)
      allocatable :: ind(:,:),iwy(:,:)
      call mma_allocate(ind,8,max_node,label='ind')
      call mma_allocate(iwy,[1,4],[0,max_node],label='iwy')
      write(6,*)' '
      write(6,*) 'now generate distinct row tableau'
      noh=0
      itm=0
      iwy=0
      ind=0
      locu=0
      jc=0

      nel =n_electron
      nm  =ns_sm
      no(1:norb_dz-1)=0
      j=0
      ja0=ja_sys
      jb0=jb_sys
      jc0=jc_sys
!  v_node
      ja(1)=ja0
      jb(1)=jb0
      jc(1)=jc0
      kk(1)=norb_dz
      do i=1,8
        ind(i,1)=0
        enddo
        jm(1)=nm
!  d_node
      imd=0
      do node=2,9
        imd=imd+1
        if(nu_ad(node).eq.0) cycle
        ja(node)=ja(1)
        jb(node)=jb0+1
        jc(node)=jc0-1
        kk(node)=norb_dz
        do i=1,8
          ind(i,node)=0
        enddo
        jm(node)=imd
      enddo
!  t_node
      imt=0
      do node=10,17
        imt=imt+1
        if(nu_ad(node).eq.0) cycle
        ja(node)=ja(1)
        jb(node)=jb0+2
        jc(node)=jc0-2
        kk(node)=norb_dz
        do i=1,8
          ind(i,node)=0
        enddo
        jm(node)=imt
      enddo
!  s_node
      ims=0
      do node=18,25
        ims=ims+1
        if(nu_ad(node).eq.0) cycle
        ja(node)=ja(1)+1
        jb(node)=jb0
        jc(node)=jc0-1
        kk(node)=norb_dz
        do i=1,8
          ind(i,node)=0
        enddo
        jm(node)=ims
      enddo
      if(jroute_sys.gt.1) then
!  d'_node
        imd=0
        do node=26,33
          imd=imd+1
          if(nu_ad(node).eq.0) cycle
          ja(node)=ja(1)+1
          jb(node)=jb0-1
          jc(node)=jc0
          kk(node)=norb_dz
          do i=1,8
           ind(i,node)=0
          enddo
          jm(node)=imd
        enddo
      endif
      if(jroute_sys.gt.2) then
!  t'_node
        imt=0
        do node=34,41
          imt=imt+1
          if(nu_ad(node).eq.0) cycle
          ja(node)=ja(1)+2
          jb(node)=jb0-2
          jc(node)=jc0
          kk(node)=norb_dz
          do i=1,8
           ind(i,node)=0
          enddo
          jm(node)=imt
        enddo
      endif
      no(norb_dz)=mxnode

      call gugadrt_ref_gfs(nel,ndj,locu,nm)


      j=0
      jk=mxnode
 7    j=j+1
      if(j.le.mxnode) then
        if(nu_ad(j).eq.0) goto 7
      endif
!      if(j.le.mxnode.and.nu_ad(j).eq.0) goto 7
      k0=kk(j)
      if(jc(j).eq.0) goto 1
      jk=jk+1
!                                            ***********
!                                            *   d=0   *
!                                            ***********
      ja(jk)=ja(j)
      jb(jk)=jb(j)
      jc(jk)=jc(j)-1
      jm(jk)=jm(j)
      do i=1,8
        ind(i,jk)=ind(i,j)
      enddo
      inb=0
      if(kk(j).eq.norb_inn-1) then
      call gugadrt_check_rcas3(jk,ind,inb,ndj,locu)
      endif
      if(inb.gt.2) then
        jk=jk-1
        jj(1,j)=0
        goto 11
      endif

      do 18 ii=no(k0)+1,jk-1
         if(ja(jk).ne.ja(ii)) goto 18
         if(jb(jk).ne.jb(ii)) goto 18
         if(jc(jk).ne.jc(ii)) goto 18
         if(jm(jk).ne.jm(ii)) goto 18
         if(kk(j).eq.norb_inn-1) goto 601
         do i=1,8
          if(ind(i,jk).ne.ind(i,ii)) goto 18
         enddo
 601     jk=jk-1
         jj(1,j)=ii
         goto 11
 18      continue
      jj(1,j)=jk
      kk(jk)=kk(j)+1
      if(kk(jk).ne.kk(jk-1)) then
      no(k0)=jk-1
      end if
 11   continue
!      v=0
 1    if(jb(j).eq.0) goto 2
      jk=jk+1
!                                            ***********
!                                            *   d=1   *
!                                            ***********
      ja(jk)=ja(j)
      jb(jk)=jb(j)-1
      jc(jk)=jc(j)
      jm(jk)=mul_tab(lsm_inn(kk(j)+1),jm(j))
      do i=1,8
        ind(i,jk)=ind(i,j)
      enddo
      ind(lsm_inn(kk(j)+1),jk)=ind(lsm_inn(kk(j)+1),j)+1
      inb=0
      if(kk(j).eq.norb_inn-1) then
      call gugadrt_check_rcas3(jk,ind,inb,ndj,locu)
      endif
      if(inb.gt.2) then
        jk=jk-1
        jj(2,j)=0
        goto 22
      endif
      do 119 ii=no(k0)+1,jk-1
         if(ja(jk).ne.ja(ii)) goto 119
         if(jb(jk).ne.jb(ii)) goto 119
         if(jc(jk).ne.jc(ii)) goto 119
         if(jm(jk).ne.jm(ii)) goto 119
         if(kk(j).eq.norb_inn-1) goto 607
           do i=1,8
            if(ind(i,jk).ne.ind(i,ii)) goto 119
          enddo
 607     jk=jk-1
         jj(2,j)=ii
         goto 22
 119  continue
      jj(2,j)=jk
      kk(jk)=kk(j)+1
      if(kk(jk).ne.kk(jk-1)) then
      no(k0)=jk-1
      end if
 22   continue
!      v=0
 2    if(ja(j).le.0) then
         goto 8
      else
         goto 5
      end if
 5    if(jc(j).eq.0) goto 3
      jk=jk+1
!                                     *************
!                                     *    d=2    *
!                                     *************
      ja(jk)=ja(j)-1
      jb(jk)=jb(j)+1
      jc(jk)=jc(j)-1
!-----------------------------------------------------------------------
      jm(jk)=mul_tab(lsm_inn(kk(j)+1),jm(j))
      do i=1,8
        ind(i,jk)=ind(i,j)
      enddo
      ind(lsm_inn(kk(j)+1),jk)=ind(lsm_inn(kk(j)+1),j)+1
      inb=0
      if(kk(j).eq.norb_inn-1) then
      call gugadrt_check_rcas3(jk,ind,inb,ndj,locu)
      endif
      if(inb.gt.2) then
        jk=jk-1
        jj(3,j)=0
        goto 44
      endif
      do 181 ii=no(k0)+1,jk-1
         if(ja(jk).ne.ja(ii)) goto 181
         if(jb(jk).ne.jb(ii)) goto 181
         if(jc(jk).ne.jc(ii)) goto 181
         if(jm(jk).ne.jm(ii)) goto 181
         if(kk(j).eq.norb_inn-1) goto 603
         do i=1,8
           if(ind(i,jk).ne.ind(i,ii)) goto 181
         enddo
 603     jk=jk-1
         jj(3,j)=ii
         goto 44
 181  continue
      jj(3,j)=jk
      kk(jk)=kk(j)+1
      if(kk(jk).ne.kk(jk-1)) then
      no(k0)=jk-1
      end if
 44   continue
!      v=0
!                                         ************
!                                         *    d=3   *
!                                         ************

 3    jk=jk+1
      ja(jk)=ja(j)-1
      jb(jk)=jb(j)
      jc(jk)=jc(j)
      jm(jk)=jm(j)
      do i=1,8
        ind(i,jk)=ind(i,j)
      enddo
      ind(lsm_inn(kk(j)+1),jk)=ind(lsm_inn(kk(j)+1),j)+2
         inb=0
      if(kk(j).eq.norb_inn-1) then
      call gugadrt_check_rcas3(jk,ind,inb,ndj,locu)
      endif
      if(inb.gt.2) then
        jk=jk-1
        jj(4,j)=0
        goto 88
      endif
      do 118 ii=no(k0)+1,jk-1
         if(ja(jk).ne.ja(ii)) goto 118
         if(jb(jk).ne.jb(ii)) goto 118
         if(jc(jk).ne.jc(ii)) goto 118
         if(jm(jk).ne.jm(ii)) goto 118
         if(kk(j).eq.norb_inn-1) goto 605
         do 606 i=1,8
           if(ind(i,jk).ne.ind(i,ii)) goto 118
 606     continue
 605     jk=jk-1
         jj(4,j)=ii
         go to 88
 118  continue
      jj(4,j)=jk
      kk(jk)=kk(j)+1
      if(kk(jk).ne.kk(jk-1)) then
      no(k0)=jk-1
      end if
 88   continue
!      v=0
      if(jk.gt.max_node) then
        write(6,*) '    the nomber of j exceeds max_node',max_node
        call abend
!       stop 777
      endif
 8    if(kk(j).le.norb_inn-1)goto 7
!  ************** external space  *************

      id=no(norb_inn)
      write(6,*)
      write(6,*)'    id=no(norb_inn)',id
      do 43 idd=no(norb_inn-1)+1,id
         if(ja(idd).eq.0.and.jb(idd).eq.0.and.jm(idd).eq.1) jv=idd
         if(ja(idd).eq.0.and.jb(idd).eq.1) then
           if(jm(idd).eq.1) jd(1)=idd
           if(jm(idd).eq.2) jd(2)=idd
           if(jm(idd).eq.3) jd(3)=idd
           if(jm(idd).eq.4) jd(4)=idd
           if(jm(idd).eq.5) jd(5)=idd
           if(jm(idd).eq.6) jd(6)=idd
           if(jm(idd).eq.7) jd(7)=idd
           if(jm(idd).eq.8) jd(8)=idd
         endif
         if(ja(idd).eq.0.and.jb(idd).eq.2) then
           if(jm(idd).eq.1) jt(1)=idd
           if(jm(idd).eq.2) jt(2)=idd
           if(jm(idd).eq.3) jt(3)=idd
           if(jm(idd).eq.4) jt(4)=idd
           if(jm(idd).eq.5) jt(5)=idd
           if(jm(idd).eq.6) jt(6)=idd
           if(jm(idd).eq.7) jt(7)=idd
           if(jm(idd).eq.8) jt(8)=idd
         endif
         if(ja(idd).eq.1.and.jb(idd).eq.0) then
           if(jm(idd).eq.1) js(1)=idd
           if(jm(idd).eq.2) js(2)=idd
           if(jm(idd).eq.3) js(3)=idd
           if(jm(idd).eq.4) js(4)=idd
           if(jm(idd).eq.5) js(5)=idd
           if(jm(idd).eq.6) js(6)=idd
           if(jm(idd).eq.7) js(7)=idd
           if(jm(idd).eq.8) js(8)=idd
      endif
 43   continue

      iwy(1,jv)=1
      do im=1,ng_sm
        if(jd(im).ne.0.and.iseg_downwei(1+im).ne.0) iwy(1,jd(im))=1
        if(jt(im).ne.0.and.iseg_downwei(9+im).ne.0) iwy(1,jt(im))=1
        if(js(im).ne.0.and.iseg_downwei(17+im).ne.0) iwy(1,js(im))=1
      enddo

      do 21 i=1,4
      iwy(i,0)=0
 21   continue
      do 20 l=norb_inn-1,norb_dz,-1
         jps=no(l-1)+1
         jpe=no(l)
      do 19 jde=jps,jpe
         j1=jj(1,jde)
         j2=jj(2,jde)
         j3=jj(3,jde)
         j4=jj(4,jde)
         iwy(1,jde)=iwy(1,j1)+iwy(1,j2)+iwy(1,j3)+iwy(1,j4)
         if(iwy(1,jde).ne.0) goto 304
           if(l.eq.norb_dz) then
             nc=0
           else
             nc=no(l-2)
           endif
           do jp0=nc+1,no(l-1)
             if(jj(1,jp0).eq.jde) jj(1,jp0)=0
             if(jj(2,jp0).eq.jde) jj(2,jp0)=0
             if(jj(3,jp0).eq.jde) jj(3,jp0)=0
             if(jj(4,jp0).eq.jde) jj(4,jp0)=0
           enddo
         goto 19
 304     if(j2.eq.0.or.iwy(1,j2).eq.0) goto 31
         iwy(2,jde)=iwy(1,j1)
 31      if(j3.eq.0.or.iwy(1,j3).eq.0) goto 32
         iwy(3,jde)=iwy(1,j1)+iwy(1,j2)
 32      if(j4.eq.0.or.iwy(1,j4).eq.0) goto 303
         iwy(4,jde)=iwy(1,jde)-iwy(1,j4)
 303     if(jde.eq.1) goto 19
         do 302 jp=jps,jde-1
           if(iwy(1,jp).ne.iwy(1,jde)) goto 302
           jq1=jj(1,jp)
           jq2=jj(2,jp)
           jq3=jj(3,jp)
           jq4=jj(4,jp)
         if(j1.ne.jq1.or.j2.ne.jq2.or.j3.ne.jq3.or.j4.ne.jq4) goto 302
           iwy(1,jde)=0
           do   jp0=no(l-2)+1,no(l-1)
             if(jj(1,jp0).eq.jde) jj(1,jp0)=jp
             if(jj(2,jp0).eq.jde) jj(2,jp0)=jp
             if(jj(3,jp0).eq.jde) jj(3,jp0)=jp
             if(jj(4,jp0).eq.jde) jj(4,jp0)=jp
           enddo
           goto 19
 302     continue
 19    continue
 20   continue
        it=mxnode
      itm(1)=1
        noh(norb_dz)=mxnode
      do 405 lr=norb_dz,norb_inn-1
        do 200 jp=no(lr)+1,no(lr+1)
          itm(jp)=0
          if(iwy(1,jp).ne.0) then
          it=it+1
          itm(jp)=it
          endif
 200    continue
      noh(lr+1)=it
 405  continue

      do 206 jpe=mxnode+1,id
        jp=itm(jpe)
      if(jp.eq.jpe) goto 206
      if(jp.eq.0) goto 206
        l=kk(jpe)
        jds=noh(l-2)+1
        jde=noh(l-1)
        do jp0=jds,jde
          if(jj(1,jp0).eq.jpe) jj(1,jp0)=jp
          if(jj(2,jp0).eq.jpe) jj(2,jp0)=jp
          if(jj(3,jp0).eq.jpe) jj(3,jp0)=jp
          if(jj(4,jp0).eq.jpe) jj(4,jp0)=jp
        enddo
      ja(jp)=ja(jpe)
      jb(jp)=jb(jpe)
      jc(jp)=jc(jpe)
      jm(jp)=jm(jpe)
      kk(jp)=kk(jpe)
      do 704 i=1,4
         iwy(i,jp)=iwy(i,jpe)
         jj(i,jp)=jj(i,jpe)
         ind(i,jp)=ind(i,jpe)
 704  continue
      do i=5,8
      ind(i,jp)=ind(i,jpe)
      enddo
 206  continue

!      open(10,file='rcas.out')
      no(norb_dz)=0
      no(norb_dz+1)=mxnode

!      write(nf10,'(2x,2i10)')norb_dz+1,no(norb_dz+1)
      do 706 lr=norb_dz,norb_inn
        no(lr+1)=noh(lr)
        write(6,'(2x,2i10)')lr+1,no(lr+1)
 706  continue

      itm(0)=0
      jv=itm(jv)
      do im=1,8
          jd(im)=itm(jd(im))
          jt(im)=itm(jt(im))
          js(im)=itm(js(im))
      enddo
      id=it
      if(it.ne.no(norb_inn+1))then
       write(6,*)'   rcas id is wrong!!   no(norb_inn)=',no(norb_inn),it
       call abend
!      stop
      end if
      iysum=0
      do j=1,mxnode
        iysum=iysum+iwy(1,j)
      enddo
!        write(6,*)'    end of rcas , node=',id,'  dimension=',iysum
      write(6,*)
      indd=no(norb_inn)
!      iprint=1
      if(iprint.eq.1) then
        write(6,*) "guga drt"
        write(6,506)
      endif
 506  format('       j    k   a  b  t jm    j0   j1   j2   j3       y1',&
     &       '    y2      y3         x   ind')

      do 541 j=1,id
        kk(j)=kk(j)+1
        if(iprint.eq.1) then
           write(6,507)j,kk(j),ja(j),jb(j),jm(j),                       &
     &             jj(1,j),jj(2,j),jj(3,j),jj(4,j),                     &
     &             iwy(2,j),iwy(3,j),iwy(4,j),iwy(1,j),(ind(i,j),i=1,8)
        endif
 541  continue
      write(6,*) 'end of rcas, drt ..........'
      write(6,*)

!      open(21,file="fort.drt",form="unformatted")
!      write(21) id
!      write(21) ja(1:id),jb(1:id),jm(1:id)
!      write(21) jj(1:4,0:id)
!      write(21) kk(0:id)
!      write(21) no(0:norb_inn+1)
!      write(21) jv,jd(1:8),jt(1:8),js(1:8)
!      close(21)

      call writedrt(id)
      call mma_deallocate(ind)
      call mma_deallocate(iwy)
 507  format(3x,2i5,1x,3i3,1x,4i5,1x,4i10,1x,8i2)
      end

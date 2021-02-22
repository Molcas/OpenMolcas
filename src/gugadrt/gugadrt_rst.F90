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
      subroutine gugadrt_rst(id,nndd)
!*******************************************
! 10 may 2007 - revised by wyb
!
#include "gendrt.fh"
#include "Sysdrt.fh"
!#ifndef _I8_
!      parameter (iintbit=32,n32int=4,n16int=2)
!#else
!      parameter (iintbit=64,n32int=2,n16int=1)
!#endif
#include "casrst_drt.fh"
#include "ref.fh"
      integer, pointer :: jabkm(:,:)
      integer, pointer :: ind(:,:)
      integer, pointer :: idjj(:,:)
      integer, pointer :: iwy(:,:)
      integer, pointer :: itm(:)
      dimension noh(max_innorb)
      dimension iextjj(n32int),iextii(n32int),jjabkm(1:n16int),         &
     &          jkabkm(1:n16int),iiabkm(1:n16int)
! estimate memory
      if(n_ref.gt.20) then
        if(norb_act.ge.10) then
          if(iintbit.eq.32) mxtnode=1000000
          if(iintbit.eq.64) mxtnode=10000000
        else
          mxtnode=1000000
        endif
      else
        if(norb_act.gt.10) then
          mxtnode=1000000
        else
          mxtnode=500000
        endif
      endif
      write(6,*)
      write(6,*) ' now generate distinct row tableau'
      allocate(jabkm(n16int,mxtnode))
      allocate(ind(n32int,0:mxtnode))
      allocate(idjj(4,0:mxtnode))
      allocate(iwy(4,0:mxtnode))
      allocate(itm(0:mxtnode))

      nm  =ns_sm
      ndj =n_ref
      no(1:norb_dz+1)=0
      jabkm(1:n16int,1:mxtnode)=0
      ind(1:n32int,0:mxtnode)=0
      idjj(1:4,0:mxtnode)=0
      iwy=0
      itm=0
      nrefbit=n32int
      jj(1:4,0:max_node)=0

! 8 bits
      iextbit=2
      nextbit=iintbit/iextbit
! 16 bits
      iabcbit=16
      nabcbit=iintbit/iabcbit

! v node
      j=0
      ja0=ja_sys
      jb0=jb_sys
      jc0=jc_sys
      write(6,"(4(1x,i4))") ja0,jb0,jc0
!  v_node
      ja(1)=ja0
      jb(1)=jb0
      jm(1)=nm
      kk(1)=norb_dz
      jaj=ja(1)
      jbj=jb(1)
      jmj=jm(1)
      kkj=kk(1)

      jjabkm(1:n16int)=0
      call wrtabkm(jjabkm,n16int,nabcbit,iabcbit,jaj,jbj,jmj,kkj)
      jabkm(1:n16int,1)=jjabkm(1:n16int)
      ind(1:nrefbit,1)=0

!  d_node
      imd=0
      do node=2,9
        imd=imd+1
        if(nu_ad(node).eq.0) cycle
        ja(node)=ja(1)
        jb(node)=jb0+1
        jm(node)=imd
        kk(node)=norb_dz
        jaj=ja(node)
        jbj=jb(node)
        jmj=jm(node)
        kkj=kk(node)

        jjabkm(1:n16int)=0
        call wrtabkm(jjabkm,n16int,nabcbit,iabcbit,jaj,jbj,jmj,kkj)
        jabkm(1:n16int,node)=jjabkm(1:n16int)
        ind(1:nrefbit,node)=0
      enddo

!  t_node
      imt=0
      do node=10,17
        imt=imt+1
        if(nu_ad(node).eq.0) cycle
        ja(node)=ja(1)
        jb(node)=jb0+2
        jm(node)=imt
        kk(node)=norb_dz
        jaj=ja(node)
        jbj=jb(node)
        jmj=jm(node)
        kkj=kk(node)

        jjabkm(1:n16int)=0
        call wrtabkm(jjabkm,n16int,nabcbit,iabcbit,jaj,jbj,jmj,kkj)
        jabkm(1:n16int,node)=jjabkm(1:n16int)
        ind(1:nrefbit,node)=0
      enddo
!  s_node
      ims=0
      do node=18,25
        ims=ims+1
        if(nu_ad(node).eq.0) cycle
        ja(node)=ja(1)+1
        jb(node)=jb0
        jm(node)=ims
        kk(node)=norb_dz
        jaj=ja(node)
        jbj=jb(node)
        jmj=jm(node)
        kkj=kk(node)
        jjabkm(1:n16int)=0
        call wrtabkm(jjabkm,n16int,nabcbit,iabcbit,jaj,jbj,jmj,kkj)
        jabkm(1:n16int,node)=jjabkm(1:n16int)
        ind(1:nrefbit,node)=0
      enddo
      if(jroute_sys.gt.1) then
!  d'_node
        imd=0
        do node=26,33
          imd=imd+1
          if(nu_ad(node).eq.0) cycle
          ja(node)=ja(1)+1
          jb(node)=jb0-1
          jm(node)=imd
          kk(node)=norb_dz
          jaj=ja(node)
          jbj=jb(node)
          jmj=jm(node)
          kkj=kk(node)
          jjabkm(1:n16int)=0
          call wrtabkm(jjabkm,n16int,nabcbit,iabcbit,jaj,jbj,jmj,kkj)
          jabkm(1:n16int,node)=jjabkm(1:n16int)
          ind(1:nrefbit,node)=0
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
          jm(node)=imt
          kk(node)=norb_dz
          jaj=ja(node)
          jbj=jb(node)
          jmj=jm(node)
          kkj=kk(node)
          jjabkm(1:n16int)=0
          call wrtabkm(jjabkm,n16int,nabcbit,iabcbit,jaj,jbj,jmj,kkj)
          jabkm(1:n16int,node)=jjabkm(1:n16int)
          ind(1:nrefbit,node)=0
        enddo
      endif
      no(norb_dz)=mxnode

      j=0
      jk=mxnode
      kk(jk)=norb_dz
!      call packnod(idkk,jk,norb_dz,nabcbit,iabcbit,mxtnode)

!**********************************************************************
!
 7    j=j+1
      if(j.le.mxnode) then
        if(nu_ad(j).eq.0) goto 7
      endif
      if(j.lt.mxnode) then
        jaj=ja(j)
        jbj=jb(j)
        jmj=jm(j)
        kkj=kk(j)
        k0=kkj
      else
        jjabkm(1:n16int)=jabkm(1:n16int,j)
        call redabkm(jjabkm,n16int,nabcbit,iabcbit,jaj,jbj,jmj,kkj)
        k0=kkj
      endif
      jkabkm(1:n16int)=jabkm(1:n16int,jk)
      call redabkm(jkabkm,n16int,nabcbit,iabcbit,jajk,jbjk,jmjk,kkjk)
      if(kkjk.ne.k0+1) no(k0)=jk
      jk=jk+1

      if(jk.gt.mxtnode) then
        write(6,*) ' the number of j exceeds max_node',mxtnode
        call abend
!       stop 777
!        call errexit(777)
      endif
! =========================================================
!                                            ***********
!                                            *   d=0   *
!                                            ***********

      jatmp=jaj
      jbtmp=jbj
      jmtmp=jmj
      kktmp=kkj+1


      iextjj(1:nrefbit)=ind(1:nrefbit,j)
      call gugadrt_njexcit(iextjj,nrefbit,iextbit,nextbit,ivalid,       &
     &                      0,kttmp,k0)
      if(ivalid.eq.0) then
        jk=jk-1
        goto 11
      endif

      if(k0.eq.norb_inn-1) then
!  act||ext
        if((jatmp.eq.1.or.jbtmp.eq.2).and.kttmp.ne.0) then
          jk=jk-1
          goto 11
        endif
        if(jbtmp.eq.1.and.kttmp.eq.2) then
          jk=jk-1
          goto 11
        endif
      end if

      call wrtabkm(jkabkm,n16int,nabcbit,iabcbit,jatmp,jbtmp,           &
     &             jmtmp,kktmp)
      do 18 ii=no(k0)+1,jk-1
        iiabkm(1:n16int)=jabkm(1:n16int,ii)
        do i=1,n16int
          if(iiabkm(i).ne.jkabkm(i)) goto 18
        enddo
        if(k0.eq.norb_inn-1) goto 601
        iextii(1:nrefbit)=ind(1:nrefbit,ii)
        do i=1,nrefbit
          if(iextii(i).ne.iextjj(i)) goto 18
        enddo
 601    jk=jk-1
        idjj(1,j)=ii
        goto 11
 18   continue

      ind(1:nrefbit,jk)=iextjj(1:nrefbit)
      idjj(1,j)=jk
      jabkm(1:n16int,jk)=jkabkm(1:n16int)

      iiabkm(1:n16int)=jabkm(1:n16int,jk-1)
      call upacknod(iiabkm,4,kj1,nabcbit,iabcbit,n16int)
      if(kktmp.ne.kj1) no(k0)=jk-1

 11   continue
! =========================================================
      if(jbj.eq.0) goto 2
      jk=jk+1
!                                            ***********
!                                            *   d=1   *
!                                            ***********
      jatmp=jaj
      jbtmp=jbj-1
      jmtmp=mul_tab(lsm_inn(k0+1),jmj)
      kktmp=kkj+1

      iextjj(1:nrefbit)=ind(1:nrefbit,j)
      call gugadrt_njexcit(iextjj,nrefbit,iextbit,nextbit,ivalid,       &
     &                     1,kttmp,k0)

      if(ivalid.eq.0) then
        jk=jk-1
        goto 22
      endif

      if(k0.eq.norb_inn-1) then
!  act||ext
        if((jatmp.eq.1.or.jbtmp.eq.2).and.kttmp.ne.0) then
          jk=jk-1
          goto 22
        endif
        if(jbtmp.eq.1.and.kttmp.eq.2) then
          jk=jk-1
          goto 22
        endif
      end if

      call wrtabkm(jkabkm,n16int,nabcbit,iabcbit,jatmp,jbtmp,           &
     &             jmtmp,kktmp)
      do 28 ii=no(k0)+1,jk-1
        iiabkm(1:n16int)=jabkm(1:n16int,ii)
        do i=1,n16int
          if(iiabkm(i).ne.jkabkm(i)) goto 28
        enddo
        if(k0.eq.norb_inn-1) goto 602
        iextii(1:nrefbit)=ind(1:nrefbit,ii)
        do i=1,nrefbit
          if(iextii(i).ne.iextjj(i)) goto 28
        enddo
 602    jk=jk-1
        idjj(2,j)=ii
        goto 22
 28   continue

      ind(1:nrefbit,jk)=iextjj(1:nrefbit)
      idjj(2,j)=jk
      jabkm(1:n16int,jk)=jkabkm(1:n16int)

      iiabkm(1:n16int)=jabkm(1:n16int,jk-1)
      call upacknod(iiabkm,4,kj1,nabcbit,iabcbit,n16int)
      if(kktmp.ne.kj1) no(k0)=jk-1

 22   continue
! =========================================================
 2    jac=norb_all-jaj-jbj
      if(jaj.le.0) then
         goto 8
      else
         goto 5
      end if
 5    if(jac.eq.0) goto 3
      jk=jk+1
!                                     *************
!                                     *    d=2    *
!                                     *************
      jatmp=jaj-1
      jbtmp=jbj+1
      jmtmp=mul_tab(lsm_inn(k0+1),jmj)
      kktmp=kkj+1

      iextjj(1:nrefbit)=ind(1:nrefbit,j)
      call gugadrt_njexcit(iextjj,nrefbit,iextbit,nextbit,ivalid,       &
     &                     2,kttmp,k0)
      if(ivalid.eq.0) then
        jk=jk-1
        goto 33
      endif

      if(k0.eq.norb_inn-1) then
!  act||ext
        if((jatmp.eq.1.or.jbtmp.eq.2).and.kttmp.ne.0) then
          jk=jk-1
          goto 33
        endif
        if(jbtmp.eq.1.and.kttmp.eq.2) then
          jk=jk-1
          goto 33
        endif
      end if

      call wrtabkm(jkabkm,n16int,nabcbit,iabcbit,jatmp,jbtmp,           &
     &             jmtmp,kktmp)
      do 38 ii=no(k0)+1,jk-1
        iiabkm(1:n16int)=jabkm(1:n16int,ii)
        do i=1,n16int
          if(iiabkm(i).ne.jkabkm(i)) goto 38
        enddo
        if(k0.eq.norb_inn-1) goto 603
        iextii(1:nrefbit)=ind(1:nrefbit,ii)
        do i=1,nrefbit
          if(iextii(i).ne.iextjj(i)) goto 38
        enddo
 603    jk=jk-1
        idjj(3,j)=ii
        goto 33
 38   continue

      ind(1:nrefbit,jk)=iextjj(1:nrefbit)
      idjj(3,j)=jk
      jabkm(1:n16int,jk)=jkabkm(1:n16int)

      iiabkm(1:n16int)=jabkm(1:n16int,jk-1)
      call upacknod(iiabkm,4,kj1,nabcbit,iabcbit,n16int)
      if(kktmp.ne.kj1) no(k0)=jk-1

 33   continue
! =========================================================
!                                         ************
!                                         *    d=3   *
!                                         ************

 3    jk=jk+1
      jatmp=jaj-1
      jbtmp=jbj
      jmtmp=jmj
      kktmp=kkj+1

      iextjj(1:nrefbit)=ind(1:nrefbit,j)
      call gugadrt_njexcit(iextjj,nrefbit,iextbit,nextbit,ivalid,       &
     &                     3,kttmp,k0)
      if(ivalid.eq.0) then
        jk=jk-1
        goto 44
      endif

      if(k0.eq.norb_inn-1) then
!  act||ext
        if((jatmp.eq.1.or.jbtmp.eq.2).and.kttmp.ne.0) then
          jk=jk-1
          goto 44
        endif
        if(jbtmp.eq.1.and.kttmp.eq.2) then
          jk=jk-1
          goto 44
        endif
      end if

      call wrtabkm(jkabkm,n16int,nabcbit,iabcbit,jatmp,jbtmp,           &
     &             jmtmp,kktmp)
      do 48 ii=no(k0)+1,jk-1
        iiabkm(1:n16int)=jabkm(1:n16int,ii)
        do i=1,n16int
          if(iiabkm(i).ne.jkabkm(i)) goto 48
        enddo
        if(k0.eq.norb_inn-1) goto 604
        iextii(1:nrefbit)=ind(1:nrefbit,ii)
        do i=1,nrefbit
          if(iextii(i).ne.iextjj(i)) goto 48
        enddo
 604    jk=jk-1
        idjj(4,j)=ii
        goto 44
 48   continue

      ind(1:nrefbit,jk)=iextjj(1:nrefbit)
      idjj(4,j)=jk
      jabkm(1:n16int,jk)=jkabkm(1:n16int)

      iiabkm(1:n16int)=jabkm(1:n16int,jk-1)
      call upacknod(iiabkm,4,kj1,nabcbit,iabcbit,n16int)
      if(kktmp.ne.kj1) no(k0)=jk-1

 44   continue
!      if(k0.eq.7) goto 999
!===========================================================
 8    if(k0.le.norb_inn-1) goto 7

!999   continue
!      if(iprint.eq.1) then
!        open(100,file="tmp.dat")
!        do i=1,jk
!          jkabkm(1:n16int)=jabkm(1:n16int,i)
!          call redabkm(jkabkm,n16int,nabcbit,iabcbit,jaj,jbj,jmj,kkj)
!          j1=idjj(1,i)
!          j2=idjj(2,i)
!          j3=idjj(3,i)
!          j4=idjj(4,i)
!          write(100,"(5(1x,i8))") i,j1,j2,j3,j4
!        enddo
!        close(100)
!      endif
!      write(6,"(10(1x,i5))") no(1:norb_all)
      write(6,*)
!  ************** external space  *************
      id=no(norb_inn)
!      write(6,508) 'befor,no=',(no(i),i=norb_dz,norb_inn)
      do idd=no(norb_inn-1)+1,id
        jjabkm(1:n16int)=jabkm(1:n16int,idd)
        call redabkm(jjabkm,n16int,nabcbit,iabcbit,jaj,jbj,jmj,kkj)
        if(jaj.eq.0.and.jbj.eq.0.and.jmj.eq.1) jv=idd
        if(jaj.eq.0.and.jbj.eq.1) then
          do i=1,ng_sm
            if(jmj.eq.i) then
              jd(i)=idd
            endif
          enddo
        endif
        if(jaj.eq.0.and.jbj.eq.2) then
          do i=1,ng_sm
            if(jmj.eq.i) then
              jt(i)=idd
            endif
          enddo
        endif
        if(jaj.eq.1.and.jbj.eq.0) then
          do i=1,ng_sm
            if(jmj.eq.i) then
              js(i)=idd
            endif
          enddo
        endif
      enddo

      iwy(1,jv)=1
      do im=1,ng_sm
        if(jd(im).ne.0.and.iseg_downwei(1+im).ne.0) then
          iwy(1,jd(im))=1
        endif
        if(jt(im).ne.0.and.iseg_downwei(9+im).ne.0) then
          iwy(1,jt(im))=1
        endif
        if(js(im).ne.0.and.iseg_downwei(17+im).ne.0) then
          iwy(1,js(im))=1
        endif
      enddo

      iwy(1:4,0)=0
      do 20 l=norb_inn-1,norb_dz,-1
        jps=no(l-1)+1
        jpe=no(l)
        do 21 jde=jps,jpe
          j1=idjj(1,jde)
          j2=idjj(2,jde)
          j3=idjj(3,jde)
          j4=idjj(4,jde)
          iwy(1,jde)=iwy(1,j1)+iwy(1,j2)+iwy(1,j3)+iwy(1,j4)
          if(iwy(1,jde).ne.0) goto 304
          do jp0=no(l-2)+1,no(l-1)
            if(idjj(1,jp0).eq.jde) idjj(1,jp0)=0
            if(idjj(2,jp0).eq.jde) idjj(2,jp0)=0
            if(idjj(3,jp0).eq.jde) idjj(3,jp0)=0
            if(idjj(4,jp0).eq.jde) idjj(4,jp0)=0
          enddo
          goto 21
 304      if(j2.eq.0.or.iwy(1,j2).eq.0) goto 31
          iwy(2,jde)=iwy(1,j1)
 31       if(j3.eq.0.or.iwy(1,j3).eq.0) goto 32
          iwy(3,jde)=iwy(1,j1)+iwy(1,j2)
 32       if(j4.eq.0.or.iwy(1,j4).eq.0) goto 303
          iwy(4,jde)=iwy(1,jde)-iwy(1,j4)
 303      if(jde.eq.1) goto 21
          do 302 jp=jps,jde-1
            if(iwy(1,jp).ne.iwy(1,jde)) goto 302
            jq1=idjj(1,jp)
            jq2=idjj(2,jp)
            jq3=idjj(3,jp)
            jq4=idjj(4,jp)
            if(j1.ne.jq1.or.j2.ne.jq2.or.j3.ne.jq3.or.j4.ne.jq4)        &
     &        goto 302
            iwy(1,jde)=0
            do jp0=no(l-2)+1,no(l-1)
              if(idjj(1,jp0).eq.jde) idjj(1,jp0)=jp
              if(idjj(2,jp0).eq.jde) idjj(2,jp0)=jp
              if(idjj(3,jp0).eq.jde) idjj(3,jp0)=jp
              if(idjj(4,jp0).eq.jde) idjj(4,jp0)=jp
            enddo
            goto 21
 302      continue
 21     continue
 20   continue

      it=mxnode
      itm(1)=1
      noh=0
      noh(norb_dz)=mxnode
      do lr=norb_dz,norb_inn-1
        do jp=no(lr)+1,no(lr+1)
          itm(jp)=0
          if(iwy(1,jp).ne.0) then
            it=it+1
            itm(jp)=it
          endif
        enddo
        noh(lr+1)=it
      enddo

      do jpe=1,mxnode
        if(nu_ad(jpe).eq.0) cycle
        jj(1,jpe)=idjj(1,jpe)
        jj(2,jpe)=idjj(2,jpe)
        jj(3,jpe)=idjj(3,jpe)
        jj(4,jpe)=idjj(4,jpe)
      enddo

      do jpe=mxnode+1,id
        jp=itm(jpe)
        if(jp.eq.0) cycle
        j1=idjj(1,jpe)
        j2=idjj(2,jpe)
        j3=idjj(3,jpe)
        j4=idjj(4,jpe)
        jjabkm(1:n16int)=jabkm(1:n16int,jpe)
        call redabkm(jjabkm,n16int,nabcbit,iabcbit,jaj,jbj,jmj,kkj)
        l=kkj
        jds=noh(l-2)+1
        jde=noh(l-1)
        do jp0=jds,jde
          if(jj(1,jp0).eq.jpe) jj(1,jp0)=jp
          if(jj(2,jp0).eq.jpe) jj(2,jp0)=jp
          if(jj(3,jp0).eq.jpe) jj(3,jp0)=jp
          if(jj(4,jp0).eq.jpe) jj(4,jp0)=jp
        enddo
        ja(jp)=jaj
        jb(jp)=jbj
        jm(jp)=jmj
        kk(jp)=kkj
        do i=1,4
          iwy(i,jp)=iwy(i,jpe)
          ji=idjj(i,jpe)
          if(itm(ji).ne.0) then
            jj(i,jp)=ji
          endif
        enddo
      enddo

!      open(10,file='rst.out')
      no(norb_dz)=0
      no(norb_dz+1)=mxnode
!      write(6,*) '   end of rst, drt ..........'
!      write(6,'(2x,2i10)')norb_dz+1,no(norb_dz+1)
      do 706 lr=norb_dz+1,norb_inn
        no(lr+1)=noh(lr)
 706  continue
      jv=itm(jv)
      itm(0)=0
      do im=1,8
        jd(im)=itm(jd(im))
        jt(im)=itm(jt(im))
        js(im)=itm(js(im))
      enddo
      iysum=0
      do j=1,mxnode
        iysum=iysum+iwy(1,j)
      enddo
      id=it
      if(it.ne.no(norb_inn+1))then
        write(6,*)'   rst id is wrong!!   no(norb_inn)=',no(norb_inn),it
        call abend
!       stop
      end if

      write(6,*)
!      iprint=1
      if(iprint.eq.1) then
        write(6,*) "guga drt"
        write(6,506)
      endif
 506  format('       j    k   a  b  t jm    j0   j1   j2   j3       y1',&
     &       '       y2      y3         x   ind')
      nndd=no(norb_inn)
      do 541 j=1,id
          kk(j)=kk(j)+1
        if(iprint.eq.1) then
          write(6,510)j,kk(j),ja(j),jb(j),jm(j),                        &
     &             jj(1,j),jj(2,j),jj(3,j),jj(4,j)
        endif
 541  continue
      write(6,*)
      write(6,*) 'end of rst, drt ..........'

      deallocate(jabkm)
      deallocate(ind)
      deallocate(idjj)
      deallocate(iwy)
      deallocate(itm)

!
!      open(21,file="fort.drt",form="unformatted")
!      write(21) id
!      write(21) ja(1:id),jb(1:id),jm(1:id)
!      write(21) jj(1:4,0:id)
!      write(21) kk(0:id)
!      write(21) no(0:norb_inn+1)
!      write(21) jv,jd(1:8),jt(1:8),js(1:8)
!      close(21)

      call writedrt(id)
      return
!508   format(3x,a10,1x,i5,1x,16i8)
 510  format(3x,10(1x,i7))
      end

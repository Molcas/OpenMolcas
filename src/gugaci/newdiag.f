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
* Copyright (C) 2008, Bingbing Suo                                     *
************************************************************************
C 04 Jul 2008 -BSUO- generalized davidson diagnolization routines
      subroutine gdavdiag(istate,mtidx,veistate)
#include "drt_h.fh"
      parameter (maxdimlu=1000,maxdimgit=10000)

c compute (h00-rouk)-1, and save the low triagular part of the
c inverse matrix at value_tmp(1:lent1-1)
      call cimtinverse(ndimt,veistate)

c h00 (h00-rouk)-1 * vector1(1:ndim)

      lent0=1
      lent1=maxdimlu*maxdimlu
      lent2=lent1+ndimt
      value_lpext(lent1+1:lent2)=vector1(mtidx+1:mtidx+ndimt)
c subroutine dsymv in blas is used
c      call dsymv('l',ndimt,1.0,ndimt,value_lpext(1),1,
c     *           vector1(mtidx+1),1)
      call matmultv(value_lpext,ndimt,maxdimlu,
     *              value_lpext(lent1+1:lent1+ndimt),
     *              vector1(mtidx+1:ndimt))
c   from h00+1 to nci_dim
      do l=ndimt+1,nci_dim
         depff=vector2(l)-veistate
c               if(abs(depff).le.depc) depff=depc
         vector1(mtidx+l)=vector1(mtidx+l)/depff
      enddo

c Avoid unused argument warnings
      if (.false.) call Unused_integer(istate)
      end

      subroutine matmultv(a,n,np,x,y)
c matrix a(n,n), vector x(n),y(n)
c ax=y
      implicit real*8 (a-h,o-z)
      dimension a(np,np),x(n),y(n)

      y(1:n)=0.d0
      do i=1,n
        do j=1,n
          y(i)=y(i)+a(i,j)*x(j)
        enddo
      enddo

      return
c...end of matmulv
      end

      subroutine cimtinverse(ndimt,vrouk)
#include "drt_h.fh"
      common /scratch/ tmpdir,len_str
      character*256 tmpdir
      character*256 filename

      parameter (maxdimlu=1000,maxdimgit=10000)

c      is=1
c      if(is.eq.1) then
c v'-v space is chosen as h0 part in gdvdd
c        ndimt=iw_sta(2,1)
c      else
c v'+d'+t'+s'-v space
c        ndimt=iw_sta(1,2)
c      endif

c      if(ndimt.le.maxdimlu) then
        filename=tmpdir(1:len_str)//"/fort.23"
        lent=len_str+8
!        write(6,*) "new diag"
!        stop 1999
!        open(nf23,file=filename(1:lent),form="unformatted")
!        read(nf23) ndimt,lent,value_lpext(1:lent)
!        close(nf23)
        neh0=maxdimlu*maxdimlu
        num=2*neh0+lent
        if(num.gt.max_tmpvalue) then
          write(6,*) "no enough scracht vectors, error 1000"
     *               ,num,max_tmpvalue,neh0,lent,ndimt
#ifdef MOLPRO
#else
      call qtrace
      call abend()
#endif
#ifdef _XIANEST_
      call qexit()
#endif
!         call abend
!          stop 1000
        endif
        value_lpext(2*neh0+1:2*neh0+lent)=value_lpext(1:lent)
        call matinverse(value_lpext,
     *                  value_lpext(neh0+1:neh0+maxdimlu**2),ndimt,
     *                  maxdimlu,value_lpext(2*neh0+1:2*neh0+lent),lent,
     *                  vrouk)
c      endif

c      stop 1000
c the inverse matrix a is also a symmetric matrix,
c so we only need to keep the low triagular part of a.
c      call savelowtra(value_lpext(1),value_lpext(neh0+1),ndimt,
c     *                maxdimlu,neh0)
      value_lpext(1:neh0)=value_lpext(neh0+1:2*neh0)
      value_lpext(neh0+1:2*neh0)=0.d0
      return
      end

      subroutine matinverse(a,y,n,np,valvec,lent,vrouk)
      implicit real*8 (a-h,o-z)
      parameter (maxdimlu=1000,maxdimgit=10000)
      integer np, indx(n)
      dimension a(np,np),y(np,np),valvec(lent)
      parameter (zero=0.d0)
*i64  parameter (zero=0.e0)

      k=0
      do i=1,n
       do j=1,i-1
          k=k+1
          a(i,j)=valvec(k)
          a(j,i)=valvec(k)
        enddo
        k=k+1
        a(i,i)=valvec(k)-vrouk
      enddo
c      write(6,*) "ndim0 v-v,",np
c      do i=1,10
c        write(6,"(10(f12.6,1x))") a(i,1:10)
c      enddo
c      stop 888

      y(1:np,1:np)=zero
      do i=1,n
        y(i,i)=1.0
      enddo
      call ludcmp(a,n,np,indx,d)
      do i=1,n
        call lubksb(a,n,np,indx,y(1,i))
      enddo

      return
      end

      subroutine savelowtra(varry,a,ndimt,maxdimlu,neh0)
      implicit real*8 (a-h,o-z)
      dimension varry(neh0),a(maxdimlu,maxdimlu)
      parameter (dzero=0.d0)
*i64  parameter (dzero=0.e0)

      varry(1:neh0)=dzero
      k=0
      do i=1,ndimt
        do j=1,i-1
          k=k+1
          varry(k)=a(i,j)
        enddo
        k=k+1
        varry(k)=a(i,i)
      enddo

      return
c...end of savelowtra
      end

      subroutine comph00()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      parameter (maxdimlu=1000)
      common /scratch/ tmpdir,len_str
      character*256 tmpdir
c     character*256 filename

      vector2(1:max_vector)=0.d0
      log_prod=2
      call readint(nf21,1,vint_ci,tmpdir,len_str)

      isend=0
      do i=1,mxnode
        ncid=iw_sta(i,1)
        if(ncid.lt.maxdimlu) then
          isend=i-1
        else
          goto 10
        endif
      enddo
10    continue
      if(iw_sta(1,2).lt.maxdimlu) isend=mxnode
      write(6,*) "isend=",isend
c      stop 888

      if(isend.lt.1) then
        write(6,*) " could not use generalized davdison diagnolization",
     *             " mehtod"
        logic_tdav=.true.
        return
      endif

c      write(13,*)' vv_drt_ci_new'
c=========== vv_one_subdrt =============================
      iml=1
      imr=1
      ipae=1
      jpae=nu_ae(ipae)
      ipael=ipae
      ipaer=ipae
      jpad_sta=1
      jpad_end=isend
      do jpad=jpad_sta,jpad_end
        call get_jpadty(jpad,jpadty,jml)
        ndime=iw_downwei(jpad,ipae)
        if(ndime.eq.0) cycle
        jpadr=jpad
        jpadl=jpad
        call seg_drt()
        call copy_to_drtl()
        jpadlr=map_jplr(jpadty,jpadty)
        jmr=jml
c        write(6,*) "jpad=",jpad,"ipae=",ipae

        call dbl_space_loop()
        call cloop_in_act()
        call ploop_in_act()
      enddo

c      stop 1000

c=========== vv_two_subdrts with same tail different head ==============
      iml=1
      imr=1
      ipae=1
      jpae=nu_ae(ipae)
      ipael=1
      jpael=nu_ae(ipae)
      ipaer=1
      jpaer=nu_ae(ipae)
      int_dd_drl=int_dd_offset(iml,iml)
      call logicg_dd(iml,iml)

      jpad1=1
      call get_jpadty(jpad1,jpad1ty,jm1)

      do jpad2=jpad1+1,jpad_end
        call get_jpadty(jpad2,jpad2ty,jm2)
        ndime=iw_downwei(jpad2,ipae)
        if(ndime.eq.0) cycle

        iwup=jpad_upwei(jpad2)
        isegdown=iseg_downwei(ipae)
        ndime=ndime*iwup*isegdown
c        vector2(1:ndime+ndim_h0)=0.0d0
!        write(6,'(2x,3i8)') jpad2,ipae,ndime
        jpad=jpad1
        call seg_drt()
        call copy_to_drtl()
        jpad=jpad2
        call seg_drt()

        jpadl=jpad1
        jml  =jm1
        jpadr=jpad2
        jmr  =jm2
        jpadlr=map_jplr(jpad1ty,jpad2ty)

c=========== do sta h_off_diagnal elements for subdrt_lr  ==============
        call cloop_in_act()
        call ploop_in_act()
c=========== do end          =============================

        call copy_to_drtl()
        jpad=jpad1
        call seg_drt()
        jpadl=jpad2
        jml  =jm2
        jpadr=jpad1
        jmr  =jm1
        jpadlr=map_jplr(jpad2ty,jpad1ty)
c========== do sta h_off_diagnal elements for subdrt_rl  ==============
        call cloop_in_act()
        call ploop_in_act()
c=========== do end          =============================
      enddo
c=========== vv_end two_subdrts with same tail =========================

c save the diagnol element into h00 matrix
      ndimh00=iw_sta(isend+1,1)
      if(iw_sta(1,2).lt.maxdimlu) ndimh00=iw_sta(1,2)
      do i=1,ndimh00
        j=i*(i+1)/2
        vector2(j)=vector1(i)
      enddo
!      ndimt=ndimh00*(ndimh00+1)/2
!      filename=tmpdir(1:len_str)//"/fort.23"
!      lent=len_str+8
!      open(nf23,file=filename(1:lent),form="unformatted")
!      write(6,*) "ndimt=",ndimh00,ndimt
      do i=1,mxnode
        write(6,*) "iw_sta(i,1)",i,iw_sta(i,1)
      enddo
      write(6,*) "iw_sta(1,2)",iw_sta(1,2),isend

!      write(nf23) ndimh00,ndimt,vector2(1:ndimt)
!      close(nf23)
c      stop 1000

      return
c...end of comph00
      end

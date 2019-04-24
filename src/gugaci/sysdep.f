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
c 11 feb 2007 -bsuo- added by suo bing
c                    this file contains subroutines which depend on plat
      subroutine version_info
#include "drt_h.fh"

      write(6,'(10x,a42)')'*****************************************'
      write(6,'(10x,a42)')'*      Xian-ci mrci program             *'
      write(6,'(10x,a42)')'*     Institute of Modern Physics       *'
      write(6,'(10x,a42)')'*        Northwest University           *'
      write(6,'(10x,a42)')'*        xian, shaanxi, china           *'
      write(6,'(10x,a42)')'*                                       *'
      write(6,'(10x,a42)')'*        report bugs and errors         *'
      write(6,'(10x,a42)')'*           wzy@nwu.edu.cn              *'
      write(6,'(10x,a42)')'*        yubin_wang@hotmail.com         *'
      write(6,'(10x,a42)')'*       bingbing_suo@hotmail.com        *'
      write(6,'(10x,a42)')'*                                       *'
      write(6,'(10x,a42)')'*****************************************'
      write(6,*)
      write(6,*)

      call get_date()

      return
      end

      subroutine get_date()
#ifdef _XIANEST_
      character*64 date
      call fdate(date)
#endif
#ifdef MOLPRO
#else
      call datimm()
#endif
      return
      end

c************************************************
      subroutine allocate_int_memory()
c************************************************
#include "drt_h.fh"
c     this subroutine is used to allocate the dynamic memory
c     for integrals
c vint_ci(:) pointer to the base address of the memory of integrals
c maxintseg  maximum length of the integral segment
c
      allocate(vint_ci(0:maxintseg+1))
      vint_ci=0.d0
      return
c...end of allocate_int_memory
      end

c************************************************
      subroutine deallocate_int_memory()
c************************************************
#include "drt_h.fh"

      deallocate(vint_ci)
      return
c...end of deallocate_int_memory
      end


      subroutine read_ml(nf,i,bv,n,m)
      implicit REAL*8 (a-h,o-z)
      dimension bv(n)
      dimension irec(64)

      idisk=0
      call idafile(nf,2,irec,64,idisk)
      idisk=irec(m)
      call ddafile(nf,2,bv,n,idisk)

c      select case (nf)
c        case(34)
c          call ddafile(nf,2,bv,n,idisk)
c        case(35)
c          call ddafile(nf,2,bv,n,idisk)
c        case default

c      end select

      return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(i)
      end

      subroutine write_ml(nf,i,bv,n,m)
      implicit REAL*8 (a-h,o-z)
      dimension bv(n)
      dimension irec(64)

      idisk=0
      if(m.eq.1) then
        irec=0
        call idafile(nf,1,irec,64,idisk)
        irec(1)=idisk
      else
        call idafile(nf,2,irec,64,idisk)
      endif

      idisk=irec(m)
      call ddafile(nf,1,bv,n,idisk)
c
c      select case (nf)
c        case(34)
c          call ddafile(nf,1,bv,n,idisk)
c        case(35)
c          call ddafile(nf,1,bv,n,idisk)
c        case default

c      end select
      irec(m+1)=idisk
      idisk=0
      call idafile(nf,1,irec,64,idisk)
c      print*,"idisk",m,idisk
      return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(i)
      end

      subroutine write_bv(nf,i,bv,n)
      implicit REAL*8 (a-h,o-z)
      dimension bv(n)

      call write_ml(nf,1,bv,n,i)

      return
      end

      subroutine read_bv(nf,i,bv,n)
      implicit REAL*8 (a-h,o-z)
      dimension bv(n)

      call read_ml(nf,1,bv,n,i)

      return
      end

c***************************************
      subroutine readint(ntyp,vintrd)
c***************************************
      implicit REAL*8 (a-h,o-z)
c
c     subroutine used to read divided integrals into
c     main memeroy
c on entry:
c -----------------
c ntyp   - integral type need to be readed
c   1 internal space integrals
c   2 integrals which have one or three indexs in external space
c   3 integrals which have two indexs in external space
c   4 external space integrals
c on out:
c -------------------
c vintrd(*) - integrals used in hamitonian matrix calculation
c
#include "files_gugaci.fh"
      dimension vintrd(*)
      dimension idx(4)
      dimension idum(1)

      idx=0
      idisk=0
      call idafile(luciint,2,idx,4,idisk)
      select case (ntyp)
        case (1)
          idisk=idx(1)
          call idafile(luciint,2,idum,1,idisk)
          lenint=idum(1)
          call ddafile(luciint,2,vintrd(2),lenint,idisk)
        case (2)
c internal external space 1 and 3 index integrals
          idisk=idx(2)
          call idafile(luciint,2,idum,1,idisk)
          lenint=idum(1)
          call ddafile(luciint,2,vintrd(2),lenint,idisk)
          !write(6,*) vintrd(1:lenint)
        case (3)
c interner external space 2 index integrals
          idisk=idx(3)
          call idafile(luciint,2,idum,1,idisk)
          lenint=idum(1)
          call ddafile(luciint,2,vintrd(2),lenint,idisk)
        case (4)
c external space intergrals
          idisk=idx(4)
          call idafile(luciint,2,idum,1,idisk)
          lenint=idum(1)
          call ddafile(luciint,2,vintrd(2),lenint,idisk)
      endselect

      !  do i=1,lenint+1
      !    write(12,"(i8,1x,f15.8)") i,vintrd(i)
      !  enddo
      return
c...end of readint
      end

      subroutine gugaciinit
      !include "comfile.fh"
#ifdef MOLPRO
      use file_qininit
      use molegeom
      use groupinfo
      use orbinfo
#endif
#include "drt_h.fh"
#include "files_gugaci.fh"
#ifdef MOLPRO
      integer :: noffset(maxrecord)
#endif

#ifdef MOLPRO
#else
      fnonemo="traone"
#endif
      fntwomo="traint"
      fndrt="cidrt"
      fnciint="ciint"
      fnloop="ciloop"
      fncidia="cidia"
      fncivec="civec"
      fncitv1="citv1"
      fncitv2="citv2"
      fnciden="ciden"
      fncimo="cimo"

      luonemo=30
      lutwomo=50
      ludrt=31
      luciint=32
      luloop=33
      lucidia=34
      lucivec=35
      lucitv1=36
      lucitv2=37
      luciden=38
      lucimo=39
#ifdef MOLPRO
#else
      call daname(ludrt,fndrt)
      call daname(luciint,fnciint)
      call daname(luloop,fnloop)
      call daname(lucidia,fncidia)
      call daname(lucivec,fncivec)
      call daname(lucitv1,fncitv1)
      call daname(lucitv2,fncitv2)
      call daname(luciden,fnciden)
      call daname(lucimo,fncimo)
#endif
#ifdef _XIANEST_
!open chkfile and drt file
      call fileopen(nfchk,anchk,11)
      call chkfil_taskctrl(2)  ! read task information

      nreps=0
!read infomation from checkfile and moint file
      call chkfil_ciorbinf(2)
      nlsm_bas(1:8)=nsbas(1:8)
      nlsm_all(1:8)=nsorb(1:8)
      ng_sm=nreps

      call fileopen(ludrt,fndrt,12)
      call fileopen(luciint,fnciint,12)
      call fileopen(luloop,fnloop,12)
      call fileopen(lucidia,fncidia,12)
      call fileopen(lucivec,fncivec,12)
      call fileopen(lucitv1,fncitv1,12)
      call fileopen(lucitv2,fncitv2,12)
      call fileopen(luciden,fnciden,12)
      call fileopen(lucimo,fncimo,12)
      call fileopen(lutwomo,fntwomo,12)
#endif

      return
      end

      subroutine gugafinalize
      implicit REAL*8 (a-h,o-z)
#include "files_gugaci.fh"
#ifdef MOLPRO
      call fileclos(ludrt,12)
      call fileclos(luciint,12)
      call fileclos(luloop,12)
      call fileclos(lucidia,12)
      call fileclos(lucivec,12)
      call fileclos(lucitv1,12)
      call fileclos(lucitv2,12)
      call fileclos(luciden,12)
      call fileclos(lucimo,12)
      call fileclos(lutwomo,12)
#else
      call daclos(ludrt)
      call daclos(luciint)
      call daclos(luloop)
      call daclos(lucidia)
      call daclos(lucivec)
      call daclos(lucitv1)
      call daclos(lucitv2)
      call daclos(luciden)
      call daclos(lucimo)
#endif
      return
      end

      subroutine memdengrad_alloc()
#include "drt_h.fh"
#include "grad_h.fh"

      ndim=norb_all*(norb_all+1)/2
      allocate(denm1(ndim))
      nc=ndim*(ndim+1)/2
      allocate(denm2(nc))

      end

      subroutine memdengrad_free()
#include "drt_h.fh"
#include "grad_h.fh"

      deallocate(denm1)
      deallocate(denm2)

      end

      subroutine memcidiag_alloc()
#include "drt_h.fh"
      integer, pointer :: jph(:),jeh(:),jwh(:)
      REAL*8, pointer :: th(:),thh(:)
      common/ptlph/jph,jeh,jwh
      common/ptlphv/th,thh

      allocate(jph(maxpl))
      allocate(jeh(maxpl))
      allocate(jwh(maxpl))
      allocate(th(maxpl))
      allocate(thh(maxpl))
      jph=0
      jeh=0
      jwh=0
      th=0
      thh=0

      return
      end

      subroutine memcidiag_dealloc()
#include "drt_h.fh"
      integer, pointer :: jph(:),jeh(:),jwh(:)
      REAL*8, pointer :: th(:),thh(:)
      common/ptlph/jph,jeh,jwh
      common/ptlphv/th,thh

      deallocate(jph)
      deallocate(jeh)
      deallocate(jwh)
      deallocate(th)
      deallocate(thh)

      return
      end

      subroutine mem_intinnindex_alloc()
#include "drt_h.fh"
#include "intsort_h.fh"

      allocate(loij(50000))
      allocate(loijk(1384150))
      allocate(loij_all(50000))
      allocate(loijk_all(1384150))
      loij(1:50000)=0
      loijk(1:1384150)=0
      loijk_all(1:50000)=0
      loijk_all(1:1384150)=0

      allocate(intind_iaqq(50000))
      allocate(intind_abkk(50000))
      lent=norb_inn*nabc+norb_all
     *    +ngw2(norb_all)+ngw3(norb_all) !2*norb_ext*norb_ext*norb_ext
      allocate(intind_iabc(lent))
      allocate(intind_ijka(50000))
      allocate(intind_ijcc(50000))
      allocate(intind_ijab(50000))
      intind_iaqq(1:50000)=0
      intind_abkk(1:50000)=0
      intind_iabc(1:lent)=0
      intind_ijka(1:50000)=0
      intind_ijcc(1:50000)=0
      intind_ijab(1:50000)=0

      allocate(intspace_iaqq(50000))
      allocate(intspace_abkk(50000))
      allocate(intspace_iabc(50000))
      allocate(intspace_ijka(50000))
      allocate(intspace_ijcc(50000))
      allocate(intspace_ijab(50000))
      intspace_iaqq(1:50000)=0
      intspace_abkk(1:50000)=0
      intspace_iabc(1:50000)=0
      intspace_ijka(1:50000)=0
      intspace_ijcc(1:50000)=0
      intspace_ijab(1:50000)=0

      voint=0.d0
      vdint=0.d0

      return
      end

      subroutine mem_intinnindex_dealloc()
#include "drt_h.fh"
#include "intsort_h.fh"

      deallocate(loij)
      deallocate(loijk)
      deallocate(loij_all)
      deallocate(loijk_all)

      deallocate(intind_iaqq)
      deallocate(intind_abkk)
      deallocate(intind_iabc)
      deallocate(intind_ijka)
      deallocate(intind_ijcc)
      deallocate(intind_ijab)

      deallocate(intspace_iaqq)
      deallocate(intspace_abkk)
      deallocate(intspace_iabc)
      deallocate(intspace_ijka)
      deallocate(intspace_ijcc)
      deallocate(intspace_ijab)

      return
      end

      subroutine allocate_casrst()
#include "drt_h.fh"
#include "pl_structure_h.fh"

      allocate(ja(max_node),jb(max_node),jm(0:max_node))
      allocate(jj(1:4,0:max_node))
      allocate(kk(0:max_node))
      ja=0; jb=0; jm=0
      jj=0
      kk=0
      return
      end

      subroutine deallocate_casrst()
#include "drt_h.fh"
#include "pl_structure_h.fh"

      deallocate(ja,jb,jm)
      deallocate(jj)
      deallocate(kk)
      return
      end

      subroutine allocate_subdrt(icase,lent)
#include "drt_h.fh"
#include "pl_structure_h.fh"

      allocate(ihy(max_wei),jj_sub(1:4,0:max_node))
      allocate(iy(1:4,0:max_node))
      if(icase.eq.1) then
        allocate(jphy(max_node))
      else
        allocate(jphy(lent))
      endif
      end

      subroutine allocate_subdrtl(icase,lent)
#include "drt_h.fh"
#include "pl_structure_h.fh"

      allocate(ihyl(max_wei),jjl_sub(1:4,0:max_node))
      allocate(iyl(1:4,0:max_node))
      if(icase.eq.1) then
        allocate(jphyl(max_node))
      else
        allocate(jphyl(lent))
      endif
      end

      subroutine deallocate_subdrt()
#include "drt_h.fh"
#include "pl_structure_h.fh"

      deallocate(ihy,jj_sub)
      deallocate(iy)
      deallocate(jphy)
      end

      subroutine deallocate_subdrtl()
#include "drt_h.fh"
#include "pl_structure_h.fh"

      deallocate(ihyl,jjl_sub)
      deallocate(iyl)
      deallocate(jphyl)
      end

#ifdef MOLPRO
#else
      real*8 function c_time()
      c_time=seconds()

      end function c_time

      integer function ipair(i,j)
      implicit none
      integer :: i,j
      if(i.ge.j) then
        ipair=i*(i-1)/2+j
      else
        ipair=j*(j-1)/2+i
      endif
      end function ipair

      subroutine trimstr(string)
! delete space character in the head and tail of the string
      implicit none
      character*(*),intent(out) :: string
      character*128  :: line
      integer :: i,j,k,mylen
      k=mylen(string)
      line(1:k)=string(1:k)
      do i=1,k
        if(string(i:i).ne." ") exit
      enddo
      string=" "
      j=k-i+1
      string(1:j)=line(i:k)

      return
      end
#endif

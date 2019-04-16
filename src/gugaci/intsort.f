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
      subroutine int_sort()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "files_gugaci.fh"
      dimension idx(4)
      data zero/0.d0/

      idx=0
      stime=c_time()

!estimate the largest block of integrals, notice we do not use symmetry
!here
      if(norb_ext.gt.norb_inn) then
! (ab|cd)
        numb=norb_ext*(norb_ext+1)/2
        numb=numb*(numb+1)/2+norb_ext*numb
      else
! (ii|ab)
        numb=norb_inn*(norb_inn+1)/2
        numb=numb*norb_ext*(norb_ext+1)/2
      endif
! (ia|bc)+(ij|ka)
      lra=norb_inn*norb_ext*(norb_ext+1)*norb_ext/2
     &   +norb_inn*(norb_inn+1)*norb_inn*norb_ext/2
      if(lra.gt.numb) numb=lra

      allocate(vint_ci(numb))
      vint_ci=0.d0
      numb=0
      write(6,900)
      call int_index(numb)
#ifdef MOLPRO
      call intrd
#else
      call intrd_molcas
#endif
c --------     vtint(*)-->vint_ci(*)    ------------
      numb=1
      call int_sort_inn(numb)         !_ext_0  (ijkl,ijkk,iijj)
      idisk=0
      call idafile(luciint,1,idx,4,idisk)
      idx(1)=idisk
      numb=numb-1
      call idafile(luciint,1,[numb],1,idisk)
      call ddafile(luciint,1,vint_ci,numb,idisk)
      write(6,902) numb
      maxintseg=numb

      vint_ci=0.d0
      numb=1
      call int_sort_inn_1(numb)       !_ext_3_1  (iabc and iaqq)
      call int_sort_inn_3(numb)       !_ext_1    (ijka)
      idx(2)=idisk
      numb=numb-1
      call idafile(luciint,1,[numb],1,idisk)
      call ddafile(luciint,1,vint_ci,numb,idisk)
      write(6,903) numb
      if(numb.gt.maxintseg) then
        maxintseg=numb
      endif

c  sum over iaqq
      viasum_0(1:norb_dz,1:norb_ext)=zero
      viasum_1(1:norb_dz,1:norb_ext)=zero
      do lri=norb_frz+1,norb_inn
        ia0=(lri-1)*norb_ext
        do lra=1,norb_ext
          ia=ia0+lra
          intposbase=intind_iaqq(ia)
          if(intposbase.eq.0) cycle
          do lrk=1,norb_dz
            idorbint=lrk*2-2
            intpos=intposbase+idorbint
            viasum_0(lri,lra)=viasum_0(lri,lra)+vint_ci(intpos)
            viasum_1(lri,lra)=viasum_1(lri,lra)+vint_ci(intpos+1)
          enddo
        enddo
      enddo

      vint_ci=0.d0
      numb=1
      call int_sort_inn_2(numb)   !_ext_2_1  (ijcc, ijab, abkk)
      idx(3)=idisk
      numb=numb-1
      call idafile(luciint,1,[numb],1,idisk)
      call ddafile(luciint,1,vint_ci,numb,idisk)
      write(6,904) numb
      if(numb.gt.maxintseg) then
        maxintseg=numb
      endif

c  sum over (abkk)
      intspace=intspace_abkk(1)
      intpos=intind_abkk(1)
      vijkk_0sum(1:intspace)=zero
      vijkk_1sum(1:intspace)=zero
      do lrk=1,norb_dz
        do ivalue=1,intspace
          vijkk_0sum(ivalue)=vijkk_0sum(ivalue)+vint_ci(intpos)
          vijkk_1sum(ivalue)=vijkk_1sum(ivalue)+vint_ci(intpos+1)
          intpos=intpos+2
        enddo
      enddo

      vint_ci=0.d0
      numb=1
      call int_sort_ext(numb)         !_ext_4_3_2  (abcd,abcc,aabb)
      idx(4)=idisk
      numb=numb-1
      call idafile(luciint,1,[numb],1,idisk)
      call ddafile(luciint,1,vint_ci,numb,idisk)
      write(6,905) numb
      if(numb.gt.maxintseg) then
        maxintseg=numb
      endif

      idisk=0
      call idafile(luciint,1,idx,4,idisk)

      write(6,910)
      deallocate(vint_ci)

      etime=c_time()
      time=etime-stime
      write(6,912) time

900   format(1x,"start integral sorting",/,
     *       1x,"start reading integral file",/)
c901   format(1x,"end reading integral file",/)
902   format(1x,"number of integrals in internal space is ",i8)
903   format(1x,"number of integrals which have one or three "
     *       ,"indexs in external space is ",i8)
904   format(1x,"number of integrals which have two indexs in "
     *       ,"external space is ",i8)
905   format(1x,"number of integrals in external space is ",i8)
910   format(1x,/,1x,"end of integral sorting")
912   format(1x,"total wall clock time for integral sorting=",
     *       f8.2," s",/)
      end

      subroutine blocks()
c**************************************
#include "drt_h.fh"
#include "intsort_h.fh"
c                              ----- 1 - electron integrals -----
      nint1=0
      do 1 ip=1,ng_sm
        nint1=nint1+((nlsm_all(ip)+1)*nlsm_all(ip))/2
  1   continue
c
c                              ----- 2 - electron integrals -----
      nint2=0
      nblock=0
      do 10 ip=1,ng_sm
        npp=(nlsm_all(ip)+1)*nlsm_all(ip)/2
        nintb=npp*(npp-1)/2
        nint2=nint2+nintb
        nblock=nblock+1
        do 11 i=1,4
  11    iblktb(i,nblock)=ip
        iblktb(5,nblock)=nint2
  10  continue
      do 20 ip=1,ng_sm
        npp=(nlsm_all(ip)+1)*nlsm_all(ip)/2
        irm=ip
      do 20 ir=1,irm-1
        nrr=(nlsm_all(ir)+1)*nlsm_all(ir)/2
        nintb=npp*nrr
        nint2=nint2+nintb
        nblock=nblock+1
        iblktb(1,nblock)=ip
        iblktb(2,nblock)=ip
        iblktb(3,nblock)=ir
        iblktb(4,nblock)=ir
        iblktb(5,nblock)=nint2
  20  continue
      do 30 ip=1,ng_sm
        iqm=ip
      do 30 iq=1,iqm-1
        ipq=nlsm_all(ip)*nlsm_all(iq)
        nintb=ipq*(ipq-1)/2
        nint2=nint2+nintb
        nblock=nblock+1
        iblktb(1,nblock)=ip
        iblktb(2,nblock)=iq
        iblktb(3,nblock)=ip
        iblktb(4,nblock)=iq
        iblktb(5,nblock)=nint2
  30  continue
      do 40 ip=4,8
        if(nlsm_all(ip).eq.0) goto 40
          iqm=ip
        do 41 iq=1,iqm-1
          if(nlsm_all(iq).eq.0) goto 41
          npq=nlsm_all(ip)*nlsm_all(iq)
          ispq=mul_tab(ip,iq)
          irm=ip
          do 42 ir=1,irm-1
            if(nlsm_all(ir).eq.0) goto 42
            ispqr=mul_tab(ispq,ir)
            ism=ir
            if(ip.eq.ir) ism=iq
            do 43 is=1,ism-1
              if(nlsm_all(is).eq.0) go to 43
              if(is.ne.ispqr) go to 43
              nrs=nlsm_all(ir)*nlsm_all(is)
              nintb=npq*nrs
              nint2=nint2+nintb
              nblock=nblock+1
              iblktb(1,nblock)=ip
              iblktb(2,nblock)=iq
              iblktb(3,nblock)=ir
              iblktb(4,nblock)=is
              iblktb(5,nblock)=nint2
  43        continue
  42      continue
  41    continue
  40  continue
      nint12=nint1+nint2

      write(6,100) nint1,nint2,nint12
c      write(12,100) nint1,nint2,nint12
  100 format(' ',1x/2x,
     -  'number of 1-electron integrals  :',i9/2x,
     -  'number of 2-electron integrals  :',i9/2x,
     -  'total number of integrals       :',i9)
      write(6,200) nint1
c      write(12,200) nint1
  200 format(' ',1x/2x,
     -  '1-electron blocks  :   1 to',i8/2x,29(1h*))
      write(6,300) (j,(iblktb(i,j),i=1,5),j=1,nblock)
c      write(12,300) (j,(iblktb(i,j),i=1,5),j=1,nblock)
  300 format(' ',1x/2x,
     -  '2-electron block description  :'/2x,
     -      40(1h*)/2x,50(3('(',i3,')',4i2,i8,3x)/2x))
      return
      end

      subroutine ff(i,j)
       i=0
       do 10 iq=1,j
       j0=iq*(iq-1)/2+iq
       if(j0.ne.j)goto 10
       i=iq
       goto 100
10     continue
100    return
       end

      subroutine int_sort_ext(ii)         !_ext_4_3_2
#include "drt_h.fh"
#include "intsort_h.fh"
      dimension lsmtmp(8)
!  _002_aa_
      ip2_aa_ext_base=ii
      do lsma=1,ng_sm
        iasta=ibsm_ext(lsma)
        ibsta=iasta+1
        ibend=iesm_ext(lsma)
        do ib=ibsta,ibend
          jb=norb_number(ib)
          do ia=iasta,ib-1
            ja=norb_number(ia)
            vint_ci(ii)=voint(jb,ja)          !(ja|jb)
!            write(10,'(2x,a20,i8,f16.8)')
!     *            '  int_sort_002',ii,vint_ci(ii)
            ii=ii+1
          enddo
        enddo
      enddo
!  _002_dd_
      ip2_dd_ext_base=ii
c     write(10,'(2x,a20,i8)')' ip2_dd_ext_base=',      ii
      do ib=2,norb_ext
        jb=norb_number(ib)
        do ia=1,ib-1
          ja=norb_number(ia)
          vint_ci(ii)=voint(ja,jb)   !(ja,jb|ja,jb)
!       write(6,'(i2,2x,i2)') ja,jb
!       write(6,'(2x,a20,i8,f16.8)')'  int_sort_002',ii,vint_ci(ii)
          ii=ii+1
        enddo
      enddo
!_003
      ip3_abd_ext_base=ii
c     write(10,'(2x,a20,i8)')' ip3_abd_ext_base=',      ii

      do ic=1,norb_ext
        lrc=norb_number(ic)
!       ip3_abd_ext_base(ic)=ii
      do lsma=1,ng_sm
        iasta=ibsm_ext(lsma)
!       iaend=iesm_ext(lsma)-1
        ibsta=iasta+1
        ibend=iesm_ext(lsma)
        do ib=ibsta,ibend
          lrb=norb_number(ib)
          do ia=iasta,ib-1
            lra=norb_number(ia)
            vint_ci(ii)=vfutei(lra,lrc,lrb,lrc)
            vint_ci(ii+1)=vfutei(lra,lrb,lrc,lrc)
!      write(10,'(2x,4i6,i8,2f16.8)')
!     :   lra,lrb,lrc,lrc,ii,vint_ci(ii),vint_ci(ii+1)
            ii=ii+2
          enddo
        enddo
      enddo
      enddo
      np3_abd_ext=0
      do isma=1,ng_sm
        na=nlsm_ext(isma)
        if ( na .le. 0 ) cycle
          np3_abd_ext=np3_abd_ext+na*(na-1)
      enddo
!_004
!      ==================================
!          la<lb<lc<ld
!      ===================================
!       1. ext _abcd_
!       ip4_abcd_ext_base=ii
      do lsmd=1,ng_sm
        idsta=ibsm_ext(lsmd)
        idend=iesm_ext(lsmd)
        do lsmc=1,lsmd
          icsta=ibsm_ext(lsmc)
          icend=iesm_ext(lsmc)
          if ( lsmc .eq. lsmd ) idsta=idsta+1
            lsmcd=mul_tab(lsmc,lsmd)
            do lsmb=1,lsmc
              lsma=mul_tab(lsmb,lsmcd)
              if ( lsma .gt. lsmb ) cycle
                ibsta=ibsm_ext(lsmb)
                ibend=iesm_ext(lsmb)
                iasta=ibsm_ext(lsma)
                iaend=iesm_ext(lsma)
              if ( lsmd .eq. lsmb ) idsta=idsta+1
              if ( lsmc .eq. lsmb ) icsta=icsta+1
              if ( lsmd .eq. lsma ) idsta=idsta+1
              if ( lsmc .eq. lsma ) icsta=icsta+1
              if ( lsmb .eq. lsma ) ibsta=ibsta+1
                lsmtmp=lsma+jp2(lsmb)+jp3(lsmc)
                ip4_abcd_ext_base(lsmtmp)=ii
c     write(10,'(2x,a20,2i8)')' ip4_abcd_ext_base=',  lsmtmp,ii
      do id=idsta,idend
        lrd=norb_number(id)
        do ic=icsta,min(id-1,icend)
          lrc=norb_number(ic)
          do ib=ibsta,min(ic-1,ibend)
            lrb=norb_number(ib)
            do ia=iasta,min(ib-1,iaend)
              lra=norb_number(ia)
              vint_ci(ii)=vfutei(lra,lrb,lrc,lrd)       !lra>lrb>lrc>lrd
              vint_ci(ii+1)=vfutei(lra,lrc,lrb,lrd)
              vint_ci(ii+2)=vfutei(lra,lrd,lrc,lrb)
              ii=ii+3
            enddo
          enddo
        enddo
      enddo

          enddo
        enddo
      enddo
      end

!===================== ext_2_1 sta ===================
      subroutine int_sort_inn_2(ii)
#include "drt_h.fh"
#include "intsort_h.fh"
c     (ijcc,ijab)
      do ismab=1,ng_sm
        ip2_ab_inn_base(ismab)=ii
        do lri=norb_frz+1,norb_inn-1
          lsmi=lsm_inn(lri)
          do lrj=lri+1,norb_inn
            lsmj=lsm_inn(lrj)
            lsmij=mul_tab(lsmi,lsmj)
            if ( ismab .ne. lsmij ) cycle
            call int_ext_2_1(lri,lrj,lsmij,ii)
          enddo
        enddo
      enddo

c     (abkk)
      do lri=1,norb_inn
        intind_abkk(lri)=ii
        intspace_abkk(lri)=0
        do lmb=1,ng_sm
          ibsta=ibsm_ext(lmb)
          ibend=iesm_ext(lmb)
          do irb=ibsta,ibend
            lrb=norb_number(irb)
            do ira=ibsta,irb-1
              lra=norb_number(ira)
              intspace_abkk(lri)=intspace_abkk(lri)+1
              vint_ci(ii)=vfutei(lra,lri,lrb,lri)
              vint_ci(ii+1)=vint_ci(ii)-2*vfutei(lra,lrb,lri,lri)
!       write(6,'(2x,4i6,i8,3f16.8)')lra,lrb,lri,lri,ii,
!     *        vint_ci(ii),vint_ci(ii+1)
              ii=ii+2
            enddo
          enddo
        enddo
      enddo

      return
      end

      subroutine int_ext_2_1(lri,lrj,lsmij,ii)
#include "drt_h.fh"
#include "intsort_h.fh"
c      write(10,*)'     start intind_ijab',lri,lrj,ii
      ij=lri-norb_frz+ngw2(lrj-norb_frz)
      if ( lsmij .eq. 1 ) then       !ijcc
c        write(10,*)'     start intind_ijcc',ii
        intind_ijcc(ij)=ii
        intspace_ijcc(ij)=0
        do ic=1,norb_ext
          intspace_ijcc(ij)=intspace_ijcc(ij)+1
          lrc=norb_number(ic)
          vint_ci(ii)=vfutei(lrj,lrc,lri,lrc)
          vint_ci(ii+1)=vfutei(lrj,lri,lrc,lrc)
c      write(10,'(2x,4i6,i8,3f16.8)')lrj,lri,lrc,lrc,ii,
c     *          vint_ci(ii),vint_ci(ii+1)
          ii=ii+2
        enddo
      endif
      intind_ijab(ij)=ii
      intspace_ijab(ij)=0
      do lsmc=1,ng_sm
        lsmd=mul_tab(lsmij,lsmc)
        if ( lsmd .gt. lsmc ) cycle
        icsta=ibsm_ext(lsmc)
        icend=iesm_ext(lsmc)
        idsta=ibsm_ext(lsmd)
        idend=iesm_ext(lsmd)
        if ( lsmc .eq. lsmd ) icsta=icsta+1
        do ic=icsta,icend
          jc=norb_number(ic)
          do id=idsta,min(idend,ic-1)
            intspace_ijab(ij)=intspace_ijab(ij)+1
            jd=norb_number(id)
            vint_ci(ii)=vfutei(jd,jc,lrj,lri)
            vint_ci(ii+1)=vfutei(jd,lrj,jc,lri)
            vint_ci(ii+2)=vfutei(jd,lri,lrj,jc)
c      write(10,'(2x,4i6,i8,3f16.8)')jd,jc,lrj,lri,ii,
c     *          vint_ci(ii),vint_ci(ii+1),vint_ci(ii+2)
            ii=ii+3
          enddo
        enddo
      enddo
      return
      end
!===================== ext_2_1 end ===================
!===================== ext_1 sta ======================
      subroutine int_sort_inn_3(ii)
#include "drt_h.fh"
#include "intsort_h.fh"
c     write(10,'(2x,a20,i8)')' int_inn_3_base=',      ii
c      write(10,*)'     start intind_ijka',ii
      do lri=norb_frz+1,norb_inn-2
        lsmi=lsm_inn(lri)
        do lrj=lri+1,norb_inn-1
          lsmj=lsm_inn(lrj)
          lsmij=mul_tab(lsmi,lsmj)
          do lrk=lrj+1,norb_inn
            lsmk=lsm_inn(lrk)
            lsmd=mul_tab(lsmij,lsmk)
            ijk=lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lrk-norb_frz)
            intind_ijka(ijk)=ii
            iasta=ibsm_ext(lsmd)
            iaend=iesm_ext(lsmd)
            do ia=iasta,iaend
              ja=norb_number(ia)
              vint_ci(ii)=vfutei(ja,lrk,lrj,lri)
              vint_ci(ii+1)=vfutei(ja,lrj,lrk,lri)
              vint_ci(ii+2)=vfutei(ja,lri,lrk,lrj)
c     write(10,'(2x,4i6,i8,3f16.8)')ja,lrk,lrj,lri,ii,
c    *       vint_ci(ii),vint_ci(ii+1),vint_ci(ii+2)
              ii=ii+3
            enddo
          enddo
        enddo
      enddo
      end

!===================== _ext_1 end ======================
      subroutine int_ext_3_2_1(lri,lsmi,ii)
#include "drt_h.fh"
#include "intsort_h.fh"
c      write(10,*)'   int_ext_3_2_1'
      iabc0=(lri-1)*nabc
c     write(10,*)'     start intind_iabc',ii
      do lsmd=1,ng_sm
        lsmbc=mul_tab(lsmi,lsmd)
        do lsmc=1,lsmd
          lsmb=mul_tab(lsmbc,lsmc)
          if ( lsmb .gt. lsmc ) cycle
            idsta=ibsm_ext(lsmd)
            idend=iesm_ext(lsmd)
            icsta=ibsm_ext(lsmc)
            icend=iesm_ext(lsmc)
            ibsta=ibsm_ext(lsmb)
            ibend=iesm_ext(lsmb)
            if ( lsmb .eq. lsmc ) icsta=icsta+1
            if ( lsmc .eq. lsmd ) idsta=idsta+1
            if ( lsmb .eq. lsmd ) idsta=idsta+1
            do ird=idsta,idend
              lrd=norb_number(ird)
              nd=ngw3(ird)
              do irc=icsta,min(icend,ird-1)
                lrc=norb_number(irc)
                nc=ngw2(irc)
                do irb=ibsta,min(ibend,irc-1)
                  lrb=norb_number(irb)
                  nb=irb
                  iabc=iabc0+nb+nc+nd
                  intind_iabc(iabc)=ii
                  vint_ci(ii)=vfutei(lrb,lrc,lrd,lri)
                  vint_ci(ii+1)=vfutei(lrb,lrd,lrc,lri)
                  vint_ci(ii+2)=vfutei(lrb,lri,lrd,lrc)
                  ii=ii+3
                enddo
              enddo
            enddo
          enddo
        enddo

c      write(10,*)'     start intind_iaqq',ii
      ib0=(lri-1)*norb_ext
      lsmb=lsmi                       !3_ibqq
      ibsta=ibsm_ext(lsmb)
      ibend=iesm_ext(lsmb)
      do irb=ibsta,ibend
        lrb=norb_number(irb)
        ib=ib0+irb
        intind_iaqq(ib)=ii
        do lrc=1,norb_inn
          vint_ci(ii)=vfutei(lrb,lrc,lri,lrc)
          vint_ci(ii+1)=vfutei(lrb,lri,lrc,lrc)
          ii=ii+2
        enddo
        do irc=norb_ext,1,-1
          lrc=norb_number(irc)
          vint_ci(ii)=vfutei(lrb,lrc,lri,lrc)
          vint_ci(ii+1)=vfutei(lrb,lri,lrc,lrc)
!      write(6,'(2x,4i6,i8,3f16.8)')lrb,lri,lrc,lrc,ii,
!     *      vint_ci(ii),vint_ci(ii+1)
          ii=ii+2
        enddo
      enddo

      end

      subroutine int_sort_inn_1(ii)
#include "drt_h.fh"
#include "intsort_h.fh"

      call determine_para_array_for_int1ind()

      do lri=norb_frz+1,norb_inn
        li=lri-norb_frz
        lsmi=lsm_inn(lri)
        call int_ext_3_2_1(lri,lsmi,ii)
      enddo
      end

!===================== ext_3_2_1 end ===================

      function list3(i,j,k)
#include "drt_h.fh"
#include "intsort_h.fh"
c            *****************
      nij   = i+ngw2(j)
      list3 = loij(nij)+2*(k-1)
c            *****************
      return
      end
c***********************************************************************
      function list4(ld,lc,lb,la)
#include "drt_h.fh"
#include "intsort_h.fh"
c                    ***************
      lra  = ncibl(la)
      njkl = ld+ngw2(lc)+ngw3(lb)
      list4= loijk(njkl)+3*(lra-1)
c                    ***************
      return
      end
************************************************************************
      subroutine int_sort_inn(numb)
#include "drt_h.fh"
#include "intsort_h.fh"
      dimension msob(120)

      msob=0
      do 12 lra=norb_inn,1,-1
        ms = lsm_inn(lra)
        msob(ms)=msob(ms)+1
        ncibl(lra) = msob(ms)
12    continue
c=======================================================================
!      write(6,'(1x,14i3)')(ncibl(i),i=1,norb_inn)

      numb = 1
      do 10  i = 1,norb_inn-1
        do 10  j = i+1,norb_inn
          if(lsm_inn(i).ne.lsm_inn(j)) goto 10

          nij=i+ngw2(j)
          loij(nij)=numb
          do 20  k = 1,norb_inn
            vint_ci(numb)=vfutei(j,k,i,k)
            vint_ci(numb+1)=vfutei(j,i,k,k)
!        write(6,'(2x,4i6,i8,3f16.8)')i,j,k,k, numb,
!     *        vint_ci(numb),vint_ci(numb+1)
            numb=numb+2
20        continue
10    continue
!      write(6,*) "num_inn",numb
!      print*, loij(1:8)
!      stop 888

c=======================================================================
c      la<lb<lc<ld
      do 30 ld = 1,norb_inn-3
        do 30 lc = ld+1,norb_inn-2
          msd  = lsm_inn(ld)
          msc  = lsm_inn(lc)
          mscd = mul_tab(msd,msc)
          do 30 lb = lc+1,norb_inn-1
            msb = lsm_inn(lb)
            msa = mul_tab(mscd,msb)

            njkl=ld+ngw2(lc)+ngw3(lb)
            loijk(njkl) = numb

            nolra=0
            do 40 la = norb_inn,lb+1,-1
              if(lsm_inn(la).ne.msa) cycle
              nolra=nolra+1
              list = loijk(njkl)+3*(nolra-1)

c        write(6,'(2x,4i3,2i7)')  la,lb,lc,ld,list,numb

              vint_ci(numb)=vfutei(la,lc,lb,ld)        !tmp stop
              vint_ci(numb+1)=vfutei(la,lb,lc,ld)
              vint_ci(numb+2)=vfutei(la,ld,lc,lb)
c     write(10,'(2x,4i6,i8,3f16.8)')la,lb,lc,ld, numb,
c    *        vint_ci(numb),vint_ci(numb+1),vint_ci(numb+2)
              numb=numb+3
40          continue
30    continue
      return
c=======================================================================
      end

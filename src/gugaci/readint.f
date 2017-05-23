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
c  16 apr 2007 - bsuo - revised to use molcas intergrals
#ifdef MOLPRO
      subroutine intrd
      use file_qininit, only : maxrecord
#include "drt_h.fh"
#include "intsort_h.fh"
#include "files_gugaci.fh"
      integer :: noffset(maxrecord)
      common /mcorb/ lsmorb(max_orb),noidx(8)
      REAL*8, pointer :: x(:)
      dimension xfock(max_orb*(max_orb+1)/2)
      dimension itratoc(ntratoc)


      nintone=0
      nmob=0
      nidx=0
      ni=0
      do i=1,ng_sm
        nism   = nlsm_all(i)
        nsmint = nism*(nism+1)/2
        itratoc(i)=nism*nism
        nmob=nmob+nlsm_bas(i)*nlsm_bas(i)
        noidx(i)=nidx
        nidx=nidx+nism
        nintone=nintone+nsmint
        lsmorb(ni+1:nidx)=i
        ni=nidx
      enddo

      noffset=0
      ni=0
      call idafile(Lutwomo,2,noffset,maxrecord,ni)
      ! 2th record in file traint, nuclear repulesive energy
      call ddafile(Lutwomo,2,vpotnuc,1,ni)
      ! 3rd record, ecore
      call ddafile(Lutwomo,2,ecor, 1,ni)  ! effective core energy
      write(6,"(a,1x,f14.8)") " Nuclear repulsive energy:",vpotnuc
      write(6,"(a,1x,f14.8)") " Frozen  core energy:     ",ecor
      vpotnuc=vpotnuc+ecor
      ! 4th record, 1e int
      call ddafile(Lutwomo,2,xfock,nintone,ni) ! 1e integrals

c write one electron fock matrix into voint
      nidx=0
      do i=1,ng_sm
        nism  = nlsm_all(i)
        idx   = noidx(i)
        nsmint = nism*(nism+1)/2
        nc=0
        do lri=1,nism
          do lrj=1,lri
            nc=nc+1
            lrcii=map_orb_order(lri+idx)
            lrcij=map_orb_order(lrj+idx)
            if(lrcii.gt.lrcij) then
              lrt=lrcij
              lrcij=lrcii
              lrcii=lrt
            endif
            voint(lrcii,lrcij)=xfock(nc+nidx)
            !write(6,'(1x,i3,1x,i3,1x,f18.9)') lrcii,lrcij,xfock(nc+nidx)
          enddo
        enddo
        nidx=nidx+nsmint
      enddo
      call readtwoeint(lutwomo,maxrecord,noffset,nlsm_all,ng_sm,mul_tab,
     *                 map_orb_order,noidx)

      !write(6,*)
      !write(6,*) "MRCI integrals"
      !do lri=1,406
      !  write(55,"(i8,f16.8)") lri,vector1(lri)
      !enddo

      return
c...end of intrd_molcas
      end

      subroutine readtwoeint(nft,maxrecord,noffset,norb,
     *                       ngsm,multab,maporb,noidx)
      implicit REAL*8 (a-h,o-z)
#include "ci_parameter.fh"
      integer :: noffset(maxrecord)
      parameter ( kbuf = ntrabuf )
      dimension norb(8),multab(8,8),maporb(max_orb),noidx(8)
      dimension itratoc(ntratoc)
      dimension buff(kbuf)

      idisk=noffset(5)
      write(6,2000)
2000  format(/7x,'symmetry',6x,' orbitals',8x,'integrals')
      do nsp=1,ngsm
        nop=norb(nsp)
        do nsq=1,nsp
          noq=norb(nsq)
          nspq=multab(nsp,nsq)
          do nsr=1,nsp
            nor=norb(nsr)
            nspqr=multab(nspq,nsr)
            nssm=nsr
            if(nsr.eq.nsp) nssm=nsq
            do nss=1,nssm
              if(nspqr.ne.nss) cycle
              nos=norb(nss)

              ityp=0
              if(nsr.eq.nss) then
                nbpq=(nop+nop**2)/2
                nbrs=(nos+nos**2)/2
                if(nsp.eq.nsr) then
c  (ii|ii) type 1 int
                  nintb=(nbpq+nbpq**2)/2
                  ityp=1
                else
c  (ii|jj) type 3 int
                  nintb=nbpq*nbrs
                  ityp=3
                endif
              else
                nbpq=nop*noq
                nbrs=nor*nos
                if(nsp.eq.nsr) then
c (ij|ij) type 2 int
                  nintb=(nbpq+nbpq**2)/2
                  ityp=2
                else
c (ij|kl) type 4 int
                  nintb=nbpq*nbrs
                  ityp=4
                endif
              endif

              if(nintb.eq.0) goto 10
              write(6,2100) nsp,nsq,nsr,nss,nop,noq,nor,nos,
     *                      nintb
2100          format(7x,4i2,1x,4i4,2x,3x,i9)

              idx=0
              iout=0
c              call ddatard(nft,buff,kbuf,idisk)
              if(nintb.gt.kbuf) then
                call ddafile(nft,2,buff,kbuf, idisk)
              else
                call ddafile(nft,2,buff,nintb,idisk)
              endif

              do li=1,nos
                ntj=nor
                if(ityp.eq.1.or.ityp.eq.3) ntj=li
                do lj=1,ntj
                  ntk=noq
                  if(ityp.eq.1.or.ityp.eq.2) ntk=li
                  do lk=1,ntk
                    ntl=nop
                    if(ityp.eq.1.or.ityp.eq.3) ntl=lk
                    if(ityp.eq.1.and.li.eq.lk) ntl=lj
                    if(ityp.eq.2.and.li.eq.lk) ntl=lj
                    do ll=1,ntl
                      iout=iout+1
                      if(iout.gt.kbuf) then
                         if(nintb-idx.lt.kbuf) then
                           call ddafile(nft,2,buff,nintb-idx,idisk)
                         else
                           call ddafile(nft,2,buff,kbuf,idisk)
                         endif
                         iout=1
                      endif
                      idx=idx+1
                      val=buff(iout)
                      lri=maporb(li+noidx(nss))
                      lrj=maporb(lj+noidx(nsr))
                      lrk=maporb(lk+noidx(nsq))
                      lrl=maporb(ll+noidx(nsp))
                      !if(ityp.ne.2) cycle
!                      write(6,"(9(1x,i2),1x,f18.9)") li,lj,lk,ll,
!     *                                               lri,lrj,lrk,lrl,
!     *                                               iout,val
                      call intrw_mol(lri,lrj,lrk,lrl,val)
                    enddo
                  enddo
                enddo
              enddo

10          continue
            enddo
          enddo
        enddo
      enddo

      return
c....end of readtwoeint
      end

#else
      subroutine intrd_molcas
#include "drt_h.fh"
#include "intsort_h.fh"
#include "files_gugaci.fh"
      common /mcorb/ lsmorb(max_orb),noidx(8)
      REAL*8, pointer :: x(:)
      dimension xfock(max_orb*(max_orb+1)/2)
      dimension itratoc(ntratoc)

      nintone=0
      nmob=0
      nidx=0
      ni=0
      do i=1,ng_sm
        nism   = nlsm_all(i)
        nsmint = nism*(nism+1)/2
        itratoc(i)=nism*nism
        nmob=nmob+nlsm_bas(i)*nlsm_bas(i)
        noidx(i)=nidx
        nidx=nidx+nism
        nintone=nintone+nsmint
        lsmorb(ni+1:nidx)=i
        ni=nidx
      enddo
      allocate(x(nmob))
      call daname(luonemo,fnonemo)
      call readtraonehead(luonemo,ecor,idisk)
      vpotnuc=ecor

c  read mo coeff, need debug, if frozen and delete orbitals are not zero
c      call ddatard(nft,x,nmob,idisk)
      call ddafile(luonemo,2,x,nmob,idisk)

c  read one electron fock matrix
      call ddafile(luonemo,2,xfock,nintone,idisk)
c      call ddatard(nft,xfock,nintone,idisk)
      call daclos(luonemo)
c  read one elctron kenetic intergrals
c      call ddatard(nft,x1e,nintone,idisk)
c      call cclose_molcas(nft)
c write one electron fock matrix into voint
      nidx=0
      do i=1,ng_sm
        nism  = nlsm_all(i)
        idx   = noidx(i)
        nsmint = nism*(nism+1)/2
        nc=0
        do lri=1,nism
          do lrj=1,lri
            nc=nc+1
            lrcii=map_orb_order(lri+idx)
            lrcij=map_orb_order(lrj+idx)
            if(lrcii.gt.lrcij) then
              lrt=lrcij
              lrcij=lrcii
              lrcii=lrt
            endif
            voint(lrcii,lrcij)=xfock(nc+nidx)
            !write(6,'(1x,i3,1x,i3,1x,f18.9)') lrcii,lrcij,xfock(nc+nidx)
          enddo
        enddo
        nidx=nidx+nsmint
      enddo

      call daname(lutwomo,fntwomo)
      call readtwoeint(lutwomo,nlsm_all,ng_sm,mul_tab,
     *                 map_orb_order,noidx)
      call daclos(lutwomo)
      write(6,*)
      deallocate(x)

      return
c...end of intrd_molcas
      end

      subroutine readtwoeint(nft,norb,ngsm,multab,maporb,noidx)
      implicit REAL*8 (a-h,o-z)
#include "ci_parameter.fh"
      parameter ( kbuf = ntrabuf )
      dimension norb(8),multab(8,8),maporb(max_orb),noidx(8)
      dimension itratoc(ntratoc)
      dimension buff(kbuf)

      idisk=0
      lenrd=ntratoc*lenintegral
      write(6,*) lenrd
c      call idatard(nft,itratoc,lenrd,idisk)
      call idafile(nft,2,itratoc,ntratoc,idisk)
c      write(6,"(10(1x,i8))") itratoc(1:10)
      write(6,2000)
2000  format(/7x,'symmetry',6x,' orbitals',8x,'integrals')
      do nsp=1,ngsm
        nop=norb(nsp)
        do nsq=1,nsp
          noq=norb(nsq)
          nspq=multab(nsp,nsq)
          do nsr=1,nsp
            nor=norb(nsr)
            nspqr=multab(nspq,nsr)
            nssm=nsr
            if(nsr.eq.nsp) nssm=nsq
            do nss=1,nssm
              if(nspqr.ne.nss) cycle
              nos=norb(nss)

              ityp=0
              if(nsr.eq.nss) then
                nbpq=(nop+nop**2)/2
                nbrs=(nos+nos**2)/2
                if(nsp.eq.nsr) then
c  (ii|ii) type 1 int
                  nintb=(nbpq+nbpq**2)/2
                  ityp=1
                else
c  (ii|jj) type 3 int
                  nintb=nbpq*nbrs
                  ityp=3
                endif
              else
                nbpq=nop*noq
                nbrs=nor*nos
                if(nsp.eq.nsr) then
c (ij|ij) type 2 int
                  nintb=(nbpq+nbpq**2)/2
                  ityp=2
                else
c (ij|kl) type 4 int
                  nintb=nbpq*nbrs
                  ityp=4
                endif
              endif

              if(nintb.eq.0) goto 10
              write(6,2100) nsp,nsq,nsr,nss,nop,noq,nor,nos,
     *                      nintb
2100          format(7x,4i2,1x,4i4,2x,3x,i9)

              idx=0
              iout=0
c              call ddatard(nft,buff,kbuf,idisk)
              call ddafile(nft,2,buff,kbuf,idisk)

              do li=1,nor
                ntj=nos
                if(nsr.eq.nss) ntj=li
                do lj=1,ntj
                  ntk=1
                  if(nsp.eq.nsr) ntk=li
                  do lk=ntk,nop
                    numin=1
                    if(nsp.eq.nsr.and.lk.eq.li)numin=lj
                    numax=noq
                    if(nsp.eq.nsq) numax=lk
                    do ll=numin,numax
                      iout=iout+1
                      if(iout.gt.kbuf) then
c                        call ddatard(nft,buff,kbuf,idisk)
                         call ddafile(nft,2,buff,kbuf,idisk)
                         iout=1
                      endif
                      idx=idx+1
                      val=buff(iout)
                      lri=maporb(lk+noidx(nsp))
                      lrj=maporb(ll+noidx(nsq))
                      lrk=maporb(li+noidx(nsr))
                      lrl=maporb(lj+noidx(nss))
!                      write(6,"(9(1x,i2),1x,f18.9)") li,lj,lk,ll,
!     *                                               lri,lrj,lrk,lrl,
!     *                                               iout,val
                      call intrw_mol(lri,lrj,lrk,lrl,val)
                    enddo
                  enddo
                enddo
              enddo

10          continue
            enddo
          enddo
        enddo
      enddo

      return
c....end of readtwoeint
      end

      subroutine readtraonehead(nft,ecor,idisk)
#include "ci_parameter.fh"
      implicit REAL*8 (a-h,o-z)
#include "maxbfn.fh"
      parameter (maxmolcasorb=maxbfn)
      dimension ncone(64),nbas(8),norb(8),nfro(8),
     *          ndel(8)
      dimension indx_idisk(64)
      character bsbl(2*4*maxmolcasorb)*1

      indx_idisk=0
      idisk=0
      call idafile(nft,2,ncone,64,idisk)
      call ddafile(nft,2,ecor,1,idisk)
      call idafile(nft,2,nsym,1,idisk)
      call idafile(nft,2,nbas,8,idisk)
      call idafile(nft,2,norb,8,idisk)
      call idafile(nft,2,nfro,8,idisk)
      call idafile(nft,2,ndel,8,idisk)
      lenrd=2*4*maxmolcasorb
      call cdafile(nft,2,bsbl,lenrd,idisk)

c#ifdef debug
c      write(6,"(a4,1x,8(2x,i8))") "ncon",ncone(1:8)
c      write(6,*) "idisk : ", idisk
c      write(6,"(a4,1x,f18.9)") "ecor",ecor
c      write(6,"(a4,1x,i8)") "nsym",nsym
c      write(6,"(a4,1x,8(2x,i8))") "nbas",nbas(1:8)
c      write(6,"(a4,1x,8(2x,i8))") "norb",norb(1:8)
c      write(6,"(a4,1x,8(2x,i8))") "nfro",nfro(1:8)
c      write(6,"(a4,1x,8(2x,i8))") "ndel",ndel(1:8)
c#endif

      return
c...end of readtraonehead
      end


#endif


c      subroutine ddatard(nft,xbuff,lend,idisk)
c      implicit REAL*8 (a-h,o-z)

c      dimension ncone(64),nbas(mxsym),norb(mxsym),nfro(mxsym),
c     *          ndel(mxsym)
c      dimension xbuff(lend)
c
c      lenrd=lend*lendbl
c      call idatard(nft,xbuff,lenrd,idisk)
c
c      return
cc...end of readtraonehead
c      end


c*************************************
c      subroutine idatard(nft,ibuf,lbuf,idisk)
cc   write by suo bing, read molcas file
cc   imul = 0, the file readed is not a multifile
cc   imul = 1, the file readed is a multifile
cc
c      implicit REAL*8 (a-h,o-z)
c      parameter (min_block_length=512)
c      character ibuf(lbuf)*1
c
c      noff=idisk*min_block_length
c      call clseek(nft,noff,nr)
c      if(nr.ne.noff) then
c        write(6,*) "error seek file : ",idisk
c        stop 888
c      endif
c      call cread_molcas(nft,ibuf,lbuf,nr)
c      if(nr.ne.lbuf) then
c        write(6,*) "error read file ",nr
c        write(6,*) "nft : ",nft
c        write(6,*) "lbuf  : ",lbuf
c        write(6,*) "idisk : ",idisk
c        stop 888
c      else
c        imo=mod(lbuf,min_block_length)
c        if(imo.eq.0) then
c          ioff=lbuf/min_block_length
c        else
c          ioff=1+lbuf/min_block_length
c        endif
c        idisk=idisk+ioff
cc        write(6,*) "ioff ",ioff
c      endif

c      return
cc...end of idatard
c      end

      subroutine intrw_mol(ik,jk,kk,lk,val)
#include "drt_h.fh"
#include "intsort_h.fh"
      dimension ind(4)
      list=0
      i=ik
      j=jk
      k=kk
      l=lk
      lri=min(i,j)
      lrj=max(i,j)
      lrk=min(k,l)
      lrl=max(k,l)
      if(lri.gt.lrk) then
        lrn=lrk
        lrk=lri
        lri=lrn
        lrn=lrl
        lrl=lrj
        lrj=lrn
      endif

      ind(1)=lri
      ind(2)=lrj
      ind(3)=lrk
      ind(4)=lrl
      do i=1,4
        do j=i+1,4
          if(ind(i).gt.ind(j)) then
            nt=ind(i)
            ind(i)=ind(j)
            ind(j)=nt
          endif
        enddo
      enddo

!      write(6,*) ind(1),ind(2),ind(3),ind(4)
!      goto 10
      if(ind(1).eq.ind(4)) then   !(iiii)
        vdint(lri,lri)=val
        goto 10
      endif

      if(lri.eq.lrj.and.lrk.eq.lrl) then   !(iikk)
        vdint(lrk,lri)=val
        goto 10
      endif

      if(lri.eq.lrk.and.lrj.eq.lrl) then   !(ijij)
!              if(lri.eq.1.and.lrj.eq.2) stop 888
        voint(lrj,lri)=val
        goto 10
      endif

      if(lri.ne.lrl.and.lrj.eq.lrk) then !(ijjl)
        list=list3_all(lri,lrl,lrj)
        vector1(list)=val
        if(lrj.eq.lrl) then
          vector1(list+1)=val
        endif
        if(lri.eq.lrj) then
          vector1(list+1)=val
        endif
        goto 10
      endif

      if(lri.ne.lrj.and.lrk.eq.lrl) then  !(ijkk)
        list=list3_all(lri,lrj,lrk)
        vector1(list+1)=val
        if(lrj.eq.lrk) then
          vector1(list)=val
        endif
        goto 10
      endif

      if(lri.eq.lrj.and.lrk.ne.lrl) then  !(iikl)
        list=list3_all(lrk,lrl,lri)
        vector1(list+1)=val
        goto 10
      endif

      if(lri.eq.lrk.and.lrj.ne.lrl) then      !(ijil)
        if(lrj.lt.lrl) then
          list=list3_all(lrj,lrl,lri)
        else
          list=list3_all(lrl,lrj,lri)
        endif
        vector1(list)=val
        goto 10
      endif

      if(lri.ne.lrk.and.lrj.eq.lrl) then   !(ijkj)
        list=list3_all(lri,lrk,lrj)
        vector1(list)=val
        if(lrj.eq.lrk) then
          vector1(list+1)=val
        endif
        goto 10
      endif

      if(lri.ne.lrj.and.lri.ne.lrk.and.
     *   lri.ne.lrl.and.lrj.ne.lrk.and.lrj.ne.lrl.and.lrk.ne.lrl) then
        list=list4_all(ind(1),ind(2),ind(3),ind(4))       !(ijkl)
        if(lrj.gt.lrk.and.lrj.gt.lrl) then
          vector1(list+2)=val
        endif
        if(lrj.gt.lrk.and.lrj.lt.lrl) then
          vector1(list)=val
        endif
        if(lrj.lt.lrk.and.lrj.lt.lrl) then
          vector1(list+1)=val
        endif

        goto 10
      endif

!10    write(6,*) "list=",list
10     return
      end

      subroutine int_index(numb)
#include "drt_h.fh"
#include "intsort_h.fh"
      dimension msob(8)

      msob=0
      do 12 la=norb_all,1,-1
        lra=norb_all-la+1
        ms = lsm(lra)
        msob(ms)=msob(ms)+1
        ncibl_all(la) = msob(ms)
12    continue

!      write(6,*)
c=======================================================================
c      write(6,'(1x,14i3)')(ncibl(i),i=1,14)

      numb = 1
      do 10  i = 1,norb_all-1
        do 10  j = i+1,norb_all
          lri=norb_number(i)
          lrj=norb_number(j)
          if(lsm(lri).ne.lsm(lrj)) goto 10
!          write(6,*) "lri=",lri,"lrj=",lrj
!          write(6,*) "lsm(lri)",lsm(lri),"lsm(lrj)",lsm(lrj)

          nij=i+ngw2(j)
c       write(6,*)'i,j,mij,nij   ',i,j,mij,nij
          loij_all(nij)=numb
!          do 20  k = 1,norb_all
!            vint_ci(numb)=vfutei(j,k,i,k)
!            vint_ci(numb+1)=vfutei(j,i,k,k)
c     write(10,'(2x,4i6,i8,3f16.8)')i,j,k,k, numb,
c    *        vint_ci(numb),vint_ci(numb+1)
            numb=numb+2*norb_all
!20        continue

10    continue
!      write(6,*) "number 3 index",numb
!      stop 777
c=======================================================================
c      la<lb<lc<ld
      do 30 ld = 1,norb_all-3
        do 30 lc = ld+1,norb_all-2
          lrd=norb_all-ld+1
          lrc=norb_all-lc+1
          msd  = lsm(lrd)
          msc  = lsm(lrc)
          mscd = mul_tab(msd,msc)
          do 30 lb = lc+1,norb_all-1
            lrb=norb_all-lb+1
            msb = lsm(lrb)
            msa = mul_tab(mscd,msb)

            njkl=ld+ngw2(lc)+ngw3(lb)
            loijk_all(njkl) = numb

            nolra=0
            do 40 la = norb_all,lb+1,-1
              lra=norb_all-la+1
              if(lsm(lra).ne.msa) cycle
!              nolra=nolra+1
!              list = loijk_all(njkl)+3*(nolra-1)
!              write(6,*) "ld,lc,lb,la",ld,lc,lb,la,list

c        write(6,'(2x,4i3,2i7)')  la,lb,lc,ld,list,numb

!             vint_ci(numb)=vfutei(la,lc,lb,ld)        !tmp stop
!             vint_ci(numb+1)=vfutei(la,lb,lc,ld)
!             vint_ci(numb+2)=vfutei(la,ld,lc,lb)
c     write(10,'(2x,4i6,i8,3f16.8)')la,lb,lc,ld, numb,
c    *        vint_ci(numb),vint_ci(numb+1),vint_ci(numb+2)
              numb=numb+3
40         continue
30    continue
!              stop 777
      return
c=======================================================================
      end

      REAL*8 function vfutei(ix,jx,kx,lx)
#include "drt_h.fh"
#include "intsort_h.fh"
      dimension ind(4)
      lri=min(ix,jx)
      lrj=max(ix,jx)
      lrk=min(kx,lx)
      lrl=max(kx,lx)
      if(lri.gt.lrk) then
        lrn=lrk
        lrk=lri
        lri=lrn
        lrn=lrl
        lrl=lrj
        lrj=lrn
      endif
      val=0.d0
!      write(6,*) ind(1),ind(2),ind(3),ind(4)
      if(lri.ne.lrl.and.lrj.eq.lrk) then !(ijjl)
        list=list3_all(lri,lrl,lrj)
        val=vector1(list)
        goto 10
      endif

      if(lri.ne.lrj.and.lrk.eq.lrl) then  !(ijkk)
        list=list3_all(lri,lrj,lrk)
        val=vector1(list+1)
        goto 10
      endif

      if(lri.eq.lrj.and.lrk.ne.lrl) then  !(iikl)
        list=list3_all(lrk,lrl,lri)
        val=vector1(list+1)
        goto 10
      endif

      if(lri.eq.lrk.and.lrj.ne.lrl) then      !(ijil)
        if(lrj.lt.lrl) then
          list=list3_all(lrj,lrl,lri)
        else
          list=list3_all(lrl,lrj,lri)
        endif
        val=vector1(list)
        goto 10
      endif

      if(lri.ne.lrk.and.lrj.eq.lrl) then   !(ijkj)
        list=list3_all(lri,lrk,lrj)
        val=vector1(list)
        goto 10
      endif

      if(lri.ne.lrj.and.lri.ne.lrk.and.
     *   lri.ne.lrl.and.lrj.ne.lrk.and.lrj.ne.lrl.and.lrk.ne.lrl) then
        ind(1)=lri
        ind(2)=lrj
        ind(3)=lrk
        ind(4)=lrl
        do i=1,4
          do j=i+1,4
            if(ind(i).gt.ind(j)) then
              nt=ind(i)
              ind(i)=ind(j)
              ind(j)=nt
            endif
          enddo
        enddo

        list=list4_all(ind(1),ind(2),ind(3),ind(4))       !(ijkl)
        if(lrj.gt.lrk.and.lrj.gt.lrl) then
          val=vector1(list+2)
        endif
        if(lrj.gt.lrk.and.lrj.lt.lrl) then
          val=vector1(list)
        endif
        if(lrj.lt.lrk.and.lrj.lt.lrl) then
          val=vector1(list+1)
        endif
        goto 10
      endif

10    vfutei=val
      end

      function list3_all(i,j,k)
#include "drt_h.fh"
#include "intsort_h.fh"
c            *****************
      nij   = i+ngw2(j)
      list3_all = loij_all(nij)+2*(k-1)
c            *****************
      return
      end
c
c***********************************************************************
      function list4_all(ld,lc,lb,la)
#include "drt_h.fh"
#include "intsort_h.fh"
c                    ***************
!      write(6,*) "ld,lc,lb,la",ld,lc,lb,la
        lra  = ncibl_all(la)
        njkl = ld+ngw2(lc)+ngw3(lb)
        list4_all= loijk_all(njkl)+3*(lra-1)
c                    ***************
      return
      end

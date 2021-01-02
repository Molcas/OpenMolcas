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
!-----------------------------------------------------------------------
!list of subroutines in this file
!      subroutine gugadrt_mole_inf()
!      subroutine gugadrt_paras_calculate()
!      subroutine arrange_orbital_molcas()
!      subroutine read_ref_state(nf)
!      subroutine gugadrt_active_drt()
!      subroutine gugadrt_ext_downwalk()
!      subroutine gugadrt_dbl_upwalk()
!      subroutine gugadrt_rst(id,nndd)
!      subroutine gugadrt_ref_gfs(nel,ndj,locu,nm)
!      subroutine gugadrt_rcas(id,indd)
!      subroutine gugadrt_check_rcas3(jk,ind,inb,ndj,locu)
!      subroutine gugadrt_dbl_downwalk()
!      subroutine gugadrt_ajphy(jp,in,jpihy)
!      subroutine gugadrt_njexcit(indjk,ljk,iextbit,nextbit,ivalid,jstep,kttmp,
!      subroutine packnod(ibuf,idx,ival,nin,nbit,lbuf)
!      subroutine packnod4(ibuf,idx1,idx2,ival,nin1,nbit1,
!      subroutine upacknod(ibuf,idx,ival,nin,nbit,lbuf)
!      subroutine redabkm(iabkm,labkm,nabcbit,iabcbit,
!      subroutine wrtabkm(iabkm,labkm,nabcbit,iabcbit,
!-----------------------------------------------------------------------
      subroutine gugainit()
! default value for performing ci calculation
#include "gendrt.fh"
#include "files_gugadrt.fh"
      parameter (maxmolcasorb=5000)
      dimension ncone(64),nbas(mxsym),norb(mxsym),nfro(mxsym),
     *          ndel(mxsym)
      dimension idx_idisk(64)
      dimension cmo(max_orb**2)
      character bsbl(2*4*maxmolcasorb)*1
      dimension dum(1),idum(1)

      fnonemo="TRAONE"
      fndrt="CIDRT"
      fncimo="CIMO"
      luonemo=30
      ludrt=31
      lucimo=32

      call daname(lucimo,fncimo)
      call daname(ludrt,fndrt)
c      call molcas_open(ludrt,fndrt)
! open file traone to read orbital informations
      call daname(luonemo,fnonemo)
c      call molcas_open(luonemo,fnonemo)
      idisk=0
      call idafile(luonemo,2,ncone,64,idisk)
      call ddafile(luonemo,2,dum,1,idisk)
c      ecor=dum(1)
      call idafile(luonemo,2,idum,1,idisk)
      nsym=idum(1)
      call idafile(luonemo,2,nbas,8,idisk)
      call idafile(luonemo,2,norb,8,idisk)
      call idafile(luonemo,2,nfro,8,idisk)
      call idafile(luonemo,2,ndel,8,idisk)
      lenrd=2*4*maxmolcasorb
      call cdafile(luonemo,2,bsbl,lenrd,idisk)
      nc=0
      do i=1,nsym
        nc=nc+nbas(i)**2
      enddo
      call ddafile(luonemo,2,cmo,nc,idisk)

      idx=0
      call idafile(lucimo,1,idx_idisk,64,idx)
      idx_idisk(2)=idx
      call cdafile(lucimo,1,bsbl,lenrd,idx)
      idx_idisk(3)=idx
      call ddafile(lucimo,1,cmo,nc,idx)
      idx_idisk(4)=idx
      idx=0
      call idafile(lucimo,1,idx_idisk,64,idx)

      call daclos(lucimo)
      call daclos(luonemo)
      nlsm_bas(1:8)=nbas(1:8)
      nlsmddel(1:8)=nfro(1:8)
      nlsmedel(1:8)=ndel(1:8)
c#ifdef debug
c      write(6,"(a4,1x,8(2x,i8))") "ncon",ncone(1:8)
c      write(6,*) "idisk : ", idisk
c      write(6,"(a4,1x,f18.9)") "ecor",ecor
c      write(6,"(a4,1x,i8)") "nsym",nsym
c      write(6,"(a4,1x,8(2x,i8))") "nbas",nbas(1:8)
c      write(6,"(a4,1x,8(2x,i8))") "norb",norb(1:8)
c      write(6,"(a4,1x,8(2x,i8))") "nfro",nfro(1:8)
c      write(6,"(a4,1x,8(2x,i8))") "ndel",ndel(1:8)

      ng_sm=nsym
      nlsm_all(1:8)=norb(1:8)
      mroot=1
      cm_cri=0.03

      return
c...end of subroutine gugadefault
      end

      subroutine gugadrt_mole_inf()
#include "gendrt.fh"
#include "files_gugadrt.fh"
#include "Sysdrt.fh"
!#ifndef _I8_
!      parameter (max_ref=64)
!#else
!      parameter (max_ref=128)
!#endif
      common /mcorb/ lsmorb(max_orb),noidx(8)
      common /refstate/ iref_occ(max_innorb,max_ref)
      dimension lsmtmp(maxgdm)
      logical log_debug
      dimension  itmpstr(72)
! copy from molcas, bsuo, jun. 30, 2009
      parameter ( ncmd=18 )
      parameter ( mxtit=10 )
      character*4 command,cmd(ncmd)
      character*72  line
      character*132 modline
      Data Cmd /'TITL','ELEC','SPIN','SYMM','ACTI',
     &          'PRIN','REFE','FIRS','INAC','CIAL',
     &          'VALE','INTE','NOCO','ONEO','EXTR',
     &          'NONI','NACT','END '/
*

      n_electron=0
      nactel=0
      spin=0.d0
      ns_sm=1
      call rdnlst(5,"GUGADRT")
      ntit=0
10    read(5,'(a)',end=991) line
      command=line(1:8)
      call upcase(command)
      if ( command(1:1).eq.'*' ) goto 10
      if (command.eq.' ') goto 10
      jcmd=0
      do icmd=1,ncmd
         if ( command.eq.cmd(icmd) ) jcmd=icmd
      end do
20    goto ( 100, 200, 300, 400, 500, 600, 700 ,800, 900,1000,
     &      1100,1200,1300,1400,1500,1600,1700,1800           ) jcmd
      write (6,*) 'input: illegal keyword'
      write (6,'(a,a)') 'command=',command
      call abend()
*
*---  process title    command ----------------------------------------*
 100  continue
      read(5,'(a)',end=991) line
      command=line(1:8)
      call upcase(command)
      if ( command(1:1).eq.'*' ) goto 100
      jcmd=0
      do icmd=1,ncmd
         if ( command.eq.cmd(icmd) ) jcmd=icmd
      end do
      if ( jcmd.ne.0 ) goto 20
      ntit=ntit+1
      goto 100
*
*---  process electron command ----------------------------------------*
 200  continue
      read(5,'(a)',end=991) line
      if ( line(1:1).eq.'*' ) goto 200
      read(line,*,err=992) n_electron
      goto 10
*
*---  process spin     command ----------------------------------------*
 300  continue
      read(5,'(a)',end=991) line
      if ( line(1:1).eq.'*' ) goto 300
      read(line,*,err=992) ispin
      spin=(ispin-1)/2.d0
      goto 10
*
*---  process symmetry command ----------------------------------------*
 400  continue
!      write (6,*)'input_guga: keyword symmetry is obsolete and ignored!
      read(5,'(a)',end=991) line
      if ( line(1:1).eq.'*' ) goto 400
      read(line,*,err=992) ns_sm
      goto 10
*
*---  process active   command ----------------------------------------*
 500  continue
      read(5,'(a)',end=991) line
      if ( line(1:1).eq.'*' ) goto 500
      modline=line//' 0 0 0 0 0 0 0 0'
      read(modline,*,err=992) (nlsm_act(i),i=1,8)
      goto 10
*
*---  process print    command ----------------------------------------*
 600  continue
      read(5,'(a)',end=991) line
      if ( line(1:1).eq.'*' ) goto 600
      read(line,*,err=992) iprint
      goto 10
*
*---  process referenc command ----------------------------------------*
 700  continue
      read(5,'(a)',end=991) line
      if ( line(1:1).eq.'*' ) goto 700
      read(line,*,err=992) n_ref,ln1
      if(n_ref.gt.max_ref) then
        write(6,*) " Warnning! Program could not deal with so much",
     *             " reference states!"
        write(6,*) " Maximum number of reference states is",max_ref
        write(6,*) " Set number of reference states into ",max_ref
        n_ref=max_ref
      endif
      if ( ln1.eq.0 ) then
        logic_mr=.true.
        goto 10
      endif
      do i=1,n_ref
        read(5,'(80i1)',end=991,err=992) iref_occ(1:ln1,i)
      enddo
      logic_mr=.true.
      goto 10
*
*---  process first    command ----------------------------------------*
 800  continue
      goto 10
*
*---  process inactive command ----------------------------------------*
 900  continue
      read(5,'(a)',end=991) line
      if ( line(1:1).eq.'*' ) goto 900
      modline=line//' 0 0 0 0 0 0 0 0'
      read(modline,*,err=992) (nlsm_dbl(i),i=1,8)
      goto 10
*
*---  process ciall    command ----------------------------------------*
1000  continue
!      read(5,'(a)',end=991) line
!      if ( line(1:1).eq.'*' ) goto 1000
!      read(line,*,err=992) ns_sm
      logic_mrelcas=.true.
      goto 10
*
*---  process valence  command ----------------------------------------*
1100  continue
      read(5,'(a)',end=991) line
      if ( line(1:1).eq.'*' ) goto 1100
      modline=line//' 0 0 0 0 0 0 0 0'
      read(modline,*,err=992) (nlsm_ext(i),i=1,8)
      goto 10
*
*---  process interact command ----------------------------------------*
1200  continue
      goto 10
*
*---  process nocorr   command ----------------------------------------*
1300  continue
      read(5,'(a)',end=991) line
      if ( line(1:1).eq.'*' ) goto 1300
      modline=line//' 0 0 0 0 0 0 0 0'
      read(modline,*,err=992)  !(ncor(i),i=1,8)
cbsuo, jun. 30, 2009 - neglect them
      goto 10
*
*---  process oneocc   command ----------------------------------------*
1400  continue
      read(5,'(a)',end=991) line
      if ( line(1:1).eq.'*' ) goto 1400
      modline=line//' 0 0 0 0 0 0 0 0'
      read(modline,*,err=992) ! (ione(i),i=1,8)
cbsuo, jun. 30, 2009 - neglect them
      goto 10
*
*---  process extract  command ----------------------------------------*
1500  write (6,*) 'input: extract option is redundant and is ignored!'
      goto 10
*
*---  process non-interact command -------------------------------------
1600  continue
      write (6,*) 'input: non-interact option is redundant and is',
     *            ' ignored!'
      goto 10
*
*---  process nactel       command -------------------------------------
1700  continue
      read(5,'(a)',end=991) line
      if ( line(1:1).eq.'*' ) goto 200
      read(line,*,err=992) nactel
      goto 10
*
*---  the end of the input is reached, print the title ----------------*
1800  continue

c frozen orbital(dbl, ext) have been delete in mo trans step, so we negl
c here.
      nlsm_frz(1:ng_sm)=0

      norb_frz=0
      norb_act=0
      norb_dz=0
      norb_all=0
      do i=1,ng_sm
        norb_frz=norb_frz+nlsm_frz(i)
        norb_dbl=norb_dbl+nlsm_dbl(i)
        norb_act=norb_act+nlsm_act(i)
        nlsm_inn(i)=nlsm_frz(i)+nlsm_dbl(i)+nlsm_act(i)
        nlsm_ext(i)=nlsm_all(i)-nlsm_inn(i)
        norb_all=norb_all+nlsm_all(i)
      enddo
      norb_inn=norb_frz+norb_dbl+norb_act
      norb_ext=norb_all-norb_inn
      norb_dz=norb_dbl+norb_frz

      nstart_act=norb_dz+1
      ngw1(1)=0
      ngw2(1)=0
      ngw2(2)=0
      ngw3(1)=0
      ngw3(2)=0
      ngw3(3)=0
      do i=1,norb_all
        ngw2(i+2) = ngw2(i+1)+i
        ngw3(i+3) = ngw3(i+2)+ngw2(i+2)
        ngw4(i+4) = ngw4(i+3)+ngw3(i+3)
      enddo
      nabc=norb_ext-2+ngw2(norb_ext-1)+ngw3(norb_ext)

      iorb=0
      do i=1,ng_sm
        do j=1,nlsm_frz(i)
          iorb=iorb+1
          lsm_inn(iorb)=i
        enddo
      enddo
      do i=1,ng_sm
        do j=1,nlsm_dbl(i)
          iorb=iorb+1
          lsm_inn(iorb)=i
        enddo
      enddo
      do i=1,ng_sm
        do j=1,nlsm_act(i)
          iorb=iorb+1
          lsm_inn(iorb)=i
        enddo
      enddo

      nact_sm=1
      do im=1,norb_inn
        if(lsm_inn(im).gt.nact_sm) nact_sm=lsm_inn(im)
      enddo

      norb_all_tmp=0
      do ngsm=1,ng_sm
        norb_all_tmp=norb_all_tmp+nlsm_all(ngsm)
      enddo
      if(norb_all_tmp.ne.norb_all) then
        write(6,*)'  input num.of orbital err! check again!'
        call abend
!        stop 777
      endif

      do l=1,norb_inn
        lr=norb_all-l+1
        lsm(lr)=lsm_inn(l)
      enddo

      lr=0
      int_dd_offset(1:8,1:8)=0
      do im=1,ng_sm
        im_lr_sta=0
        do iml=1,ng_sm
          imr=mul_tab(im,iml)
          if ( imr .gt. iml ) cycle
          int_dd_offset(iml,imr)=im_lr_sta
          int_dd_offset(imr,iml)=im_lr_sta
          if(iml.eq.imr)im_lr_sta=im_lr_sta+
     :           nlsm_ext(iml)*(nlsm_ext(iml)-1)/2
          if(iml.ne.imr)im_lr_sta=im_lr_sta+
     :           nlsm_ext(iml)*nlsm_ext(imr)
        enddo
        do l=1,nlsm_ext(im)
          lr=lr+1
          lsm(lr)=im
        enddo
      enddo
      lsm_inn(norb_inn+1)=lsm(norb_ext)

      lsmtmp(1:8)=0
      if (logic_mr) then
        do i=1,n_ref
          itmpstr(1:norb_act)=iref_occ(1:norb_act,i)
          do j=1,norb_dz
             iref_occ(j,i)=2
          enddo
          iref_occ(norb_dz+1:norb_dz+norb_act,i)=itmpstr(1:norb_act)
        enddo
      endif

      do i=1,ng_sm   !norb_inn
         lsmtmp(i)=0
      enddo
      do i=1,norb_inn
         lsmtmp(lsm_inn(i))=lsmtmp(lsm_inn(i))+1
      enddo

      itmp=0
      do i=1,ng_sm
         ibsm_ext(i)=itmp+1
         itmp=itmp+nlsm_ext(i)
         iesm_ext(i)=itmp
      enddo
!============ chck input data block ====================
      if(norb_frz.gt.norb_dz) then
        write(6,*)' check input data: norb_frz  error'
        call abend
!        stop 777
      endif
! check number of electrons
      if(nactel.ne.0) then
        if(n_electron.ne.0) then
          ne_act=n_electron-2*norb_dz
          if(nactel.ne.ne_act) then
            write(6,*) "Input error, Error in checking ELECtron",
     &                 " and NACTel!"
            call abend
          endif
        else
          ne_act=nactel
          n_electron=2*norb_dz+nactel
        endif
      else
        if(n_electron.ne.0) then
          ne_act=n_electron-2*norb_dz
          nactel=ne_act
        else
          ne_act=0
          write(6,*) "Input error, you must input ELECtron or NACTel!"
          call abend
        endif
      endif
      nde=2*norb_dz
      if(nde.gt.n_electron) then
        write(6,"(1x,42a)") 'check input date: number of elctrons error'
        write(6,"(1x,a20,1x,i4)") 'number of electrons ',n_electron
        write(6,"(1x,a36,1x,i4)") "number of doubly occupied electrons "
     *          ,nde
        call abend
!       stop 777
      endif

      norb1=norb_inn+norb_ext
      norb2=0
      do i=1,ng_sm
        norb2=norb2+nlsm_all(i)
      enddo
      if(norb1.ne.norb2) then
        write(6,*)' check input data: orb_number error'
        write(6,*) norb1,norb2
        call abend
!        stop 777
      endif

      if ( logic_mr) then
        write(6,*) " Refrence states"
!        ne_act=n_electron-2*norb_dz
        do i=1,n_ref
          write(6,"(80i1)") iref_occ(1:norb_inn,i)
          ms_ref=1
          neact=0
          do j=norb_dz+1,norb_inn
            if(iref_occ(j,i).eq.1) ms_ref=mul_tab(ms_ref,lsm_inn(j))
            if(iref_occ(j,i).eq.1) neact=neact+1
            if(iref_occ(j,i).eq.2) neact=neact+2
          enddo
          if(ms_ref.ne.ns_sm.or.neact.ne.ne_act) then
            write(6,*)'  input ref_conf  err! check again! ',i
            call abend
!           stop 777
          endif
        enddo
      endif
!============ block end =============================
c****************************************************
      log_debug=.false.
      if(log_debug) then
        write(6,1001)
        write(6,1002) norb_all,norb_frz,norb_dz,norb_inn,norb_ext
        write(6,1002) nlsm_all(1:ng_sm)
        write(6,1002) nlsm_frz(1:ng_sm)
        write(6,1002) nlsm_dbl(1:ng_sm)
        write(6,1002) nlsm_inn(1:ng_sm)
        write(6,*) "lsm inn"
        write(6,1002) lsm_inn(1:norb_inn)
        write(6,*) "lsm all",norb_all
        write(6,1002) (lsm(i),i=norb_all,1,-1)
      endif
c*****************************************************

! merge into molcas
c write date into cidrt for ci calculation
      noidx=0
      idisk=0
      call idafile(ludrt,1,noidx,2,idisk)
! group symmetry
      call idafile(ludrt,1,[ng_sm],1,idisk)
! state symmetry
      call idafile(ludrt,1,[ns_sm],1,idisk)
! number of roots to be cal
      call idafile(ludrt,1,[mroot],1,idisk)
! number of corelation electrons
      call idafile(ludrt,1,[n_electron],1,idisk)
! number of active electrons
      call idafile(ludrt,1,[nactel],1,idisk)
! spin symmetry of the state, 2s+1
      call idafile(ludrt,1,[ispin],1,idisk)
! dbl orb
      call idafile(ludrt,1,nlsm_dbl,8,idisk)
! act orb
      call idafile(ludrt,1,nlsm_act,8,idisk)
! all correlated orb
      call idafile(ludrt,1,nlsm_all,8,idisk)
! num. basis
      call idafile(ludrt,1,nlsm_bas,8,idisk)
! method to choose ref. state, 4 for cas, 2 for rst
      if(logic_mr) then
        call idafile(ludrt,1,[2],1,idisk)
        call idafile(ludrt,1,[n_ref],1,idisk)
! reference states
!       iref_occ(norb_dz+1:norb_dz+norb_act,i)=itmpstr(1:norb_act)
        do i=1,n_ref
           call idafile(ludrt,1,iref_occ(1,i),norb_inn,idisk)
        enddo
      endif
      if(logic_mrelcas) then
        call idafile(ludrt,1,[4],1,idisk)
      endif
! number of ref. states, if logic_mr = .true.
      noidx(2)=idisk
      idisk=0
      call idafile(ludrt,1,noidx,2,idisk)
      noidx=0
      do i=1,ng_sm
        noidx(i)=i
      enddo
      mul=nint(2*spin)+1
      write(6,*)'-----------------------------------------------'
      write(6,*)'    ci orbitals information'
      ndisk=0
      do i=1,8
        ndisk=ndisk+nlsm_bas(i)
      enddo
      write(6,*)'    num. of basis:          ',ndisk
      write(6,*)'    num. of orbitals:       ',norb_all
      write(6,*)'    num. of active-orbitals:',norb_act
      write(6,*)'    num. of electrons:      ',n_electron
      write(6,*)'    multiplicity:           ',mul
      write(6,*)'    symmetry:               ',ns_sm
      write(6,*)
      write(6,*)'    oribtials per-symmtry'
      write(6,1003) noidx(1:ng_sm)
      write(6,1004) nlsmddel(1:ng_sm)
      write(6,1005) nlsm_dbl(1:ng_sm)
      write(6,1006) nlsm_act(1:ng_sm)
      write(6,1007) nlsm_ext(1:ng_sm)
      write(6,1008) nlsmedel(1:ng_sm)
      write(6,*)'-----------------------------------------------'

      return
! following information will be drts
! notice we do not close file ludrt here, it will be close after drt
! information is written into this file
!      print*, "in subroutine gugadrt_mole_inf"
!      stop 1000
!---------------------------------------------------------------------
991   write (6,*) 'input: end of input file encountered'
      write (6,'(a,a)') 'last command: ',command
      call abend()
992   write (6,*) 'input: error while reading input!'
      write (6,'(a,a)') 'last command: ',command
      call abend()
1001  format(1x,"norb all group sm")
1002  format(8(1x,i3))
1003  format(16x,8(i4))
1004  format(8x,"frozen  ",8(i4))
1005  format(8x,"double  ",8(i4))
1006  format(8x,"active  ",8(i4))
1007  format(8x,"virtual ",8(i4))
1008  format(8x,"frz ext ",8(i4))
      end

      subroutine gugadrt_paras_calculate()
#include "gendrt.fh"
!      include "paraconstants_h.for"
!      include "intsort_h.for"
c      data inlptb_new/
!         1  2  3  4  5  6  7  8  9  10  11  12 13
c a^r=1
c     *    -1, 0, 0, 0, 0, 0,-7, 0,-9,-10,  0,-12, 0,
c a^l=2
c     *     0, 0, 0, 0, 0,-6, 0,-8, 0,  0,-11,  0, 0,
c b_r=3
c     *    4, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0,
c b_l=4
c     *    5, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0,
c b^r=5
c     *    0, 7, 8,10,11, 0, 0, 0, 0,  0,  0,  0, 0,
c b^l=6
c     *    0, 0, 9, 0,12, 0, 0, 0, 0,  0,  0,  0, 9,
c c^'=7
c     *    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 3,
c c^"=8
c     *    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,13,
c d_l^r=9
c     *    6, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0,
c d^r^r=10
c     *    0,-2, 0,-4, 0, 0, 0, 0, 0,  0,  0,  0, 0,
c d^r^l=11
c     *    0, 0,-3, 0,-5, 0, 0, 0, 0,  0,  0,  0, 0,
c c^'" =12
c     *    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 0/
!
!****************************************************************
!   ar      =1 (+a^r)     drr     =2 (+d^rr)   drl     =3 (+d^rl)
!   arbr    =4 (+d^rr)    arbl    =5 (+d^rl)   ard_l^r =6 (+a^l)
!   drrb^r  =7 (+a^r)     drlb^r  =8 (+a^l)    drlb^l  =9 (+a^r)
!   arbrb^r =10 (+a^r)    arblb^r =11 (+a^l)   arblb^l =12 (+a^r)
!   drl     =13 (*)
!****************************************************************

!     data inlptb_new
!    */ -1,  0,  4,  5,  0,  0,  1,  1,  6,  0,  0,  1,
!    *   0,  0,  0,  0,  7,  0,  2,  2,  0, -2,  0,  2,
!    *   0,  0,  0,  0,  8,  9,  3,  3,  0,  0, -3,  3,
!    *   0,  0,  0,  0, 10,  0,  4,  4,  0, -4,  0,  4,
!    *   0,  0,  0,  0, 11, 12,  5,  5,  0,  0, -5,  5,
!    *   0, -6,  0,  0,  0,  0,  6,  6,  0,  0,  0,  6,
!    *  -7,  0,  0,  0,  0,  0,  7,  7,  0,  0,  0,  7,
!    *   0, -8,  0,  0,  0,  0,  8,  8,  0,  0,  0,  8,
!    *  -9,  0,  0,  0,  0,  0,  9,  9,  0,  0,  0,  9,
!    * -10,  0,  0,  0,  0,  0, 10, 10,  0,  0,  0, 10,
!    *   0,-11,  0,  0,  0,  0, 11, 11,  0,  0,  0, 11,
!    * -12,  0,  0,  0,  0,  0, 12, 12,  0,  0,  0, 12,
!    *   0,  0,  0,  0,  0,  9,  3, 13,  0,  0,  0,  0/

!     data indord_plptype/0,0,0,1,1,3,3,1,1,2,2,2/   !severe_new_error_1
      ja_sys=int(n_electron*0.5d0-spin)-norb_dz
      jb_sys=int(spin+spin)
      jc_sys=norb_all-ja_sys-jb_sys

      if ( jb_sys .eq.0 ) jroute_sys=1
      if ( jb_sys .eq.1 ) jroute_sys=2
      if ( jb_sys .gt.1 ) jroute_sys=3

      return
      end

c****************************************************
      subroutine arrange_orbital_molcas()
c****************************************************
c  arrange orbital for ci calculation, in meld and
c  molcas program, the orbital are arranged as the symmetry
c  block. we transfer them to ci order
c -----
c map_order_orbital    ab ---> ci
c
#include "gendrt.fh"
!      include "intsort_h.for"
      common /mcorb/ lsmorb(max_orb),noidx(8)
      dimension lsmorbcount(ng_sm),map_tmp(max_orb)
      logical logi_norb_inn(norb_all)

      logi_norb_inn(1:norb_all)=.false.
      iorb=norb_all
      do la=1,norb_all
         norb_number(la)=iorb
         iorb=iorb-1
      enddo

      nim=0
      lsmorbcount(1)=nim
      do im=2,ng_sm
        nim=nim+nlsm_all(im-1)
        lsmorbcount(im)=nim
      enddo

      if ( logic_assign_actorb ) then
        do lr=1,norb_inn
          lr_scf=map_orb_order(lr)
          logi_norb_inn(lr_scf)=.true.
        enddo
        goto 200
      else
        do lr=1,norb_inn
          lsmr=lsm_inn(lr)
          lsmorbcount(lsmr)=lsmorbcount(lsmr)+1
          lr_scf=lsmorbcount(lsmr)
          map_orb_order(lr)=lr_scf
          logi_norb_inn(lr_scf)=.true.
        enddo
      endif
200   lr_scf0=norb_all
      la=norb_inn+1
      do ms=ng_sm,1,-1
        lr_scf0=lr_scf0-nlsm_all(ms)
        lr_scf=lr_scf0
        do lra=1,nlsm_all(ms)
          lr_scf=lr_scf+1
          if(logi_norb_inn(lr_scf)) cycle
          map_orb_order(la)=lr_scf
          la=la+1
        enddo
      enddo

      isum2=0
      isum3=0
      do i=1,ng_sm
        jp2(i)=isum2
        isum2=isum2+i
        jp3(i)=isum3
        isum3=isum3+isum2
      enddo

      iccount=1
      do lrd = 1,norb_inn
!        ipwt(lrd) = iccount
        iccount   = iccount+2
      enddo
      lsmorbcount =0
      do lrd=norb_dz,1,-1
        lsmid=lsm_inn(lrd)
        lsmorbcount(lsmid)=lsmorbcount(lsmid)+1
!        ipws(lrd)=(lsmorbcount(lsmid)-1)*3+1
      enddo

      map_tmp(1:norb_all)=map_orb_order(1:norb_all)
      do 40 i=1,norb_all
        do j=1,norb_all
          if(map_tmp(j).eq.i) then
            map_orb_order(i)=j
            goto 40
          endif
        enddo
40    continue

c      write(6,*) "map_order_orbit"
c      write(6,1001) map_orb_order(1:norb_all)
      return
c1001  format(20(1x,i3))
c...end of arrange_orbital_molcas
      end

      subroutine gugadrt_active_drt()
#include "gendrt.fh"
      common/casrst/ja(max_node),jb(max_node),jm(0:max_node)
     :    ,jj(4,0:max_node),kk(0:max_node),no(0:max_innorb)
     :    ,jv,jd(8),jt(8),js(8)
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
100   if(logic_mr)  call gugadrt_rst(ndd,indd)
      if(logic_mrelcas) call gugadrt_rcas(ndd,indd)   !npp=3

      nu_ae(1)=jv
      do im=1,ng_sm
        nu_ae(1+im)=jd(im)
        nu_ae(9+im)=jt(im)
        nu_ae(17+im)=js(im)
      enddo
      jds=1
      jde=mxnode
      jpe=no(norb_inn+1)

      iin(1:max_node)=0
      jp=jv
      iin(1:jpe)=0
      iin(0)=0
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
      enddo
      nci_dim=iseg_sta(25)+iseg_dim(25)
      iseg_sta(26)=nci_dim
      do jp=1,25
        if(iseg_downwei(jp).eq.0) cycle
        iseg_upwei(jp)=iseg_dim(jp)/iseg_downwei(jp)
      enddo
! to the end of dbl,act,ext parts
200   call gugadrt_dbl_downwalk()
c      write(6,*)'  end of drt,nci_dim= ',nci_dim
c      write(6,*)'number of cfss: ',nci_dim
      write(6,*)
      write(6,*)'-----------------------------------------------'
      write(6,*)'    csf information'
      write(6,*)'    num. of configurations:        ',nci_dim
      ndim=iseg_dim(1)
      write(6,*)'    num. of valence states:        ',ndim
      ndim=0
      do i=2,9
        ndim=ndim+iseg_dim(i)
      enddo
      write(6,"(5x,a32,1x,i12)")'num. of doublet couple singles: ',ndim
      ndim=0
      do i=10,17
        ndim=ndim+iseg_dim(i)
      enddo
      write(6,"(5x,a32,1x,i12)")'num. of triplet couple doubles: ',ndim
       ndim=0
      do i=18,25
        ndim=ndim+iseg_dim(i)
      enddo
      write(6,"(5x,a32,1x,i12)")'num. of singlet couple doubles: ',ndim

      write(6,*)'-----------------------------------------------'
      return
      end

      subroutine gugadrt_ext_downwalk()
#include "gendrt.fh"
      common/casrst/ja(max_node),jb(max_node),jm(0:max_node)
     :    ,jj(4,0:max_node),kk(0:max_node),no(0:max_innorb)
     :    ,jv,jd(8),jt(8),js(8)
      dimension iwmij(8)
      nu_ae(1)=1
      do im=1,ng_sm
        nu_ae(1+im)=1+im
        nu_ae(9+im)=9+im
        nu_ae(17+im)=17+im
      enddo

      iwmij=0
      iseg_downwei(nu_ae(1))=1
      do imi=1,ng_sm
        iseg_downwei(nu_ae(1+imi))=nlsm_ext(imi)
        do imj=imi,ng_sm
          imij=mul_tab(imi,imj)
          if(imij.ne.1) then
           iwmij(imij)=iwmij(imij)+nlsm_ext(imi)*nlsm_ext(imj)
             cycle
          endif
         iwmij(1)=iwmij(1)+nlsm_ext(imi)*(nlsm_ext(imi)-1)/2
        enddo
      enddo
      do im=1,ng_sm
        iseg_downwei(nu_ae(9+im))=iwmij(im)
        iseg_downwei(nu_ae(17+im))=iwmij(im)
      enddo
      iseg_downwei(nu_ae(18))=iseg_downwei(nu_ae(18))+norb_ext
      return
      end

      subroutine gugadrt_dbl_upwalk()
#include "gendrt.fh"
      if(norb_dbl.eq.1) then
c     v(1),d(2-9),s(18-25)           for s=0
c     v(1),d(2-9),s(18-25),d'(26-33)   for s<>0

        mxnode=17+ng_sm
        lri=norb_frz+1
        lsmi=lsm_inn(lri)
        lsmid=mul_tab(lsmi,ns_sm)
!for node_v
        nu_ad(1)=1
        jpad_upwei(1)=1
!for node_d
        nu_ad(1+lsmid)=1+lsmid
        jpad_upwei(1+lsmid)=1
!for node_s
        nu_ad(17+ns_sm)=17+ns_sm
        jpad_upwei(17+ns_sm)=1

        if(jroute_sys.eq.1) return
        mxnode=25+ng_sm
!for node_d'
        nu_ad(25+lsmid)=25+lsmid
        jpad_upwei(25+lsmid)=1
        return
      endif
      nu_ad=0
      jpad_upwei=0

      nu_ad(1)=1
      jpad_upwei(1)=1
      if(norb_dbl.eq.0) then
        mxnode=1
        return
      endif
      do lri=norb_frz+1,norb_dz
        lsmi=lsm_inn(lri)
        lsmid=mul_tab(lsmi,ns_sm)
        no_d=lsmid+1
        jpad_upwei(no_d)=jpad_upwei(no_d)+1
        do lrj=lri+1,norb_dz
          lsmj=lsm_inn(lrj)
           lsmij=mul_tab(lsmi,lsmj)
          lsmit=mul_tab(lsmij,ns_sm)
          no_t =lsmit+9
          jpad_upwei(no_t)=jpad_upwei(no_t)+1
        enddo
      enddo
c     v(1),d(2-9),t(10-17),s(18-25),d'(26-33),t'(34-41)
      select case (jroute_sys)
        case(1)
          goto 100
        case(2)
          goto 200
        case(3)
          goto 300
      end select
100     mxnode=25                     !v,d,t,s
        jpad_upwei(18:25)=jpad_upwei(10:17)
        jpad_upwei(17+ns_sm)=jpad_upwei(17+ns_sm)+norb_dbl
        goto 500
200     mxnode=25+8
        jpad_upwei(18:25)=jpad_upwei(10:17)+jpad_upwei(10:17)
        jpad_upwei(17+ns_sm)=jpad_upwei(17+ns_sm)+norb_dbl
        jpad_upwei(26:33)=jpad_upwei(2:9)
        goto 500
300     mxnode=25+8+8
        jpad_upwei(18:25)=jpad_upwei(10:17)+jpad_upwei(10:17)
        jpad_upwei(17+ns_sm)=jpad_upwei(17+ns_sm)+norb_dbl
        jpad_upwei(26:33)=jpad_upwei(2:9)
        jpad_upwei(34:41)=jpad_upwei(10:17)

500   do node=2,mxnode
        iw=jpad_upwei(node)
        if(iw.eq.0) cycle
        nu_ad(node)=node
      enddo
      return
      end

c*******************************************
      subroutine gugadrt_rst(id,nndd)
c*******************************************
c 10 may 2007 - revised by wyb
c
#include "gendrt.fh"
#include "Sysdrt.fh"
!#ifndef _I8_
!      parameter (iintbit=32,n32int=4,n16int=2)
!#else
!      parameter (iintbit=64,n32int=2,n16int=1)
!#endif
      common/casrst/ja(max_node),jb(max_node),jm(0:max_node)
     :    ,jj(4,0:max_node),kk(0:max_node),no(0:max_innorb)
     :    ,jv,jd(8),jt(8),js(8)
      common/ref/ndj,ndjgrop,ndjmod
      integer, pointer :: jabkm(:,:)
      integer, pointer :: ind(:,:)
      integer, pointer :: idjj(:,:)
      integer, pointer :: iwy(:,:)
      integer, pointer :: itm(:)
      dimension noh(max_innorb)
      dimension iextjj(n32int),iextii(n32int),jjabkm(1:n16int),
     *          jkabkm(1:n16int),iiabkm(1:n16int)
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

c 8 bits
      iextbit=2
      nextbit=iintbit/iextbit
c 16 bits
      iabcbit=16
      nabcbit=iintbit/iabcbit

c v node
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
c      call packnod(idkk,jk,norb_dz,nabcbit,iabcbit,mxtnode)

c**********************************************************************
c
7     j=j+1
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
c =========================================================
c                                            ***********
c                                            *   d=0   *
c                                            ***********

      jatmp=jaj
      jbtmp=jbj
      jmtmp=jmj
      kktmp=kkj+1


      iextjj(1:nrefbit)=ind(1:nrefbit,j)
      call gugadrt_njexcit(iextjj,nrefbit,iextbit,nextbit,ivalid,
     *                      0,kttmp,k0)
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

      call wrtabkm(jkabkm,n16int,nabcbit,iabcbit,jatmp,jbtmp,
     *             jmtmp,kktmp)
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
601     jk=jk-1
        idjj(1,j)=ii
        goto 11
18    continue

      ind(1:nrefbit,jk)=iextjj(1:nrefbit)
      idjj(1,j)=jk
      jabkm(1:n16int,jk)=jkabkm(1:n16int)

      iiabkm(1:n16int)=jabkm(1:n16int,jk-1)
      call upacknod(iiabkm,4,kj1,nabcbit,iabcbit,n16int)
      if(kktmp.ne.kj1) no(k0)=jk-1

11    continue
c =========================================================
      if(jbj.eq.0) goto 2
      jk=jk+1
c                                            ***********
c                                            *   d=1   *
c                                            ***********
      jatmp=jaj
      jbtmp=jbj-1
      jmtmp=mul_tab(lsm_inn(k0+1),jmj)
      kktmp=kkj+1

      iextjj(1:nrefbit)=ind(1:nrefbit,j)
      call gugadrt_njexcit(iextjj,nrefbit,iextbit,nextbit,ivalid,
     *                     1,kttmp,k0)

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

      call wrtabkm(jkabkm,n16int,nabcbit,iabcbit,jatmp,jbtmp,
     *             jmtmp,kktmp)
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
602     jk=jk-1
        idjj(2,j)=ii
        goto 22
28    continue

      ind(1:nrefbit,jk)=iextjj(1:nrefbit)
      idjj(2,j)=jk
      jabkm(1:n16int,jk)=jkabkm(1:n16int)

      iiabkm(1:n16int)=jabkm(1:n16int,jk-1)
      call upacknod(iiabkm,4,kj1,nabcbit,iabcbit,n16int)
      if(kktmp.ne.kj1) no(k0)=jk-1

22    continue
c =========================================================
2     jac=norb_all-jaj-jbj
      if(jaj.le.0) then
         goto 8
      else
         goto 5
      end if
5     if(jac.eq.0) goto 3
      jk=jk+1
c                                     *************
c                                     *    d=2    *
c                                     *************
      jatmp=jaj-1
      jbtmp=jbj+1
      jmtmp=mul_tab(lsm_inn(k0+1),jmj)
      kktmp=kkj+1

      iextjj(1:nrefbit)=ind(1:nrefbit,j)
      call gugadrt_njexcit(iextjj,nrefbit,iextbit,nextbit,ivalid,
     *                     2,kttmp,k0)
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

      call wrtabkm(jkabkm,n16int,nabcbit,iabcbit,jatmp,jbtmp,
     *             jmtmp,kktmp)
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
603     jk=jk-1
        idjj(3,j)=ii
        goto 33
38    continue

      ind(1:nrefbit,jk)=iextjj(1:nrefbit)
      idjj(3,j)=jk
      jabkm(1:n16int,jk)=jkabkm(1:n16int)

      iiabkm(1:n16int)=jabkm(1:n16int,jk-1)
      call upacknod(iiabkm,4,kj1,nabcbit,iabcbit,n16int)
      if(kktmp.ne.kj1) no(k0)=jk-1

33    continue
c =========================================================
c                                         ************
c                                         *    d=3   *
c                                         ************

3     jk=jk+1
      jatmp=jaj-1
      jbtmp=jbj
      jmtmp=jmj
      kktmp=kkj+1

      iextjj(1:nrefbit)=ind(1:nrefbit,j)
      call gugadrt_njexcit(iextjj,nrefbit,iextbit,nextbit,ivalid,
     *                     3,kttmp,k0)
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

      call wrtabkm(jkabkm,n16int,nabcbit,iabcbit,jatmp,jbtmp,
     *             jmtmp,kktmp)
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
604     jk=jk-1
        idjj(4,j)=ii
        goto 44
48    continue

      ind(1:nrefbit,jk)=iextjj(1:nrefbit)
      idjj(4,j)=jk
      jabkm(1:n16int,jk)=jkabkm(1:n16int)

      iiabkm(1:n16int)=jabkm(1:n16int,jk-1)
      call upacknod(iiabkm,4,kj1,nabcbit,iabcbit,n16int)
      if(kktmp.ne.kj1) no(k0)=jk-1

44    continue
c      if(k0.eq.7) goto 999
c===========================================================
8     if(k0.le.norb_inn-1) goto 7

c999   continue
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
c      write(6,"(10(1x,i5))") no(1:norb_all)
      write(6,*)
c  ************** external space  *************
      id=no(norb_inn)
c      write(6,508) 'befor,no=',(no(i),i=norb_dz,norb_inn)
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
304       if(j2.eq.0.or.iwy(1,j2).eq.0) goto 31
          iwy(2,jde)=iwy(1,j1)
31        if(j3.eq.0.or.iwy(1,j3).eq.0) goto 32
          iwy(3,jde)=iwy(1,j1)+iwy(1,j2)
32        if(j4.eq.0.or.iwy(1,j4).eq.0) goto 303
          iwy(4,jde)=iwy(1,jde)-iwy(1,j4)
303       if(jde.eq.1) goto 21
          do 302 jp=jps,jde-1
            if(iwy(1,jp).ne.iwy(1,jde)) goto 302
            jq1=idjj(1,jp)
            jq2=idjj(2,jp)
            jq3=idjj(3,jp)
            jq4=idjj(4,jp)
            if(j1.ne.jq1.or.j2.ne.jq2.or.j3.ne.jq3.or.j4.ne.jq4)
     *        goto 302
            iwy(1,jde)=0
            do jp0=no(l-2)+1,no(l-1)
              if(idjj(1,jp0).eq.jde) idjj(1,jp0)=jp
              if(idjj(2,jp0).eq.jde) idjj(2,jp0)=jp
              if(idjj(3,jp0).eq.jde) idjj(3,jp0)=jp
              if(idjj(4,jp0).eq.jde) idjj(4,jp0)=jp
            enddo
            goto 21
302       continue
21      continue
20    continue

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

c      open(10,file='rst.out')
      no(norb_dz)=0
      no(norb_dz+1)=mxnode
c      write(6,*) '   end of rst, drt ..........'
c      write(6,'(2x,2i10)')norb_dz+1,no(norb_dz+1)
      do 706 lr=norb_dz+1,norb_inn
        no(lr+1)=noh(lr)
706   continue
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
506   format('       j    k   a  b  t jm    j0   j1   j2   j3       y1',
     :       '       y2      y3         x   ind')
      nndd=no(norb_inn)
      do 541 j=1,id
          kk(j)=kk(j)+1
        if(iprint.eq.1) then
          write(6,510)j,kk(j),ja(j),jb(j),jm(j),
     :             jj(1,j),jj(2,j),jj(3,j),jj(4,j)
        endif
541   continue
      write(6,*)
      write(6,*) 'end of rst, drt ..........'

      deallocate(jabkm)
      deallocate(ind)
      deallocate(idjj)
      deallocate(iwy)
      deallocate(itm)

c
c      open(21,file="fort.drt",form="unformatted")
c      write(21) id
c      write(21) ja(1:id),jb(1:id),jm(1:id)
c      write(21) jj(1:4,0:id)
c      write(21) kk(0:id)
c      write(21) no(0:norb_inn+1)
c      write(21) jv,jd(1:8),jt(1:8),js(1:8)
c      close(21)

      call writedrt(id)
      return
c508   format(3x,a10,1x,i5,1x,16i8)
510   format(3x,10(1x,i7))
      end

      subroutine gugadrt_ref_gfs(nel,ndj,locu,nm)
#include "gendrt.fh"
#include "Sysdrt.fh"
      common/casrst/ja(max_node),jb(max_node),jm(0:max_node)
     :    ,jj(4,0:max_node),kk(0:max_node),no(0:max_innorb)
     :    ,jv,jd(8),jt(8),js(8)
      dimension lhsm(8),locu(8,max_ref),lscu(0:8,max_ref)
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

      subroutine gugadrt_rcas(id,indd)
#include "gendrt.fh"
#include "Sysdrt.fh"
#include "stdalloc.fh"
      common/casrst/ja(max_node),jb(max_node),jm(0:max_node)
     :    ,jj(4,0:max_node),kk(0:max_node),no(0:max_innorb)
     :    ,jv,jd(8),jt(8),js(8)
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
7     j=j+1
      if(j.le.mxnode) then
        if(nu_ad(j).eq.0) goto 7
      endif
!      if(j.le.mxnode.and.nu_ad(j).eq.0) goto 7
      k0=kk(j)
      if(jc(j).eq.0) goto 1
      jk=jk+1
c                                            ***********
c                                            *   d=0   *
c                                            ***********
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
18       continue
      jj(1,j)=jk
      kk(jk)=kk(j)+1
      if(kk(jk).ne.kk(jk-1)) then
      no(k0)=jk-1
      end if
11    continue
c      v=0
1     if(jb(j).eq.0) goto 2
      jk=jk+1
c                                            ***********
c                                            *   d=1   *
c                                            ***********
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
607      jk=jk-1
         jj(2,j)=ii
         goto 22
119   continue
      jj(2,j)=jk
      kk(jk)=kk(j)+1
      if(kk(jk).ne.kk(jk-1)) then
      no(k0)=jk-1
      end if
22    continue
c      v=0
2     if(ja(j).le.0) then
         goto 8
      else
         goto 5
      end if
5     if(jc(j).eq.0) goto 3
      jk=jk+1
c                                     *************
c                                     *    d=2    *
c                                     *************
      ja(jk)=ja(j)-1
      jb(jk)=jb(j)+1
      jc(jk)=jc(j)-1
c-----------------------------------------------------------------------
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
603      jk=jk-1
         jj(3,j)=ii
         goto 44
181   continue
      jj(3,j)=jk
      kk(jk)=kk(j)+1
      if(kk(jk).ne.kk(jk-1)) then
      no(k0)=jk-1
      end if
44    continue
c      v=0
c                                         ************
c                                         *    d=3   *
c                                         ************

3     jk=jk+1
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
606      continue
605      jk=jk-1
         jj(4,j)=ii
         go to 88
118   continue
      jj(4,j)=jk
      kk(jk)=kk(j)+1
      if(kk(jk).ne.kk(jk-1)) then
      no(k0)=jk-1
      end if
88    continue
c      v=0
      if(jk.gt.max_node) then
        write(6,*) '    the nomber of j exceeds max_node',max_node
        call abend
!       stop 777
      endif
8     if(kk(j).le.norb_inn-1)goto 7
c  ************** external space  *************

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
43    continue

      iwy(1,jv)=1
      do im=1,ng_sm
        if(jd(im).ne.0.and.iseg_downwei(1+im).ne.0) iwy(1,jd(im))=1
        if(jt(im).ne.0.and.iseg_downwei(9+im).ne.0) iwy(1,jt(im))=1
        if(js(im).ne.0.and.iseg_downwei(17+im).ne.0) iwy(1,js(im))=1
      enddo

      do 21 i=1,4
      iwy(i,0)=0
21    continue
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
304      if(j2.eq.0.or.iwy(1,j2).eq.0) goto 31
         iwy(2,jde)=iwy(1,j1)
31       if(j3.eq.0.or.iwy(1,j3).eq.0) goto 32
         iwy(3,jde)=iwy(1,j1)+iwy(1,j2)
32       if(j4.eq.0.or.iwy(1,j4).eq.0) goto 303
         iwy(4,jde)=iwy(1,jde)-iwy(1,j4)
303      if(jde.eq.1) goto 19
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
302      continue
19     continue
20    continue
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
200     continue
      noh(lr+1)=it
405   continue

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
704   continue
      do i=5,8
      ind(i,jp)=ind(i,jpe)
      enddo
206   continue

c      open(10,file='rcas.out')
      no(norb_dz)=0
      no(norb_dz+1)=mxnode

c      write(nf10,'(2x,2i10)')norb_dz+1,no(norb_dz+1)
      do 706 lr=norb_dz,norb_inn
        no(lr+1)=noh(lr)
        write(6,'(2x,2i10)')lr+1,no(lr+1)
706   continue

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
c        write(6,*)'    end of rcas , node=',id,'  dimension=',iysum
      write(6,*)
      indd=no(norb_inn)
!      iprint=1
      if(iprint.eq.1) then
        write(6,*) "guga drt"
        write(6,506)
      endif
506   format('       j    k   a  b  t jm    j0   j1   j2   j3       y1',
     :       '    y2      y3         x   ind')

      do 541 j=1,id
        kk(j)=kk(j)+1
        if(iprint.eq.1) then
           write(6,507)j,kk(j),ja(j),jb(j),jm(j),
     :             jj(1,j),jj(2,j),jj(3,j),jj(4,j),
     :             iwy(2,j),iwy(3,j),iwy(4,j),iwy(1,j),(ind(i,j),i=1,8)
        endif
541   continue
      write(6,*) 'end of rcas, drt ..........'
      write(6,*)

c      open(21,file="fort.drt",form="unformatted")
c      write(21) id
c      write(21) ja(1:id),jb(1:id),jm(1:id)
c      write(21) jj(1:4,0:id)
c      write(21) kk(0:id)
c      write(21) no(0:norb_inn+1)
c      write(21) jv,jd(1:8),jt(1:8),js(1:8)
c      close(21)

      call writedrt(id)
      call mma_deallocate(ind)
      call mma_deallocate(iwy)
507   format(3x,2i5,1x,3i3,1x,4i5,1x,4i10,1x,8i2)
      end

      subroutine gugadrt_check_rcas3(jk,ind,inb,ndj,locu)
#include "gendrt.fh"
      common/casrst/ja(max_node),jb(max_node),jm(0:max_node)
     :    ,jj(4,0:max_node),kk(0:max_node),no(0:max_innorb)
     :    ,jv,jd(8),jt(8),js(8)
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

c     juv,just(nost,nost),jud(nost)
c     |  \  1         |
c     | d,dd,s(i=i)   |
c     |    \ s,t,tt(i<j)|
c     |     \       1 2 |     deal with inner of dbl_space
c     |ss(i>j)\       |
c     |  2 1  \       |
      subroutine gugadrt_dbl_downwalk()
#include "gendrt.fh"
c     integer lsml(10,10)       !to del
      if(norb_dbl.ne.0) goto 200
c----------- norb_dbl=0 ------------------------------------------------
      do im=1,ng_sm
        nnd=iseg_sta(1+im)
        nnt=iseg_sta(9+im)
        nns=iseg_sta(17+im)
        do lri=norb_dz,norb_frz+1,-1
          ismi=lsm_inn(lri)
         if(ismi.ne.im) cycle
          jud(lri)=nnd
          nnd=nnd+iseg_downwei(1+im)
        enddo
        do lrj=norb_dz,norb_frz+1,-1
         ismj=lsm_inn(lrj)
          do lri=lrj,1,-1
           ismi=lsm_inn(lri)
           ismij=mul_tab(ismi,ismj)
            if(ismij.ne.im) cycle
            just(lri,lrj)=nns
            nns=nns+iseg_downwei(17+im)
            if(lri.eq.lrj) cycle
           just(lrj,lri)=nnt
            nnt=nnt+iseg_downwei(9+im)
           enddo
        enddo
      enddo
c----------- norb_dbl=0 ------------------------------------------------
c----------- norb_dbl<>0 -----------------------------------------------
200   continue
      do im=1,ng_sm
        nnd=0
        nns=0
        do lri=norb_frz+1,norb_dz
          ismi=mul_tab(lsm_inn(lri),ns_sm)
         if(ismi.ne.im) cycle
          jud(lri)=nnd
          nnd=nnd+1
        enddo
        do lri=norb_frz+1,norb_dz-1
          ismi=mul_tab(lsm_inn(lri),ns_sm)
          do lrj=lri+1,norb_dz      !tmp
           ismj=lsm_inn(lrj)
           ismij=mul_tab(ismi,ismj)
            if(ismij.ne.im) cycle
            just(lri,lrj)=nns
            nns=nns+1
         enddo
        enddo
        if(im.eq.ns_sm) then
          do lr0=norb_frz+1,norb_dz
            just(lr0,lr0)=nns
            nns=nns+1
          enddo
        endif
        do lri=norb_frz+1,norb_dz-1
          ismi=mul_tab(lsm_inn(lri),ns_sm)
          do lrj=lri+1,norb_dz      !tmp
           ismj=lsm_inn(lrj)
           ismij=mul_tab(ismi,ismj)
            if(ismij.ne.im) cycle
            just(lrj,lri)=nns
            nns=nns+1
         enddo
        enddo
      enddo

      return
      end

      function gugadrt_iwalk_ad(jdbl,jext,iwa,iwd)
#include "gendrt.fh"
      iwup=jpad_upwei(jdbl)
      isegdown=iseg_downwei(jext)
      gugadrt_iwalk_ad=(iwa*iwup+iwd)*isegdown+iw_sta(jdbl,jext)
      return
      end

      subroutine gugadrt_ajphy(jp,in,jpihy)
#include "gendrt.fh"
      common/casrst/ja(max_node),jb(max_node),jm(0:max_node)
     :    ,jj(4,0:max_node),kk(0:max_node),no(0:max_innorb)
     :    ,jv,jd(8),jt(8),js(8)
      common/sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),
     :     jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)
      dimension iin(0:max_node),jpihy(max_wei)

      iin(0)=0
      if(jp.eq.jpad) then
      in=1
      jpihy(1)=0
      return
      endif
      lr=kk(jp)
c     write(6,*)'  ajphy,jp,start,end',jp,no(nst-lr)+1,no(nst-lr+1)
      jpe=no(lr+1)
      iin(jpad:jpe)=0
      iin(jp)=1
      do 10 jpn=no(lr-1),jpad,-1
      iin(jpn)=iin(jj_sub(1,jpn))+iin(jj_sub(2,jpn))+iin(jj_sub(3,jpn))
     :        +iin(jj_sub(4,jpn))
10    continue
      in=iin(jpad)

      do 30 l=1,in
        jpihy(l)=0
        jy=l
        nn=jpad
        do 349 i=norb_dz+1,lr-1
          idr=0
          do 40 j=1,4
            if(jj(j,nn).eq.0) goto 40
            if(jy.gt.iin(jj_sub(j,nn))) goto 353
            idr=j
            goto 350
353         jy=jy-iin(jj_sub(j,nn))
40        continue
350       if(idr.ne.1)jpihy(l)=jpihy(l)+iy(idr,nn)
          nn=jj_sub(idr,nn)
349     continue
30    continue
      return
      end

      subroutine gugadrt_njexcit(indjk,ljk,iextbit,nextbit,ivalid,
     *                           jstep,kttmp,k0)
#include "gendrt.fh"
#include "Sysdrt.fh"
      common /refstate/ iref_occ(max_innorb,max_ref)
      common/ref/ndj,ndjgrop,ndjmod
      dimension indjk(ljk),itexcit(n_ref)

      kp=k0
      do idxref=1,n_ref
        call upacknod(indjk,idxref,ival,nextbit,iextbit,ljk)
        if(jstep.eq.1.or.jstep.eq.2) then
          if(iref_occ(kp+1,idxref).eq.0) ival=ival+1
        endif
        if(jstep.eq.3) then
          if(iref_occ(kp+1,idxref).eq.0) ival=ival+2
          if(iref_occ(kp+1,idxref).eq.1) ival=ival+1
        endif
        if(ival.gt.2) ival=3
        itexcit(idxref)=ival
      enddo
      inm=itexcit(1)
      do idxref=2,n_ref
        inm=min(inm,itexcit(idxref))
      enddo

      if(inm.gt.2) then
        ivalid=0
      else
        kttmp=inm
        ivalid=1
        if(jstep.ne.0) then
          do i=1,n_ref
            ival=itexcit(i)
            call packnod(indjk,i,ival,nextbit,iextbit,ljk)
          enddo
        endif
      endif

      return
      end

c**************************************************
      subroutine packnod(ibuf,idx,ival,nin,nbit,lbuf)
c**************************************************
c  pack integral ival into ibuf on bit mode
c  ibuf() integral buffer array
c  idx    index
c  nin    number of integrals in one integral
c  lbuf   length of ibuf
      implicit real*8 (a-h,o-z)
      dimension ibuf(lbuf)
      integer*4, parameter :: one4=1
      integer, parameter :: i4=kind(one4)

      inv=ival
      nimod=mod(idx,nin)
      if(nimod.eq.0) then
        ngrp=idx/nin
        nidbit=0
      else
        ngrp=idx/nin+1
        nidbit=(nin-nimod)*nbit
      endif

c      write(6,"(b64.64)") inv
#ifdef _AIX_
      call abend
#else
c IFG: changed to avoid compiler warnings, although it is
c      probably the compiler's fault
c     call mvbits(inv,0,nbit,ibuf(ngrp),nidbit)
      call mvbits(inv,0,int(nbit,i4),ibuf(ngrp),int(nidbit,i4))
#endif
c      write(6,"(b64.64)") ibuf(ngrp)

      return
c...end of packnod
      end

c**************************************************
      subroutine packnod4(ibuf,idx1,idx2,ival,nin1,nbit1,
     *                    lbuf1,lbuf2)
c**************************************************
c  pack integral ival into ibuf on bit mode
c  ibuf() integral buffer array
c  idx1   index of the ibuf(i,*)
c  idx2   index of the ibuf(*,i)
c  ival   the value will packed
c  nin1   number of integrals in one integral in ibuf(i,*)
c  nbit1  number of bits packed in ibuf(i,*)
c  nin2   number of integrals in one integral in ibuf(*,i)
c  nbit2  number of bits packed in ibuf(*,i)
c  lbuf   length of ibuf
      implicit real*8 (a-h,o-z)
      dimension ibuf(lbuf1,0:lbuf2)
      integer*4, parameter :: one4=1
      integer, parameter :: i4=kind(one4)

      inv=ival
      nimod1=mod(idx1,nin1)
      if(nimod1.eq.0) then
        ngrp1=idx1/nin1
        nidbit1=0
      else
        ngrp1=idx1/nin1+1
        nidbit1=(nin1-nimod1)*nbit1
      endif
#ifdef _AIX_
      call abend
#else

c IFG: changed to avoid compiler warnings, although it is
c      probably the compiler's fault
c     call mvbits(inv,0,nbit1,ibuf(ngrp1,idx2),nidbit1)
      call mvbits(inv,0,int(nbit1,i4),ibuf(ngrp1,idx2),int(nidbit1,i4))
#endif
      return
c...end of packnod4
      end

c**************************************************
      subroutine upacknod(ibuf,idx,ival,nin,nbit,lbuf)
c**************************************************
c  pack integral ival into ibuf on bit mode
c  ibuf() integral buffer array
c  idx    index
c  nin    number of integrals in one integral
c  lbuf   length of ibuf
      implicit real*8 (a-h,o-z)
      dimension ibuf(lbuf)
      integer*4, parameter :: one4=1
      integer, parameter :: i4=kind(one4)

      nimod=mod(idx,nin)
      if(nimod.eq.0) then
        ngrp=idx/nin
        nidbit=0
      else
        ngrp=idx/nin+1
        nidbit=(nin-nimod)*nbit
      endif

      isp=nidbit
c      inv=ibuf(ngrp)
c      ival=ibits(inv,isp,nbit-1)
      ival=0
#ifdef _AIX_
      call abend
#else

c IFG: changed to avoid compiler warnings, although it is
c      probably the compiler's fault
c     call mvbits(ibuf(ngrp),isp,nbit,ival,0)
      call mvbits(ibuf(ngrp),int(isp,i4),int(nbit,i4),ival,0)
#endif
c      write(6,*) isp,nbit-1
c      write(6,"(b64.64)") ival

      return
c...end of upacknod
      end


      subroutine redabkm(iabkm,labkm,nabcbit,iabcbit,
     *                   jatmp,jbtmp,jmtmp,kktmp)
c this subroutine unpack ja jb jm kt from compress arrays
      implicit real*8 (a-h,o-z)
      dimension iabkm(1:labkm)

      call upacknod(iabkm,1,jatmp,nabcbit,iabcbit,labkm)
      call upacknod(iabkm,2,jbtmp,nabcbit,iabcbit,labkm)
      call upacknod(iabkm,3,jmtmp,nabcbit,iabcbit,labkm)
      call upacknod(iabkm,4,kktmp,nabcbit,iabcbit,labkm)

      return
c...end of reabtm
      end

      subroutine wrtabkm(iabkm,labkm,nabcbit,iabcbit,
     *                   jatmp,jbtmp,jmtmp,kktmp)
c this subroutine unpack ja jb jm kt from compress arrays
      implicit real*8 (a-h,o-z)
      dimension iabkm(1:labkm)

      call packnod(iabkm,1,jatmp,nabcbit,iabcbit,labkm)
      call packnod(iabkm,2,jbtmp,nabcbit,iabcbit,labkm)
      call packnod(iabkm,3,jmtmp,nabcbit,iabcbit,labkm)
      call packnod(iabkm,4,kktmp,nabcbit,iabcbit,labkm)

      return
c...end of reabtm
      end

      subroutine writedrt(id)
#include "gendrt.fh"
#include "files_gugadrt.fh"
      common/casrst/ja(max_node),jb(max_node),jm(0:max_node)
     :    ,jj(4,0:max_node),kk(0:max_node),no(0:max_innorb)
     :    ,jv,jd(8),jt(8),js(8)
      dimension jbuf(4*(id+1)),idx(2)

      nc=1
      do i=0,id
        jbuf(nc:nc+3)=jj(1:4,i)
        nc=nc+4
      enddo
      idisk=0
!  number of nodes
      call idafile(ludrt,2,idx,2,idisk)
      idisk=idx(2)
      call idafile(ludrt,1,[id],1,idisk)
      call idafile(ludrt,1,ja,id,idisk)
      call idafile(ludrt,1,jb,id,idisk)
      call idafile(ludrt,1,jm,id,idisk)
      call idafile(ludrt,1,jbuf,4*(id+1),idisk)
      call idafile(ludrt,1,kk(0),1+id,idisk)
      call idafile(ludrt,1,no(0),norb_inn+2,idisk)
      call idafile(ludrt,1,[jv],1,idisk)
      call idafile(ludrt,1,jd,8,idisk)
      call idafile(ludrt,1,jt,8,idisk)
      call idafile(ludrt,1,js,8,idisk)

      return
      end

      subroutine gugadrt_gugafinalize()
! default value for performing ci calculation
#include "gendrt.fh"
#include "files_gugadrt.fh"

      call daclos(ludrt)
      end

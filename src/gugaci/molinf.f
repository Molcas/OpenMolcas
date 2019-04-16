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
c 28 dec 2007 -bsuo- the initial vectors for davidson diagnolazation are
c                    the vectors obtained by diagnol h0 space will be us
c                    initial vector for mrci calculation.
c 18 jul 2007 -bsuo- multi-root calculation is revised to perform calcul
c                    nci_dim*mroot>max_civector, the roots are calculated
c                    by one, old subroutine for davidson diagonalization
c                    deleted. the subroutines to init the initial vector
c                    revised.
c 15 mar 2007 -bsuo- new method for multi-root calculation finished
c 23 feb 2007 -bsuo- revised by suo bing for multi-root calculation
c
c this module consains some routines to calculate the segmentation
c value of the partial loops,information of the molecular orbital
c and the davidson diagonalization calculation
c
c*************************************************
      subroutine mole_inf()
c*************************************************
#include "drt_h.fh"
      common /thresh/ vthreen,vthrealp,vthreresid
      character*72  title(mxtit)
      parameter ( ncmd=9 )
      character*4 command,cmd(ncmd)
      character*72  line
      Data Cmd /'TITL','NRRO','MAXI','CPRO','PTHR',
     &          'CONV','PROR','REST',
#ifdef MOLPRO
#else
     &          'END '/
#endif
#ifdef _XIANEST_
     &         '$END'/
#endif
      logical logic_restart
#ifdef MOLPRO
#else
      call qenter("INPUT")
      call rdnlst(5,"GUGACI")
#endif

c init some program control logic variables
c trandional davidson diagnolization method is used
      logic_tdav=.true.
      logic_inivec_read=.false.
      logic_calpro=.false.
      cm_cri=0.05
      pror=0.0001
c set the default convergence threshhold
      vthreen=1.d-8
      vthrealp=1.d-6
      vthreresid=1.d-8
      mroot=1
      maxciiter=30
      ntit=0
10    read(5,'(a)',end=991) line
      command=line(1:8)
      call upcase(command)
      if ( command(1:1).eq.'*' ) goto 10
#ifdef _XIANEST_
      if(command(1:4).eq."$MRC") goto 10
#endif
      if (command.eq.' ') goto 10
      jcmd=0
      do icmd=1,ncmd
         if ( command.eq.cmd(icmd) ) jcmd=icmd
      end do
20    goto ( 100, 200, 300, 400, 500, 600,
     &       700, 800, 900
     &       ) jcmd
      write (6,*) 'input: illegal keyword'
      write (6,'(a,a)') 'command=',command
#ifdef MOLPRO
#else
      call qtrace
      call abend()
#endif
#ifdef _XIANEST_
      call qexit()
#endif
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
      if ( ntit.le.mxtit ) title(ntit)=line
      goto 100
*
*---  process nrroot command ----------------------------------------*
 200  continue
      read(5,'(a)',end=991) line
      if(line(1:1).eq."*") goto 200
      read(line,*,err=992) mroot
      goto 10

*
*---  process Maxiterations command ------------------------------------
 300  continue
      read(5,'(a)',end=991) line
      if(line(1:1).eq."*") goto 300
      read(line,*,err=992) maxciiter
      goto 10

*
*---  process Calculate property command -------------------------------
 400  continue
      logic_calpro=.true.
      goto 10

*
*---  process Thresh print command -------------------------------------
 500  continue
      read(5,'(a)',end=991) line
      if(line(1:1).eq."*") goto 500
      read(line,*,err=992) cm_cri
      goto 10
*
*---  process Convergence threshold command ----------------------------
 600  continue
      read(5,'(a)',end=991) line
      if(line(1:1).eq."*") goto 600
      read(line,*,err=992) vthreen,vthrealp,vthreresid
      goto 10
*
*---  process print orbital command ------------------------------------
 700  continue
      read(5,'(a)',end=991)  line
      if(line(1:1).eq."*") goto 700
      read(line,*,err=992) pror
      goto 10


*---  process restart command ------------------------------------------
 800  continue
      logic_restart=.true.
      goto 10

*--- End of GUGACI input ----------------------*
 900  continue

      call mole_inf_molcas()
      return
991   write (6,*) 'input: end of input file encountered'
      write (6,'(a,a)') 'last command: ',command
#ifdef _XIANEST_
      call qexit()
#endif
#ifdef MOLPRO
#else
      call qtrace
      call abend()
#endif
992   write (6,*) 'input: error while reading input!'
      write (6,'(a,a)') 'last command: ',command
#ifdef _XIANEST_
      call qexit()
#endif
#ifdef MOLPRO
#else
      call qtrace
      call abend()
#endif
      end

      subroutine mole_inf_molcas()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "files_gugaci.fh"
      common /thresh/ vthreen,vthrealp,vthreresid
      common /mcorb/ lsmorb(max_orb),noidx(8)
      dimension lsmtmp(maxgdm)
      dimension idum(1)

c      open(nf1,file="drt.inp")
c      read(nf1,*)
c      read(nf1,*) mroot,mth_eigen,kin
c      read(nf1,*) n_electron, spin, ng_sm, ns_sm, cm_cri
c      read(nf1,*) nlsm_all(1:ng_sm)
c      read(nf1,*) nlsm_frz(1:ng_sm)
c      read(nf1,*) nlsm_dbl(1:ng_sm)
c      read(nf1,*) nlsm_act(1:ng_sm)
c      read(nf1,*) log_thre
! merge into molcas
c write date into cidrt for ci calculation
      noidx=0
      idisk=0
      call idafile(ludrt,2,noidx,2,idisk)
! group symmetry
      call idafile(ludrt,2,idum,1,idisk)
      ng_sm=idum(1)
! state symmetry
      call idafile(ludrt,2,idum,1,idisk)
      ns_sm=idum(1)
! number of roots to be cal
      call idafile(ludrt,2,idum,1,idisk)
      nroot=idum(1)
! number of corelation electrons
      call idafile(ludrt,2,idum,1,idisk)
      n_electron=idum(1)
! number of active electrons
      call idafile(ludrt,2,idum,1,idisk)
      nactel=idum(1)
! spin symmetry of the state, 2s+1
      call idafile(ludrt,2,idum,1,idisk)
      ispin=idum(1)
! dbl orb
      call idafile(ludrt,2,nlsm_dbl,8,idisk)
! act orb
      call idafile(ludrt,2,nlsm_act,8,idisk)
! all correlated orb
      call idafile(ludrt,2,nlsm_all,8,idisk)
! num. basis
      call idafile(ludrt,2,nlsm_bas,8,idisk)
      spin=(ispin-1)/2.d0

      norb_frz=0
      norb_act=0
      norb_dz=0
      norb_all=0
      noidx=0
      idx=0
      ni=0
      do i=1,ng_sm
        norb_frz=norb_frz+nlsm_frz(i)
        norb_dbl=norb_dbl+nlsm_dbl(i)
        norb_act=norb_act+nlsm_act(i)
        nlsm_inn(i)=nlsm_frz(i)+nlsm_dbl(i)+nlsm_act(i)
        nlsm_ext(i)=nlsm_all(i)-nlsm_inn(i)
        norb_all=norb_all+nlsm_all(i)
        noidx(i)=idx
        idx=idx+nlsm_all(i)
        lsmorb(ni+1:idx)=i
        ni=idx
      enddo
      norb_inn=norb_frz+norb_dbl+norb_act
      norb_ext=norb_all-norb_inn
      norb_dz=norb_dbl+norb_frz

      nstart_act=norb_dz+1
      ngw1=0
      ngw2=0
      ngw3=0
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
#ifdef MOLPRO
#else
      call qtrace
      call abend()
#endif
#ifdef _XIANEST_
      call qexit()
#endif
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

      logic_mr=.false.
      logic_mrelcas=.false.
      logic_assign_actorb=.false.

      call idafile(ludrt,2,idum,1,idisk)
      imrcas_case=idum(1)
      if(imrcas_case.eq.2) then
        logic_mr=.true.
        call idafile(ludrt,2,idum,1,idisk)
        n_ref=idum(1)
        do i=1,n_ref
          call idafile(ludrt,2,iref_occ(1,i),norb_inn,idisk)
        enddo
      endif
      if(imrcas_case.eq.4) then
        logic_mrelcas=.true.
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

c****************************************************
#ifdef _DEBUG
      idebug=1
      if(idebug.eq.1) then
        write(6,1001)
        write(6,1002) norb_all,norb_frz,norb_dz,norb_inn,norb_ext
        write(6,1002) nlsm_all(1:ng_sm)
        write(6,1002) nlsm_frz(1:ng_sm)
        write(6,1002) nlsm_dbl(1:ng_sm)
        write(6,1002) nlsm_inn(1:ng_sm)
        write(6,1002) nlsm_ext(1:ng_sm)
        write(6,*) "lsm inn, imrcas_case",imrcas_case
        write(6,1002) lsm_inn(1:norb_inn)
        write(6,*) "lsm all",norb_all
        write(6,1002) (lsm(i),i=norb_all,1,-1)
      endif
1001  format(1x,"all frz dz inn ext")
1002  format(8(1x,i3))
#endif
c*****************************************************

      return
      end

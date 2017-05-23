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
* Copyright (C) 1996-2006, T. Thorsteinsson and D. L. Cooper           *
************************************************************************
      subroutine mkguess2_cvb(orbs,cvb,irdorbs,orbsao)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
#include "mo_cvb.fh"
      dimension orbs(norb,norb),cvb(nvb)
      dimension irdorbs(norb),orbsao(nbas_mo,norb)
      save thresh
      data thresh/1d-10/

      call izero(irdorbs,norb)
c  -- transfer from orbs if applicable --
c  (Newly assigned memory => orbs will be zero)
      do 100 iorb=1,norb
      if(dnrm2_(norb,orbs(1,iorb),1).gt.thresh)then
        irdorbs(iorb)=1
        call fmove(orbs(1,iorb),orbsao(1,iorb),norb)
      endif
100   continue
c  -- restore from previous optim --
      if(.not.up2date_cvb('RESTGS'))then
        if(up2date_cvb('WRITEGS'))then
          call rdi_cvb(ndetvb1,1,recn_tmp04,0)
          i1=mstacki_cvb(ndetvb1)
          i2=mstackr_cvb(ndetvb1)
          call mkrestgs_cvb(orbsao,irdorbs,cvb,
     >      w(lw(9)),iw(ll(11)),iw(ll(12)),iw(i1),w(i2))
          call mfreei_cvb(i1)
        endif
        call untouch_cvb('RESTGS')
      endif
c  -- read from file --
      if(.not.up2date_cvb('STRTGS'))then
        call setstrtvb_cvb(strtvb)
        if(tstfile_cvb(strtvb))call mkstrtgs_cvb(orbsao,irdorbs,cvb,
     >    strtvb,kbasiscvb)
        call untouch_cvb('STRTGS')
      endif
c  -- input --
      if(.not.up2date_cvb('INPGS'))then
        i1 = mstacki_cvb(norb)
        call rdioff_cvb(6,recinp,ioffs)
        call rdi_cvb(iw(i1),norb,recinp,ioffs)
        call rdioff_cvb(5,recinp,ioffs)
        do 200 iorb=1,norb
        if(iw(iorb+i1-1).eq.1)then
c MO basis ...
          irdorbs(iorb)=1
          call rdr_cvb(orbsao(1,iorb),norb,recinp,ioffs)
        elseif(iw(iorb+i1-1).eq.2)then
c AO basis ...
          irdorbs(iorb)=2
          call rdr_cvb(orbsao(1,iorb),nbas_mo,recinp,ioffs)
        endif
200     ioffs=ioffs+mxaobf
        call mfreei_cvb(i1)

        i1 = mstackr_cvb(nvbinp)
        call rdioff_cvb(7,recinp,ioffs)
        call rdrs_cvb(w(i1),nvbinp,recinp,ioffs)
        if(dnrm2_(nvbinp,w(i1),1).gt.thresh)then
          call rdioff_cvb(3,recinp,ioffs)
          call rdis_cvb(kbasiscvb,1,recinp,ioffs)
          call fmove(w(i1),w(lv(2)),nvbinp)
        endif
        call mfreer_cvb(i1)

        call untouch_cvb('INPGS')
      endif
c  -- semi-random --
c  Leading diagonal, random but positive orbital overlaps :
      dum=rand_cvb(.777d0)
      c=1d-1
      do 300 iorb=1,norb
      if(irdorbs(iorb).eq.0)then
        irdorbs(iorb)=1
        do 400 ii=1,norb
        orbsao(ii,iorb)=c*rand_cvb(zero)
400     if(ii.eq.iorb)orbsao(ii,iorb)=one
      else
c  Dummy calls to RAND to get consistent guesses :
        do 500 ii=1,norb
500     dum=rand_cvb(zero)
      endif
300   continue

c  Collect orbitals and transform AO -> MO :
      norb_ao=0
      do iorb=1,norb
      if(irdorbs(iorb).eq.1)then
        call fmove(orbsao(1,iorb),orbs(1,iorb),norb)
      elseif(irdorbs(iorb).eq.2)then
        norb_ao=norb_ao+1
        if(norb_ao.ne.iorb)call fmove(orbsao(1,iorb),
     >    orbsao(1,norb_ao),nbas_mo)
      endif
      enddo
      i1=mstackr_cvb(norb*norb_ao)
      call ao2mo_cvb(orbsao,w(i1),norb_ao)
      iorb_ao=0
      do iorb=1,norb
      if(irdorbs(iorb).eq.2)then
        iorb_ao=iorb_ao+1
        call fmove(w((iorb_ao-1)*norb+i1),orbs(1,iorb),norb)
      endif
      enddo
      call mfreer_cvb(i1)

      call nize_cvb(orbs,norb,dum,norb,0,0)

      if(abs(detm_cvb(orbs,norb)).lt.1d-8)then
        dum=rand_cvb(.777d0)
        c=1d-1
        do 600 iorb=1,norb
        do 600 ii=1,norb
600     orbs(ii,iorb)=orbs(ii,iorb)+c*(one-two*rand_cvb(zero))
        if(abs(detm_cvb(orbs,norb)).lt.1d-8)then
          if(ip(1).ge.0)
     >      write(6,'(2a)')' Starting orbital guess was near-singular',
     >        ' - using semi-random guess instead.'
          dum=rand_cvb(.777d0)
          c=1d-1
          do 700 iorb=1,norb
          do 700 ii=1,norb
          orbs(ii,iorb)=c*rand_cvb(zero)
700       if(ii.eq.iorb)orbs(ii,iorb)=one
        else
          if(ip(1).ge.0)
     >      write(6,'(2a)')' Starting orbital guess was near-singular',
     >        ' - scrambling orbital coefficients.'
        endif
        call nize_cvb(orbs,norb,dum,norb,0,0)
      endif

      call nize_cvb(orbs,norb,dum,norb,0,0)

c  Perfect-pairing spin function(s) :
      if(dnrm2_(nvb,cvb,1).lt.thresh)then
        kbasiscvb=kbasis
        call ppgs_cvb(cvb)
      endif

      cnrm=dnrm2_(nvb,cvb,1)
      if(cnrm.lt.thresh)then
        write(6,*)' Fatal error - starting structure coefficients ',
     >    'all zero !'
        call abend_cvb()
      endif

      if(kbasiscvb.ne.kbasis)then
        call mktrnspn_cvb()
        call untouch_cvb('TRNSPN')
      endif

      if(ploc)call rtransf_plc(orbs,cvb)

      if(ip(1).ge.2.and..not.endvar)then
        write(6,'(/,a)')' Wavefunction guess :'
        call report_cvb(orbs,norb)
        write(6,'(/,a)')' Structure coefficients :'
        write(6,'(a)')' ------------------------'
        call vecprint_cvb(cvb,nvb)
      endif

      return
      end

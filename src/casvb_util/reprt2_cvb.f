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
      subroutine reprt2_cvb(orbs,cvb,
     >  civec,civb,civbs,civbh,citmp,
     >  sstruc,sstruc2,
     >  orbinv,sorbs,owrk,gjorb,gjorb2,gjorb3,
     >  cvbstot,cvbsspn,cvbdet,dvbdet,evbdet,
     >  dmat,occ)
      implicit real*8 (a-h,o-z)
      logical make_sstruc
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "formats_cvb.fh"
      dimension orbs(norb,norb),cvb(nvb)
      dimension civec(ndet),civb(ndet)
      dimension civbs(ndet),civbh(ndet),citmp(ndet)
      dimension sstruc(nvb,nvb),sstruc2(nvb,nvb)
      dimension orbinv(norb,norb),sorbs(norb,norb),owrk(norb,norb)
      dimension gjorb(*),gjorb2(*),gjorb3(*)
      dimension cvbstot(nvb),cvbsspn(nvb)
      dimension cvbdet(ndetvb),dvbdet(ndetvb),evbdet(ndetvb)
      dimension dmat(norb,norb),occ(norb)

      if(ip(5).ge.1)then
        call report_cvb(orbs,norb)
        write(6,'(/,a)')' Structure coefficients :'
        write(6,'(a)')' ------------------------'
        call vecprint_cvb(cvb,nvb)
      endif

c First save CI vector
      call str2vbc_cvb(cvb,cvbdet)
      call vb2cic_cvb(cvbdet,civb)
      call gaussj_cvb(orbs,gjorb)
      call applyt_cvb(civb,gjorb)
      call proj_cvb(civb)
      call cinorm2_cvb(civb,cnrm)
      call ciscale_cvb(civb,one/cnrm)
      call symweight_cvb(civb,civb,wsym)

c Save VB wavefunction
      call setsavvb_cvb(savvb)
      if(valid_cvb(savvb))then
        if(ip(5).ge.1)then
          write(6,'(a)')' '
          call prtfid_cvb(' Saving VB wavefunction to ',savvb)
        endif
        call putguess_cvb(orbs,cvb,savvb)
      endif
      call putci_cvb(civb)
      call fzero(dmat,norb*norb)
      call onedens_cvb(civb,civb,dmat,.true.,0)
c  Before overwriting CIVB/CIVEC get SVB :
      if(lcalcsvb)then
        call cidot_cvb(civec,civb,svb)
        call untouch_cvb('SVB')
      endif
      if(lcalcevb)then
        call cicopy_cvb(civb,civbh)
        call applyh_cvb(civbh)
        call cidot_cvb(civb,civbh,evb)
        evb=evb+corenrg
        call untouch_cvb('EVB')
c  ESYM needed for variational calculation
        call symweight_cvb(civb,civbh,esym)
        do 100 i=1,nirrep
        if(wsym(i).gt.1d-10)then
          esym(i)=esym(i)/wsym(i)+corenrg
        else
          esym(i)=zero
        endif
100     continue
      endif
      if(lcalcsvb.and.ip(5).ge.1)then
        write(6,'(a)')' '
        write(6,formE)' Svb :      ',svb
        if((lcalcevb.and.ip(5).ge.1).or.
     >    (icrit.eq.2.and.ip(3).lt.0.and.ip(5).ge.1))
     >    write(6,formE)' Evb :      ',evb
      else
        if((lcalcevb.and.ip(5).ge.1).or.
     >    (icrit.eq.2.and.ip(3).lt.0.and.ip(5).ge.1))then
          write(6,'(a)')' '
          write(6,formE)' Evb :      ',evb
        endif
      endif
      if(variat.and.(.not.endvar))return

      call makegjorbs_cvb(orbs,gjorb,gjorb2,gjorb3)
      call makecivbs_cvb(civbs,orbs,gjorb,gjorb2,gjorb3,cvbdet)

      call ciscale_cvb(civbs,one/cnrm)
      call dscal_(nvb,one/cnrm,cvb,1)
      call dscal_(ndetvb,one/cnrm,cvbdet,1)
      call finalresult_cvb()
c -- Analysis
      if(.not.lcalccivbs)then
c  CIVBS has been evaluated previously (NB. based on unnormalized CVB) :
        call ciscale_cvb(civbs,one/cnrm)
        call ci2vbg_cvb(civbs,dvbdet)
      else
        call cicopy_cvb(civb,civbs)
        call transp_cvb(orbs,owrk,norb,norb)
        call gaussj_cvb(owrk,gjorb2)
        call applyt_cvb(civbs,gjorb2)
        call ci2vbg_cvb(civbs,dvbdet)
      endif
      call vb2strg_cvb(dvbdet,cvbstot)
      sum1=ddot_(nvb,cvb,1,cvbstot,1)
      call vb2strg_cvb(cvbdet,cvbsspn)
      sum2=ddot_(nvb,cvb,1,cvbsspn,1)
      if((ip(5).ge.1.and.ivbweights.lt.0)
     >  .or.(ivbweights.ge.0.and.mod(ivbweights,2).eq.1))then
        write(6,'(/,a)')' Chirgwin-Coulson weights of structures :'
        write(6,'(a)')' ----------------------------------------'
        write(6,formVBWnorm) ' VB spin+space (norm ',sum1,') :'
        do i=1,nvb
        dvbdet(i)=cvb(i)*cvbstot(i)/sum1
        enddo
        call vecprint_cvb(dvbdet,nvb)
        if(sum2.lt.1d3)then
          write(6,formVBWnorm)' VB spin only  (norm ',sum2,') :'
        else
          formVBWnorm(4:4)='e'
          write(6,formVBWnorm)' VB spin only  (norm ',sum2,') :'
          formVBWnorm(4:4)='f'
        endif
        do i=1,nvb
        dvbdet(i)=cvb(i)*cvbsspn(i)/sum2
        enddo
        call vecprint_cvb(dvbdet,nvb)
      endif
      make_sstruc=(ivbweights.gt.1.or.ishstruc.eq.1)
      if(make_sstruc)then
        if(.not.proj)then
          call mxattb_cvb(orbs,orbs,norb,norb,norb,sorbs)
          call gaussj_cvb(sorbs,gjorb3)
        else
          call gaussj_cvb(orbs,gjorb)
          call transp_cvb(orbs,owrk,norb,norb)
          call gaussj_cvb(owrk,gjorb2)
        endif
        do 7200 k=1,nvb
        call fzero(sstruc(1,k),nvb)
        sstruc(k,k)=1d0
        call str2vbf_cvb(sstruc(1,k),dvbdet)
        call vb2cif_cvb(dvbdet,civb)
        if(.not.proj)then
          call applyt_cvb(civb,gjorb3)
        else
          call applyt_cvb(civb,gjorb)
          call proj_cvb(civb)
          call applyt_cvb(civb,gjorb2)
        endif
        call ci2vbg_cvb(civb,dvbdet)
7200    call vb2strg_cvb(dvbdet,sstruc(1,k))
        do 7220 k=1,nvb
        do 7220 l=k+1,nvb
        sstruc(k,l)=.5d0*(sstruc(k,l)+sstruc(l,k))
7220    sstruc(l,k)=sstruc(k,l)
      endif

      if(ivbweights.gt.1)then
c  --  VB analysis Lowdin and Inverse  - begin  --
        if(mod(ivbweights,8).gt.3)then
          write(6,'(/,a)')
     >      ' Inverse-overlap weights of structures :'
          write(6,'(a)')
     >      ' ---------------------------------------'
          call fmove_cvb(sstruc,sstruc2,nvb*nvb)
          call mxinv_cvb(sstruc2,nvb)
c  Use DVBDET for weights :
          sum=zero
          do 7300 k=1,nvb
          dvbdet(k)=cvb(k)*cvb(k)/sstruc2(k,k)
7300      sum=sum+dvbdet(k)
          write(6,formVBWnorm)' VB spin+space (norm ',sum,') :'
          call dscal_(nvb,one/sum,dvbdet,1)
          call vecprint_cvb(dvbdet,nvb)
        endif

        if(mod(ivbweights,4).gt.1)then
          write(6,'(/,a)')
     >      ' Weights of Lowdin-orthogonalized structures :'
          write(6,'(a)')
     >      ' ---------------------------------------------'
          call fmove_cvb(sstruc,sstruc2,nvb*nvb)
c  Normalise overlap matrix before square root :
c  Use CVBDET for normalized structure coefficients :
          call fmove_cvb(cvb,cvbdet,nvb)
          do 7400 k=1,nvb
          fac=sqrt(sstruc2(k,k))
          cvbdet(k)=fac*cvbdet(k)
          fac=one/fac
          call dscal_(nvb,fac,sstruc2(1,k),1)
          call dscal_(nvb,fac,sstruc2(k,1),nvb)
7400      sstruc2(k,k)=one

          call mxsqrt_cvb(sstruc2,nvb,1)
c  Use DVBDET for weights :
          call mxatb_cvb(sstruc2,cvbdet,nvb,nvb,1,dvbdet)
          sum=zero
          do 7500 k=1,nvb
          dvbdet(k)=dvbdet(k)*dvbdet(k)
7500      sum=sum+dvbdet(k)
          write(6,formVBWnorm) ' VB spin+space (norm ',sum,') :'
          call dscal_(nvb,one/sum,dvbdet,1)
          call vecprint_cvb(dvbdet,nvb)
        endif
      endif
c  --  VB analysis Lowdin and Inverse  - end    --

c  Spin correlation analysis
c  Transform spin -> det to cater for non-orthogonal spin basis
      if(sij.and.sc)then
        call str2vbc_cvb(cvb,cvbdet)
        call str2vbg_cvb(cvbstot,dvbdet)
        call str2vbg_cvb(cvbsspn,evbdet)
        call scorr_cvb(cvbdet,dvbdet,evbdet)
      endif

c  Weights of CASSCF vector in VB basis
      if(lciweights)then
        if(.not.ifcasci_cvb())then
          write(6,'(2a)')' Warning - no CIWEIGHTS without CASSCF',
     >      ' vector.'
        else
          call str2vbc_cvb(cvb,cvbdet)
          call vb2cic_cvb(cvbdet,civb)
          call ciweight_cvb(civec,civbs,civb,citmp,civbh,
     >      orbs,sorbs,orbinv,owrk,gjorb,gjorb2,gjorb3)
        endif
      endif

      if(ishstruc.eq.1)then
        write(6,'(/,a)')' Overlap matrix between structures :'
        write(6,'(a)')  ' -----------------------------------'
        call mxprintd_cvb(sstruc,nvb,nvb,0)
        if(ifhamil_cvb())then
          call gaussj_cvb(orbs,gjorb)
          call transp_cvb(orbs,owrk,norb,norb)
          call gaussj_cvb(owrk,gjorb2)
          do 7250 k=1,nvb
          call fzero(sstruc2(1,k),nvb)
          sstruc2(k,k)=1d0
          call str2vbf_cvb(sstruc2(1,k),dvbdet)
          call vb2cif_cvb(dvbdet,civb)
          call applyt_cvb(civb,gjorb)
          call proj_cvb(civb)
          call applyh_cvb(civb)
          call proj_cvb(civb)
          call applyt_cvb(civb,gjorb2)
          call ci2vbg_cvb(civb,dvbdet)
7250      call vb2strg_cvb(dvbdet,sstruc2(1,k))
          do 7270 k=1,nvb
          do 7270 l=k+1,nvb
          sstruc2(k,l)=.5d0*(sstruc2(k,l)+sstruc2(l,k))
7270      sstruc2(l,k)=sstruc2(k,l)
          write(6,'(/,a)')' Hamiltonian matrix between structures :'
          write(6,'(a)')  ' ---------------------------------------'
          call mxprintd_cvb(sstruc2,nvb,nvb,0)
          call mxgendiag_cvb(sstruc2,sstruc,dvbdet,nvb)
          nr_print=min(20,nvb)
          if(nr_print.lt.nvb)
     >      write(6,'(/,a,i4,a)')' Printing',nr_print,' lowest roots:'
          do 7290 iroot=1,nr_print
          write(6,formroot)' Root no.',iroot,' energy=',dvbdet(iroot),
     >      ' :'
7290      call vecprint_cvb(sstruc2(1,iroot),nvb)
        endif
      endif

      call dscal_(nvb,cnrm,cvb,1)

      if(ip(5).ge.1.and.nirrep.gt.1)then
        write(6,'(/,a)')
     >    ' Symmetry contributions to total VB wavefunction :'
        write(6,'(a)')
     >    ' -------------------------------------------------'
        iimx=min(4,nirrep)
        write(6,formSymW)' Irreps 1 to',iimx,' :',(wsym(ii),ii=1,iimx)
        if(nirrep.gt.4)then
          iimx=min(8,nirrep)
          write(6,formSymW)' Irreps 5 to',iimx,' :',(wsym(ii),ii=5,iimx)
        endif
        if(lcalcevb)then
          write(6,'(/,a)')
     >      ' Energies for components > 1d-10 :'
          write(6,'(a)')
     >      ' ---------------------------------'
          iimx=min(4,nirrep)
          write(6,formSymW)' Irreps 1 to',iimx,' :',(esym(ii),ii=1,iimx)
          if(nirrep.gt.4)then
            iimx=min(8,nirrep)
            write(6,formSymW)' Irreps 5 to',iimx,' :',
     >        (esym(ii),ii=5,iimx)
          endif
        endif
      endif
      if(ip(5).ge.1)then
        write(6,'(/,a)')' One-electron density :'
        write(6,'(a)')' ----------------------'
        call mxprint_cvb(dmat,norb,norb,0)
        call mxdiag_cvb(dmat,occ,norb)
c  Sort NOs in order of increasing occ. numbers :
        do 8000 iorb=1,norb/2
        call dswap_(norb,dmat(1,iorb),1,dmat(1,norb+1-iorb),1)
        swp=occ(iorb)
        occ(iorb)=occ(norb+1-iorb)
8000    occ(norb+1-iorb)=swp
        write(6,'(/,a)')' Natural orbitals :'
        write(6,'(a)')' ------------------'
        call mxprint_cvb(dmat,norb,norb,0)
        write(6,'(/,a)')' Occupation numbers :'
        write(6,'(a)')' --------------------'
        call mxprint_cvb(occ,1,norb,0)
        occ_nel=0d0
        do i=1,norb
        occ_nel=occ_nel+occ(i)
        enddo
        if(abs(occ_nel-DBLE(nel)).gt..1d0)then
          write(6,*)' Error, sum of occupation numbers not equal to ',
     >      'number of electrons :',occ_nel,nel
          call abend_cvb()
        endif
      endif
      call molden_cvb()
      return
      end

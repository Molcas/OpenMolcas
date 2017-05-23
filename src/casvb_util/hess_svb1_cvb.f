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
      subroutine hess_svb1_cvb(orbs,
     >   civecp,civbs,civb,citmp,
     >   orbinv,sorbs,owrk,
     >   gjorb,gjorb2,gjorb3,
     >   dvbdet,
     >   grad1,grad2,hessorb,
     >   vec1,iorts,
     >   hessinp,hessout,owrk2,owrk3)
      implicit real*8 (a-h,o-z)
      logical orbopt2,strucopt2
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "frag_cvb.fh"
#include "fx_cvb.fh"
#include "malloc_cvb.fh"
      dimension orbs(norb,norb)
      dimension civecp(ndet),civbs(ndet),civb(ndet),citmp(ndet)
      dimension orbinv(norb,norb),sorbs(norb,norb),owrk(norb,norb)
      dimension gjorb(*),gjorb2(*),gjorb3(*)
      dimension dvbdet(ndetvb)
      dimension grad1(npr),grad2(npr)
      dimension hessorb(nprorb,nprorb)
      dimension iorts(2,nort)

c  VEC1 dimension is MAX(NPRORB,NDETVB)
      dimension vec1(*)
      dimension hessinp(npr),hessout(npr)

      hess_orb_nrm=dnrm2_(nprorb,hessinp,1)
      orbopt2=hess_orb_nrm.gt.1d-10
      if(nprvb.gt.0) then
        hess_ci_nrm=dnrm2_(nprvb,hessinp(nprorb+1),1)
      else
        hess_ci_nrm=zero
      endif
      strucopt2=strucopt.and.hess_ci_nrm.gt.1d-10
      if(orbopt2.and..not.strucopt2)n_orbhess=n_orbhess+1
      if(strucopt2.and..not.orbopt2)n_cihess=n_cihess+1

      call transp_cvb(orbs,owrk,norb,norb)
      call mxattb_cvb(orbs,orbs,norb,norb,norb,sorbs)

      call fzero(hessout,npr)
      if(orbopt2)call mxatb_cvb(hessorb,hessinp,
     >  nprorb,nprorb,1,hessout)
c  Combinations of gradients :
      g1f=ddot_(npr,grad1,1,hessinp,1)
      g2f=ddot_(npr,grad2,1,hessinp,1)
      fac1=g1f*oaa3+g2f*aa2
      fac2=g1f*aa2
      call daxpy_(npr,fac1,grad1,1,hessout,1)
      call daxpy_(npr,fac2,grad2,1,hessout,1)

      if(orbopt2.and.strucopt)then
        call mxunfold_cvb(hessinp(1),owrk,norb)
        call transp_cvb(owrk,owrk,norb,norb)
        call cizero_cvb(citmp)
        call oneexc_cvb(civecp,citmp,owrk,.true.,2)
        call mkgrd_cvb(civb,citmp,vec1,dvbdet,npr,.false.)
        call daxpy_(nprvb,aa1,vec1(nprorb+1),1,
     >    hessout(nprorb+1),1)
        call transp_cvb(owrk,owrk2,norb,norb)
        if(.not.(proj.or.projcas))then
          call mxatb_cvb(sorbs,owrk2,norb,norb,norb,owrk3)
          call mxatb_cvb(owrk3,orbinv,norb,norb,norb,owrk2)
          call addvec(owrk,owrk2,owrk,norb*norb)
        endif
        call cizero_cvb(citmp)
        call oneexc_cvb(civbs,citmp,owrk,.true.,2)
        call mkgrd_cvb(civb,citmp,vec1,dvbdet,npr,.false.)
        call daxpy_(nprvb,oaa2,vec1(nprorb+1),1,
     >    hessout(nprorb+1),1)
      endif
      if(strucopt2)then
        call str2vbf_cvb(hessinp(1+nprorb),dvbdet)
        call vb2cif_cvb(dvbdet,citmp)
c  Structure coeff. <-> orbital
        call mkgrd_cvb(citmp,civbs,vec1,dvbdet,nprorb,.true.)
        call daxpy_(nprorb,oaa2,vec1,1,hessout,1)
        call mkgrd_cvb(citmp,civecp,vec1,dvbdet,nprorb,.true.)
        call daxpy_(nprorb,aa1,vec1,1,hessout,1)
        if(proj.or.projcas)
     >    call oneexc_cvb(civb,citmp,hessinp,.false.,1)
c  Structure coeff. <-> all
        call applyts_cvb(citmp,orbs,gjorb,gjorb2,gjorb3)
        call mkgrd_cvb(civb,citmp,vec1,dvbdet,npr,.true.)
        call daxpy_(npr,oaa2,vec1,1,hessout,1)
c  2nd-order term for structure coefficients
        if(nfrag.gt.1)then
          call str2vbf_cvb(hessinp(1+nprorb),dvbdet)
          i1 = mstackr_cvb(ndetvb)
          i2 = mstackr_cvb(nvb)
          call ci2ordr_cvb(civbs,dvbdet,w(i1))
          call vb2strg_cvb(w(i1),w(i2))
          call daxpy_(nvb,oaa2,w(i2),1,hessout(1+nprorb),1)
          call ci2ordr_cvb(civecp,dvbdet,w(i1))
          call vb2strg_cvb(w(i1),w(i2))
          call daxpy_(nvb,aa1,w(i2),1,hessout(1+nprorb),1)
          call mfreer_cvb(i1)
        endif
      elseif(proj.or.projcas)then
        call cizero_cvb(citmp)
        call oneexc_cvb(civb,citmp,hessinp,.false.,1)
c  Structure coeff. <-> all
        call applyts_cvb(citmp,orbs,gjorb,gjorb2,gjorb3)
        call mkgrd_cvb(civb,citmp,vec1,dvbdet,npr,.true.)
        call daxpy_(npr,oaa2,vec1,1,hessout,1)
      endif

      if(orbopt2.and.nort.gt.0)then
c  Non-linear correction for orthogonality constraints :
        call fmove(sorbs,owrk,norb*norb)
        call mxinv_cvb(owrk,norb)
        do 100 iort=1,nort
        iorb=iorts(1,iort)
        jorb=iorts(2,iort)
        corr1=zero
        do 200 korb=1,norb
        ki=korb+(iorb-1)*(norb-1)
        if(korb.gt.iorb)ki=ki-1
        kj=korb+(jorb-1)*(norb-1)
        if(korb.gt.jorb)kj=kj-1
        if(korb.ne.iorb)corr1=corr1+owrk(jorb,korb)*
     >    (aa1*grad2(ki)+.5d0*oaa2*grad1(ki))
200     if(korb.ne.jorb)corr1=corr1+owrk(iorb,korb)*
     >    (aa1*grad2(kj)+.5d0*oaa2*grad1(kj))
        corr1=-.5d0*corr1
        do 300 korb=1,norb
        if(korb.eq.iorb)goto 300
        ki=korb+(iorb-1)*(norb-1)
        if(korb.gt.iorb)ki=ki-1
        do 400 lorb=1,norb
        if(lorb.eq.jorb)goto 400
        lj=lorb+(jorb-1)*(norb-1)
        if(lorb.gt.jorb)lj=lj-1
        hessout(ki)=hessout(ki)+sorbs(korb,lorb)*corr1*hessinp(lj)
        hessout(lj)=hessout(lj)+sorbs(korb,lorb)*corr1*hessinp(ki)
400     continue
300     continue
100     continue
      endif
      return
      end

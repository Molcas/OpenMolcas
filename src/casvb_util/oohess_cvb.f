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
      subroutine oohess_cvb(orbs,civecp,civbs,civb,
     >   orbinv,sorbs,owrk,
     >   grad2,gradx,hessorb,hesst)
c  Evaluate "cheap" orbital <-> orbital part of hessian :
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "fx_cvb.fh"
      dimension orbs(norb,norb)
      dimension civecp(ndet),civbs(ndet),civb(ndet)
      dimension orbinv(norb,norb),sorbs(norb,norb),owrk(norb,norb)
      dimension grad2(npr),gradx(norb,norb)
      dimension hessorb(nprorb,nprorb)
      dimension hesst(norb*norb,norb*norb)

      if(icrit.eq.1)then
        oaa2_use=oaa2
        aa1_use=aa1
      elseif(icrit.eq.2)then
        oaa2_use=f2
        aa1_use=f1
      endif

      call fzero(hessorb,nprorb*nprorb)
      if(icrit.eq.1.and..not.(proj.or.projcas))then
        call dev2b_cvb(civbs,civecp,civb,hessorb,hesst,
     >    oaa2_use,aa1_use,gradx,grad2)

        call mxattb_cvb(orbs,orbs,norb,norb,norb,sorbs)
        call fmove(sorbs,orbinv,norb*norb)
        call mxinv_cvb(orbinv,norb)

        do 100 jorb=1,norb
        do 100 iorb=1,norb
        iprm=iorb+(jorb-1)*norb
        call mxatb_cvb(orbinv,hesst(1,iprm),norb,norb,norb,owrk)
100     call mxatb_cvb(owrk,sorbs,norb,norb,norb,hesst(1,iprm))
        iprm1=0
        do 200 iorb=1,norb
        do 200 jorb=1,norb
        if(jorb.eq.iorb)goto 200
        iprm1=iprm1+1
        ifr1=jorb+(iorb-1)*norb
        iprm2=0
        do 300 korb=1,norb
        do 300 lorb=1,norb
        if(lorb.eq.korb)goto 300
        iprm2=iprm2+1
        ifr2=korb+(lorb-1)*norb
        if(iprm2.le.iprm1)then
          hessorb(iprm2,iprm1)=hessorb(iprm2,iprm1)+
     >      oaa2_use*hesst(ifr2,ifr1)
          hessorb(iprm1,iprm2)=hessorb(iprm2,iprm1)
        endif
300     continue
200     continue
      elseif(icrit.eq.1)then
        call dev2a_cvb(civbs,civecp,civb,hessorb,oaa2_use,aa1_use)
      else
        call cidaxpy_cvb(-ww/ovraa,civbs,civecp)
        call cizero_cvb(civbs)
        call dev2c_cvb(civecp,civb,hessorb,aa1_use)
      endif
      return
      end

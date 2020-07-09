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
      subroutine mkorbfree2_cvb(orbs,north,corth,
     >  irels,relorb,ifxorb,iorts,irots,
     >  trprm,owrk,owrk2,orbinv,idel)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

      dimension orbs(norb,norb)
      dimension north(norb),corth(norb,niorth)
      dimension irels(2,nijrel),relorb(norb,norb,nijrel)
      dimension ifxorb(norb),iorts(2,nort),irots(2,ndrot)
      dimension trprm(nprorb,nprorb)
      dimension owrk(norb,norb),owrk2(norb,norb)
      dimension orbinv(norb,norb),idel(nprorb)
      dimension dum(1)
      save thresh
      data thresh/1.d-7/

      if(orbfr_is_unit)then
        nfrorb=nprorb
        goto 1000
      endif
      call fzero(trprm,nprorb*nprorb)
      call izero(idel,nprorb)
      call fmove_cvb(orbs,orbinv,norb*norb)
      call mxinv_cvb(orbinv,norb)

      nc=0
      ioffs=0
      do 100 iorb=1,norb
      if(north(iorb).gt.0.and.ifxorb(iorb).ne.1)then
c  Transform simple constraints to basis of VB orbitals :
        call mxattb_cvb(orbs,corth(1,1+ioffs),norb,norb,north(iorb),
     >    owrk)
        call span_cvb(owrk,north(iorb),ncon,dum,norb,0)
        ishift=(iorb-1)*(norb-1)
        do 300 icon=1,ncon
        nc=nc+1
        i2=0
        do 301 i=1,norb
        if(i.ne.iorb)then
          i2=i2+1
          trprm(i2+ishift,nc)=owrk(i,icon)
        endif
301     continue
300     continue
      elseif(ifxorb(iorb).eq.1)then
        ishift=(iorb-1)*(norb-1)
        do 350 icon=1,norb-1
        nc=nc+1
        trprm(icon+ishift,nc)=one
350     continue
      endif
      ioffs=ioffs+north(iorb)
100   continue

      call mxattb_cvb(orbs,orbs,norb,norb,norb,owrk)
      call ortelim_cvb(trprm,iorts,irots,owrk,
     >  nc,nprorb,norb*(norb-1),nrem)
      call izero(idel,nprorb)
      do 500 i=1,nrem
      idel(i)=1
500   continue

      do 600 irel=1,nijrel
      iorb=irels(1,irel)
      jorb=irels(2,irel)
      call mxatb_cvb(relorb(1,1,irel),orbs,norb,norb,norb,owrk)
      call mxatb_cvb(orbinv,owrk,norb,norb,norb,owrk2)
      if(abs(abs(owrk2(iorb,jorb))-one).gt.thresh)then
        write(6,*)' Transformation matrix cannot be correct !'
        call mxprint_cvb(owrk2,norb,norb,0)
        call abend_cvb()
      endif
      ishift=(jorb-1)*(norb-1)
      ishift2=(iorb-1)*(norb-1)
      i2=0
      do 700 i=1,norb
      if(i.eq.iorb)goto 700
      i2=i2+1
      j2=0
      do 720 j=1,norb
      if(j.eq.jorb)goto 720
      j2=j2+1
      do 740 iprm=1,nprorb
      trprm(i2+ishift2,iprm)=trprm(i2+ishift2,iprm)
     >  +owrk2(i,j)*trprm(j2+ishift,iprm)
740   continue
720   continue
700   continue
      ioff=1+(iorb-1)*(norb-1)
      ioff2=1+iorb*(norb-1)
      nl1=(iorb-1)*(norb-1)
      nl2=(norb-iorb)*(norb-1)
      do 800 i=nrem+1,nprorb
      sum1=ddot_(norb-1,trprm(ioff,i),1,trprm(ioff,i),1)
      sum2=ddot_(nl1,trprm(1,i),1,trprm(1,i),1)
      if(nl2.gt.0) sum2=sum2
     >    +ddot_(nl2,trprm(ioff2,i),1,trprm(ioff2,i),1)
      if(sum1.gt.thresh.and.sum2.lt.thresh)idel(i)=1
800   continue
600   continue
      nfrorb=0
      do 900 i=1,norb*(norb-1)
      if(idel(i).ne.1)then
        nfrorb=nfrorb+1
        call fmove_cvb(trprm(1,i),trprm(1,nfrorb),nprorb)
      endif
900   continue
      call fzero(trprm(1,nfrorb+1),(nprorb-nfrorb)*nprorb)
      call nize_cvb(trprm,nfrorb,dum,nprorb,0,0)

1000  nfr=nfrorb+nfrvb
      orbopt=(nfrorb.ne.0)
      return
      end

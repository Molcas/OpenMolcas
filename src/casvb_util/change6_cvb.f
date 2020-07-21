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
      subroutine change6_cvb()
      implicit real*8 (a-h,o-z)
      logical changed
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "rls_cvb.fh"
#include "davtune_cvb.fh"
      save icase

      changed=.false.
      if(changed)call touch_cvb('CIFREE')

      nprorb=norb*(norb-1)
      if(strucopt)then
        nprvb=nvb
        npr=nprorb+nvb
      elseif(.not.(imethod.eq.4.or.imethod.eq.6))then
        nprvb=0
        npr=nprorb
      else
        npr=0
      endif
      if(chpcmp_cvb(npr))changed=.true.
      if((.not.(imethod.eq.4.or.imethod.eq.6)).and.ifinish.eq.0)then
c  Standard 2nd-order procedure :
        icase=1
      elseif(imethod.eq.4.and.icrit.eq.1.and.ifinish.eq.0)then
c  Overlap-based Davidson
        icase=2
      elseif(imethod.eq.4.and.icrit.eq.2.and.ifinish.eq.0)then
c  Energy-based Davidson
        icase=3
      elseif(imethod.eq.6.or.ifinish.eq.1.or.ifinish.eq.2)then
c  No arrays needed
        icase=4
      else
        icase=5
      endif
      if(chpcmp_cvb(icase))changed=.true.
      if(changed)call touch_cvb('MEM6')
      return

      entry chop6_cvb()
      if(release(6))call mfreer_cvb(lp(1))
      release(6)=.true.
      release(7)=.false.
      lp(1) = mstackr_cvb(0)

      call setcnt2_cvb(6,0)
      if(icase.eq.1)then
c  Standard non-linear optimization procedure :
        lp(1)= mstackr_cvb(norb*norb+nvb+1+mxirrep)
        lp(2)= mstackr_cvb(npr)
        lq(3)= mstackr_cvb(nprorb*nprorb)
        lq(4)= mstackr_cvb(norb**4)
        lp(5)= mstackr_cvb(npr)
        lp(6)= mstackr_cvb(npr)
        lq(7)= mstackr_cvb(npr)
        lq(8)= mstackr_cvb(npr)
        lq(9)= mstackr_cvb(norb*norb)
c  Vec1 work array
        lq(10)= mstackr_cvb(max(npr,ndetvb))
      elseif(icase.eq.2)then
c  Overlap-based Davidson optimization :
        iremain=mavailr_cvb()
        maxdav=min(mxiter,nvb,mxdav)

        memwrk=ndetvb+5*norb*norb+3*ihlf_cvb(norb+2*norb*norb)
        do 1 idav=maxdav,1,-1
c  NEED is approx req. memory :
        need=2*nvb*idav+2*nvb+idav+1000+memwrk
        if(need.lt.iremain)goto 2
1       continue
        idav=0
        if(nvb.eq.0)then
          need=1000+memwrk
          if(need.lt.iremain)goto 2
        endif
        write(6,*)' Not enough memory for Davidson!',need,iremain
        call abend_cvb()
2       maxdav=idav

      elseif(icase.eq.3)then
c  Energy-based Davidson optimization :
        iremain=mavailr_cvb()
        maxdav=min(mxiter,nvb,mxdav)

        mem_applyh=ndet+neread
        ncimx=0
        do ir=1,nirrep
        ncimx=max(ncimx,ncivb(ir))
        enddo
        if(ncimx.ne.ndet)mem_applyh=mem_applyh+ncimx
        memwrk=ndetvb+3*norb*norb+2*ihlf_cvb(norb+2*norb*norb)

        do 11 idav=maxdav,1,-1
c  NEED is approx req. memory :
        need=3*nvb*idav+nvb+idav*(2*idav+3)+1000+mem_applyh+memwrk
        if(need.lt.iremain)goto 12
11      continue
        idav=0
        if(nvb.eq.0)then
          need=1000+memwrk
          if(need.lt.iremain)goto 12
        endif
        write(6,*)' Not enough memory for Davidson!',need,iremain
        call abend_cvb()
12      maxdav=idav

      elseif(icase.eq.4)then
c  Wavefunction analysis :
        mstackr_cvb0=mstackr_cvb(0)
        if(((.not.variat).or.endvar).and.
     >    (ivbweights.gt.1.or.ishstruc.eq.1))then
          lp(1)= mstackr_cvb(nvb*nvb)
          lp(2)= mstackr_cvb(nvb*nvb)
        else
          lp(1)= mstackr_cvb0
          lp(2)= mstackr_cvb0
        endif
        do i=3,11
        lp(i)= mstackr_cvb0
        enddo
      endif
      return
      end

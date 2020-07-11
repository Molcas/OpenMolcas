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
c  *********************************************************************
c  *                                                                   *
c  *  CNFPRT   := Print configurations.                                *
c  *                                                                   *
c  *********************************************************************
      subroutine cnfprt_cvb(iconfs,nconf1,nel1)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"


#include "malloc_cvb.fh"
      dimension iconfs(noe,nconf1)

      i1 = mstacki_cvb(noe)
c  Main loop over configurations :
      do 100 iconf=1,nconf1
c  Prepare iw(i1) for print
      ioffs=i1-1
      do 200 iorb=1,norb
      if(iconfs(iorb,iconf).eq.2)then
        iw(1+ioffs)=iorb
        iw(2+ioffs)=iorb
        ioffs=ioffs+2
      endif
200   continue
      do 300 iorb=1,norb
      if(iconfs(iorb,iconf).eq.1)then
        iw(1+ioffs)=iorb
        ioffs=ioffs+1
      endif
300   continue
      write(6,'(i8,a,20i3)')iconf,'   =>  ',(iw(ii+i1-1),ii=1,nel1)
100   continue
      call mfreei_cvb(i1)
      return
      end
      function nvb_cvb(kbasis_loc)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"


#include "frag_cvb.fh"

      ndetvb=0
      ndetvb2=0
      nvbr=0
      do 100 ifrag=1,nfrag
      if(kbasis_loc.ne.6)then
        nvb_fr(ifrag)=nvbr_fr(ifrag)
      else
        if(absym(1))then
          nvb_fr(ifrag)=ndetvb2_fr(ifrag)
        else
          nvb_fr(ifrag)=ndetvb_fr(ifrag)
        endif
      endif
      ndetvb=ndetvb+ndetvb_fr(ifrag)
      ndetvb2=ndetvb2+ndetvb2_fr(ifrag)
      nvbr=nvbr+nvbr_fr(ifrag)
100   continue

      if(kbasis_loc.ne.6)then
        nvb_loc=nvbr
      else
        if(absym(1))then
          nvb_loc=ndetvb2
        else
          nvb_loc=ndetvb
        endif
      endif
      nvb_cvb=nvb_loc
      return
      end
      function ifns_cvb(nel1,nalf1,kbasis1)
      implicit real*8 (a-h,o-z)

      nbet1=nel1-nalf1
      if(nbet1.gt.nalf1)then
        nsw=nalf1
        nalf1=nbet1
        nbet1=nsw
      endif
      if(kbasis1.ne.6)then
        call icomb_cvb(nel1,nbet1,iretval1)
        call icomb_cvb(nel1,nbet1-1,iretval2)
        ifn=iretval1-iretval2
      else
        call icomb_cvb(nel1,nalf1,ifn)
        if(nalf1.eq.nbet1)ifn=(ifn+1)/2
      endif
      ifns_cvb=ifn
      return
      end
      function ndet_cvb(nel1,nalf1)
      implicit real*8 (a-h,o-z)

      call icomb_cvb(nel1,nalf1,nd)
      ndet_cvb=nd
      return
      end

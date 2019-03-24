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
      subroutine ortelim_cvb(trprm,iorts,irots,sorbs,
     >  nc,npr1,norbprm,nrem)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension trprm(npr1,npr1),iorts(2,nort),irots(2,ndrot)
      dimension sorbs(norb,norb)

      i1 = mstackrz_cvb(norbprm*max(nc+nort+ndrot,norbprm))
      i1ff=i1-1
      do 100 i=1,nc
100   call fmove_cvb(trprm(1,i),w(1+(i-1)*norbprm+i1ff),norbprm)
      do 200 iort=1,nort
      iorb=iorts(1,iort)
      jorb=iorts(2,iort)
      do 200 korb=1,norb
      ki=korb+(iorb-1)*(norb-1)
      if(korb.gt.iorb)ki=ki-1
      kj=korb+(jorb-1)*(norb-1)
      if(korb.gt.jorb)kj=kj-1
      if(korb.ne.iorb)w(ki+(iort+nc-1)*norbprm+i1ff)=sorbs(korb,jorb)
200   if(korb.ne.jorb)w(kj+(iort+nc-1)*norbprm+i1ff)=sorbs(korb,iorb)
      do 300 irot=1,ndrot
      iorb=irots(1,irot)
      jorb=irots(2,irot)
      do 300 korb=1,norb
      ki=korb+(iorb-1)*(norb-1)
      if(korb.gt.iorb)ki=ki-1
      kj=korb+(jorb-1)*(norb-1)
      if(korb.gt.jorb)kj=kj-1
      if(korb.ne.iorb)w(ki+(irot+nc+nort-1)*norbprm+i1ff)=
     >  sorbs(korb,jorb)
300   if(korb.ne.jorb)w(kj+(irot+nc+nort-1)*norbprm+i1ff)=
     >  -sorbs(korb,iorb)
      call span_cvb(w(i1),nc+nort+ndrot,nrem,dum,norbprm,0)
      call compl_cvb(w(i1),nrem,norbprm)

      call fzero(trprm,npr1*npr1)
      do 400 i=1,norbprm
400   call fmove_cvb(w(1+(i-1)*norbprm+i1ff),trprm(1,i),norbprm)

      call mfreer_cvb(i1)
      return
      end

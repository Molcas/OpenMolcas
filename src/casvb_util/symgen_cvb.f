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
      subroutine symgen_cvb(nalf1,nbet1,nda1,ndb1,
     >  isymalf,isymbet,iasyind,ibsyind,
     >  ialfsym,ibetsym,irpdet,irpalf,irpbet,
     >  mingrph,maxgrph,nk,locc,lunocc,xalf,xbet,icount)
      implicit real*8 (a-h,o-w,y-z),integer(x)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

      dimension isymalf(nda1),isymbet(ndb1)
      dimension iasyind(0:mxirrep),ibsyind(0:mxirrep)
      dimension ialfsym(nda1),ibetsym(ndb1)
      dimension irpdet(mxirrep),irpalf(mxirrep),irpbet(mxirrep)
      dimension mingrph(0:norb),maxgrph(0:norb)
      dimension nk(0:norb),locc(norb+1),lunocc(norb+1)
      dimension xalf(0:norb,0:nalf1),xbet(0:norb,0:nbet1)
      dimension icount(mxirrep)

c Alpha loop:
      call izero(irpalf,mxirrep)
      do 100 iorb=0,norb
      mingrph(iorb)=max(iorb-norb+nalf1,0)
100   maxgrph(iorb)=min(iorb,nalf1)
      call weight_cvb(xalf,mingrph,maxgrph,nalf1,norb)
      call imove_cvb(maxgrph,nk,norb+1)
      call occupy_cvb(nk,norb,locc,lunocc)
      index=1
200   continue
      irp=1
      do 250 ia=1,nalf1
250   irp=md2h(irp,ityp(locc(ia)))
      irpalf(irp)=irpalf(irp)+1
      ialfsym(index)=irp
      call loind_cvb(norb,nalf1,nk,mingrph,maxgrph,
     >                       locc,lunocc,index,xalf,*200)
      iasyind(0)=0
      do 260 irp=1,mxirrep
260   iasyind(irp)=iasyind(irp-1)+irpalf(irp)
      call izero(icount,mxirrep)
      do 275 ida=1,nda1
      irrep=ialfsym(ida)
      icount(irrep)=icount(irrep)+1
275   isymalf(icount(irrep)+iasyind(irrep-1))=ida

c Beta loop:
      call izero(irpbet,mxirrep)
      do 300 iorb=0,norb
      mingrph(iorb)=max(iorb-norb+nbet1,0)
300   maxgrph(iorb)=min(iorb,nbet1)
      call weight_cvb(xbet,mingrph,maxgrph,nbet1,norb)
      call imove_cvb(maxgrph,nk,norb+1)
      call occupy_cvb(nk,norb,locc,lunocc)
      index=1
400   continue
      irp=1
      do 450 ib=1,nbet1
450   irp=md2h(irp,ityp(locc(ib)))
      irpbet(irp)=irpbet(irp)+1
      ibetsym(index)=irp
      call loind_cvb(norb,nbet1,nk,mingrph,maxgrph,
     >                       locc,lunocc,index,xbet,*400)
      ibsyind(0)=0
      do 460 irp=1,mxirrep
460   ibsyind(irp)=ibsyind(irp-1)+irpbet(irp)
      call izero(icount,mxirrep)
      do 475 idb=1,ndb1
      irrep=ibetsym(idb)
      icount(irrep)=icount(irrep)+1
475   isymbet(icount(irrep)+ibsyind(irrep-1))=idb

      do 500 irp=1,mxirrep
      irpdet(irp)=0
      do 500 jrp=1,mxirrep
500   irpdet(irp)=irpdet(irp)+irpalf(jrp)*irpbet(md2h(irp,jrp))
      return
      end

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
* Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
*               1996-2006, David L. Cooper                             *
************************************************************************
      subroutine symgen_cvb(nalf1,nbet1,nda1,ndb1,
     >  isymalf,isymbet,iasyind,ibsyind,
     >  ialfsym,ibetsym,irpdet,irpalf,irpbet,
     >  mingrph,maxgrph,nk,locc,lunocc,xalf,xbet,icount)
      implicit real*8 (a-h,o-w,y-z),integer(x)
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
      maxgrph(iorb)=min(iorb,nalf1)
100   continue
      call weight_cvb(xalf,mingrph,maxgrph,nalf1,norb)
      call imove_cvb(maxgrph,nk,norb+1)
      call occupy_cvb(nk,norb,locc,lunocc)
      index=1
200   continue
      irp=1
      do 250 ia=1,nalf1
      irp=md2h(irp,ityp(locc(ia)))
250   continue
      irpalf(irp)=irpalf(irp)+1
      ialfsym(index)=irp
      call loind_cvb(norb,nalf1,nk,mingrph,maxgrph,
     >                       locc,lunocc,index,xalf,*200)
      iasyind(0)=0
      do 260 irp=1,mxirrep
      iasyind(irp)=iasyind(irp-1)+irpalf(irp)
260   continue
      call izero(icount,mxirrep)
      do 275 ida=1,nda1
      irrep=ialfsym(ida)
      icount(irrep)=icount(irrep)+1
      isymalf(icount(irrep)+iasyind(irrep-1))=ida
275   continue

c Beta loop:
      call izero(irpbet,mxirrep)
      do 300 iorb=0,norb
      mingrph(iorb)=max(iorb-norb+nbet1,0)
      maxgrph(iorb)=min(iorb,nbet1)
300   continue
      call weight_cvb(xbet,mingrph,maxgrph,nbet1,norb)
      call imove_cvb(maxgrph,nk,norb+1)
      call occupy_cvb(nk,norb,locc,lunocc)
      index=1
400   continue
      irp=1
      do 450 ib=1,nbet1
      irp=md2h(irp,ityp(locc(ib)))
450   continue
      irpbet(irp)=irpbet(irp)+1
      ibetsym(index)=irp
      call loind_cvb(norb,nbet1,nk,mingrph,maxgrph,
     >                       locc,lunocc,index,xbet,*400)
      ibsyind(0)=0
      do 460 irp=1,mxirrep
      ibsyind(irp)=ibsyind(irp-1)+irpbet(irp)
460   continue
      call izero(icount,mxirrep)
      do 475 idb=1,ndb1
      irrep=ibetsym(idb)
      icount(irrep)=icount(irrep)+1
      isymbet(icount(irrep)+ibsyind(irrep-1))=idb
475   continue

      do 500 irp=1,mxirrep
      irpdet(irp)=0
      do 501 jrp=1,mxirrep
      irpdet(irp)=irpdet(irp)+irpalf(jrp)*irpbet(md2h(irp,jrp))
501   continue
500   continue
      return
      end

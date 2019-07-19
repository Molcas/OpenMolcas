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
      subroutine permvb2_cvb(v1,iperm,vb,iapr,ixapr,
     >  xalf,xbet,mingrph,maxgrph,
     >  nk,locc,lunocc,inewocc,inocc2,negs,
     >  inda,phsa,indb,phsb,v2,ialg)
      implicit real*8 (a-h,o-w,y-z),integer(x)
      logical vb
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

      dimension iperm(norb)
      dimension iapr(ndetvb),ixapr(nda+1)
      dimension xalf(0:norb,0:nalf),xbet(0:norb,0:nbet)
      dimension mingrph(0:norb),maxgrph(0:norb)
      dimension nk(0:norb),locc(norb+1),lunocc(norb+1)
      dimension inewocc(norb),inocc2(norb),negs(norb)
      dimension inda(nda),phsa(nda),indb(ndb),phsb(ndb)
c  V1 is dimensioned either NDET or NDETVB according to CI/VB
c  V2 is dimensioned NDET/NDA or NDETVB according to CI/VB
      dimension v1(*),v2(*)

c  Some tests of permutation
c  Valid ?
      call izero(negs,norb)
      do 10 i=1,norb
      iprm=abs(iperm(i))
      if(iprm.lt.1.or.iprm.gt.norb)then
        write(6,*)' Illegal orbital permutation!'
        call abend_cvb()
      endif
10    negs(iprm)=negs(iprm)+1
      do 20 iorb=1,norb
      if(negs(iorb).ne.1)then
        write(6,*)' Illegal orbital permutation!'
        call abend_cvb()
      endif
20    continue
c  Return if identity
      do 30 iorb=1,norb
30    if(iperm(iorb).ne.iorb)goto 35
      return
35    continue
c  Use IALG=2 if only phase changes
      do 40 iorb=1,norb
40    if(abs(iperm(iorb)).ne.iorb)goto 45
      ialg=2
45    continue
      call izero(negs,norb)
      do 50 i=1,norb
50    if(iperm(i).lt.0)negs(abs(iperm(i)))=1
c Alpha loop:
      call izero(inocc2,norb)
      do 100 iorb=0,norb
      mingrph(iorb)=max(iorb-norb+nalf,0)
100   maxgrph(iorb)=min(iorb,nalf)
      call weight_cvb(xalf,mingrph,maxgrph,nalf,norb)
      call imove_cvb(maxgrph,nk,norb+1)
      call occupy_cvb(nk,norb,locc,lunocc)
      index=1
200   continue
      call izero(inewocc,norb)
      do 225 ialf=1,nalf
225   inewocc(abs(iperm(locc(ialf))))=ialf
      ineg=0
      ia=0
      do 250 iorb=1,norb
      if(inewocc(iorb).ne.0)then
        ia=ia+1
        inocc2(ia)=inewocc(iorb)
        inewocc(iorb)=1
        if(negs(iorb).eq.1)ineg=ineg+1
      endif
250   continue
      if(mod(ineg,2).eq.0)then
        phsa(index)= party_cvb(inocc2,nalf)
      else
        phsa(index)=-party_cvb(inocc2,nalf)
      endif
      inda(index)=indget_cvb(inewocc,nalf,norb,xalf)

      call loind_cvb(norb,nalf,nk,mingrph,maxgrph,
     >                       locc,lunocc,index,xalf,*200)
c Beta loop:
      call izero(inocc2,norb)
      do 400 iorb=0,norb
      mingrph(iorb)=max(iorb-norb+nbet,0)
400   maxgrph(iorb)=min(iorb,nbet)
      call weight_cvb(xbet,mingrph,maxgrph,nbet,norb)
      call imove_cvb(maxgrph,nk,norb+1)
      call occupy_cvb(nk,norb,locc,lunocc)
      index=1
500   continue
      call izero(inewocc,norb)
      do 525 ibet=1,nbet
525   inewocc(abs(iperm(locc(ibet))))=ibet
      ineg=0
      ib=0
      do 550 iorb=1,norb
      if(inewocc(iorb).ne.0)then
        ib=ib+1
        inocc2(ib)=inewocc(iorb)
        inewocc(iorb)=1
        if(negs(iorb).eq.1)ineg=ineg+1
      endif
550   continue
      if(mod(ineg,2).eq.0)then
        phsb(index)= party_cvb(inocc2,nbet)
      else
        phsb(index)=-party_cvb(inocc2,nbet)
      endif
      indb(index)=indget_cvb(inewocc,nbet,norb,xbet)

      call loind_cvb(norb,nbet,nk,mingrph,maxgrph,
     >                       locc,lunocc,index,xbet,*500)

      if(vb)then
        call fzero(v2,ndetvb)
        do 1000 ia=1,nda
        iato=inda(ia)
        do 1000 ixa=ixapr(ia),ixapr(ia+1)-1
        ib=iapr(ixa)
        ibto=indb(ib)
        do 1100 ixato=ixapr(iato),ixapr(iato+1)-1
        if(iapr(ixato).eq.ibto)goto 1200
1100    continue
c  Shouldn't get here ...
        write(6,'(a,100i3)')
     >    ' Error, VB determinants not closed under permutation :',iperm
        call abend_cvb()
1200    continue
1000    v2(ixa)=phsa(ia)*phsb(ib)*v1(ixato)
        call fmove_cvb(v2,v1,ndetvb)
      elseif(ialg.eq.1)then
c  Brute force strategy if enough memory (x1.5 faster) :
        do 2000 ib=1,ndb
        iboff=(ib-1)*nda
        inboff=(indb(ib)-1)*nda
        do 2000 ia=1,nda
2000    v2(ia+iboff)=phsa(ia)*phsb(ib)*v1(inda(ia)+inboff)
        call fmove_cvb(v2,v1,ndet)
      elseif(ialg.eq.2)then
c  More-or-less in-place update of V1 :
        do 3000 ia=1,nda
        if(ia.eq.inda(ia))then
          if(phsa(ia).eq.-one)then
            ioffs=ia-nda
            do 3100 ib=1,ndb
3100        v1(ib*nda+ioffs)=-v1(ib*nda+ioffs)
          endif
        elseif(inda(ia).ne.0)then
c  Cyclic permutation involving IA :
          ioffs=ia-nda
          do 3300 ib=1,ndb
3300      v2(ib)=v1(ib*nda+ioffs)
          iat=ia
3400      continue
          if(phsa(iat).eq.one)then
            ioffs1=iat-nda
            ioffs2=inda(iat)-nda
            do 3500 ib=1,ndb
3500        v1(ib*nda+ioffs1)=v1(ib*nda+ioffs2)
          else
            ioffs1=iat-nda
            ioffs2=inda(iat)-nda
            do 3600 ib=1,ndb
3600        v1(ib*nda+ioffs1)=-v1(ib*nda+ioffs2)
          endif
          iatold=iat
          iat=inda(iat)
          inda(iatold)=0
          if(inda(iat).ne.ia)goto 3400
          if(phsa(iat).eq.one)then
            ioffs=iat-nda
            do 3700 ib=1,ndb
3700        v1(ib*nda+ioffs)=v2(ib)
          else
            ioffs=iat-nda
            do 3800 ib=1,ndb
3800        v1(ib*nda+ioffs)=-v2(ib)
          endif
          inda(iat)=0
        endif
3000    continue
        do 4000 ib=1,ndb
        if(ib.eq.indb(ib))then
          if(phsb(ib).eq.-one)then
            ioffs=(ib-1)*nda
            do 4100 ia=1,nda
4100        v1(ia+ioffs)=-v1(ia+ioffs)
          endif
        elseif(indb(ib).ne.0)then
c  Cyclic permutation involving IB :
          ioffs=(ib-1)*nda
          do 4300 ia=1,nda
4300      v2(ia)=v1(ia+ioffs)
          ibt=ib
4400      continue
          if(phsb(ibt).eq.one)then
            ioffs1=(ibt-1)*nda
            ioffs2=(indb(ibt)-1)*nda
            do 4500 ia=1,nda
4500        v1(ia+ioffs1)=v1(ia+ioffs2)
          else
            ioffs1=(ibt-1)*nda
            ioffs2=(indb(ibt)-1)*nda
            do 4600 ia=1,nda
4600        v1(ia+ioffs1)=-v1(ia+ioffs2)
          endif
          ibtold=ibt
          ibt=indb(ibt)
          indb(ibtold)=0
          if(indb(ibt).ne.ib)goto 4400
          if(phsb(ibt).eq.one)then
            ioffs=(ibt-1)*nda
            do 4700 ia=1,nda
4700        v1(ia+ioffs)=v2(ia)
          else
            ioffs=(ibt-1)*nda
            do 4800 ia=1,nda
4800        v1(ia+ioffs)=-v2(ia)
          endif
          indb(ibt)=0
        endif
4000    continue
      endif
      return
      end
c  **********************************
c  ** Routines involving CI and VB **
c  **********************************

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
      subroutine vbgenabdet2_cvb(
     >  idetavb,idetbvb,
     >  iconfs,nconf,nconfion,
     >  ndetvb,nel,noe,
     >  nalf,nbet,norb,
     >  xalf,xbet,mingrph,maxgrph,
     >  inewocc,iaccm)
      implicit real*8 (a-h,o-w,y-z),integer(x)
#include "malloc_cvb.fh"
      dimension idetavb(ndetvb),idetbvb(ndetvb)
      dimension iconfs(noe,nconf),nconfion(0:nel)
      dimension xalf(0:norb,0:nalf),xbet(0:norb,0:nbet)
      dimension mingrph(0:norb),maxgrph(0:norb)
      dimension inewocc(norb),iaccm(norb)
      logical debug
      data debug/.false./

c Set xalf and xbet for indget :
      do 100 iorb=0,norb
      mingrph(iorb)=max(iorb-norb+nalf,0)
      maxgrph(iorb)=min(iorb,nalf)
100   continue
      call weight_cvb(xalf,mingrph,maxgrph,nalf,norb)
      do 200 iorb=0,norb
      mingrph(iorb)=max(iorb-norb+nbet,0)
      maxgrph(iorb)=min(iorb,nbet)
200   continue
      call weight_cvb(xbet,mingrph,maxgrph,nbet,norb)

c  Transformation matrix VB structures  <->  determinants
      incrdet=0
      ioff_nconf=0
      do 1000 ion=0,nel/2
      nelsing=nel-2*ion
      nalfsing=nalf-ion
      nbetsing=nbet-ion
c  Skip if configurations are incompatible with Ms :
      if(nalfsing.lt.0.or.nbetsing.lt.0.or.nelsing.lt.0)goto 1001

c  Generate alpha/beta strings for singly occupied electrons :
      call icomb_cvb(nelsing,nalfsing,nstring)
      iastr = mstacki_cvb(nalfsing*nstring)
      ibstr = mstacki_cvb(nbetsing*nstring)
      call stringen_cvb(nelsing,nalfsing,iw(iastr),iw(ibstr),nstring)
      if(debug)then
        write(6,*)' ionicity=',ion,' nconf=',nconfion(ion)
        write(6,*)' check alpha strings :'
        do i=1,nstring
        write(6,*)i,' => ',(iw(ii+iastr-1+(i-1)*nalfsing),ii=1,nalfsing)
        enddo
        write(6,*)' check beta strings :'
        do i=1,nstring
        write(6,*)i,' => ',(iw(ii+ibstr-1+(i-1)*nbetsing),ii=1,nbetsing)
        enddo
      endif

      do 1100 iconf=ioff_nconf+1,ioff_nconf+nconfion(ion)
      call imove_cvb(iconfs(1,iconf),inewocc,norb)
      incr=0
      do 1200 iorb=1,norb
      if(inewocc(iorb).eq.1)then
        incr=incr+1
        iaccm(incr)=iorb
      endif
      inewocc(iorb)=max(0,inewocc(iorb)-1)
1200  continue

c  Spin string loop :
      do 1300 index=1,nstring

c  Alpha index in full string space ...
      do 1400 i=1,nalfsing
      iaocc=iaccm(iw(i+(index-1)*nalfsing+iastr-1))
      inewocc(iaocc)=inewocc(iaocc)+1
1400  continue
      iaind=indget_cvb(inewocc,nalf,norb,xalf)
      do 1500 i=1,nalfsing
      iaocc=iaccm(iw(i+(index-1)*nalfsing+iastr-1))
      inewocc(iaocc)=inewocc(iaocc)-1
1500  continue

c  Beta index in full string space ...
      do 1600 i=1,nbetsing
      ibocc=iaccm(iw(i+(index-1)*nbetsing+ibstr-1))
      inewocc(ibocc)=inewocc(ibocc)+1
1600  continue
      ibind=indget_cvb(inewocc,nbet,norb,xbet)
      do 1700 i=1,nbetsing
      ibocc=iaccm(iw(i+(index-1)*nbetsing+ibstr-1))
      inewocc(ibocc)=inewocc(ibocc)-1
1700  continue

      incrdet=incrdet+1
      idetavb(incrdet)=iaind
      idetbvb(incrdet)=ibind
1300  continue
1100  continue
      call mfreei_cvb(iastr)
1001  ioff_nconf=ioff_nconf+nconfion(ion)
1000  continue
      if(debug)then
        write(6,*)' idetavb='
        write(6,'(10i6)')idetavb
        write(6,*)' idetbvb='
        write(6,'(10i6)')idetbvb
      endif
      return
      end

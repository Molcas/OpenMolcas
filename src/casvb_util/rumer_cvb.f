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
      subroutine rumer_cvb(bikcof,
     > nel,nalf,nbet,ndet,ifns,kbasis,iprint,nswpdim,
     > minspn,maxspn,nkspn,
     > minswp,maxswp,nkswp,ioccswp,locca,lnocca,
     > xdet,xspin,iwork,
     > ialfs,ibets)
      implicit real*8 (a-h,o-w,y-z),integer(x)
      dimension bikcof(ndet,ifns)
      dimension minspn(0:nel), maxspn(0:nel), nkspn(0:nel)
      dimension minswp(0:2*nbet), maxswp(0:2*nbet)
      dimension nkswp(0:2*nbet)
      dimension ioccswp(nbet,nswpdim)
      dimension locca(nel), lnocca(nel)
      dimension xdet(0:nel,0:nalf),xspin((nel+1)*(nalf+1)),
     > iwork(nel)
      dimension ialfs(nalf),ibets(nbet)
      save one
      data one/1.0d0/

      call fzero(bikcof,ifns*ndet)

c Determinant weight array (index from alpha spin string):
      call weightfl_cvb(xdet,nalf,nel)
      if(ndet.ne.xdet(nel,nalf))then
        write(6,*) ' Discrepancy in NDET:',ndet,xdet(nel,nalf)
        call abend_cvb()
      endif

c
c Rumer spin functions
c
      if(iprint.ge.2.and.(kbasis.eq.3.or.kbasis.eq.4))
     >  write(6,6100) 2**nbet
      abphase=DBLE(1-2*mod(nbet*(nbet-1)/2,2))

c Prepare NKs for a<->b interchanges in (ab-ba) terms :
      nbet2=nbet+nbet
      do 1000 iorb=0,nbet2
      minswp(iorb)=iorb/2
      maxswp(iorb)=(iorb+1)/2
1000  continue
      call imove_cvb(maxswp,nkswp,nbet2+1)
      iswp=0
1100  iswp=iswp+1
      do 1200 ia=1,nbet
      ioccswp(ia,iswp)=nkswp(2*ia)-nkswp(2*ia-1)
1200  continue
      call loop_cvb(nbet2,nkswp,minswp,maxswp,*1100)

c Spin function weight arrays:
      do 2000 iorb=0,nel
      minspn(iorb)=max(iorb-nalf,0)
      maxspn(iorb)=min(iorb/2,nbet)
2000  continue
      call weight_cvb(xspin,minspn,maxspn,nbet,nel)
      if(ifns.ne.xspin((nel+1)*(nbet+1)).and.kbasis.ne.6)then
        write(6,*) ' Discrepancy in IFNS:',ifns,xspin((nel+1)*(nbet+1))
        call abend_cvb()
      endif
      call imove_cvb(maxspn,nkspn,nel+1)
      call occupy_cvb(nkspn,nel,locca,lnocca)

c Loop:
      index=1
c Determine pairings
2100  continue
      if(kbasis.eq.6.and.index.gt.ifns)goto 3300
      do 2200 ib=1,nbet
      ibets(ib)=locca(ib)
      do 2300 ia=nalf,1,-1
      ialfs(ib)=lnocca(ia)
      if(ialfs(ib).lt.ibets(ib))then
        do 2400 iachek=1,ib-1
        if(ialfs(iachek).eq.ialfs(ib))goto 2300
2400    continue
        goto 2500
      endif
2300  continue
2500  continue
2200  continue
      do 2600 ib=nbet+1,nalf
      do 2700 ia=1,nalf
      ialfs(ib)=lnocca(ia)
      do 2800 iachek=1,ib-1
      if(ialfs(iachek).eq.ialfs(ib))goto 2700
2800  continue
      goto 2600
2700  continue
2600  continue

      if(iprint.ge.2.and.(kbasis.eq.3.or.kbasis.eq.4))
     >  write(6,6200) index,(ialfs(ii),ibets(ii),ii=1,nbet)

      if(kbasis.eq.4)then
        bikvalue=one
      else
        bikvalue=abphase*party_cvb(ialfs,nalf)*party_cvb(ibets,nbet)
      endif

      do 3000 iswp=1,nswpdim
      call izero(iwork,nel)
      do 3100 ia=1,nalf
      iwork(ialfs(ia))=1
3100  continue
      do 3200 ia=1,nbet
      if(ioccswp(ia,iswp).ne.0)then
c  IA => alpha electron in position number 2 (= swap).
        iwork(ialfs(ia))=0
        iwork(ibets(ia))=1
      endif
3200  continue
      bikcof(indget_cvb(iwork,nalf,nel,xdet),index)=bikvalue
3000  continue
      call loind_cvb(nel,nbet,nkspn,minspn,maxspn,
     >                       locca,lnocca,index,xspin,*2100)
3300  continue

c Normalise
      scale=one/sqrt(dble(2**nbet))
      call dscal_(ndet*ifns,scale,bikcof,1)
      return
6100  format(/,' Number of determinants per structure:',i4)
6200  format(2x,i3,' ==> ',8(i2,'-',i2,3x))
      end

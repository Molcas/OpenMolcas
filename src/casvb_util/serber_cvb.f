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
      subroutine serber_cvb(bikcof,
     > nel,nalf,nbet,ndet,ifns,
     > minspn,maxspn,nkspn,locca,lnocca,xspin,
     > ialfs,ibets,ianti)
      implicit real*8 (a-h,o-w,y-z),integer(x)
      dimension bikcof(ndet,ifns)
      dimension minspn(0:nel), maxspn(0:nel), nkspn(0:nel)
      dimension locca(nel), lnocca(nel)
      dimension xspin((nel+1)*(nalf+1))
      dimension ialfs(nalf),ibets(nbet),ianti(ifns)
      dimension dum(1)

c
c Serber spin functions
c

c For each Rumer spin function we determine (minus) the "antisymmetry"
c number, as explained in [Spin Eigenfunctions, R. Pauncz, Sec. 5.5].

c Spin function weight arrays:
      do 2000 iorb=0,nel
      minspn(iorb)=max(iorb-nalf,0)
      maxspn(iorb)=min(iorb/2,nbet)
2000  continue
      call weight_cvb(xspin,minspn,maxspn,nbet,nel)
      if(ifns.ne.xspin((nel+1)*(nbet+1)))then
        write(6,*) ' Discrepancy in IFNS:',ifns,xspin((nel+1)*(nbet+1))
        call abend_cvb()
      endif
      call imove_cvb(maxspn,nkspn,nel+1)
      call occupy_cvb(nkspn,nel,locca,lnocca)
c Loop:
      index=1
c Determine pairings
2100  continue
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

      ianti(index)=0
      do 2600 ib=1,nbet
      if(mod(ialfs(ib),2).eq.1.and.ialfs(ib).eq.ibets(ib)-1)
     >  ianti(index)=ianti(index)-1
2600  continue
      call loind_cvb(nel,nbet,nkspn,minspn,maxspn,
     >                       locca,lnocca,index,xspin,*2100)

c Sort according to decreasing values of IANTI :
      nc=0
      do 4000 iantival=nbet,0,-1
      do 4001 k=1,ifns
      if(ianti(k).eq.-iantival)then
        nc=nc+1
        ianti(k)=nc
      endif
4001  continue
4000  continue

      do 5000 k=1,ifns
      if(ianti(k).ne.k)then
        do 5100 l=1,ifns
        if(ianti(l).eq.k)goto 5200
5100    continue
        write(6,*)' Error - swap function not found!',k,ianti(k)
        call abend_cvb()
5200    call dswap_(ndet,bikcof(1,k),1,bikcof(1,l),1)
        ianti(l)=ianti(k)
        ianti(k)=k
      endif
5000  continue

c Orthonormalize :
      call schmidtn_cvb(bikcof,ifns,dum,ndet,0)
      return
      end

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
      subroutine gaussj2_cvb(a,lrow,lcol,ibook,irows,ijs,oijs,n)
      implicit real*8 (a-h,o-z)
      dimension a(n,n)
      dimension lrow(n),lcol(n),ibook(n)
      dimension irows(n),ijs(2,n*n),oijs(n*n)
      save zero,one,thresh
      data zero/0.d0/,one/1.d0/,thresh/1.d-10/

c  initialize imx & jmx to suppress compiler warnings ...
      imx=0
      jmx=0

      nij=n*n
      do 100 i=1,n
      irows(i)=i
100   ibook(i)=0
      do 1000 imain=1,n
      amx=zero
      do 1100 i=1,n
      if(ibook(i).ne.1)then
        do 1200 j=1,n
        if(ibook(j).eq.0)then
          if(abs(a(i,j)).ge.amx)then
            amx=abs(a(i,j))
            imx=i
            jmx=j
          endif
        elseif(ibook(j).gt.1)then
          write(6,*)' Singular matrix in GAUSSJ !'
          call abend_cvb()
        endif
1200    continue
      endif
1100  continue
      ibook(jmx)=ibook(jmx)+1
      if(imx.ne.jmx)then
        do 1300 ii=1,n
        dum=a(imx,ii)
        a(imx,ii)=a(jmx,ii)
1300    a(jmx,ii)=dum
        idum=irows(imx)
        irows(imx)=irows(jmx)
        irows(jmx)=idum
      endif
      lrow(imain)=imx
      lcol(imain)=jmx
      if(abs(a(jmx,jmx)).lt.thresh)then
        ihad=imain-1
        goto 3000
      endif
      ijs(1,nij)=irows(jmx)
      ijs(2,nij)=irows(jmx)
      oijs(nij)=a(jmx,jmx)
      nij=nij-1
      oneovamx=one/a(jmx,jmx)
      a(jmx,jmx)=one
      do 1500 ii=1,n
1500  a(jmx,ii)=oneovamx*a(jmx,ii)
      do 1700 ii2=1,n
      if(ii2.ne.jmx)then
        dum=a(ii2,jmx)
        a(ii2,jmx)=zero
        do 1800 ii=1,n
1800    a(ii2,ii)=a(ii2,ii)-dum*a(jmx,ii)
        ijs(1,nij)=irows(ii2)
        ijs(2,nij)=irows(jmx)
        oijs(nij)=dum
        nij=nij-1
      endif
1700  continue
1000  continue
      do 2000 ii=n,1,-1
      if(lrow(ii).ne.lcol(ii))then
        do 2100 ii2=1,n
        dum=a(ii2,lrow(ii))
        a(ii2,lrow(ii))=a(ii2,lcol(ii))
2100    a(ii2,lcol(ii))=dum
      endif
2000  continue
      return
3000  continue
      do 3100 i=1,n
      do 3200 j=1,ihad
3200  if(lcol(j).eq.i)goto 3100
      ijs(1,nij)=irows(i)
      ijs(2,nij)=irows(i)
      oijs(nij)=zero
      nij=nij-1
      do 3300 j=1,n
      if(j.ne.i)then
        ijs(1,nij)=irows(j)
        ijs(2,nij)=irows(i)
        oijs(nij)=a(j,i)
        nij=nij-1
      endif
3300  continue
3100  continue
      return
      end

************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
       subroutine grc32y (mapda,mapdb,mapdc,mapia,mapib,mapic,
     & mvec,ssa,ssb,possc0,ix)
c
#include "ccsd1.fh"
c
       integer mapda(0:512,1:6)
       integer mapdb(0:512,1:6)
       integer mapdc(0:512,1:6)
c
       integer mapia(1:8,1:8,1:8)
       integer mapib(1:8,1:8,1:8)
       integer mapic(1:8,1:8,1:8)
c
       integer mvec(1:4096,1:7)
       integer possc0
       integer ssa,ssb
c
c     help variables
c
       integer nhelp1,nhelp2,nhelp4
       integer nhelp41,nhelp42
       integer ntest2
       integer sa1,sa2,sa3,sb1,sb2,sa23
       integer ia,ib,iy,ix
       integer possct
c
c     sctructure A(p,qr)*B(qr)=YC(p)
c
c1.0  prepare mapdc,mapic
c
       call grc0 (1,0,mapda(0,1),0,0,0,mmul(ssa,ssb),
     & possc0,possct,mapdc,mapic)
c
c1.1  define limitations - A p,q>r must be tested - ntest1
c     - B p>q must be tested - ntest2
c
c      if (mapda(0,6).eq.2) then
c      ntest1=1
c      else
c      ntest1=0
c      end if
c
       if (mapdb(0,6).eq.1) then
       ntest2=1
       else
       ntest2=0
       end if
c
c1.2  def symm states and test the limitations
c
       ix=1
       do 100 sb1=1,nsym
       sa2=sb1
c
       sb2=mmul(ssb,sb1)
       sa3=sb2
       sa23=mmul(sa2,sa3)
       if ((ntest2.eq.1).and.(sb1.lt.sb2)) then
c     Meggie out
       goto 100
       end if
c
       sa1=mmul(ssa,sa23)
c
c1.3  def mvec,mapdc and mapdi
c
       ia=mapia(sa1,sa2,1)
       ib=mapib(sb1,1,1)
       iy=mapic(1,1,1)
c
c     yes/no
       if ((mapda(ia,2).gt.0).and.(mapdb(ib,2).gt.0)) then
       nhelp1=1
       else
       goto 100
       end if
c
c     rowA
       nhelp2=dimm(mapda(0,1),sa1)
c
c     sum
       nhelp41=dimm(mapda(0,2),sa2)
       nhelp42=dimm(mapda(0,3),sa3)
       if ((ntest2.eq.1).and.(sa2.eq.sa3)) then
       nhelp4=nhelp41*(nhelp41-1)/2
       else
       nhelp4=nhelp41*nhelp42
       end if
c
       mvec(ix,1)=nhelp1
       mvec(ix,2)=mapda(ia,1)
       mvec(ix,3)=mapdb(ib,1)
       mvec(ix,4)=mapdc(iy,1)
       mvec(ix,5)=nhelp2
       mvec(ix,6)=nhelp4
       mvec(ix,7)=0
c
       ix=ix+1
c
 100    continue
       ix=ix-1
c
       return
       end

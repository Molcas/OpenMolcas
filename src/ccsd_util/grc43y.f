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
       subroutine grc43y (mapda,mapdb,mapdc,mapia,mapib,mapic,
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
       integer nhelp41,nhelp42,nhelp43
       integer ntest1,ntest2
       integer sa1,sa2,sa3,sa4,sb1,sb2,sb3,sa23,sa234,sb12
       integer nsymb2
       integer ia,ib,iy,ix
       integer possct
c
c     sctructure A(p,qrs)*B(qrs)=YC(p)
c
c1.0  prepare mapdc,mapic
c
       call grc0 (1,0,mapda(0,1),0,0,0,mmul(ssa,ssb),
     & possc0,possct,mapdc,mapic)
c
c1.1  define limitations - p,q>r,s must be tested - ntest1
c     p,q,r>s must be tested - ntest2
c
       if (mapda(0,6).eq.2) then
       ntest1=1
       else
       ntest1=0
       end if
c
       if (mapda(0,6).eq.3) then
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
       if (ntest1.eq.1) then
       nsymb2=sb1
       else
       nsymb2=nsym
       end if
c
       do 101 sb2=1,nsymb2
       sb12=mmul(sb1,sb2)
       sa3=sb2
       sa23=mmul(sa2,sa3)
c
       sb3=mmul(ssb,sb12)
       sa4=sb3
       sa234=mmul(sa23,sa4)
       if ((ntest2.eq.1).and.(sb2.lt.sb3)) then
c     Meggie out
       goto 101
       end if
c
       sa1=mmul(ssa,sa234)
       if ((ntest1.eq.1).and.(sb1.lt.sb2)) then
c     Meggie out
       goto 101
       end if
c
c1.3  def mvec,mapdc and mapdi
c
       ia=mapia(sa1,sa2,sa3)
       ib=mapib(sb1,sb2,1)
       iy=mapic(1,1,1)
c
c     yes/no
       if ((mapda(ia,2).gt.0).and.(mapdb(ib,2).gt.0)) then
       nhelp1=1
       else
       goto 101
       end if
c
c     rowA
       nhelp2=dimm(mapda(0,1),sa1)
c
c     sum
       nhelp41=dimm(mapda(0,2),sa2)
       nhelp42=dimm(mapda(0,3),sa3)
       nhelp43=dimm(mapda(0,4),sa4)
       if ((ntest1.eq.1).and.(sa2.eq.sa3)) then
       nhelp4=nhelp41*(nhelp41-1)*nhelp43/2
       else if ((ntest2.eq.1).and.(sa3.eq.sa4)) then
       nhelp4=nhelp41*nhelp42*(nhelp42-1)/2
       else
       nhelp4=nhelp41*nhelp42*nhelp43
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
 101    continue
 100    continue
       ix=ix-1
c
       return
       end

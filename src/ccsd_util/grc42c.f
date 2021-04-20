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
       subroutine grc42C (mapda,mapdb,mapdc,mapia,mapib,mapic,
     & mvec,ssa,ssb,pbar,possc0,ix)
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
       integer pbar,possc0
       integer ssa,ssb
c
c     help variables
c
       integer nhelp1,nhelp2,nhelp3,nhelp4
       integer nhelp21,nhelp22,nhelp23
       integer ntest1,ntest2
       integer sa1,sa2,sa3,sa4,sb1,sb2,sa12,sa123
       integer nsyma2,nsyma3
       integer ia,ib,ic,ix
       integer possct
c
c1*
c
       if (pbar.eq.1) then
c
c     not possible
c
       else if (pbar.eq.2) then
c
c     sctructure A(pq,rs)*B(rs)=C(pq)
c     implemented in grc42y
c
       else if (pbar.eq.3) then
c
c     sctructure A(pqr,s)*B(s,t)=C(pqr,t)
c
c1.1  define limitations -  p>q,r,s must be tested - ntest1
c     p,q>r,s must be tested - ntest2
c
c1.0  prepare mapdc,mapic
c
       call grc0 (4,mapda(0,6),mapda(0,1),mapda(0,2),mapda(0,3),mapdb(0,
     &            2),mmul(ssa,ssb),
     & possc0,possct,mapdc,mapic)
c
       if (mapda(0,6).eq.1) then
       ntest1=1
       else
       ntest1=0
       end if
c
       if (mapda(0,6).eq.2) then
       ntest2=1
       else
       ntest2=0
       end if
c
c1.2  def symm states and test the limitations
c
       ix=1
       do 100 sa1=1,nsym
       if (ntest1.eq.1) then
       nsyma2=sa1
       else
       nsyma2=nsym
       end if
c
       do 101 sa2=1,nsyma2
       sa12=mmul(sa1,sa2)
       if (ntest2.eq.1) then
       nsyma3=sa2
       else
       nsyma3=nsym
       end if
c
       do 102 sa3=1,nsyma3
       sa123=mmul(sa12,sa3)
c
       sa4=mmul(ssa,sa123)
       sb1=sa4
c
       sb2=mmul(ssb,sb1)
c
c1.3  def mvec,mapdc and mapdi
c
       ia=mapia(sa1,sa2,sa3)
       ib=mapib(sb1,1,1)
       ic=mapic(sa1,sa2,sa3)
c
c     yes/no
       if ((mapda(ia,2).gt.0).and.(mapdb(ib,2).gt.0)) then
       nhelp1=1
       else
       goto 102
       end if
c
c     rowA
       nhelp21=dimm(mapda(0,1),sa1)
       nhelp22=dimm(mapda(0,2),sa2)
       nhelp23=dimm(mapda(0,3),sa3)
       if ((ntest1.eq.1).and.(sa1.eq.sa2)) then
       nhelp2=nhelp21*(nhelp21-1)*nhelp23/2
       else if ((ntest2.eq.1).and.(sa2.eq.sa3)) then
       nhelp2=nhelp21*nhelp22*(nhelp22-1)/2
       else
       nhelp2=nhelp21*nhelp22*nhelp23
       end if
c
c     colB
       nhelp3=dimm(mapdb(0,2),sb2)
c
c     sum
       nhelp4=dimm(mapda(0,4),sa4)
c
       mvec(ix,1)=nhelp1
       mvec(ix,2)=mapda(ia,1)
       mvec(ix,3)=mapdb(ib,1)
       mvec(ix,4)=mapdc(ic,1)
       mvec(ix,5)=nhelp2
       mvec(ix,6)=nhelp4
       mvec(ix,7)=nhelp3
c
       ix=ix+1
c
 102    continue
 101    continue
 100    continue
c
c
       end if
       ix=ix-1
c
       return
       end

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
       subroutine grc34C (mapda,mapdb,mapdc,mapia,mapib,mapic,
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
       integer nhelp31,nhelp32
       integer nhelp41,nhelp42
       integer ntest1,ntest2
       integer sa1,sa2,sa3,sb1,sb2,sb3,sb4,sa12,sb12,sb123
       integer ia,ib,ic,ix
       integer possct
c
c1*
c
       if (pbar.eq.1) then
c
c     sctructure A(p,qr)*B(qr,st)=C(p,st)
c
c1.1  define limitations - q>r,s,t must be tested - ntest1
c     - q,r,s>t must be tested - ntest2
c
       if ((mapdb(0,6).eq.1).or.(mapdb(0,6).eq.4)) then
       ntest1=1
       else
       ntest1=0
       end if
c
       if ((mapdb(0,6).eq.3).or.(mapdb(0,6).eq.4)) then
       ntest2=1
       else
       ntest2=0
       end if
c
c1.0  prepare mapdc,mapic
c
       if (ntest2.eq.1) then
       nhelp1=2
       else
       nhelp1=0
       end if
c
       call grc0 (3,nhelp1,mapda(0,1),mapdb(0,3),mapdb(0,4),0,mmul(ssa,
     &            ssb),
     & possc0,possct,mapdc,mapic)
c
c1.2  def symm states and test the limitations
c
       ix=1
       do 130 sa1=1,nsym
c
       do 120 sa2=1,nsym
       sa12=mmul(sa1,sa2)
       sb1=sa2
c
       sa3=mmul(ssa,sa12)
       sb2=sa3
       sb12=mmul(sb1,sb2)
       if ((ntest1.eq.1).and.(sa2.lt.sa3)) then
c     Meggie out
       goto 120
       end if
c
       do 100 sb3=1,nsym
       sb123=mmul(sb12,sb3)
c
       sb4=mmul(ssb,sb123)
       if ((ntest2.eq.1).and.(sb3.lt.sb4)) then
c     Meggie out
       goto 100
       end if
c
c1.3  def mvec,mapdc and mapdi
c
       ia=mapia(sa1,sa2,1)
       ib=mapib(sb1,sb2,sb3)
       ic=mapic(sa1,sb3,1)
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
c     colB
       nhelp31=dimm(mapdb(0,3),sb3)
       nhelp32=dimm(mapdb(0,4),sb4)
       if ((ntest2.eq.1).and.(sb3.eq.sb4)) then
       nhelp3=nhelp31*(nhelp31-1)/2
       else
       nhelp3=nhelp31*nhelp32
       end if
c
c     sum
       nhelp41=dimm(mapda(0,2),sa2)
       nhelp42=dimm(mapda(0,3),sa3)
       if ((ntest1.eq.1).and.(sa2.eq.sa3)) then
       nhelp4=nhelp41*(nhelp41-1)/2
       else
       nhelp4=nhelp41*nhelp42
       end if
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
 100    continue
 120    continue
 130    continue
c
       else if (pbar.eq.3) then
c
c     sctructure A(pq,r)*B(r,stu)=C(pq,stu)
c     not used
c
       end if
       ix=ix-1
c
       return
       end

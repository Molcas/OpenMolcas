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
       subroutine grc22C (mapda,mapdb,mapdc,mapia,mapib,mapic,
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
       integer sa1,sa2,sb1,sb2
       integer ia,ib,ic,ix
       integer possct
c
c1*
c
       if (pbar.eq.1) then
c
c     sctructure A(p,q)*B(q,r)=C(p,r)
c
c1.0  prepare mapdc,mapic
c
       call grc0 (2,0,mapda(0,1),mapdb(0,2),0,0,mmul(ssa,ssb),
     & possc0,possct,mapdc,mapic)
c
c1.1  define limitations - no limitations
c
c1.2  def symm states and test the limitations
c
       ix=1
       do 100 sa1=1,nsym
c
       sa2=mmul(ssa,sa1)
       sb1=sa2
c
       sb2=mmul(ssb,sb1)
c
c1.3  def mvec,mapdc and mapdi
c
       ia=mapia(sa1,1,1)
       ib=mapib(sb1,1,1)
       ic=mapic(sa1,1,1)
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
       nhelp3=dimm(mapdb(0,2),sb2)
c
c     sum
       nhelp4=dimm(mapda(0,2),sa2)
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
c
c
       end if
       ix=ix-1
c
       return
       end

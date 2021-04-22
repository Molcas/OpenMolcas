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
* Copyright (C) 2006, Pavel Neogrady                                   *
************************************************************************
       subroutine multstack (wrk,wrksize,
     &                       mapda,mapdb,mapdc,mapia,mapib,mapic,
     &                       ssa,ssb,possc0,bsize)
c
c        This is a special routine for multiplying of the:
c        C(ij,Bp) = A(ij,cd) . B(cd,Bp)
c        where Bp is a limited (partial) sumation over
c        b index (#b - bsize), namely those, hich are stacked
c
c        This routine is used only in stacking in sumoverab
c        process and is a modification of grc42y routine.
c        Type of index Bp is registered as for standard b, but
c        all lengths of blocks are calculated with
c        bsize, instead of dimm(typb,symb) and the symmetry
c       of b is ignored, B and C are treated as 2index
c
c        P.N. 17.02.06
c
c
#include "ccsd1.fh"
#include "wrk.fh"
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
       integer possc0,bsize
       integer ssa,ssb
c
c     help variables
c
       integer nhelp1,nhelp2,nhelp3,nhelp4
       integer nhelp21,nhelp22,nhelp41,nhelp42
       integer ntest1,ntest2
       integer sa1,sa2,sa3,sa4,sb1,sb2,sa34,sa134
CLD    integer sa1,sa2,sa3,sa4,sb1,sb2,sb3,sa12,sa34,sa134,sb12
CLD    integer nsyma2
       integer ia,ib,ix,iy
CLD    integer ia,ib,ic,ix,iy
       integer possct
c
c1*
c
c     sctructure A(pq,rs)*B(rs,t)=C(pq,t)
c     sctructure A(pq,rs)*B(rs)=YC(pq)
c
c1.1  define limitations -  p>q,r,s must be tested - ntest1
c     p,q,r>s must be tested - ntest2
c
       if ((mapda(0,6).eq.1).or.(mapda(0,6).eq.4)) then
       ntest1=1
       else
       ntest1=0
       end if
c
       if ((mapda(0,6).eq.3).or.(mapda(0,6).eq.4)) then
       ntest2=1
       else
       ntest2=0
       end if
c
c1.0  prepare mapdc,mapic
c
       call grc0stack (bsize,ntest1,mapda(0,1),mapda(0,2),mapdb(0,3),
     &                 0,mmul(ssa,ssb),possc0,possct,mapdc,mapic)

c
c
c1.2  def symm states and test the limitations
c
       ix=1
       do 100 sb1=1,nsym
       sa3=sb1
c
       sb2=mmul(ssb,sb1)
       sa4=sb2
       sa34=mmul(sa3,sa4)
       if ((ntest2.eq.1).and.(sb1.lt.sb2)) then
c     Meggie out
       goto 100
       end if
c
       do 50 sa1=1,nsym
       sa134=mmul(sa1,sa34)
c
       sa2=mmul(ssa,sa134)
       if ((ntest1.eq.1).and.(sa1.lt.sa2)) then
c     Meggie out
       goto 50
       end if
c
c1.3  def mvec,mapdc and mapdi
c
       ia=mapia(sa1,sa2,sa3)
       ib=mapib(sb1,1,1)
       iy=mapic(sa1,1,1)
c
c     yes/no
       if ((mapda(ia,2).gt.0).and.(mapdb(ib,2).gt.0)) then
       nhelp1=1
       else
       goto 50
       end if
c
c     rowA
       nhelp21=dimm(mapda(0,1),sa1)
       nhelp22=dimm(mapda(0,2),sa2)
       if ((ntest1.eq.1).and.(sa1.eq.sa2)) then
       nhelp2=nhelp21*(nhelp21-1)/2
       else
       nhelp2=nhelp21*nhelp22
       end if
c
c     sum
       nhelp41=dimm(mapda(0,3),sa3)
       nhelp42=dimm(mapda(0,4),sa4)
       if ((ntest2.eq.1).and.(sa3.eq.sa4)) then
       nhelp4=nhelp41*(nhelp41-1)/2
       else
       nhelp4=nhelp41*nhelp42
       end if
c
c     colBp
       nhelp3=bsize
c
       mvec(ix,1)=nhelp1
       mvec(ix,2)=mapda(ia,1)
       mvec(ix,3)=mapdb(ib,1)
       mvec(ix,4)=mapdc(iy,1)
       mvec(ix,5)=nhelp2
       mvec(ix,6)=nhelp4
       mvec(ix,7)=nhelp3
c
       ix=ix+1
c
 50     continue
 100    continue
       ix=ix-1
c
c
c*        multiplying
c
        call multc0 (wrk,wrksize,
     &               mvec,ix,mapdc,1)
c
       return
       end

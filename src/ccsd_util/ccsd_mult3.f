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

c
c     this file contains following routines:
c     grc34c
c     grc33c
c     grc32c
c     grc24c
c     grc23c
c     grc22c
c
c     -------------------------------------------------------
c
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
c
c     -----------------------------
c
       subroutine grc33C (mapda,mapdb,mapdc,mapia,mapib,mapic,
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
       integer nhelp21,nhelp22,nhelp31,nhelp32
       integer nhelp41,nhelp42
       integer ntest1,ntest2
       integer sa1,sa2,sa3,sb1,sb2,sb3,sa12,sb12
       integer nsyma2
       integer ia,ib,ic,ix
       integer possct
c
c1*
c
       if (pbar.eq.1) then
c
c     sctructure A(p,qr)*B(qr,s)=C(p,s)
c
c1.0  prepare mapdc,mapic
c
       call grc0 (2,0,mapda(0,1),mapdb(0,3),0,0,mmul(ssa,ssb),
     & possc0,possct,mapdc,mapic)
c
c1.1  define limitations - p,q>r must be tested - ntest1
c
       if (mapda(0,6).eq.2) then
       ntest1=1
       else
       ntest1=0
       end if
c
c1.2  def symm states and test the limitations
c
       ix=1
       do 100 sa1=1,nsym
c
       do 101 sa2=1,nsym
       sa12=mmul(sa1,sa2)
       sb1=sa2
c
       sa3=mmul(ssa,sa12)
       if ((ntest1.eq.1).and.(sa2.lt.sa3)) then
c     Meggie out
       goto 101
       end if
       sb2=sa3
       sb12=mmul(sb1,sb2)
c
       sb3=mmul(ssb,sb12)
c
c1.3  def mvec,mapdc and mapdi
c
       ia=mapia(sa1,sa2,1)
       ib=mapib(sb1,sb2,1)
       ic=mapic(sa1,sb3,1)
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
c     colB
       nhelp3=dimm(mapdb(0,3),sb3)
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
 101    continue
 100    continue
c
       else if (pbar.eq.2) then
c
c     sctructure A(pq,r)*B(r,st)=C(pq,st)
c
c2.1  define limitations - A p>q,r must be tested - ntest1
c     B p,q>r must be tested - ntest2
c
       if (mapda(0,6).eq.1) then
       ntest1=1
       else
       ntest1=0
       end if
c
       if (mapdb(0,6).eq.2) then
       ntest2=1
       else
       ntest2=0
       end if
c
c2.0  prepare mapdc,mapic
c
       if ((ntest1.eq.1).and.(ntest2.eq.1)) then
       nhelp1=4
       else if (ntest1.eq.1) then
       nhelp1=1
       else if (ntest2.eq.1) then
       nhelp1=3
       else
       nhelp1=0
       end if
c
       call grc0 (4,nhelp1,mapda(0,1),mapda(0,2),mapdb(0,2),mapdb(0,3),
     &            mmul(ssa,ssb),
     & possc0,possct,mapdc,mapic)
c
c2.2  def symm states and test the limitations
c
       ix=1
       do 240 sa1=1,nsym
       if (ntest1.eq.1) then
       nsyma2=sa1
       else
       nsyma2=nsym
       end if
c
       do 230 sa2=1,nsyma2
       sa12=mmul(sa1,sa2)
c
       sa3=mmul(ssa,sa12)
       sb1=sa3
c
       do 220 sb2=1,nsym
       sb12=mmul(sb1,sb2)
c
       sb3=mmul(ssb,sb12)
       if ((ntest2.eq.1).and.(sb2.lt.sb3)) then
c     Meggie out
       goto 220
       end if
c
c2.3  def mvec,mapdc and mapdi
c
       ia=mapia(sa1,sa2,sa3)
       ib=mapib(sb1,sb2,sb3)
       ic=mapic(sa1,sa2,sb2)
c
c     yes/no
       if ((mapda(ia,2).gt.0).and.(mapdb(ib,2).gt.0)) then
       nhelp1=1
       else
       goto 220
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
c     colB
       nhelp31=dimm(mapdb(0,2),sb2)
       nhelp32=dimm(mapdb(0,3),sb3)
       if ((ntest2.eq.1).and.(sb2.eq.sb3)) then
       nhelp3=nhelp31*(nhelp31-1)/2
       else
       nhelp3=nhelp31*nhelp32
       end if
c
c     sum
       nhelp4=dimm(mapda(0,3),sa3)
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
 220    continue
 230    continue
 240    continue
c
c
       end if
       ix=ix-1
c
       return
       end
c
c     -----------------------------
c
       subroutine grc32C (mapda,mapdb,mapdc,mapia,mapib,mapic,
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
       integer nhelp21,nhelp22
       integer ntest1
       integer sa1,sa2,sa3,sb1,sb2,sa12
       integer nsyma2
       integer ia,ib,ic,ix
       integer possct
c
c1*
c
       if (pbar.eq.1) then
c
c     sctructure A(p,qr)*B(qr)=y(p)
c     implemented in grc32y
c
       else if (pbar.eq.2) then
c
c     sctructure A(pq,r)*B(r,s)=C(pq,s)
c
c1.0  prepare mapdc,mapic
c
       call grc0 (3,mapdb(0,6),mapda(0,1),mapda(0,2),mapdb(0,2),0,
     &            mmul(ssa,ssb),
     & possc0,possct,mapdc,mapic)
c
c1.1  define limitations - p>q,r must be tested - ntest1
c
       if (mapdb(0,6).eq.1) then
       ntest1=1
       else
       ntest1=0
       end if
c
c1.2  def symm states and test the limitations
c
       ix=1
       do 130 sa1=1,nsym
       if (ntest1.eq.1) then
       nsyma2=sa1
       else
       nsyma2=nsym
       end if
c
       do 120 sa2=1,nsyma2
       sa12=mmul(sa1,sa2)
c
       sa3=mmul(ssa,sa12)
       sb1=sa3
c
       sb2=mmul(ssb,sb1)
c
c1.3  def mvec,mapdc and mapdi
c
       ia=mapia(sa1,sa2,1)
       ib=mapib(sb1,1,1)
       ic=mapic(sa1,sa2,1)
c
c     yes/no
       if ((mapda(ia,2).gt.0).and.(mapdb(ib,2).gt.0)) then
       nhelp1=1
       else
       goto 120
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
c     colB
       nhelp3=dimm(mapdb(0,2),sb2)
c
c     sum
       nhelp4=dimm(mapda(0,3),sa3)
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
 120    continue
 130    continue
c
c
       end if
       ix=ix-1
c
       return
       end
c
c     -----------------------------
c
       subroutine grc24C (mapda,mapdb,mapdc,mapia,mapib,mapic,
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
       integer nhelp31,nhelp32,nhelp33
       integer ntest1,ntest2
       integer sa1,sa2,sb1,sb2,sb3,sb4,sb12,sb123
       integer nsymb3
       integer ia,ib,ic,ix
       integer possct
c
c1*
c
       if (pbar.eq.1) then
c
c     sctructure A(p,q)*B(q,rst)=C(p,rst)
c
c1.0  prepare mapdc,mapic
c
       call grc0 (4,mapdb(0,6),mapda(0,1),mapdb(0,2),mapdb(0,3),mapdb(0,
     &            4),mmul(ssa,ssb),
     & possc0,possct,mapdc,mapic)
c
c1.1  define limitations - q,r>s,t must be tested - ntest1
c     - q,r,s>t must be tested - ntest2
c
       if (mapdb(0,6).eq.2) then
       ntest1=1
       else
       ntest1=0
       end if
c
       if (mapdb(0,6).eq.3) then
       ntest2=1
       else
       ntest2=0
       end if
c
c1.2  def symm states and test the limitations
c
       ix=1
       do 100 sa1=1,nsym
c
       sa2=mmul(ssa,sa1)
       sb1=sa2
c
       do 101 sb2=1,nsym
       sb12=mmul(sb1,sb2)
       if (ntest1.eq.1) then
       nsymb3=sb2
       else
       nsymb3=nsym
       end if
c
       do 102 sb3=1,nsymb3
       sb123=mmul(sb12,sb3)
c
       sb4=mmul(ssb,sb123)
       if ((ntest2.eq.1).and.(sb3.lt.sb4)) then
c     Meggie out
       goto 102
       end if
c
c1.3  def mvec,mapdc and mapdi
c
       ia=mapia(sa1,1,1)
       ib=mapib(sb1,sb2,sb3)
       ic=mapic(sa1,sb2,sb3)
c
c     yes/no
       if ((mapda(ia,2).gt.0).and.(mapdb(ib,2).gt.0)) then
       nhelp1=1
       else
       goto 102
       end if
c
c     rowA
       nhelp2=dimm(mapda(0,1),sa1)
c
c     colB
       nhelp31=dimm(mapdb(0,2),sb2)
       nhelp32=dimm(mapdb(0,3),sb3)
       nhelp33=dimm(mapdb(0,4),sb4)
       if ((ntest1.eq.1).and.(sb2.eq.sb3)) then
       nhelp3=nhelp31*(nhelp31-1)*nhelp33/2
       else if ((ntest2.eq.1).and.(sb3.eq.sb4)) then
       nhelp3=nhelp31*nhelp32*(nhelp32-1)/2
       else
       nhelp3=nhelp31*nhelp32*nhelp33
       end if
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
c
c     -----------------------------
c
       subroutine grc23C (mapda,mapdb,mapdc,mapia,mapib,mapic,
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
       integer ntest1
       integer sa1,sa2,sb1,sb2,sb3,sb12
       integer ia,ib,ic,ix
       integer possct
c
c1*
c
       if (pbar.eq.1) then
c
c     sctructure A(p,q)*B(q,rs)=C(p,rs)
c
c1.0  prepare mapdc,mapic
c
       call grc0 (3,mapdb(0,6),mapda(0,1),mapdb(0,2),mapdb(0,3),0,
     &            mmul(ssa,ssb),
     & possc0,possct,mapdc,mapic)
c
c1.1  define limitations - q,r>s must be tested - ntest1
c
       if (mapdb(0,6).eq.2) then
       ntest1=1
       else
       ntest1=0
       end if
c
c1.2  def symm states and test the limitations
c
       ix=1
       do 100 sa1=1,nsym
c
       sa2=mmul(ssa,sa1)
       sb1=sa2
c
       do 101 sb2=1,nsym
       sb12=mmul(sb1,sb2)
c
       sb3=mmul(ssb,sb12)
       if ((ntest1.eq.1).and.(sb2.lt.sb3)) then
c     Meggie out
       goto 101
       end if
c
c1.3  def mvec,mapdc and mapdi
c
       ia=mapia(sa1,1,1)
       ib=mapib(sb1,sb2,1)
       ic=mapic(sa1,sb2,1)
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
c     colB
       nhelp31=dimm(mapdb(0,2),sb2)
       nhelp32=dimm(mapdb(0,3),sb3)
       if ((ntest1.eq.1).and.(sb2.eq.sb3)) then
       nhelp3=nhelp31*(nhelp31-1)/2
       else
       nhelp3=nhelp31*nhelp32
       end if
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
 101    continue
 100    continue
c
c
       end if
       ix=ix-1
c
       return
       end
c
c     -----------------------------
c
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
c
c     -----------------------------
c

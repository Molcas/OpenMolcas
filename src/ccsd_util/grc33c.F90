!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
       subroutine grc33C (mapda,mapdb,mapdc,mapia,mapib,mapic,          &
     & mvec,ssa,ssb,pbar,possc0,ix)
!
#include "ccsd1.fh"
!
       integer mapda(0:512,1:6)
       integer mapdb(0:512,1:6)
       integer mapdc(0:512,1:6)
!
       integer mapia(1:8,1:8,1:8)
       integer mapib(1:8,1:8,1:8)
       integer mapic(1:8,1:8,1:8)
!
       integer mvec(1:4096,1:7)
       integer pbar,possc0
       integer ssa,ssb
!
!     help variables
!
       integer nhelp1,nhelp2,nhelp3,nhelp4
       integer nhelp21,nhelp22,nhelp31,nhelp32
       integer nhelp41,nhelp42
       integer ntest1,ntest2
       integer sa1,sa2,sa3,sb1,sb2,sb3,sa12,sb12
       integer nsyma2
       integer ia,ib,ic,ix
       integer possct
!
!1*
!
       if (pbar.eq.1) then
!
!     sctructure A(p,qr)*B(qr,s)=C(p,s)
!
!1.0  prepare mapdc,mapic
!
       call grc0 (2,0,mapda(0,1),mapdb(0,3),0,0,mmul(ssa,ssb),          &
     & possc0,possct,mapdc,mapic)
!
!1.1  define limitations - p,q>r must be tested - ntest1
!
       if (mapda(0,6).eq.2) then
       ntest1=1
       else
       ntest1=0
       end if
!
!1.2  def symm states and test the limitations
!
       ix=1
       do 100 sa1=1,nsym
!
       do 101 sa2=1,nsym
       sa12=mmul(sa1,sa2)
       sb1=sa2
!
       sa3=mmul(ssa,sa12)
       if ((ntest1.eq.1).and.(sa2.lt.sa3)) then
!     Meggie out
       goto 101
       end if
       sb2=sa3
       sb12=mmul(sb1,sb2)
!
       sb3=mmul(ssb,sb12)
!
!1.3  def mvec,mapdc and mapdi
!
       ia=mapia(sa1,sa2,1)
       ib=mapib(sb1,sb2,1)
       ic=mapic(sa1,sb3,1)
!
!     yes/no
       if ((mapda(ia,2).gt.0).and.(mapdb(ib,2).gt.0)) then
       nhelp1=1
       else
       goto 101
       end if
!
!     rowA
       nhelp2=dimm(mapda(0,1),sa1)
!
!     colB
       nhelp3=dimm(mapdb(0,3),sb3)
!
!     sum
       nhelp41=dimm(mapda(0,2),sa2)
       nhelp42=dimm(mapda(0,3),sa3)
       if ((ntest1.eq.1).and.(sa2.eq.sa3)) then
       nhelp4=nhelp41*(nhelp41-1)/2
       else
       nhelp4=nhelp41*nhelp42
       end if
!
       mvec(ix,1)=nhelp1
       mvec(ix,2)=mapda(ia,1)
       mvec(ix,3)=mapdb(ib,1)
       mvec(ix,4)=mapdc(ic,1)
       mvec(ix,5)=nhelp2
       mvec(ix,6)=nhelp4
       mvec(ix,7)=nhelp3
!
       ix=ix+1
!
 101    continue
 100    continue
!
       else if (pbar.eq.2) then
!
!     sctructure A(pq,r)*B(r,st)=C(pq,st)
!
!2.1  define limitations - A p>q,r must be tested - ntest1
!     B p,q>r must be tested - ntest2
!
       if (mapda(0,6).eq.1) then
       ntest1=1
       else
       ntest1=0
       end if
!
       if (mapdb(0,6).eq.2) then
       ntest2=1
       else
       ntest2=0
       end if
!
!2.0  prepare mapdc,mapic
!
       if ((ntest1.eq.1).and.(ntest2.eq.1)) then
       nhelp1=4
       else if (ntest1.eq.1) then
       nhelp1=1
       else if (ntest2.eq.1) then
       nhelp1=3
       else
       nhelp1=0
       end if
!
       call grc0 (4,nhelp1,mapda(0,1),mapda(0,2),mapdb(0,2),mapdb(0,3), &
     &            mmul(ssa,ssb),                                        &
     & possc0,possct,mapdc,mapic)
!
!2.2  def symm states and test the limitations
!
       ix=1
       do 240 sa1=1,nsym
       if (ntest1.eq.1) then
       nsyma2=sa1
       else
       nsyma2=nsym
       end if
!
       do 230 sa2=1,nsyma2
       sa12=mmul(sa1,sa2)
!
       sa3=mmul(ssa,sa12)
       sb1=sa3
!
       do 220 sb2=1,nsym
       sb12=mmul(sb1,sb2)
!
       sb3=mmul(ssb,sb12)
       if ((ntest2.eq.1).and.(sb2.lt.sb3)) then
!     Meggie out
       goto 220
       end if
!
!2.3  def mvec,mapdc and mapdi
!
       ia=mapia(sa1,sa2,sa3)
       ib=mapib(sb1,sb2,sb3)
       ic=mapic(sa1,sa2,sb2)
!
!     yes/no
       if ((mapda(ia,2).gt.0).and.(mapdb(ib,2).gt.0)) then
       nhelp1=1
       else
       goto 220
       end if
!
!     rowA
       nhelp21=dimm(mapda(0,1),sa1)
       nhelp22=dimm(mapda(0,2),sa2)
       if ((ntest1.eq.1).and.(sa1.eq.sa2)) then
       nhelp2=nhelp21*(nhelp21-1)/2
       else
       nhelp2=nhelp21*nhelp22
       end if
!
!     colB
       nhelp31=dimm(mapdb(0,2),sb2)
       nhelp32=dimm(mapdb(0,3),sb3)
       if ((ntest2.eq.1).and.(sb2.eq.sb3)) then
       nhelp3=nhelp31*(nhelp31-1)/2
       else
       nhelp3=nhelp31*nhelp32
       end if
!
!     sum
       nhelp4=dimm(mapda(0,3),sa3)
!
       mvec(ix,1)=nhelp1
       mvec(ix,2)=mapda(ia,1)
       mvec(ix,3)=mapdb(ib,1)
       mvec(ix,4)=mapdc(ic,1)
       mvec(ix,5)=nhelp2
       mvec(ix,6)=nhelp4
       mvec(ix,7)=nhelp3
!
       ix=ix+1
!
 220    continue
 230    continue
 240    continue
!
!
       end if
       ix=ix-1
!
       return
       end

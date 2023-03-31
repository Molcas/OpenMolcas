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
       subroutine unpackab2 (wrk,wrksize,                               &
     & mapdn,mapin,mapdr1,mapir1,mapdr2,mapir2,                         &
     &                       mapdr3,mapir3,                             &
     & mapdr4,mapir4,mapdr5,mapir5,mapdr6,mapir6,ssn,key,aeqb)
!
!     mapdn - direct map matrix corresponding to N (Input)
!     mapin - direct map matrix corresponding to N (Input)
!     mapdr1- inverse map matrix corresponding to R1 (Input)
!     mapir1- inverse map matrix corresponding to R1 (Input)
!     mapdr2- inverse map matrix corresponding to R2 (Input)
!     mapir2- inverse map matrix corresponding to R2 (Input)
!     mapdr3- inverse map matrix corresponding to R3 (Input)
!     mapir3- inverse map matrix corresponding to R3 (Input)
!     mapdr4- inverse map matrix corresponding to R4 (Input)
!     mapir4- inverse map matrix corresponding to R4 (Input)
!     mapdr5- inverse map matrix corresponding to R4 (Input)
!     mapir5- inverse map matrix corresponding to R4 (Input)
!     mapdr6- inverse map matrix corresponding to R4 (Input)
!     mapir6- inverse map matrix corresponding to R4 (Input)
!     ssn   - overall symmtry state of N (Input)
!     key   - define state of _a,_b  key   a, b (Input)
!     1    S  S
!     2    V  S
!     3    S  V
!     4    V  V
!     aqeb  - aqeb=0 if a.ne.b; =1 if a=b (Input)
!
!     this routine unpack: N _a_b (p,q) = <ab|pq>              key   aeqb
!     N _a_b (p,q) a>=b-aa  -> R1 _a,_b(j,e)aaaa = <ab||je>     4     0
!     -> R2 _a,_b(j,e)bbbb = <ab||je>    1-4    0
!     -> R3 _a,_b(j,e)abab = <ab||je>    3,4   0,1
!     -> R4 _b,_a(j,e)abab = <ba||je>    2,4   0,1
!     -> R5 _a,_b(j,e)abba = <ab||je>    3,4   0,1
!     -> R6 _b,_a(j,e)abba = <ba||je>    2,4   0,1
!
#include "ccsd1.fh"
#include "wrk.fh"
!
       integer mapdn(0:512,1:6)
       integer mapdr1(0:512,1:6)
       integer mapdr2(0:512,1:6)
       integer mapdr3(0:512,1:6)
       integer mapdr4(0:512,1:6)
       integer mapdr5(0:512,1:6)
       integer mapdr6(0:512,1:6)
!
       integer mapin(1:8,1:8,1:8)
       integer mapir1(1:8,1:8,1:8)
       integer mapir2(1:8,1:8,1:8)
       integer mapir3(1:8,1:8,1:8)
       integer mapir4(1:8,1:8,1:8)
       integer mapir5(1:8,1:8,1:8)
       integer mapir6(1:8,1:8,1:8)
!
       integer ssn,key,aeqb
!
!     help variables
!
       integer symp,symq,dimp,dimq,dime,dimj
       integer in,inm,inp,ir1,ir2,ir3,ir4,ir5,ir6
       integer possn,possnp,possnm,possr1,possr2,possr3,possr4,possr5,  &
     & possr6
       integer lengthn
!
       do 1000 symp=1,nsym
       symq=mmul(ssn,symp)
!
       in=mapin(symp,1,1)
       lengthn=mapdn(in,2)
!
       if (lengthn.eq.0) goto 1000
!
       if (symp.eq.symq) then
!     symp=symq
!
       possn=mapdn(in,1)
       dimp=dimm(5,symp)
!
!I.1  def R1 - _a_b(j,e)aaaa
!
       if ((key.eq.4).and.(aeqb.eq.0)) then
       ir1=mapir1(symp,1,1)
       possr1=mapdr1(ir1,1)
       dimj=noa(symp)
       dime=nva(symq)
       if (mapdr1(ir1,2).gt.0) then
       call unpckhelp5 (wrk(possn),wrk(possr1),dimp,dimj,dime,          &
     & 0,noa(symp),noa(symp),nva(symp))
       end if
       end if
!
!I.2  def R2 - _a_b(j,e)bbbb
!
       if (aeqb.eq.0) then
       ir2=mapir2(symp,1,1)
       possr2=mapdr2(ir2,1)
       dimj=nob(symp)
       dime=nvb(symq)
       if (mapdr2(ir2,2).gt.0) then
       call unpckhelp5 (wrk(possn),wrk(possr2),dimp,dimj,dime,          &
     & 0,nob(symp),nob(symp),nvb(symp))
       end if
       end if
!
!I.3  def R3 - _a_b(j,e)abba
!
       if ((key.eq.2).or.(key.eq.4)) then
       ir3=mapir3(symp,1,1)
       possr3=mapdr3(ir3,1)
       dimj=nob(symp)
       dime=nva(symq)
       if (mapdr3(ir3,2).gt.0) then
       call unpckhelp7 (wrk(possn),wrk(possr3),dimp,dimp,dimj,          &
     & dime,0,nob(symp),noa(symq),nva(symq))
       end if
       end if
!
!I.4  def R4 - _b_a(j,e)abba
!
       if (((key.eq.3).or.(key.eq.4)).and.(aeqb.eq.0)) then
       ir4=mapir4(symp,1,1)
       possr4=mapdr4(ir4,1)
       dimj=nob(symp)
       dime=nva(symq)
       if (mapdr4(ir4,2).gt.0) then
       call unpckhelp6 (wrk(possn),wrk(possr4),dimp,dimp,dimj,          &
     & dime,0,nob(symp),noa(symq),nva(symq))
       end if
       end if
!
!I.5  def R5 - _a_b (j,e)abab
!
       if ((key.eq.2).or.(key.eq.4)) then
       ir5=mapir5(symp,1,1)
       possr5=mapdr5(ir5,1)
       dimj=noa(symp)
       dime=nvb(symq)
       if (mapdr5(ir5,2).gt.0) then
       call unpckhelp3 (wrk(possn),wrk(possr5),dimp,dimp,dimj,          &
     & dime,0,noa(symp),nob(symq),nvb(symq))
       end if
       end if
!
!I.6  def R6 - _b_a(j,e)abab
!
       if (((key.eq.3).or.(key.eq.4)).and.(aeqb.eq.0)) then
       ir6=mapir6(symp,1,1)
       possr6=mapdr6(ir6,1)
       dimj=noa(symp)
       dime=nvb(symq)
       if (mapdr6(ir6,2).gt.0) then
       call unpckhelp4 (wrk(possn),wrk(possr6),dimp,dimp,dimj,          &
     & dime,0,noa(symp),nob(symq),nvb(symq))
       end if
       end if
!
       else
!     symp.ne.symq
!
       inp=mapin(symp,1,1)
       inm=mapin(symq,1,1)
       possnp=mapdn(inp,1)
       possnm=mapdn(inm,1)
       dimp=dimm(5,symp)
       dimq=dimm(5,symq)
!
!II.1 def R1 - _a_b(j,e)aaaa
!
       if (key.eq.4) then
       ir1=mapir1(symp,1,1)
       possr1=mapdr1(ir1,1)
       dimj=noa(symp)
       dime=nva(symq)
       if (mapdr1(ir1,2).gt.0) then
       call unpckhelp2 (wrk(possnp),wrk(possnm),wrk(possr1),dimp,dimq,  &
     &                  dimj,dime,                                      &
     & 0,noa(symp),noa(symq),nva(symq))
       end if
       end if
!
!II.2 def R2 - _a_b(j,e)bbbb
!
       ir2=mapir2(symp,1,1)
       possr2=mapdr2(ir2,1)
       dimj=nob(symp)
       dime=nvb(symq)
       if (mapdr2(ir2,2).gt.0) then
       call unpckhelp2 (wrk(possnp),wrk(possnm),wrk(possr2),dimp,dimq,  &
     &                  dimj,dime,                                      &
     & 0,nob(symp),nob(symq),nvb(symq))
       end if
!
!II.3 def R3 - _a_b(j,e)abba
!
       if ((key.eq.2).or.(key.eq.4)) then
       ir3=mapir3(symq,1,1)
       possr3=mapdr3(ir3,1)
       dimj=nob(symq)
       dime=nva(symp)
       if (mapdr3(ir3,2).gt.0) then
       call unpckhelp7 (wrk(possnp),wrk(possr3),dimp,dimq,dimj,dime,    &
     & 0,nob(symq),noa(symp),nva(symp))
       end if
       end if
!
!II.4 def R4 - _b_a(j,e)abba
!
       if ((key.eq.3).or.(key.eq.4)) then
       ir4=mapir4(symp,1,1)
       possr4=mapdr4(ir4,1)
       dimj=nob(symp)
       dime=nva(symq)
       if (mapdr4(ir4,2).gt.0) then
       call unpckhelp6 (wrk(possnp),wrk(possr4),dimp,dimq,dimj,dime,    &
     & 0,nob(symp),noa(symq),nva(symq))
       end if
       end if
!
!II.5 def R5 - _a_b(j,e)abab
!
       if ((key.eq.2).or.(key.eq.4)) then
       ir5=mapir5(symp,1,1)
       possr5=mapdr5(ir5,1)
       dimj=noa(symp)
       dime=nvb(symq)
       if (mapdr5(ir5,2).gt.0) then
       call unpckhelp3 (wrk(possnp),wrk(possr5),dimp,dimq,dimj,dime,    &
     & 0,noa(symp),nob(symq),nvb(symq))
       end if
       end if
!
!II.6 def R6 - _b_a(j,e)abab
!
       if ((key.eq.3).or.(key.eq.4)) then
       ir6=mapir6(symq,1,1)
       possr6=mapdr6(ir6,1)
       dimj=noa(symq)
       dime=nvb(symp)
       if (mapdr6(ir6,2).gt.0) then
       call unpckhelp4 (wrk(possnp),wrk(possr6),dimp,dimq,dimj,dime,    &
     & 0,noa(symq),nob(symp),nvb(symp))
       end if
       end if
!
       end if
!
 1000   continue
!
       return
       end

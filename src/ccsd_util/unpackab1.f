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
       subroutine unpackab1 (wrk,wrksize,
     & mapdn,mapin,mapdr1,mapir1,
     & mapdr2,mapir2,mapdr3,mapir3,mapdr4,mapir4,ssn,key,aeqb)
c
c     mapdn - direct map matrix corresponding to N (Input)
c     mapin - direct map matrix corresponding to N (Input)
c     mapdr1- inverse map matrix corresponding to R1 (Input)
c     mapir1- inverse map matrix corresponding to R1 (Input)
c     mapdr2- inverse map matrix corresponding to R2 (Input)
c     mapir2- inverse map matrix corresponding to R2 (Input)
c     mapdr3- inverse map matrix corresponding to R3 (Input)
c     mapir3- inverse map matrix corresponding to R3 (Input)
c     mapdr4- inverse map matrix corresponding to R4 (Input)
c     mapir4- inverse map matrix corresponding to R4 (Input)
c     ssn   - overall symmtry state of N (Input)
c     key   - define state of _a,_b  key   a, b (Input)
c     1    S  S
c     2    V  S
c     3    S  V
c     4    V  V
c     aqeb  - aqeb=0 if a.ne.b; =1 if a=b (Input)
c
c     this routine unpack: N _a_b (p,q) = <ab|pq>              key  aeqb
c     N _a_b (p,q) a>=b-bb  -> R1 _a,_b(ef)aaaa = <ab||ef>      4    0
c     -> R2 _a,_b(ef)bbbb = <ab||ef>     1-4   0
c     -> R3 _a,_b(e,f)abab = <ab||ef>    3,4   0,1
c     -> R4 _b,_a(e,f)abab = <ba||fe>    2,4   0,1
c
c     !N.B. mylim, ze pri II.3;II.4;III.3;III.4 maju byt possn +--+ a nie ++++
c     ako su terazky
c
#include "ccsd1.fh"
#include "wrk.fh"
c
       integer mapdn(0:512,1:6)
       integer mapdr1(0:512,1:6)
       integer mapdr2(0:512,1:6)
       integer mapdr3(0:512,1:6)
       integer mapdr4(0:512,1:6)
c
       integer mapin(1:8,1:8,1:8)
       integer mapir1(1:8,1:8,1:8)
       integer mapir2(1:8,1:8,1:8)
       integer mapir3(1:8,1:8,1:8)
       integer mapir4(1:8,1:8,1:8)
c
       integer ssn,key,aeqb
c
c     help variables
c
       integer symp,symq,dimp,dimq,dime,dimf,dimef,lengthn
       integer ir1,ir2,ir3,ir4,possr1,possr2,possr3,possr4
       integer possn,in,inp,inm,possnp,possnm
c
       do 1000 symp=1,nsym
       symq=mmul(ssn,symp)
c
       in=mapin(symp,1,1)
       lengthn=mapdn(in,2)
c
       if (lengthn.eq.0) goto 1000
c
       if (symp.eq.symq) then
c     symp=symq
c
       possn=mapdn(in,1)
       dimp=dimm(5,symp)
c
cI.1  def R1 - _a_b(ef)aaaa
c
       if ((key.eq.4).and.(aeqb.eq.0)) then
       ir1=mapir1(symp,1,1)
       possr1=mapdr1(ir1,1)
       dimef=mapdr1(ir1,2)
       if (mapdr1(ir1,2).gt.0) then
       call unpckhelp1 (wrk(possn),wrk(possr1),dimp,dimef,noa(symp),
     &                  nva(symp))
       end if
       end if
c
cI.2  def R2 - _a_b(ef)bbbb
c
       if (aeqb.eq.0) then
       ir2=mapir2(symp,1,1)
       possr2=mapdr2(ir2,1)
       dimef=mapdr2(ir2,2)
       if (mapdr2(ir2,2).gt.0) then
       call unpckhelp1 (wrk(possn),wrk(possr2),dimp,dimef,nob(symp),
     &                  nvb(symp))
       end if
       end if
c
cI.3  def R3 - _a_b(e,f)abab
c
       if ((key.eq.2).or.(key.eq.4)) then
       ir3=mapir3(symp,1,1)
       possr3=mapdr3(ir3,1)
       dime=dimm(3,symp)
       dimf=dimm(4,symq)
       if (mapdr3(ir3,2).gt.0) then
       call unpckhelp3 (wrk(possn),wrk(possr3),dimp,dimp,dime,
     & dimf,noa(symp),nva(symp),nob(symq),nvb(symq))
       end if
       end if
c
cI.4  def R4 - _b_a(e,f)abab
c
       if (((key.eq.3).or.(key.eq.4)).and.(aeqb.eq.0)) then
       ir4=mapir4(symp,1,1)
       possr4=mapdr4(ir4,1)
       dime=dimm(3,symp)
       dimf=dimm(4,symq)
       if (mapdr4(ir4,2).gt.0) then
       call unpckhelp4 (wrk(possn),wrk(possr4),dimp,dimp,dime,
     & dimf,noa(symp),nva(symp),nob(symq),nvb(symq))
       end if
       end if
c
       else if (symp.gt.symq) then
c     symp > symq
c
       inp=mapin(symp,1,1)
       inm=mapin(symq,1,1)
       possnp=mapdn(inp,1)
       possnm=mapdn(inm,1)
       dimp=dimm(5,symp)
       dimq=dimm(5,symq)
c
cII.1 def R1 - _a_b(ef)aaaa
c
       if (key.eq.4) then
       ir1=mapir1(symp,1,1)
       possr1=mapdr1(ir1,1)
       dime=dimm(3,symp)
       dimf=dimm(3,symq)
       if (mapdr1(ir1,2).gt.0) then
       call unpckhelp2 (wrk(possnp),wrk(possnm),wrk(possr1),dimp,dimq,
     &                  dime,dimf,
     & noa(symp),nva(symp),noa(symq),nva(symq))
       end if
       end if
c
cII.2 def R2 - _a_b(ef)bbbb
c
       ir2=mapir2(symp,1,1)
       possr2=mapdr2(ir2,1)
       dime=dimm(4,symp)
       dimf=dimm(4,symq)
       if (mapdr2(ir2,2).gt.0) then
       call unpckhelp2 (wrk(possnp),wrk(possnm),wrk(possr2),dimp,dimq,
     &                  dime,dimf,
     & nob(symp),nvb(symp),nob(symq),nvb(symq))
       end if
c
cII.3 def R3 - _a_b(e,f)abab
c
       if ((key.eq.2).or.(key.eq.4)) then
       ir3=mapir3(symp,1,1)
       possr3=mapdr3(ir3,1)
       dime=dimm(3,symp)
       dimf=dimm(4,symq)
       if (mapdr3(ir3,2).gt.0) then
       call unpckhelp3 (wrk(possnp),wrk(possr3),dimp,dimq,dime,dimf,
     & noa(symp),nva(symp),nob(symq),nvb(symq))
       end if
       end if
c
cII.4 def R4 - _b_a(e,f)abab
c
       if ((key.eq.3).or.(key.eq.4)) then
       ir4=mapir4(symq,1,1)
       possr4=mapdr4(ir4,1)
       dime=dimm(3,symq)
       dimf=dimm(4,symp)
       if (mapdr4(ir4,2).gt.0) then
       call unpckhelp4 (wrk(possnp),wrk(possr4),dimp,dimq,dime,dimf,
     & noa(symq),nva(symq),nob(symp),nvb(symp))
       end if
       end if
c
       else
c     symp < symq
c%
       inp=mapin(symp,1,1)
       inm=mapin(symq,1,1)
       possnp=mapdn(inp,1)
       possnm=mapdn(inm,1)
       dimp=dimm(5,symp)
       dimq=dimm(5,symq)
c
cIII.1def R1 - _a_b(ef)aaaa -no operation
cIII.2def R2 - _a_b(ef)bbbb -no operation
c
cIII.3def R3 - _a_b(e,f)abab
c
       if ((key.eq.2).or.(key.eq.4)) then
       ir3=mapir3(symp,1,1)
       possr3=mapdr3(ir3,1)
       dime=dimm(3,symp)
       dimf=dimm(4,symq)
       if (mapdr3(ir3,2).gt.0) then
       call unpckhelp3 (wrk(possnp),wrk(possr3),dimp,dimq,dime,dimf,
     & noa(symp),nva(symp),nob(symq),nvb(symq))
       end if
       end if
c
cIII.4def R4 - _b_a(e,f)abab
c
       if ((key.eq.3).or.(key.eq.4)) then
       ir4=mapir4(symq,1,1)
       possr4=mapdr4(ir4,1)
       dime=dimm(3,symq)
       dimf=dimm(4,symp)
       if (mapdr4(ir4,2).gt.0) then
       call unpckhelp4 (wrk(possnp),wrk(possr4),dimp,dimq,dime,dimf,
     & noa(symq),nva(symq),nob(symp),nvb(symp))
       end if
       end if
c
c%
       end if
c
 1000   continue
c
       return
       end

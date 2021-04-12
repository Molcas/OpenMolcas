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
       subroutine unpackab3 (wrk,wrksize,
     & mapdn,mapin,mapdr1,mapir1,ssn,nabnow,possab0,lentotab,key)
c
c     mapdn - direct map matrix corresponding to N (Input)
c     mapin - direct map matrix corresponding to N (Input)
c     mapdr1- inverse map matrix corresponding to R1 (Input)
c     mapir1- inverse map matrix corresponding to R1 (Input)
c     ssn   - overall symmtry state of N (Input)
c     nabnow- #of b indexes stored in stack (dimension of stack)
c     possab0 - real possition of AB_stack (I)
c     lentotab- lenrgth of one record in stack (per one b ind) (I)
c     key   - key, specifying which integrals need to be defined (I)
c       this routine unpack: N _a_b (p,q) = <ab|pq>              key
c       N _a_b (p,q) a>=b-bb  -> R1 _a,_b(ef)aaaa = <ab||ef>      1
c                             -> R1 _a,_b(ef)bbbb = <ab||ef>      2
c                             -> R1 _a,_b(e,f)abab = <ab||ef>     3
c                             -> R1 _b,_a(e,f)abab = <ba||fe>     4
c
c     !N.B. mylim, ze pri II.3;II.4;III.3;III.4 maju byt possn +--+ a nie ++++
c     ako su terazky
c
#include "ccsd1.fh"
#include "wrk.fh"
c
       integer mapdn(0:512,1:6)
       integer mapin(1:8,1:8,1:8)
       integer mapdr1(0:512,1:6)
       integer mapir1(1:8,1:8,1:8)
c
       integer ssn,key,nabnow,possab0,lentotab
c
c     help variables
c
       integer symp,symq,dimp,dimq,dime,dimf,dimef,lengthn
       integer ir1,possr1,bb
       integer possn,in,inp,inm,possnp,possnm
c
c
       do 2000 bb=1,nabnow
c
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
cold   possn=mapdn(in,1)
       possn=possab0+(mapdn(in,1)-mapdn(1,1))+(bb-1)*lentotab
       dimp=dimm(5,symp)
c
cI.1  def R1 - _a_b(ef)aaaa
c
       if (key.eq.1) then
       ir1=mapir1(symp,1,1)
       possr1=mapdr1(ir1,1)
       dimef=(dimm(3,symp)*(dimm(3,symp)-1))/2
       if (mapdr1(ir1,2).gt.0) then
       call unpckhelp8 (wrk(possn),wrk(possr1),dimp,dimef,noa(symp),
     &                  nva(symp),bb,nabnow)
       end if
       end if
c
cI.2  def R1 - _a_b(ef)bbbb
c
       if (key.eq.2) then
       ir1=mapir1(symp,1,1)
       possr1=mapdr1(ir1,1)
       dimef=(dimm(4,symp)*(dimm(4,symp)-1))/2
       if (mapdr1(ir1,2).gt.0) then
       call unpckhelp8 (wrk(possn),wrk(possr1),dimp,dimef,nob(symp),
     &                  nvb(symp),bb,nabnow)
       end if
       end if
c
cI.3  def R1 - _a_b(e,f)abab
c
       if (key.eq.3) then
       ir1=mapir1(symp,1,1)
       possr1=mapdr1(ir1,1)
       dime=dimm(3,symp)
       dimf=dimm(4,symq)
       if (mapdr1(ir1,2).gt.0) then
       call unpckhelp10(wrk(possn),wrk(possr1),dimp,dimp,dime,
     & dimf,noa(symp),nva(symp),nob(symq),nvb(symq),bb,nabnow)
       end if
       end if
c
cI.4  def R1 - _b_a(e,f)abab
c
       if (key.eq.4) then
       ir1=mapir1(symp,1,1)
       possr1=mapdr1(ir1,1)
       dime=dimm(3,symp)
       dimf=dimm(4,symq)
       if (mapdr1(ir1,2).gt.0) then
       call unpckhelp11(wrk(possn),wrk(possr1),dimp,dimp,dime,
     & dimf,noa(symp),nva(symp),nob(symq),nvb(symq),bb,nabnow)
       end if
       end if
c
       else if (symp.gt.symq) then
c     symp > symq
c
       inp=mapin(symp,1,1)
       inm=mapin(symq,1,1)
cold   possnp=mapdn(inp,1)
       possnp=possab0+(mapdn(inp,1)-mapdn(1,1))+(bb-1)*lentotab
cold   possnm=mapdn(inm,1)
       possnm=possab0+(mapdn(inm,1)-mapdn(1,1))+(bb-1)*lentotab
       dimp=dimm(5,symp)
       dimq=dimm(5,symq)
c
cII.1 def R1 - _a_b(ef)aaaa
c
       if (key.eq.1) then
       ir1=mapir1(symp,1,1)
       possr1=mapdr1(ir1,1)
       dime=dimm(3,symp)
       dimf=dimm(3,symq)
       if (mapdr1(ir1,2).gt.0) then
       call unpckhelp9 (wrk(possnp),wrk(possnm),wrk(possr1),dimp,dimq,
     &                  dime,dimf,
     & noa(symp),nva(symp),noa(symq),nva(symq),bb,nabnow)
       end if
       end if
c
cII.2 def R1 - _a_b(ef)bbbb
c
       if (key.eq.2) then
       ir1=mapir1(symp,1,1)
       possr1=mapdr1(ir1,1)
       dime=dimm(4,symp)
       dimf=dimm(4,symq)
       if (mapdr1(ir1,2).gt.0) then
       call unpckhelp9 (wrk(possnp),wrk(possnm),wrk(possr1),dimp,dimq,
     &                  dime,dimf,
     & nob(symp),nvb(symp),nob(symq),nvb(symq),bb,nabnow)
       end if
       end if
c
cII.3 def R1 - _a_b(e,f)abab
c
       if (key.eq.3) then
       ir1=mapir1(symp,1,1)
       possr1=mapdr1(ir1,1)
       dime=dimm(3,symp)
       dimf=dimm(4,symq)
       if (mapdr1(ir1,2).gt.0) then
       call unpckhelp10(wrk(possnp),wrk(possr1),dimp,dimq,dime,dimf,
     & noa(symp),nva(symp),nob(symq),nvb(symq),bb,nabnow)
       end if
       end if
c
cII.4 def R4 - _b_a(e,f)abab
c
       if (key.eq.4) then
       ir1=mapir1(symq,1,1)
       possr1=mapdr1(ir1,1)
       dime=dimm(3,symq)
       dimf=dimm(4,symp)
       if (mapdr1(ir1,2).gt.0) then
       call unpckhelp11(wrk(possnp),wrk(possr1),dimp,dimq,dime,dimf,
     & noa(symq),nva(symq),nob(symp),nvb(symp),bb,nabnow)
       end if
       end if
c
       else
c     symp < symq
c%

       inp=mapin(symp,1,1)
       inm=mapin(symq,1,1)
cold   possnp=mapdn(inp,1)
cold   possnm=mapdn(inm,1)
       possnp=possab0+(mapdn(inp,1)-mapdn(1,1))+(bb-1)*lentotab
       possnm=possab0+(mapdn(inm,1)-mapdn(1,1))+(bb-1)*lentotab
       dimp=dimm(5,symp)
       dimq=dimm(5,symq)
c
cIII.1def R1 - _a_b(ef)aaaa -no operation
cIII.2def R1 - _a_b(ef)bbbb -no operation
c
cIII.3def R1 - _a_b(e,f)abab
c
       if (key.eq.3) then
       ir1=mapir1(symp,1,1)
       possr1=mapdr1(ir1,1)
       dime=dimm(3,symp)
       dimf=dimm(4,symq)
       if (mapdr1(ir1,2).gt.0) then
       call unpckhelp10(wrk(possnp),wrk(possr1),dimp,dimq,dime,dimf,
     & noa(symp),nva(symp),nob(symq),nvb(symq),bb,nabnow)
       end if
       end if
c
cIII.4def R3 - _b_a(e,f)abab
c
       if (key.eq.4) then
       ir1=mapir1(symq,1,1)
       possr1=mapdr1(ir1,1)
       dime=dimm(3,symq)
       dimf=dimm(4,symp)
       if (mapdr1(ir1,2).gt.0) then
       call unpckhelp11(wrk(possnp),wrk(possr1),dimp,dimq,dime,dimf,
     & noa(symq),nva(symq),nob(symp),nvb(symp),bb,nabnow)
       end if
       end if
c
c%
       end if
c
 1000   continue
 2000   continue
c
       return
       end

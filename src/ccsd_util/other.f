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
c       This file contains following routines:
c     unpackab1
c     unpackab2
c       unpckhelp1
c       unpckhelp2
c       unpckhelp3
c       unpckhelp4
c       unpckhelp5
c       unpckhelp6
c       unpckhelp7
c
c     saverest1
c     saverest2
c
c       files for stacking:
c     unpackab3
c       unpckhelp8
c       unpckhelp9
c       unpckhelp10
c       unpckhelp11
c
c     --------------------
c     --------------------
c
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
       integer symp,symq,dimp,dimq,dime,dimf,dimef,lenghtn
       integer ir1,ir2,ir3,ir4,possr1,possr2,possr3,possr4
       integer possn,in,inp,inm,possnp,possnm
c
       do 1000 symp=1,nsym
       symq=mmul(ssn,symp)
c
       in=mapin(symp,1,1)
       lenghtn=mapdn(in,2)
c
       if (lenghtn.eq.0) goto 1000
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
c
c     --------------------
c
       subroutine unpackab2 (wrk,wrksize,
     & mapdn,mapin,mapdr1,mapir1,mapdr2,mapir2,
     &                       mapdr3,mapir3,
     & mapdr4,mapir4,mapdr5,mapir5,mapdr6,mapir6,ssn,key,aeqb)
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
c     mapdr5- inverse map matrix corresponding to R4 (Input)
c     mapir5- inverse map matrix corresponding to R4 (Input)
c     mapdr6- inverse map matrix corresponding to R4 (Input)
c     mapir6- inverse map matrix corresponding to R4 (Input)
c     ssn   - overall symmtry state of N (Input)
c     key   - define state of _a,_b  key   a, b (Input)
c     1    S  S
c     2    V  S
c     3    S  V
c     4    V  V
c     aqeb  - aqeb=0 if a.ne.b; =1 if a=b (Input)
c
c     this routine unpack: N _a_b (p,q) = <ab|pq>              key   aeqb
c     N _a_b (p,q) a>=b-aa  -> R1 _a,_b(j,e)aaaa = <ab||je>     4     0
c     -> R2 _a,_b(j,e)bbbb = <ab||je>    1-4    0
c     -> R3 _a,_b(j,e)abab = <ab||je>    3,4   0,1
c     -> R4 _b,_a(j,e)abab = <ba||je>    2,4   0,1
c     -> R5 _a,_b(j,e)abba = <ab||je>    3,4   0,1
c     -> R6 _b,_a(j,e)abba = <ba||je>    2,4   0,1
c
#include "ccsd1.fh"
#include "wrk.fh"
c
       integer mapdn(0:512,1:6)
       integer mapdr1(0:512,1:6)
       integer mapdr2(0:512,1:6)
       integer mapdr3(0:512,1:6)
       integer mapdr4(0:512,1:6)
       integer mapdr5(0:512,1:6)
       integer mapdr6(0:512,1:6)
c
       integer mapin(1:8,1:8,1:8)
       integer mapir1(1:8,1:8,1:8)
       integer mapir2(1:8,1:8,1:8)
       integer mapir3(1:8,1:8,1:8)
       integer mapir4(1:8,1:8,1:8)
       integer mapir5(1:8,1:8,1:8)
       integer mapir6(1:8,1:8,1:8)
c
       integer ssn,key,aeqb
c
c     help variables
c
       integer symp,symq,dimp,dimq,dime,dimj
       integer in,inm,inp,ir1,ir2,ir3,ir4,ir5,ir6
       integer possn,possnp,possnm,possr1,possr2,possr3,possr4,possr5,
     & possr6
       integer lenghtn
c
       do 1000 symp=1,nsym
       symq=mmul(ssn,symp)
c
       in=mapin(symp,1,1)
       lenghtn=mapdn(in,2)
c
       if (lenghtn.eq.0) goto 1000
c
       if (symp.eq.symq) then
c     symp=symq
c
       possn=mapdn(in,1)
       dimp=dimm(5,symp)
c
cI.1  def R1 - _a_b(j,e)aaaa
c
       if ((key.eq.4).and.(aeqb.eq.0)) then
       ir1=mapir1(symp,1,1)
       possr1=mapdr1(ir1,1)
       dimj=noa(symp)
       dime=nva(symq)
       if (mapdr1(ir1,2).gt.0) then
       call unpckhelp5 (wrk(possn),wrk(possr1),dimp,dimj,dime,
     & 0,noa(symp),noa(symp),nva(symp))
       end if
       end if
c
cI.2  def R2 - _a_b(j,e)bbbb
c
       if (aeqb.eq.0) then
       ir2=mapir2(symp,1,1)
       possr2=mapdr2(ir2,1)
       dimj=nob(symp)
       dime=nvb(symq)
       if (mapdr2(ir2,2).gt.0) then
       call unpckhelp5 (wrk(possn),wrk(possr2),dimp,dimj,dime,
     & 0,nob(symp),nob(symp),nvb(symp))
       end if
       end if
c
cI.3  def R3 - _a_b(j,e)abba
c
       if ((key.eq.2).or.(key.eq.4)) then
       ir3=mapir3(symp,1,1)
       possr3=mapdr3(ir3,1)
       dimj=nob(symp)
       dime=nva(symq)
       if (mapdr3(ir3,2).gt.0) then
       call unpckhelp7 (wrk(possn),wrk(possr3),dimp,dimp,dimj,
     & dime,0,nob(symp),noa(symq),nva(symq))
       end if
       end if
c
cI.4  def R4 - _b_a(j,e)abba
c
       if (((key.eq.3).or.(key.eq.4)).and.(aeqb.eq.0)) then
       ir4=mapir4(symp,1,1)
       possr4=mapdr4(ir4,1)
       dimj=nob(symp)
       dime=nva(symq)
       if (mapdr4(ir4,2).gt.0) then
       call unpckhelp6 (wrk(possn),wrk(possr4),dimp,dimp,dimj,
     & dime,0,nob(symp),noa(symq),nva(symq))
       end if
       end if
c
cI.5  def R5 - _a_b (j,e)abab
c
       if ((key.eq.2).or.(key.eq.4)) then
       ir5=mapir5(symp,1,1)
       possr5=mapdr5(ir5,1)
       dimj=noa(symp)
       dime=nvb(symq)
       if (mapdr5(ir5,2).gt.0) then
       call unpckhelp3 (wrk(possn),wrk(possr5),dimp,dimp,dimj,
     & dime,0,noa(symp),nob(symq),nvb(symq))
       end if
       end if
c
cI.6  def R6 - _b_a(j,e)abab
c
       if (((key.eq.3).or.(key.eq.4)).and.(aeqb.eq.0)) then
       ir6=mapir6(symp,1,1)
       possr6=mapdr6(ir6,1)
       dimj=noa(symp)
       dime=nvb(symq)
       if (mapdr6(ir6,2).gt.0) then
       call unpckhelp4 (wrk(possn),wrk(possr6),dimp,dimp,dimj,
     & dime,0,noa(symp),nob(symq),nvb(symq))
       end if
       end if
c
       else
c     symp.ne.symq
c
       inp=mapin(symp,1,1)
       inm=mapin(symq,1,1)
       possnp=mapdn(inp,1)
       possnm=mapdn(inm,1)
       dimp=dimm(5,symp)
       dimq=dimm(5,symq)
c
cII.1 def R1 - _a_b(j,e)aaaa
c
       if (key.eq.4) then
       ir1=mapir1(symp,1,1)
       possr1=mapdr1(ir1,1)
       dimj=noa(symp)
       dime=nva(symq)
       if (mapdr1(ir1,2).gt.0) then
       call unpckhelp2 (wrk(possnp),wrk(possnm),wrk(possr1),dimp,dimq,
     &                  dimj,dime,
     & 0,noa(symp),noa(symq),nva(symq))
       end if
       end if
c
cII.2 def R2 - _a_b(j,e)bbbb
c
       ir2=mapir2(symp,1,1)
       possr2=mapdr2(ir2,1)
       dimj=nob(symp)
       dime=nvb(symq)
       if (mapdr2(ir2,2).gt.0) then
       call unpckhelp2 (wrk(possnp),wrk(possnm),wrk(possr2),dimp,dimq,
     &                  dimj,dime,
     & 0,nob(symp),nob(symq),nvb(symq))
       end if
c
cII.3 def R3 - _a_b(j,e)abba
c
       if ((key.eq.2).or.(key.eq.4)) then
       ir3=mapir3(symq,1,1)
       possr3=mapdr3(ir3,1)
       dimj=nob(symq)
       dime=nva(symp)
       if (mapdr3(ir3,2).gt.0) then
       call unpckhelp7 (wrk(possnp),wrk(possr3),dimp,dimq,dimj,dime,
     & 0,nob(symq),noa(symp),nva(symp))
       end if
       end if
c
cII.4 def R4 - _b_a(j,e)abba
c
       if ((key.eq.3).or.(key.eq.4)) then
       ir4=mapir4(symp,1,1)
       possr4=mapdr4(ir4,1)
       dimj=nob(symp)
       dime=nva(symq)
       if (mapdr4(ir4,2).gt.0) then
       call unpckhelp6 (wrk(possnp),wrk(possr4),dimp,dimq,dimj,dime,
     & 0,nob(symp),noa(symq),nva(symq))
       end if
       end if
c
cII.5 def R5 - _a_b(j,e)abab
c
       if ((key.eq.2).or.(key.eq.4)) then
       ir5=mapir5(symp,1,1)
       possr5=mapdr5(ir5,1)
       dimj=noa(symp)
       dime=nvb(symq)
       if (mapdr5(ir5,2).gt.0) then
       call unpckhelp3 (wrk(possnp),wrk(possr5),dimp,dimq,dimj,dime,
     & 0,noa(symp),nob(symq),nvb(symq))
       end if
       end if
c
cII.6 def R6 - _b_a(j,e)abab
c
       if ((key.eq.3).or.(key.eq.4)) then
       ir6=mapir6(symq,1,1)
       possr6=mapdr6(ir6,1)
       dimj=noa(symq)
       dime=nvb(symp)
       if (mapdr6(ir6,2).gt.0) then
       call unpckhelp4 (wrk(possnp),wrk(possr6),dimp,dimq,dimj,dime,
     & 0,noa(symq),nob(symp),nvb(symp))
       end if
       end if
c
       end if
c
 1000   continue
c
       return
       end
c
c     --------------------
c
       subroutine unpckhelp1 (a,b,dimp,dimef,eadd,noe)
c
c     this routine do:
c     b(ef) = a(pe,qf)-a(qf,pe) for symp=symq
c
       integer dimp,dimef,eadd,noe
       real*8 a(1:dimp,1:dimp)
       real*8 b(1:dimef)
c
c     help variables
       integer pe,qf,ef
c
       ef=0
       do 100 pe=eadd+2,eadd+noe
       do 100 qf=eadd+1,pe-1
       ef=ef+1
       b(ef)=a(pe,qf)-a(qf,pe)
 100    continue
c
       return
       end
c
c     ----------
c
       subroutine unpckhelp2 (ap,am,b,dimp,dimq,dime,dimf,eadd,noe,fadd,
     &                        nof)
c
c     this routine do:
c     b(e,f) = a(pe,qf)-a(qf,pe) for symp>symq
c
       integer dimp,dimq,dime,dimf,eadd,noe,fadd,nof
       real*8 ap(1:dimp,1:dimq)
       real*8 am(1:dimq,1:dimp)
       real*8 b(1:dime,1:dimf)
c
c     help variables
       integer pe,qf,f
c
       do 100 qf=fadd+1,fadd+nof
       f=qf-fadd
       do 100 pe=eadd+1,eadd+noe
       b(pe-eadd,f)=ap(pe,qf)-am(qf,pe)
 100    continue
c
       return
       end
c
c     ----------
c
       subroutine unpckhelp3 (a,b,dimp,dimq,dime,dimf,eadd,noe,fadd,nof)
c
c     this routine do:
c     b(e,f) = a(pe,qf)
c
       integer dimp,dimq,dime,dimf,eadd,noe,fadd,nof
       real*8 a(1:dimp,1:dimq)
       real*8 b(1:dime,1:dimf)
c
c     help variables
       integer pe,qf,f
c
       do 100 qf=fadd+1,fadd+nof
       f=qf-fadd
       do 100 pe=eadd+1,eadd+noe
       b(pe-eadd,f)=a(pe,qf)
 100    continue
c
       return
       end
c
c     ----------
c
       subroutine unpckhelp4 (a,b,dimp,dimq,dime,dimf,eadd,noe,fadd,nof)
c
c     this routine do:
c     b(e,f) =  a(pf,qe)
c
       integer dimp,dimq,dime,dimf,eadd,noe,fadd,nof
       real*8 a(1:dimp,1:dimq)
       real*8 b(1:dime,1:dimf)
c
c     help variables
       integer qe,pf,f
c
       do 100 pf=fadd+1,fadd+nof
       f=pf-fadd
       do 100 qe=eadd+1,eadd+noe
       b(qe-eadd,f)=a(pf,qe)
 100    continue
c
       return
       end
c
c     ----------
c
       subroutine unpckhelp5 (a,b,dimp,dimj,dime,jadd,noj,eadd,noe)
c
c     this routine do:
c     b(j,e) = a(pj,qe)-a(qe,pj) for symp=symq
c
       integer dimp,dime,dimj,eadd,noe,jadd,noj
       real*8 a(1:dimp,1:dimp)
       real*8 b(1:dimj,1:dime)
c
c     help variables
       integer pj,qe,e
c
       do 100 qe=eadd+1,eadd+noe
       e=qe-eadd
       do 100 pj=jadd+1,jadd+noj
       b(pj-jadd,e)=a(pj,qe)-a(qe,pj)
 100    continue
c
       return
       end
c
c     ----------
c
       subroutine unpckhelp6 (a,b,dimp,dimq,dime,dimf,eadd,noe,fadd,nof)
c
c     this routine do:
c     b(e,f) = -a(pe,qf)
c
       integer dimp,dimq,dime,dimf,eadd,noe,fadd,nof
       real*8 a(1:dimp,1:dimq)
       real*8 b(1:dime,1:dimf)
c
c     help variables
       integer pe,qf,f
c
       do 100 qf=fadd+1,fadd+nof
       f=qf-fadd
       do 100 pe=eadd+1,eadd+noe
       b(pe-eadd,f)=-a(pe,qf)
 100    continue
c
       return
       end
c
c     ----------
c
       subroutine unpckhelp7 (a,b,dimp,dimq,dime,dimf,eadd,noe,fadd,nof)
c
c     this routine do:
c     b(e,f) =  -a(pf,qe)
c
       integer dimp,dimq,dime,dimf,eadd,noe,fadd,nof
       real*8 a(1:dimp,1:dimq)
       real*8 b(1:dime,1:dimf)
c
c     help variables
       integer qe,pf,f
c
       do 100 pf=fadd+1,fadd+nof
       f=pf-fadd
       do 100 qe=eadd+1,eadd+noe
       b(qe-eadd,f)=-a(pf,qe)
 100    continue
c
       return
       end
c
c     -------------------------
c
       subroutine saverest1 (wrk,wrksize,
     & lunrst)
c
c     this routine save restart informations:
c     t13,t14,t21,t22,t23
c
#include "ccsd1.fh"
#include "ccsd2.fh"
#include "wrk.fh"
c
       integer lunrst
c
c     help variables
c
       integer rc
c
c0    return if need
       if (keyrst.eq.0) return
c
c1    rewind tape
       call filemanager (2,lunrst,rc)
c
c2    write T1aa
       call wrtmediate (wrk,wrksize,
     & lunrst,mapdt13,mapit13,rc)
c
c3    write T1bb
       call wrtmediate (wrk,wrksize,
     & lunrst,mapdt14,mapit14,rc)
c
c4    write T2aaaa
       call wrtmediate (wrk,wrksize,
     & lunrst,mapdt21,mapit21,rc)
c
c5    write T2bbbb
       call wrtmediate (wrk,wrksize,
     & lunrst,mapdt22,mapit22,rc)
c
c6    write T2abab
       call wrtmediate (wrk,wrksize,
     & lunrst,mapdt23,mapit23,rc)
c
       return
       end
c
c     -------------------------
c
       subroutine saverest2 (lunrst,energy,niter,iokey,daddr)
c
c     this routine save restart informations:
c     energy, niter
c     to prepaired possition in lunrst
c

#include "SysDef.fh"
       integer lunrst,niter,iokey,daddr
       real*8 energy
c
c1    write energy,niter
       if (iokey.eq.1) then
c      Fortran IO
       write (lunrst) energy,niter
c
       else
c      MOLCAS IO
       call ddafile (lunrst,1,[energy],1,daddr)
       call idafile (lunrst,1,[niter],1,daddr)
       end if
c
       return
       end
c
c     -------------------------
c
c       files for stacking:
c     unpackab3
c       unpckhelp8
c       unpckhelp9
c       unpckhelp10
c       unpckhelp11
c
c     --------------------
c
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
       integer symp,symq,dimp,dimq,dime,dimf,dimef,lenghtn
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
       lenghtn=mapdn(in,2)
c
       if (lenghtn.eq.0) goto 1000
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
c
c     --------------------
c
       subroutine unpckhelp8 (a,b,dimp,dimef,eadd,noe,bb,dimb)
c
c     this routine do:
c     b(ef,_Bb) = a(pe,qf)-a(qf,pe) for symp=symq
c
       integer dimp,dimef,eadd,noe,bb,dimb
       real*8 a(1:dimp,1:dimp)
       real*8 b(1:dimef,1:dimb)
c
c     help variables
       integer pe,qf,ef
c
       ef=0
       do 100 pe=eadd+2,eadd+noe
       do 100 qf=eadd+1,pe-1
       ef=ef+1
       b(ef,bb)=a(pe,qf)-a(qf,pe)
 100    continue
c
       return
       end
c
c     ----------
c
       subroutine unpckhelp9 (ap,am,b,dimp,dimq,dime,dimf,eadd,noe,fadd,
     &                        nof,bb,dimb)
c
c     this routine do:
c     b(e,f,_Bb) = a(pe,qf)-a(qf,pe) for symp>symq
c
       integer dimp,dimq,dime,dimf,eadd,noe,fadd,nof,bb,dimb
       real*8 ap(1:dimp,1:dimq)
       real*8 am(1:dimq,1:dimp)
       real*8 b(1:dime,1:dimf,1:dimb)
c
c     help variables
       integer pe,qf,f
c
       do 100 qf=fadd+1,fadd+nof
       f=qf-fadd
       do 100 pe=eadd+1,eadd+noe
       b(pe-eadd,f,bb)=ap(pe,qf)-am(qf,pe)
 100    continue
c
       return
       end
c
c     ----------
c
       subroutine unpckhelp10(a,b,dimp,dimq,dime,dimf,eadd,noe,fadd,nof,
     c                        bb,dimb)
c
c     this routine do:
c     b(e,f,_Bb) = a(pe,qf)
c
       integer dimp,dimq,dime,dimf,eadd,noe,fadd,nof,bb,dimb
       real*8 a(1:dimp,1:dimq)
       real*8 b(1:dime,1:dimf,1:dimb)
c
c     help variables
       integer pe,qf,f
c
       do 100 qf=fadd+1,fadd+nof
       f=qf-fadd
       do 100 pe=eadd+1,eadd+noe
       b(pe-eadd,f,bb)=a(pe,qf)
 100    continue
c
       return
       end
c
c     ----------
c
       subroutine unpckhelp11(a,b,dimp,dimq,dime,dimf,eadd,noe,fadd,nof,
     c                        bb,dimb)
c
c     this routine do:
c     b(e,f,_Bb) =  a(pf,qe)
c
       integer dimp,dimq,dime,dimf,eadd,noe,fadd,nof,bb,dimb
       real*8 a(1:dimp,1:dimq)
       real*8 b(1:dime,1:dimf,1:dimb)
c
c     help variables
       integer qe,pf,f
c
       do 100 pf=fadd+1,fadd+nof
       f=pf-fadd
       do 100 qe=eadd+1,eadd+noe
       b(qe-eadd,f,bb)=a(pf,qe)
 100    continue
c
       return
       end
c
c     -------------------------
c
        subroutine prmap (mapd,mapi)
c

       integer mapd(0:512,1:6)
       integer mapi(1:8,1:8,1:8)
       integer ii,i,j
       ii=mapd(0,5)
       do i=0,ii
       write (6,99) i,(mapd(i,j),j=1,6)
99      format (i3,6x,i10,1x,5(i6,1x))
       end do
c
        write (6,*) mapi(1,1,1),mapi(2,1,1)
       return
       end

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
c     this file contains folloeing routines:
c     add
c     add10
c     add21
c     add22
c     add32
c     add41
c     add42
c     add43
c     add44
c
c     -----------------------------------
c
       subroutine cct3_add (wrk,wrksize,
     & ninda,nindb,nindext,typext,u,v,ssu,ssv,factor,
     & mapda,ssa,mapdb,mapib,ssb,rc)
c
c     this routine do:
c     B(indb) = B(indb) + factor * A(inda)
c
c     ninda   - # of indexes in matrix A (1-4)
c     nindb   - # of indexes in matrix B (1-4)
c     nindext - # of external (frozen, fixed) indexes (now:0-2)
c     typext  - characterize external indexes as follows:
c     0 - no frozen index
c     1 - frozen index p
c     2 - frozen index q
c     3 - frozen index r
c     4 - frozen index s
c     5 - frozen indexes p,q
c     6 - frozen indexes r,s
c     u       - value of first external index (if any, else 0)
c     v       - value of second external index (if any, else 0)
c     ssu     - symmetry of u (if any, else 1)
c     ssv     - symmetry of v (if any, else 1)
c     factor  - multiplicative factor (see def)
c     mapda   - direct map matrix corresponding to A (see docc.txt)
c     ssa     - overall spin state of matrix A
c     mapdb   - direct map matrix corresponding to B (see docc.txt)
c     mapib   - inverse map matrix corresponding to B (see docc.txt)
c     ssb     - overall spin state of matrix B
c     rc      - return (error) code
c
c     Table of present implementations:
c
c     nindB  nindxet  typext  typB  =>   Implemented
c     >4                                   No
c
c     4       0        0    0-4            Yes
c     4       1       1-4   0,4            Yes
c     4       1       1-4   2,3            No
c     4       2        5    0,4            Yes
c     4       2        5    2,3            No
c     4       2       6-n                  No
c     4       3                            No
c
c     3       0       0     0-2            Yes
c     3       1       1-3    0             Yes
c     3       1       1-3   1,2            No
c     3       2                            No
c
c     2       0       0     0,1            Yes
c     2       1       1-2    0             Yes
c     2       1       1-2    1             No
c     2       2                            No
c
c     1                                    No
c
c
c     !N.B. oprav co je oznacene c@!
c
#include "t31.fh"
#include "wrk.fh"
c
       integer ninda,nindb,nindext,typext,u,v,ssu,ssv,ssa,ssb,rc
       real*8 factor
       integer mapda(0:512,1:6)
       integer mapdb(0:512,1:6)
c
       integer mapib(1:8,1:8,1:8)
c
c     help variables
c
       integer sa1,sa2,sa3,sa4,ssp,ssq,pq
       integer nhelp1,nhelp2,nhelp3,nhelp4,nhelp5
       integer nhelp6,nhelp7,nhelp8,nhelp9,nhelp10
       integer ia,ib,ibm
       integer typa,typb,p,q
       real*8 fact
c      To fix some 'uninitialized' warnings
       p=0
       q=0
c     general tests
c
       nhelp1=nindA+nindext
c
       if (nhelp1.ne.nindb) then
c     RC=1  : incompatible (nindA, nindB and nindext, Stup)
       rc=1
       return
       end if
c
       nhelp1=mmul(ssu,ssv)
       nhelp1=mmul(ssa,nhelp1)
       if (nhelp1.ne.ssb) then
c     RC=2  : incompatible (ssa, ssb ,ssu and ssv, Stup)
       rc=2
       return
       end if
c
       typa=mapda(0,6)
       typb=mapdb(0,6)
       fact=factor
c
       if (nindext.gt.0) then
       if ((typb.ge.1).and.(typb.le.3)) then
c     RC=3 : nindext>0, typB is 1,2 or 3 (NCI)
       rc=3
       return
       end if
       end if

c
       if ((nindext.eq.2).and.(typb.eq.4)) then
c
c     def p,q,ssp,ssq,fact(new)
c
       if (ssu.gt.ssv) then
c     ssu>ssv
       p=u
       q=v
       ssp=ssu
       ssq=ssv
       fact=factor
       else if (ssu.eq.ssv) then
c     ssu=ssv
       if (u.ge.v) then
       p=u
       q=v
       ssp=ssu
       ssq=ssv
       fact=factor
       else
       p=v
       q=u
       ssp=ssv
       ssq=ssu
       fact=-factor
       end if
       else
c     ssu<ssv
       p=v
       q=u
       ssp=ssv
       ssq=ssu
       fact=-factor
       end if
c
       end if
c
       if (nindb.eq.4) then
c
c     **********  -> B(pqrs) **********
c
       if (nindext.eq.0) then
c
c400  case B(pqrs) <-- A(pqrs) or B(p,q,r,s) <-- A(p,q,r,s)
c
c     tests
c
       if (typa.ne.typb) then
c     RC=4 : nindB=4, nindext=0 (TypA incompatible with TypB ,Stup)
       rc=4
       return
       end if
c
       do 400 ia=1,mapda(0,5)
c
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
       sa3=mapda(ia,5)
       sa4=mapda(ia,6)
c
       ib=mapib(sa1,sa2,sa3)
c
c     def lenght
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 400
c
c     def possA,possB
       nhelp2=mapda(ia,1)
       nhelp3=mapdb(ib,1)
c
       call cct3_add10 (wrk(nhelp2),wrk(nhelp3),nhelp1,fact)
c
 400    continue
c
       else if (nindext.eq.1) then
c
       if (typext.eq.1) then
c
       if (typb.eq.0) then
c
c4110 case B(p,q,r,s) <-- A(q,r,s)
c
c     tsets
c
       if (typa.ne.0) then
c     RC=5 : nindB=4, nindeext=1, typext=1, typB=0 (typA is not 0, Stup)
       rc=5
       return
       end if
c
       do 4110 ia=1,mapda(0,5)
c
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
       sa3=mapda(ia,5)
c
       ib=mapib(ssu,sa1,sa2)
c
c     def lenght
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 4110
c
c     def possA,possB
       nhelp2=mapda(ia,1)
       nhelp3=mapdb(ib,1)
c
c     def dimp,dimq,dimr,dims
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
       nhelp7=dimm(mapdb(0,4),mapdb(ib,6))
c
c     def fictive dimensions
       nhelp8=nhelp5*nhelp6*nhelp7
c
       call cct3_add21 (wrk(nhelp2),wrk(nhelp3),u,nhelp4,nhelp8,fact)
c
 4110   continue
c
       else if (typb.eq.4) then
c
c4114 case B(pq,rs) <-- A(q,rs)
c
c     tsets
c
       if (typa.ne.2) then
c     RC=6  : nindB=4, nindeext=1, typext=1, typB=4 (typA is not 2, Stup)
       rc=6
       return
       end if
c
       do 4114 ia=1,mapda(0,5)
c
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
       sa3=mapda(ia,5)
c
       ib=mapib(ssu,sa1,sa2)
       ibm=mapib(sa1,ssu,sa2)
c
c     def lenght
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 4114
c
c     def possA
       nhelp2=mapda(ia,1)
c
c
       if (ssu.gt.sa1) then
c
c     def dimp,dimq,dimr,dims
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
       nhelp7=dimm(mapdb(0,4),mapdb(ib,6))
c
c     def fictive dimensions
       if (sa2.eq.sa3) then
       nhelp8=nhelp6*(nhelp6-1)/2
       else
       nhelp8=nhelp6*nhelp7
       end if
c
c     def possB
       nhelp3=mapdb(ib,1)
c     def fictive dimensions
       nhelp9=nhelp5*nhelp8
       call cct3_add21 (wrk(nhelp2),wrk(nhelp3),u,nhelp4,nhelp9,fact)
c
       else if (ssu.eq.sa1) then
c
c     def dimp,dimq,dimr,dims
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
       nhelp7=dimm(mapdb(0,4),mapdb(ib,6))
c
c     def fictive dimensions
       if (sa2.eq.sa3) then
       nhelp8=nhelp6*(nhelp6-1)/2
       else
       nhelp8=nhelp6*nhelp7
       end if
c
c     def possB
       nhelp3=mapdb(ib,1)
c     def fictive dimensions
       nhelp9=nhelp4*(nhelp4-1)/2
       call cct3_add41 (wrk(nhelp2),wrk(nhelp3),
     &  u,nhelp4,nhelp9,nhelp8,fact)
c
       else
c     ssu<sa1  B(qp,rs) <-- -A_p (q,rs)
c
c     def dimp,dimq,dimr,dims
       nhelp4=dimm(mapdb(0,1),mapdb(ibm,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ibm,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ibm,5))
       nhelp7=dimm(mapdb(0,4),mapdb(ibm,6))
c
c     def fictive dimensions
       if (sa2.eq.sa3) then
       nhelp8=nhelp6*(nhelp6-1)/2
       else
       nhelp8=nhelp6*nhelp7
       end if
c
c     def possB-
       nhelp3=mapdb(ibm,1)
       call cct3_add32 (wrk(nhelp2),wrk(nhelp3),
     &  u,nhelp4,nhelp5,nhelp8,-fact)
c
       end if
c
 4114   continue
c
       else
c     RC=7 : nindB=4, nindext=1, typext=1, (typA is not 0 or 4, (NCI))
       rc=7
       return
       end if
c
       else if (typext.eq.2) then
c
       if (typb.eq.0) then
c
c4120 case B(p,q,r,s) <-- A(p,r,s)
c
c     tsets
c
       if (typa.ne.0) then
c     RC=8 : nindB=4, nindeext=1, typext=2, typB=0 (typA is not 0, Stup)
       rc=8
       return
       end if
c
       do 4120 ia=1,mapda(0,5)
c
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
       sa3=mapda(ia,5)
c
       ib=mapib(sa1,ssu,sa2)
c
c     def lenght
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 4120
c
c     def possA,possB
       nhelp2=mapda(ia,1)
       nhelp3=mapdb(ib,1)
c
c     def dimp,dimq,dimr,dims
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
       nhelp7=dimm(mapdb(0,4),mapdb(ib,6))
c
c     def fictive dimensions
       nhelp8=nhelp6*nhelp7
c
       call cct3_add32 (wrk(nhelp2),wrk(nhelp3),
     &  u,nhelp4,nhelp5,nhelp8,fact)
c
 4120   continue
c
       else if (typb.eq.4) then
c
c4124 case B(pq,rs) <-- A(p,rs)
c
c     tsets
c
       if (typa.ne.2) then
c     RC=9 : nindB=4, nindeext=1, typext=2, typB=4 (typA is not 2, Stup)
       rc=9
       return
       end if
c
       do 4124 ia=1,mapda(0,5)
c
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
       sa3=mapda(ia,5)
c
       ib=mapib(sa1,ssu,sa2)
       ibm=mapib(ssu,sa1,sa2)
c
c     def lenght
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 4124
c
c     def possA
       nhelp2=mapda(ia,1)
c
c
       if (sa1.gt.ssu) then
c
c     def dimp,dimq,dimr,dims
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
       nhelp7=dimm(mapdb(0,4),mapdb(ib,6))
c
c     def fictive dimensions
       if (sa2.eq.sa3) then
       nhelp8=nhelp6*(nhelp6-1)/2
       else
       nhelp8=nhelp6*nhelp7
       end if
c
c     def possB
       nhelp3=mapdb(ib,1)
       call cct3_add32 (wrk(nhelp2),
     &  wrk(nhelp3),u,nhelp4,nhelp5,nhelp8,fact)
c
c
       else if (sa1.eq.ssu) then
c
c     def dimp,dimq,dimr,dims
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
       nhelp7=dimm(mapdb(0,4),mapdb(ib,6))
c
c     def fictive dimensions
       if (sa2.eq.sa3) then
       nhelp8=nhelp6*(nhelp6-1)/2
       else
       nhelp8=nhelp6*nhelp7
       end if
c
c     def possB
       nhelp3=mapdb(ib,1)
c     def fictive dimensions
       nhelp9=nhelp4*(nhelp4-1)/2
       call cct3_add42 (wrk(nhelp2),wrk(nhelp3),
     & u,nhelp4,nhelp9,nhelp8,fact)
c
c
       else
c     sa1<ssu  B(qp,rs) <-- -A_q (p,rs)
c
c     def dimp,dimq,dimr,dims
       nhelp4=dimm(mapdb(0,1),mapdb(ibm,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ibm,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ibm,5))
       nhelp7=dimm(mapdb(0,4),mapdb(ibm,6))
c
c     def fictive dimensions
       if (sa2.eq.sa3) then
       nhelp8=nhelp6*(nhelp6-1)/2
       else
       nhelp8=nhelp6*nhelp7
       end if
c
c     def possB-
       nhelp3=mapdb(ibm,1)
c     def fictive index
       nhelp9=nhelp8*nhelp5
       call cct3_add21 (wrk(nhelp2),wrk(nhelp3),u,nhelp4,nhelp9,-fact)
c
       end if
c
 4124   continue
c
       else
c     RC=10: nindB=4, nindext=1, typext=2, (typA is not 0 or 4, (NCI))
       rc=10
       return
       end if
c
       else if (typext.eq.3) then
c
       if (typb.eq.0) then
c
c4130 case B(p,q,r,s) <-- A(p,q,s)
c
c     tsets
c
       if (typa.ne.0) then
c     RC=11: nindB=4, nindeext=1, typext=3, typB=0 (typA is not 0, Stup)
       rc=11
       return
       end if
c
       do 4130 ia=1,mapda(0,5)
c
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
       sa3=mapda(ia,5)
c
       ib=mapib(sa1,sa2,ssu)
c
c     def lenght
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 4130
c
c     def possA,possB
       nhelp2=mapda(ia,1)
       nhelp3=mapdb(ib,1)
c
c     def dimp,dimq,dimr,dims
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
       nhelp7=dimm(mapdb(0,4),mapdb(ib,6))
c
c     def fictive dimensions
       nhelp8=nhelp4*nhelp5
c
       call cct3_add32 (wrk(nhelp2),wrk(nhelp3),
     &  u,nhelp8,nhelp6,nhelp7,fact)
c
 4130   continue
c
       else if (typb.eq.4) then
c
c4134 case B(pq,rs) <-- A(pq,s)
c@!   oprav to tak ako v typext 1 a 2
c
c     tsets
c
       if (typa.ne.1) then
c     RC=12: nindB=4, nindeext=1, typext=3, typB=4 (typA is not 1, Stup)
       rc=12
       return
       end if
c
       do 4134 ia=1,mapda(0,5)
c
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
       sa3=mapda(ia,5)
c
       ib=mapib(sa1,sa2,ssu)
       ibm=mapib(sa1,sa2,sa3)
c
c     def lenght
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 4134
c
c     def possA
       nhelp2=mapda(ia,1)
c
c     def dimp,dimq,dimr,dims
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
       nhelp7=dimm(mapdb(0,4),mapdb(ib,6))
c
c     def fictive dimensions
       if (sa1.eq.sa2) then
       nhelp8=nhelp4*(nhelp4-1)/2
       else
       nhelp8=nhelp4*nhelp5
       end if
c
       if (ssu.gt.sa3) then
c
c     def possB
       nhelp3=mapdb(ib,1)
       call cct3_add32 (wrk(nhelp2),wrk(nhelp3),
     &  u,nhelp8,nhelp6,nhelp7,fact)
c
       else if (ssu.eq.sa3) then
c
c     def possB
       nhelp3=mapdb(ib,1)
c     def fictive dimensions
       nhelp9=nhelp6*(nhelp6-1)/2
       call cct3_add43 (wrk(nhelp2),wrk(nhelp3),
     &  u,nhelp8,nhelp9,nhelp6,fact)
c
       else
c     ssu<sa3  B(pq,sr) <-- -A_r (pq,s)
c     def possB-
       nhelp3=mapdb(ibm,1)
c     def fictive dimension
       nhelp9=nhelp8*nhelp7
       call cct3_add22 (wrk(nhelp2),wrk(nhelp3),u,nhelp9,nhelp6,-fact)
c
       end if
c
 4134   continue
c
       else
c     RC=13: nindB=4, nindext=1, typext=3, (typA is not 0 or 4, (NCI))
       rc=13
       return
       end if
c
       else if (typext.eq.4) then
c
       if (typb.eq.0) then
c
c4140 case B(p,q,r,s) <-- A(p,q,r)
c
c     tsets
c
       if (typa.ne.0) then
c     RC=14: nindB=4, nindeext=1, typext=4, typB=0 (typA is not 0, Stup)
       rc=14
       return
       end if
c
       do 4140 ia=1,mapda(0,5)
c
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
       sa3=mapda(ia,5)
c
       ib=mapib(sa1,sa2,sa3)
c
c     def lenght
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 4140
c
c     def possA,possB
       nhelp2=mapda(ia,1)
       nhelp3=mapdb(ib,1)
c
c     def dimp,dimq,dimr,dims
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
       nhelp7=dimm(mapdb(0,4),mapdb(ib,6))
c
c     def fictive dimensions
       nhelp8=nhelp4*nhelp5*nhelp6
c
       call cct3_add22 (wrk(nhelp2),wrk(nhelp3),u,nhelp8,nhelp7,fact)
c
 4140   continue
c

       else if (typb.eq.4) then
c
c4144 case B(pq,rs) <-- A(pq,r)
c@!   oprav to tak ako v typext 1 a 2
c
c     tsets
c
       if (typa.ne.1) then
c     RC=15: nindB=4, nindeext=1, typext=4, typB=4 (typA is not 1, Stup)
       rc=15
       return
       end if
c
       do 4144 ia=1,mapda(0,5)
c
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
       sa3=mapda(ia,5)
c
       ib=mapib(sa1,sa2,sa3)
       ibm=mapib(sa1,sa2,ssu)
c
c     def lenght
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 4144
c
c     def possA
       nhelp2=mapda(ia,1)
c
c     def dimp,dimq,dimr,dims
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
       nhelp7=dimm(mapdb(0,4),mapdb(ib,6))
c
c     def fictive dimensions
       if (sa1.eq.sa2) then
       nhelp8=nhelp4*(nhelp4-1)/2
       else
       nhelp8=nhelp4*nhelp5
       end if
c
       if (sa3.gt.ssu) then
c
c     def possB
       nhelp3=mapdb(ib,1)
c     def fictive dimension
       nhelp9=nhelp8*nhelp6
       call cct3_add22 (wrk(nhelp2),wrk(nhelp3),u,nhelp8,nhelp7,fact)
c
       else if (sa3.eq.ssu) then
c
c     def possB
       nhelp3=mapdb(ib,1)
c     def fictive dimensions
       nhelp9=nhelp6*(nhelp6-1)/2
       call cct3_add44 (wrk(nhelp2),
     &  wrk(nhelp3),u,nhelp8,nhelp9,nhelp6,fact)
c
       else
c     sa3<ssu  B(pq,sr) <-- -A_s (pq,r)
c     def possB-
       nhelp3=mapdb(ibm,1)
       call cct3_add32 (wrk(nhelp2),wrk(nhelp3),
     &  u,nhelp8,nhelp7,nhelp6,-fact)
c
       end if
c
 4144   continue
c
       else
c     RC=16: nindB=4, nindext=1, typext=4, (typA is not 0 or 4, (NCI))
       rc=16
       return
       end if
c
       else
c     RC=17: nindB=4, nindext=1, typext=@  (Stup)
       rc=17
       return
       end if
c
       else if (nindext.eq.2) then
c
       if (typext.eq.5) then
c
       if (typb.eq.0) then
c
c4250 case B(p,q,r,s) <-- A(r,s)
c
c     tests
c
       if ((typb.eq.0).and.(typa.ne.0)) then
c     RC=18: nindB=4, nindext=2, typext=5, typB=0 (typA is not 0, Stup)
       rc=18
       return
       end if
c
       do 4250 ia=1,mapda(0,5)
c
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
c@
       nhelp1=mmul(ssu,ssv)
       nhelp1=mmul(nhelp1,sa1)
       nhelp1=mmul(nhelp1,ssb)
       if (nhelp1.ne.sa2) then
       write(6,*) ' Add Bpqrs <- Ars incorrect',ssp,ssq,sa1,nhelp1,sa1,
     & sa2
       goto 4250
       end if
c@@
c
       ib=mapib(ssu,ssv,sa1)
c
c     def lenght
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 4250
c
c     def possA,possB
       nhelp2=mapda(ia,1)
       nhelp3=mapdb(ib,1)
c
c     def dimp,dimq,dimr,dims
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
       nhelp7=dimm(mapdb(0,4),mapdb(ib,6))
c
c     calc joined pq index
       pq=(v-1)*nhelp4+u
c
c     calc fictive lenghts
       nhelp9=nhelp6*nhelp7
       nhelp10=nhelp4*nhelp5
c
       call cct3_add21 (wrk(nhelp2),wrk(nhelp3),pq,nhelp10,nhelp9,fact)
c
 4250   continue
c
       else if (typb.eq.4) then
c
c4254 case B(pq,rs) <-- A(rs)
c
c     tests
c
       if ((typb.eq.4).and.(typa.ne.1)) then
c     RC=19: nindB=4, nindext=2, typext=5, typB=4 (typA is not 1, Stup)
       rc=19
       return
       end if
c
       do 4254 ia=1,mapda(0,5)
c
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
c
c@
       nhelp1=mmul(ssp,ssq)
       nhelp1=mmul(nhelp1,sa1)
       nhelp1=mmul(nhelp1,ssb)
       if (nhelp1.ne.sa2) then
       write(6,*) ' Add Bpqrs <- Ars incorrect'
       goto 4254
       end if
c@@
       ib=mapib(ssp,ssq,sa1)
c
c     def lenght
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 4254
c
c     def possA,possB
       nhelp2=mapda(ia,1)
       nhelp3=mapdb(ib,1)
c
c     def dimp,dimq,dimr,dims
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
       nhelp7=dimm(mapdb(0,4),mapdb(ib,6))
c
c     calc joined pq index and fictive lenght of pq pair
       if (ssp.eq.ssq) then
       pq=(p-1)*(p-2)/2+q
       nhelp10=nhelp4*(nhelp4-1)/2
       else
       pq=(q-1)*nhelp4+p
       nhelp10=nhelp4*nhelp5
       end if
c
c     calc fictive lenghts
       if (sa1.eq.sa2) then
       nhelp9=nhelp6*(nhelp6-1)/2
       else
       nhelp9=nhelp6*nhelp7
       end if
c
       call cct3_add21 (wrk(nhelp2),wrk(nhelp3),pq,nhelp10,nhelp9,fact)
c
 4254   continue
c
       else
c     RC=20: nindB=4, nindext=2, typext=5 (typB is not 0 or 4, NCI)
       rc=20
       return
       end if
c
       else if (typext.eq.6) then
c
c426  case B(p,q,r,s) <-- A(p,q) and B(pq,rs) <-- A(pq)
c
c     RC=21: nindB=4, nindext=2, typext=6, NCI)
       rc=21
       return
c
       else
c     RC=22: nindB=4, nindext=2, (typext is not 5 or 6, NCI)
       rc=22
       return
       end if
c
       else
c     RC=23: nindB=4, nindext>2 (NCI)
       rc=23
       return
       end if
c
       else if (nindb.eq.3) then
c
c     **********  -> B(pqr) **********
c
       if (nindext.eq.0) then
c
c300  case B(pqr) <-- A(pqr)
c
c     tests
c
       if (typa.ne.typb) then
c     RC=24: nindB=3, nindext=0 (TypA incompatible with TypB ,Stup)
       rc=24
       return
       end if
c
       do 300 ia=1,mapda(0,5)
c
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
       sa3=mapda(ia,5)
c
       ib=mapib(sa1,sa2,1)
c
c     def lenght
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 300
c
c     def possA,possB
       nhelp2=mapda(ia,1)
       nhelp3=mapdb(ib,1)
c
       call cct3_add10 (wrk(nhelp2),wrk(nhelp3),nhelp1,fact)
c
 300    continue
c
       else if (nindext.eq.1) then
c
       if (typext.eq.1) then
c
c311  case B(p,q,r) <-- A(q,r)
c
       if ((typa.eq.0).and.(typb.eq.0)) then
c
c311  case B(p,q,r) <-- A(q,r)
c
       do 311 ia=1,mapda(0,5)
c
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
c
       ib=mapib(ssu,sa1,1)
c
c     def lenght
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 311
c
c     def possA,possB
       nhelp2=mapda(ia,1)
       nhelp3=mapdb(ib,1)
c
c     def dimp,dimq,dimr
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
c
c     def fictive dimensions
       nhelp7=nhelp5*nhelp6
c
       call cct3_add21 (wrk(nhelp2),wrk(nhelp3),u,nhelp4,nhelp7,fact)
c
 311    continue
c
       else
c     RC=25: nindB=3, nindext=1, typext=1 (tybA,B is not 0, NCI)
       rc=25
       return
       end if
c
       else if (typext.eq.2) then
c
       if ((typa.eq.0).and.(typb.eq.0)) then
c
c312  case B(p,q,r) <-- A(p,r)
c
       do 312 ia=1,mapda(0,5)
c
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
c
       ib=mapib(sa1,ssu,1)
c
c     def lenght
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 312
c
c     def possA,possB
       nhelp2=mapda(ia,1)
       nhelp3=mapdb(ib,1)
c
c     def dimp,dimq,dimr
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
c
       call cct3_add32 (wrk(nhelp2),wrk(nhelp3),
     &  u,nhelp4,nhelp6,nhelp7,fact)
c
 312    continue
c
       else
c     RC=26: nindB=3, nindext=1, typext=2 (tybA,B is not 0, NCI)
       rc=26
       return
       end if
c
       else if (typext.eq.3) then
c
       if ((typa.eq.0).and.(typb.eq.0)) then
c
c313  case B(p,q,r) <-- A(p,q)
c
       do 313 ia=1,mapda(0,5)
c
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
c
       ib=mapib(sa1,sa2,1)
c
c     def lenght
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 313
c
c     def possA,possB
       nhelp2=mapda(ia,1)
       nhelp3=mapdb(ib,1)
c
c     def dimp,dimq,dimr
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
       nhelp6=dimm(mapdb(0,3),mapdb(ib,5))
c
c     def fictive dimensions
       nhelp7=nhelp4*nhelp5
c
       call cct3_add21 (wrk(nhelp2),wrk(nhelp3),u,nhelp7,nhelp6,fact)
c
 313    continue
c
       else
c     RC=27: nindB=3, nindext=1, typext=3 (tybA,B is not 0, NCI)
       rc=27
       return
       end if
c
       else
c     RC=28: nindB=3 , typext=@ (Stup)
       rc=28
       return
       end if
c
       else
c     RC=29: nindB=3, nindext>1 (NCI)
       rc=29
       return
       end if
c
       else if (nindb.eq.2) then
c
c     **********  -> B(pq) **********
c
       if (nindext.eq.0) then
c
c200  case B(p,q) <-- A(p,q)
c
       do 200 ia=1,mapda(0,5)
c
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
c
       ib=mapib(sa1,1,1)
c
c     def lenght
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 200
c
c     def possA,possB
       nhelp2=mapda(ia,1)
       nhelp3=mapdb(ib,1)
c
       call cct3_add10 (wrk(nhelp2),wrk(nhelp3),nhelp1,fact)
c
 200    continue
c
       else if (nindext.eq.1) then
c
       if (typext.eq.1) then
c
c211  case B(p,q) <-- A(q)
c
       if ((typa.eq.0).and.(typb.eq.0)) then
c
       do 211 ia=1,mapda(0,5)
c
       ib=mapib(ssu,1,1)
c
c     def lenght
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 211
c
c     def possA,possB
       nhelp2=mapda(ia,1)
       nhelp3=mapdb(ib,1)
c
c     def dimp,dimq
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
c
       call cct3_add21 (wrk(nhelp2),wrk(nhelp3),u,nhelp4,nhelp5,fact)
c
 211    continue
c
       else
c     RC=30: nindB=2, nindext=1, typext=1, (typA,B is not 0, NCI)
       rc=30
       return
       end if
c
       else if (typext.eq.2) then
c
c212  case B(p,q) <-- A(p)
c
       if ((typa.eq.0).and.(typb.eq.0)) then
c
       do 212 ia=1,mapda(0,5)
c
       sa1=mapda(ia,3)
       ib=mapib(sa1,1,1)
c
c     def lenght
       nhelp1=mapda(ia,2)
       if (nhelp1.eq.0) goto 212
c
c     def possA,possB
       nhelp2=mapda(ia,1)
       nhelp3=mapdb(ib,1)
c
c     def dimp,dimq
       nhelp4=dimm(mapdb(0,1),mapdb(ib,3))
       nhelp5=dimm(mapdb(0,2),mapdb(ib,4))
c
       call cct3_add22 (wrk(nhelp2),wrk(nhelp3),u,nhelp4,nhelp5,fact)
c
 212    continue
c
       else
c     RC=31: nindB=2, nindext=1, typext=2, (typA,B is not 0, NCI)
       rc=31
       return
       end if
c
       else
c     RC=32: nindB=2, nindext=1, typext=@ (Stup)
       rc=32
       return
       end if
c
       else
c     RC=33: nindB=2, ininext>1 (NCI)
       rc=33
       return
       end if
c
       else
c     RC=34: nindb less then 2 (NCI/Stup)
       rc=34
       return
       end if
c
       return
       end
c
c     -------------------------
c
       subroutine cct3_add10 (a,b,dimp,fact)
c
c     this routine do:
c     B(p) <-- fact * A(p)
c
       integer dimp
       real*8 fact
       real*8 b(1:dimp)
       real*8 a(1:dimp)
c
c     help variable
c
       integer p
c
       do 100 p=1,dimp
       b(p)=b(p)+fact*a(p)
 100    continue
c
       return
       end
c
c     -------------------------
c
       subroutine cct3_add21 (a,b,p,dimp,dimq,fact)
c
c     this routine do:
c     B(p,q) <-- fact * A(q) for given p
c
       integer dimp,dimq,p
       real*8 fact
       real*8 b(1:dimp,1:dimq)
       real*8 a(1:dimq)
c
c     help variable
c
       integer q
c
       do 100 q=1,dimq
       b(p,q)=b(p,q)+fact*a(q)
 100    continue
c
       return
       end
c
c     -------------------------
c
       subroutine cct3_add22 (a,b,q,dimp,dimq,fact)
c
c     this routine do:
c     B(p,q) <-- fact * A(p) for given q
c
       integer dimp,dimq,q
       real*8 fact
       real*8 b(1:dimp,1:dimq)
       real*8 a(1:dimp)
c
c     help variable
c
       integer p
c
       do 100 p=1,dimp
       b(p,q)=b(p,q)+fact*a(p)
 100    continue
c
       return
       end
c
c     -------------------------
c
       subroutine cct3_add32 (a,b,q,dimp,dimq,dimr,fact)
c
c     this routine do:
c     B(p,q,r) <-- fact * A(p,r) for given q
c
       integer dimp,dimq,dimr,q
       real*8 fact
       real*8 b(1:dimp,1:dimq,1:dimr)
       real*8 a(1:dimp,1:dimr)
c
c     help variable
c
       integer p,r
c
       do 100 r=1,dimr
       do 101 p=1,dimp
       b(p,q,r)=b(p,q,r)+fact*a(p,r)
 101    continue
 100    continue
c
       return
       end
c
c     -------------------------
c
       subroutine cct3_add41 (a,b,p,dimp,dimpq,dimr,fact)

c     this routine do:
c     B(pq,r) <-- fact * A(q,r) for given p
c
#include "t31.fh"
       integer dimp,dimpq,dimr,p
       real*8 fact
       real*8 b(1:dimpq,1:dimr)
       real*8 a(1:dimp,1:dimr)
c
c     help variable
c
       integer q,r,pq,qp
c
       if (p.eq.1) goto 101
c
       do 100 r=1,dimr
       pq=nshf(p)
c
       do 50 q=1,p-1
       pq=pq+1
       b(pq,r)=b(pq,r)+fact*a(q,r)
 50     continue
c
 100    continue
c
 101    if (p.eq.dimp) then
       return
       end if
c
       do 200 r=1,dimr
c
       do 150 q=p+1,dimp
       qp=nshf(q)+p
       b(qp,r)=b(qp,r)-fact*a(q,r)
 150    continue
c
 200    continue
c
       return
       end
c
c     -------------------------
c
       subroutine cct3_add42 (a,b,q,dimq,dimpq,dimr,fact)

c     this routine do:
c     B(pq,r) <-- fact * A(p,r) for given q
c
#include "t31.fh"
       integer dimq,dimpq,dimr,q
       real*8 fact
       real*8 b(1:dimpq,1:dimr)
       real*8 a(1:dimq,1:dimr)
c
c     help variable
c
       integer pq,qp,r,p
c
       if (q.eq.1) goto 101
c
       do 100 r=1,dimr
       qp=nshf(q)
c
       do 50 p=1,q-1
       qp=qp+1
       b(qp,r)=b(qp,r)-fact*a(p,r)
 50     continue
c
 100    continue
c
 101    if (q.eq.dimq) then
       return
       end if
c
       do 200 r=1,dimr
c
       do 150 p=q+1,dimq
       pq=nshf(p)+q
       b(pq,r)=b(pq,r)+fact*a(p,r)
 150    continue
c
 200    continue
c
       return
       end
c
c     -------------------------
c
       subroutine cct3_add43 (a,b,q,dimp,dimqr,dimr,fact)

c     this routine do:
c     B(p,qr) <-- fact * A(p,r) for given q
c
#include "t31.fh"
       integer dimp,dimqr,dimr,q
       real*8 fact
       real*8 b(1:dimp,1:dimqr)
       real*8 a(1:dimp,1:dimr)
c
c     help variable
c
       integer p,qr,rq,r
c
       if (q.eq.1) goto 101
c
       qr=nshf(q)
       do 100 r=1,q-1
       qr=qr+1
c
       do 50 p=1,dimp
       b(p,qr)=b(p,qr)+fact*a(p,r)
 50     continue
c
 100    continue
c
 101    if (q.eq.dimr) then
       return
       end if
c
c
       do 200 r=q+1,dimr
       rq=nshf(r)+q
       do 150 p=1,dimp
       b(p,rq)=b(p,rq)-fact*a(p,r)
 150    continue
c
 200    continue
c
       return
       end
c
c     -------------------------
c
       subroutine cct3_add44 (a,b,r,dimp,dimqr,dimq,fact)

c     this routine do:
c     B(p,qr) <-- fact * A(p,q) for given r
c
#include "t31.fh"
       integer dimp,dimqr,dimq,r
       real*8 fact
       real*8 b(1:dimp,1:dimqr)
       real*8 a(1:dimp,1:dimq)
c
c     help variable
c
       integer p,qr,rq,q
c
       if (r.eq.1) goto 101
c
       rq=nshf(r)
       do 100 q=1,r-1
       rq=rq+1
c
       do 50 p=1,dimp
       b(p,rq)=b(p,rq)-fact*a(p,q)
 50     continue
c
 100    continue
c
 101    if (r.eq.dimq) then
       return
       end if
c
c
       do 200 q=r+1,dimq
       qr=nshf(q)+r
       do 150 p=1,dimp
       b(p,qr)=b(p,qr)+fact*a(p,q)
 150    continue
c
 200    continue
c
       return
       end
c

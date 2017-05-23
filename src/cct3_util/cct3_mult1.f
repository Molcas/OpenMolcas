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
c     mult
c     multc0
c     multy0
c
c     -------------------------------------------------------
c
       subroutine cct3_mult (wrk,wrksize,
     & ninda,nindb,nindc,nindsum,
     & mapda,mapia,ssa,mapdb,mapib,ssb,mapdc,mapic,ssc,
     & possc0,rc)
c
c     ninda   - # of indexes in matrix A (Input)
c     nindb   - # of indexes in matrix B (Input)
c     nindc   - # of indexes in matrix C (Input, for test)
c     nindsum - # of sumation indexes    (Input)
c     mapda   - direct map matrix corresponding to A  (Input)
c     mapia   - inverse map matrix corresponding to A  (Input)
c     ssa     - overall symmetry state of matrix A (Input)
c     mapdb   - direct map matrix corresponding to B  (Input)
c     mapib   - inverse map matrix corresponding to B  (Input)
c     ssb     - overall symmetry state of matrix B (Input)
c     mapdc   - direct map matrix corresponding to C  (Output)
c     mapic   - inverse map matrix corresponding to C  (Output)
c     ssc     - overall symmetry state of matrix C  (Output)
c     possc0  - initial possition of matrix C  (Input)
c
c
c     This routine realize matrix-matrix and matrix-vector multiplications
c     A(indA)*B(indB)=C(indC) or A(indA)*B(indB)=Y(indC)
c
c     N.B.
c     ninda, nindb, and nindsum are independent,
c     nindc is tested if it is in agreement with previous ones - for test
c     typA, typB are tested if they are in mutual agreement
c
c
c     Table of implemented processes
c
c     nindA   nindB  nindsum =>nindC        Operation               Implemented
c     4       4       4       0     A(pqrs)*B(pqrs)  =S             Not yet
c     4       4       3       2     A(p,qrs)*B(qrs,t)=C(p,t)          Yes
c     4       4       2       4     A(pq,rs)*B(rs,tu)=C(pq,tu)        Yes
c     4       4       1       6     A(pqr,s)*B(s,tuv)=C(pqr,tuv)      No
c     4       4       0       8     A(pqrs)*B(tuvz)=C(pqrs,tuvz)      No
c
c     4       3       3       1     A(p,qrs)*B(qrs)  =Y(p)            Yes
c     4       3       2       3     A(pq,rs)*B(rs,t) =C(pq,t)         Yes
c     4       3       1       5     A(pqr,s)*B(s,tu) =C(pqr,tu)       No
c     4       3       0       7     A(pqrs)*B(tuv)   =C(pqrs,tuv)     No
c
c     4       2       2       2     A(pq,rs)*B(pq)   =Y(pq)           Yes
c     4       2       1       4     A(pqr,s)*B(s,t)  =C(pqr,t)        Yes
c     4       2       0       6     A(pqrs)*B(tu)    =C(pqrs,tu)      No
c
c     4       1       1       3     A(pqr,s)*B(s)    =Y(pqr)        Not yet
c     4       1       0       5     A(pqrs)*B(t)     =C(pqrs,t)       No
c
c     3       4       3       1     A(pqr) *B(pqr,s) =Y(s)          Not need
c     3       4       2       3     A(p,qr)*B(qr,st) =C(p,st)         Yes
c     3       4       1       5     A(pq,r)*B(r,stu) =C(pq,stu)       No
c     3       4       0       7     A(pqr)*B(stuv) =C(pqr,stuv)       No
c
c     3       3       3       0     A(pqr) *B(pqr)   =S             Not yet
c     3       3       2       2     A(p,qr)*B(qr,s)  =C(p,s)          Yes
c     3       3       1       4     A(pq,r)*B(r,st)  =C(pq,st)        Yes
c     3       3       0       6     A(pqr) *B(stu)   =C(pqrmstu)      No
c
c     3       2       2       1     A(p,qr)*B(qr)    =Y(p)            Yes
c     3       2       1       3     A(pq,r)*B(r,s)   =C(pq,s)         Yes
c     3       2       0       5     A(pqr) *B(st)    =C(pqr,st)       No
c
c     3       1       1       2     A(pq,r)*B(r)     =Y(pq)         Not yet
c     3       1       0       4     A(pqr) *B(s)     =C(pqr,s)        No
c
c     2       4       2       2     A(pq) *B(pq,rs)  =Y(rs)         Not need
c     2       4       1       4     A(p,q)*B(q,rst)  =C(p,rst)        Yes
c     2       4       0       6     A(pq) *B(rstu)   =C(pq,rstu)      No
c
c     2       3       2       1     A(pq) *B(pq,r)   =Y(r)          Not need
c     2       3       1       3     A(p,q)*B(q,rs)   =C(p,rs)         Yes
c     2       3       0       5     A(pq) *B(rst)    =C(pq,rst)       No
c
c     2       2       2       0     A(pq) *B(pq)     =S             Not yet
c     2       2       1       2     A(p,q)*B(q,r)    =C(p,r)          Yes
c     2       2       0       4     A(pq) *B(rs)     =C(pq,rs)        No
c
c     2       1       1       1     A(p,q)*B(q)      =Y(p)          Not yet
c     2       1       0       3     A(pq) *B(r)      =C(pq,r)         No
c
c     1       4       1       3     A(p)*B(p,qrs)    =Y(qrs)        Not need
c     1       4       0       5     A(p)*B(qrst)     =C(p,qrst)       No
c
c     1       3       1       2     A(p)*B(p,qr)     =Y(qr)         Not need
c     1       3       0       4     A(p)*B(qrs)      =C(p,qrs)        No
c
c     1       2       1       1     A(p)*B(p,q)      =Y(q)          Not need
c     1       2       0       3     A(p)*B(qr)       =C(p,qr)         No
c
c     1       1       1       0     A(p)*B(p)        =S             Not yet
c     1       1       0       2     A(p)*B(q)        =C(p,q)          No
c
c     Legend used in error messeges:
c     NCI = Not Currently Implemented
c     Stup= Stupidity
c     @   = inproper value
c
c
#include "t31.fh"
#include "wrk.fh"
c
       integer ninda,nindb,nindc,nindsum,ssa,ssb,ssc,possc0,rc
c
       integer mapda(0:512,1:6)
       integer mapdb(0:512,1:6)
       integer mapdc(0:512,1:6)
c
       integer mapia(1:8,1:8,1:8)
       integer mapib(1:8,1:8,1:8)
       integer mapic(1:8,1:8,1:8)
c
c     help variables
c
       integer mvec(1:4096,1:7)
       integer typa,typb,ix

       rc=0
       ssc=mmul(ssa,ssb)
       typa=mapda(0,6)
       typb=mapdb(0,6)
c
c
c
       if (ninda.eq.4) then
c
c     ****** nind A = 4  ******
c
       if (nindb.eq.4) then
c
c44   *** nind A =4, B=4
c
       if (nindsum.eq.1) then
c     RC=1 : nindA=4, nindB=4, nindsum=1 (too large metiate, NCI)
       rc=1
       return
c
       else if (nindsum.eq.2) then
c
c     case A(pqrs)*B(rstu)=C(pqtu)
c
c     tests
c
       if (nindc.ne.4) then
c     RC=2 : nindA=4, nindB=4, nindsum=2 (# of ind. in C is not 4, Stup)
       rc=2
       return
       end if
c
       if (typa.eq.2) then
c     RC=3 : nindA=4, nindB=4, nindsum=2 (typA is 2, Stup)
       rc=3
       return
       end if
c
       if (typb.eq.2) then
c     RC=4 : nindA=4, nindB=4, nindsum=2 (typB is 2, Stup)
       rc=4
       return
       end if
c
       if ((typa.eq.3).or.(typa.eq.4)) then
c     there is  r>s in A
       if ((typb.eq.4).or.(typb.eq.1)) then
c     OK
       else
c     RC=5 : nindA=4, nindB=4, nindsum=2 (typA incompatible with typB , Stup)
       rc=5
       return
       end if
       else
c     there is no r>s in A
       if ((typb.eq.4).or.(typb.eq.1)) then
c     RC=5 : nindA=4, nindB=4, nindsum=2 (typA incompatible with typB , Stup)
       rc=5
       return
       else
c     OK
       end if
       end if
c
c     call cct3_grc44C and multc0
c
       call cct3_grc44C (mapda,mapdb,mapdc,mapia,mapib,mapic,
     & mvec,ssa,ssb,2,possc0,ix)
c
       call cct3_multc0 (wrk,wrksize,
     & mvec,ix,mapdc,1)
c
       else if (nindsum.eq.3) then
c
c     case A(pqrs)*B(qrst)=C(pt)
c
c     tests
c
       if (nindc.ne.2) then
c     RC=6 : nindA=4, nindB=4, nindsum=3 (# of ind. in C is not 2, Stup)
       rc=6
       return
       end if
c
       if (typa.eq.1) then
c     RC=7 : nindA=4, nindB=4, nindsum=3 (typA is 1, Stup)
       rc=7
       return
       end if
c
       if (typb.eq.3) then
c     RC=8 : nindA=4, nindB=4, nindsum=3 (typB is 3, Stup)
       rc=8
       return
       end if
c
       if ((typa.eq.2).and.(typb.ne.1)) then
c     RC=9 : nindA=4, nindB=4, nindsum=3 (typA incomp. with typB, Stup)
       rc=9
       return
       end if
c
       if ((typa.eq.3).and.(typb.ne.2)) then
c     RC=9 : nindA=4, nindB=4, nindsum=3 (typA incomp. with typB, Stup)
c     like in previous case
       rc=9
       return
       end if
c
c     call cct3_grc44C and multc0
c
       call cct3_grc44C (mapda,mapdb,mapdc,mapia,mapib,mapic,
     & mvec,ssa,ssb,1,possc0,ix)
c
       call cct3_multc0 (wrk,wrksize,
     & mvec,ix,mapdc,1)
c
       else if (nindsum.eq.4) then
c     RC=10: nindA=4, nindB=4, nindsum=4 (NCI)
       rc=10
       return
c
       else
c     RC=11: nindA=4, nindB=4, nindsum=@ (Stup)
       rc=11
       return
       end if
c
       else if (nindb.eq.3) then
c
c43   *** nind A =4, B=3
c
       if (nindsum.eq.3) then
c
c     case A(pqrs)*B(qrs)=Y(p)
c
c     tests
c
       if (nindc.ne.1) then
c     RC=12: nindA=4, nindB=3, nindsum=3 (# of index in C is not 1, Stup)
       rc=12
       return
       end if
c
       if ((typa.eq.1).or.(typa.eq.4)) then
c     RC=13: nindA=4, nindB=3, nindsum=3 (typA is 1 or 4, Stup)
       rc=13
       return
       end if
c
       if ((typa.eq.2).and.(typb.ne.1)) then
c     RC=14: nindA=4, nindB=3, nindsum=3 (typA incomp. with typB, Stup)
       rc=14
       return
       end if
c
       if ((typa.eq.3).and.(typb.ne.2)) then
c     RC=14: nindA=4, nindB=3, nindsum=3 (typA incomp. with typB, Stup)
c     like in previous case
       rc=14
       return
       end if
c
c     call cct3_grc43y and multy0
c
       call cct3_grc43y (mapda,mapdb,mapdc,mapia,mapib,mapic,
     & mvec,ssa,ssb,possc0,ix)

       call cct3_multy0 (wrk,wrksize,
     & mvec,ix,mapdc,1)
c
       else if (nindsum.eq.2) then
c
c     case A(pqrs)*B(rst)=C(pqt)
c
c     tests
c
       if (nindc.ne.3) then
c     RC=15: nindA=4, nindB=3, nindsum=2 (# of index in C is not 3, Stup)
       rc=15
       return
       end if
c
       if (typa.eq.2) then
c     RC=16: nindA=4, nindB=3, nindsum=2 (typA is 2, Stup)
       rc=16
       return
       end if
c
       if (typb.eq.2) then
c     RC=17: nindA=4, nindB=3, nindsum=2 (typB is 2, Stup)
       rc=17
       return
       end if
c
       if (((typa.eq.3).or.(typa.eq.4)).and.(typb.ne.1)) then
c     RC=18: nindA=4, nindB=3, nindsum=2 (typA incomp. with typB, Stup)
       rc=18
       return
       end if
c
c     call cct3_grc43c and multc0
c
       call cct3_grc43c (mapda,mapdb,mapdc,mapia,mapib,mapic,
     & mvec,ssa,ssb,2,possc0,ix)

       call cct3_multc0 (wrk,wrksize,
     & mvec,ix,mapdc,1)
c
       else if (nindsum.eq.1) then
c     RC=19: nindA=4, nindB=3, nindsum=1 (too large mediate, NCI)
       rc=19
       return
c
       else
c     RC=20: nindA=4, nindB=3, nindsum=@ (Stup)
       rc=20
       return
       end if
c
       else if (nindb.eq.2) then
c
c42   *** nind A =4, B=2
c
       if (nindsum.eq.1) then
c
c     case A(pqr,s)*B(s,t)=C(pqr,t)
c
c     tests
c
       if (nindc.ne.4) then
c     RC=21: nindA=4, nindB=2, nindsum=1 (# of index in C is not 4, Stup)
       rc=21
       return
       end if
c
       if ((typa.eq.3).or.(typa.eq.4)) then
c     RC=22: nindA=4, nindB=2, nindsum=1 (typA is 3 or 4, Stup)
       rc=22
       return
       end if
c
       if (typb.eq.1) then
c     RC=23: nindA=4, nindB=2, nindsum=1 (typA is 1, Stup)
       rc=23
       return
       end if
c
c     call cct3_grc42c and multc0
c
       call cct3_grc42c (mapda,mapdb,mapdc,mapia,mapib,mapic,
     & mvec,ssa,ssb,3,possc0,ix)
c
       call cct3_multc0 (wrk,wrksize,
     & mvec,ix,mapdc,1)
c
       else if (nindsum.eq.2) then
c
c     case A(pqrs)*B(rs)=Y(pq)
c
c     tests
c
       if (nindc.ne.2) then
c     RC=24: nindA=4, nindB=2, nindsum=2 (# of index in C is not 2, Stup)
       rc=24
       return
       end if
c
       if (typa.eq.2) then
c     RC=25: nindA=4, nindB=2, nindsum=2 (typA is 2, Stup)
       rc=25
       return
       end if
c
       if (((typa.eq.3).or.(typa.eq.4)).and.(typb.ne.1)) then
c     RC=26: nindA=4, nindB=2, nindsum=2 (typA incomp. with typB, Stup)
       rc=26
       return
       end if
c
c     call cct3_grc42y and multy0
c
       call cct3_grc42y (mapda,mapdb,mapdc,mapia,mapib,mapic,
     & mvec,ssa,ssb,possc0,ix)
c
       call cct3_multy0 (wrk,wrksize,
     & mvec,ix,mapdc,1)
c
       else
c     RC=27: nindA=4, nindB=2, nindsum=@ (Stup)
       rc=27
       return
       end if
c
       else if (nindb.eq.1) then
c
c41   *** nind A =4, B=1
c
c     RC=28: nindA=4, nindB=1 (NCI)
       rc=28
       return
c
       else
c
c4@   *** nind A =4, B=@
c
c     RC=29  : nindA=4, nindB=@ (Stup)
       rc=29
       return
       end if
c
       else if (ninda.eq.3) then
c
c     ****** nind A = 3  ******
c
       if (nindb.eq.4) then
c
c34   *** nind A =3, B=4
c
       if (nindsum.eq.1) then
c     RC=30: nindA=3, nindB=4, nindsum=1 (too large mediate, NCI)
       rc=30
       return
c
       else if (nindsum.eq.2) then
c
c     case A(p,qr)*B(qr,st)=C(p,st)
c
c     tests
c
       if (nindc.ne.3) then
c     RC=31: nindA=3, nindB=4, nindsum=2 (# of ind. in C is not 3, Stup)
       rc=31
       return
       end if
c
       if (typa.eq.1) then
c     RC=32: nindA=3, nindB=4, nindsum=2 (typA is 1, Stup)
       rc=32
       return
       end if
c
       if (typb.eq.2) then
c     RC=33: nindA=3, nindB=4, nindsum=2 (typB is 2, Stup)
       rc=33
       return
       end if
c
       if (((typb.eq.1).or.(typb.eq.4)).and.(typa.ne.2)) then
c     RC=34: nindA=3, nindB=4, nindsum=2 (typA incomp. with typB, Stup)
       rc=34
       return
       end if
c
c     call cct3_grc34c and multc0
c
       call cct3_grc34c (mapda,mapdb,mapdc,mapia,mapib,mapic,
     & mvec,ssa,ssb,1,possc0,ix)
c
       call cct3_multc0 (wrk,wrksize,
     & mvec,ix,mapdc,1)
c
       else if (nindsum.eq.3) then
c     RC=35: nindA=3, nindB=4, nindsum=3 (NCI)
       rc=35
       return
       end if
c
       else if (nindb.eq.3) then
c
c33   *** nind A =3, B=3
c
       if (nindsum.eq.1) then
c
c     case A(pq,r)*B(r,st)=A(pq,st)
c
c     tests
c
       if (nindc.ne.4) then
c     RC=36: nindA=3, nindB=3, nindsum=1 (# of ind. in C is not 4, Stup)
       rc=36
       return
       end if
c
       if (typa.eq.2) then
c     RC=37: nindA=3, nindB=3, nindsum=1 (typA is 2, Stup)
       rc=37
       return
       end if
c
       if (typb.eq.1) then
c     RC=38: nindA=3, nindB=3, nindsum=1 (typB is 1, Stup)
       rc=38
       return
       end if
c
c     call cct3_grc34c and multc0
c
       call cct3_grc34c (mapda,mapdb,mapdc,mapia,mapib,mapic,
     & mvec,ssa,ssb,2,possc0,ix)

       call cct3_multc0 (wrk,wrksize,
     & mvec,ix,mapdc,1)
c
       else if (nindsum.eq.2) then
c
c     case A(p,qr)*B(rq,s)=C(p,s)
c
c     tests
c
       if (nindc.ne.2) then
c     RC=39: nindA=3, nindB=3, nindsum=2 (# of ind. in C is not 2, Stup)
       rc=39
       return
       end if
c
       if (typa.eq.1) then
c     RC=40: nindA=3, nindB=3, nindsum=2 (typA is 1, Stup)
       rc=40
       return
       end if
c
       if (typb.eq.2) then
c     RC=41: nindA=3, nindB=3, nindsum=2 (typB is 2, Stup)
       rc=41
       return
       end if
c
       if ((typa.eq.2).and.(typb.ne.1)) then
c     RC=42: nindA=3, nindB=3, nindsum=2 (typA incomp. with typB, Stup)
       rc=42
       return
       end if
c
c     call cct3_grc33c and multc0
c
       call cct3_grc33c (mapda,mapdb,mapdc,mapia,mapib,mapic,
     & mvec,ssa,ssb,1,possc0,ix)

       call cct3_multc0 (wrk,wrksize,
     & mvec,ix,mapdc,1)
c
       else if (nindsum.eq.3) then
c     RC=43: nindA=3, nindB=3, nindsum=3 (NCI)
       rc=43
       return
       end if
c
       else if (nindb.eq.2) then
c
c32   *** nind A =3, B=2
c
       if (nindsum.eq.1) then
c
c     case A(pq,r)*B(r,s)=C(pq,s)
c
c     tests
c
       if (nindc.ne.3) then
c     RC=44: nindA=3, nindB=2, nindsum=1 (# of ind. in C is not 3, Stup)
       rc=44
       return
       end if
c
       if (typa.eq.2) then
c     RC=45: nindA=3, nindB=2, nindsum=1 (typA is 2, Stup)
       rc=45
       return
       end if
c
       if (typb.eq.1) then
c     RC=46: nindA=3, nindB=2, nindsum=1 (typB is 1, Stup)
       rc=46
       return
       end if
c
c     call cct3_grc32c and multc0
c
       call cct3_grc32c (mapda,mapdb,mapdc,mapia,mapib,mapic,
     & mvec,ssa,ssb,2,possc0,ix)

       call cct3_multc0 (wrk,wrksize,
     & mvec,ix,mapdc,1)
c
       else if (nindsum.eq.2) then
c
c     case A(p,qr)*B(qr)=C(p)
c
c     tests
c
       if (nindc.ne.1) then
c     RC=47: nindA=3, nindB=2, nindsum=2 (# of ind. in C is not 1, Stup)
       rc=47
       return
       end if
c
       if (typa.eq.1) then
c     RC=48: nindA=3, nindB=2, nindsum=2 (typA is 1, Stup)
       rc=48
       return
       end if
c
       if ((typb.eq.1).and.(typa.ne.2)) then
c     RC=49: nindA=3, nindB=2, nindsum=2 (typA incomp. with typB, Stup)
       rc=49
       return
       end if
c
c     call cct3_grc32y and multy0
c
       call cct3_grc32y (mapda,mapdb,mapdc,mapia,mapib,mapic,
     & mvec,ssa,ssb,possc0,ix)
c
       call cct3_multy0 (wrk,wrksize,
     & mvec,ix,mapdc,1)
c
       else
c     RC=50: nindA=3, nindB=2, nindsum=@ (Stup)
       rc=50
       return
       end if
c
       else if (nindb.eq.1) then
c
c31   *** nind A =3, B=1
c
c     RC=51: nindA=3, nindB=1 (NCI)
       rc=51
       return
c
       else
c
c3@   *** nind A =3, B=@
c
c     RC=52: nindA=3, nindB=@ (Stup)
       rc=52
       return
       end if
c
       else if (ninda.eq.2) then
c
c     ****** nind A = 2  ******
c
       if (nindb.eq.4) then
c
c24   *** nind A =2, B=4
c
       if (nindsum.eq.1) then
c
c     case A(p,q)*B(q,rst)=C(p,rst)
c
c     tests
c
       if (nindc.ne.4) then
c     RC=53: nindA=2, nindB=4, nindsum=1 (# of ind. in C is not 4, Stup)
       rc=53
       return
       end if
c
       if (typa.eq.1) then
c     RC=54: nindA=2, nindB=4, nindsum=1 (typA is 1, Stup)
       rc=54
       return
       end if
c
       if ((typb.eq.1).or.(typb.eq.4)) then
c     RC=55: nindA=2, nindB=4, nindsum=1 (typB is 1 or 4, Stup)
       rc=55
       return
       end if
c
c     call cct3_grc24c and multc0
c
       call cct3_grc24C (mapda,mapdb,mapdc,mapia,mapib,mapic,
     & mvec,ssa,ssb,1,possc0,ix)
c
       call cct3_multc0 (wrk,wrksize,
     & mvec,ix,mapdc,1)
c
       else if (nindsum.eq.2) then
c     RC=56: nindA=2, nindB=4, nindsum=2 (NCI)
       rc=56
       return
c
       else
c     RC=57: nindA=2, nindB=4, nindsum=@ (Stup)
       rc=57
       return
       end if
c
       else if (nindb.eq.3) then
c
c23   *** nind A =2, B=3
c
       if (nindsum.eq.1) then
c
c     case A(p,q)*B(q,rs)=C(p,rs)
c
c     tests
c
       if (nindc.ne.3) then
c     RC=58: nindA=2, nindB=3, nindsum=1 (# of ind. in C is not 3, Stup)
       rc=58
       return
       end if
c
       if (typa.eq.1) then
c     RC=59: nindA=2, nindB=3, nindsum=1 (typA is 1, Stup)
       rc=59
       return
       end if
c
       if (typb.eq.1) then
c     RC=60: nindA=2, nindB=3, nindsum=1 (typB is 1, Stup)
       rc=60
       return
       end if
c
c     call cct3_grc23c and multc0
c
       call cct3_grc23C (mapda,mapdb,mapdc,mapia,mapib,mapic,
     & mvec,ssa,ssb,1,possc0,ix)

       call cct3_multc0 (wrk,wrksize,
     & mvec,ix,mapdc,1)
c
       else if (nindsum.eq.2) then
c     RC=61: nindA=2, nindB=3, nindsum=2 (NCI)
       rc=61
       return
c
       else
c     RC=62: nindA=2, nindB=3, nindsum=@ (Stup)
       rc=62
       return
       end if
c
       else if (nindb.eq.2) then
c
c22   *** nind A =2, B=2
c
       if (nindsum.eq.1) then
c
c     case A(p,q)*B(q,r)=C(p,r)
c
c     tests
c
       if (nindc.ne.2) then
c     RC=63: nindA=2, nindB=2, nindsum=1 (# of ind. in C is not 2, Stup)
       rc=63
       return
       end if
c
       if (typa.eq.1) then
c     RC=64: nindA=2, nindB=2, nindsum=1 (typA is 1, Stup)
       rc=64
       return
       end if
c
       if (typb.eq.1) then
c     RC=65: nindA=2, nindB=2, nindsum=1 (typB is 1, Stup)
       rc=65
       return
       end if
c
c     call cct3_grc22c and multc0
c
       call cct3_grc22C (mapda,mapdb,mapdc,mapia,mapib,mapic,
     & mvec,ssa,ssb,1,possc0,ix)
c
       call cct3_multc0 (wrk,wrksize,
     & mvec,ix,mapdc,1)
c
       else if (nindsum.eq.2) then
c     RC=66: nindA=2, nindB=2, nindsum=2 (NCI)
       rc=66
       return
c
       else
c     RC=67: nindA=2, nindB=2, nindsum=@ (Stup)
       rc=67
       return
       end if
c
       else if (nindb.eq.1) then
c
c21   *** nind A =2, B=1
c
c     RC=68: nindA=2, nindB=1 (NCI)
       rc=68
       return
c
       else
c
c2@   *** nind A =2, B=@
c
c     RC=69: nindA=2, nindB=@
       rc=69
       return
       end if
c
       else if (ninda.eq.1) then
c
c     ****** nind A = 1  ******
c
c1x   RC=70: nind=1 (NCI)
       rc=70
       return
c
       else
c
c     ****** nind A = @  ******
c
c@x   RC=71: inproper nindA
       rc=71
       return
       end if
c
c
       return
       end
c
c     ------------------
c
       subroutine cct3_multc0 (wrk,wrksize,
     & mvec,ix,mapdc,key)
c
c     This routine realize multiplying according mvec
c     for C=A*B
c     N.B. if key=0, C file is not vanished (ie can be used for
c     adding to some existing file)
c
c     If C=A*B process is faster or comparambe with C=AT*B then mchntyp should be set to 1.
c     If C=AT*B is significantly faster than C=A*T (more than 20%), than mchntyp should be set
c     to 2. (default is 1)
c     if mchntyp is 2, than
c     1) proceses with scale(A)/scale(B) > scalelim will be calculated as C=A*B
c     2) processes with scale(A)/scale(B) < scalelim will be calculated as C=AT*B
c     Note, that for mchntyp =2 more memory is required, due to requirement of
c     aditional o2v2 help file possd0        (parameter possd0 is transported through t31.fh, not
c     through ccsd2.fh)
c
c
#include "t31.fh"
#include "wrk.fh"
       integer mvec(1:4096,1:7)
       integer ix,key
       integer mapdc(0:512,1:6)

c     help variables
c
       integer nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,nhelp6
       integer iix,ic
       real*8 scale
c
c1    set C=0
c
       if (key.eq.1) then
c
c     C matrix must be vanished
c
       do 100 ic=1,mapdc(0,5)
       nhelp1=mapdc(ic,1)
       nhelp2=mapdc(ic,2)
       call cct3_mv0zero (nhelp2,nhelp2,wrk(nhelp1))
 100    continue
c
       end if
c
c2    C=C+A*B
c
       if (ix.eq.0) then
       return
       end if
c
       do 200 iix=1,ix
c
c     skip this sumation if yes/no=0
       if (mvec(iix,1).eq.0) goto 200
c
c     realize individial sumation
c
c     def possitions of A,B,C
       nhelp1=mvec(iix,2)
       nhelp2=mvec(iix,3)
       nhelp3=mvec(iix,4)
c
c     def rowA(rowC), colA(rowB,sum), colB(colC)
       nhelp4=mvec(iix,5)
       nhelp5=mvec(iix,6)
       nhelp6=mvec(iix,7)
c
       if (mchntyp.eq.1) then
c
c*    Typ 1
       call cct3_mc0c1a3b (nhelp4,nhelp5,nhelp5,nhelp6,nhelp4,nhelp6,
     & nhelp4,nhelp5,nhelp6,wrk(nhelp1),wrk(nhelp2),wrk(nhelp3))
c
       else
c
c*    Typ2
       scale=(1.0d0*nhelp4)/(1.0d0*nhelp6)
       if (scale.gt.slim) then
       call cct3_mc0c1a3b (nhelp4,nhelp5,nhelp5,nhelp6,nhelp4,nhelp6,
     & nhelp4,nhelp5,nhelp6,wrk(nhelp1),wrk(nhelp2),wrk(nhelp3))

       else
c     map D=AT
       call cct3_map21 (wrk(nhelp1),wrk(possd0),nhelp4,nhelp5,2,1,1)
c     calc C=DT*B
       call cct3_mc0c1at3b (nhelp5,nhelp4,nhelp5,nhelp6,nhelp4,nhelp6,
     & nhelp4,nhelp5,nhelp6,wrk(possd0),wrk(nhelp2),wrk(nhelp3))
       end if
c
       end if
c
 200    continue
c
       return
       end
c
c     -------------------------
c
       subroutine cct3_multy0 (wrk,wrksize,
     & mvec,ix,mapdy,key)
c
c     This routine realize multiplying according mvec
c     for Y=A*B
c     N.B. if key=0, Y file is not vanished (ie can be used for
c     adding to some existing file)
c
#include "t31.fh"
#include "wrk.fh"
       integer mvec(1:4096,1:7)
       integer ix,key
       integer mapdy(0:512,1:6)

c     help variables
c
       integer nhelp1,nhelp2,nhelp3,nhelp4,nhelp5
       integer iix,iy
c
c1    set C=0
c
       if (key.eq.1) then
c
c     Y vector must be vanished
c
       do 100 iy=1,mapdy(0,5)
       nhelp1=mapdy(iy,1)
       nhelp2=mapdy(iy,2)
       call cct3_mv0zero (nhelp2,nhelp2,wrk(nhelp1))
 100    continue
c
       end if
c
c2    Y=Y+A*B
c
       if (ix.eq.0) then
       return
       end if
c
       do 200 iix=1,ix
c
c     skip this sumation if yes/no=0
       if (mvec(iix,1).eq.0) goto 200
c
c     realize individial sumation
c
c     def possitions of A,B,Y
       nhelp1=mvec(iix,2)
       nhelp2=mvec(iix,3)
       nhelp3=mvec(iix,4)
c
c     def rowA(rowY), colA(sum)
       nhelp4=mvec(iix,5)
       nhelp5=mvec(iix,6)
c
       call cct3_mv0v1a3u (nhelp4,nhelp5,nhelp5,nhelp4,
     & nhelp4,nhelp5,1,1,
     & wrk(nhelp1),wrk(nhelp2),wrk(nhelp3))
c
 200    continue
c
       return
       end
c
c     -------------------------------
c

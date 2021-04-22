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
       subroutine fack (wrk,wrksize,
     & nind,newtyp,mapda,ssa,mapia,mapdb,mapib,possb0,
     &                  rc)
c
c     nind   - # of indexes in matrices A,B (Input)
c     newtyp - typ of final matrix B (Input) (i.e. mapdb(0,6), see docc.txt)
c     mapda  - direct map matrix corresponding to A (Input) (see docc.txt)
c     ssa    - overall symmetry state of matrix A (Input)
c     mapia  - inverse map matrix corresponding to A (Input) (see docc.txt)
c     mapdb  - direct map matrix corresponding to B (Output) (see docc.txt)
c     mapib  - inverse map matrix corresponding to B (Output) (see docc.txt)
c     possb0 - initial possition of B matrix in WRK (Input)
c     rc     - return (error) code (Output)
c
c     this routine realize following packings :

c     Operation                   Nind    TypA    TypB
c
c     B(pq)     = A(p,q)     - A(q,p)          2       0       1
c     B(pq,r)   = A(p,q,r)   - A(q,p,r)        3       0       1
c     B(p,qr)   = A(p,q,r)   - A(p,r,q)        3       0       2
c     B(pq,r,s) = A(p,q,r,s) - A(q,p,r,s)      4       0       1
c     B(p,q,rs) = A(p,q,r,s) - A(q,p,s,r)      4       0       3
c     B(pq,rs)  = NCI                          4       0       4
c     B(pq,rs)  = A(pq,r,s)  - A(pq,s,r)       4       1       4
c     B(pq,rs)  = A(p,q,rs)  - A(q,p,rs)       4       3       4
c
#include "ccsd1.fh"
#include "wrk.fh"
c
       integer nind,newtyp,ssa,possb0,rc
       integer mapda(0:512,1:6)
       integer mapdb(0:512,1:6)
       integer mapia(1:8,1:8,1:8)
       integer mapib(1:8,1:8,1:8)
c
c     help variables
c
       integer nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,nhelp6,nhelp7,nhelp8,
     & nhelp9
       integer sb1,sb2,sb3,sb4
       integer rc1,ib,iap,iam,posst
c
       rc=0
c
c     get mapdb,mapib
c
       call grc0 (nind,newtyp,mapda(0,1),mapda(0,2),mapda(0,3),mapda(0,
     &            4),ssa,
     & possb0,posst,mapdb,mapib)
c
       if (nind.lt.2) then
c
c     *********  less than 2 ind ***
c
c     RC=1 less than 2 indexes
       rc=1
       return
c
       else if (nind.eq.2) then
c
c     *********  2 indexes *********
c
c2.1  A(p,q) -> B(pq)
c
       if (mapda(0,6).ne.0) then
c     RC=2 bad type for: nind=2
       rc=2
       return
       end if
c
       do 210 ib=1,mapdb(0,5)
c
       sb1=mapdb(ib,3)
       sb2=mapdb(ib,4)
       if (mapdb(ib,2).eq.0) goto 210
c
       if (sb1.eq.sb2) then
c     sym b1 = sym b2
c
c     def ia+,-
       iap=mapia(sb1,1,1)
c
c     poss B,A+,A-
       nhelp1=mapdb(ib,1)
       nhelp2=mapda(iap,1)
c
c     dimp-s
       nhelp4=dimm(mapdb(0,1),sb1)
c
c     def size
       nhelp9=nhelp4*(nhelp4-1)/2
c
       call pack210 (wrk(nhelp2),wrk(nhelp1),nhelp9,nhelp4,rc1)
c
       else
c     sym b1 > sym b2
c
c     def ia+,-
       iap=mapia(sb1,1,1)
       iam=mapia(sb2,1,1)
c
c     poss B,A+,A-
       nhelp1=mapdb(ib,1)
       nhelp2=mapda(iap,1)
       nhelp3=mapda(iam,1)
c
c     dimp-s
       nhelp4=dimm(mapdb(0,1),sb1)
       nhelp5=dimm(mapdb(0,2),sb2)
c
       call pack211 (wrk(nhelp2),wrk(nhelp3),wrk(nhelp1),nhelp4,nhelp5,
     &               rc1)
c

       end if
c
 210    continue
c
       else if (nind.eq.3) then
c
c     *********  3 indexes *********
c
       if (mapda(0,6).ne.0) then
c     RC=3 bad type for: nind=3
       rc=3
       return
       end if
c
       if (newtyp.eq.1) then
c
c3.1  A(p,q,r) -> B(pq,r)
c
       do 310 ib=1,mapdb(0,5)
c
       sb1=mapdb(ib,3)
       sb2=mapdb(ib,4)
       sb3=mapdb(ib,5)
       if (mapdb(ib,2).eq.0) goto 310
c
       if (sb1.eq.sb2) then
c     sym b1 = sym b2
c
c     def ia+,-
       iap=mapia(sb1,sb2,1)
c
c     poss B,A+,A-
       nhelp1=mapdb(ib,1)
       nhelp2=mapda(iap,1)
c
c     dimp-s
       nhelp4=dimm(mapdb(0,1),sb1)
       nhelp6=dimm(mapdb(0,3),sb3)
c
c     def size
       nhelp9=nhelp4*(nhelp4-1)/2
c
       call pack310 (wrk(nhelp2),wrk(nhelp1),nhelp9,nhelp6,nhelp4,rc1)
c
       else
c     sym b1 > sym b2
c
c     def ia+,-
       iap=mapia(sb1,sb2,1)
       iam=mapia(sb2,sb1,1)
c
c     poss B,A+,A-
       nhelp1=mapdb(ib,1)
       nhelp2=mapda(iap,1)
       nhelp3=mapda(iam,1)
c
c     dimp-s
       nhelp4=dimm(mapdb(0,1),sb1)
       nhelp5=dimm(mapdb(0,2),sb2)
       nhelp6=dimm(mapdb(0,3),sb3)
c
       call pack311 (wrk(nhelp2),wrk(nhelp3),wrk(nhelp1),nhelp4,nhelp5,
     &               nhelp6,rc1)
c
       end if
c
 310    continue
c
       else if (newtyp.eq.2) then
c
c3.2  A(p,q,r) -> B(p,qr)
c
       do 320 ib=1,mapdb(0,5)
c
       sb1=mapdb(ib,3)
       sb2=mapdb(ib,4)
       sb3=mapdb(ib,5)
       if (mapdb(ib,2).eq.0) goto 320
c
       if (sb2.eq.sb3) then
c     sym b2 = sym b3
c
c     def ia+,-
       iap=mapia(sb1,sb2,1)
c
c     poss B,A+,A-
       nhelp1=mapdb(ib,1)
       nhelp2=mapda(iap,1)
c
c     dimp-s
       nhelp4=dimm(mapdb(0,1),sb1)
       nhelp5=dimm(mapdb(0,2),sb2)
c
c     def size
       nhelp9=nhelp5*(nhelp5-1)/2
c
       call pack320 (wrk(nhelp2),wrk(nhelp1),nhelp4,nhelp9,nhelp5,rc1)
c
       else
c     sym b2 > sym b3
c
c     def ia+,-
       iap=mapia(sb1,sb2,1)
       iam=mapia(sb1,sb3,1)
c
c     poss B,A+,A-
       nhelp1=mapdb(ib,1)
       nhelp2=mapda(iap,1)
       nhelp3=mapda(iam,1)
c
c     dimp-s
       nhelp4=dimm(mapdb(0,1),sb1)
       nhelp5=dimm(mapdb(0,2),sb2)
       nhelp6=dimm(mapdb(0,3),sb3)
c
       call pack321 (wrk(nhelp2),wrk(nhelp3),wrk(nhelp1),nhelp4,nhelp5,
     &               nhelp6,rc1)
c
       end if
c
 320    continue
c
       else
c     RC=4 : bad newtyp for: nind=3
       rc=4
       return
       end if
c
       else if (nind.eq.4) then
c
c     *********  4 indexes *********
c
       if (mapda(0,6).eq.0) then
c
c*    case A(p,q,r,s) => newtyp can be 1,3,4
c
       if (newtyp.eq.1) then
c
c4.1  A(p,q,r,s) -> B(pq,r,s)
c
       do 410 ib=1,mapdb(0,5)
c
       sb1=mapdb(ib,3)
       sb2=mapdb(ib,4)
       sb3=mapdb(ib,5)
       sb4=mapdb(ib,6)
       if (mapdb(ib,2).eq.0) goto 410
c
       if (sb1.eq.sb2) then
c     sym b1 = sym b2
c
c     def ia+,-
       iap=mapia(sb1,sb2,sb3)
c
c     poss B,A+,A-
       nhelp1=mapdb(ib,1)
       nhelp2=mapda(iap,1)
c
c     dimp-s
       nhelp4=dimm(mapdb(0,1),sb1)
       nhelp6=dimm(mapdb(0,3),sb3)
       nhelp7=dimm(mapdb(0,4),sb4)
c
c     def size
       nhelp8=nhelp4*(nhelp4-1)/2
       nhelp9=nhelp6*nhelp7
c
       call pack310 (wrk(nhelp2),wrk(nhelp1),nhelp8,nhelp9,nhelp4,rc1)
c
       else
c     sym b1 > sym b2
c
c     def ia+,-
       iap=mapia(sb1,sb2,sb3)
       iam=mapia(sb2,sb1,sb3)
c
c     poss B,A+,A-
       nhelp1=mapdb(ib,1)
       nhelp2=mapda(iap,1)
       nhelp3=mapda(iam,1)
c
c     dimp-s
       nhelp4=dimm(mapdb(0,1),sb1)
       nhelp5=dimm(mapdb(0,2),sb2)
       nhelp6=dimm(mapdb(0,3),sb3)
       nhelp7=dimm(mapdb(0,4),sb4)
c
c     def size
       nhelp9=nhelp6*nhelp7
c
       call pack311 (wrk(nhelp2),wrk(nhelp3),wrk(nhelp1),nhelp4,nhelp5,
     &               nhelp9,rc1)
c
       end if
c
 410    continue
c
       else if (newtyp.eq.3) then
c
c4.2  A(p,q,r,s) -> B(p,q,rs)
c
       do 420 ib=1,mapdb(0,5)
c
       sb1=mapdb(ib,3)
       sb2=mapdb(ib,4)
       sb3=mapdb(ib,5)
       sb4=mapdb(ib,6)
       if (mapdb(ib,2).eq.0) goto 420
c
       if (sb3.eq.sb4) then
c     sym b3 = sym b4
c
c     def ia+,-
       iap=mapia(sb1,sb2,sb3)
c
c     poss B,A+,A-
       nhelp1=mapdb(ib,1)
       nhelp2=mapda(iap,1)
c
c     dimp-s
       nhelp4=dimm(mapdb(0,1),sb1)
       nhelp5=dimm(mapdb(0,2),sb2)
       nhelp6=dimm(mapdb(0,3),sb3)
c
c     def size
       nhelp8=nhelp6*(nhelp6-1)/2
       nhelp9=nhelp4*nhelp5
c
       call pack320 (wrk(nhelp2),wrk(nhelp1),nhelp9,nhelp8,nhelp6,rc1)
c
       else
c     sym b3 > sym b4
c
c     def ia+,-
       iap=mapia(sb1,sb2,sb3)
       iam=mapia(sb1,sb2,sb4)
c
c     poss B,A+,A-
       nhelp1=mapdb(ib,1)
       nhelp2=mapda(iap,1)
       nhelp3=mapda(iam,1)
c
c     dimp-s
       nhelp4=dimm(mapdb(0,1),sb1)
       nhelp5=dimm(mapdb(0,2),sb2)
       nhelp6=dimm(mapdb(0,3),sb3)
       nhelp7=dimm(mapdb(0,4),sb4)
c
c     def size
       nhelp9=nhelp4*nhelp5
c
       call pack321 (wrk(nhelp2),wrk(nhelp3),wrk(nhelp1),nhelp9,nhelp6,
     &               nhelp7,rc1)
c
       end if
c
 420    continue
c
       else if (newtyp.eq.4) then
c
c4.3  A(p,q,r,s) -> B(pq,rs)
c
c     RC=5 : Not Currently Implemented : nind=4, oldtyp=0, newtyp=4
       rc=5
       return
c
       else
c     RC=6 : incompatible newtyp for: nind=4, oldtyp=0
       rc=6
       return
       end if
c
       else if (mapda(0,6).eq.1) then
c
c*    case A(pq,r,s) => newtyp can be 4
c
       if (newtyp.eq.4) then
c
c4.4  A(pq,r,s) -> B(pq,rs)
c
       do 440 ib=1,mapdb(0,5)
c
       sb1=mapdb(ib,3)
       sb2=mapdb(ib,4)
       sb3=mapdb(ib,5)
       sb4=mapdb(ib,6)
       if (mapdb(ib,2).eq.0) goto 440
c
       if (sb3.eq.sb4) then
c     sym b3 = sym b4
c
c     def ia+,-
       iap=mapia(sb1,sb2,sb3)
c
c     poss B,A+,A-
       nhelp1=mapdb(ib,1)
       nhelp2=mapda(iap,1)
c
c     dimp-s
       nhelp4=dimm(mapdb(0,1),sb1)
       nhelp5=dimm(mapdb(0,2),sb2)
       nhelp6=dimm(mapdb(0,3),sb3)
c
c     def size
       nhelp8=nhelp6*(nhelp6-1)/2
       if (sb1.eq.sb2) then
       nhelp9=nhelp4*(nhelp4-1)/2
       else
       nhelp9=nhelp4*nhelp5
       end if
c
       call pack320 (wrk(nhelp2),wrk(nhelp1),nhelp9,nhelp8,nhelp6,rc1)
c
       else
c     sym b3 > sym b4
c
c     def ia+,-
       iap=mapia(sb1,sb2,sb3)
       iam=mapia(sb1,sb2,sb4)
c
c     poss B,A+,A-
       nhelp1=mapdb(ib,1)
       nhelp2=mapda(iap,1)
       nhelp3=mapda(iam,1)
c
c     dimp-s
       nhelp4=dimm(mapdb(0,1),sb1)
       nhelp5=dimm(mapdb(0,2),sb2)
       nhelp6=dimm(mapdb(0,3),sb3)
       nhelp7=dimm(mapdb(0,4),sb4)
c
c     def size
       if (sb1.eq.sb2) then
       nhelp9=nhelp4*(nhelp4-1)/2
       else
       nhelp9=nhelp4*nhelp5
       end if
c
       call pack321 (wrk(nhelp2),wrk(nhelp3),wrk(nhelp1),nhelp9,nhelp6,
     &               nhelp7,rc1)
c
       end if
c
 440    continue
c
       else
c     RC=7 : incompatible newtyp for: nind=4, oldtyp=1
       rc=7
       return
       end if
c
       else if (mapda(0,6).eq.3) then
c
c*    case A(p,q,rs) => newtyp can be 4
c
       if (newtyp.eq.4) then
c
c4.5  A(p,q,rs) -> B(pq,rs)
c
       do 450 ib=1,mapdb(0,5)
c
       sb1=mapdb(ib,3)
       sb2=mapdb(ib,4)
       sb3=mapdb(ib,5)
       sb4=mapdb(ib,6)
       if (mapdb(ib,2).eq.0) goto 450
c
       if (sb1.eq.sb2) then
c     sym b1 = sym b2
c
c     def ia+,-
       iap=mapia(sb1,sb2,sb3)
c
c     poss B,A+,A-
       nhelp1=mapdb(ib,1)
       nhelp2=mapda(iap,1)
c
c     dimp-s
       nhelp4=dimm(mapdb(0,1),sb1)
       nhelp6=dimm(mapdb(0,3),sb3)
       nhelp7=dimm(mapdb(0,4),sb4)
c
c     def size
       nhelp8=nhelp4*(nhelp4-1)/2
       if (sb3.eq.sb4) then
       nhelp9=nhelp6*(nhelp6-1)/2
       else
       nhelp9=nhelp6*nhelp7
       end if
c
       call pack310 (wrk(nhelp2),wrk(nhelp1),nhelp8,nhelp9,nhelp4,rc1)
c
       else
c     sym b1 > sym b2
c
c     def ia+,-
       iap=mapia(sb1,sb2,sb3)
       iam=mapia(sb2,sb1,sb3)
c
c     poss B,A+,A-
       nhelp1=mapdb(ib,1)
       nhelp2=mapda(iap,1)
       nhelp3=mapda(iam,1)
c
c     dimp-s
       nhelp4=dimm(mapdb(0,1),sb1)
       nhelp5=dimm(mapdb(0,2),sb2)
       nhelp6=dimm(mapdb(0,3),sb3)
       nhelp7=dimm(mapdb(0,4),sb4)
c
c     def size
       if (sb3.eq.sb4) then
       nhelp9=nhelp6*(nhelp6-1)/2
       else
       nhelp9=nhelp6*nhelp7
       end if
c
       call pack311 (wrk(nhelp2),wrk(nhelp3),wrk(nhelp1),nhelp4,nhelp5,
     &               nhelp9,rc1)
c
       end if
c
 450    continue
c
       else
c     RC=8 : incompatible newtyp for: nind=4, oldtyp=3
       rc=8
       return
       end if
c
       else
c     RC=9 : incompatible typold for: nind=4
       rc=9
       return
       end if
c
       else
c
c     *********  more than 4 ind ***
c
c     RC=10: more than 4 indexes
       rc=10
       return
       end if
c
       return
       end

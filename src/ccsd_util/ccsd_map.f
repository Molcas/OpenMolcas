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
c     map
c     noperm
c     map41
c     map42
c     map31
c     map32
c     map21
c     map22
c     map11
c
c     -----------------------------------------------------------
c
       subroutine map (wrk,wrksize,
     & nind,p,q,r,s,
     & mapda,mapia,ssa,mapdb,mapib,possb0,posst,rc)
c
c     this routine realize mappings
c
c     B(indb) <-- A(inda)
c     where inda are order of indexes in mtx A and indb = Perm(inda)
c
c     nind   - number of indexes in matrix A (and B)  (Input)
c     p      - possition of 1-st index od mtx A in mtx B  (Input)
c     q      - possition of 2-nd index od mtx A in mtx B  (Input)
c     r      - possition of 3-rd index od mtx A in mtx B  (Input)
c     s      - possition of 4-th index od mtx A in mtx B  (Input)
c     mapda  - direct map matrix corresponding to A  (Input)
c     mapia  - inverse map matrix corresponding to A  (Input)
c     ssa    - overall symetry state  of matrix A  (Input)
c     mapdb  - direct map matrix corresponding to B  (Output)
c     mapib  - inverse map matrix corresponding to B  (Output)
c     possb0 - initial possition of matrix B in WRK  (Input)
c     posst  - final possition of matrix B in WRK (Output, not used yet)
c     rc     - return (error) code  (Output)
c
c
c     The table of implemented permutations
c
c     ninda  typA   p  q  r  s            Operation                 Implemented
c
c     4     0     all 24 comb.    A(1,2,3,4) -> B(p,q,r,s)             Yes
c
c     4     1     1  2  3  4      A(12,3,4)  -> B(12,3,4)              Yes
c     4     1     1  2  4  3      A(12,3,4)  -> B(12,4,3)              Yes
c     4     1     2  3  1  4      A(12,3,4)  -> B(3,12,4)              Yes
c     4     1           4  1      A(12,3,4)  -> B(4,12,3)              Yes
c     4     1     3  4  1  2      A(12,3,4)  -> B(3,4,12)              Yes
c     4     1           2  1      A(12,3,4)  -> B(4,3,12)              Yes
c     4     1     other comb.                                          No
c
c     4     2     1  2  3  4      A(1,23,4)  -> B(1,23,4)              Yes
c     4     2     4  2  3  1      A(1,23,4)  -> B(4,23,1)              Yes
c     4     2     3  1  2  4      A(1,23,4)  -> B(23,1,4)              Yes
c     4     2     4        3      A(1,23,4)  -> B(23,4,1)              Yes
c     4     2     1  3  4  2      A(1,23,4)  -> B(1,4,23)              Yes
c     4     2     2        1      A(1,23,4)  -> B(4,1,23)              Yes
c     4     2     other comb.                                          No
c
c     4     3     1  2  3  4      A(1,2,34)  -> B(1,2,34)              Yes
c     4     3     2  1  3  4      A(1,2,34)  -> B(2,1,34)              Yes
c     4     3     3  4  1  2      A(1,2,34)  -> B(34,1,2)              Yes
c     4     3     4  4            A(1,2,34)  -> B(34,2,1)              Yes
c     4     3     1  4  2  3      A(1,2,34)  -> B(1,34,2)              Yes
c     4     3     4  1            A(1,2,34)  -> B(2,34,1)              Yes
c     4     3     other comb.                                          No
c
c     4     4     1  2  3  4      A(12,34)   -> B(12,34)               Yes
c     4     4     3  4  1  2      A(12,34)   -> B(32,12)               Yes
c     4     4     other comb.                                          No
c
c     3     0     all 6 comb.     A(1,2,3)   -> B(p,q,r)               Yes
c
c     3     1     1  2  3  -      A(12,3)    -> B(12,3)                Yes
c     3     1     2  3  1  -      A(12,3)    -> B(3,12)                Yes
c     3     1     other comb.                                          No
c
c     3     2     1  2  3         A(1,23)    -> B(1,23)                Yes
c     3     2     3  1  2         A(1,23)    -> B(23,1)                Yes
c     3     2     other comb.                                          No
c
c     2     0     all 2 comb.     A(1,2)     -> B(p,q)                 Yes
c
c     2     1     1  2  -  -      A(12)      -> B(12)                  Yes
c     2     1     other comb.                                          No
c
c
#include "ccsd1.fh"
#include "wrk.fh"

c
       integer nind,p,q,r,s,ssa,possb0,posst,rc
       integer mapda(0:512,1:6),mapdb(0:512,1:6)
       integer mapia(1:8,1:8,1:8),mapib(1:8,1:8,1:8)
c
       integer type(1:4),dl(1:4),sa(1:4)
       integer nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,nhelp6,nhelp7
       integer ia,ib,typ,newtyp
c
       rc=0
c
c     some test of correctness
c
       nhelp1=p+q+r+s
       if ((nind.eq.1).and.(nhelp1.eq.1)) then
       else if ((nind.eq.2).and.(nhelp1.eq.3)) then
       else if ((nind.eq.3).and.(nhelp1.eq.6)) then
       else if ((nind.eq.4).and.(nhelp1.eq.10)) then
       else
c     RC=1 - uncorrect indexes
       rc=1
       return
       end if
c
c     ******** No permutation *******
c
       if (nind.eq.1) then
       call noperm (wrk,wrksize,
     & mapda,mapia,mapdb,mapib,possb0,posst)
       return
       end if
c
       if ((nind.eq.2).and.(p.eq.1).and.(q.eq.2)) then
       call noperm (wrk,wrksize,
     & mapda,mapia,mapdb,mapib,possb0,posst)
       return
       end if
c
       if ((nind.eq.3).and.(p.eq.1).and.(q.eq.2).and.(r.eq.3)) then
       call noperm (wrk,wrksize,
     & mapda,mapia,mapdb,mapib,possb0,posst)
       return
       end if
c
       if ((nind.eq.4).and.(p.eq.1).and.(q.eq.2).and.(r.eq.3)
     &     .and.(s.eq.4)) then
       call noperm (wrk,wrksize,
     & mapda,mapia,mapdb,mapib,possb0,posst)
       return
       end if
c
       if (nind.eq.2) then
c
c     *********** 2 index ***********
c
       typ=mapda(0,6)
c
       if (typ.eq.0) then
c
c2.1  map A(p,q) -> B(q,p)
c
c     get mapdb,mapib
c
       type(p)=mapda(0,1)
       type(q)=mapda(0,2)
       call grc0 (nind,0,type(1),type(2),0,0,ssa,
     & possb0,posst,mapdb,mapib)
c
       do 210 ia=1,mapda(0,5)
       if (mapda(ia,2).eq.0) goto 210
c
       sa(p)=mapda(ia,3)
       sa(q)=mapda(ia,4)
       ib=mapib(sa(1),1,1)
c
c     possitions of A,B
       nhelp1=mapda(ia,1)
       nhelp2=mapdb(ib,1)
c
c     dimp,dimq
       nhelp3=dimm(mapda(0,1),sa(p))
       nhelp4=dimm(mapda(0,2),sa(q))
       call map21 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,p,q,1)
c
 210    continue
c
       else
c     RC=2 - bad typ for: nind=2
       rc=2
       end if
c
       else if (nind.eq.3) then
c
c     *********** 3 index ***********
c
       typ=mapda(0,6)
c
       if (typ.eq.0) then
c
c3.1  map A(p,q,r) -> B(p1,q1,r1)
c
c     get mapdb,mapib
c
       type(p)=mapda(0,1)
       type(q)=mapda(0,2)
       type(r)=mapda(0,3)
       call grc0 (nind,0,type(1),type(2),type(3),0,ssa,
     & possb0,posst,mapdb,mapib)
c
       do 310 ia=1,mapda(0,5)
       if (mapda(ia,2).eq.0) goto 310
c
       sa(p)=mapda(ia,3)
       sa(q)=mapda(ia,4)
       sa(r)=mapda(ia,5)
       ib=mapib(sa(1),sa(2),1)
c
c     possitions of A,B
       nhelp1=mapda(ia,1)
       nhelp2=mapdb(ib,1)
c
c     dimp,dimq,r
       nhelp3=dimm(mapda(0,1),sa(p))
       nhelp4=dimm(mapda(0,2),sa(q))
       nhelp5=dimm(mapda(0,3),sa(r))
       call map31 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,p,q,r,1)
c
 310    continue
c
       else if (typ.eq.1) then
c
c3.2  map A(pq,r)
c     => only sophystical order of p,q,r is 2,3,1 i.e. B(r,pq)
c     => newtyp=2
c
       if ((p.eq.2).and.(q.eq.3).and.(r.eq.1)) then
c
c     get mapdb,mapib
c
       type(p)=mapda(0,1)
       type(q)=mapda(0,2)
       type(r)=mapda(0,3)
       call grc0 (nind,2,type(1),type(2),type(3),0,ssa,
     & possb0,posst,mapdb,mapib)
c
       do 320 ia=1,mapda(0,5)
       if (mapda(ia,2).eq.0) goto 320
c
       sa(p)=mapda(ia,3)
       sa(q)=mapda(ia,4)
       sa(r)=mapda(ia,5)
       ib=mapib(sa(1),sa(2),1)
c
c     possitions of A,B
       nhelp1=mapda(ia,1)
       nhelp2=mapdb(ib,1)
c
c     dimp,dimq,dimr
       nhelp3=dimm(mapda(0,1),sa(p))
       nhelp4=dimm(mapda(0,2),sa(q))
       nhelp5=dimm(mapda(0,3),sa(r))
c
c     dimpq
       if (sa(1).eq.sa(2)) then
       nhelp7=nhelp3*(nhelp3-1)/2
       else
       nhelp7=nhelp3*nhelp4
       end if
       call map21 (wrk(nhelp1),wrk(nhelp2),nhelp7,nhelp5,2,1,1)
c
 320    continue
c
       else
c     RC=3 - bad order of p,q,r for: nind=3,typ=1
       rc=3
       end if
c
       else if (typ.eq.2) then
c
c3.3  map A(p,qr)
c     => only sophystical order of p,q,r is 3,1,2 i.e. B(qr,p)
c     => typb=1
c
       if ((p.eq.3).and.(q.eq.1).and.(r.eq.2)) then
c
c     get mapdb,mapib
c
       type(p)=mapda(0,1)
       type(q)=mapda(0,2)
       type(r)=mapda(0,3)
       call grc0 (nind,1,type(1),type(2),type(3),0,ssa,
     & possb0,posst,mapdb,mapib)
c
       do 330 ia=1,mapda(0,5)
       if (mapda(ia,2).eq.0) goto 330
c
       sa(p)=mapda(ia,3)
       sa(q)=mapda(ia,4)
       sa(r)=mapda(ia,5)
       ib=mapib(sa(1),sa(2),1)
c
c     possitions of A,B
       nhelp1=mapda(ia,1)
       nhelp2=mapdb(ib,1)
c
c     dimp,dimq,dimr
       nhelp3=dimm(mapda(0,1),sa(p))
       nhelp4=dimm(mapda(0,2),sa(q))
       nhelp5=dimm(mapda(0,3),sa(r))
c
c     dimpq
       if (sa(1).eq.sa(2)) then
       nhelp7=nhelp4*(nhelp4-1)/2
       else
       nhelp7=nhelp4*nhelp5
       end if
       call map21 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp7,2,1,1)
c
 330    continue
c
       else
c     RC=4 - bad order of p,q,r for: nind=3,typ=2
       rc=4
       end if
c
       else
c     RC=5 - bad type for: nind=3
       rc=5
       end if
c
       else if (nind.eq.4) then
c
c     *********** 4 index ***********
c
       typ=mapda(0,6)
c
       if (typ.eq.0) then
c
c4.1  map A(p,q,r,s)
c
c     get mapdb,mapib
c
       type(p)=mapda(0,1)
       type(q)=mapda(0,2)
       type(r)=mapda(0,3)
       type(s)=mapda(0,4)
       call grc0 (nind,0,type(1),type(2),type(3),type(4),ssa,
     & possb0,posst,mapdb,mapib)
c
       do 410 ia=1,mapda(0,5)
       if (mapda(ia,2).eq.0) goto 410
c
       sa(p)=mapda(ia,3)
       sa(q)=mapda(ia,4)
       sa(r)=mapda(ia,5)
       sa(s)=mapda(ia,6)
       ib=mapib(sa(1),sa(2),sa(3))
c
c     possitions of A,B
       nhelp1=mapda(ia,1)
       nhelp2=mapdb(ib,1)
c
c     dimp,dimq,dimr,dims
       nhelp3=dimm(mapda(0,1),sa(p))
       nhelp4=dimm(mapda(0,2),sa(q))
       nhelp5=dimm(mapda(0,3),sa(r))
       nhelp6=dimm(mapda(0,4),sa(s))
       call map41 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,nhelp6,
     &             p,q,r,s,1)
c
 410    continue
c
       else if ((typ.gt.0).and.(typ.lt.4)) then
c
c4.2  map A(p,q,r,s) with one limitation
c
c     def newtype and do corresponding tests
c
       if (typ.eq.1) then
       if ((p.gt.3).or.((q-p).ne.1)) then
c     RC=6 bad order of p,q,r,s, for: nind=4,typ=1
       rc=6
       return
       else
       nhelp1=p
       end if
       else if (typ.eq.2) then
       if ((q.gt.3).or.((r-q).ne.1)) then
c     RC=7 bad order of p,q,r,s, for: nind=4,typ=2
       rc=7
       return
       else
       nhelp1=q
       end if
       else if (typ.eq.3) then
       if ((r.gt.3).or.((s-r).ne.1)) then
c     RC=8 bad order of p,q,r,s, for: nind=4,typ=3
       rc=8
       return
       else
       nhelp1=r
       end if
       end if
       newtyp=nhelp1
c
c     get mapdb,mapib
c
       type(p)=mapda(0,1)
       type(q)=mapda(0,2)
       type(r)=mapda(0,3)
       type(s)=mapda(0,4)
       call grc0 (nind,newtyp,type(1),type(2),type(3),type(4),ssa,
     & possb0,posst,mapdb,mapib)

       do 420 ia=1,mapda(0,5)
       if (mapda(ia,2).eq.0) goto 420
c
       sa(p)=mapda(ia,3)
       sa(q)=mapda(ia,4)
       sa(r)=mapda(ia,5)
       sa(s)=mapda(ia,6)
       ib=mapib(sa(1),sa(2),sa(3))
c
c     possitions of A,B
       nhelp1=mapda(ia,1)
       nhelp2=mapdb(ib,1)
c
c     dimp,dimq,dimr,dims
       dl(p)=dimm(mapda(0,1),sa(p))
       dl(q)=dimm(mapda(0,2),sa(q))
       dl(r)=dimm(mapda(0,3),sa(r))
       dl(s)=dimm(mapda(0,4),sa(s))
c
       if ((newtyp.eq.1).and.(sa(1).eq.sa(2))) then
c     B(p1q1,r1,s1) case => oldtyp can be 2,3
       nhelp7=dl(1)*(dl(1)-1)/2
c
       if (typ.eq.1) then
c     A(pq,r,s) case
c     A(pq,r,s) -> A(pq,s,r)
       call map31 (wrk(nhelp1),wrk(nhelp2),nhelp7,dl(r),dl(s),1,3,2,1)
       else if (typ.eq.2) then
c     A(p,qr,s) case
       if (p.eq.3) then
c     A (p,qr,s) -> B(qr,p,s)
       call map31 (wrk(nhelp1),wrk(nhelp2),dl(p),nhelp7,dl(s),2,1,3,1)
       else
c     A (p,qr,s) -> B(qr,s,p)
       call map31 (wrk(nhelp1),wrk(nhelp2),dl(p),nhelp7,dl(s),3,1,2,1)
       end if
       else
c     A(p,q,rs) case
       if (p.eq.3) then
c     A (p,q,rs) -> B(rs,p,q)
       call map31 (wrk(nhelp1),wrk(nhelp2),dl(p),dl(q),nhelp7,2,3,1,1)
       else
c     A (p,q,rs) -> B(rs,q,p)
       call map31 (wrk(nhelp1),wrk(nhelp2),dl(p),dl(q),nhelp7,3,2,1,1)
       end if
       end if
c
       else if ((newtyp.eq.2).and.(sa(2).eq.sa(3))) then
c     B(p1,q1r1,s1) case => oldtyp can be 1,3
       nhelp7=dl(2)*(dl(2)-1)/2
       if (typ.eq.1) then
c     A(pq,r,s) case
       if (r.eq.1) then
c     A (pq,r,s) -> B(r,pq,s)
       call map31 (wrk(nhelp1),wrk(nhelp2),nhelp7,dl(r),dl(s),2,1,3,1)
       else
c     A (pq,r,s) -> B(s,pq,p)
       call map31 (wrk(nhelp1),wrk(nhelp2),nhelp7,dl(r),dl(s),2,3,1,1)
       end if
       else if (typ.eq.2) then
c     A(p,qr,s) case
c     A(p,qr,s) -> B(s,qr,p)
       call map31 (wrk(nhelp1),wrk(nhelp2),dl(p),nhelp7,dl(s),3,2,1,1)
       else
c     A(p,q,rs) case
       if (p.eq.1) then
c     A (p,q,rs) -> B(p,rs,q)
       call map31 (wrk(nhelp1),wrk(nhelp2),dl(p),dl(q),nhelp7,1,3,2,1)
       else
c     A (p,q,rs) -> B(q,rs,p)
       call map31 (wrk(nhelp1),wrk(nhelp2),dl(p),dl(q),nhelp7,3,1,2,1)
       end if
       end if
c
       else if ((newtyp.eq.3).and.(sa(3).eq.sa(4))) then
c     B(p1,q1,r1s1) case => oldtyp can be 1,2
       nhelp7=dl(3)*(dl(3)-1)/2
       if (typ.eq.1) then
c     A(pq,r,s) case
       if (r.eq.2) then
c     A (pq,r,s) -> B(r,s,pq)
       call map31 (wrk(nhelp1),wrk(nhelp2),nhelp7,dl(r),dl(s),3,1,2,1)
       else
c     A (pq,r,s) -> B(s,r,pq)
       call map31 (wrk(nhelp1),wrk(nhelp2),nhelp7,dl(r),dl(s),3,2,1,1)
       end if
       else if (typ.eq.2) then
c     A(p,qr,s) case
       if (p.eq.1) then
c     A (p,qr,s) -> B(p,s,qr)
       call map31 (wrk(nhelp1),wrk(nhelp2),dl(p),nhelp7,dl(s),1,3,2,1)
       else
c     A (p,qr,s) -> B(s,p,qr)
       call map31 (wrk(nhelp1),wrk(nhelp2),dl(p),nhelp7,dl(s),2,3,1,1)
       end if
       else
c     A(p,q,rs) case
c     A(p,q,rs) -> A(q,p,rs)
       call map31 (wrk(nhelp1),wrk(nhelp2),dl(p),dl(q),nhelp7,2,1,3,1)
       end if
c
       else
c     B(p1,q1,r1,s1) case
       call map41 (wrk(nhelp1),wrk(nhelp2),dl(p),dl(q),dl(r),dl(s),p,q,
     &             r,s,1)
       end if
c
 420    continue
c
c
       else if (typ.eq.4) then
c
c4.2  map A(pq,rs) -> B(rs,pq)
c
c     tests
c
       if ((p.ne.3).and.(q.ne.4).and.(r.ne.1).and.(s.ne.2)) then
c     RC=9 : nind=4, typA=4 (bad order of p,q,r,s, Stup)
       rc=9
       return
       end if
c
c     get mapdb,mapib
c
       type(p)=mapda(0,1)
       type(q)=mapda(0,2)
       type(r)=mapda(0,3)
       type(s)=mapda(0,4)
       call grc0 (nind,4,type(1),type(2),type(3),type(4),ssa,
     & possb0,posst,mapdb,mapib)
c
       do 430 ia=1,mapda(0,5)
       if (mapda(ia,2).eq.0) goto 430
c
       sa(p)=mapda(ia,3)
       sa(q)=mapda(ia,4)
       sa(r)=mapda(ia,5)
       sa(s)=mapda(ia,6)
       ib=mapib(sa(1),sa(2),sa(3))
c
c     possitions of A,B
       nhelp1=mapda(ia,1)
       nhelp2=mapdb(ib,1)
c
c     dimp,dimq,dimr,dims
       nhelp3=dimm(mapda(0,1),sa(p))
       nhelp4=dimm(mapda(0,2),sa(q))
       nhelp5=dimm(mapda(0,3),sa(r))
       nhelp6=dimm(mapda(0,4),sa(s))
c
       if ((sa(p).eq.sa(q)).and.(sa(r).eq.sa(s))) then
c     A(pq,rs) -> B(rs,pq)
       nhelp3=nhelp3*(nhelp3-1)/2
       nhelp5=nhelp5*(nhelp5-1)/2
       call map21 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp5,2,1,1)
c
       else if (sa(p).eq.sa(q)) then
c     A(pq,r,s) -> B(r,s,pq)
       nhelp3=nhelp3*(nhelp3-1)/2
       call map31 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp5,nhelp6,3,1,2,1)
c
       else if (sa(r).eq.sa(s)) then
c     A(p,q,rs) -> B(rs,p,q)
       nhelp5=nhelp5*(nhelp5-1)/2
       call map31 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,2,3,1,1)
c
       else
c     A(p,q,r,s) -> B(r,s,p,q)
       call map41 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,nhelp6,
     &             3,4,1,2,1)
       end if
c
 430    continue
c
       else
c     RC=10 bad type for: nind=4
       rc=10
       end if
c
       else
c
c     *********  more than 4 indexes
c
c     RC=11 bad number of indexes
       rc=11
       end if
c
c
       return
       end
c
c     -----------
c
       subroutine noperm (wrk,wrksize,
     & mapda,mapia,mapdb,mapib,poss0,posst)
c
c     realize mapping without permutation
c     define mapd,mapi
c
#include "ccsd1.fh"
#include "wrk.fh"
c
       integer poss0,posst
       integer mapda(0:512,1:6),mapdb(0:512,1:6)
       integer mapia(1:8,1:8,1:8),mapib(1:8,1:8,1:8)
c
c     help variables
c
       integer ib,nhelp,i,j,k
c
c     def mapib
c
       do 10 k=1,nsym
       do 11 j=1,nsym
       do 12 i=1,nsym
       mapib(i,j,k)=mapia(i,j,k)
 12     continue
 11     continue
 10     continue
c
c     def initial values
c
       do nhelp=1,6
       mapdb(0,nhelp)=mapda(0,nhelp)
       end do
c
       posst=poss0
       do 100 ib=1,mapda(0,5)
       do nhelp=2,6
       mapdb(ib,nhelp)=mapda(ib,nhelp)
       end do
       mapdb(ib,1)=posst
       posst=posst+mapdb(ib,2)
c
       call map11 (wrk(mapda(ib,1)),wrk(mapdb(ib,1)),mapda(ib,2),1)
c
 100    continue
c
       return
       end
c
c     --------------------------
c
       subroutine map41 (a,b,dimp,dimq,dimr,dims,p,q,r,s,nfact)
c
c     maping A(p1,q1,r1,s1) -> nfact* B(p2,q2,r2,s2)
c
       real*8 a(*)
       real*8 b(*)
       integer dim(4)
       integer dimp,dimq,dimr,dims,p,q,r,s,nfact
       dim(p)=dimp
       dim(q)=dimq
       dim(r)=dimr
       dim(s)=dims
       call map42 (a,b,dimp,dimq,dimr,dims,dim(1),dim(2),dim(3),dim(4),
     &             p,q,r,s,nfact)
       return
       end
c
       subroutine map42 (a,b,dimp,dimq,dimr,dims,dim1,dim2,dim3,dim4,p,
     &                   q,r,s,nfact)
c
       integer dimp,dimq,dimr,dims,dim1,dim2,dim3,dim4,p,q,r,s,nfact
       real*8 a(1:dimp,1:dimq,1:dimr,1:dims)
       real*8 b(1:dim1,1:dim2,1:dim3,1:dim4)
c     integer index(1:4)
c
       integer pp,qq,rr,ss
c
       if (nfact.eq.1) then
c
c     factor + 1
c
c     do 100 ss=1,dims
c     index(s)=ss
c     do 100 rr=1,dimr
c     index(r)=rr
c     do 100 qq=1,dimq
c     index(q)=qq
c     do 100 pp=1,dimp
c     index(p)=pp
c     b(index(1),index(2),index(3),index(4))=a(pp,qq,rr,ss)
c100  continue
c
       if (p.eq.1) then
c     1***
       if (q.eq.2) then
c     12**
       if (r.eq.3) then
c     123* (4)
       do 111 ss=1,dims
       do 1110 rr=1,dimr
       do 1111 qq=1,dimq
       do 1112 pp=1,dimp
       b(pp,qq,rr,ss)=a(pp,qq,rr,ss)
 1112   continue
 1111   continue
 1110   continue
 111    continue
       else
c     124* (3)
       do 112 ss=1,dims
       do 1120 rr=1,dimr
       do 1121 qq=1,dimq
       do 1122 pp=1,dimp
       b(pp,qq,ss,rr)=a(pp,qq,rr,ss)
 1122   continue
 1121   continue
 1120   continue
 112    continue
       end if
       else if (q.eq.3) then
c     13**
       if (r.eq.2) then
c     132* (4)
       do 113 ss=1,dims
       do 1130 rr=1,dimr
       do 1131 qq=1,dimq
       do 1132 pp=1,dimp
       b(pp,rr,qq,ss)=a(pp,qq,rr,ss)
 1132   continue
 1131   continue
 1130   continue
 113    continue
       else
c     134* (2)
       do 114 ss=1,dims
       do 1140 rr=1,dimr
       do 1141 qq=1,dimq
       do 1142 pp=1,dimp
       b(pp,ss,qq,rr)=a(pp,qq,rr,ss)
 1142   continue
 1141   continue
 1140   continue
 114    continue
       end if
       else if (q.eq.4) then
c     14**
       if (r.eq.2) then
c     142* (3)
       do 115 ss=1,dims
       do 1150 rr=1,dimr
       do 1151 qq=1,dimq
       do 1152 pp=1,dimp
       b(pp,rr,ss,qq)=a(pp,qq,rr,ss)
 1152   continue
 1151   continue
 1150   continue
 115    continue
       else
c     143* (2)
       do 116 ss=1,dims
       do 1160 rr=1,dimr
       do 1161 qq=1,dimq
       do 1162 pp=1,dimp
       b(pp,ss,rr,qq)=a(pp,qq,rr,ss)
 1162   continue
 1161   continue
 1160   continue
 116    continue
       end if
       end if
c
       else if (p.eq.2) then
c     2***
       if (q.eq.1) then
c     21**
       if (r.eq.3) then
c     213* (4)
       do 121 ss=1,dims
       do 1210 rr=1,dimr
       do 1211 qq=1,dimq
       do 1212 pp=1,dimp
       b(qq,pp,rr,ss)=a(pp,qq,rr,ss)
 1212   continue
 1211   continue
 1210   continue
 121    continue
       else
c     214* (3)
       do 122 ss=1,dims
       do 1220 rr=1,dimr
       do 1221 qq=1,dimq
       do 1222 pp=1,dimp
       b(qq,pp,ss,rr)=a(pp,qq,rr,ss)
 1222   continue
 1221   continue
 1220   continue
 122    continue
       end if
       else if (q.eq.3) then
c     23**
       if (r.eq.1) then
c     231* (4)
       do 123 ss=1,dims
       do 1230 rr=1,dimr
       do 1231 qq=1,dimq
       do 1232 pp=1,dimp
       b(rr,pp,qq,ss)=a(pp,qq,rr,ss)
 1232   continue
 1231   continue
 1230   continue
 123    continue
       else
c     234* (1)
       do 124 ss=1,dims
       do 1240 rr=1,dimr
       do 1241 qq=1,dimq
       do 1242 pp=1,dimp
       b(ss,pp,qq,rr)=a(pp,qq,rr,ss)
 1242   continue
 1241   continue
 1240   continue
 124    continue
       end if
       else if (q.eq.4) then
c     24**
       if (r.eq.3) then
c     243* (1)
       do 125 ss=1,dims
       do 1250 rr=1,dimr
       do 1251 qq=1,dimq
       do 1252 pp=1,dimp
       b(ss,pp,rr,qq)=a(pp,qq,rr,ss)
 1252   continue
 1251   continue
 1250   continue
 125    continue
       else
c     241* (3)
       do 126 ss=1,dims
       do 1260 rr=1,dimr
       do 1261 qq=1,dimq
       do 1262 pp=1,dimp
       b(rr,pp,ss,qq)=a(pp,qq,rr,ss)
 1262   continue
 1261   continue
 1260   continue
 126    continue
       end if
       end if
c
       else if (p.eq.3) then
c     3***
       if (q.eq.1) then
c     31**
       if (r.eq.2) then
c     312* (4)
       do 131 ss=1,dims
       do 1310 rr=1,dimr
       do 1311 qq=1,dimq
       do 1312 pp=1,dimp
       b(qq,rr,pp,ss)=a(pp,qq,rr,ss)
 1312   continue
 1311   continue
 1310   continue
 131    continue
       else
c     314* (2)
       do 132 ss=1,dims
       do 1320 rr=1,dimr
       do 1321 qq=1,dimq
       do 1322 pp=1,dimp
       b(qq,ss,pp,rr)=a(pp,qq,rr,ss)
 1322   continue
 1321   continue
 1320   continue
 132    continue
       end if
       else if (q.eq.2) then
c     32**
       if (r.eq.1) then
c     321* (4)
       do 133 ss=1,dims
       do 1330 rr=1,dimr
       do 1331 qq=1,dimq
       do 1332 pp=1,dimp
       b(rr,qq,pp,ss)=a(pp,qq,rr,ss)
 1332   continue
 1331   continue
 1330   continue
 133    continue
       else
c     324* (1)
       do 134 ss=1,dims
       do 1340 rr=1,dimr
       do 1341 qq=1,dimq
       do 1342 pp=1,dimp
       b(ss,qq,pp,rr)=a(pp,qq,rr,ss)
 1342   continue
 1341   continue
 1340   continue
 134    continue
       end if
       else if (q.eq.4) then
c     34**
       if (r.eq.1) then
c     341* (2)
       do 135 ss=1,dims
       do 1350 rr=1,dimr
       do 1351 qq=1,dimq
       do 1352 pp=1,dimp
       b(rr,ss,pp,qq)=a(pp,qq,rr,ss)
 1352   continue
 1351   continue
 1350   continue
 135    continue
       else
c     342* (1)
       do 136 ss=1,dims
       do 1360 rr=1,dimr
       do 1361 qq=1,dimq
       do 1362 pp=1,dimp
       b(ss,rr,pp,qq)=a(pp,qq,rr,ss)
 1362   continue
 1361   continue
 1360   continue
 136    continue
       end if
       end if
c
       else if (p.eq.4) then
c     4***
       if (q.eq.1) then
c     41**
       if (r.eq.3) then
c     413* (2)
       do 141 ss=1,dims
       do 1410 rr=1,dimr
       do 1411 qq=1,dimq
       do 1412 pp=1,dimp
       b(qq,ss,rr,pp)=a(pp,qq,rr,ss)
 1412   continue
 1411   continue
 1410   continue
 141    continue
       else
c     412* (3)
       do 142 ss=1,dims
       do 1420 rr=1,dimr
       do 1421 qq=1,dimq
       do 1422 pp=1,dimp
       b(qq,rr,ss,pp)=a(pp,qq,rr,ss)
 1422   continue
 1421   continue
 1420   continue
 142    continue
       end if
       else if (q.eq.2) then
c     42**
       if (r.eq.1) then
c     421* (3)
       do 143 ss=1,dims
       do 1430 rr=1,dimr
       do 1431 qq=1,dimq
       do 1432 pp=1,dimp
       b(rr,qq,ss,pp)=a(pp,qq,rr,ss)
 1432   continue
 1431   continue
 1430   continue
 143    continue
       else
c     423* (1)
       do 144 ss=1,dims
       do 1440 rr=1,dimr
       do 1441 qq=1,dimq
       do 1442 pp=1,dimp
       b(ss,qq,rr,pp)=a(pp,qq,rr,ss)
 1442   continue
 1441   continue
 1440   continue
 144    continue
       end if
       else if (q.eq.3) then
c     43**
       if (r.eq.1) then
c     431* (2)
       do 145 ss=1,dims
       do 1450 rr=1,dimr
       do 1451 qq=1,dimq
       do 1452 pp=1,dimp
       b(rr,ss,qq,pp)=a(pp,qq,rr,ss)
 1452   continue
 1451   continue
 1450   continue
 145    continue
       else
c     432* (1)
       do 146 ss=1,dims
       do 1460 rr=1,dimr
       do 1461 qq=1,dimq
       do 1462 pp=1,dimp
       b(ss,rr,qq,pp)=a(pp,qq,rr,ss)
 1462   continue
 1461   continue
 1460   continue
 146    continue
       end if
       end if
c
       end if
c
       else
c
c     factor = -1
c
c     do 200 ss=1,dims
c     index(s)=ss
c     do 200 rr=1,dimr
c     index(r)=rr
c     do 200 qq=1,dimq
c     index(q)=qq
c     do 200 pp=1,dimp
c     index(p)=pp
c     b(index(1),index(2),index(3),index(4))=-a(pp,qq,rr,ss)
c200  continue
c
       if (p.eq.1) then
c     1***
       if (q.eq.2) then
c     12**
       if (r.eq.3) then
c     123* (4)
       do 211 ss=1,dims
       do 2110 rr=1,dimr
       do 2111 qq=1,dimq
       do 2112 pp=1,dimp
       b(pp,qq,rr,ss)=-a(pp,qq,rr,ss)
 2112   continue
 2111   continue
 2110   continue
 211    continue
       else
c     124* (3)
       do 212 ss=1,dims
       do 2120 rr=1,dimr
       do 2121 qq=1,dimq
       do 2122 pp=1,dimp
       b(pp,qq,ss,rr)=-a(pp,qq,rr,ss)
 2122   continue
 2121   continue
 2120   continue
 212    continue
       end if
       else if (q.eq.3) then
c     13**
       if (r.eq.2) then
c     132* (4)
       do 213 ss=1,dims
       do 2130 rr=1,dimr
       do 2131 qq=1,dimq
       do 2132 pp=1,dimp
       b(pp,rr,qq,ss)=-a(pp,qq,rr,ss)
 2132   continue
 2131   continue
 2130   continue
 213    continue
       else
c     134* (2)
       do 214 ss=1,dims
       do 2140 rr=1,dimr
       do 2141 qq=1,dimq
       do 2142 pp=1,dimp
       b(pp,ss,qq,rr)=-a(pp,qq,rr,ss)
 2142   continue
 2141   continue
 2140   continue
 214    continue
       end if
       else if (q.eq.4) then
c     14**
       if (r.eq.2) then
c     142* (3)
       do 215 ss=1,dims
       do 2150 rr=1,dimr
       do 2151 qq=1,dimq
       do 2152 pp=1,dimp
       b(pp,rr,ss,qq)=-a(pp,qq,rr,ss)
 2152   continue
 2151   continue
 2150   continue
 215    continue
       else
c     143* (2)
       do 216 ss=1,dims
       do 2160 rr=1,dimr
       do 2161 qq=1,dimq
       do 2162 pp=1,dimp
       b(pp,ss,rr,qq)=-a(pp,qq,rr,ss)
 2162   continue
 2161   continue
 2160   continue
 216    continue
       end if
       end if
c
       else if (p.eq.2) then
c     2***
       if (q.eq.1) then
c     21**
       if (r.eq.3) then
c     213* (4)
       do 221 ss=1,dims
       do 2210 rr=1,dimr
       do 2211 qq=1,dimq
       do 2212 pp=1,dimp
       b(qq,pp,rr,ss)=-a(pp,qq,rr,ss)
 2212   continue
 2211   continue
 2210   continue
 221    continue
       else
c     214* (3)
       do 222 ss=1,dims
       do 2220 rr=1,dimr
       do 2221 qq=1,dimq
       do 2222 pp=1,dimp
       b(qq,pp,ss,rr)=-a(pp,qq,rr,ss)
 2222   continue
 2221   continue
 2220   continue
 222    continue
       end if
       else if (q.eq.3) then
c     23**
       if (r.eq.1) then
c     231* (4)
       do 223 ss=1,dims
       do 2230 rr=1,dimr
       do 2231 qq=1,dimq
       do 2232 pp=1,dimp
       b(rr,pp,qq,ss)=-a(pp,qq,rr,ss)
 2232   continue
 2231   continue
 2230   continue
 223    continue
       else
c     234* (1)
       do 224 ss=1,dims
       do 2240 rr=1,dimr
       do 2241 qq=1,dimq
       do 2242 pp=1,dimp
       b(ss,pp,qq,rr)=-a(pp,qq,rr,ss)
 2242   continue
 2241   continue
 2240   continue
 224    continue
       end if
       else if (q.eq.4) then
c     24**
       if (r.eq.3) then
c     243* (1)
       do 225 ss=1,dims
       do 2250 rr=1,dimr
       do 2251 qq=1,dimq
       do 2252 pp=1,dimp
       b(ss,pp,rr,qq)=-a(pp,qq,rr,ss)
 2252   continue
 2251   continue
 2250   continue
 225    continue
       else
c     241* (3)
       do 226 ss=1,dims
       do 2260 rr=1,dimr
       do 2261 qq=1,dimq
       do 2262 pp=1,dimp
       b(rr,pp,ss,qq)=-a(pp,qq,rr,ss)
 2262   continue
 2261   continue
 2260   continue
 226    continue
       end if
       end if
c
       else if (p.eq.3) then
c     3***
       if (q.eq.1) then
c     31**
       if (r.eq.2) then
c     312* (4)
       do 231 ss=1,dims
       do 2310 rr=1,dimr
       do 2311 qq=1,dimq
       do 2312 pp=1,dimp
       b(qq,rr,pp,ss)=-a(pp,qq,rr,ss)
 2312   continue
 2311   continue
 2310   continue
 231    continue
       else
c     314* (2)
       do 232 ss=1,dims
       do 2320 rr=1,dimr
       do 2321 qq=1,dimq
       do 2322 pp=1,dimp
       b(qq,ss,pp,rr)=-a(pp,qq,rr,ss)
 2322   continue
 2321   continue
 2320   continue
 232    continue
       end if
       else if (q.eq.2) then
c     32**
       if (r.eq.1) then
c     321* (4)
       do 233 ss=1,dims
       do 2330 rr=1,dimr
       do 2331 qq=1,dimq
       do 2332 pp=1,dimp
       b(rr,qq,pp,ss)=-a(pp,qq,rr,ss)
 2332   continue
 2331   continue
 2330   continue
 233    continue
       else
c     324* (1)
       do 234 ss=1,dims
       do 2340 rr=1,dimr
       do 2341 qq=1,dimq
       do 2342 pp=1,dimp
       b(ss,qq,pp,rr)=-a(pp,qq,rr,ss)
 2342   continue
 2341   continue
 2340   continue
 234    continue
       end if
       else if (q.eq.4) then
c     34**
       if (r.eq.1) then
c     341* (2)
       do 235 ss=1,dims
       do 2350 rr=1,dimr
       do 2351 qq=1,dimq
       do 2352 pp=1,dimp
       b(rr,ss,pp,qq)=-a(pp,qq,rr,ss)
 2352   continue
 2351   continue
 2350   continue
 235    continue
       else
c     342* (1)
       do 236 ss=1,dims
       do 2360 rr=1,dimr
       do 2361 qq=1,dimq
       do 2362 pp=1,dimp
       b(ss,rr,pp,qq)=-a(pp,qq,rr,ss)
 2362   continue
 2361   continue
 2360   continue
 236    continue
       end if
       end if
c
       else if (p.eq.4) then
c     4***
       if (q.eq.1) then
c     41**
       if (r.eq.3) then
c     413* (2)
       do 241 ss=1,dims
       do 2410 rr=1,dimr
       do 2411 qq=1,dimq
       do 2412 pp=1,dimp
       b(qq,ss,rr,pp)=-a(pp,qq,rr,ss)
 2412   continue
 2411   continue
 2410   continue
 241    continue
       else
c     412* (3)
       do 242 ss=1,dims
       do 2420 rr=1,dimr
       do 2421 qq=1,dimq
       do 2422 pp=1,dimp
       b(qq,rr,ss,pp)=-a(pp,qq,rr,ss)
 2422   continue
 2421   continue
 2420   continue
 242    continue
       end if
       else if (q.eq.2) then
c     42**
       if (r.eq.1) then
c     421* (3)
       do 243 ss=1,dims
       do 2430 rr=1,dimr
       do 2431 qq=1,dimq
       do 2432 pp=1,dimp
       b(rr,qq,ss,pp)=-a(pp,qq,rr,ss)
 2432   continue
 2431   continue
 2430   continue
 243    continue
       else
c     423* (1)
       do 244 ss=1,dims
       do 2440 rr=1,dimr
       do 2441 qq=1,dimq
       do 2442 pp=1,dimp
       b(ss,qq,rr,pp)=-a(pp,qq,rr,ss)
 2442   continue
 2441   continue
 2440   continue
 244    continue
       end if
       else if (q.eq.3) then
c     43**
       if (r.eq.1) then
c     431* (2)
       do 245 ss=1,dims
       do 2450 rr=1,dimr
       do 2451 qq=1,dimq
       do 2452 pp=1,dimp
       b(rr,ss,qq,pp)=-a(pp,qq,rr,ss)
 2452   continue
 2451   continue
 2450   continue
 245    continue
       else
c     432* (1)
       do 246 ss=1,dims
       do 2460 rr=1,dimr
       do 2461 qq=1,dimq
       do 2462 pp=1,dimp
       b(ss,rr,qq,pp)=-a(pp,qq,rr,ss)
 2462   continue
 2461   continue
 2460   continue
 246    continue
       end if
       end if
c
       end if
c
       end if
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(s)
       end
c
c     --------------------------
c
       subroutine map31 (a,b,dimp,dimq,dimr,p,q,r,nfact)
c
c     maping A(p1,q1,r1) -> nfact*B(p2,q2,r2)
c
       real*8 a(*)
       real*8 b(*)
       integer dim(3)
       integer dimp,dimq,dimr,p,q,r,nfact
       dim(p)=dimp
       dim(q)=dimq
       dim(r)=dimr
       call map32 (a,b,dimp,dimq,dimr,dim(1),dim(2),dim(3),p,q,r,nfact)
       return
       end
c
       subroutine map32 (a,b,dimp,dimq,dimr,dim1,dim2,dim3,p,q,r,nfact)
c
       integer dimp,dimq,dimr,dim1,dim2,dim3,p,q,r,nfact
       real*8 a(1:dimp,1:dimq,1:dimr)
       real*8 b(1:dim1,1:dim2,1:dim3)
c     integer index(1:3)
c
       integer pp,qq,rr
c
       if (nfact.eq.1) then
c
c     fact = + 1
c
c     do 100 rr=1,dimr
c     index(r)=rr
c     do 100 qq=1,dimq
c     index(q)=qq
c     do 100 pp=1,dimp
c     index(p)=pp
c     b(index(1),index(2),index(3))=a(pp,qq,rr)
c100  continue
c
       if (p.eq.1) then
c     1**
       if (q.eq.2) then
c     12* (3)
       do 110 rr=1,dimr
       do 111 qq=1,dimq
       do 112 pp=1,dimp
       b(pp,qq,rr)=a(pp,qq,rr)
 112    continue
 111    continue
 110    continue
       else
c     13* (2)
       do 120 rr=1,dimr
       do 121 qq=1,dimq
       do 122 pp=1,dimp
       b(pp,rr,qq)=a(pp,qq,rr)
 122    continue
 121    continue
 120    continue
       end if
       else if (p.eq.2) then
c     2**
       if (q.eq.1) then
c     21* (3)
       do 130 rr=1,dimr
       do 131 qq=1,dimq
       do 132 pp=1,dimp
       b(qq,pp,rr)=a(pp,qq,rr)
 132    continue
 131    continue
 130    continue
       else
c     23* (1)
       do 140 rr=1,dimr
       do 141 qq=1,dimq
       do 142 pp=1,dimp
       b(rr,pp,qq)=a(pp,qq,rr)
 142    continue
 141    continue
 140    continue
       end if
       else if (p.eq.3) then
c     3**
       if (q.eq.1) then
c     31* (2)
       do 150 rr=1,dimr
       do 151 qq=1,dimq
       do 152 pp=1,dimp
       b(qq,rr,pp)=a(pp,qq,rr)
 152    continue
 151    continue
 150    continue
       else
c     32* (1)
       do 160 rr=1,dimr
       do 161 qq=1,dimq
       do 162 pp=1,dimp
       b(rr,qq,pp)=a(pp,qq,rr)
 162    continue
 161    continue
 160    continue
       end if
       end if
c
c
       else
c
c     factor = -1
c
c     do 200 rr=1,dimr
c     index(r)=rr
c     do 200 qq=1,dimq
c     index(q)=qq
c     do 200 pp=1,dimp
c     index(p)=pp
c     b(index(1),index(2),index(3))=-a(pp,qq,rr)
c     200        continue
c
       if (p.eq.1) then
c     1**
       if (q.eq.2) then
c     12* (3)
       do 210 rr=1,dimr
       do 211 qq=1,dimq
       do 212 pp=1,dimp
       b(pp,qq,rr)=-a(pp,qq,rr)
 212    continue
 211    continue
 210    continue
       else
c     13* (2)
       do 220 rr=1,dimr
       do 221 qq=1,dimq
       do 222 pp=1,dimp
       b(pp,rr,qq)=-a(pp,qq,rr)
 222    continue
 221    continue
 220    continue
       end if
       else if (p.eq.2) then
c     2**
       if (q.eq.1) then
c     21* (3)
       do 230 rr=1,dimr
       do 231 qq=1,dimq
       do 232 pp=1,dimp
       b(qq,pp,rr)=-a(pp,qq,rr)
 232    continue
 231    continue
 230    continue
       else
c     23* (1)
       do 240 rr=1,dimr
       do 241 qq=1,dimq
       do 242 pp=1,dimp
       b(rr,pp,qq)=-a(pp,qq,rr)
 242    continue
 241    continue
 240    continue
       end if
       else if (p.eq.3) then
c     3**
       if (q.eq.1) then
c     31* (2)
       do 250 rr=1,dimr
       do 251 qq=1,dimq
       do 252 pp=1,dimp
       b(qq,rr,pp)=-a(pp,qq,rr)
 252    continue
 251    continue
 250    continue
       else
c     32* (1)
       do 260 rr=1,dimr
       do 261 qq=1,dimq
       do 262 pp=1,dimp
       b(rr,qq,pp)=-a(pp,qq,rr)
 262    continue
 261    continue
 260    continue
       end if
       end if
c
       end if
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(r)
       end
c
c     --------------------------
c
       subroutine map21 (a,b,dimp,dimq,p,q,nfact)
c
c     maping A(p1,q1) -> nfact*B(p2,q2)
c
       real*8 a(*)
       real*8 b(*)
       integer dim(2)
       integer dimp,dimq,p,q,nfact
       dim(p)=dimp
       dim(q)=dimq
       call map22 (a,b,dimp,dimq,dim(1),dim(2),p,q,nfact)
       return
       end
c
       subroutine map22 (a,b,dimp,dimq,dim1,dim2,p,q,nfact)
c
       integer dimp,dimq,dim1,dim2,p,q,nfact
       real*8 a(1:dimp,1:dimq)
       real*8 b(1:dim1,1:dim2)
c       integer index(1:2)
c
       integer pp,qq
c
       if (nfact.eq.1) then
c
c     nfact = + 1
c
c     do 100 qq=1,dimq
c     index(q)=qq
c     do 100 pp=1,dimp
c     index(p)=pp
c     b(index(1),index(2))=a(pp,qq)
c100  continue
c
       if (p.eq.1) then
c     1* (2)
       do 110 qq=1,dimq
       do 111 pp=1,dimp
       b(pp,qq)=a(pp,qq)
 111    continue
 110    continue
       else
c     2* (1)
       do 120 qq=1,dimq
       do 121 pp=1,dimp
       b(qq,pp)=a(pp,qq)
 121    continue
 120    continue
       end if
c
c
       else
c
c     nfact = -1
c
c     do 200 qq=1,dimq
c     index(q)=qq
c     do 200 pp=1,dimp
c     index(p)=pp
c     b(index(1),index(2))=-a(pp,qq)
c200  continue
c
       if (p.eq.1) then
c     1* (2)
       do 210 qq=1,dimq
       do 211 pp=1,dimp
       b(pp,qq)=-a(pp,qq)
 211    continue
 210    continue
       else
c     2* (1)
       do 220 qq=1,dimq
       do 221 pp=1,dimp
       b(qq,pp)=-a(pp,qq)
 221    continue
 220    continue
       end if
c
       end if
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(q)
       end
c
c     --------------------------
c
       subroutine map11 (a,b,dimp,nfact)
c
c     maping A(p) -> nfact*B(p)
c
       integer dimp,nfact
       real*8 a(1:dimp)
       real*8 b(1:dimp)
c
       integer pp
c
       if (nfact.eq.1) then
c
       do 100 pp=1,dimp
       b(pp)=a(pp)
 100    continue
c
       else
c
       do 200 pp=1,dimp
       b(pp)=-a(pp)
 200    continue
c
       end if
c
       return
       end
c

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

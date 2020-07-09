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

c     this file contains following routines:
c
c     t3addpck
c     t3aphlp1
c     t3aphlp2
c     t3aphlp3
c     t3aphlp4
c     t3aphlp5
c     t3aphlp6
c     t3aphlp7
c     t3aphlp8
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       subroutine t3addpck (wrk,wrksize,
     & nind,typap,mapda,mapia,ssa,mapdb,mapib,
     & ns,szkey,rc)
c
c     nind   - # of indexes in matrices A,B (Input)
c     typap  - typ of operation (see Table) (Input)
c     mapda  - direct map matrix corresponding to A (Input)
c     mapia  - inverse map matrix corresponding to A (Input)
c     ssa    - overall symmetry state of matrix A (Input)
c     mapdb  - direct map matrix corresponding to B (Input)
c     mapib  - inverse map matrix corresponding to B (Input)
c     ns     - signum of the operation (+-1) (Input)
c     szkey  - zet zero key (I)
c     = 0 no vanishing
c     = 1 set B=0 at the beggining
c     rc     - return (error) code (Output)
c
c     this routine is an addition to CCSD add/pack routine. It do not
c     realize operations, which is done there (add/pack), but only
c     additional ones, required in T3 code, namely:
c
c     Operation                        Nind TypA TypB typap
c
c     B(pqr)  = B + ns*(A(qr,p)-A(pr,q)+A(pq,r))    3    1    5    1
c     B(pq,r) = B + ns*(A(q,r,p)-A(p,r,q))          3    1    1    2
c     B(p,qr) = B + ns*(     -A(p,r,q)+A(p,q,r))    3    1    2    3
c
c     N.B. typab is redundant, it can be determined from mapd's
c     N.B. mapib is redundant
c
#include "t31.fh"
#include "wrk.fh"
c
       integer nind,typap,ssa,ns,szkey,rc
       integer mapda(0:512,1:6)
       integer mapdb(0:512,1:6)
       integer mapia(1:8,1:8,1:8)
       integer mapib(1:8,1:8,1:8)
c
c     help variables
c
       integer ib,possb
       integer ia1,ia2,ia3,possa1,possa2,possa3
       integer syma,symb,symc,dima,dimb,dimc
       integer nhelp1,nhelp2
c
c
c0    some tests
c
       if (nind.ne.3) then
c     RC=1 : nind is not 3 (NCI)
       rc=1
       return
       end if
c
       if ((typap.eq.1).and.
     & ((mapda(0,6).ne.1).and.(mapdb(0,6).ne.5))) then
c     RC=2 : typap=1, typA=1, and typB is not 5 (Stup)
       rc=2
       return
       end if
c
       if ((typap.eq.2).and.
     & ((mapda(0,6).ne.0).and.(mapdb(0,6).ne.1))) then
c     RC=3 : typap=2, typA=0, and typB is not 1 (Stup)
       rc=3
       return
       end if
c
       if ((typap.eq.3).and.
     & ((mapda(0,6).ne.0).and.(mapdb(0,6).ne.2))) then
c     RC=4 : typap=3, typA=0, and typB is not 2 (Stup)
       rc=4
       return
       end if
c
c
       if (typap.eq.1) then
c1    case B(pqr) = B(pqr) + ns*(A(qr,p)-A(pr,q)+A(pq,r))
c
       do 100 ib=1,mapdb(0,5)
c
c1.*  def possition of b
       possb=mapdb(ib,1)
c
c1.*  def symmetry status
       syma=mapdb(ib,3)
       symb=mapdb(ib,4)
       symc=mapdb(ib,5)
c
c1.*  def dimensions
       dima=dimm(mapdb(0,1),syma)
       dimb=dimm(mapdb(0,2),symb)
       dimc=dimm(mapdb(0,3),symc)
c
c1.*  realize packing
c
       if (syma.eq.symc) then
c1.a  case syma=symb=symc
c
c1.a.*find address for a1,2,3 (the same one)
       ia1=mapia(symb,symc,1)
c
c1.a.*def possition of a1,2,3 (the same one)
       possa1=mapda(ia1,1)
c
c1.a.*def additional dimensions ab,abc
       nhelp1=dima*(dima-1)/2
       nhelp2=dima*(dima-1)*(dima-2)/6
c
c1.a.*do packing
       call t3aphlp4 (wrk(possa1),wrk(possb),dima,
     & nhelp1,nhelp2,ns,szkey)
c
       else if (syma.eq.symb) then
c1.b  case syma=symb.ne.symc
c
c1.b.*find address for a1,2,3
       ia1=mapia(symb,symc,1)
       ia2=mapia(syma,symc,1)
       ia3=mapia(syma,symb,1)
c
c1.b.*def possition of a1,2,3
       possa1=mapda(ia1,1)
       possa2=mapda(ia2,1)
       possa3=mapda(ia3,1)
c
c1.b.*def additional dimensions ab
       nhelp1=dima*(dima-1)/2
c
c1.b.*do packing
       call t3aphlp2 (wrk(possa1),wrk(possa2),wrk(possa3),
     & wrk(possb),
     & dima,dimb,dimc,nhelp1,ns,szkey)
c
       else if (symb.eq.symc) then
c1.c  case syma.ne.symb=symc
c
c1.c.*find address for a1,2,3
       ia1=mapia(symb,symc,1)
       ia2=mapia(syma,symc,1)
       ia3=mapia(syma,symb,1)
c
c1.c.*def possition of a1,2,3
       possa1=mapda(ia1,1)
       possa2=mapda(ia2,1)
       possa3=mapda(ia3,1)
c
c1.c.*def additional dimensions bc
       nhelp1=dimb*(dimb-1)/2
c
c1.c.*do packing
       call t3aphlp3 (wrk(possa1),wrk(possa2),wrk(possa3),
     & wrk(possb),
     & dima,dimb,dimc,nhelp1,ns,szkey)
c
       else
c1.d  case syma.ne.symb.ne.symc
c
c1.d.*find address for a1,2,3
       ia1=mapia(symb,symc,1)
       ia2=mapia(syma,symc,1)
       ia3=mapia(syma,symb,1)
c
c1.d.*def possition of a1,2,3
       possa1=mapda(ia1,1)
       possa2=mapda(ia2,1)
       possa3=mapda(ia3,1)
c
c1.d.*do packing
       call t3aphlp1 (wrk(possa1),wrk(possa2),wrk(possa3),
     & wrk(possb),
     & dima,dimb,dimc,ns,szkey)
c
       end if
c
c
 100    continue

       else if (typap.eq.2) then
c2    case B(pq,r) = B(pq,r) + ns*(A(q,r,p)-A(p,r,q))
c
       do 200 ib=1,mapdb(0,5)
c
c2.*  def possition of b
       possb=mapdb(ib,1)
c
c2.*  def symmetry status
       syma=mapdb(ib,3)
       symb=mapdb(ib,4)
       symc=mapdb(ib,5)
c
c2.*  def dimensions
       dima=dimm(mapdb(0,1),syma)
       dimb=dimm(mapdb(0,2),symb)
       dimc=dimm(mapdb(0,3),symc)
c
c2.*  realize packing
c
       if (syma.eq.symb) then
c2.a  case syma=symb,symc
c
c2.a.*find address for a1,2
       ia1=mapia(symb,symc,1)
       ia2=mapia(syma,symc,1)
c
c2.a.*def possition of a1,2,3
       possa1=mapda(ia1,1)
       possa2=mapda(ia2,1)
c
c2.a.*def additional dimensions ab
       nhelp1=dima*(dima-1)/2
c
c2.a.*do packing
       call t3aphlp6 (wrk(possa1),wrk(possa2),wrk(possb),
     & dima,dimb,dimc,nhelp1,ns,szkey)
c
       else
c2.b  case syma.ne.symb,symc
c
c2.b.*find address for a1,2
       ia1=mapia(symb,symc,1)
       ia2=mapia(syma,symc,1)
c
c2.b.*def possition of a1,2,3
       possa1=mapda(ia1,1)
       possa2=mapda(ia2,1)
c
c2.b.*do packing
       call t3aphlp5 (wrk(possa1),wrk(possa2),wrk(possb),
     & dima,dimb,dimc,ns,szkey)
c
       end if
c
 200    continue
c
c
       else if (typap.eq.3) then
c3    case B(p,qr) = B(p,qr) + ns*(-A(p,r,q)+A(p,q,r))
c
       do 300 ib=1,mapdb(0,5)
c
c3.*  def possition of b
       possb=mapdb(ib,1)
c
c3.*  def symmetry status
       syma=mapdb(ib,3)
       symb=mapdb(ib,4)
       symc=mapdb(ib,5)
c
c3.*  def dimensions
       dima=dimm(mapdb(0,1),syma)
       dimb=dimm(mapdb(0,2),symb)
       dimc=dimm(mapdb(0,3),symc)
c
c3.*  realize packing
c
       if (symb.eq.symc) then
c3.a  case syma,symb=symc
c
c3.a.*find address for a2,3
       ia2=mapia(syma,symc,1)
       ia3=mapia(syma,symb,1)
c
c3.a.*def possition of a2,3
       possa2=mapda(ia2,1)
       possa3=mapda(ia3,1)
c
c3.a.*def additional dimensions bc
       nhelp1=dimb*(dimb-1)/2
c
c3.a.*do packing
       call t3aphlp8 (wrk(possa2),wrk(possa3),wrk(possb),
     & dima,dimb,nhelp1,ns,szkey)
c
       else
c3.b  case syma,symb.ne.symc
c
c3.b.*find address for a2,3
       ia2=mapia(syma,symc,1)
       ia3=mapia(syma,symb,1)
c
c3.b.*def possition of a2,3
       possa2=mapda(ia2,1)
       possa3=mapda(ia3,1)
c
c3.b.*do packing
       call t3aphlp7 (wrk(possa2),wrk(possa3),wrk(possb),
     & dima,dimb,dimc,ns,szkey)
c
       end if
c
 300    continue
c
c
       else
c     RC=5 , typap is not 1,2,3 (NCI)
       rc=5
       return
c
       end if
c
       return
c Avoid unused argument warnings
       if (.false.) then
         call Unused_integer(ssa)
         call Unused_integer_array(mapib)
       end if
       end
c
c     -----------------------------
c
       subroutine t3aphlp1 (a1,a2,a3,b,dimp,dimq,dimr,ns,szkey)
c
c     this routine realize following procedure
c     for symp.ne.symq.ne.symr
c
c     B(p,q,r)=B(p,q,r)+ns.(A1(q,r,p)-A2(p,r,q)+A3(p,q,r))
c
c     a1     - A1 matrix (I)
c     a2     - A2 matrix (I)
c     a3     - A3 matrix (I)
c     b      - B  matrix (I/O)
c     dimp   - dimension of p index (I)
c     dimq   - dimension of q index (I)
c     dimr   - dimension of r index (I)
c     ns     - singum of the permutation (+-1) (I)
c     szkey  - set zero key (I)
c     = 0 no vanishing
c     = 1 set B=0 at the beggining
c
       integer dimp,dimq,dimr,ns,szkey
       real*8 a1(1:dimq,1:dimr,1:dimp)
       real*8 a2(1:dimp,1:dimr,1:dimq)
       real*8 a3(1:dimp,1:dimq,1:dimr)
       real*8 b(1:dimp,1:dimq,1:dimr)
c
c     help variables
c
       integer p,q,r
       integer nhelp
c
c
       if (szkey.eq.1) then
       nhelp=dimp*dimq*dimr
       call cct3_mv0zero (nhelp,nhelp,b)
       end if
c
       if (ns.eq.1) then
c     phase +1
c
       do 102 r=1,dimr
       do 1020 q=1,dimq
       do 1021 p=1,dimp
       b(p,q,r)=b(p,q,r)+a3(p,q,r)
 1021   continue
 1020   continue
 102    continue
c
       do 104 r=1,dimr
       do 1040 q=1,dimq
       do 1041 p=1,dimp
       b(p,q,r)=b(p,q,r)-a2(p,r,q)
 1041   continue
 1040   continue
 104    continue
c
       do 106 r=1,dimr
       do 1060 q=1,dimq
       do 1061 p=1,dimp
       b(p,q,r)=b(p,q,r)+a1(q,r,p)
 1061   continue
 1060   continue
 106    continue
c
       else
c     phase -1
c
       do 202 r=1,dimr
       do 2020 q=1,dimq
       do 2021 p=1,dimp
       b(p,q,r)=b(p,q,r)-a3(p,q,r)
 2021   continue
 2020   continue
 202    continue
c
       do 204 r=1,dimr
       do 2040 q=1,dimq
       do 2041 p=1,dimp
       b(p,q,r)=b(p,q,r)+a2(p,r,q)
 2041   continue
 2040   continue
 204    continue
c
       do 206 r=1,dimr
       do 2060 q=1,dimq
       do 2061 p=1,dimp
       b(p,q,r)=b(p,q,r)-a1(q,r,p)
 2061   continue
 2060   continue
 206    continue
c
       end if
c
       return
       end
c
c     ----------------
c
       subroutine t3aphlp2 (a1,a2,a3,b,dimp,dimq,dimr,dimpq,ns,szkey)
c
c     this routine realize following procedure
c     for symp=symq.ne.symr
c
c     B(pq,r)=B(pq,r)+ns.(A1(q,r,p)-A2(p,r,q)+A3(pq,r))
c
c     a1     - A1 matrix (I)
c     a2     - A2 matrix (I)
c     a3     - A3 matrix (I)
c     b      - B  matrix (I/O)
c     dimp   - dimension of p index (I)
c     dimq   - dimension of q index (I)
c     dimr   - dimension of r index (I)
c     dimpq  - dimension of pq (I)
c     ns     - singum of the permutation (+-1) (I)
c     szkey  - set zero key (I)
c     = 0 no vanishing
c     = 1 set B=0 at the beggining
c
#include "t31.fh"
       integer dimp,dimq,dimr,dimpq,ns,szkey
       real*8 a1(1:dimq,1:dimr,1:dimp)
       real*8 a2(1:dimp,1:dimr,1:dimq)
       real*8 a3(1:dimpq,1:dimr)
       real*8 b(1:dimpq,1:dimr)
c
c     help variables
c
       integer p,q,r,pq,pq0
       integer nhelp
c
c
       if (szkey.eq.1) then
       nhelp=dimpq*dimr
       call cct3_mv0zero (nhelp,nhelp,b)
       end if
c
       if (ns.eq.1) then
c     phase +1
c
       do 102 r=1,dimr
       do 1020 pq=1,dimpq
       b(pq,r)=b(pq,r)+a3(pq,r)
 1020   continue
 102    continue
c
       do 104 r=1,dimr
       do 1040 p=2,dimp
       pq0=nshf(p)
       do 1041 q=1,p-1
       b(pq0+q,r)=b(pq0+q,r)-a2(p,r,q)
 1041   continue
 1040   continue
 104    continue
c
       do 106 r=1,dimr
       do 1060 p=2,dimp
       pq0=nshf(p)
       do 1061 q=1,p-1
       b(pq0+q,r)=b(pq0+q,r)+a1(q,r,p)
 1061   continue
 1060   continue
 106    continue
c
       else
c     phase -1
c
       do 202 r=1,dimr
       do 2020 pq=1,dimpq
       b(pq,r)=b(pq,r)-a3(pq,r)
 2020   continue
 202    continue
c
       do 204 r=1,dimr
       do 2040 p=2,dimp
       pq0=nshf(p)
       do 2041 q=1,p-1
       b(pq0+q,r)=b(pq0+q,r)+a2(p,r,q)
 2041   continue
 2040   continue
 204    continue
c
       do 206 r=1,dimr
       do 2060 p=2,dimp
       pq0=nshf(p)
       do 2061 q=1,p-1
       b(pq0+q,r)=b(pq0+q,r)-a1(q,r,p)
 2061   continue
 2060   continue
 206    continue
c
       end if
c
       return
       end
c
c     ----------------
c
       subroutine t3aphlp3 (a1,a2,a3,b,dimp,dimq,dimr,dimqr,ns,szkey)
c
c     this routine realize following procedure
c     for symp.ne.symq=symr
c
c     B(p,qr)=B(p,qr)+ns.(A1(qr,p)-A2(p,r,q)+A3(p,q,r))
c
c     a1     - A1 matrix (I)
c     a2     - A2 matrix (I)
c     a3     - A3 matrix (I)
c     b      - B  matrix (I/O)
c     dimp   - dimension of p index (I)
c     dimq   - dimension of q index (I)
c     dimr   - dimension of r index (I)
c     dimqr  - dimension of qr (I)
c     ns     - singum of the permutation (+-1) (I)
c     szkey  - set zero key (I)
c     = 0 no vanishing
c     = 1 set B=0 at the beggining
c
#include "t31.fh"
       integer dimp,dimq,dimr,dimqr,ns,szkey
       real*8 a1(1:dimqr,1:dimp)
       real*8 a2(1:dimp,1:dimr,1:dimq)
       real*8 a3(1:dimp,1:dimq,1:dimr)
       real*8 b(1:dimp,1:dimqr)
c
c     help variables
c
       integer p,q,r,qr,qr0
       integer nhelp
c
c
       if (szkey.eq.1) then
       nhelp=dimp*dimqr
       call cct3_mv0zero (nhelp,nhelp,b)
       end if
c
       if (ns.eq.1) then
c     phase +1
c
       do 102 q=2,dimq
       qr0=nshf(q)
       do 1020 r=1,q-1
       qr=qr0+r
       do 1021 p=1,dimp
       b(p,qr)=b(p,qr)+a3(p,q,r)
 1021   continue
 1020   continue
 102    continue
c
       do 104 q=2,dimq
       qr0=nshf(q)
       do 1040 r=1,q-1
       qr=qr0+r
       do 1041 p=1,dimp
       b(p,qr)=b(p,qr)-a2(p,r,q)
 1041   continue
 1040   continue
 104    continue
c
       do 106 qr=1,dimqr
       do 1060 p=1,dimp
       b(p,qr)=b(p,qr)+a1(qr,p)
 1060   continue
 106    continue
c
       else
c     phase -1
c
       do 202 q=2,dimq
       qr0=nshf(q)
       do 2020 r=1,q-1
       qr=qr0+r
       do 2021 p=1,dimp
       b(p,qr)=b(p,qr)-a3(p,q,r)
 2021   continue
 2020   continue
 202    continue
c
       do 204 q=2,dimq
       qr0=nshf(q)
       do 2040 r=1,q-1
       qr=qr0+r
       do 2041 p=1,dimp
       b(p,qr)=b(p,qr)+a2(p,r,q)
 2041   continue
 2040   continue
 204    continue
c
       do 206 qr=1,dimqr
       do 2060 p=1,dimp
       b(p,qr)=b(p,qr)-a1(qr,p)
 2060   continue
 206    continue
c
       end if
c
       return
       end
c
c     ----------------
c
       subroutine t3aphlp4 (a,b,dimp,dimpq,dimpqr,ns,szkey)
c
c     this routine realize following procedure
c     for symp=symq=symr
c
c     B(pqr)=B(pqr)+ns.(A1(qr,p)-A2(pr,q)+A3(pq,r))
c
c     a      - A  matrix (I)
c     b      - B  matrix (I/O)
c     dimp   - dimension of p (q,r) index (I)
c     dimpq  - dimension of pq index (I)
c     dimpqr - dimension of pqr (I)
c     ns     - singum of the permutation (+-1) (I)
c     szkey  - set zero key (I)
c     = 0 no vanishing
c     = 1 set B=0 at the beggining
c
c     N.B. tie sumacie by sa mozno dali trochu vylepsit
c
       integer dimp,dimpq,dimpqr,ns,szkey
       real*8 a(1:dimpq,1:dimp)
       real*8 b(1:dimpqr)
c
c     help variables
c
       integer p,q,r,pqr,pq,pr,qr,pq0
c
c
       if (szkey.eq.1) then
       call cct3_mv0zero (dimpqr,dimpqr,b)
       end if
c
       if (ns.eq.1) then
c     phase +1
c
       pqr=0
c
       do 100 p=3,dimp
       pq0=(p-1)*(p-2)/2
       qr=0
       do 101 q=2,p-1
       pq=pq0+q
       pr=(p-1)*(p-2)/2
       do 102 r=1,q-1
       pr=pr+1
       qr=qr+1
       pqr=pqr+1
       b(pqr)=b(pqr)+a(qr,p)-a(pr,q)+a(pq,r)
 102    continue
 101    continue
 100    continue
c
       else
c     phase -1
c
       pqr=0
c
       do 200 p=3,dimp
       pq0=(p-1)*(p-2)/2
       qr=0
       do 201 q=2,p-1
       pq=pq0+q
       pr=(p-1)*(p-2)/2
       do 202 r=1,q-1
       pr=pr+1
       qr=qr+1
       pqr=pqr+1
       b(pqr)=b(pqr)-a(qr,p)+a(pr,q)-a(pq,r)
 202    continue
 201    continue
 200    continue
c
c
       end if
c
       return
       end
c
c     ----------------
c
       subroutine t3aphlp5 (a1,a2,b,dimp,dimq,dimr,ns,szkey)
c
c     this routine realize following procedure
c     for symp.ne.symq.ne.symr
c
c     B(p,q,r)=B(p,q,r)+ns*(A1(q,r,p)-A2(p,r,q))
c
c     a1     - A1 matrix (I)
c     a2     - A2 matrix (I)
c     b      - B  matrix (I/O)
c     dimp   - dimension of p index (I)
c     dimq   - dimension of q index (I)
c     dimr   - dimension of r index (I)
c     ns     - singum of the permutation (+-1) (I)
c     szkey  - set zero key (I)
c     = 0 no vanishing
c     = 1 set B=0 at the beggining
c
       integer dimp,dimq,dimr,ns,szkey
       real*8 a1(1:dimq,1:dimr,1:dimp)
       real*8 a2(1:dimp,1:dimr,1:dimq)
       real*8 b(1:dimp,1:dimq,1:dimr)
c
c     help variables
c
       integer p,q,r
       integer nhelp
c
c
       if (szkey.eq.1) then
       nhelp=dimp*dimq*dimr
       call cct3_mv0zero (nhelp,nhelp,b)
       end if
c
       if (ns.eq.1) then
c     phase +1
c
       do 104 r=1,dimr
       do 1040 q=1,dimq
       do 1041 p=1,dimp
       b(p,q,r)=b(p,q,r)-a2(p,r,q)
 1041   continue
 1040   continue
 104    continue
c
       do 106 r=1,dimr
       do 1060 q=1,dimq
       do 1061 p=1,dimp
       b(p,q,r)=b(p,q,r)+a1(q,r,p)
 1061   continue
 1060   continue
 106    continue
c
       else
c     phase -1
c
       do 204 r=1,dimr
       do 2040 q=1,dimq
       do 2041 p=1,dimp
       b(p,q,r)=b(p,q,r)+a2(p,r,q)
 2041   continue
 2040   continue
 204    continue
c
       do 206 r=1,dimr
       do 2060 q=1,dimq
       do 2061 p=1,dimp
       b(p,q,r)=b(p,q,r)-a1(q,r,p)
 2061   continue
 2060   continue
 206    continue
c
       end if
c
       return
       end
c
c     ----------------
c
       subroutine t3aphlp6 (a1,a2,b,dimp,dimq,dimr,dimpq,ns,szkey)
c
c     this routine realize following procedure
c     for symp=symq.ne.symr
c
c     B(pq,r)=B(pq,r)+ns*(A1(q,r,p)-A2(p,r,q))
c
c     a1     - A1 matrix (I)
c     a2     - A2 matrix (I)
c     b      - B  matrix (I/O)
c     dimp   - dimension of p index (I)
c     dimq   - dimension of q index (I)
c     dimr   - dimension of r index (I)
c     dimpq  - dimension of pq (I)
c     ns     - singum of the permutation (+-1) (I)
c     szkey  - set zero key (I)
c     = 0 no vanishing
c     = 1 set B=0 at the beggining
c
#include "t31.fh"
       integer dimp,dimq,dimr,dimpq,ns,szkey
       real*8 a1(1:dimq,1:dimr,1:dimp)
       real*8 a2(1:dimp,1:dimr,1:dimq)
       real*8 b(1:dimpq,1:dimr)
c
c     help variables
c
       integer p,q,r,pq0
       integer nhelp
c
c
       if (szkey.eq.1) then
       nhelp=dimpq*dimr
       call cct3_mv0zero (nhelp,nhelp,b)
       end if
c
       if (ns.eq.1) then
c     phase +1
c
       do 104 r=1,dimr
       do 1040 p=2,dimp
       pq0=nshf(p)
       do 1041 q=1,p-1
       b(pq0+q,r)=b(pq0+q,r)-a2(p,r,q)
 1041   continue
 1040   continue
 104    continue
c
       do 106 r=1,dimr
       do 1060 p=2,dimp
       pq0=nshf(p)
       do 1061 q=1,p-1
       b(pq0+q,r)=b(pq0+q,r)+a1(q,r,p)
 1061   continue
 1060   continue
 106    continue
c
       else
c     phase -1
c
       do 204 r=1,dimr
       do 2040 p=2,dimp
       pq0=nshf(p)
       do 2041 q=1,p-1
       b(pq0+q,r)=b(pq0+q,r)+a2(p,r,q)
 2041   continue
 2040   continue
 204    continue
c
       do 206 r=1,dimr
       do 2060 p=2,dimp
       pq0=nshf(p)
       do 2061 q=1,p-1
       b(pq0+q,r)=b(pq0+q,r)-a1(q,r,p)
 2061   continue
 2060   continue
 206    continue
c
       end if
c
       return
       end
c
c     ----------------
c
       subroutine t3aphlp7 (a2,a3,b,dimp,dimq,dimr,ns,szkey)
c
c     this routine realize following procedure
c     for symp.ne.symq.ne.symr
c
c     B(p,q,r)=B(p,q,r)+ns*(-A2(p,r,q)+A3(p,q,r))
c
c     a2     - A2 matrix (I)
c     a3     - A3 matrix (I)
c     b      - B  matrix (I/O)
c     dimp   - dimension of p index (I)
c     dimq   - dimension of q index (I)
c     dimr   - dimension of r index (I)
c     ns     - singum of the permutation (+-1) (I)
c     szkey  - set zero key (I)
c     = 0 no vanishing
c     = 1 set B=0 at the beggining
c
       integer dimp,dimq,dimr,ns,szkey
       real*8 a2(1:dimp,1:dimr,1:dimq)
       real*8 a3(1:dimp,1:dimq,1:dimr)
       real*8 b(1:dimp,1:dimq,1:dimr)
c
c     help variables
c
       integer p,q,r
       integer nhelp
c
c
       if (szkey.eq.1) then
       nhelp=dimp*dimq*dimr
       call cct3_mv0zero (nhelp,nhelp,b)
       end if
c
       if (ns.eq.1) then
c     phase +1
c
       do 102 r=1,dimr
       do 1020 q=1,dimq
       do 1021 p=1,dimp
       b(p,q,r)=b(p,q,r)+a3(p,q,r)
 1021   continue
 1020   continue
 102    continue
c
       do 104 r=1,dimr
       do 1040 q=1,dimq
       do 1041 p=1,dimp
       b(p,q,r)=b(p,q,r)-a2(p,r,q)
 1041   continue
 1040   continue
 104    continue
c
       else
c     phase -1
c
       do 202 r=1,dimr
       do 2020 q=1,dimq
       do 2021 p=1,dimp
       b(p,q,r)=b(p,q,r)-a3(p,q,r)
 2021   continue
 2020   continue
 202    continue
c
       do 204 r=1,dimr
       do 2040 q=1,dimq
       do 2041 p=1,dimp
       b(p,q,r)=b(p,q,r)+a2(p,r,q)
 2041   continue
 2040   continue
 204    continue
c
       end if
c
       return
       end
c
c     ----------------
c
       subroutine t3aphlp8 (a2,a3,b,dimp,dimq,dimqr,ns,szkey)
c
c     this routine realize following procedure
c     for symp.ne.symq=symr
c
c     B(p,qr)=B(p,qr)+ns*(-A2(p,r,q)+A3(p,q,r))
c
c     a2     - A2 matrix (I)
c     a3     - A3 matrix (I)
c     b      - B  matrix (I/O)
c     dimp   - dimension of p index (I)
c     dimq   - dimension of q index (I)
c     dimqr  - dimension of qr (I)
c     ns     - singum of the permutation (+-1) (I)
c     szkey  - set zero key (I)
c     = 0 no vanishing
c     = 1 set B=0 at the beggining
c
#include "t31.fh"
       integer dimp,dimq,dimqr,ns,szkey
       real*8 a2(1:dimp,1:dimq,1:dimq)
       real*8 a3(1:dimp,1:dimq,1:dimq)
       real*8 b(1:dimp,1:dimqr)
c
c     help variables
c
       integer p,q,r,qr0,qr
       integer nhelp
c
c
       if (szkey.eq.1) then
       nhelp=dimp*dimqr
       call cct3_mv0zero (nhelp,nhelp,b)
       end if
c
       if (ns.eq.1) then
c     phase +1
c
       do 102 q=2,dimq
       qr0=nshf(q)
       do 1020 r=1,q-1
       qr=qr0+r
       do 1021 p=1,dimp
       b(p,qr)=b(p,qr)+a2(p,q,r)
 1021   continue
 1020   continue
 102    continue
c
       do 104 q=2,dimq
       qr0=nshf(q)
       do 1040 r=1,q-1
       qr=qr0+r
       do 1041 p=1,dimp
       b(p,qr)=b(p,qr)-a2(p,r,q)
 1041   continue
 1040   continue
 104    continue
c
       else
c     phase -1
c
       do 202 q=2,dimq
       qr0=nshf(q)
       do 2020 r=1,q-1
       qr=qr0+r
       do 2021 p=1,dimp
       b(p,qr)=b(p,qr)-a2(p,q,r)
 2021   continue
 2020   continue
 202    continue
c
       do 204 q=2,dimq
       qr0=nshf(q)
       do 2040 r=1,q-1
       qr=qr0+r
       do 2041 p=1,dimp
       b(p,qr)=b(p,qr)+a2(p,r,q)
 2041   continue
 2040   continue
 204    continue
c
       end if
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_real_array(a3)
       end
c
c     ----------------
c

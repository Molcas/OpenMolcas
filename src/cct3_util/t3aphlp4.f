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

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

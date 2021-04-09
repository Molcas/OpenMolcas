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
       subroutine diish1 (wrk,wrksize,
     & nind,rdiis1,mapd1,mapd2,mapd3,mapd4,
     & mapi1,mapi2,mapi3,mapi4,ndiis,szkey)
c
c     this routine upgrade rdiis1(p,q) maptix
c     rdiis1(p,q) = szkey*rdiis1(p,q) + (xp|xq)
c
c     nind   - number of indexes in vectors (i)
c     rdiis1 - overlap matrix of 1-5 vectors (O)
c     mapd1  - direct map of vector 1 (i)
c     mapd2  - direct map of vector 2 (i)
c     mapd3  - direct map of vector 3 (i)
c     mapd4  - direct map of vector 4 (i)
c     mapi1  - inverse map of vector 1 (i)
c     mapi2  - inverse map of vector 2 (i)
c     mapi3  - inverse map of vector 3 (i)
c     mapi4  - inverse map of vector 4 (i)
c     if tere is less than 5 vectors, use any map
c     ndiis  - dimension of DIIS (2-4) (I)
c     szkey  - 0 - no vanishing rdiis1 at the beggining
c     1 - vanishing rdiis1 at the beggining
c
#include "wrk.fh"
#include "ccsd1.fh"
c
       real*8 rdiis1(1:4,1:4)
       integer mapd1(0:512,1:6)
       integer mapd2(0:512,1:6)
       integer mapd3(0:512,1:6)
       integer mapd4(0:512,1:6)
       integer mapi1(1:8,1:8,1:8)
       integer mapi2(1:8,1:8,1:8)
       integer mapi3(1:8,1:8,1:8)
       integer mapi4(1:8,1:8,1:8)
c
       integer nind,ndiis,szkey
c
c     help variables
c
       integer nhelp
       real*8 scalar
       integer rc,num
c
c
       num=ndiis+1
c
       if (szkey.eq.1) then
       nhelp=4*4
       call mv0zero(nhelp,nhelp,rdiis1)
       end if
c
       if (num.gt.0) then
c     calc X11
       call multdot (wrk,wrksize,
     & nind,mapd1,mapi1,1,mapd1,mapi1,1,scalar,rc)
       rdiis1(1,1)=rdiis1(1,1)+scalar
       end if
c
       if (num.gt.1) then
c     X21
       call multdot (wrk,wrksize,
     & nind,mapd2,mapi2,1,mapd1,mapi1,1,scalar,rc)
       rdiis1(2,1)=rdiis1(2,1)+scalar
       rdiis1(1,2)=rdiis1(1,2)+scalar
c     X22
       call multdot (wrk,wrksize,
     & nind,mapd2,mapi2,1,mapd2,mapi2,1,scalar,rc)
       rdiis1(2,2)=rdiis1(2,2)+scalar
       end if
c
       if (num.gt.2) then
c     X31
       call multdot (wrk,wrksize,
     & nind,mapd3,mapi3,1,mapd1,mapi1,1,scalar,rc)
       rdiis1(3,1)=rdiis1(3,1)+scalar
       rdiis1(1,3)=rdiis1(1,3)+scalar
c     X32
       call multdot (wrk,wrksize,
     & nind,mapd3,mapi3,1,mapd2,mapi2,1,scalar,rc)
       rdiis1(3,2)=rdiis1(3,2)+scalar
       rdiis1(2,3)=rdiis1(2,3)+scalar
c     X33
       call multdot (wrk,wrksize,
     & nind,mapd3,mapi3,1,mapd3,mapi3,1,scalar,rc)
       rdiis1(3,3)=rdiis1(3,3)+scalar
       end if
c
       if (num.gt.3) then
c     X41
       call multdot (wrk,wrksize,
     & nind,mapd4,mapi4,1,mapd1,mapi1,1,scalar,rc)
       rdiis1(4,1)=rdiis1(4,1)+scalar
       rdiis1(1,4)=rdiis1(1,4)+scalar
c     X42
       call multdot (wrk,wrksize,
     & nind,mapd4,mapi4,1,mapd2,mapi2,1,scalar,rc)
       rdiis1(4,2)=rdiis1(4,2)+scalar
       rdiis1(2,4)=rdiis1(2,4)+scalar
c     X43
       call multdot (wrk,wrksize,
     & nind,mapd4,mapi4,1,mapd3,mapi3,1,scalar,rc)
       rdiis1(4,3)=rdiis1(4,3)+scalar
       rdiis1(3,4)=rdiis1(3,4)+scalar
c     X44
       call multdot (wrk,wrksize,
     & nind,mapd4,mapi4,1,mapd4,mapi4,1,scalar,rc)
       rdiis1(4,4)=rdiis1(4,4)+scalar
       end if
c
       return
       end

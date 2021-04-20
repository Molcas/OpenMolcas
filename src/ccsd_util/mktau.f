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
       subroutine mktau (wrk,wrksize,
     & mapdt2,mapit2,mapdt1a,mapit1a,mapdt1b,mapit1b,
     &                   fact,rc)
c
c     this routine do:
c     t2(abij) = t2(abij) + fact. (t1(ai).t1(bj)-t1(bi).t1(aj))
c     N.B. T24a,4b must be of type 4, T2abab of type 0
c
c     mapdt2  - direct map of T2 (I)
c     mapit2        - inverse map of T2 (I)
c     mapdt1a        - direct map of T1aa (I)
c     mapit1a        - inverse map of T1aa (I)
c     mapdt1b        - direct map of T1bb (I)
c     mapit1b        - inverse map of T1bb (I)
c     fact    - numerical factor (I)
c     rc        - return (error) code
c
c
#include "ccsd1.fh"
#include "wrk.fh"
       integer rc
       real*8 fact
c
       integer mapdt2(0:512,1:6)
       integer mapit2(1:8,1:8,1:8)
c
       integer mapdt1a(0:512,1:6)
       integer mapit1a(1:8,1:8,1:8)
c
       integer mapdt1b(0:512,1:6)
       integer mapit1b(1:8,1:8,1:8)
c
c     help variables
c
       integer posst2,posst1a,posst1b,posst11,posst12
       integer dimi,dimj,dima,dimb,dimab,dimij,syma,symb,symi,symj
       integer iit2,iit1a,iit1b,iit11,iit12
c
       rc=0
c
       if (mapdt2(0,6).eq.0) then
cI.1  T2abab case

       do 100 iit2=1,mapdt2(0,5)
c
       posst2=mapdt2(iit2,1)
       syma=mapdt2(iit2,3)
       symb=mapdt2(iit2,4)
       symi=mapdt2(iit2,5)
       symj=mapdt2(iit2,6)
       dima=nva(syma)
       dimb=nvb(symb)
       dimi=noa(symi)
       dimj=nob(symj)
       iit1a=mapit1a(syma,1,1)
       iit1b=mapit1b(symb,1,1)
       posst1a=mapdt1a(iit1a,1)
       posst1b=mapdt1b(iit1b,1)
c
       if ((syma.eq.symi).and.(symb.eq.symj).and.(mapdt2(iit2,2).gt.0))
     & then
       call mktauhelp1 (wrk(posst2),wrk(posst1a),wrk(posst1b),
     & dima,dimb,dimi,dimj,noa(symi),nob(symj),fact)
       end if
c
 100    continue
c
       else if ((mapdt2(0,6).eq.4).and.(mapdt2(0,1).eq.3)) then
cI.2  T2aaaa case
c
       do 200 iit2=1,mapdt2(0,5)
c
       posst2=mapdt2(iit2,1)
       syma=mapdt2(iit2,3)
       symb=mapdt2(iit2,4)
       symi=mapdt2(iit2,5)
       symj=mapdt2(iit2,6)
       dima=nva(syma)
       dimb=nva(symb)
       dimi=noa(symi)
       dimj=noa(symj)
       iit11=mapit1a(syma,1,1)
       iit12=mapit1a(symb,1,1)
       posst11=mapdt1a(iit11,1)
       posst12=mapdt1a(iit12,1)
c
       if ((syma.eq.symi).and.(symb.eq.symj).
     & and.(syma.ne.symj).and.(mapdt2(iit2,2).gt.0)) then
cI.2.*case T2(sym1,sym2,sym1,sym2)
c
       call mktauhelp1 (wrk(posst2),wrk(posst11),wrk(posst12),
     & dima,dimb,dimi,dimj,noa(syma),noa(symb),fact)
c
       else if ((syma.eq.symi).and.(symb.eq.symj).
     & and.(syma.eq.symj).and.(mapdt2(iit2,2).gt.0)) then
cI.2.*case T2(sym1,sym1,sym1,sym1)
c
       dimab=(dima*(dima-1))/2
       dimij=(dimi*(dimi-1))/2
       call mktauhelp2 (wrk(posst2),wrk(posst11),
     & dimab,dimij,dima,dimi,noa(syma),fact)
c
       end if
c
 200    continue
c
       else if ((mapdt2(0,6).eq.4).and.(mapdt2(0,1).eq.4)) then
cI.3  T2bbbb case
c
       do 300 iit2=1,mapdt2(0,5)
c
       posst2=mapdt2(iit2,1)
       syma=mapdt2(iit2,3)
       symb=mapdt2(iit2,4)
       symi=mapdt2(iit2,5)
       symj=mapdt2(iit2,6)
       dima=nvb(syma)
       dimb=nvb(symb)
       dimi=nob(symi)
       dimj=nob(symj)
       iit11=mapit1b(syma,1,1)
       iit12=mapit1b(symb,1,1)
       posst11=mapdt1b(iit11,1)
       posst12=mapdt1b(iit12,1)
c
       if ((syma.eq.symi).and.(symb.eq.symj).
     & and.(syma.ne.symj).and.(mapdt2(iit2,2).gt.0)) then
cI.3.*case T2(sym1,sym2,sym1,sym2)
c
       call mktauhelp1 (wrk(posst2),wrk(posst11),wrk(posst12),
     & dima,dimb,dimi,dimj,nob(syma),nob(symb),fact)
c
       else if ((syma.eq.symi).and.(symb.eq.symj).
     & and.(syma.eq.symj).and.(mapdt2(iit2,2).gt.0)) then
cI.3.*case T2(sym1,sym1,sym1,sym1)
c
       dimab=(dima*(dima-1))/2
       dimij=(dimi*(dimi-1))/2
       call mktauhelp2 (wrk(posst2),wrk(posst11),
     & dimab,dimij,dima,dimi,nob(syma),fact)
c
       end if
c
 300    continue
c
       else
cI.4  RC=1 : incorrect mapdt for T2
       rc=1
       return
       end if
c
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer_array(mapit2)
       end

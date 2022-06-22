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
       subroutine contt147 (wrk,wrksize,
     & lunt2o1,lunt2o2,lunt2o3)
c
c     this routine do contributions T14 and T17
c     T14: t1n(a,i) <- sum(me) [ T2o(ae,im) . FIII(e,m)]
c     T17: t1n(a,i) <- sum(e,m>n) [ T2o(ae,mn) . <ie||mn> ]
c
c     1. T1n(a,i)aa <- sum(m,e-aa) [ T2o(a,e,i,m)aaaa . FIII(e,m)aa ]
c     2. T1n(a,i)aa <- sum(m,e-bb) [ T2o(a,e,i,m)abab . FIII(e,m)bb ]
c     3. T1n(a,i)bb <- sum(m,e-bb) [ T2o(a,e,i,m)bbbb . FIII(e,m)bb ]
c     4. T1n(a,i)bb <- sum(m,e-aa) [ T2o(e,a,m,i)abab . FIII(e,m)aa ]
c     5. T1n(a,i)aa <- - sum(e,m>n-aaa) [ T2o(a,e,mn)aaaa  . <mn||ie>aaaa ]
c     6. T1n(a,i)aa <- - sum(e,m,n-bab) [ T2o(a,e,m,n)abab . <mn||ie>abab ]
c     7. T1n(a,i)bb <- - sum(e,m>n-bbb) [ T2o(a,e,mn)bbbb  . <mn||ie>bbbb ]
c     8. T1n(a,i)bb <- + sum(e,m,n-aab) [ T2o(e,a,m,n)abab . <mn||ie>abba ]
c
c     N.B. use and destroy : V1,V2,V3,M1
c     N.B. # of read : 3

       use Para_Info, only: MyRank
#include "ccsd2.fh"
#include "parallel.fh"
#include "wrk.fh"
       integer lunt2o1,lunt2o2,lunt2o3
c
c     help variables
c
       integer posst,rc,ssc
c
c
cpar
      if ((myRank.eq.idbaab).or.(myRank.eq.idaabb).or.
     &    (myRank.eq.idfin)) then
c15.1 read V1(cd,kl) <= T2o(cd,kl)aaaa
       call filemanager (2,lunt2o1,rc)
       call getmediate (wrk,wrksize,
     & lunt2o1,possv10,mapdv1,mapiv1,rc)
      end if
c
c
c
c1    T1n(a,i)aa <- sum(m,e-aa) [ T2o(a,e,i,m)aaaa . FIII(e,m)aa ]
c
cpar
      if (myRank.eq.idbaab) then
c
c1.1  expand V2(a,e,i,m) <= V1(ae,im)
       call expand (wrk,wrksize,
     & 4,4,mapdv1,mapiv1,1,possv20,mapdv2,mapiv2,rc)
c
c1.2  map V3(a,i,e,m) <= V2(a,e,i,m)
       call map (wrk,wrksize,
     & 4,1,3,2,4,mapdv2,mapiv2,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
c
c1.3  mult M1(a,i) <= V3(a,i,e,m) . FIII(e,m)aa
       call mult (wrk,wrksize,
     & 4,2,2,2,mapdv3,mapiv3,1,mapdf31,mapif31,1,mapdm1,
     &            mapim1,ssc,possm10,rc)
c
c1.4  add t1n(a,i)aa <- M1(a,i)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,1.0d0,mapdm1,1,mapdt13,mapit13,1,rc)
c
        end if
c
c
c5    T1n(a,i)aa <-  - sum(e,m>n-aaa) [ T2o(a,e,mn)aaaa  . <mn||ie>aaaa ]
c
cpar
      if (myRank.eq.idfin) then
c
c5.1  expand V2(a,e,mn) <= V1(ae,mn)
       call expand (wrk,wrksize,
     & 4,5,mapdv1,mapiv1,1,possv20,mapdv2,mapiv2,rc)
c
c5.2  map V3(e,mn,i) <= <ie||mn>aaaa
       call map (wrk,wrksize,
     & 4,4,1,2,3,mapdw11,mapiw11,1,mapdv3,mapiv3,possv30,
     &           posst,rc)
c
c5.3  mult M1(a,i) <= V2(a,e,mn) . V3(e,mn,i)
       call mult (wrk,wrksize,
     & 4,4,2,3,mapdv2,mapiv2,1,mapdv3,mapiv3,1,mapdm1,mapim1,
     &            ssc,possm10,rc)
c
c5.4  add t1n(a,i)aa <-  - M1(a,i)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,-1.0d0,mapdm1,1,mapdt13,mapit13,1,rc)
c
        end if
c
c
c
cpar
      if ((myRank.eq.idbaab).or.(myRank.eq.idaabb).or.
     &    (myRank.eq.idfin)) then
c37.1 read V1(cd,kl) <= T2o(cd,kl)bbbb
       call filemanager (2,lunt2o2,rc)
       call getmediate (wrk,wrksize,
     & lunt2o2,possv10,mapdv1,mapiv1,rc)
       end if
c
c
c
c3    T1n(a,i)bb <- sum(m,e-bb) [ T2o(a,e,i,m)bbbb . FIII(e,m)bb ]
c
cpar
      if (myRank.eq.idaabb) then
c
c3.1  expand V2(a,e,i,m) <= V1(ae,im)
       call expand (wrk,wrksize,
     & 4,4,mapdv1,mapiv1,1,possv20,mapdv2,mapiv2,rc)
c
c3.2  map V3(a,i,e,m) <= V2(a,e,i,m)
       call map (wrk,wrksize,
     & 4,1,3,2,4,mapdv2,mapiv2,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
c
c3.3  mult M1(a,i) <= V3(a,i,e,m) . FIII(e,m)bb
       call mult (wrk,wrksize,
     & 4,2,2,2,mapdv3,mapiv3,1,mapdf32,mapif32,1,mapdm1,
     &            mapim1,ssc,possm10,rc)
c
c3.4  add t1n(a,i)bb <- M1(a,i)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,1.0d0,mapdm1,1,mapdt14,mapit14,1,rc)
c
        end if
c
c
c7    T1n(a,i)bb <-  - sum(e,m>n-bbb) [ T2o(a,e,mn)bbbb  . <mn||ie>bbbb ]
c
cpar
      if (myRank.eq.idfin) then
c
c7.1  expand V2(a,e,mn) <= V1(ae,mn)
       call expand (wrk,wrksize,
     & 4,5,mapdv1,mapiv1,1,possv20,mapdv2,mapiv2,rc)
c
c7.2  map V3(e,mn,i) <= <ie||mn>bbbb
       call map (wrk,wrksize,
     & 4,4,1,2,3,mapdw12,mapiw12,1,mapdv3,mapiv3,possv30,
     &           posst,rc)
c
c7.3  mult M1(a,i) <= V2(a,e,mn) . V3(e,mn,i)
       call mult (wrk,wrksize,
     & 4,4,2,3,mapdv2,mapiv2,1,mapdv3,mapiv3,1,mapdm1,mapim1,
     &            ssc,possm10,rc)
c
c7.4  add t1n(a,i)bb <- - M1(a,i)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,-1.0d0,mapdm1,1,mapdt14,mapit14,1,rc)
c
        end if
c
c
c
cpar
      if ((myRank.eq.idbaab).or.(myRank.eq.idaabb).or.
     &    (myRank.eq.idfin)) then
c2468.1 read V1(c,d,k,l) <= T2o(c,d,k,l)abab
       call filemanager (2,lunt2o3,rc)
       call getmediate (wrk,wrksize,
     & lunt2o3,possv10,mapdv1,mapiv1,rc)
      end if
c
c
c
c2    T1n(a,i)aa <- sum(m,e-bb) [ T2o(a,e,i,m)abab . FIII(e,m)bb ]
c
cpar
      if (myRank.eq.idaabb) then
c
c2.1  map V3(a,i,e,m) <= V1(a,e,i,m)
       call map (wrk,wrksize,
     & 4,1,3,2,4,mapdv1,mapiv1,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
c
c2.2  mult M1(a,i) <= V3(a,i,e,m) . FIII(e,m)bb
       call mult (wrk,wrksize,
     & 4,2,2,2,mapdv3,mapiv3,1,mapdf32,mapif32,1,mapdm1,
     &            mapim1,ssc,possm10,rc)
c
c2.3  add t1n(a,i)aa <- M1(a,i)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,1.0d0,mapdm1,1,mapdt13,mapit13,1,rc)
c
        end if
c
c
c4    T1n(a,i)bb <- sum(m,e-aa) [ T2o(e,a,m,i)abab . FIII(e,m)aa ]
c
cpar
      if (myRank.eq.idbaab) then
c
c4.1  map V3(a,i,e,m) <= V1(e,a,m,i)
       call map (wrk,wrksize,
     & 4,3,1,4,2,mapdv1,mapiv1,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
c
c4.2  mult M1(a,i) <= V3(a,i,e,m) . FIII(e,m)aa
       call mult (wrk,wrksize,
     & 4,2,2,2,mapdv3,mapiv3,1,mapdf31,mapif31,1,mapdm1,
     &            mapim1,ssc,possm10,rc)
c
c4.3  add t1n(a,i)bb <- M1(a,i)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,1.0d0,mapdm1,1,mapdt14,mapit14,1,rc)
c
        end if
c
c
cpar
      if (myRank.eq.idfin) then
c
c6    T1n(a,i)aa <-  - sum(e,m,n-bab) [ T2o(a,e,m,n)abab . <mn||ie>abab ]
c
c6.1  map V3(e,m,n,i) <= <ie||mn>abab
       call map (wrk,wrksize,
     & 4,4,1,2,3,mapdw13,mapiw13,1,mapdv3,mapiv3,possv30,
     &           posst,rc)
c
c6.2  mult M1(a,i) <= V1(a,e,m,n) . V3(e,m,n,i)
       call mult (wrk,wrksize,
     & 4,4,2,3,mapdv1,mapiv1,1,mapdv3,mapiv3,1,mapdm1,mapim1,
     &            ssc,possm10,rc)
c
c6.3  add t1n(a,i)aa <- - M1(a,i)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,-1.0d0,mapdm1,1,mapdt13,mapit13,1,rc)
c
c
c
c8    T1n(a,i)bb <- + sum(e,m,n-aab) [ T2o(e,a,m,n)abab . <mn||ie>abba ]
c
c8.1  map V2(a,e,m,n) <= V1(e,a,m,n)
       call map (wrk,wrksize,
     & 4,2,1,3,4,mapdv1,mapiv1,1,mapdv2,mapiv2,possv20,posst,
     &           rc)
c
c8.2  map V3(e,m,n,i) <= <ie||mn>baab
       call map (wrk,wrksize,
     & 4,4,1,2,3,mapdw14,mapiw14,1,mapdv3,mapiv3,possv30,
     &           posst,rc)
c
c8.3  mult M1(a,i) <= V2(a,e,m,n) . V3(e,m,n,i)
       call mult (wrk,wrksize,
     & 4,4,2,3,mapdv2,mapiv2,1,mapdv3,mapiv3,1,mapdm1,mapim1,
     &            ssc,possm10,rc)
c
c8.4  add t1n(a,i)bb <-  M1(a,i)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,1.0d0,mapdm1,1,mapdt14,mapit14,1,rc)
c
cpar
        end if
c
       return
       end

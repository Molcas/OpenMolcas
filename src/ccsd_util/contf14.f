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
       subroutine contf14 (wrk,wrksize,
     & lunabij1,lunabij2,lunabij3,lunt2o1,lunt2o2,
     &                     lunt2o3)
c
c     this routine do:
c     2)  FI4   f1(a,e) <- -sum(m>n,f) [Tap(af,mn) . <ef||mn>]
c
c     N.B. use and destroy : V1,V2,V3,V4,M1
c     N.B. # of get v2o2 : 6
c     N.B. possible fusion with f24 graph
c
       use Para_Info, only: MyRank
       implicit none
#include "ccsd2.fh"
#include "parallel.fh"
#include "wrk.fh"
       integer lunabij1,lunabij2,lunabij3,lunt2o1,lunt2o2,lunt2o3
c
c     help variables
c
       integer posst,rc,ssc
c
c1    f1(a,e)aa <- -sum(m>n,f-aaa) [ Tap(a,f,mn)aaaa . <ef||mn>aaaa ]
c
cpar
      if (myRank.eq.idbaab) then
c
c1.1  read V1(af,mn) <= T2o(af,mn)aaaa
       call filemanager (2,lunt2o1,rc)
       call getmediate (wrk,wrksize,
     & lunt2o1,possv10,mapdv1,mapiv1,rc)
c
c1.2  make Tap V1(af,mn) from V1(af,mn)
       call mktau (wrk,wrksize,
     & mapdv1,mapiv1,mapdt11,mapit11,mapdt12,mapit12,0.5d0,
     &             rc)
c
c1.3  expand V2(a,f,mn) <= V1(af,mn)
       call expand (wrk,wrksize,
     & 4,5,mapdv1,mapiv1,1,possv20,mapdv2,mapiv2,rc)
c
c1.4  read V1(ef,mn) <= <ef||mn>aaaa
       call filemanager (2,lunabij1,rc)
       call getmediate (wrk,wrksize,
     & lunabij1,possv10,mapdv1,mapiv1,rc)
c
c1.5  expand V3(e,f,mn) <= V1(ef,mn)
       call expand (wrk,wrksize,
     & 4,5,mapdv1,mapiv1,1,possv30,mapdv3,mapiv3,rc)
c
c1.6  map V1(f,mn,e) <= V3(e,f,mn)
       call map (wrk,wrksize,
     & 4,4,1,2,3,mapdv3,mapiv3,1,mapdv1,mapiv1,possv10,posst,
     &           rc)
c
c1.7  mult M1(a,e) <= V2(a,f,mn) . V1(f,mn,e)
       call mult (wrk,wrksize,
     & 4,4,2,3,mapdv2,mapiv2,1,mapdv1,mapiv1,1,mapdm1,mapim1,
     &            ssc,possm10,rc)
c
c1.8  add f1(a,e)aa <- M1(a,e)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,-1.0d0,mapdm1,1,mapdf11,mapif11,1,rc)
c
       end if
c
c
c2    f1(a,e)bb <- -sum(m>n,f-bbb) [ Tap(a,f,mn)bbbb . <ef||mn>bbbb ]
c
cpar
      if (myRank.eq.idaabb) then
c
c2.1  read V1(af,mn) <= T2o(af,mn)bbbb
       call filemanager (2,lunt2o2,rc)
       call getmediate (wrk,wrksize,
     & lunt2o2,possv10,mapdv1,mapiv1,rc)
c
c2.2  make Tap V1(af,mn) from V1(af,mn)
       call mktau (wrk,wrksize,
     & mapdv1,mapiv1,mapdt11,mapit11,mapdt12,mapit12,0.5d0,
     &             rc)
c
c2.3  expand V2(a,f,mn) <= V1(af,mn)
       call expand (wrk,wrksize,
     & 4,5,mapdv1,mapiv1,1,possv20,mapdv2,mapiv2,rc)
c
c2.4  read V1(ef,mn) <= <ef||mn>bbbb
       call filemanager (2,lunabij2,rc)
       call getmediate (wrk,wrksize,
     & lunabij2,possv10,mapdv1,mapiv1,rc)
c
c2.5  expand V3(e,f,mn) <= V1(ef,mn)
       call expand (wrk,wrksize,
     & 4,5,mapdv1,mapiv1,1,possv30,mapdv3,mapiv3,rc)
c
c2.6  map V1(f,mn,e) <= V3(e,f,mn)
       call map (wrk,wrksize,
     & 4,4,1,2,3,mapdv3,mapiv3,1,mapdv1,mapiv1,possv10,posst,
     &           rc)
c
c2.7  mult M1(a,e) <= V2(a,f,mn) . V1(f,mn,e)
       call mult (wrk,wrksize,
     & 4,4,2,3,mapdv2,mapiv2,1,mapdv1,mapiv1,1,mapdm1,mapim1,
     &            ssc,possm10,rc)
c
c2.8  add f1(a,e)bb <- M1(a,e)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,-1.0d0,mapdm1,1,mapdf12,mapif12,1,rc)
c
       end if
c
c
c3    f1(a,e)aa <- -sum(m,n,f-abb) [ Tap(a,f,m,n)abab . <ef||mn>abab
c4    f1(a,e)bb <- -sum(m,n,f-aba) [ Tap(f,a,m,n)abab . <fe||mn>abab
c
cpar
      if ((myRank.eq.idbaab).or.(myRank.eq.idaabb)) then
c
c34.1 read V1(c,d,m,n) <= T2o(c,d,m,n)abab
       call filemanager (2,lunt2o3,rc)
       call getmediate (wrk,wrksize,
     & lunt2o3,possv10,mapdv1,mapiv1,rc)
c
c34.2 read V2(c,d,m,n) <= <cd||mn>abab
       call filemanager (2,lunabij3,rc)
       call getmediate (wrk,wrksize,
     & lunabij3,possv20,mapdv2,mapiv2,rc)
c
c34.3 make Tap V1(c,d,m,n) from V1(c,d,m,n)
       call mktau (wrk,wrksize,
     & mapdv1,mapiv1,mapdt11,mapit11,mapdt12,mapit12,0.5d0,
     &             rc)
c
c
c3.1  map V3(f,m,n,e) <= V2(e,f,m,n)
       call map (wrk,wrksize,
     & 4,4,1,2,3,mapdv2,mapiv2,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
c
      end if
c
cpar
      if (myRank.eq.idbaab) then
c
c3.2  mult M1(a,e) <= V1(a,f,m,n) . V3(f,m,n,e)
       call mult (wrk,wrksize,
     & 4,4,2,3,mapdv1,mapiv1,1,mapdv3,mapiv3,1,mapdm1,mapim1,
     &            ssc,possm10,rc)
c
c3.3  add f1(a,e))aa <- - M1(a,e)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,-1.0d0,mapdm1,1,mapdf11,mapif11,1,rc)
      end if
c
c
cpar
      if (myRank.eq.idaabb) then
c
c4.1  map V4(a,f,m,n) <= V1(f,a,m,n)
       call map (wrk,wrksize,
     & 4,2,1,3,4,mapdv1,mapiv1,1,mapdv4,mapiv4,possv40,posst,
     &           rc)
c
c4.2  map V3(f,m,n,e) <= V2 (f,e,m,n)
       call map (wrk,wrksize,
     & 4,1,4,2,3,mapdv2,mapiv2,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
c
c4.3  mult M1(a,e) <= V4(a,f,m,n) . V3(f,m,n,e)
       call mult (wrk,wrksize,
     & 4,4,2,3,mapdv4,mapiv4,1,mapdv3,mapiv3,1,mapdm1,mapim1,
     &            ssc,possm10,rc)
c
c4.4  add f1(a,e))bb <- - M1(a,e)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,-1.0d0,mapdm1,1,mapdf12,mapif12,1,rc)
c
        end if
c
       return
       end

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
       subroutine contf32 (wrk,wrksize,
     & lunabij1,lunabij2,lunabij3)
c
c     this routine do:
c     FIII2 f3(e,m) <- sum(n,f) [ <ef||mn> . T1o(n,f) ]
c
c     N.B use and destroy : V1,V2,M1
c     N.B # of read : 3
c
       use Para_Info, only: MyRank
#include "ccsd2.fh"
#include "parallel.fh"
#include "wrk.fh"
       integer lunabij1,lunabij2,lunabij3
c
c     help variables
c
       integer posst,rc,ssc
c
c1    f3(e,m)aa <- sum(n,f-aa) [ <ef||mn>aaaa . T1o(n,f)aa ]
c
cpar
      if (myRank.eq.idbaab) then
c
c1.1  read V1(ef,mn) <= <ef||mn>aaaa
       call filemanager (2,lunabij1,rc)
       call getmediate (wrk,wrksize,
     & lunabij1,possv10,mapdv1,mapiv1,rc)
c
c1.2  expand V2(e,f,m,n) <= V1(ef,mn)
       call expand (wrk,wrksize,
     & 4,4,mapdv1,mapiv1,1,possv20,mapdv2,mapiv2,rc)
c
c1.3  map V1(e,m,f,n) <= V2(e,f,m,n)
       call map (wrk,wrksize,
     & 4,1,3,2,4,mapdv2,mapiv2,1,mapdv1,mapiv1,possv10,posst,
     &           rc)
c
c1.4  mult M1(e,m) = V1(e,m,f,n) . T1o(f,n)aa
       call mult (wrk,wrksize,
     & 4,2,2,2,mapdv1,mapiv1,1,mapdt11,mapit11,1,mapdm1,
     &            mapim1,ssc,possm10,rc)
c
c1.5  add f3(e,m)aa <- M1(e,m)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,1.0d0,mapdm1,1,mapdf31,mapif31,1,rc)
c
        end if
c
c
c2    f3(e,m)bb <- sum(n,f-bb) [ <ef||mn>bbbb . T1o(n,f)bb ]
c
cpar
      if (myRank.eq.idaabb) then
c
c2.1  read V1(ef,mn) <= <ef||mn>bbbb
       call filemanager (2,lunabij2,rc)
       call getmediate (wrk,wrksize,
     & lunabij2,possv10,mapdv1,mapiv1,rc)
c
c2.2  expand V2(e,f,m,n) <= V1(ef,mn)
       call expand (wrk,wrksize,
     & 4,4,mapdv1,mapiv1,1,possv20,mapdv2,mapiv2,rc)
c
c2.3  map V1(e,m,f,n) <= V2(e,f,m,n)
       call map (wrk,wrksize,
     & 4,1,3,2,4,mapdv2,mapiv2,1,mapdv1,mapiv1,possv10,posst,
     &           rc)
c
c2.4  mult M1(e,m) = V1(e,m,f,n) . T1o(f,n)bb
       call mult (wrk,wrksize,
     & 4,2,2,2,mapdv1,mapiv1,1,mapdt12,mapit12,1,mapdm1,
     &            mapim1,ssc,possm10,rc)
c
c2.5  add f3(e,m)bb <- M1(e,m)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,1.0d0,mapdm1,1,mapdf32,mapif32,1,rc)
c
       end if
c
c
c3    f3(e,m)aa <- sum(n,f-bb) [ <ef||mn>abab . T1o(n,f)bb ]
c4    f3(e,m)bb <- sum(n,f-aa) [ <fe||nm>abab . T1o(n,f)aa ]
c
cpar
      if ((myRank.eq.idbaab).or.(myRank.eq.idaabb)) then
c34.1 read V1(c,d,k,l) <= <cd||kl>abab
       call filemanager (2,lunabij3,rc)
       call getmediate (wrk,wrksize,
     & lunabij3,possv10,mapdv1,mapiv1,rc)
       end if
c
c
cpar
      if (myRank.eq.idbaab) then
c
c3.1  map V2(e,m,f,n) <= V1 (e,f,m,n)
       call map (wrk,wrksize,
     & 4,1,3,2,4,mapdv1,mapiv1,1,mapdv2,mapiv2,possv20,posst,
     &           rc)
c
c3.2  mult M1(e,m) = V2(e,m,f,n) . T1o(f,n)bb
       call mult (wrk,wrksize,
     & 4,2,2,2,mapdv2,mapiv2,1,mapdt12,mapit12,1,mapdm1,
     &            mapim1,ssc,possm10,rc)
c
c3.3  add f3(e,m)aa <- M1(e,m)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,1.0d0,mapdm1,1,mapdf31,mapif31,1,rc)
c
       end if
c
c
cpar
      if (myRank.eq.idaabb) then
c
c4.1  map V2(e,m,f,n) <= V1 (f,e,n,m)
       call map (wrk,wrksize,
     & 4,3,1,4,2,mapdv1,mapiv1,1,mapdv2,mapiv2,possv20,posst,
     &           rc)
c
c4.2  mult M1(e,m) = V2(e,m,f,n) . T1o(f,n)aa
       call mult (wrk,wrksize,
     & 4,2,2,2,mapdv2,mapiv2,1,mapdt11,mapit11,1,mapdm1,
     &            mapim1,ssc,possm10,rc)
c
c4.3  add f3(e,m)bb <- M1(e,m)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,1.0d0,mapdm1,1,mapdf32,mapif32,1,rc)
c
       end if
c
       return
       end

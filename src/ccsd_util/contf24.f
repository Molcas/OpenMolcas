!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
       subroutine contf24 (wrk,wrksize,                                 &
     & lunabij1,lunabij2,lunabij3,lunt2o1,lunt2o2,                      &
     &                     lunt2o3)
!
!     this routine do:
!     5)  FII4  f2(m,i) <- sum(n,e>f) [ <ef||mn> . Tap(in,ef)]
!
!     N.B. use and destroy : V1,V2,V3,V4,M1
!     N.B. # of get v2o2 : 6
!     N.B. possible fusion with f14 graph
!
       use Para_Info, only: MyRank
#include "ccsd2.fh"
#include "parallel.fh"
#include "wrk.fh"
       integer lunabij1,lunabij2,lunabij3,lunt2o1,lunt2o2,lunt2o3
!
!     help variables
!
       integer posst,rc,ssc
!
!1    f2(m,i)aa <- sum(n,e>f-aaa) [ <ef||mn>aaaa . Tap(ef,in)aaaa ]
!
!par
      if (myRank.eq.idbaab) then
!
!1.1  read V1(af,mn) <= T2o(ef,in)aaaa
       call filemanager (2,lunt2o1,rc)
       call getmediate (wrk,wrksize,                                    &
     & lunt2o1,possv10,mapdv1,mapiv1,rc)
!
!1.2  make Tap V1(ef,in) from V1(ef,in)
       call mktau (wrk,wrksize,                                         &
     & mapdv1,mapiv1,mapdt11,mapit11,mapdt12,mapit12,0.5d0,             &
     &             rc)
!
!1.3  expand V2(ef,i,n) <= V1(ef,in)
       call expand (wrk,wrksize,                                        &
     & 4,6,mapdv1,mapiv1,1,possv20,mapdv2,mapiv2,rc)
!
!1.4  map V4(ef,n,i) <= V2(ef,i,n)
       call map (wrk,wrksize,                                           &
     & 4,1,2,4,3,mapdv2,mapiv2,1,mapdv4,mapiv4,possv40,posst,           &
     &           rc)
!
!1.5  read V1(ef,mn) <= <ef||mn>aaaa
       call filemanager (2,lunabij1,rc)
       call getmediate (wrk,wrksize,                                    &
     & lunabij1,possv10,mapdv1,mapiv1,rc)
!
!1.6  expand V3(ef,m,n) <= V1(ef,mn)
       call expand (wrk,wrksize,                                        &
     & 4,6,mapdv1,mapiv1,1,possv30,mapdv3,mapiv3,rc)
!
!1.7  map V1(m,ef,n) <= V3(ef,m,n)
       call map (wrk,wrksize,                                           &
     & 4,2,3,1,4,mapdv3,mapiv3,1,mapdv1,mapiv1,possv10,posst,           &
     &           rc)
!
!1.8  mult M1(m,i) <= V1(m,ef,n) . V4(ef,n,i)
       call mult (wrk,wrksize,                                          &
     & 4,4,2,3,mapdv1,mapiv1,1,mapdv4,mapiv4,1,mapdm1,mapim1,           &
     &            ssc,possm10,rc)
!
!1.9  add f2(m,i)aa <- M1(m,i)
       call add (wrk,wrksize,                                           &
     & 2,2,0,0,0,0,1,1,1.0d0,mapdm1,1,mapdf21,mapif21,1,rc)
!
       end if
!
!
!
!2    f2(m,i)bb <- sum(n,e>f-bbb) [ <ef||mn>bbbb . Tap(ef,in)bbbb ]
!
!par
      if (myRank.eq.idaabb) then
!
!2.1  read V1(af,mn) <= T2o(ef,in)bbbb
       call filemanager (2,lunt2o2,rc)
       call getmediate (wrk,wrksize,                                    &
     & lunt2o2,possv10,mapdv1,mapiv1,rc)
!
!2.2  make Tap V1(ef,in) from V1(ef,in)
       call mktau (wrk,wrksize,                                         &
     & mapdv1,mapiv1,mapdt11,mapit11,mapdt12,mapit12,0.5d0,             &
     &             rc)
!
!2.3  expand V2(ef,i,n) <= V1(ef,in)
       call expand (wrk,wrksize,                                        &
     & 4,6,mapdv1,mapiv1,1,possv20,mapdv2,mapiv2,rc)
!
!2.4  map V4(ef,n,i) <= V2(ef,i,n)
       call map (wrk,wrksize,                                           &
     & 4,1,2,4,3,mapdv2,mapiv2,1,mapdv4,mapiv4,possv40,posst,           &
     &           rc)
!
!2.5  read V1(ef,mn) <= <ef||mn>bbbb
       call filemanager (2,lunabij2,rc)
       call getmediate (wrk,wrksize,                                    &
     & lunabij2,possv10,mapdv1,mapiv1,rc)
!
!2.6  expand V3(ef,m,n) <= V1(ef,mn)
       call expand (wrk,wrksize,                                        &
     & 4,6,mapdv1,mapiv1,1,possv30,mapdv3,mapiv3,rc)
!
!2.7  map V1(m,ef,n) <= V3(ef,m,n)
       call map (wrk,wrksize,                                           &
     & 4,2,3,1,4,mapdv3,mapiv3,1,mapdv1,mapiv1,possv10,posst,           &
     &           rc)
!
!2.8  mult M1(m,i) <= V1(m,ef,n) . V4(ef,n,i)
       call mult (wrk,wrksize,                                          &
     & 4,4,2,3,mapdv1,mapiv1,1,mapdv4,mapiv4,1,mapdm1,mapim1,           &
     &            ssc,possm10,rc)
!
!2.9  add f2(m,i)bb <- M1(m,i)
       call add (wrk,wrksize,                                           &
     & 2,2,0,0,0,0,1,1,1.0d0,mapdm1,1,mapdf22,mapif22,1,rc)
!
       end if
!
!
!3*   f2(m,i)aa <- sum(nef-bab) [ <ef||mn>abab . Tap(e,f,i,n)abab ]
!4*   f2(m,i)bb <- sum(nef-aab) [ <ef||nm>abab . Tap(e,f,n,i)abab ]
!
!par
      if ((myRank.eq.idbaab).or.(myRank.eq.idaabb)) then
!
!34.1 read V1(e,f,k,l) <= T2o(e,f,k,l)abab
       call filemanager (2,lunt2o3,rc)
       call getmediate (wrk,wrksize,                                    &
     & lunt2o3,possv10,mapdv1,mapiv1,rc)
!
!34.2 read V2(e,f,k,l) <= <ef||kl>abab
       call filemanager (2,lunabij3,rc)
       call getmediate (wrk,wrksize,                                    &
     & lunabij3,possv20,mapdv2,mapiv2,rc)
!
!34.3 make Tap V1(e,f,k,l) from V1(e,f,k,l)
       call mktau (wrk,wrksize,                                         &
     & mapdv1,mapiv1,mapdt11,mapit11,mapdt12,mapit12,0.5d0,             &
     &             rc)
       end if
!
!par
      if (myRank.eq.idbaab) then
!
!3.1  map V3(m,e,f,n) <= V2(e,f,m,n)
       call map (wrk,wrksize,                                           &
     & 4,2,3,1,4,mapdv2,mapiv2,1,mapdv3,mapiv3,possv30,posst,           &
     &           rc)
!
!3.2  map V4(e,f,n,i) <= V1(e,f,i,n)
       call map (wrk,wrksize,                                           &
     & 4,1,2,4,3,mapdv1,mapiv1,1,mapdv4,mapiv4,possv40,posst,           &
     &           rc)
!
!3.3  mult M1(m,i) <= V3(m,e,f,n) . V4(e,f,n,i)
       call mult (wrk,wrksize,                                          &
     & 4,4,2,3,mapdv3,mapiv3,1,mapdv4,mapiv4,1,mapdm1,mapim1,           &
     &            ssc,possm10,rc)
!
!3.4  add f2(m,i)aa <- M1(m,i)
       call add (wrk,wrksize,                                           &
     & 2,2,0,0,0,0,1,1,1.0d0,mapdm1,1,mapdf21,mapif21,1,rc)
!
       end if
!
!par
      if (myRank.eq.idaabb) then
!
!4.1  map V3(m,e,f,n) <= V2(e,f,n,m)
       call map (wrk,wrksize,                                           &
     & 4,2,3,4,1,mapdv2,mapiv2,1,mapdv3,mapiv3,possv30,posst,           &
     &           rc)
!
!4.3  mult M1(m,i) <= V3(m,e,f,n) . V1(e,f,n,i)
       call mult (wrk,wrksize,                                          &
     & 4,4,2,3,mapdv3,mapiv3,1,mapdv1,mapiv1,1,mapdm1,mapim1,           &
     &            ssc,possm10,rc)
!
!4.4  add f2(m,i)bb <- M1(m,i)
       call add (wrk,wrksize,                                           &
     & 2,2,0,0,0,0,1,1,1.0d0,mapdm1,1,mapdf22,mapif22,1,rc)
!
       end if
!
       return
       end

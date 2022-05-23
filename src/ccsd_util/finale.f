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
c
c     this file contains following routines :
c     finale
c     contf12
c     contf14
c     contf22
c     contf23
c     contf24
c     contf32
c     contw1
c     contt12
c     contt13
c     contt147
c     contt29
c     contf4
c     contf5
c
c     -----------------------------------------------------------
c
       subroutine finale (wrk,wrksize,
     & lunabij1,lunabij2,lunabij3,
     & lunt2o1,lunt2o2,lunt2o3)
c
c     this routine do
c     1)  FI2   f1(a,e) <- -0.5 sum(m) [t1o(a,m) . fok(e,m)]
c     2)  FI4   f1(a,e) <- -sum(m>n,f) [Tap(af,mn) . <ef||mn>]
c
c     3)  FII2  f2(m,i) <- -0.5 sum(m) [fok(m,e) . T1o(e,i)]
c     4)  FII3  f2(m,i) <- sum(e,n) [ <ie||mn> . T1o(e,n)]
c     5)  FII4  f2(m,i) <- sum(n,e>f) [ <ef||mn> . Tap(in,ef)]
c
c     6)  FIII2 f3(e,m) <- sum(n,f) [ <ef||mn> . T1o(n,f) ]
c
c     7)  FIV1 FIV(b,e) <= FI(b,e)
c     8)  FIV2 FIV(b,e) <- -0.5 sum(m) [ T1o(b,m) . FIII(e,m)]
c     9)  T22 T2n(ab,ij) <- P(ab) sum(e) [t2o(ae,ij) . F4(b.e)]
c
c     10) FV1 FV(m,j) <= FII(m,j)
c     11) FV2 FV(m,j)aa  <- -0.5 sum (e) [ FIII(e,m) . T1o(e,j) ]
c     12) T23 T2n(ab,ij) <- - P(ij) sum(m) [T2o(ab,im) . F5(m,j) ]
c
c     13) WI1 w1(mnij) <= <mn||ij>
c     14) W12 w1(mnij) <- P(ij) sum(e) [ <ie||mn> . T1o(e,j) ]
c     15) WI3 w1(mnij) <- sum(e>f) [ <ef||mn> . Tau(ef,ij) ]
c     16) T24 t2n(ab,ij) <- sum(m>n) [ Tau(ab,mn) . W1(mn,ij) ]
c
c     17) T12 t1n(a,i) <- sum(e) [FI(a,e) . T1o(e,i)]
c     18) T13 t1n(a,i) <- - sum(m) [T1o(a,m) . FII(m,i)]
c     19) T14 t1n(a,i) <- sum(me) [ T2o(ae,im) . FIII(e,m)]
c     20) T17 t1n(a,i) <- sum(e,m>n) [ T2o(ae,mn) . <ie||mn> ]
c
c     ..) T22 done ad 9)
c     ..) T23 done ad 12)
c     ..) T24 done as 16)
c     21) T29 t2n(ab,ij) <- -P(a,b) sum(m) [T1o(a,m) . <mb||ij>]
c
c
#include "ccsd1.fh"
#include "ccsd2.fh"
#include "wrk.fh"
c
       integer lunabij1,lunabij2,lunabij3
       integer lunt2o1,lunt2o2,lunt2o3
c
c1
       call contf12 (wrk,wrksize)
c2
       call contf14 (wrk,wrksize,
     & lunabij1,lunabij2,lunabij3,lunt2o1,lunt2o2,lunt2o3)
c3
       call contf22 (wrk,wrksize)
c4
       call contf23 (wrk,wrksize)
c5
       call contf24 (wrk,wrksize,
     & lunabij1,lunabij2,lunabij3,lunt2o1,lunt2o2,lunt2o3)
c6
       call contf32 (wrk,wrksize,
     & lunabij1,lunabij2,lunabij3)
c7,8,9
       call contf4 (wrk,wrksize,
     & lunt2o1,lunt2o2,lunt2o3)
c10,11,12
       call contf5 (wrk,wrksize,
     & lunt2o1,lunt2o2,lunt2o3)
c13,14,15,16
       call contw1 (wrk,wrksize,
     & lunabij1,lunabij2,lunabij3,lunt2o1,lunt2o2,lunt2o3)
c17
       call contt12 (wrk,wrksize)
c18
       call contt13 (wrk,wrksize)
c19,20
       call contt147 (wrk,wrksize,
     & lunt2o1,lunt2o2,lunt2o3)
c21
       call contt29 (wrk,wrksize)
c
c
       return
       end
c
c     -----------------------------------
c
       subroutine contf12 (wrk,wrksize)
c
c     this routine do
c     FI2   f1(a,e) <- -0.5 sum(m) [t1o(a,m) . fok(e,m)]
c
c     N.B. use and destroy: M1,M2
c
        use Para_Info, only: MyRank
#include "ccsd1.fh"
#include "ccsd2.fh"
#include "parallel.fh"
#include "wrk.fh"
c
c     help variables
c
       integer rc,posst,ssc
c
c1    f1(a,e)aa <- sum(m-a) [T1o(a,m)aa . fok(e,m)aa]
c
cpar
      if (myRank.eq.idbaab) then
c
c1.1  map M1(m,e) <- fok(e,m)aa
       call map (wrk,wrksize,
     & 2,2,1,0,0,mapdfk3,mapifk3,1,mapdm1,mapim1,possm10,
     &           posst,rc)
c1.2  mult M2(a,e) = t1o(a,m)aa . M1(m,e)
       call mult (wrk,wrksize,
     & 2,2,2,1,mapdt11,mapit11,1,mapdm1,mapim1,1,mapdm2,
     &            mapim2,ssc,possm20,rc)
c1.3  add f1(a,e)aa <- -0.5 M2(a,e)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,-0.5d0,mapdm2,1,mapdf11,mapif11,1,rc)
c
       end if
c
c
c
c2    f1(a,e)bb <- sum(m-b) [T1o(a,m)bb . fok(e,m)bb]
c
cpar
      if (myRank.eq.idaabb) then
c
c2.1  map M1(m,e) <- fok(e,m)bb
       call map (wrk,wrksize,
     & 2,2,1,0,0,mapdfk4,mapifk4,1,mapdm1,mapim1,possm10,
     &           posst,rc)
c2.2  mult M2(a,e) = t1o(a,m)bb . M1(m,e)
       call mult (wrk,wrksize,
     & 2,2,2,1,mapdt12,mapit12,1,mapdm1,mapim1,1,mapdm2,
     &            mapim2,ssc,possm20,rc)
c2.3  add f1(a,e)bb <- -0.5 M2(a,e)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,-0.5d0,mapdm2,1,mapdf12,mapif12,1,rc)
c
        end if
c
       return
       end
c
c     -----------------------------------
c
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
c
c     -----------------------------------
c
       subroutine contf22 (wrk,wrksize)
c
c     this routine do:
c     f2(m,i) <- 0.5 sum(m) [fok(m,e) . T1o(e,i)]
c
c     N.B. use and destroy : M1,M2
c
       use Para_Info, only: MyRank
#include "ccsd2.fh"
#include "parallel.fh"
#include "wrk.fh"
c
c     help variables
c
       integer posst,rc,ssc
c
c1    f2(m,i)aa <- 0.5 sum(e-a) [ fok(e,m)aa . t1o(e,i)aa]
c
cpar
      if (myRank.eq.idbaab) then
c
c1.1  map M1(m,e) <= fok(e,m)aa
       call map (wrk,wrksize,
     & 2,2,1,0,0,mapdfk3,mapifk3,1,mapdm1,mapim1,possm10,
     &           posst,rc)
c
c1.2  mult M2(m,i) <= M1(m,e) . T1o(e,i)aa
       call mult (wrk,wrksize,
     & 2,2,2,1,mapdm1,mapim1,1,mapdt11,mapit11,1,mapdm2,
     &            mapim2,ssc,possm20,rc)
c
c1.3  add f2(m,i)aa <- 0.5 M2(m,i)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,0.5d0,mapdm2,1,mapdf21,mapif21,1,rc)
c
       end if
c
c
c2    f2(m,i)aa <- 0.5 sum(e-b) [ fok(e,m)bb . t1o(e,i)bb]
c
cpar
      if (myRank.eq.idaabb) then
c
c2.1  map M1(m,e) <= fok(e,m)bb
       call map (wrk,wrksize,
     & 2,2,1,0,0,mapdfk4,mapifk4,1,mapdm1,mapim1,possm10,
     &           posst,rc)
c
c2.2  mult M2(m,i) <= M1(m,e) . T1o(e,i)bb
       call mult (wrk,wrksize,
     & 2,2,2,1,mapdm1,mapim1,1,mapdt12,mapit12,1,mapdm2,
     &            mapim2,ssc,possm20,rc)
c
c2.3  add f2(m,i)bb <- 0.5 M2(m,i)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,0.5d0,mapdm2,1,mapdf22,mapif22,1,rc)
c
        end if
c
       return
       end
c
c     -----------------------------------
c
       subroutine contf23 (wrk,wrksize)
c
c     this routine do:
c     f2(m,i) <- sum(e,n) [ <ie||mn> . T1o(e,n)]
c
       use Para_Info, only: MyRank
#include "ccsd2.fh"
#include "parallel.fh"
#include "wrk.fh"
c
c     help variables
c
       integer posst,rc,ssc
c
c1    f2(m,i)aa <- sum(e,n-aa) [ <ie||mn>aaaa . t1o(e,n)aa ]
c
cpar
      if (myRank.eq.idbaab) then
c
c1.1  expand V1(i,e,m,n) <= <ie||mn>aaaa
       call expand (wrk,wrksize,
     & 4,3,mapdw11,mapiw11,1,possv10,mapdv1,mapiv1,rc)
c
c1.2  map V2(m,i,e,n) <= V1(i,e,m,n)
       call map (wrk,wrksize,
     & 4,2,3,1,4,mapdv1,mapiv1,1,mapdv2,mapiv2,possv20,posst,
     &           rc)
c
c1.3  mult M1(m,i) = V2(m,i,e,n) . T2o(e,n)aa
       call mult (wrk,wrksize,
     & 4,2,2,2,mapdv2,mapiv2,1,mapdt11,mapit11,1,mapdm1,
     &            mapim1,ssc,possm10,rc)
c
c1.4  add f2(m,i)aa <- M1(m,i)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,1.0d0,mapdm1,1,mapdf21,mapif21,1,rc)
c
c
c
c2    f2(m,i)aa <- sum(e,n-bb) [ <ie||mn>abab . t1o(e,n)bb ]
c
c2.1  map V2(m,i,e,n) <= <ie||mn>abab
       call map (wrk,wrksize,
     & 4,2,3,1,4,mapdw13,mapiw13,1,mapdv2,mapiv2,possv20,
     &           posst,rc)
c
c2.2  mult M1(m,i) = V2(m,i,e,n) . T2o(e,n)bb
       call mult (wrk,wrksize,
     & 4,2,2,2,mapdv2,mapiv2,1,mapdt12,mapit12,1,mapdm1,
     &            mapim1,ssc,possm10,rc)
c
c2.3  add f2(m,i)aa <- M1(m,i)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,1.0d0,mapdm1,1,mapdf21,mapif21,1,rc)
c
       end if
c
c
c3    f2(m,i)bb <- sum(e,n-bb) [ <ie||mn>bbbb . t1o(e,n)bb ]
c
cpar
      if (myRank.eq.idaabb) then
c
c3.1  expand V1(i,e,m,n) <= <ie||mn>bbbb
       call expand (wrk,wrksize,
     & 4,3,mapdw12,mapiw12,1,possv10,mapdv1,mapiv1,rc)
c
c3.2  map V2(m,i,e,n) <= V1(i,e,m,n)
       call map (wrk,wrksize,
     & 4,2,3,1,4,mapdv1,mapiv1,1,mapdv2,mapiv2,possv20,posst,
     &           rc)
c
c3.3  mult M1(m,i) =  V2(m,i,e,n) . T2o(e,n)bb
       call mult (wrk,wrksize,
     & 4,2,2,2,mapdv2,mapiv2,1,mapdt12,mapit12,1,mapdm1,
     &            mapim1,ssc,possm10,rc)
c
c3.4  add f2(m,i)bb <- M1(m,i)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,1.0d0,mapdm1,1,mapdf22,mapif22,1,rc)
c
c
c
c4    f2(m,i)bb <- - sum(e,n-aa) [ <ie||nm>baab . t1o(e,n)aa ]
c
c4.1  map V2(m,i,e,n) <= <ie||nm>baab
       call map (wrk,wrksize,
     & 4,2,3,4,1,mapdw14,mapiw14,1,mapdv2,mapiv2,possv20,
     &           posst,rc)
c
c4.2  mult M1(m,i) = V2(m,i,e,n) . T2o(e,n)aa
       call mult (wrk,wrksize,
     & 4,2,2,2,mapdv2,mapiv2,1,mapdt11,mapit11,1,mapdm1,
     &            mapim1,ssc,possm10,rc)
c
c4.3  add f2(m,i)aa <- M1(m,i)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,-1.0d0,mapdm1,1,mapdf22,mapif22,1,rc)
c
       end if
c
       return
       end
c
c     -----------------------------------
c
       subroutine contf24 (wrk,wrksize,
     & lunabij1,lunabij2,lunabij3,lunt2o1,lunt2o2,
     &                     lunt2o3)
c
c     this routine do:
c     5)  FII4  f2(m,i) <- sum(n,e>f) [ <ef||mn> . Tap(in,ef)]
c
c     N.B. use and destroy : V1,V2,V3,V4,M1
c     N.B. # of get v2o2 : 6
c     N.B. possible fusion with f14 graph
c
       use Para_Info, only: MyRank
#include "ccsd2.fh"
#include "parallel.fh"
#include "wrk.fh"
       integer lunabij1,lunabij2,lunabij3,lunt2o1,lunt2o2,lunt2o3
c
c     help variables
c
       integer posst,rc,ssc
c
c1    f2(m,i)aa <- sum(n,e>f-aaa) [ <ef||mn>aaaa . Tap(ef,in)aaaa ]
c
cpar
      if (myRank.eq.idbaab) then
c
c1.1  read V1(af,mn) <= T2o(ef,in)aaaa
       call filemanager (2,lunt2o1,rc)
       call getmediate (wrk,wrksize,
     & lunt2o1,possv10,mapdv1,mapiv1,rc)
c
c1.2  make Tap V1(ef,in) from V1(ef,in)
       call mktau (wrk,wrksize,
     & mapdv1,mapiv1,mapdt11,mapit11,mapdt12,mapit12,0.5d0,
     &             rc)
c
c1.3  expand V2(ef,i,n) <= V1(ef,in)
       call expand (wrk,wrksize,
     & 4,6,mapdv1,mapiv1,1,possv20,mapdv2,mapiv2,rc)
c
c1.4  map V4(ef,n,i) <= V2(ef,i,n)
       call map (wrk,wrksize,
     & 4,1,2,4,3,mapdv2,mapiv2,1,mapdv4,mapiv4,possv40,posst,
     &           rc)
c
c1.5  read V1(ef,mn) <= <ef||mn>aaaa
       call filemanager (2,lunabij1,rc)
       call getmediate (wrk,wrksize,
     & lunabij1,possv10,mapdv1,mapiv1,rc)
c
c1.6  expand V3(ef,m,n) <= V1(ef,mn)
       call expand (wrk,wrksize,
     & 4,6,mapdv1,mapiv1,1,possv30,mapdv3,mapiv3,rc)
c
c1.7  map V1(m,ef,n) <= V3(ef,m,n)
       call map (wrk,wrksize,
     & 4,2,3,1,4,mapdv3,mapiv3,1,mapdv1,mapiv1,possv10,posst,
     &           rc)
c
c1.8  mult M1(m,i) <= V1(m,ef,n) . V4(ef,n,i)
       call mult (wrk,wrksize,
     & 4,4,2,3,mapdv1,mapiv1,1,mapdv4,mapiv4,1,mapdm1,mapim1,
     &            ssc,possm10,rc)
c
c1.9  add f2(m,i)aa <- M1(m,i)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,1.0d0,mapdm1,1,mapdf21,mapif21,1,rc)
c
       end if
c
c
c
c2    f2(m,i)bb <- sum(n,e>f-bbb) [ <ef||mn>bbbb . Tap(ef,in)bbbb ]
c
cpar
      if (myRank.eq.idaabb) then
c
c2.1  read V1(af,mn) <= T2o(ef,in)bbbb
       call filemanager (2,lunt2o2,rc)
       call getmediate (wrk,wrksize,
     & lunt2o2,possv10,mapdv1,mapiv1,rc)
c
c2.2  make Tap V1(ef,in) from V1(ef,in)
       call mktau (wrk,wrksize,
     & mapdv1,mapiv1,mapdt11,mapit11,mapdt12,mapit12,0.5d0,
     &             rc)
c
c2.3  expand V2(ef,i,n) <= V1(ef,in)
       call expand (wrk,wrksize,
     & 4,6,mapdv1,mapiv1,1,possv20,mapdv2,mapiv2,rc)
c
c2.4  map V4(ef,n,i) <= V2(ef,i,n)
       call map (wrk,wrksize,
     & 4,1,2,4,3,mapdv2,mapiv2,1,mapdv4,mapiv4,possv40,posst,
     &           rc)
c
c2.5  read V1(ef,mn) <= <ef||mn>bbbb
       call filemanager (2,lunabij2,rc)
       call getmediate (wrk,wrksize,
     & lunabij2,possv10,mapdv1,mapiv1,rc)
c
c2.6  expand V3(ef,m,n) <= V1(ef,mn)
       call expand (wrk,wrksize,
     & 4,6,mapdv1,mapiv1,1,possv30,mapdv3,mapiv3,rc)
c
c2.7  map V1(m,ef,n) <= V3(ef,m,n)
       call map (wrk,wrksize,
     & 4,2,3,1,4,mapdv3,mapiv3,1,mapdv1,mapiv1,possv10,posst,
     &           rc)
c
c2.8  mult M1(m,i) <= V1(m,ef,n) . V4(ef,n,i)
       call mult (wrk,wrksize,
     & 4,4,2,3,mapdv1,mapiv1,1,mapdv4,mapiv4,1,mapdm1,mapim1,
     &            ssc,possm10,rc)
c
c2.9  add f2(m,i)bb <- M1(m,i)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,1.0d0,mapdm1,1,mapdf22,mapif22,1,rc)
c
       end if
c
c
c3*   f2(m,i)aa <- sum(nef-bab) [ <ef||mn>abab . Tap(e,f,i,n)abab ]
c4*   f2(m,i)bb <- sum(nef-aab) [ <ef||nm>abab . Tap(e,f,n,i)abab ]
c
cpar
      if ((myRank.eq.idbaab).or.(myRank.eq.idaabb)) then
c
c34.1 read V1(e,f,k,l) <= T2o(e,f,k,l)abab
       call filemanager (2,lunt2o3,rc)
       call getmediate (wrk,wrksize,
     & lunt2o3,possv10,mapdv1,mapiv1,rc)
c
c34.2 read V2(e,f,k,l) <= <ef||kl>abab
       call filemanager (2,lunabij3,rc)
       call getmediate (wrk,wrksize,
     & lunabij3,possv20,mapdv2,mapiv2,rc)
c
c34.3 make Tap V1(e,f,k,l) from V1(e,f,k,l)
       call mktau (wrk,wrksize,
     & mapdv1,mapiv1,mapdt11,mapit11,mapdt12,mapit12,0.5d0,
     &             rc)
       end if
c
cpar
      if (myRank.eq.idbaab) then
c
c3.1  map V3(m,e,f,n) <= V2(e,f,m,n)
       call map (wrk,wrksize,
     & 4,2,3,1,4,mapdv2,mapiv2,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
c
c3.2  map V4(e,f,n,i) <= V1(e,f,i,n)
       call map (wrk,wrksize,
     & 4,1,2,4,3,mapdv1,mapiv1,1,mapdv4,mapiv4,possv40,posst,
     &           rc)
c
c3.3  mult M1(m,i) <= V3(m,e,f,n) . V4(e,f,n,i)
       call mult (wrk,wrksize,
     & 4,4,2,3,mapdv3,mapiv3,1,mapdv4,mapiv4,1,mapdm1,mapim1,
     &            ssc,possm10,rc)
c
c3.4  add f2(m,i)aa <- M1(m,i)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,1.0d0,mapdm1,1,mapdf21,mapif21,1,rc)
c
       end if
c
cpar
      if (myRank.eq.idaabb) then
c
c4.1  map V3(m,e,f,n) <= V2(e,f,n,m)
       call map (wrk,wrksize,
     & 4,2,3,4,1,mapdv2,mapiv2,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
c
c4.3  mult M1(m,i) <= V3(m,e,f,n) . V1(e,f,n,i)
       call mult (wrk,wrksize,
     & 4,4,2,3,mapdv3,mapiv3,1,mapdv1,mapiv1,1,mapdm1,mapim1,
     &            ssc,possm10,rc)
c
c4.4  add f2(m,i)bb <- M1(m,i)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,1.0d0,mapdm1,1,mapdf22,mapif22,1,rc)
c
       end if
c
       return
       end
c
c     -----------------------------------
c
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
c
c     -----------------------------------
c
       subroutine contw1 (wrk,wrksize,
     & lunabij1,lunabij2,lunabij3,lunt2o1,lunt2o2,
     &                    lunt2o3)
c
c     this routine calc contributions to W1 i.e. W11,W12,W13
c     and controibution frim W1 to T2 i.e. T24
c
c
c     Mediate : WI(m,n,i,j)
c     Spin states : WI(mn,ij)aaaa, WI(mn,ij)bbbb, WI(m,n,i,j)abab
c     Contributions:
c
c     WI1
cI.1  WI(mn,ij)aaaa   <= <mn||ij>aaaa
cJ.1  WI(mn,ij)bbbb   <= <mn||ij>bbbb
cK.1  WI(m,n,i,j)abab <= <mn||ij>abab
c
c     WI2
cI.2  Q(mn,i,j)aaaa  <= sum(e-a) [ <ie||mn>aaaa . T1o(e,j)aa ]
c     WI(mn,ij)aaaa   <- Q(mn,i,j)aaaa - Q(mn,j,i)aaaa  i>j
cJ.2  Q(mn,i,j)bbbb  <= sum(e-b) [ <ie||mn>bbbb . T1o(e,j)bb ]
c     WI(mn,ij)bbbb   <- Q(mn,i,j)bbbb - Q(mn,j,i)bbbb  i>j
cK.2  WI(m,n,i,j)abab <- - sum(e-a) [ <je||mn>baab . T1o(e,i)aa ]
c     <- sum(e-b) [ <ie||mn>abab . T1o(e,j)bb ]
c
c     WI3
cI.3  WI(mn,ij)aaaa   <- sum(e>f-aa) [ <ef||mn>aaaa . Tau(ef,ij)aaaa ]
cJ.3  WI(mn,ij)bbbb   <- sum(e>f-bb) [ <ef||mn>bbbb . Tau(ef,ij)bbbb ]
cK.3  WI(m,n,i,j)abab <- sum(e,f-ab) [ <ef||mn>abab . Tau(e,f,i,j)abab ]
c
c     T24
cI.4  T2n(ab,ij)aaaa   <- sum(m>n-aa) [ Tau(ab,mn)aaaa . W1(mn,ij)aaaa ]
cJ.4  T2n(ab,ij)bbbb   <- sum(m>n-bb) [ Tau(ab,mn)bbbb . W1(mn,ij)bbbb ]
cK.4  T2n(a,b,i,j)abab <- sum(m,n-ab) [ Tau(a,b,m,n)abab . W1(m,n,i,j)abab ]
c
c     N.B. use and destroy : V1,V2,V3,V4
c     N.B. number of read  : 6
c
       use Para_Info, only: MyRank
#include "ccsd2.fh"
#include "parallel.fh"
#include "wrk.fh"
       integer lunabij1,lunabij2,lunabij3,lunt2o1,lunt2o2,lunt2o3
c
c     help variables
c
       integer posst,rc,ssc
c
c
cpar
      if (myRank.eq.idfin) then
c
cI    case W1(mn,ij)aaaa
c
cI.1.1map V1(mn,ij) <= <mn||ij>aaaa
       call map (wrk,wrksize,
     & 4,1,2,3,4,mapdw01,mapiw01,1,mapdv1,mapiv1,possv10,
     &           posst,rc)
c
c
c
cI.2.1map V2(mn,i,e) <= <ie||mn>aaaa
       call map (wrk,wrksize,
     & 4,3,4,1,2,mapdw11,mapiw11,1,mapdv2,mapiv2,possv20,
     &           posst,rc)
c
cI.2.2mult V3(mn,i,j) <= V2(mn,i,e) . T1o(e,j)aa
       call mult (wrk,wrksize,
     & 4,2,4,1,mapdv2,mapiv2,1,mapdt11,mapit11,1,mapdv3,
     &            mapiv3,ssc,possv30,rc)
c
cI.2.3pack V2(mn,ij) <= V3(mn,i,j)
       call fack (wrk,wrksize,
     & 4,4,mapdv3,1,mapiv3,mapdv2,mapiv2,possv20,rc)
c
cI.2.4add V1(mn,ij) <- V2(mn,ij)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv2,1,mapdv1,mapiv1,1,rc)
c
c
c
cI.3.1read V2(ef,ij) <= T2o(ef,ij)aaaa
       call filemanager (2,lunt2o1,rc)
       call getmediate (wrk,wrksize,
     & lunt2o1,possv20,mapdv2,mapiv2,rc)
c
cI.3.2make Tau V2(ef,ij) from V2(ef,ij)
       call mktau (wrk,wrksize,
     & mapdv2,mapiv2,mapdt11,mapit11,mapdt12,mapit12,1.0d0,
     &             rc)
c
cI.3.3read V3(ef,mn) <= <ef||mn>aaaa
       call filemanager (2,lunabij1,rc)
       call getmediate (wrk,wrksize,
     & lunabij1,possv30,mapdv3,mapiv3,rc)
c
cI.3.4map V4(mn,ef) <= V3(ef,mn)
       call map (wrk,wrksize,
     & 4,3,4,1,2,mapdv3,mapiv3,1,mapdv4,mapiv4,possv40,posst,
     &           rc)
c
cI.3.5mult V3(mn,ij) = V4(mn,ef) . V2(ef,ij)
       call mult (wrk,wrksize,
     & 4,4,4,2,mapdv4,mapiv4,1,mapdv2,mapiv2,1,mapdv3,mapiv3,
     &            ssc,possv30,rc)
c
cI.3.6add V1(mn,ij) <- V3(mn,ij)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv3,1,mapdv1,mapiv1,1,rc)
c
c
c
cI.4.0Tau(ab,mn) are in V2(ab,mn) from I.3.2
c     W1(mn,ij)  are in V1(mn,ij)
c
cI.4.1V3(ab,ij) = V2(ab,mn) . V1(mn,ij)
       call mult (wrk,wrksize,
     & 4,4,4,2,mapdv2,mapiv2,1,mapdv1,mapiv1,1,mapdv3,mapiv3,
     &            ssc,possv30,rc)
c
cI.4.2add t2n(ab,ij)aaaa <- V3(ab,ij)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv3,1,mapdt21,mapit21,1,rc)
c
c
c
cJ    case W1(mn,ij)bbbb
c
cJ.1.1map V1(mn,ij) <= <mn||ij>bbbb
       call map (wrk,wrksize,
     & 4,1,2,3,4,mapdw02,mapiw02,1,mapdv1,mapiv1,possv10,
     &           posst,rc)
c
c
c
cJ.2.1map V2(mn,i,e) <= <ie||mn>bbbb
       call map (wrk,wrksize,
     & 4,3,4,1,2,mapdw12,mapiw12,1,mapdv2,mapiv2,possv20,
     &           posst,rc)
c
cJ.2.2mult V3(mn,i,j) <= V2(mn,i,e) . T1o(e,j)bb
       call mult (wrk,wrksize,
     & 4,2,4,1,mapdv2,mapiv2,1,mapdt12,mapit12,1,mapdv3,
     &            mapiv3,ssc,possv30,rc)
c
cJ.2.3pack V2(mn,ij) <= V3(mn,i,j)
       call fack (wrk,wrksize,
     & 4,4,mapdv3,1,mapiv3,mapdv2,mapiv2,possv20,rc)
c
cJ.2.4add V1(mn,ij) <- V2(mn,ij)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv2,1,mapdv1,mapiv1,1,rc)
c
c
c
cJ.3.1read V2(ef,ij) <= T2o(ef,ij)bbbb
       call filemanager (2,lunt2o2,rc)
       call getmediate (wrk,wrksize,
     & lunt2o2,possv20,mapdv2,mapiv2,rc)
c
cJ.3.2make Tau V2(ef,ij) from V2(ef,ij)
       call mktau (wrk,wrksize,
     & mapdv2,mapiv2,mapdt11,mapit11,mapdt12,mapit12,1.0d0,
     &             rc)
c
cJ.3.3read V3(ef,mn) <= <ef||mn>bbbb
       call filemanager (2,lunabij2,rc)
       call getmediate (wrk,wrksize,
     & lunabij2,possv30,mapdv3,mapiv3,rc)
c
cJ.3.4map V4(mn,ef) <= V3(ef,mn)
       call map (wrk,wrksize,
     & 4,3,4,1,2,mapdv3,mapiv3,1,mapdv4,mapiv4,possv40,posst,
     &           rc)
c
cJ.3.5mult V3(mn,ij) = V4(mn,ef) . V2(ef,ij)
       call mult (wrk,wrksize,
     & 4,4,4,2,mapdv4,mapiv4,1,mapdv2,mapiv2,1,mapdv3,mapiv3,
     &            ssc,possv30,rc)
c
cJ.3.6add V1(mn,ij) <- V3(mn,ij)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv3,1,mapdv1,mapiv1,1,rc)
c
c
c
cJ.4.0Tau(ab,mn) are in V2(ab,mn) from J.3.2
c     W1(mn,ij)  are in V1(mn,ij)
c
cJ.4.1mult V3(ab,ij) = V2(ab,mn) . V1(mn,ij)
       call mult (wrk,wrksize,
     & 4,4,4,2,mapdv2,mapiv2,1,mapdv1,mapiv1,1,mapdv3,mapiv3,
     &            ssc,possv30,rc)
c
cJ.4.2add t2n(ab,ij)bbbb <- V3(ab,ij)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv3,1,mapdt22,mapit22,1,rc)
c
c
c
cK    case W1(mn,ij)abab
c
cK.1.1map V1(m,n,i,j) <= <mn||ij>abab
       call map (wrk,wrksize,
     & 4,1,2,3,4,mapdw03,mapiw03,1,mapdv1,mapiv1,possv10,
     &           posst,rc)
c
c
c
cK.2.1map V2(m,n,j,e) <= <je||mn>baab
       call map (wrk,wrksize,
     & 4,3,4,1,2,mapdw14,mapiw14,1,mapdv2,mapiv2,possv20,
     &           posst,rc)
c
cK.2.2mult V3(m,n,j,i) <= V2(m,n,j,e) . T1o(e,i)aa
       call mult (wrk,wrksize,
     & 4,2,4,1,mapdv2,mapiv2,1,mapdt11,mapit11,1,mapdv3,
     &            mapiv3,ssc,possv30,rc)
c
cK.2.3map V2(m,n,i,j) <= V3(m,n,j,i)
       call map (wrk,wrksize,
     & 4,1,2,4,3,mapdv3,mapiv3,1,mapdv2,mapiv2,possv20,posst,
     &           rc)
c
cK.2.4add V1(m,n,i,j) <- - V2(m,n,i,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,-1.0d0,mapdv2,1,mapdv1,mapiv1,1,rc)
c
c
cK.2.5map V2(m,n,i,e) <= <ie||mn>abab
       call map (wrk,wrksize,
     & 4,3,4,1,2,mapdw13,mapiw13,1,mapdv2,mapiv2,possv20,
     &           posst,rc)
c
cK.2.6mult V3(m,n,i,j) <= V2(m,n,i,e) . T1o(e,j)bb
       call mult (wrk,wrksize,
     & 4,2,4,1,mapdv2,mapiv2,1,mapdt12,mapit12,1,mapdv3,
     &            mapiv3,ssc,possv30,rc)
c
cK.2.7add V1(m,n,i,j) <- V3(m,n,i,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv3,1,mapdv1,mapiv1,1,rc)
c
c
c
cK.3.1read V2(e,f,i,j) <= T2o(e,f,i,j)abab
       call filemanager (2,lunt2o3,rc)
       call getmediate (wrk,wrksize,
     & lunt2o3,possv20,mapdv2,mapiv2,rc)
c
cK.3.2make Tau V2(e,f,i,j) from V2(e,f,i,j)
       call mktau (wrk,wrksize,
     & mapdv2,mapiv2,mapdt11,mapit11,mapdt12,mapit12,1.0d0,
     &             rc)
c
cK.3.3read V3(e,f,m,n) <= <ef||mn>abab
       call filemanager (2,lunabij3,rc)
       call getmediate (wrk,wrksize,
     & lunabij3,possv30,mapdv3,mapiv3,rc)
c
cK.3.4map V4(m,n,e,f) <= V3(e,f,m,n)
       call map (wrk,wrksize,
     & 4,3,4,1,2,mapdv3,mapiv3,1,mapdv4,mapiv4,possv40,posst,
     &           rc)
c
cK.3.5mult V3(m,n,i,j) = V4(m,n,e,f) . V2(e,f,i,j)
       call mult (wrk,wrksize,
     & 4,4,4,2,mapdv4,mapiv4,1,mapdv2,mapiv2,1,mapdv3,mapiv3,
     &            ssc,possv30,rc)
c
cK.3.6add V1(m,n,i,j) <- V3(m,n,i,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv3,1,mapdv1,mapiv1,1,rc)
c
c
c
cK.4.0Tau(a,b,m,n) are in V2(ab,mn) from K.3.2
c     W1(m,n,i,j)  are in V1(m,n,i,j)
c
cK.4.1mult V3(a,b,i,j) = V2(a,b,m,n) . V1(m,n,i,j)
       call mult (wrk,wrksize,
     & 4,4,4,2,mapdv2,mapiv2,1,mapdv1,mapiv1,1,mapdv3,mapiv3,
     &            ssc,possv30,rc)
c
cK.4.2add t2n(a,b,i,j)abab <- V3(a,b,i,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv3,1,mapdt23,mapit23,1,rc)
c
cpar
       end if
c
       return
       end
c
c     -----------------------------------
c
       subroutine contt12 (wrk,wrksize)
c
c     this routine do T12 contribution:
c     t1n(a,i) <- sum(e) [FI(a,e) . T1o(e,i)
c
c     N.B. use and destroy : M1
c     N.B. Parallel : in the case where idaaaa.ne.idbaab and
c          idaaaa.ne.idbaab this routine runs contributions to
c          T1n also on idaaaa and idbbbb nodes, since on these
c          nodes there is a specific
c          part of contributions F13 (see notes in sumoverb routine)
c           which are not presented on pilot nodes
c
       use Para_Info, only: MyRank
#include "ccsd2.fh"
#include "parallel.fh"
#include "wrk.fh"
c
c     help variables
c
       integer rc,ssc
c
c1    t1n(a,i)aa <- sum(e-a) [F1(a,e)aa . T1o(e,i)aa ]
c
cpar
      if ((myRank.eq.idbaab).or.(myRank.eq.idaaaa)) then
c1.1  mult M1(a,i) <= F1(a,e)aa . T1o(e,i)aa
       call mult (wrk,wrksize,
     & 2,2,2,1,mapdf11,mapif11,1,mapdt11,mapit11,1,mapdm1,
     &            mapim1,ssc,possm10,rc)
c
c1.2  add t1n(a,i)aa <- M1(a,i)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,1.0d0,mapdm1,1,mapdt13,mapit13,1,rc)
      end if
c
c
c2    t1n(a,i)bb <- sum(e-a) [F1(a,e)bb . T1o(e,i)bb ]
c
cpar
      if ((myRank.eq.idaabb).or.(myRank.eq.idbbbb)) then
c2.1  mult M1(a,i) <= F1(a,e)bb . T1o(e,i)bb
       call mult (wrk,wrksize,
     & 2,2,2,1,mapdf12,mapif12,1,mapdt12,mapit12,1,mapdm1,
     &            mapim1,ssc,possm10,rc)
c
c2.2  add t1n(a,i)bb <- M1(a,i)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,1.0d0,mapdm1,1,mapdt14,mapit14,1,rc)
      end if
c
       return
       end
c
c     -----------------------------------
c
       subroutine contt13 (wrk,wrksize)
c
c     this routine do T13 contribution:
c     t1n(a,i) <- sum(m) [T1o(a,m) . FII(m,i)]
c
c     N.B. use and destroy : M1
c
       use Para_Info, only: MyRank
#include "ccsd2.fh"
#include "parallel.fh"
#include "wrk.fh"
c
c     help variables
c
       integer rc,ssc
c
c1    t1n(a,i)aa <- - sum(m-a) [T1o(a,m)aa  . FII(m,i)aa]
c
cpar
      if (myRank.eq.idbaab) then
c1.1  mult M1(a,i) <= T1o(a,m)aa  . FII(m,i)aa
       call mult (wrk,wrksize,
     & 2,2,2,1,mapdt11,mapit11,1,mapdf21,mapif21,1,mapdm1,
     &            mapim1,ssc,possm10,rc)
c
c1.2  add t1n(a,i)aa <- M1(a,i)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,-1.0d0,mapdm1,1,mapdt13,mapit13,1,rc)
      end if
c
c
c
c2    t1n(a,i)bb <- - sum(m-b) [T1o(a,m)bb  . FII(m,i)bb]
c
cpar
      if (myRank.eq.idaabb) then
c2.1  mult M1(a,i) <= T1o(a,m)bb  . FII(m,i)bb
       call mult (wrk,wrksize,
     & 2,2,2,1,mapdt12,mapit12,1,mapdf22,mapif22,1,mapdm1,
     &            mapim1,ssc,possm10,rc)
c
c2.2  add t1n(a,i)bb <- M1(a,i)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,-1.0d0,mapdm1,1,mapdt14,mapit14,1,rc)
      end if
c
       return
       end
c
c     -----------------------------------
c
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
c
c     -----------------------------------
c
       subroutine contt29 (wrk,wrksize)
c
c     this routine do:
c
c     T29
c1    Q(c,d,ij)aaaa   <= - sum(m-a)  [ T1o(c,m)aa . <ij||md>aaaa ]
c     T2n(ab,ij)aaaa   <- Q(a,b,ij)aaaa - Q(b,a,ij)aaaa
c2    Q(c,d,ij)bbbb   <= - sum(m-b)  [ T1o(c,m)bb . <ij||md>bbbb ]
c     T2n(ab,ij)bbbb   <- Q(a,b,ij)bbbb - Q(b,a,ij)bbbb
c3    T2n(a,b,i,j)abab <- - sum(m-a)  [ T1o(a,m)aa . <ij||mb)abab ]
c4    T2n(a,b,i,j)abab <- + sum(m-b)  [ T1o(b,m)bb . <ij||ma>abba ]
c
c     N.B. use and destroy : V1,V2
c
       use Para_Info, only: MyRank
#include "ccsd2.fh"
#include "parallel.fh"
#include "wrk.fh"
c
c     help variables
c
       integer posst,rc,ssc
c
cpar
        if (myRank.eq.idfin) then
c
c1    T2n(ab,ij)aaaa <- -P(a,b) sum(m-a) [T1o(a,m)aa . <md||ij>aaaa]
c
c1.1  mult V1(c,d,ij) <= T1o(c,m)aa . <md||ij>aaaa
       call mult (wrk,wrksize,
     & 2,4,4,1,mapdt11,mapit11,1,mapdw11,mapiw11,1,mapdv1,
     &            mapiv1,ssc,possv10,rc)
c
c1.2  pack V2(cd,ij) <= V1(c,d,ij) - V1(d,c,ij)
       call fack (wrk,wrksize,
     & 4,4,mapdv1,1,mapiv1,mapdv2,mapiv2,possv20,rc)
c
c1.3  add t2n(ab,ij)aaaa <- - V2(ab,ij)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,-1.0d0,mapdv2,1,mapdt21,mapit21,1,rc)
c
c
c
c2    T2n(ab,ij)bbbb <- -P(a,b) sum(m-b) [T1o(a,m)bb . <md||ij>bbbb]
c
c2.1  mult V1(c,d,ij) <= T1o(c,m)bb . <md||ij>bbbb
       call mult (wrk,wrksize,
     & 2,4,4,1,mapdt12,mapit12,1,mapdw12,mapiw12,1,mapdv1,
     &            mapiv1,ssc,possv10,rc)
c
c2.2  pack V2(cd,ij) <= V1(c,d,ij) - V1(d,c,ij)
       call fack (wrk,wrksize,
     & 4,4,mapdv1,1,mapiv1,mapdv2,mapiv2,possv20,rc)
c
c2.3  add t2n(ab,ij)bbbb <- - V2(ab,ij)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,-1.0d0,mapdv2,1,mapdt22,mapit22,1,rc)
c
c
c
c3    T2n(a,b,i,j)abab <- - sum(m-a)  [ T1o(a,m)aa . <ij||mb)abab ]
c
c3.1  mult V1(a,b,i,j) <= T1o(a,m)aa . <mb||ij>abab
       call mult (wrk,wrksize,
     & 2,4,4,1,mapdt11,mapit11,1,mapdw13,mapiw13,1,mapdv1,
     &            mapiv1,ssc,possv10,rc)
c
c3.2  add t2n(a,b,i,j)abab <- - V1(a,b,i,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,-1.0d0,mapdv1,1,mapdt23,mapit23,1,rc)
c
c
c
c4    T2n(a,b,i,j)abab <- + sum(m-b)  [ T1o(b,m)bb . <ij||ma>abba ]
c
c4.1  mult V1(b,a,i,j) <= T1o(b,m)bb . <ma||ij>baab
       call mult (wrk,wrksize,
     & 2,4,4,1,mapdt12,mapit12,1,mapdw14,mapiw14,1,mapdv1,
     &            mapiv1,ssc,possv10,rc)
c
c4.2  map V2(a,b,i,j) <= V1(b,a,i,j)
       call map (wrk,wrksize,
     & 4,2,1,3,4,mapdv1,mapiv1,1,mapdv2,mapiv2,possv20,posst,
     &           rc)
c
c4.3  add t2n(a,b,i,j)abab <-  V2(a,b,i,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv2,1,mapdt23,mapit23,1,rc)
c
cpar
        end if
c
       return
       end
c
c     -----------------------------------
c
       subroutine contf4 (wrk,wrksize,
     & lunt2o1,lunt2o2,lunt2o3)
c
c     this routine do FIV1,FIV2,T22
c
c     FIV1
c1    FIV(b,e)aa <= FI(b,e)aa
c2    FIV(b,e)bb <= FI(b,e)bb
c
c     FIV2
c3    FIV(b,e)aa <- -0.5 sum(m-a) [ T1o(b,m)aa . FIII(e,m)aa ]
c4    FIV(b,e)bb <- -0.5 sum(m-b) [ T1o(b,m)bb . FIII(e,m)bb ]
c
c     T22
c5    Q(c,d,ij)aaaa   <- - sum(e-a) [ FIV(c,e)aa . T2o(e,d,ij)aaaa ]
c5    T2n(ab,ij)aaaa   <- Q(b,a,ij)aaaa - Q(a,b,ij)aaaa
c6    Q(c,d,ij)bbbb   <- - sum(e-b) [ FIV(c,e)bb . T2o(e,d,ij)bbbb ]
c6    T2n(ab,ij)bbbb   <- Q(b,a,ij)bbbb - Q(a,b,ij)bbbb
c7    T2n(a,b,i,j)abab <- sum (e-a)  [ FIV(a,e)aa . T2o(e,b,i,j)abab ]
c8    T2n(a,b,i,j)abab <- sum (e-b)  [ FIV(b,e)bb . T2o(a,e,i,j)abab ]
c
c     N.B. use and destroy : V1,V2,V3,M1,M2,M3
c     N.B. # of read       : 3
c     N.B. Parallel : in the case where idaaaa.ne.idbaab and
c          idaaaa.ne.idbaab this routine runs contributions to
c          T2n also on idaaaa and idbbbb nodes, but only with
c          FIV = FI, since on these nodes there is a specific
c          part of contributions F13 (see notes in sumoverb routine)
c           which are not presented on pilot nodes
c
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
       if ((myRank.eq.idbaab).or.(myRank.eq.idaabb)
     & .or.(myRank.eq.idaaaa).or.(myRank.eq.idbbbb)) then
c78.0 read V1(c,d,i,j) <= T2o(c,d,i,j)abab
       call filemanager (2,lunt2o3,rc)
       call getmediate (wrk,wrksize,
     & lunt2o3,possv10,mapdv1,mapiv1,rc)
       end if
c
c
c
c
c1    FIV(b,e)aa <= FI(b,e)aa
c
cpar
       if ((myRank.eq.idbaab).or.(myRank.eq.idaaaa)) then
c1.1  map M1(b,e) <= F1(b,e)aa
       call map (wrk,wrksize,
     & 2,1,2,0,0,mapdf11,mapif11,1,mapdm1,mapim1,possm10,
     &           posst,rc)
       end if
c
c
c
c3    FIV(b,e)aa <- -0.5 sum(m-a) [ T1o(b,m)aa . FIII(e,m)aa ]
c
cpar
       if (myRank.eq.idbaab) then
c
c3.1  map M2(m,e) <= f3(e,m)aa
       call map (wrk,wrksize,
     & 2,2,1,0,0,mapdf31,mapif31,1,mapdm2,mapim2,possm20,
     &           posst,rc)
c
c3.2  mult M3(b,e) <= T1o(b,m)aa . M2(m,e)
       call mult (wrk,wrksize,
     & 2,2,2,1,mapdt11,mapit11,1,mapdm2,mapim2,1,mapdm3,
     &            mapim3,ssc,possm30,rc)
c
c3.3  add M1(b,e) <- -0.5 M3(b,e)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,-0.5d0,mapdm3,1,mapdm1,mapim1,1,rc)
cparend
        end if
c
c
cpar
       if ((myRank.eq.idbaab).or.(myRank.eq.idaaaa)) then
c
c5    Q(c,d,ij)aaaa   <- - sum(e-a) [ FIV(c,e)aa . T2o(e,d,ij)aaaa ]
c5    T2n(ab,ij)aaaa   <- Q(b,a,ij)aaaa - Q(a,b,ij)aaaa
c
c5.1  read V2(ed,ij) <= T2o(ed,ij)aaaa
       call filemanager (2,lunt2o1,rc)
       call getmediate (wrk,wrksize,
     & lunt2o1,possv20,mapdv2,mapiv2,rc)
c
c5.2  expand V3(e,d,ij) <= V2 (ed,ij)
       call expand (wrk,wrksize,
     & 4,5,mapdv2,mapiv2,1,possv30,mapdv3,mapiv3,rc)
c
c5.3  mult V2(c,d,ij) <= M1(c,e) . V3(e,d,ij)
       call mult (wrk,wrksize,
     & 2,4,4,1,mapdm1,mapim1,1,mapdv3,mapiv3,1,mapdv2,mapiv2,
     &            ssc,possv20,rc)
c
c5.4  pack V3(ab,ij) <= V2(a,b,ij) - V2(b,a,ij)
       call fack (wrk,wrksize,
     & 4,4,mapdv2,1,mapiv2,mapdv3,mapiv3,possv30,rc)
c
c5.5  add T2n(ab,ij)aaaa <- V3(ab,ij)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv3,1,mapdt21,mapit21,1,rc)
c
c
c
c7    T2n(a,b,i,j)abab <- sum (e-a)  [ FIV(a,e)aa . T2o(e,b,i,j)abab ]
c
c7.1  mult V2(a,b,i,j) <= M1(a,e) . V1(e,b,i,j)
       call mult (wrk,wrksize,
     & 2,4,4,1,mapdm1,mapim1,1,mapdv1,mapiv1,1,mapdv2,mapiv2,
     &            ssc,possv20,rc)
c
c7.2  add T2n(a,b,i,j)abab <- V2(a,b,i,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv2,1,mapdt23,mapit23,1,rc)
cparend
        end if
c
c
c
c2    FIV(b,e)bb <= FI(b,e)bb
c
cpar
       if ((myRank.eq.idaabb).or.(myRank.eq.idbbbb)) then
c2.1  map M1(b,e) <= F1(b,e)bb
       call map (wrk,wrksize,
     & 2,1,2,0,0,mapdf12,mapif12,1,mapdm1,mapim1,possm10,
     &           posst,rc)
         end if
c
c
c
c4    FIV(b,e)bb <- -0.5 sum(m-b) [ T1o(b,m)bb . FIII(e,m)bb ]
c
cpar
       if (myRank.eq.idaabb) then
c4.1  map M2(m,e) <= f3(e,m)bb
       call map (wrk,wrksize,
     & 2,2,1,0,0,mapdf32,mapif32,1,mapdm2,mapim2,possm20,
     &           posst,rc)
c
c4.2  mult M3(b,e) <= T1o(b,m)bb . M2(m,e)
       call mult (wrk,wrksize,
     & 2,2,2,1,mapdt12,mapit12,1,mapdm2,mapim2,1,mapdm3,
     &            mapim3,ssc,possm30,rc)
c
c4.3  add M1(b,e) <- -0.5 M3(b,e)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,-0.5d0,mapdm3,1,mapdm1,mapim1,1,rc)
cparend
        end if
c
c
c
cpar
       if ((myRank.eq.idaabb).or.(myRank.eq.idbbbb)) then
c
c6    Q(c,d,ij)bbbb   <- - sum(e-b) [ FIV(c,e)bb . T2o(e,d,ij)bbbb ]
c6    T2n(ab,ij)bbbb   <- Q(b,a,ij)bbbb - Q(a,b,ij)bbbb
c
c6.1  read V2(cd,ij) <= T2o(cd,ij)bbbb
       call filemanager (2,lunt2o2,rc)
       call getmediate (wrk,wrksize,
     & lunt2o2,possv20,mapdv2,mapiv2,rc)
c
c6.2  expand V3(e,d,ij) <= V2 (cd,ij)
       call expand (wrk,wrksize,
     & 4,5,mapdv2,mapiv2,1,possv30,mapdv3,mapiv3,rc)
c
c6.3  mult V2(c,d,ij) <= M1(c,d) . V3(e,d,ij)
       call mult (wrk,wrksize,
     & 2,4,4,1,mapdm1,mapim1,1,mapdv3,mapiv3,1,mapdv2,mapiv2,
     &            ssc,possv20,rc)
c
c6.4  pack V3(ab,ij) <= V2(a,b,ij) - V2(b,a,ij)
       call fack (wrk,wrksize,
     & 4,4,mapdv2,1,mapiv2,mapdv3,mapiv3,possv30,rc)
c
c6.5  add T2n(ab,ij)bbbb <- V3(ab,ij)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv3,1,mapdt22,mapit22,1,rc)
c
c
c
c8    T2n(a,b,i,j)abab <- sum (e-b)  [ FIV(b,e)bb . T2o(a,e,i,j)abab ]
c
c8.1  map V3(e,a,i,j) <= V1(a,e,i,j)
       call map (wrk,wrksize,
     & 4,2,1,3,4,mapdv1,mapiv1,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
c
c8.2  mult V2(b,a,i,j) <= M1(b,e) . V3(e,a,i,j)
       call mult (wrk,wrksize,
     & 2,4,4,1,mapdm1,mapim1,1,mapdv3,mapiv3,1,mapdv2,mapiv2,
     &            ssc,possv20,rc)
c
c8.3  map V3(a,b,i,j) <= V2(b,a,i,j)
       call map (wrk,wrksize,
     & 4,2,1,3,4,mapdv2,mapiv2,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
c
c8.4  add T2n(a,b,i,j)abab <- V3(a,b,i,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv3,1,mapdt23,mapit23,1,rc)
cparend
        end if
c
       return
       end
c
c     -----------------------------------
c
       subroutine contf5 (wrk,wrksize,
     & lunt2o1,lunt2o2,lunt2o3)
c
c     this routine do : FV1,FV2,T23
c
c     FV1
c1    FV(m,j)aa  <= FII(m,j)aa
c2    FV(m,j)bb  <= FII(m,j)bb
c
c     FV2
c3    FV(m,j)aa  <- -0.5 sum (e-a) [ FIII(e,m)aa . T1o(e,j)aa ]
c4    FV(m,j)bb  <- -0.5 sum (e-b) [ FIII(e,m)bb . T1o(e,j)bb ]
c
c     T23
c5    Q(ab,k,l)aaaa   <= sum(m-a)   [ T2o(ab,k,m)aaaa . FV(m,l)aa ]
c5    T2n(ab,ij)aaaa   <- Q(ab,j,i)aaaa - Q(ab,i,j)aaaa
c6    Q(ab,k,l)bbbb   <= sum(m-b)   [ T2o(ab,k,m)bbbb . FV(m,l)bb ]
c6    T2n(ab,ij)bbbb   <- Q(ab,j,i)bbbb - Q(ab,i,j)bbbb
c7    T2n(a,b,i,j)abab <- - sum(m-b) [ T2o(a,b,i,m)abab . FV(m,j)bb ]
c8    T2n(a,b,i,j)abab <- - sum(m-a) [ T2o(a,b,m,j)abab . FV(m,i)aa ]
c
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
c
cpar
       if ((myRank.eq.idbaab).or.(myRank.eq.idaabb)) then
c78.0 read V1(a,b,k,l) <= T2o(a,b,k,l)abab
       call filemanager (2,lunt2o3,rc)
       call getmediate (wrk,wrksize,
     & lunt2o3,possv10,mapdv1,mapiv1,rc)
       end if
c
c
cpar
       if (myRank.eq.idbaab) then
c
c1    FV(m,j)aa  <= FII(m,j)aa
c
c1.1  map M1(m,j) <= F2(m,j)aa
       call map (wrk,wrksize,
     & 2,1,2,0,0,mapdf21,mapif21,1,mapdm1,mapim1,possm10,
     &           posst,rc)
c
c
c
c3    FV(m,j)aa  <- 0.5 sum (e-a) [ FIII(e,m)aa . T1o(e,j)aa ]
c
c3.1  map M2(m,e) <= F3(e,m)aa
       call map (wrk,wrksize,
     & 2,2,1,0,0,mapdf31,mapif31,1,mapdm2,mapim2,possm20,
     &           posst,rc)
c
c3.2  mult M3(m,j) <= M2(m,e) . T1o(e,j)aa
       call mult (wrk,wrksize,
     & 2,2,2,1,mapdm2,mapim2,1,mapdt11,mapit11,1,mapdm3,
     &            mapim3,ssc,possm30,rc)
c
c3.3  add M1(m,j) <- 0.5 M3(m,j)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,0.5d0,mapdm3,1,mapdm1,mapim1,1,rc)
c
c
c
c5    Q(ab,k,l)aaaa   <= sum(m-a)   [ T2o(ab,k,m)aaaa . FV(m,l)aa ]
c5    T2n(ab,ij)aaaa   <- Q(ab,j,i)aaaa - Q(ab,i,j)aaaa
c
c5.1  read V3(ab,kl) <= T2o(ab,kl)aaaa
       call filemanager (2,lunt2o1,rc)
       call getmediate (wrk,wrksize,
     & lunt2o1,possv30,mapdv3,mapiv3,rc)
c
c5.2  expand V2(ab,k,l) <- V3(ab,kl)
       call expand (wrk,wrksize,
     & 4,6,mapdv3,mapiv3,1,possv20,mapdv2,mapiv2,rc)
c
c5.3  mult V3(ab,k,l) <= V2(a,b,k,m) . M1(m,l)
       call mult (wrk,wrksize,
     & 4,2,4,1,mapdv2,mapiv2,1,mapdm1,mapim1,1,mapdv3,mapiv3,
     &            ssc,possv30,rc)
c
c5.4  pack V2(ab,ij) <= V3(ab,i,j) - V3(ab,j,i)
       call fack (wrk,wrksize,
     & 4,4,mapdv3,1,mapiv3,mapdv2,mapiv2,possv20,rc)
c
c5.5  add T2n(ab,ij)aaaa <- - V2(ab,ij)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,-1.0d0,mapdv2,1,mapdt21,mapit21,1,rc)
c
c
c
c8    T2n(a,b,i,j)abab <- - sum(m-a) [ T2o(a,b,m,j)abab . FV(m,i)aa ]
c
c8.1  map V2(a,b,j,m) <= V1(a,b,m,j)
       call map (wrk,wrksize,
     & 4,1,2,4,3,mapdv1,mapiv1,1,mapdv2,mapiv2,possv20,posst,
     &           rc)
c
c8.2  mult V3(a,b,j,i) <= V2(a,b,j,m) . M1(m,i)
       call mult (wrk,wrksize,
     & 4,2,4,1,mapdv2,mapiv2,1,mapdm1,mapim1,1,mapdv3,mapiv3,
     &            ssc,possv30,rc)
c
c8.3  map V2(a,b,i,j) <= V3(a,b,j,i)
       call map (wrk,wrksize,
     & 4,1,2,4,3,mapdv3,mapiv3,1,mapdv2,mapiv2,possv20,posst,
     &           rc)
c
c8.4  add T2n(a,b,i,j)abab <- - V2(a,b,i,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,-1.0d0,mapdv2,1,mapdt23,mapit23,1,rc)
cpar
        end if
c
c
cpar
       if (myRank.eq.idaabb) then
c
c2    FV(m,j)bb  <= FII(m,j)bb
c
c2.1  map M1(m,j) <= F2(m,j)bb
       call map (wrk,wrksize,
     & 2,1,2,0,0,mapdf22,mapif22,1,mapdm1,mapim1,possm10,
     &           posst,rc)
c
c
c
c4    FV(m,j)bb  <- 0.5 sum (e-b) [ FIII(e,m)bb . T1o(e,j)bb ]
c
c4.1  map M2(m,e) <= F3(e,m)bb
       call map (wrk,wrksize,
     & 2,2,1,0,0,mapdf32,mapif32,1,mapdm2,mapim2,possm20,
     &           posst,rc)
c
c4.2  mult M3(m,j) <= M2(m,e) . T1o(e,j)bb
       call mult (wrk,wrksize,
     & 2,2,2,1,mapdm2,mapim2,1,mapdt12,mapit12,1,mapdm3,
     &            mapim3,ssc,possm30,rc)
c
c4.3  add M1(m,j) <- 0.5 M3(m,j)
       call add (wrk,wrksize,
     & 2,2,0,0,0,0,1,1,0.5d0,mapdm3,1,mapdm1,mapim1,1,rc)
c
c
c
c6    Q(ab,k,l)bbbb   <= sum(m-b)   [ T2o(ab,k,m)bbbb . FV(m,l)bb ]
c6    T2n(ab,ij)bbbb   <- Q(ab,j,i)bbbb - Q(ab,i,j)bbbb
c
c6.1  read V3(ab,kl) <= T2o(ab,kl)bbbb
       call filemanager (2,lunt2o2,rc)
       call getmediate (wrk,wrksize,
     & lunt2o2,possv30,mapdv3,mapiv3,rc)
c
c6.2  expand V2(ab,k,l) <- V3(ab,kl)
       call expand (wrk,wrksize,
     & 4,6,mapdv3,mapiv3,1,possv20,mapdv2,mapiv2,rc)
c
c6.3  mult V3(ab,k,l) <= V2(a,b,k,m) . M1(m,l)
       call mult (wrk,wrksize,
     & 4,2,4,1,mapdv2,mapiv2,1,mapdm1,mapim1,1,mapdv3,mapiv3,
     &            ssc,possv30,rc)
c
c6.4  pack V2(ab,ij) <= V3(ab,i,j) - V3(ab,j,i)
       call fack (wrk,wrksize,
     & 4,4,mapdv3,1,mapiv3,mapdv2,mapiv2,possv20,rc)
c
c6.5  add T2n(ab,ij)bbbb <- - V2(ab,ij)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,-1.0d0,mapdv2,1,mapdt22,mapit22,1,rc)
c
c
c
c7    T2n(a,b,i,j)abab <- - sum(m-b) [ T2o(a,b,i,m)abab . FV(m,j)bb ]
c
c7.1  mult V3(a,b,i,j) <= V1(a,b,i,m) . M1(m,j)
       call mult (wrk,wrksize,
     & 4,2,4,1,mapdv1,mapiv1,1,mapdm1,mapim1,1,mapdv3,mapiv3,
     &            ssc,possv30,rc)
c
c7.2  add T2n(a,b,i,j)abab <- - V3(a,b,i,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,-1.0d0,mapdv3,1,mapdt23,mapit23,1,rc)
cpar
        end if
c
       return
       end

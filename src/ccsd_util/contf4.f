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

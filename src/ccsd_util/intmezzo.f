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
c     intermezzo
c     getw3
c
c     ---------------------------------------------------------
c
       subroutine intermezzo (wrk,wrksize,
     & lunw3aaaa,lunw3bbbb,lunw3abba,
     & lunw3baab,lunw3aabb,lunw3bbaa,lunt2o1,lunt2o2,lunt2o3,
     & lunabij1,lunabij2,lunabij3)
c
c     this routine calculate contributions
c     WIII3, WIII4 and T26
c
c     assigment of spin combinations:
c
c     WIII(m,e,b,j)aaaa  - I
c     WIII(m,e,b,j)bbbb  - K
c     WIII(m,e,b,j)aabb  - L
c     WIII(m,e,b,j)abba  - M
c     WIII(m,e,b,j)baab  - N
c     WIII(m,e,b,j)bbaa  - J
c
c     WIII3
c     WIII(m,e,b,j)aaaa <- sum(n-a) [ <mn||je>aaaa . T1o(n,b)aa ]
c     WIII(m,e,b,j)bbbb <- sum(n-b) [ <mn||je>bbbb . T1o(n,b)bb ]
c     WIII(m,e,b,j)aabb <- sum(n-b) [ <mn||je>abba . T1o(n,b)bb ]
c     WIII(m,e,b,j)abba <- sum(n-b) [ <mn||je>abab . T1o(n,b)bb ]
c     WIII(m,e,b,j)baab <- -sum(n-a) [ <nm||je>abba . T1o(n,b)aa ]
c     WIII(m,e,b,j)bbaa <- -sum(n-a) [ <je||nm>abab . T1o(n,b)aa ]
c
c     WIII4
c     Q(f,b,j,n)aaaa   <= 0.5 T2o(f,b,j,n)aaaa + T1o(f,j)aa . T1o(b,n)aa
c     WIII(m,e,b,j)aaaa <- -sum(n,f-aa)     [ Q(f,b,j,n)aaaa   . <ef||mn>aaaa ]
c     <- 0.5 sum(n,f-bb)  [ T2o(b,f,j,n)abab . <ef||mn>abab ]
c     WIII(m,e,b,j)aabb <- 0.5 sum(n,f-aa)  [ T2o(f,b,n,j)abab . <ef||mn>aaaa ]
c     <- -sum(n,f-bb)     [ Q(f,b,j,n)bbbb   . <ef||mn>abab ]
c     Q(f,b,j,n)bbbb   <= 0.5 T2o(f,b,j,n)bbbb + T1o(f,j)nn . T1o(b,n)bb
c     WIII(m,e,b,j)bbbb <- -sum(n,f-bb)     [ Q(f,b,j,n)bbbb   . <ef||mn>bbbb ]
c     <- 0.5 sum(n,f-aa)  [ T2o(f,b,n,j)abab . <fe||nm>abab ]
c     WIII(m,e,b,j)bbaa <- 0.5 sum(n,f-bb)  [ T2o(b,f,j,n)abab . <ef||mn>bbbb ]
c     <- - sum(n,f-aa)    [ Q(f,b,j,n)aaaa   . <fe||nm>abab ]
c     Q(f,b,j,n)abab   <= 0.5 T2o(f,b,j,n)abab + T1o(f,j)aa . T1o(b,n)bb
c     WIII(m,e,b,j)abba <- sum(n,f-ba)      [ Q(f,b,j,n)abab   . <fe||mn>abab ]
c     Q(b,f,n,j)abab   <= 0.5 T2o(b,f,n,j)abab + T1o(b,n)aa . T1o(fmj)bb
c     WIII(m,e,b,j)baab <- sum(n,f-ab)      [ Q(b,f,n,j)abab   . <ef||nm>abab ]
c
c
c     T26
c     R1(a,i,b,j)aaaa <= sum(m,e-aa) [ T2o(a,e,i,m)aaaa . WIII(m,e,b,j)aaaa ]
c     <- sum(m,e-bb) [ T2o(a,e,i,m)abab . WIII(m,e,b,j)bbaa ]
c     T2n(ab,ij)aaaa   <= {1(a,i,b,j)-R1(b,i,a,j)-R1(a,j,b,i)+R1(b,j,a,i)}aaaa
c     R1(a,i,b,j)bbbb <= sum(m,e-bb) [ T2o(a,e,i,m)bbbb . WIII(m,e,b,j)bbbb ]
c     <- sum(m,e-aa) [ T2o(e,a,m,i)abab . WIII(m,e,b,j)aabb ]
c     T2n(ab,ij)bbbb   <= {1(a,i,b,j)-R1(b,i,a,j)-R1(a,j,b,i)+R1(b,j,a,i)}bbbb
c     T2n(a,b,i,j)abab <- sum(m,e-aa) [ T2o(a,e,i,m)aaaa . WIII(m,e,b,j)aabb ]
c     <- sum(m,e-aa) [ T2o(e,b,m,j)abab . WIII(m,e,a,i)aaaa ]
c     <- sum(m,e-bb) [ T2o(a,e,i,m)abab . WIII(m,e,b,j)bbbb ]
c     <- sum(m,e-bb) [ T2o(b,e,j,m)bbbb . WIII(m,e,a,i)bbaa ]
c     <- sum(m,e-ab) [ T2o(a,e,m,j)abab . WIII(m,e,b,i)abba ]
c     <- sum(m,e-ab) [ T2o(e,b,i,m)abab . WIII(m,e,a,j)baab ]
c
c     N.B. use and destry : V1,V2,V3,V4,M1,M2
c     N.B. # of read      : 30 + 6
c     # of write     : 2
c
        use Para_Info, only: MyRank
        implicit none
#include "ccsd2.fh"
#include "wrk.fh"
#include "parallel.fh"
c
       integer lunw3aaaa,lunw3bbbb,lunw3abba,lunw3baab,lunw3aabb,
     & lunw3bbaa
       integer lunt2o1,lunt2o2,lunt2o3,lunabij1,lunabij2,lunabij3
c
c     help variableas
c
       integer posst,rc,ssc,lunqaaaa,lunqbbbb
c
c
cA.1  rewind nw3 files
cpar
        if (myRank.eq.idaaaa) then
       call filemanager (2,lunw3aaaa,rc)
        end if
c
        if (myRank.eq.idbaab) then
       call filemanager (2,lunw3baab,rc)
        end if
c
        if (myRank.eq.idbbaa) then
       call filemanager (2,lunw3bbaa,rc)
        end if
c
        if (myRank.eq.idbbbb) then
       call filemanager (2,lunw3bbbb,rc)
        end if
c
        if (myRank.eq.idabba) then
       call filemanager (2,lunw3abba,rc)
        end if
c
        if (myRank.eq.idaabb) then
       call filemanager (2,lunw3aabb,rc)
        end if
c
c0.1  map M1(i,a) <- T1o(a,i)aa
       call map (wrk,wrksize,
     & 2,2,1,0,0,mapdt11,mapit11,1,mapdm1,mapim1,possm10,
     &           posst,rc)
c0.2  map M2(i,a) <- T1o(a,i)bb
       call map (wrk,wrksize,
     & 2,2,1,0,0,mapdt12,mapit12,1,mapdm2,mapim2,possm20,
     &           posst,rc)
c
c
c     part I - W3aaaa
cpar  all contributions are calculated for idaaaa only,
c     just Qaaaa is calculated both on idaaaa and idbbaa
c
c     par
      if (idaaaa.eq.myRank) then
c
cI.1  get V1(m,e,a,j) <- W3aaaa(m,e,a,j)
       call getw3 (wrk,wrksize,
     & lunw3aaaa,1)
c
cI.2  WIII(m,e,b,j)aaaa <- sum(n-a) [ <mn||je>aaaa . T1o(n,b)aa ]
c
cI.2.1expand V2(j,e,m,n) <- <je||mn>aaaa
       call expand (wrk,wrksize,
     & 4,3,mapdw11,mapiw11,1,possv20,mapdv2,mapiv2,rc)
cI.2.2map V3(m,e,j,n) <- V2(j,e,m,n)
       call map (wrk,wrksize,
     & 4,3,2,1,4,mapdv2,mapiv2,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
cI.2.3mult V2(m,e,j,b) <- V3(m,e,j,n) . M1(n,b)
       call mult (wrk,wrksize,
     & 4,2,4,1,mapdv3,mapiv3,1,mapdm1,mapim1,1,mapdv2,mapiv2,
     &            ssc,possv20,rc)
cI.2.4map V3(m,e,b,j) <- V2(m,e,j,b)
       call map (wrk,wrksize,
     & 4,1,2,4,3,mapdv2,mapiv2,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
cI.2.5add V1(m,e,b,j) <- 1.0d0 . V3(m,e,b,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv3,1,mapdv1,mapiv1,1,rc)
cparend
        end if
c
cI.3  WIII4
c     Q(f,b,j,n)aaaa   <= 0.5 T2o(f,b,j,n)aaaa + T1o(f,j)aa . T1o(b,n)aa
c     WIII(m,e,b,j)aaaa <- -sum(n,f-aa)     [ Q(f,b,j,n)aaaa   . <ef||mn>aaaa ]
c     <- 0.5 sum(n,f-bb)  [ T2o(b,f,j,n)abab . <ef||mn>abab ]
c
cpar
        if ((myRank.eq.idaaaa).or.(myRank.eq.idbbaa)) then
cI.3.1get V3(fb,jn) <- T2o(fb,jn)aaaa
       call filemanager (2,lunt2o1,rc)
       call getmediate (wrk,wrksize,
     & lunt2o1,possv30,mapdv3,mapiv3,rc)
cI.3.2expand V2(f,b,j,n) <- V3(fb,jn)
       call expand (wrk,wrksize,
     & 4,4,mapdv3,mapiv3,1,possv20,mapdv2,mapiv2,rc)
cI.3.3mkQ V2(f,b,j,n) <- 0.5 V2(f,b,j,n) + T1o(f,j)aa . T1o(b,n)aa
       call mkq(wrk,wrksize,
     & mapdv2,mapiv2,mapdt11,mapit11,mapdt11,mapit11,0.5d0,rc)
cI.3.4map V3(f,n,b,j) <- V2(f,b,j,n)
       call map (wrk,wrksize,
     & 4,1,3,4,2,mapdv2,mapiv2,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
cparend
        end if
cpar
        if (myRank.eq.idbbaa) then
cI.3.5write V3(f,n,b,j) to lunqaaaa
       call filemanager (1,lunqaaaa,rc)
       call wrtmediate (wrk,wrksize,
     & lunqaaaa,mapdv3,mapiv3,rc)
        end if
c
cpar
        if (myRank.eq.idaaaa) then
cI.3.6get V2(ef,mn) = <ef||mn>aaaa from luna file
       call filemanager (2,lunabij1,rc)
       call getmediate (wrk,wrksize,
     & lunabij1,possv20,mapdv2,mapiv2,rc)
cI.3.7expand V4(e,f,m,n) <- V2(ef,mn)
       call expand (wrk,wrksize,
     & 4,4,mapdv2,mapiv2,1,possv40,mapdv4,mapiv4,rc)
cI.3.8map V2(m,e,f,n) <- V4(e,f,m,n)
       call map (wrk,wrksize,
     & 4,2,3,1,4,mapdv4,mapiv4,1,mapdv2,mapiv2,possv20,posst,
     &           rc)
cI.3.9mult V4(m,e,b,j) <- V2(m,e,f,n) . V3(f,n,b,j)
       call mult (wrk,wrksize,
     & 4,4,4,2,mapdv2,mapiv2,1,mapdv3,mapiv3,1,mapdv4,mapiv4,
     &            ssc,possv40,rc)
cI.3.10 add V1(m,e,b,j) <- -1.0d0 . V4(m,e,b,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,-1.0d0,mapdv4,1,mapdv1,mapiv1,1,rc)
cI.3.11 get V2(e,f,m,n) <- <ef||mn>abab
       call filemanager (2,lunabij3,rc)
       call getmediate (wrk,wrksize,
     & lunabij3,possv20,mapdv2,mapiv2,rc)
cI.3.12 map V3(m,e,f,n) <- V2(e,f,m,n)
       call map (wrk,wrksize,
     & 4,2,3,1,4,mapdv2,mapiv2,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
cI.3.13 get V2(b,f,j,n) <- T2o(b,f,j,n)abab
       call filemanager (2,lunt2o3,rc)
       call getmediate (wrk,wrksize,
     & lunt2o3,possv20,mapdv2,mapiv2,rc)
cI.3.14 map V4(f,n,b,j) <- V2(b,f,j,n)
       call map (wrk,wrksize,
     & 4,3,1,4,2,mapdv2,mapiv2,1,mapdv4,mapiv4,possv40,posst,
     &           rc)
cI.3.15 mult V2(m,e,b,j) <- V3(m,e,f,n) . V4(f,n,b,j)
       call mult (wrk,wrksize,
     & 4,4,4,2,mapdv3,mapiv3,1,mapdv4,mapiv4,1,mapdv2,mapiv2,
     &            ssc,possv20,rc)
cI.3.16 add V1(m,e,b,j) <- 0.5d0 . V2(m,e,b,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,0.5d0,mapdv2,1,mapdv1,mapiv1,1,rc)
c
cI.4  R1(a,i,b,j)aaaa <= sum(m,e-aa) [ T2o(a,e,i,m)aaaa . WIII(m,e,b,j)aaaa ]
c     T2n(ab,ij)aaaa   <= {1(a,i,b,j)-R1(b,i,a,j)-R1(a,j,b,i)+R1(b,j,a,i)}aaaa
cI.4.1get V2(ae,im) <- T2o(ae,im)aaaa
       call filemanager (2,lunt2o1,rc)
       call getmediate (wrk,wrksize,
     & lunt2o1,possv20,mapdv2,mapiv2,rc)
cI.4.2expand V3(a,e,i,m) <- V2(ae,im)
       call expand (wrk,wrksize,
     & 4,4,mapdv2,mapiv2,1,possv30,mapdv3,mapiv3,rc)
cI.4.3map V2(a,i,m,e) <- V3(a,e,i,m)
       call map (wrk,wrksize,
     & 4,1,4,2,3,mapdv3,mapiv3,1,mapdv2,mapiv2,possv20,posst,
     &           rc)
cI.4.4mult V4(a,i,b,j) <- V2(a,i,m,e) . V1(m,e,b,j)
       call mult (wrk,wrksize,
     & 4,4,4,2,mapdv2,mapiv2,1,mapdv1,mapiv1,1,mapdv4,mapiv4,
     &            ssc,possv40,rc)
cI.4.5map V3(a,b,i,j) <- V4(a,i,b,j)
       call map (wrk,wrksize,
     & 4,1,3,2,4,mapdv4,mapiv4,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
cI.4.6pack V3(ab,ij) <- V2(ab,i,j) <- V3(a,b,i,j)
       call fack (wrk,wrksize,
     & 4,1,mapdv3,1,mapiv3,mapdv2,mapiv2,possv20,rc)
       call fack (wrk,wrksize,
     & 4,4,mapdv2,1,mapiv2,mapdv3,mapiv3,possv30,rc)
cI.4.7add T2n(ab,ij)aaaa <- 1.0d0 V3(ab,ij)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv3,1,mapdt21,mapit21,1,rc)
c
cI.5  T2n(a,b,i,j)abab <- sum(m,e-aa) [ T2o(e,b,m,j)abab . WIII(m,e,a,i)aaaa ]
cI.5.1get V4(e,b,m,j) <- T2o(e,b,m,j)abab
       call filemanager (2,lunt2o3,rc)
       call getmediate (wrk,wrksize,
     & lunt2o3,possv40,mapdv4,mapiv4,rc)
cI.5.2map V2(b,j,m,e) <- V4(e,b,m,j)
       call map (wrk,wrksize,
     & 4,4,1,3,2,mapdv4,mapiv4,1,mapdv2,mapiv2,possv20,posst,
     &           rc)
cI.5.3mult V4(b,j,a,i) <- V2(b,j,m,e) . V1(m,e,a,i)
       call mult (wrk,wrksize,
     & 4,4,4,2,mapdv2,mapiv2,1,mapdv1,mapiv1,1,mapdv4,mapiv4,
     &            ssc,possv40,rc)
cI.5.4map V3(a,b,i,j) <- V4(b,j,a,i)
       call map (wrk,wrksize,
     & 4,2,4,1,3,mapdv4,mapiv4,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
cI.5.5add T2n(a,b,i,j)abab <- 1.0d0 . V3(a,b,i,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv3,1,mapdt23,mapit23,1,rc)
cparend
       end if
c
c
c     J part W3bbaa
c
c     par
      if (myRank.eq.idbbaa) then
c
cJ.1  get V1(m,e,a,j) <- W3bbaa(m,e,a,j)
       call getw3 (wrk,wrksize,
     & lunw3bbaa,6)
c
cJ.2  WIII(m,e,b,j)bbaa <- -sum(n-a) [ <je||nm>abab . T1o(n,b)aa ]
cJ.2.1map V3(m,e,j,n) <- <j,e||n,m>abab
       call map (wrk,wrksize,
     & 4,3,2,4,1,mapdw13,mapiw13,1,mapdv3,mapiv3,possv30,
     &           posst,rc)
cJ.2.2mult V4(m,e,j,b) <- V3(m,e,j,n) . M1(n,b)
       call mult (wrk,wrksize,
     & 4,2,4,1,mapdv3,mapiv3,1,mapdm1,mapim1,1,mapdv4,mapiv4,
     &            ssc,possv40,rc)
cJ.2.3map V3(m,e,b,j) <- V4(m,e,j,b)
       call map (wrk,wrksize,
     & 4,1,2,4,3,mapdv4,mapiv4,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
cJ.2.4add V1(m,e,b,j) <- -1.0d0 V3(m,e,b,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,-1.0d0,mapdv3,1,mapdv1,mapiv1,1,rc)
c
cJ.2  WIII(m,e,b,j)bbaa <- 0.5 sum(n,f-bb)  [ T2o(b,f,j,n)abab . <ef||mn>bbbb ]
c     <- - sum(n,f-aa)    [ Q(f,b,j,n)aaaa   . <fe||nm>abab ]
cJ.2.1get V2(f,n,b,j) from lunqaaaa (produced in I step) and close it
       call filemanager (2,lunqaaaa,rc)
       call getmediate (wrk,wrksize,
     & lunqaaaa,possv20,mapdv2,mapiv2,rc)
       call filemanager (3,lunqaaaa,rc)
cJ.2.2get V3(f,e,n,m) <- <fe||nm>abab
       call filemanager (2,lunabij3,rc)
       call getmediate (wrk,wrksize,
     & lunabij3,possv30,mapdv3,mapiv3,rc)
cJ.2.3map V4(m,e,f,n) <- V3(f,e,n,m)
       call map (wrk,wrksize,
     & 4,3,2,4,1,mapdv3,mapiv3,1,mapdv4,mapiv4,possv40,posst,
     &           rc)
cJ.2.4mult V3(m,e,b,j) <- V4(m,e,f,n) . V2(f,n,b,j)
       call mult (wrk,wrksize,
     & 4,4,4,2,mapdv4,mapiv4,1,mapdv2,mapiv2,1,mapdv3,mapiv3,
     &            ssc,possv30,rc)
cJ.2.5add V1(m,e,b,j) <- -1.0d0 . V3(m,e,b,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,-1.0d0,mapdv3,1,mapdv1,mapiv1,1,rc)
cJ.2.6get V2(b,f,j,n) <- T2o(b,f,j,n)abab
       call filemanager (2,lunt2o3,rc)
       call getmediate (wrk,wrksize,
     & lunt2o3,possv20,mapdv2,mapiv2,rc)
cJ.2.7map V3(f,n,b,j) <- V2(b,f,j,n)
       call map (wrk,wrksize,
     & 4,3,1,4,2,mapdv2,mapiv2,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
cJ.3.8get V2(ef,mn) = <ef||mn>bbbb from lunb file
       call filemanager (2,lunabij2,rc)
       call getmediate (wrk,wrksize,
     & lunabij2,possv20,mapdv2,mapiv2,rc)
cJ.3.9exp V4(e,f,m,n) <- V2(ef,mn)
       call expand (wrk,wrksize,
     & 4,4,mapdv2,mapiv2,1,possv40,mapdv4,mapiv4,rc)
cJ.3.10        map V2(m,e,f,n) <- V4(e,f,m,n)
       call map (wrk,wrksize,
     & 4,2,3,1,4,mapdv4,mapiv4,1,mapdv2,mapiv2,possv20,posst,
     &           rc)
cJ.3.11 mult V4(m,e,b,j) <- V2(m,e,f,n) . V3(f,n,b,j)
       call mult (wrk,wrksize,
     & 4,4,4,2,mapdv2,mapiv2,1,mapdv3,mapiv3,1,mapdv4,mapiv4,
     &            ssc,possv40,rc)
cJ.3.12        add V1(m,e,b,j) <- 0.5d0 . V4(m,e,b,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,0.5d0,mapdv4,1,mapdv1,mapiv1,1,rc)
c
cJ.4  R1(a,i,b,j)aaaa <= sum(m,e-bb) [ T2o(a,e,i,m)abab . WIII(m,e,b,j)bbaa ]
c     T2n(ab,ij)aaaa   <= {1(a,i,b,j)-R1(b,i,a,j)-R1(a,j,b,i)+R1(b,j,a,i)}aaaa
cJ.4.1get V3(a,e,i,m) <- T2o(a,e,i,m)abab
       call filemanager (2,lunt2o3,rc)
       call getmediate (wrk,wrksize,
     & lunt2o3,possv30,mapdv3,mapiv3,rc)
cJ.4.2map V2(a,i,m,e) <- V3(a,e,i,m)
       call map (wrk,wrksize,
     & 4,1,4,2,3,mapdv3,mapiv3,1,mapdv2,mapiv2,possv20,posst,
     &           rc)
cJ.4.3mult V4(a,i,b,j) <- V2(a,i,m,e) . V1(m,e,b,j)
       call mult (wrk,wrksize,
     & 4,4,4,2,mapdv2,mapiv2,1,mapdv1,mapiv1,1,mapdv4,mapiv4,
     &            ssc,possv40,rc)
cJ.4.4map V3(a,b,i,j) <- V4(a,i,b,j)
       call map (wrk,wrksize,
     & 4,1,3,2,4,mapdv4,mapiv4,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
cJ.4.5pack V3(ab,ij) <- V2(ab,i,j) <- V3(a,b,i,j)
       call fack (wrk,wrksize,
     & 4,1,mapdv3,1,mapiv3,mapdv2,mapiv2,possv20,rc)
       call fack (wrk,wrksize,
     & 4,4,mapdv2,1,mapiv2,mapdv3,mapiv3,possv30,rc)
cJ.4.6add T2n(ab,ij)aaaa <- 1.0d0 V2(ab,ij)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv3,1,mapdt21,mapit21,1,rc)
c
cJ.5  T2n(a,b,i,j)abab <- sum(m,e-aa) [ T2o(b,e,j,m)bbbb . WIII(m,e,a,i)bbaa ]
cJ.5.1get V4(be,jm) <- T2o(be,jm)bbbb
       call filemanager (2,lunt2o2,rc)
       call getmediate (wrk,wrksize,
     & lunt2o2,possv40,mapdv4,mapiv4,rc)
cJ.5.2expand V3(b,e,j,m) <-  V4(be,jm)
       call expand (wrk,wrksize,
     & 4,4,mapdv4,mapiv4,1,possv30,mapdv3,mapiv3,rc)
cJ.5.2map V2(b,j,m,e) <- V3(b,e,j,m)
       call map (wrk,wrksize,
     & 4,1,4,2,3,mapdv3,mapiv3,1,mapdv2,mapiv2,possv20,posst,
     &           rc)
cJ.5.3mult V3(b,j,a,i) <- V2(b,j,m,e) . V1(m,e,a,i)
       call mult (wrk,wrksize,
     & 4,4,4,2,mapdv2,mapiv2,1,mapdv1,mapiv1,1,mapdv3,mapiv3,
     &            ssc,possv30,rc)
cJ.5.4map V2(a,b,i,j) <- V3(b,j,a,i)
       call map (wrk,wrksize,
     & 4,2,4,1,3,mapdv3,mapiv3,1,mapdv2,mapiv2,possv20,posst,
     &           rc)
cJ.5.5add T2n(a,b,i,j)abab <- 1.0d0 . V2(a,b,i,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv2,1,mapdt23,mapit23,1,rc)
cparend
       end if
c
c
c     part K - W3bbbb
cpar  all contributions are calculated for idbbbb only,
c     just Qbbbb is calculated both on idbbbb and idaabb
c
cpar
      if (myRank.eq.idbbbb) then
c
cK.1  get V1(m,e,a,j) <- W3bbbb(m,e,a,j)
       call getw3 (wrk,wrksize,
     & lunw3bbbb,2)
c
cK.2  WIII(m,e,b,j)bbbb <- sum(n-b) [ <mn||je>bbbb . T1o(n,b)bb ]
cK.2.1expand V2(j,e,m,n) <- <je||mn>bbbb
       call expand (wrk,wrksize,
     & 4,3,mapdw12,mapiw12,1,possv20,mapdv2,mapiv2,rc)
cK.2.2map V3(m,e,j,n) <- V2(j,e,m,n)
       call map (wrk,wrksize,
     & 4,3,2,1,4,mapdv2,mapiv2,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
cK.2.3mult V2(m,e,j,b) <- V3(m,e,j,n) . M2(n,b)
       call mult (wrk,wrksize,
     & 4,2,4,1,mapdv3,mapiv3,1,mapdm2,mapim2,1,mapdv2,mapiv2,
     &            ssc,possv20,rc)
cK.2.4map V3(m,e,b,j) <- V2(m,e,j,b)
       call map (wrk,wrksize,
     & 4,1,2,4,3,mapdv2,mapiv2,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
cK.2.5add V1(m,e,b,j) <- 1.0d0 . V3(m,e,b,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv3,1,mapdv1,mapiv1,1,rc)
cparend
        end if
c
cK.3  WIII4
c     Q(f,b,j,n)bbbb   <= 0.5 T2o(f,b,j,n)bbbb + T1o(f,j)bb . T1o(b,n)bb
c     WIII(m,e,b,j)bbbb <- -sum(n,f-bb)     [ Q(f,b,j,n)bbbb   . <ef||mn>bbbb ]
c     <- 0.5 sum(n,f-aa)  [ T2o(f,b,n,j)abab . <fe||nm>abab ]
c
cpar
        if ((myRank.eq.idbbbb).or.(myRank.eq.idaabb)) then
cK.3.1get V3(fb,jn) <- T2o(fb,jn)bbbb
       call filemanager (2,lunt2o2,rc)
       call getmediate (wrk,wrksize,
     & lunt2o2,possv30,mapdv3,mapiv3,rc)
cK.3.2expand V2(f,b,j,n) <- V3(fb,jn)
       call expand (wrk,wrksize,
     & 4,4,mapdv3,mapiv3,1,possv20,mapdv2,mapiv2,rc)
cK.3.3mkQ V2(f,b,j,n) <- 0.5 V2(f,b,j,n) + T1o(f,j)bb . T1o(b,n)bb
       call mkq(wrk,wrksize,
     & mapdv2,mapiv2,mapdt12,mapit12,mapdt12,mapit12,0.5d0,rc)
cK.3.4map V3(f,n,b,j) <- V2(f,b,j,n)
       call map (wrk,wrksize,
     & 4,1,3,4,2,mapdv2,mapiv2,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
cparend
        end if
cpar
        if (myRank.eq.idaabb) then
cK.3.5write V3(f,n,b,j) to lunqbbbb
       call filemanager (1,lunqbbbb,rc)
       call wrtmediate (wrk,wrksize,
     & lunqbbbb,mapdv3,mapiv3,rc)
        end if
cpar
        if (myRank.eq.idbbbb) then
cK.3.6get V2(ef,mn) = <ef||mn>bbbb from luna file
       call filemanager (2,lunabij2,rc)
       call getmediate (wrk,wrksize,
     & lunabij2,possv20,mapdv2,mapiv2,rc)
cK.3.7expand V4(e,f,m,n) <- V2(ef,mn)
       call expand (wrk,wrksize,
     & 4,4,mapdv2,mapiv2,1,possv40,mapdv4,mapiv4,rc)
cK.3.8map V2(m,e,f,n) <- V4(e,f,m,n)
       call map (wrk,wrksize,
     & 4,2,3,1,4,mapdv4,mapiv4,1,mapdv2,mapiv2,possv20,posst,
     &           rc)
cK.3.9mult V4(m,e,b,j) <- V2(m,e,f,n) . V3(f,n,b,j)
       call mult (wrk,wrksize,
     & 4,4,4,2,mapdv2,mapiv2,1,mapdv3,mapiv3,1,mapdv4,mapiv4,
     &            ssc,possv40,rc)
cK.3.10 add V1(m,e,b,j) <- -1.0d0 . V4(m,e,b,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,-1.0d0,mapdv4,1,mapdv1,mapiv1,1,rc)
cK.3.11 get V2(f,e,n,m) <- <fe||nm>abab
       call filemanager (2,lunabij3,rc)
       call getmediate (wrk,wrksize,
     & lunabij3,possv20,mapdv2,mapiv2,rc)
cK.3.12 map V3(m,e,f,n) <- V2(f,e,n,m)
       call map (wrk,wrksize,
     & 4,3,2,4,1,mapdv2,mapiv2,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
cK.3.13 get V2(f,b,n,j) <- T2o(f,b,n,j)abab
       call filemanager (2,lunt2o3,rc)
       call getmediate (wrk,wrksize,
     & lunt2o3,possv20,mapdv2,mapiv2,rc)
cK.3.14 map V4(f,n,b,j) <- V2(f,b,n,j)
       call map (wrk,wrksize,
     & 4,1,3,2,4,mapdv2,mapiv2,1,mapdv4,mapiv4,possv40,posst,
     &           rc)
cK.3.15 mult V2(m,e,b,j) <- V3(m,e,f,n) . V4(f,n,b,j)
       call mult (wrk,wrksize,
     & 4,4,4,2,mapdv3,mapiv3,1,mapdv4,mapiv4,1,mapdv2,mapiv2,
     &            ssc,possv20,rc)
cK.3.16 add V1(m,e,b,j) <- 0.5d0 . V2(m,e,b,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,0.5d0,mapdv2,1,mapdv1,mapiv1,1,rc)
c
cK.4  R1(a,i,b,j)bbbb <= sum(m,e-bb) [ T2o(a,e,i,m)bbbb . WIII(m,e,b,j)bbbb ]
c     T2n(ab,ij)bbbb   <= {1(a,i,b,j)-R1(b,i,a,j)-R1(a,j,b,i)+R1(b,j,a,i)}bbbb
cK.4.1get V2(ae,im) <- T2o(ae,im)bbbb
       call filemanager (2,lunt2o2,rc)
       call getmediate (wrk,wrksize,
     & lunt2o2,possv20,mapdv2,mapiv2,rc)
cK.4.2expand V3(a,e,i,m) <- V2(ae,im)
       call expand (wrk,wrksize,
     & 4,4,mapdv2,mapiv2,1,possv30,mapdv3,mapiv3,rc)
cK.4.3map V2(a,i,m,e) <- V3(a,e,i,m)
       call map (wrk,wrksize,
     & 4,1,4,2,3,mapdv3,mapiv3,1,mapdv2,mapiv2,possv20,posst,
     &           rc)
cK.4.4mult V4(a,i,b,j) <- V2(a,i,m,e) . V1(m,e,b,j)
       call mult (wrk,wrksize,
     & 4,4,4,2,mapdv2,mapiv2,1,mapdv1,mapiv1,1,mapdv4,mapiv4,
     &            ssc,possv40,rc)
cK.4.5map V3(a,b,i,j) <- V4(a,i,b,j)
       call map (wrk,wrksize,
     & 4,1,3,2,4,mapdv4,mapiv4,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
cK.4.6pack V3(ab,ij) <- V2(ab,i,j) <- V3(a,b,i,j)
       call fack (wrk,wrksize,
     & 4,1,mapdv3,1,mapiv3,mapdv2,mapiv2,possv20,rc)
       call fack (wrk,wrksize,
     & 4,4,mapdv2,1,mapiv2,mapdv3,mapiv3,possv30,rc)
cK.4.7add T2n(ab,ij)bbbb <- 1.0d0 V3(ab,ij)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv3,1,mapdt22,mapit22,1,rc)
c
c
cK.5  T2n(a,b,i,j)abab <- sum(m,e-bb) [ T2o(a,e,i,m)abab . WIII(m,e,b,j)bbbb ]
cK.5.1get V4(a,e,i,m) <- T2o(a,e,i,m)abab
       call filemanager (2,lunt2o3,rc)
       call getmediate (wrk,wrksize,
     & lunt2o3,possv40,mapdv4,mapiv4,rc)
cK.5.2map V2(a,i,m,e) <- V4(a,e,i,m)
       call map (wrk,wrksize,
     & 4,1,4,2,3,mapdv4,mapiv4,1,mapdv2,mapiv2,possv20,posst,
     &           rc)
cK.5.3mult V4(a,i,b,j) <- V2(a,i,m,e) . V1(m,e,b,j)
       call mult (wrk,wrksize,
     & 4,4,4,2,mapdv2,mapiv2,1,mapdv1,mapiv1,1,mapdv4,mapiv4,
     &            ssc,possv40,rc)
cK.5.4map V3(a,b,i,j) <- V4(a,i,b,j)
       call map (wrk,wrksize,
     & 4,1,3,2,4,mapdv4,mapiv4,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
cK.5.5add T2n(a,b,i,j)abab <- 1.0d0 . V3(a,b,i,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv3,1,mapdt23,mapit23,1,rc)
cparend
       end if
c
c
c     L part W3aabb
c
cpar
      if (myRank.eq.idaabb) then
c
cL.1  get V1(m,e,a,j) <- W3aabb(m,e,a,j)
       call getw3 (wrk,wrksize,
     & lunw3aabb,3)
c
cL.2  WIII(m,e,b,j)aabb <- sum(n-b) [ <mn||je>abba . T1o(n,b)bb ]
cL.2.1map V3(m,e,j,n) <- <j,e||m,n>baab
       call map (wrk,wrksize,
     & 4,3,2,1,4,mapdw14,mapiw14,1,mapdv3,mapiv3,possv30,
     &           posst,rc)
cL.2.2mult V4(m,e,j,b) <- V3(m,e,j,n) . M2(n,b)
       call mult (wrk,wrksize,
     & 4,2,4,1,mapdv3,mapiv3,1,mapdm2,mapim2,1,mapdv4,mapiv4,
     &            ssc,possv40,rc)
cL.2.3map V3(m,e,b,j) <- V4(m,e,j,b)
       call map (wrk,wrksize,
     & 4,1,2,4,3,mapdv4,mapiv4,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
cL.2.4add V1(m,e,b,j) <- 1.0d0 V3(m,e,b,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv3,1,mapdv1,mapiv1,1,rc)
c
cL.2  WIII(m,e,b,j)aabb <- 0.5 sum(n,f-aa)  [ T2o(f,b,n,j)abab . <ef||mn>aaaa ]
c     <- - sum(n,f-bb)    [ Q(f,b,j,n)bbbb   . <ef||mn>abab ]
cL.2.1get V2(f,n,b,j) from lunqbbbb (produced in K step) and close it
       call filemanager (2,lunqbbbb,rc)
       call getmediate (wrk,wrksize,
     & lunqbbbb,possv20,mapdv2,mapiv2,rc)
       call filemanager (3,lunqbbbb,rc)
cL.2.2get V3(e,f,m,n) <- <ef||mn>abab
       call filemanager (2,lunabij3,rc)
       call getmediate (wrk,wrksize,
     & lunabij3,possv30,mapdv3,mapiv3,rc)
cL.2.3map V4(m,e,f,n) <- V3(e,f,m,n)
       call map (wrk,wrksize,
     & 4,2,3,1,4,mapdv3,mapiv3,1,mapdv4,mapiv4,possv40,posst,
     &           rc)
cL.2.4mult V3(m,e,b,j) <- V4(m,e,f,n) . V2(f,n,b,j)
       call mult (wrk,wrksize,
     & 4,4,4,2,mapdv4,mapiv4,1,mapdv2,mapiv2,1,mapdv3,mapiv3,
     &            ssc,possv30,rc)
cL.2.5add V1(m,e,b,j) <- -1.0d0 . V3(m,e,b,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,-1.0d0,mapdv3,1,mapdv1,mapiv1,1,rc)
cL.2.6get V2(f,b,n,j) <- T2o(f,b,n,j)abab
       call filemanager (2,lunt2o3,rc)
       call getmediate (wrk,wrksize,
     & lunt2o3,possv20,mapdv2,mapiv2,rc)
cL.2.7map V3(f,n,b,j) <- V2(f,b,n,j)
       call map (wrk,wrksize,
     & 4,1,3,2,4,mapdv2,mapiv2,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
cL.3.8get V2(ef,mn) = <ef||mn>aaaa from lunb file
       call filemanager (2,lunabij1,rc)
       call getmediate (wrk,wrksize,
     & lunabij1,possv20,mapdv2,mapiv2,rc)
cL.3.9exp V4(e,f,m,n) <- V2(ef,mn)
       call expand (wrk,wrksize,
     & 4,4,mapdv2,mapiv2,1,possv40,mapdv4,mapiv4,rc)
cL.3.10        map V2(m,e,f,n) <- V4(e,f,m,n)
       call map (wrk,wrksize,
     & 4,2,3,1,4,mapdv4,mapiv4,1,mapdv2,mapiv2,possv20,posst,
     &           rc)
cL.3.11 mult V4(m,e,b,j) <- V2(m,e,f,n) . V3(f,n,b,j)
       call mult (wrk,wrksize,
     & 4,4,4,2,mapdv2,mapiv2,1,mapdv3,mapiv3,1,mapdv4,mapiv4,
     &            ssc,possv40,rc)
cL.3.12        add V1(m,e,b,j) <- 0.5d0 . V4(m,e,b,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,0.5d0,mapdv4,1,mapdv1,mapiv1,1,rc)
c
cL.4  R1(a,i,b,j)bbbb <- sum(m,e-aa) [ T2o(e,a,m,i)abab . WIII(m,e,b,j)aabb ]
c     T2n(ab,ij)bbbb   <= {1(a,i,b,j)-R1(b,i,a,j)-R1(a,j,b,i)+R1(b,j,a,i)}bbbb
cL.4.1get V3(e,a,m,i) <- T2o(e,a,m,i)abab
       call filemanager (2,lunt2o3,rc)
       call getmediate (wrk,wrksize,
     & lunt2o3,possv30,mapdv3,mapiv3,rc)
cL.4.2map V2(a,i,m,e) <- V3(e,a,m,i)
       call map (wrk,wrksize,
     & 4,4,1,3,2,mapdv3,mapiv3,1,mapdv2,mapiv2,possv20,posst,
     &           rc)
cL.4.3mult V4(a,i,b,j) <- V2(a,i,m,e) . V1(m,e,b,j)
       call mult (wrk,wrksize,
     & 4,4,4,2,mapdv2,mapiv2,1,mapdv1,mapiv1,1,mapdv4,mapiv4,
     &            ssc,possv40,rc)
cL.4.4map V3(a,b,i,j) <- V4(a,i,b,j)
       call map (wrk,wrksize,
     & 4,1,3,2,4,mapdv4,mapiv4,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
cL.4.5pack V3(ab,ij) <- V2(ab,i,j) <- V3(a,b,i,j)
       call fack (wrk,wrksize,
     & 4,1,mapdv3,1,mapiv3,mapdv2,mapiv2,possv20,rc)
       call fack (wrk,wrksize,
     & 4,4,mapdv2,1,mapiv2,mapdv3,mapiv3,possv30,rc)
cL.4.6add T2n(ab,ij)aaaa <- 1.0d0 V2(ab,ij)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv3,1,mapdt22,mapit22,1,rc)
c
cL.5  T2n(a,b,i,j)abab <- sum(m,e-aa) [ T2o(a,e,i,m)aaaa . WIII(m,e,b,j)aabb ]
cL.5.1get V4(ae,im) <- T2o(ae,im)aaaa
       call filemanager (2,lunt2o1,rc)
       call getmediate (wrk,wrksize,
     & lunt2o1,possv40,mapdv4,mapiv4,rc)
cL.5.2expand V3(a,e,i,m) <-  V4(ae,im)
       call expand (wrk,wrksize,
     & 4,4,mapdv4,mapiv4,1,possv30,mapdv3,mapiv3,rc)
cL.5.2map V2(a,i,m,e) <- V3(a,e,i,m)
       call map (wrk,wrksize,
     & 4,1,4,2,3,mapdv3,mapiv3,1,mapdv2,mapiv2,possv20,posst,
     &           rc)
cL.5.3mult V3(a,i,b,j) <- V2(a,i,m,e) . V1(m,e,b,j)
       call mult (wrk,wrksize,
     & 4,4,4,2,mapdv2,mapiv2,1,mapdv1,mapiv1,1,mapdv3,mapiv3,
     &            ssc,possv30,rc)
cL.5.4map V2(a,b,i,j) <- V3(a,i,b,j)
       call map (wrk,wrksize,
     & 4,1,3,2,4,mapdv3,mapiv3,1,mapdv2,mapiv2,possv20,posst,
     &           rc)
cL.5.5add T2n(a,b,i,j)abab <- 1.0d0 . V2(a,b,i,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv2,1,mapdt23,mapit23,1,rc)
cparend
       end if
c
c
c     M part W3abba
c
cpar
      if (myRank.eq.idabba) then
c
cM.1  get V1(m,e,a,j) <- W3abba
       call getw3 (wrk,wrksize,
     & lunw3abba,4)
c
cM.2  WIII(m,e,b,j)abba <- sum(n-b) [ <mn||je>abab . T1o(n,b)bb ]
cM.2.1map V3(m,e,j,n) <- <j,e||m,n>abab
       call map (wrk,wrksize,
     & 4,3,2,1,4,mapdw13,mapiw13,1,mapdv3,mapiv3,possv30,
     &           posst,rc)
cM.2.2mult V4(m,e,j,b) <- V3(m,e,j,n) . M2(n,b)
       call mult (wrk,wrksize,
     & 4,2,4,1,mapdv3,mapiv3,1,mapdm2,mapim2,1,mapdv4,mapiv4,
     &            ssc,possv40,rc)
cM.2.3map V3(m,e,b,j) <- V4(m,e,j,b)
       call map (wrk,wrksize,
     & 4,1,2,4,3,mapdv4,mapiv4,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
cM.2.4add V1(m,e,b,j) <- 1.0d0 V3(m,e,b,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv3,1,mapdv1,mapiv1,1,rc)
c
cM.3  Q(f,b,j,n)abab   <= 0.5 T2o(f,b,j,n)abab + T1o(f,j)aa . T1o(b,n)bb
c     WIII(m,e,b,j)abba <- sum(n,f-ba)      [ Q(f,b,j,n)abab   . <fe||mn>abab ]
cM.3.1get V4(f,b,j,n) <- T2o(f,b,j,n)abab
       call filemanager (2,lunt2o3,rc)
       call getmediate (wrk,wrksize,
     & lunt2o3,possv40,mapdv4,mapiv4,rc)
cM.3.2mkQ V4(f,b,j,n) <- 0.5 V4(f,b,j,n) + T1o(f,j)aa . T1o(b,n)bb
       call mkq(wrk,wrksize,
     & mapdv4,mapiv4,mapdt11,mapit11,mapdt12,mapit12,0.5d0,rc)
cM.3.3map V3(f,n,b,j) <- V4(f,b,j,n)
       call map (wrk,wrksize,
     & 4,1,3,4,2,mapdv4,mapiv4,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
cM.3.4get V2(f,e,m,n) <- <fe||mn>abab
       call filemanager (2,lunabij3,rc)
       call getmediate (wrk,wrksize,
     & lunabij3,possv20,mapdv2,mapiv2,rc)
cM.3.5map V4(m,e,f,n) <- V2(f,e,m,n)
       call map (wrk,wrksize,
     & 4,3,2,1,4,mapdv2,mapiv2,1,mapdv4,mapiv4,possv40,posst,
     &           rc)
cM.3.6mult V2(m,e,b,j) <- V4(m,e,f,n) . V3(f,n,b,j)
       call mult (wrk,wrksize,
     & 4,4,4,2,mapdv4,mapiv4,1,mapdv3,mapiv3,1,mapdv2,mapiv2,
     &            ssc,possv20,rc)
cM.3.7add V1(m,e,b,j) <- 1.0d0 V2(m,e,b,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv2,1,mapdv1,mapiv1,1,rc)
c
cM.4  T2n(a,b,i,j)abab <- sum(m,e-ab) [ T2o(a,e,m,j)abab . WIII(m,e,b,i)abba ]
cM.4.1get V2(a,e,m,j) <- T2o(a,e,m,j)abab
       call filemanager (2,lunt2o3,rc)
       call getmediate (wrk,wrksize,
     & lunt2o3,possv20,mapdv2,mapiv2,rc)
cM.4.2map V3(a,j,m,e) <- V2(a,e,m,j)
       call map (wrk,wrksize,
     & 4,1,4,3,2,mapdv2,mapiv2,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
cM.4.3mult V4(a,j,b,i) <- V3(a,j,m,e) . V1(m,e,b,i)
       call mult (wrk,wrksize,
     & 4,4,4,2,mapdv3,mapiv3,1,mapdv1,mapiv1,1,mapdv4,mapiv4,
     &            ssc,possv40,rc)
cM.4.4map V3(a,b,i,j) <- V4(a,j,b,i)
       call map (wrk,wrksize,
     & 4,1,4,2,3,mapdv4,mapiv4,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
cM.4.5add T2n(a,b,i,j)abab <- 1.0d0 V3(a,b,i,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv3,1,mapdt23,mapit23,1,rc)
cparend
       end if
c
c
c     N part W3baab
c
cpar
      if (myRank.eq.idbaab) then
c
cN.1  get V1(m,e,a,j) <- W3baab
       call getw3 (wrk,wrksize,
     & lunw3baab,5)
c
cN.2  WIII(m,e,b,j)baab <- - sum(n-a) [ <je||nm>baab . T1o(n,b)aa ]
cN.2.1map V3(m,e,j,n) <- <j,e||n,m>baab
       call map (wrk,wrksize,
     & 4,3,2,4,1,mapdw14,mapiw14,1,mapdv3,mapiv3,possv30,
     &           posst,rc)
cN.2.2mult V4(m,e,j,b) <- V3(m,e,j,n) . M1(n,b)
       call mult (wrk,wrksize,
     & 4,2,4,1,mapdv3,mapiv3,1,mapdm1,mapim1,1,mapdv4,mapiv4,
     &            ssc,possv40,rc)
cN.2.3map V3(m,e,b,j) <- V4(m,e,j,b)
       call map (wrk,wrksize,
     & 4,1,2,4,3,mapdv4,mapiv4,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
cN.2.4add V1(m,e,b,j) <- -1.0d0 V3(m,e,b,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,-1.0d0,mapdv3,1,mapdv1,mapiv1,1,rc)
c
cN.3  Q(b,f,n,j)abab   <= 0.5 T2o(b,f,n,j)abab + T1o(b,n)aa . T1o(f,j)bb
c     WIII(m,e,b,j)baab <- sum(n,f-ab)      [ Q(b,f,n,j)abab   . <ef||nm>abab ]
cN.3.1get V4(b,f,n,j) <- T2o(b,f,n,j)abab
       call filemanager (2,lunt2o3,rc)
       call getmediate (wrk,wrksize,
     & lunt2o3,possv40,mapdv4,mapiv4,rc)
cN.3.2mkQ V4(b,f,n,j) <- 0.5 V4(b,f,n,j) + T1o(b,n)aa . T1o(f,j)bb
       call mkq(wrk,wrksize,
     & mapdv4,mapiv4,mapdt11,mapit11,mapdt12,mapit12,0.5d0,rc)
cN.3.3map V3(f,n,b,j) <- V4(b,f,n,j)
       call map (wrk,wrksize,
     & 4,3,1,2,4,mapdv4,mapiv4,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
cN.3.4get V2(e,f,n,m) <- <ef||nm>abab
       call filemanager (2,lunabij3,rc)
       call getmediate (wrk,wrksize,
     & lunabij3,possv20,mapdv2,mapiv2,rc)
cN.3.5map V4(m,e,f,n) <- V2(e,f,n,m)
       call map (wrk,wrksize,
     & 4,2,3,4,1,mapdv2,mapiv2,1,mapdv4,mapiv4,possv40,posst,
     &           rc)
cN.3.6mult V2(m,e,b,j) <- V4(m,e,f,n) . V3(f,n,b,j)
       call mult (wrk,wrksize,
     & 4,4,4,2,mapdv4,mapiv4,1,mapdv3,mapiv3,1,mapdv2,mapiv2,
     &            ssc,possv20,rc)
cN.3.7add V1(m,e,b,j) <- 1.0d0 V2(m,e,b,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv2,1,mapdv1,mapiv1,1,rc)
c
cN.4  T2n(a,b,i,j)abab <- sum(m,e-ab) [ T2o(e,b,i,m)abab . WIII(m,e,a,j)baab ]
cN.4.1get V2(e,b,i,m) <- T2o(e,b,i,m)abab
       call filemanager (2,lunt2o3,rc)
       call getmediate (wrk,wrksize,
     & lunt2o3,possv20,mapdv2,mapiv2,rc)
cN.4.2map V3(b,i,m,e) <- V2(e,b,i,m)
       call map (wrk,wrksize,
     & 4,4,1,2,3,mapdv2,mapiv2,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
cN.4.3mult V4(b,i,a,j) <- V3(b,i,m,e) . V1(m,e,a,j)
       call mult (wrk,wrksize,
     & 4,4,4,2,mapdv3,mapiv3,1,mapdv1,mapiv1,1,mapdv4,mapiv4,
     &            ssc,possv40,rc)
cN.4.4map V3(a,b,i,j) <- V4(b,i,a,j)
       call map (wrk,wrksize,
     & 4,2,3,1,4,mapdv4,mapiv4,1,mapdv3,mapiv3,possv30,posst,
     &           rc)
cN.4.5add T2n(a,b,i,j)abab <- 1.0d0 V3(a,b,i,j)
       call add (wrk,wrksize,
     & 4,4,0,0,0,0,1,1,1.0d0,mapdv3,1,mapdt23,mapit23,1,rc)
cparend
       end if
c
       return
       end
c
c     --------------------------
c
       subroutine getw3 (wrk,wrksize,
     & lunw3xxxx,nxxxx)
c
c     This routine reconstruct W3(m,e,a,j)xxxx from lunw3xxxx file
c     to V1(m,e,a,j) and define corresponding mapdv1 and mapiv1
c     This routine also close lunw3xxxx file
c
c     lunw3xxxx - lun of opened w3xxxx file
c     nxxxx     - xxxx identifier
c     1 - aaaa
c     2 - bbbb
c     3 - aabb
c     4 - abba
c     5 - baab
c     6 - bbaa
c
c     the structure of lunw3xxxx is:
c
c     do syma=1,nsym
c     map of H _a(m,e,j)bbbb (W3(m,e,a,j))
c     skip cycle over a if length of all files is zero
c     do a=1,nvx(syma) [ x is a or b ]
c     if (h1length.gt.0) then
c     write H(m,e,j)
c     end if
c     end do
c     end do
c
#include "ccsd1.fh"
#include "ccsd2.fh"
#include "wrk.fh"
c
       integer lunw3xxxx,nxxxx
c
c     help variables
c
       integer rc,syma,a,h1length,posst,aup,aalfayes,iiv1,v1length
c
c0    def aalfayes
c
       if ((nxxxx.eq.1).or.(nxxxx.eq.5).or.(nxxxx.eq.6)) then
       aalfayes=1
       else
       aalfayes=0
       end if
c
c1.1  def maps of V1(m,e,a,j)
c
       if (nxxxx.eq.1) then
       call grc0 (4,0,1,3,3,1,1,possv10,posst,mapdv1,mapiv1)
       else if (nxxxx.eq.2) then
       call grc0 (4,0,2,4,4,2,1,possv10,posst,mapdv1,mapiv1)
       else if (nxxxx.eq.3) then
       call grc0 (4,0,1,3,4,2,1,possv10,posst,mapdv1,mapiv1)
       else if (nxxxx.eq.4) then
       call grc0 (4,0,1,4,4,1,1,possv10,posst,mapdv1,mapiv1)
       else if (nxxxx.eq.5) then
       call grc0 (4,0,2,3,3,2,1,possv10,posst,mapdv1,mapiv1)
       else if (nxxxx.eq.6) then
       call grc0 (4,0,2,4,3,1,1,possv10,posst,mapdv1,mapiv1)
       end if
c
c1.2  vanish V1
       iiv1=mapdv1(0,5)
       v1length=mapdv1(iiv1,1)+mapdv1(iiv1,2)-possv10
       call mv0zero (v1length,v1length,wrk(possv10))
c
c2    rewind tape lunw3xxxx
       call filemanager (2,lunw3xxxx,rc)
c
c3    loop over symA
c
       do 3000 syma=1,nsym
c
c3.1  get map of H _a(m,e,j) to mapd,i H1
       call getmap (lunw3xxxx,possh10,h1length,mapdh1,mapih1,rc)
c
c3.2  skip cycle over a if length of H1 is 0
       if (h1length.eq.0) goto 3000
c
c3.3  loop over all a in this symmetry
c
       if (aalfayes.eq.1) then
       aup=nva(syma)
       else
       aup=nvb(syma)
       end if
c
       do 2500 a=1,aup
c
       if (h1length.gt.0) then
c
c3.3.1read H1 if any
       call rea (lunw3xxxx,h1length,wrk(possh10))
c
c3.3.2insert H1 into V1 for given a and syma
       call add (wrk,wrksize,
     & 3,4,1,3,a,0,syma,1,1.0d0,mapdh1,syma,mapdv1,mapiv1,1,
     &           rc)
c
       end if
c
 2500   continue
c
 3000   continue
c
c4    close lunw3xxxx file
       call filemanager (3,lunw3xxxx,rc)
c
       return
       end

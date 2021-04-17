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

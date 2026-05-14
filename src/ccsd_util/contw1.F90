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

subroutine contw1(wrk,wrksize,lunabij1,lunabij2,lunabij3,lunt2o1,lunt2o2,lunt2o3)
! this routine calcs contributions to W1 i.e. W11,W12,W13
! and contribution from W1 to T2 i.e. T24
!
! Mediate : WI(m,n,i,j)
! Spin states : WI(mn,ij)aaaa, WI(mn,ij)bbbb, WI(m,n,i,j)abab
! Contributions:
!
!     WI1
!I.1  WI(mn,ij)aaaa   <= <mn||ij>aaaa
!J.1  WI(mn,ij)bbbb   <= <mn||ij>bbbb
!K.1  WI(m,n,i,j)abab <= <mn||ij>abab
!
!     WI2
!I.2  Q(mn,i,j)aaaa  <= sum(e-a) [ <ie||mn>aaaa . T1o(e,j)aa ]
!     WI(mn,ij)aaaa   <- Q(mn,i,j)aaaa - Q(mn,j,i)aaaa  i>j
!J.2  Q(mn,i,j)bbbb  <= sum(e-b) [ <ie||mn>bbbb . T1o(e,j)bb ]
!     WI(mn,ij)bbbb   <- Q(mn,i,j)bbbb - Q(mn,j,i)bbbb  i>j
!K.2  WI(m,n,i,j)abab <- - sum(e-a) [ <je||mn>baab . T1o(e,i)aa ]
!     <- sum(e-b) [ <ie||mn>abab . T1o(e,j)bb ]
!
!     WI3
!I.3  WI(mn,ij)aaaa   <- sum(e>f-aa) [ <ef||mn>aaaa . Tau(ef,ij)aaaa ]
!J.3  WI(mn,ij)bbbb   <- sum(e>f-bb) [ <ef||mn>bbbb . Tau(ef,ij)bbbb ]
!K.3  WI(m,n,i,j)abab <- sum(e,f-ab) [ <ef||mn>abab . Tau(e,f,i,j)abab ]
!
!     T24
!I.4  T2n(ab,ij)aaaa   <- sum(m>n-aa) [ Tau(ab,mn)aaaa . W1(mn,ij)aaaa ]
!J.4  T2n(ab,ij)bbbb   <- sum(m>n-bb) [ Tau(ab,mn)bbbb . W1(mn,ij)bbbb ]
!K.4  T2n(a,b,i,j)abab <- sum(m,n-ab) [ Tau(a,b,m,n)abab . W1(m,n,i,j)abab ]
!
! N.B. use and destroy : V1,V2,V3,V4
! N.B. number of read  : 6

use ccsd_global, only: idfin, t11, t12, t21, t22, t23, v1, v2, v3, v4, w01, w02, w03, w11, w12, w13, w14
use Para_Info, only: MyRank
use Constants, only: One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize
real(kind=wp), intent(inout) :: wrk(wrksize)
integer(kind=iwp), intent(_IN_) :: lunabij1, lunabij2, lunabij3, lunt2o1, lunt2o2, lunt2o3
integer(kind=iwp) :: post, rc, ssc

!par
if (myRank == idfin) then

  !I case W1(mn,ij)aaaa

  !I.1.1 map V1(mn,ij) <= <mn||ij>aaaa
  call map(wrk,wrksize,4,1,2,3,4,w01,1,v1,post,rc)

  !I.2.1 map V2(mn,i,e) <= <ie||mn>aaaa
  call map(wrk,wrksize,4,3,4,1,2,w11,1,v2,post,rc)

  !I.2.2 mult V3(mn,i,j) <= V2(mn,i,e) . T1o(e,j)aa
  call ccmult(wrk,wrksize,4,2,4,1,v2,1,t11,1,v3,ssc,rc)

  !I.2.3 pack V2(mn,ij) <= V3(mn,i,j)
  call fack(wrk,wrksize,4,4,v3,1,v2,rc)

  !I.2.4 add V1(mn,ij) <- V2(mn,ij)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,One,v2,1,v1,1,rc)

  !I.3.1 read V2(ef,ij) <= T2o(ef,ij)aaaa
  call filemanager(2,lunt2o1,rc)
  call getmediate(wrk,wrksize,lunt2o1,v2,rc)

  !I.3.2 make Tau V2(ef,ij) from V2(ef,ij)
  call mktau(wrk,wrksize,v2,t11,t12,One,rc)

  !I.3.3 read V3(ef,mn) <= <ef||mn>aaaa
  call filemanager(2,lunabij1,rc)
  call getmediate(wrk,wrksize,lunabij1,v3,rc)

  !I.3.4 map V4(mn,ef) <= V3(ef,mn)
  call map(wrk,wrksize,4,3,4,1,2,v3,1,v4,post,rc)

  !I.3.5 mult V3(mn,ij) = V4(mn,ef) . V2(ef,ij)
  call ccmult(wrk,wrksize,4,4,4,2,v4,1,v2,1,v3,ssc,rc)

  !I.3.6 add V1(mn,ij) <- V3(mn,ij)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,One,v3,1,v1,1,rc)

  !I.4.0 Tau(ab,mn) are in V2(ab,mn) from I.3.2
  !      W1(mn,ij)  are in V1(mn,ij)

  !I.4.1 V3(ab,ij) = V2(ab,mn) . V1(mn,ij)
  call ccmult(wrk,wrksize,4,4,4,2,v2,1,v1,1,v3,ssc,rc)

  !I.4.2 add t2n(ab,ij)aaaa <- V3(ab,ij)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,One,v3,1,t21,1,rc)

  !J case W1(mn,ij)bbbb

  !J.1.1 map V1(mn,ij) <= <mn||ij>bbbb
  call map(wrk,wrksize,4,1,2,3,4,w02,1,v1,post,rc)

  !J.2.1 map V2(mn,i,e) <= <ie||mn>bbbb
  call map(wrk,wrksize,4,3,4,1,2,w12,1,v2,post,rc)

  !J.2.2 mult V3(mn,i,j) <= V2(mn,i,e) . T1o(e,j)bb
  call ccmult(wrk,wrksize,4,2,4,1,v2,1,t12,1,v3,ssc,rc)

  !J.2.3 pack V2(mn,ij) <= V3(mn,i,j)
  call fack(wrk,wrksize,4,4,v3,1,v2,rc)

  !J.2.4 add V1(mn,ij) <- V2(mn,ij)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,One,v2,1,v1,1,rc)

  !J.3.1 read V2(ef,ij) <= T2o(ef,ij)bbbb
  call filemanager(2,lunt2o2,rc)
  call getmediate(wrk,wrksize,lunt2o2,v2,rc)

  !J.3.2 make Tau V2(ef,ij) from V2(ef,ij)
  call mktau(wrk,wrksize,v2,t11,t12,One,rc)

  !J.3.3 read V3(ef,mn) <= <ef||mn>bbbb
  call filemanager(2,lunabij2,rc)
  call getmediate(wrk,wrksize,lunabij2,v3,rc)

  !J.3.4 map V4(mn,ef) <= V3(ef,mn)
  call map(wrk,wrksize,4,3,4,1,2,v3,1,v4,post,rc)

  !J.3.5 mult V3(mn,ij) = V4(mn,ef) . V2(ef,ij)
  call ccmult(wrk,wrksize,4,4,4,2,v4,1,v2,1,v3,ssc,rc)

  !J.3.6 add V1(mn,ij) <- V3(mn,ij)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,One,v3,1,v1,1,rc)

  !J.4.0 Tau(ab,mn) are in V2(ab,mn) from J.3.2
  !      W1(mn,ij)  are in V1(mn,ij)

  !J.4.1 mult V3(ab,ij) = V2(ab,mn) . V1(mn,ij)
  call ccmult(wrk,wrksize,4,4,4,2,v2,1,v1,1,v3,ssc,rc)

  !J.4.2 add t2n(ab,ij)bbbb <- V3(ab,ij)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,One,v3,1,t22,1,rc)

  !K case W1(mn,ij)abab

  !K.1.1 map V1(m,n,i,j) <= <mn||ij>abab
  call map(wrk,wrksize,4,1,2,3,4,w03,1,v1,post,rc)

  !K.2.1 map V2(m,n,j,e) <= <je||mn>baab
  call map(wrk,wrksize,4,3,4,1,2,w14,1,v2,post,rc)

  !K.2.2 mult V3(m,n,j,i) <= V2(m,n,j,e) . T1o(e,i)aa
  call ccmult(wrk,wrksize,4,2,4,1,v2,1,t11,1,v3,ssc,rc)

  !K.2.3 map V2(m,n,i,j) <= V3(m,n,j,i)
  call map(wrk,wrksize,4,1,2,4,3,v3,1,v2,post,rc)

  !K.2.4 add V1(m,n,i,j) <- - V2(m,n,i,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,-One,v2,1,v1,1,rc)

  !K.2.5 map V2(m,n,i,e) <= <ie||mn>abab
  call map(wrk,wrksize,4,3,4,1,2,w13,1,v2,post,rc)

  !K.2.6 mult V3(m,n,i,j) <= V2(m,n,i,e) . T1o(e,j)bb
  call ccmult(wrk,wrksize,4,2,4,1,v2,1,t12,1,v3,ssc,rc)

  !K.2.7 add V1(m,n,i,j) <- V3(m,n,i,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,One,v3,1,v1,1,rc)

  !K.3.1 read V2(e,f,i,j) <= T2o(e,f,i,j)abab
  call filemanager(2,lunt2o3,rc)
  call getmediate(wrk,wrksize,lunt2o3,v2,rc)

  !K.3.2 make Tau V2(e,f,i,j) from V2(e,f,i,j)
  call mktau(wrk,wrksize,v2,t11,t12,One,rc)

  !K.3.3 read V3(e,f,m,n) <= <ef||mn>abab
  call filemanager(2,lunabij3,rc)
  call getmediate(wrk,wrksize,lunabij3,v3,rc)

  !K.3.4 map V4(m,n,e,f) <= V3(e,f,m,n)
  call map(wrk,wrksize,4,3,4,1,2,v3,1,v4,post,rc)

  !K.3.5 mult V3(m,n,i,j) = V4(m,n,e,f) . V2(e,f,i,j)
  call ccmult(wrk,wrksize,4,4,4,2,v4,1,v2,1,v3,ssc,rc)

  !K.3.6 add V1(m,n,i,j) <- V3(m,n,i,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,One,v3,1,v1,1,rc)

  !K.4.0 Tau(a,b,m,n) are in V2(ab,mn) from K.3.2
  !      W1(m,n,i,j)  are in V1(m,n,i,j)

  !K.4.1 mult V3(a,b,i,j) = V2(a,b,m,n) . V1(m,n,i,j)
  call ccmult(wrk,wrksize,4,4,4,2,v2,1,v1,1,v3,ssc,rc)

  !K.4.2 add t2n(a,b,i,j)abab <- V3(a,b,i,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,One,v3,1,t23,1,rc)

end if

return

end subroutine contw1

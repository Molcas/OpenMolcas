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

subroutine contf24(wrk,wrksize,lunabij1,lunabij2,lunabij3,lunt2o1,lunt2o2,lunt2o3)
! this routine does:
! 5)  FII4  f2(m,i) <- sum(n,e>f) [ <ef||mn> . Tap(in,ef)]
!
! N.B. use and destroy : V1,V2,V3,V4,M1
! N.B. # of get v2o2 : 6
! N.B. possible fusion with f14 graph

use ccsd_global, only: f21, f22, idaabb, idbaab, m1, t11, t12, v1, v2, v3, v4
use Para_Info, only: MyRank
use Constants, only: One, Half
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize
real(kind=wp), intent(inout) :: wrk(wrksize)
integer(kind=iwp), intent(_IN_) :: lunabij1, lunabij2, lunabij3, lunt2o1, lunt2o2, lunt2o3
integer(kind=iwp) :: post, rc, ssc

!1 f2(m,i)aa <- sum(n,e>f-aaa) [ <ef||mn>aaaa . Tap(ef,in)aaaa ]

!par
if (myRank == idbaab) then

  !1.1 read V1(af,mn) <= T2o(ef,in)aaaa
  call filemanager(2,lunt2o1,rc)
  call getmediate(wrk,wrksize,lunt2o1,v1,rc)

  !1.2 make Tap V1(ef,in) from V1(ef,in)
  call mktau(wrk,wrksize,v1,t11,t12,Half,rc)

  !1.3 expand V2(ef,i,n) <= V1(ef,in)
  call expand(wrk,wrksize,4,6,v1,1,v2,rc)

  !1.4 map V4(ef,n,i) <= V2(ef,i,n)
  call map(wrk,wrksize,4,1,2,4,3,v2,1,v4,post,rc)

  !1.5 read V1(ef,mn) <= <ef||mn>aaaa
  call filemanager(2,lunabij1,rc)
  call getmediate(wrk,wrksize,lunabij1,v1,rc)

  !1.6 expand V3(ef,m,n) <= V1(ef,mn)
  call expand(wrk,wrksize,4,6,v1,1,v3,rc)

  !1.7 map V1(m,ef,n) <= V3(ef,m,n)
  call map(wrk,wrksize,4,2,3,1,4,v3,1,v1,post,rc)

  !1.8 mult M1(m,i) <= V1(m,ef,n) . V4(ef,n,i)
  call ccmult(wrk,wrksize,4,4,2,3,v1,1,v4,1,m1,ssc,rc)

  !1.9 add f2(m,i)aa <- M1(m,i)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,One,m1,1,f21,1,rc)

end if

!2 f2(m,i)bb <- sum(n,e>f-bbb) [ <ef||mn>bbbb . Tap(ef,in)bbbb ]

!par
if (myRank == idaabb) then

  !2.1 read V1(af,mn) <= T2o(ef,in)bbbb
  call filemanager(2,lunt2o2,rc)
  call getmediate(wrk,wrksize,lunt2o2,v1,rc)

  !2.2 make Tap V1(ef,in) from V1(ef,in)
  call mktau(wrk,wrksize,v1,t11,t12,Half,rc)

  !2.3 expand V2(ef,i,n) <= V1(ef,in)
  call expand(wrk,wrksize,4,6,v1,1,v2,rc)

  !2.4 map V4(ef,n,i) <= V2(ef,i,n)
  call map(wrk,wrksize,4,1,2,4,3,v2,1,v4,post,rc)

  !2.5 read V1(ef,mn) <= <ef||mn>bbbb
  call filemanager(2,lunabij2,rc)
  call getmediate(wrk,wrksize,lunabij2,v1,rc)

  !2.6 expand V3(ef,m,n) <= V1(ef,mn)
  call expand(wrk,wrksize,4,6,v1,1,v3,rc)

  !2.7 map V1(m,ef,n) <= V3(ef,m,n)
  call map(wrk,wrksize,4,2,3,1,4,v3,1,v1,post,rc)

  !2.8 mult M1(m,i) <= V1(m,ef,n) . V4(ef,n,i)
  call ccmult(wrk,wrksize,4,4,2,3,v1,1,v4,1,m1,ssc,rc)

  !2.9 add f2(m,i)bb <- M1(m,i)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,One,m1,1,f22,1,rc)

end if

!3* f2(m,i)aa <- sum(nef-bab) [ <ef||mn>abab . Tap(e,f,i,n)abab ]
!4* f2(m,i)bb <- sum(nef-aab) [ <ef||nm>abab . Tap(e,f,n,i)abab ]

!par
if ((myRank == idbaab) .or. (myRank == idaabb)) then

  !34.1 read V1(e,f,k,l) <= T2o(e,f,k,l)abab
  call filemanager(2,lunt2o3,rc)
  call getmediate(wrk,wrksize,lunt2o3,v1,rc)

  !34.2 read V2(e,f,k,l) <= <ef||kl>abab
  call filemanager(2,lunabij3,rc)
  call getmediate(wrk,wrksize,lunabij3,v2,rc)

  !34.3 make Tap V1(e,f,k,l) from V1(e,f,k,l)
  call mktau(wrk,wrksize,v1,t11,t12,Half,rc)
end if

!par
if (myRank == idbaab) then

  !3.1 map V3(m,e,f,n) <= V2(e,f,m,n)
  call map(wrk,wrksize,4,2,3,1,4,v2,1,v3,post,rc)

  !3.2 map V4(e,f,n,i) <= V1(e,f,i,n)
  call map(wrk,wrksize,4,1,2,4,3,v1,1,v4,post,rc)

  !3.3 mult M1(m,i) <= V3(m,e,f,n) . V4(e,f,n,i)
  call ccmult(wrk,wrksize,4,4,2,3,v3,1,v4,1,m1,ssc,rc)

  !3.4 add f2(m,i)aa <- M1(m,i)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,One,m1,1,f21,1,rc)

end if

!par
if (myRank == idaabb) then

  !4.1 map V3(m,e,f,n) <= V2(e,f,n,m)
  call map(wrk,wrksize,4,2,3,4,1,v2,1,v3,post,rc)

  !4.3 mult M1(m,i) <= V3(m,e,f,n) . V1(e,f,n,i)
  call ccmult(wrk,wrksize,4,4,2,3,v3,1,v1,1,m1,ssc,rc)

  !4.4 add f2(m,i)bb <- M1(m,i)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,One,m1,1,f22,1,rc)

end if

return

end subroutine contf24

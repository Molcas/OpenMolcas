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

subroutine contf14(wrk,wrksize,lunabij1,lunabij2,lunabij3,lunt2o1,lunt2o2,lunt2o3)
! this routine does:
! 2)  FI4   f1(a,e) <- -sum(m>n,f) [Tap(af,mn) . <ef||mn>]
!
! N.B. use and destroy : V1,V2,V3,V4,M1
! N.B. # of get v2o2 : 6
! N.B. possible fusion with f24 graph

use ccsd_global, only: f11, f12, idaabb, idbaab, m1, t11, t12, v1, v2, v3, v4
use Para_Info, only: MyRank
use Constants, only: One, Half
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize
real(kind=wp), intent(inout) :: wrk(wrksize)
integer(kind=iwp), intent(_IN_) :: lunabij1, lunabij2, lunabij3, lunt2o1, lunt2o2, lunt2o3
integer(kind=iwp) :: post, rc, ssc

!1 f1(a,e)aa <- -sum(m>n,f-aaa) [ Tap(a,f,mn)aaaa . <ef||mn>aaaa ]

!par
if (myRank == idbaab) then

  !1.1 read V1(af,mn) <= T2o(af,mn)aaaa
  call filemanager(2,lunt2o1,rc)
  call getmediate(wrk,wrksize,lunt2o1,v1,rc)

  !1.2 make Tap V1(af,mn) from V1(af,mn)
  call mktau(wrk,wrksize,v1,t11,t12,Half,rc)

  !1.3 expand V2(a,f,mn) <= V1(af,mn)
  call expand(wrk,wrksize,4,5,v1,1,v2,rc)

  !1.4 read V1(ef,mn) <= <ef||mn>aaaa
  call filemanager(2,lunabij1,rc)
  call getmediate(wrk,wrksize,lunabij1,v1,rc)

  !1.5 expand V3(e,f,mn) <= V1(ef,mn)
  call expand(wrk,wrksize,4,5,v1,1,v3,rc)

  !1.6 map V1(f,mn,e) <= V3(e,f,mn)
  call map(wrk,wrksize,4,4,1,2,3,v3,1,v1,post,rc)

  !1.7 mult M1(a,e) <= V2(a,f,mn) . V1(f,mn,e)
  call ccmult(wrk,wrksize,4,4,2,3,v2,1,v1,1,m1,ssc,rc)

  !1.8 add f1(a,e)aa <- M1(a,e)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,-One,m1,1,f11,1,rc)

end if

!2 f1(a,e)bb <- -sum(m>n,f-bbb) [ Tap(a,f,mn)bbbb . <ef||mn>bbbb ]

!par
if (myRank == idaabb) then

  !2.1 read V1(af,mn) <= T2o(af,mn)bbbb
  call filemanager(2,lunt2o2,rc)
  call getmediate(wrk,wrksize,lunt2o2,v1,rc)

  !2.2 make Tap V1(af,mn) from V1(af,mn)
  call mktau(wrk,wrksize,v1,t11,t12,Half,rc)

  !2.3 expand V2(a,f,mn) <= V1(af,mn)
  call expand(wrk,wrksize,4,5,v1,1,v2,rc)

  !2.4 read V1(ef,mn) <= <ef||mn>bbbb
  call filemanager(2,lunabij2,rc)
  call getmediate(wrk,wrksize,lunabij2,v1,rc)

  !2.5 expand V3(e,f,mn) <= V1(ef,mn)
  call expand(wrk,wrksize,4,5,v1,1,v3,rc)

  !2.6 map V1(f,mn,e) <= V3(e,f,mn)
  call map(wrk,wrksize,4,4,1,2,3,v3,1,v1,post,rc)

  !2.7 mult M1(a,e) <= V2(a,f,mn) . V1(f,mn,e)
  call ccmult(wrk,wrksize,4,4,2,3,v2,1,v1,1,m1,ssc,rc)

  !2.8 add f1(a,e)bb <- M1(a,e)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,-One,m1,1,f12,1,rc)

end if

!3 f1(a,e)aa <- -sum(m,n,f-abb) [ Tap(a,f,m,n)abab . <ef||mn>abab
!4 f1(a,e)bb <- -sum(m,n,f-aba) [ Tap(f,a,m,n)abab . <fe||mn>abab

!par
if ((myRank == idbaab) .or. (myRank == idaabb)) then

  !34.1 read V1(c,d,m,n) <= T2o(c,d,m,n)abab
  call filemanager(2,lunt2o3,rc)
  call getmediate(wrk,wrksize,lunt2o3,v1,rc)

  !34.2 read V2(c,d,m,n) <= <cd||mn>abab
  call filemanager(2,lunabij3,rc)
  call getmediate(wrk,wrksize,lunabij3,v2,rc)

  !34.3 make Tap V1(c,d,m,n) from V1(c,d,m,n)
  call mktau(wrk,wrksize,v1,t11,t12,Half,rc)

  !3.1 map V3(f,m,n,e) <= V2(e,f,m,n)
  call map(wrk,wrksize,4,4,1,2,3,v2,1,v3,post,rc)

end if

!par
if (myRank == idbaab) then

  !3.2 mult M1(a,e) <= V1(a,f,m,n) . V3(f,m,n,e)
  call ccmult(wrk,wrksize,4,4,2,3,v1,1,v3,1,m1,ssc,rc)

  !3.3 add f1(a,e))aa <- - M1(a,e)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,-One,m1,1,f11,1,rc)
end if

!par
if (myRank == idaabb) then

  !4.1 map V4(a,f,m,n) <= V1(f,a,m,n)
  call map(wrk,wrksize,4,2,1,3,4,v1,1,v4,post,rc)

  !4.2 map V3(f,m,n,e) <= V2 (f,e,m,n)
  call map(wrk,wrksize,4,1,4,2,3,v2,1,v3,post,rc)

  !4.3 mult M1(a,e) <= V4(a,f,m,n) . V3(f,m,n,e)
  call ccmult(wrk,wrksize,4,4,2,3,v4,1,v3,1,m1,ssc,rc)

  !4.4 add f1(a,e))bb <- - M1(a,e)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,-One,m1,1,f12,1,rc)

end if

return

end subroutine contf14

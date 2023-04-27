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

subroutine contf32(wrk,wrksize,lunabij1,lunabij2,lunabij3)
! this routine does:
! FIII2 f3(e,m) <- sum(n,f) [ <ef||mn> . T1o(n,f) ]
!
! N.B use and destroy : V1,V2,M1
! N.B # of read : 3

use ccsd_global, only: f31, f32, idaabb, idbaab, m1, t11, t12, v1, v2
use Para_Info, only: MyRank
use Constants, only: One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize
real(kind=wp), intent(inout) :: wrk(wrksize)
integer(kind=iwp), intent(_IN_) :: lunabij1, lunabij2, lunabij3
integer(kind=iwp) :: post, rc, ssc

!1 f3(e,m)aa <- sum(n,f-aa) [ <ef||mn>aaaa . T1o(n,f)aa ]

!par
if (myRank == idbaab) then

  !1.1 read V1(ef,mn) <= <ef||mn>aaaa
  call filemanager(2,lunabij1,rc)
  call getmediate(wrk,wrksize,lunabij1,v1,rc)

  !1.2 expand V2(e,f,m,n) <= V1(ef,mn)
  call expand(wrk,wrksize,4,4,v1,1,v2,rc)

  !1.3 map V1(e,m,f,n) <= V2(e,f,m,n)
  call map(wrk,wrksize,4,1,3,2,4,v2,1,v1,post,rc)

  !1.4 mult M1(e,m) = V1(e,m,f,n) . T1o(f,n)aa
  call ccmult(wrk,wrksize,4,2,2,2,v1,1,t11,1,m1,ssc,rc)

  !1.5 add f3(e,m)aa <- M1(e,m)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,One,m1,1,f31,1,rc)

end if

!2 f3(e,m)bb <- sum(n,f-bb) [ <ef||mn>bbbb . T1o(n,f)bb ]

!par
if (myRank == idaabb) then

  !2.1 read V1(ef,mn) <= <ef||mn>bbbb
  call filemanager(2,lunabij2,rc)
  call getmediate(wrk,wrksize,lunabij2,v1,rc)

  !2.2 expand V2(e,f,m,n) <= V1(ef,mn)
  call expand(wrk,wrksize,4,4,v1,1,v2,rc)

  !2.3 map V1(e,m,f,n) <= V2(e,f,m,n)
  call map(wrk,wrksize,4,1,3,2,4,v2,1,v1,post,rc)

  !2.4 mult M1(e,m) = V1(e,m,f,n) . T1o(f,n)bb
  call ccmult(wrk,wrksize,4,2,2,2,v1,1,t12,1,m1,ssc,rc)

  !2.5 add f3(e,m)bb <- M1(e,m)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,One,m1,1,f32,1,rc)

end if

!3 f3(e,m)aa <- sum(n,f-bb) [ <ef||mn>abab . T1o(n,f)bb ]
!4 f3(e,m)bb <- sum(n,f-aa) [ <fe||nm>abab . T1o(n,f)aa ]

!par
if ((myRank == idbaab) .or. (myRank == idaabb)) then

  !34.1 read V1(c,d,k,l) <= <cd||kl>abab
  call filemanager(2,lunabij3,rc)
  call getmediate(wrk,wrksize,lunabij3,v1,rc)

end if

!par
if (myRank == idbaab) then

  !3.1 map V2(e,m,f,n) <= V1 (e,f,m,n)
  call map(wrk,wrksize,4,1,3,2,4,v1,1,v2,post,rc)

  !3.2 mult M1(e,m) = V2(e,m,f,n) . T1o(f,n)bb
  call ccmult(wrk,wrksize,4,2,2,2,v2,1,t12,1,m1,ssc,rc)

  !3.3 add f3(e,m)aa <- M1(e,m)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,One,m1,1,f31,1,rc)

end if

!par
if (myRank == idaabb) then

  !4.1 map V2(e,m,f,n) <= V1 (f,e,n,m)
  call map(wrk,wrksize,4,3,1,4,2,v1,1,v2,post,rc)

  !4.2 mult M1(e,m) = V2(e,m,f,n) . T1o(f,n)aa
  call ccmult(wrk,wrksize,4,2,2,2,v2,1,t11,1,m1,ssc,rc)

  !4.3 add f3(e,m)bb <- M1(e,m)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,One,m1,1,f32,1,rc)

end if

return

end subroutine contf32

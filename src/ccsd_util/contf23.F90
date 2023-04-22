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

subroutine contf23(wrk,wrksize)
! this routine does:
! f2(m,i) <- sum(e,n) [ <ie||mn> . T1o(e,n)]

use ccsd_global, only: f21, f22, idaabb, idbaab, m1, t11, t12, v1, v2, w11, w12, w13, w14
use Para_Info, only: MyRank
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize
real(kind=wp), intent(inout) :: wrk(wrksize)
integer(kind=iwp) :: post, rc, ssc

!1 f2(m,i)aa <- sum(e,n-aa) [ <ie||mn>aaaa . t1o(e,n)aa ]

!par
if (myRank == idbaab) then

  !1.1 expand V1(i,e,m,n) <= <ie||mn>aaaa
  call expand(wrk,wrksize,4,3,w11,1,v1,rc)

  !1.2 map V2(m,i,e,n) <= V1(i,e,m,n)
  call map(wrk,wrksize,4,2,3,1,4,v1,1,v2,post,rc)

  !1.3 mult M1(m,i) = V2(m,i,e,n) . T2o(e,n)aa
  call ccmult(wrk,wrksize,4,2,2,2,v2,1,t11,1,m1,ssc,rc)

  !1.4 add f2(m,i)aa <- M1(m,i)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,One,m1,1,f21,1,rc)

  !2 f2(m,i)aa <- sum(e,n-bb) [ <ie||mn>abab . t1o(e,n)bb ]

  !2.1 map V2(m,i,e,n) <= <ie||mn>abab
  call map(wrk,wrksize,4,2,3,1,4,w13,1,v2,post,rc)

  !2.2 mult M1(m,i) = V2(m,i,e,n) . T2o(e,n)bb
  call ccmult(wrk,wrksize,4,2,2,2,v2,1,t12,1,m1,ssc,rc)

  !2.3 add f2(m,i)aa <- M1(m,i)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,One,m1,1,f21,1,rc)

end if

!3 f2(m,i)bb <- sum(e,n-bb) [ <ie||mn>bbbb . t1o(e,n)bb ]

!par
if (myRank == idaabb) then

  !3.1 expand V1(i,e,m,n) <= <ie||mn>bbbb
  call expand(wrk,wrksize,4,3,w12,1,v1,rc)

  !3.2 map V2(m,i,e,n) <= V1(i,e,m,n)
  call map(wrk,wrksize,4,2,3,1,4,v1,1,v2,post,rc)

  !3.3 mult M1(m,i) =  V2(m,i,e,n) . T2o(e,n)bb
  call ccmult(wrk,wrksize,4,2,2,2,v2,1,t12,1,m1,ssc,rc)

  !3.4 add f2(m,i)bb <- M1(m,i)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,One,m1,1,f22,1,rc)

  !4 f2(m,i)bb <- - sum(e,n-aa) [ <ie||nm>baab . t1o(e,n)aa ]

  !4.1 map V2(m,i,e,n) <= <ie||nm>baab
  call map(wrk,wrksize,4,2,3,4,1,w14,1,v2,post,rc)

  !4.2 mult M1(m,i) = V2(m,i,e,n) . T2o(e,n)aa
  call ccmult(wrk,wrksize,4,2,2,2,v2,1,t11,1,m1,ssc,rc)

  !4.3 add f2(m,i)aa <- M1(m,i)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,-One,m1,1,f22,1,rc)

end if

return

end subroutine contf23

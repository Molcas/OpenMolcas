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

subroutine contt147(wrk,wrksize,lunt2o1,lunt2o2,lunt2o3)
! this routine does contributions T14 and T17
! T14: t1n(a,i) <- sum(me) [ T2o(ae,im) . FIII(e,m)]
! T17: t1n(a,i) <- sum(e,m>n) [ T2o(ae,mn) . <ie||mn> ]
!
! 1. T1n(a,i)aa <- sum(m,e-aa) [ T2o(a,e,i,m)aaaa . FIII(e,m)aa ]
! 2. T1n(a,i)aa <- sum(m,e-bb) [ T2o(a,e,i,m)abab . FIII(e,m)bb ]
! 3. T1n(a,i)bb <- sum(m,e-bb) [ T2o(a,e,i,m)bbbb . FIII(e,m)bb ]
! 4. T1n(a,i)bb <- sum(m,e-aa) [ T2o(e,a,m,i)abab . FIII(e,m)aa ]
! 5. T1n(a,i)aa <- - sum(e,m>n-aaa) [ T2o(a,e,mn)aaaa  . <mn||ie>aaaa ]
! 6. T1n(a,i)aa <- - sum(e,m,n-bab) [ T2o(a,e,m,n)abab . <mn||ie>abab ]
! 7. T1n(a,i)bb <- - sum(e,m>n-bbb) [ T2o(a,e,mn)bbbb  . <mn||ie>bbbb ]
! 8. T1n(a,i)bb <- + sum(e,m,n-aab) [ T2o(e,a,m,n)abab . <mn||ie>abba ]
!
! N.B. use and destroy : V1,V2,V3,M1
! N.B. # of read : 3

use ccsd_global, only: f31, f32, idaabb, idbaab, idfin, m1, t13, t14, v1, v2, v3, w11, w12, w13, w14
use Para_Info, only: MyRank
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: wrksize, lunt2o1, lunt2o2, lunt2o3
real(kind=wp) :: wrk(wrksize)
integer(kind=iwp) :: posst, rc, ssc

!par
if ((myRank == idbaab) .or. (myRank == idaabb) .or. (myRank == idfin)) then
  !15.1 read V1(cd,kl) <= T2o(cd,kl)aaaa
  call filemanager(2,lunt2o1,rc)
  call getmediate(wrk,wrksize,lunt2o1,v1%pos0,v1%d,v1%i,rc)
end if

!1 T1n(a,i)aa <- sum(m,e-aa) [ T2o(a,e,i,m)aaaa . FIII(e,m)aa ]

!par
if (myRank == idbaab) then

  !1.1 expand V2(a,e,i,m) <= V1(ae,im)
  call expand(wrk,wrksize,4,4,v1%d,v1%i,1,v2%pos0,v2%d,v2%i,rc)

  !1.2 map V3(a,i,e,m) <= V2(a,e,i,m)
  call map(wrk,wrksize,4,1,3,2,4,v2%d,v2%i,1,v3%d,v3%i,v3%pos0,posst,rc)

  !1.3 mult M1(a,i) <= V3(a,i,e,m) . FIII(e,m)aa
  call mult(wrk,wrksize,4,2,2,2,v3%d,v3%i,1,f31%d,f31%i,1,m1%d,m1%i,ssc,m1%pos0,rc)

  !1.4 add t1n(a,i)aa <- M1(a,i)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,One,m1%d,1,t13%d,t13%i,1,rc)

end if

!5 T1n(a,i)aa <-  - sum(e,m>n-aaa) [ T2o(a,e,mn)aaaa  . <mn||ie>aaaa ]

!par
if (myRank == idfin) then

  !5.1 expand V2(a,e,mn) <= V1(ae,mn)
  call expand(wrk,wrksize,4,5,v1%d,v1%i,1,v2%pos0,v2%d,v2%i,rc)

  !5.2 map V3(e,mn,i) <= <ie||mn>aaaa
  call map(wrk,wrksize,4,4,1,2,3,w11%d,w11%i,1,v3%d,v3%i,v3%pos0,posst,rc)

  !5.3 mult M1(a,i) <= V2(a,e,mn) . V3(e,mn,i)
  call mult(wrk,wrksize,4,4,2,3,v2%d,v2%i,1,v3%d,v3%i,1,m1%d,m1%i,ssc,m1%pos0,rc)

  !5.4 add t1n(a,i)aa <-  - M1(a,i)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,-One,m1%d,1,t13%d,t13%i,1,rc)

end if

!par
if ((myRank == idbaab) .or. (myRank == idaabb) .or. (myRank == idfin)) then

  !37.1 read V1(cd,kl) <= T2o(cd,kl)bbbb
  call filemanager(2,lunt2o2,rc)
  call getmediate(wrk,wrksize,lunt2o2,v1%pos0,v1%d,v1%i,rc)

end if

!3 T1n(a,i)bb <- sum(m,e-bb) [ T2o(a,e,i,m)bbbb . FIII(e,m)bb ]

!par
if (myRank == idaabb) then

  !3.1 expand V2(a,e,i,m) <= V1(ae,im)
  call expand(wrk,wrksize,4,4,v1%d,v1%i,1,v2%pos0,v2%d,v2%i,rc)

  !3.2 map V3(a,i,e,m) <= V2(a,e,i,m)
  call map(wrk,wrksize,4,1,3,2,4,v2%d,v2%i,1,v3%d,v3%i,v3%pos0,posst,rc)

  !3.3 mult M1(a,i) <= V3(a,i,e,m) . FIII(e,m)bb
  call mult(wrk,wrksize,4,2,2,2,v3%d,v3%i,1,f32%d,f32%i,1,m1%d,m1%i,ssc,m1%pos0,rc)

  !3.4 add t1n(a,i)bb <- M1(a,i)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,One,m1%d,1,t14%d,t14%i,1,rc)

end if

!7 T1n(a,i)bb <-  - sum(e,m>n-bbb) [ T2o(a,e,mn)bbbb  . <mn||ie>bbbb ]

!par
if (myRank == idfin) then

  !7.1 expand V2(a,e,mn) <= V1(ae,mn)
  call expand(wrk,wrksize,4,5,v1%d,v1%i,1,v2%pos0,v2%d,v2%i,rc)

  !7.2 map V3(e,mn,i) <= <ie||mn>bbbb
  call map(wrk,wrksize,4,4,1,2,3,w12%d,w12%i,1,v3%d,v3%i,v3%pos0,posst,rc)

  !7.3 mult M1(a,i) <= V2(a,e,mn) . V3(e,mn,i)
  call mult(wrk,wrksize,4,4,2,3,v2%d,v2%i,1,v3%d,v3%i,1,m1%d,m1%i,ssc,m1%pos0,rc)

  !7.4 add t1n(a,i)bb <- - M1(a,i)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,-One,m1%d,1,t14%d,t14%i,1,rc)

end if

!par
if ((myRank == idbaab) .or. (myRank == idaabb) .or. (myRank == idfin)) then

  !2468.1 read V1(c,d,k,l) <= T2o(c,d,k,l)abab
  call filemanager(2,lunt2o3,rc)
  call getmediate(wrk,wrksize,lunt2o3,v1%pos0,v1%d,v1%i,rc)

end if

!2 T1n(a,i)aa <- sum(m,e-bb) [ T2o(a,e,i,m)abab . FIII(e,m)bb ]

!par
if (myRank == idaabb) then

  !2.1 map V3(a,i,e,m) <= V1(a,e,i,m)
  call map(wrk,wrksize,4,1,3,2,4,v1%d,v1%i,1,v3%d,v3%i,v3%pos0,posst,rc)

  !2.2 mult M1(a,i) <= V3(a,i,e,m) . FIII(e,m)bb
  call mult(wrk,wrksize,4,2,2,2,v3%d,v3%i,1,f32%d,f32%i,1,m1%d,m1%i,ssc,m1%pos0,rc)

  !2.3 add t1n(a,i)aa <- M1(a,i)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,One,m1%d,1,t13%d,t13%i,1,rc)

end if

!4 T1n(a,i)bb <- sum(m,e-aa) [ T2o(e,a,m,i)abab . FIII(e,m)aa ]

!par
if (myRank == idbaab) then

  !4.1 map V3(a,i,e,m) <= V1(e,a,m,i)
  call map(wrk,wrksize,4,3,1,4,2,v1%d,v1%i,1,v3%d,v3%i,v3%pos0,posst,rc)

  !4.2 mult M1(a,i) <= V3(a,i,e,m) . FIII(e,m)aa
  call mult(wrk,wrksize,4,2,2,2,v3%d,v3%i,1,f31%d,f31%i,1,m1%d,m1%i,ssc,m1%pos0,rc)

  !4.3 add t1n(a,i)bb <- M1(a,i)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,One,m1%d,1,t14%d,t14%i,1,rc)

end if

!par
if (myRank == idfin) then

  !6 T1n(a,i)aa <-  - sum(e,m,n-bab) [ T2o(a,e,m,n)abab . <mn||ie>abab ]

  !6.1 map V3(e,m,n,i) <= <ie||mn>abab
  call map(wrk,wrksize,4,4,1,2,3,w13%d,w13%i,1,v3%d,v3%i,v3%pos0,posst,rc)

  !6.2 mult M1(a,i) <= V1(a,e,m,n) . V3(e,m,n,i)
  call mult(wrk,wrksize,4,4,2,3,v1%d,v1%i,1,v3%d,v3%i,1,m1%d,m1%i,ssc,m1%pos0,rc)

  !6.3 add t1n(a,i)aa <- - M1(a,i)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,-One,m1%d,1,t13%d,t13%i,1,rc)

  !8 T1n(a,i)bb <- + sum(e,m,n-aab) [ T2o(e,a,m,n)abab . <mn||ie>abba ]

  !8.1 map V2(a,e,m,n) <= V1(e,a,m,n)
  call map(wrk,wrksize,4,2,1,3,4,v1%d,v1%i,1,v2%d,v2%i,v2%pos0,posst,rc)

  !8.2 map V3(e,m,n,i) <= <ie||mn>baab
  call map(wrk,wrksize,4,4,1,2,3,w14%d,w14%i,1,v3%d,v3%i,v3%pos0,posst,rc)

  !8.3 mult M1(a,i) <= V2(a,e,m,n) . V3(e,m,n,i)
  call mult(wrk,wrksize,4,4,2,3,v2%d,v2%i,1,v3%d,v3%i,1,m1%d,m1%i,ssc,m1%pos0,rc)

  !8.4 add t1n(a,i)bb <-  M1(a,i)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,One,m1%d,1,t14%d,t14%i,1,rc)

end if

return

end subroutine contt147

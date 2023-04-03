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

subroutine contf5(wrk,wrksize,lunt2o1,lunt2o2,lunt2o3)
! this routine does: FV1,FV2,T23
!
!     FV1
!1    FV(m,j)aa  <= FII(m,j)aa
!2    FV(m,j)bb  <= FII(m,j)bb
!
!     FV2
!3    FV(m,j)aa  <- -0.5 sum (e-a) [ FIII(e,m)aa . T1o(e,j)aa ]
!4    FV(m,j)bb  <- -0.5 sum (e-b) [ FIII(e,m)bb . T1o(e,j)bb ]
!
!     T23
!5    Q(ab,k,l)aaaa   <= sum(m-a)   [ T2o(ab,k,m)aaaa . FV(m,l)aa ]
!5    T2n(ab,ij)aaaa   <- Q(ab,j,i)aaaa - Q(ab,i,j)aaaa
!6    Q(ab,k,l)bbbb   <= sum(m-b)   [ T2o(ab,k,m)bbbb . FV(m,l)bb ]
!6    T2n(ab,ij)bbbb   <- Q(ab,j,i)bbbb - Q(ab,i,j)bbbb
!7    T2n(a,b,i,j)abab <- - sum(m-b) [ T2o(a,b,i,m)abab . FV(m,j)bb ]
!8    T2n(a,b,i,j)abab <- - sum(m-a) [ T2o(a,b,m,j)abab . FV(m,i)aa ]

use Para_Info, only: MyRank
use Constants, only: One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: wrksize, lunt2o1, lunt2o2, lunt2o3
real(kind=wp) :: wrk(wrksize)
#include "ccsd2.fh"
#include "parallel.fh"
integer(kind=iwp) :: posst, rc, ssc

!par
if ((myRank == idbaab) .or. (myRank == idaabb)) then

  !78.0 read V1(a,b,k,l) <= T2o(a,b,k,l)abab
  call filemanager(2,lunt2o3,rc)
  call getmediate(wrk,wrksize,lunt2o3,possv10,mapdv1,mapiv1,rc)

end if

!par
if (myRank == idbaab) then

  !1 FV(m,j)aa  <= FII(m,j)aa

  !1.1 map M1(m,j) <= F2(m,j)aa
  call map(wrk,wrksize,2,1,2,0,0,mapdf21,mapif21,1,mapdm1,mapim1,possm10,posst,rc)

  !3 FV(m,j)aa  <- 0.5 sum (e-a) [ FIII(e,m)aa . T1o(e,j)aa ]

  !3.1 map M2(m,e) <= F3(e,m)aa
  call map(wrk,wrksize,2,2,1,0,0,mapdf31,mapif31,1,mapdm2,mapim2,possm20,posst,rc)

  !3.2 mult M3(m,j) <= M2(m,e) . T1o(e,j)aa
  call mult(wrk,wrksize,2,2,2,1,mapdm2,mapim2,1,mapdt11,mapit11,1,mapdm3,mapim3,ssc,possm30,rc)

  !3.3 add M1(m,j) <- 0.5 M3(m,j)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,Half,mapdm3,1,mapdm1,mapim1,1,rc)

  !5 Q(ab,k,l)aaaa   <= sum(m-a)   [ T2o(ab,k,m)aaaa . FV(m,l)aa ]
  !5 T2n(ab,ij)aaaa   <- Q(ab,j,i)aaaa - Q(ab,i,j)aaaa

  !5.1 read V3(ab,kl) <= T2o(ab,kl)aaaa
  call filemanager(2,lunt2o1,rc)
  call getmediate(wrk,wrksize,lunt2o1,possv30,mapdv3,mapiv3,rc)

  !5.2 expand V2(ab,k,l) <- V3(ab,kl)
  call expand(wrk,wrksize,4,6,mapdv3,mapiv3,1,possv20,mapdv2,mapiv2,rc)

  !5.3 mult V3(ab,k,l) <= V2(a,b,k,m) . M1(m,l)
  call mult(wrk,wrksize,4,2,4,1,mapdv2,mapiv2,1,mapdm1,mapim1,1,mapdv3,mapiv3,ssc,possv30,rc)

  !5.4 pack V2(ab,ij) <= V3(ab,i,j) - V3(ab,j,i)
  call fack(wrk,wrksize,4,4,mapdv3,1,mapiv3,mapdv2,mapiv2,possv20,rc)

  !5.5 add T2n(ab,ij)aaaa <- - V2(ab,ij)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,-One,mapdv2,1,mapdt21,mapit21,1,rc)

  !8 T2n(a,b,i,j)abab <- - sum(m-a) [ T2o(a,b,m,j)abab . FV(m,i)aa ]

  !8.1 map V2(a,b,j,m) <= V1(a,b,m,j)
  call map(wrk,wrksize,4,1,2,4,3,mapdv1,mapiv1,1,mapdv2,mapiv2,possv20,posst,rc)

  !8.2 mult V3(a,b,j,i) <= V2(a,b,j,m) . M1(m,i)
  call mult(wrk,wrksize,4,2,4,1,mapdv2,mapiv2,1,mapdm1,mapim1,1,mapdv3,mapiv3,ssc,possv30,rc)

  !8.3 map V2(a,b,i,j) <= V3(a,b,j,i)
  call map(wrk,wrksize,4,1,2,4,3,mapdv3,mapiv3,1,mapdv2,mapiv2,possv20,posst,rc)

  !8.4 add T2n(a,b,i,j)abab <- - V2(a,b,i,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,-One,mapdv2,1,mapdt23,mapit23,1,rc)
  !par
end if

!par
if (myRank == idaabb) then

  !2 FV(m,j)bb  <= FII(m,j)bb

  !2.1 map M1(m,j) <= F2(m,j)bb
  call map(wrk,wrksize,2,1,2,0,0,mapdf22,mapif22,1,mapdm1,mapim1,possm10,posst,rc)

  !4 FV(m,j)bb  <- 0.5 sum (e-b) [ FIII(e,m)bb . T1o(e,j)bb ]

  !4.1 map M2(m,e) <= F3(e,m)bb
  call map(wrk,wrksize,2,2,1,0,0,mapdf32,mapif32,1,mapdm2,mapim2,possm20,posst,rc)

  !4.2 mult M3(m,j) <= M2(m,e) . T1o(e,j)bb
  call mult(wrk,wrksize,2,2,2,1,mapdm2,mapim2,1,mapdt12,mapit12,1,mapdm3,mapim3,ssc,possm30,rc)

  !4.3 add M1(m,j) <- 0.5 M3(m,j)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,Half,mapdm3,1,mapdm1,mapim1,1,rc)

  !6 Q(ab,k,l)bbbb   <= sum(m-b)   [ T2o(ab,k,m)bbbb . FV(m,l)bb ]
  !6 T2n(ab,ij)bbbb   <- Q(ab,j,i)bbbb - Q(ab,i,j)bbbb

  !6.1 read V3(ab,kl) <= T2o(ab,kl)bbbb
  call filemanager(2,lunt2o2,rc)
  call getmediate(wrk,wrksize,lunt2o2,possv30,mapdv3,mapiv3,rc)

  !6.2 expand V2(ab,k,l) <- V3(ab,kl)
  call expand(wrk,wrksize,4,6,mapdv3,mapiv3,1,possv20,mapdv2,mapiv2,rc)

  !6.3 mult V3(ab,k,l) <= V2(a,b,k,m) . M1(m,l)
  call mult(wrk,wrksize,4,2,4,1,mapdv2,mapiv2,1,mapdm1,mapim1,1,mapdv3,mapiv3,ssc,possv30,rc)

  !6.4 pack V2(ab,ij) <= V3(ab,i,j) - V3(ab,j,i)
  call fack(wrk,wrksize,4,4,mapdv3,1,mapiv3,mapdv2,mapiv2,possv20,rc)

  !6.5 add T2n(ab,ij)bbbb <- - V2(ab,ij)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,-One,mapdv2,1,mapdt22,mapit22,1,rc)

  !7 T2n(a,b,i,j)abab <- - sum(m-b) [ T2o(a,b,i,m)abab . FV(m,j)bb ]

  !7.1 mult V3(a,b,i,j) <= V1(a,b,i,m) . M1(m,j)
  call mult(wrk,wrksize,4,2,4,1,mapdv1,mapiv1,1,mapdm1,mapim1,1,mapdv3,mapiv3,ssc,possv30,rc)

  !7.2 add T2n(a,b,i,j)abab <- - V3(a,b,i,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,-One,mapdv3,1,mapdt23,mapit23,1,rc)
  !par
end if

return

end subroutine contf5

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

subroutine contf4(wrk,wrksize,lunt2o1,lunt2o2,lunt2o3)
! this routine does: FIV1,FIV2,T22

!     FIV1
!1    FIV(b,e)aa <= FI(b,e)aa
!2    FIV(b,e)bb <= FI(b,e)bb
!
!     FIV2
!3    FIV(b,e)aa <- -0.5 sum(m-a) [ T1o(b,m)aa . FIII(e,m)aa ]
!4    FIV(b,e)bb <- -0.5 sum(m-b) [ T1o(b,m)bb . FIII(e,m)bb ]
!
!     T22
!5    Q(c,d,ij)aaaa   <- - sum(e-a) [ FIV(c,e)aa . T2o(e,d,ij)aaaa ]
!5    T2n(ab,ij)aaaa   <- Q(b,a,ij)aaaa - Q(a,b,ij)aaaa
!6    Q(c,d,ij)bbbb   <- - sum(e-b) [ FIV(c,e)bb . T2o(e,d,ij)bbbb ]
!6    T2n(ab,ij)bbbb   <- Q(b,a,ij)bbbb - Q(a,b,ij)bbbb
!7    T2n(a,b,i,j)abab <- sum (e-a)  [ FIV(a,e)aa . T2o(e,b,i,j)abab ]
!8    T2n(a,b,i,j)abab <- sum (e-b)  [ FIV(b,e)bb . T2o(a,e,i,j)abab ]
!
! N.B. use and destroy : V1,V2,V3,M1,M2,M3
! N.B. # of read       : 3
! N.B. Parallel : in the case where idaaaa /= idbaab and
!      idaaaa /= idbaab this routine runs contributions to
!      T2n also on idaaaa and idbbbb nodes, but only with
!      FIV = FI, since on these nodes there is a specific
!      part of contributions F13 (see notes in sumoverb routine)
!      which are not presented on pilot nodes

use ccsd_global, only: f11, f12, f31, f32, idaaaa, idaabb, idbaab, idbbbb, m1, m2, m3, t11, t12, t21, t22, t23, v1, v2, v3
use Para_Info, only: MyRank
use Constants, only: One, Half
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize
real(kind=wp), intent(inout) :: wrk(wrksize)
integer(kind=iwp), intent(_IN_) :: lunt2o1, lunt2o2, lunt2o3
integer(kind=iwp) :: post, rc, ssc

!par
if ((myRank == idbaab) .or. (myRank == idaabb) .or. (myRank == idaaaa) .or. (myRank == idbbbb)) then
  !78.0 read V1(c,d,i,j) <= T2o(c,d,i,j)abab
  call filemanager(2,lunt2o3,rc)
  call getmediate(wrk,wrksize,lunt2o3,v1,rc)
end if

!1 FIV(b,e)aa <= FI(b,e)aa

!par
if ((myRank == idbaab) .or. (myRank == idaaaa)) then

  !1.1 map M1(b,e) <= F1(b,e)aa
  call map(wrk,wrksize,2,1,2,0,0,f11,1,m1,post,rc)

end if

!3 FIV(b,e)aa <- -0.5 sum(m-a) [ T1o(b,m)aa . FIII(e,m)aa ]

!par
if (myRank == idbaab) then

  !3.1 map M2(m,e) <= f3(e,m)aa
  call map(wrk,wrksize,2,2,1,0,0,f31,1,m2,post,rc)

  !3.2 mult M3(b,e) <= T1o(b,m)aa . M2(m,e)
  call ccmult(wrk,wrksize,2,2,2,1,t11,1,m2,1,m3,ssc,rc)

  !3.3 add M1(b,e) <- -0.5 M3(b,e)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,-Half,m3,1,m1,1,rc)

end if
!parend

!par
if ((myRank == idbaab) .or. (myRank == idaaaa)) then

  !5 Q(c,d,ij)aaaa   <- - sum(e-a) [ FIV(c,e)aa . T2o(e,d,ij)aaaa ]
  !5 T2n(ab,ij)aaaa   <- Q(b,a,ij)aaaa - Q(a,b,ij)aaaa

  !5.1 read V2(ed,ij) <= T2o(ed,ij)aaaa
  call filemanager(2,lunt2o1,rc)
  call getmediate(wrk,wrksize,lunt2o1,v2,rc)

  !5.2 expand V3(e,d,ij) <= V2 (ed,ij)
  call expand(wrk,wrksize,4,5,v2,1,v3,rc)

  !5.3 mult V2(c,d,ij) <= M1(c,e) . V3(e,d,ij)
  call ccmult(wrk,wrksize,2,4,4,1,m1,1,v3,1,v2,ssc,rc)

  !5.4 pack V3(ab,ij) <= V2(a,b,ij) - V2(b,a,ij)
  call fack(wrk,wrksize,4,4,v2,1,v3,rc)

  !5.5 add T2n(ab,ij)aaaa <- V3(ab,ij)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,One,v3,1,t21,1,rc)

  !7 T2n(a,b,i,j)abab <- sum (e-a)  [ FIV(a,e)aa . T2o(e,b,i,j)abab ]

  !7.1 mult V2(a,b,i,j) <= M1(a,e) . V1(e,b,i,j)
  call ccmult(wrk,wrksize,2,4,4,1,m1,1,v1,1,v2,ssc,rc)

  !7.2 add T2n(a,b,i,j)abab <- V2(a,b,i,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,One,v2,1,t23,1,rc)

end if
!parend

!2 FIV(b,e)bb <= FI(b,e)bb

!par
if ((myRank == idaabb) .or. (myRank == idbbbb)) then

  !2.1 map M1(b,e) <= F1(b,e)bb
  call map(wrk,wrksize,2,1,2,0,0,f12,1,m1,post,rc)

end if

!4 FIV(b,e)bb <- -0.5 sum(m-b) [ T1o(b,m)bb . FIII(e,m)bb ]

!par
if (myRank == idaabb) then

  !4.1 map M2(m,e) <= f3(e,m)bb
  call map(wrk,wrksize,2,2,1,0,0,f32,1,m2,post,rc)

  !4.2 mult M3(b,e) <= T1o(b,m)bb . M2(m,e)
  call ccmult(wrk,wrksize,2,2,2,1,t12,1,m2,1,m3,ssc,rc)

  !4.3 add M1(b,e) <- -0.5 M3(b,e)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,-Half,m3,1,m1,1,rc)

end if
!parend

!par
if ((myRank == idaabb) .or. (myRank == idbbbb)) then

  !6 Q(c,d,ij)bbbb   <- - sum(e-b) [ FIV(c,e)bb . T2o(e,d,ij)bbbb ]
  !6 T2n(ab,ij)bbbb   <- Q(b,a,ij)bbbb - Q(a,b,ij)bbbb

  !6.1 read V2(cd,ij) <= T2o(cd,ij)bbbb
  call filemanager(2,lunt2o2,rc)
  call getmediate(wrk,wrksize,lunt2o2,v2,rc)

  !6.2 expand V3(e,d,ij) <= V2 (cd,ij)
  call expand(wrk,wrksize,4,5,v2,1,v3,rc)

  !6.3 mult V2(c,d,ij) <= M1(c,d) . V3(e,d,ij)
  call ccmult(wrk,wrksize,2,4,4,1,m1,1,v3,1,v2,ssc,rc)

  !6.4 pack V3(ab,ij) <= V2(a,b,ij) - V2(b,a,ij)
  call fack(wrk,wrksize,4,4,v2,1,v3,rc)

  !6.5 add T2n(ab,ij)bbbb <- V3(ab,ij)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,One,v3,1,t22,1,rc)

  !8 T2n(a,b,i,j)abab <- sum (e-b)  [ FIV(b,e)bb . T2o(a,e,i,j)abab ]

  !8.1 map V3(e,a,i,j) <= V1(a,e,i,j)
  call map(wrk,wrksize,4,2,1,3,4,v1,1,v3,post,rc)

  !8.2 mult V2(b,a,i,j) <= M1(b,e) . V3(e,a,i,j)
  call ccmult(wrk,wrksize,2,4,4,1,m1,1,v3,1,v2,ssc,rc)

  !8.3 map V3(a,b,i,j) <= V2(b,a,i,j)
  call map(wrk,wrksize,4,2,1,3,4,v2,1,v3,post,rc)

  !8.4 add T2n(a,b,i,j)abab <- V3(a,b,i,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,One,v3,1,t23,1,rc)

end if
!parend

return

end subroutine contf4

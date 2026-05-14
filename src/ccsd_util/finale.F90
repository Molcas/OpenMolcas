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

subroutine finale(wrk,wrksize,lunabij1,lunabij2,lunabij3,lunt2o1,lunt2o2,lunt2o3)
! this routine does:
! 1)  FI2   f1(a,e) <- -0.5 sum(m) [t1o(a,m) . fok(e,m)]
! 2)  FI4   f1(a,e) <- -sum(m>n,f) [Tap(af,mn) . <ef||mn>]
!
! 3)  FII2  f2(m,i) <- -0.5 sum(m) [fok(m,e) . T1o(e,i)]
! 4)  FII3  f2(m,i) <- sum(e,n) [ <ie||mn> . T1o(e,n)]
! 5)  FII4  f2(m,i) <- sum(n,e>f) [ <ef||mn> . Tap(in,ef)]
!
! 6)  FIII2 f3(e,m) <- sum(n,f) [ <ef||mn> . T1o(n,f) ]
!
! 7)  FIV1 FIV(b,e) <= FI(b,e)
! 8)  FIV2 FIV(b,e) <- -0.5 sum(m) [ T1o(b,m) . FIII(e,m)]
! 9)  T22 T2n(ab,ij) <- P(ab) sum(e) [t2o(ae,ij) . F4(b.e)]
!
! 10) FV1 FV(m,j) <= FII(m,j)
! 11) FV2 FV(m,j)aa  <- -0.5 sum (e) [ FIII(e,m) . T1o(e,j) ]
! 12) T23 T2n(ab,ij) <- - P(ij) sum(m) [T2o(ab,im) . F5(m,j) ]
!
! 13) WI1 w1(mnij) <= <mn||ij>
! 14) W12 w1(mnij) <- P(ij) sum(e) [ <ie||mn> . T1o(e,j) ]
! 15) WI3 w1(mnij) <- sum(e>f) [ <ef||mn> . Tau(ef,ij) ]
! 16) T24 t2n(ab,ij) <- sum(m>n) [ Tau(ab,mn) . W1(mn,ij) ]
!
! 17) T12 t1n(a,i) <- sum(e) [FI(a,e) . T1o(e,i)]
! 18) T13 t1n(a,i) <- - sum(m) [T1o(a,m) . FII(m,i)]
! 19) T14 t1n(a,i) <- sum(me) [ T2o(ae,im) . FIII(e,m)]
! 20) T17 t1n(a,i) <- sum(e,m>n) [ T2o(ae,mn) . <ie||mn> ]
!
! ..) T22 done as 9)
! ..) T23 done as 12)
! ..) T24 done as 16)
! 21) T29 t2n(ab,ij) <- -P(a,b) sum(m) [T1o(a,m) . <mb||ij>]

use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize
real(kind=wp), intent(inout) :: wrk(wrksize)
integer(kind=iwp), intent(_IN_) :: lunabij1, lunabij2, lunabij3, lunt2o1, lunt2o2, lunt2o3

!1
call contf12(wrk,wrksize)
!2
call contf14(wrk,wrksize,lunabij1,lunabij2,lunabij3,lunt2o1,lunt2o2,lunt2o3)
!3
call contf22(wrk,wrksize)
!4
call contf23(wrk,wrksize)
!5
call contf24(wrk,wrksize,lunabij1,lunabij2,lunabij3,lunt2o1,lunt2o2,lunt2o3)
!6
call contf32(wrk,wrksize,lunabij1,lunabij2,lunabij3)
!7,8,9
call contf4(wrk,wrksize,lunt2o1,lunt2o2,lunt2o3)
!10,11,12
call contf5(wrk,wrksize,lunt2o1,lunt2o2,lunt2o3)
!13,14,15,16
call contw1(wrk,wrksize,lunabij1,lunabij2,lunabij3,lunt2o1,lunt2o2,lunt2o3)
!17
call contt12(wrk,wrksize)
!18
call contt13(wrk,wrksize)
!19,20
call contt147(wrk,wrksize,lunt2o1,lunt2o2,lunt2o3)
!21
call contt29(wrk,wrksize)

return

end subroutine finale

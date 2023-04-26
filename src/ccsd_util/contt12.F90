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

subroutine contt12(wrk,wrksize)
! this routine does T12 contribution:
! t1n(a,i) <- sum(e) [FI(a,e) . T1o(e,i)
!
! N.B. use and destroy : M1
! N.B. Parallel : in the case where idaaaa /= idbaab and
!      idaaaa /= idbaab this routine runs contributions to
!      T1n also on idaaaa and idbbbb nodes, since on these
!      nodes there is a specific
!      part of contributions F13 (see notes in sumoverb routine)
!      which are not presented on pilot nodes

use ccsd_global, only: f11, f12, idaaaa, idaabb, idbaab, idbbbb, m1, t11, t12, t13, t14
use Para_Info, only: MyRank
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize
real(kind=wp), intent(inout) :: wrk(wrksize)
integer(kind=iwp) :: rc, ssc

!1 t1n(a,i)aa <- sum(e-a) [F1(a,e)aa . T1o(e,i)aa ]

!par
if ((myRank == idbaab) .or. (myRank == idaaaa)) then
  !1.1 mult M1(a,i) <= F1(a,e)aa . T1o(e,i)aa
  call ccmult(wrk,wrksize,2,2,2,1,f11,1,t11,1,m1,ssc,rc)

  !1.2 add t1n(a,i)aa <- M1(a,i)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,One,m1,1,t13,1,rc)
end if

!2 t1n(a,i)bb <- sum(e-a) [F1(a,e)bb . T1o(e,i)bb ]

!par
if ((myRank == idaabb) .or. (myRank == idbbbb)) then
  !2.1 mult M1(a,i) <= F1(a,e)bb . T1o(e,i)bb
  call ccmult(wrk,wrksize,2,2,2,1,f12,1,t12,1,m1,ssc,rc)

  !2.2 add t1n(a,i)bb <- M1(a,i)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,One,m1,1,t14,1,rc)
end if

return

end subroutine contt12

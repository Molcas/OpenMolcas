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

subroutine contt13(wrk,wrksize)
! this routine does T13 contribution:
! t1n(a,i) <- sum(m) [T1o(a,m) . FII(m,i)]
!
! N.B. use and destroy : M1

use ccsd_global, only: f21, f22, idaabb, idbaab, m1, t11, t12, t13, t14
use Para_Info, only: MyRank
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize
real(kind=wp), intent(inout) :: wrk(wrksize)
integer(kind=iwp) :: rc, ssc

!1 t1n(a,i)aa <- - sum(m-a) [T1o(a,m)aa  . FII(m,i)aa]

!par
if (myRank == idbaab) then
  !1.1 mult M1(a,i) <= T1o(a,m)aa  . FII(m,i)aa
  call ccmult(wrk,wrksize,2,2,2,1,t11,1,f21,1,m1,ssc,rc)

  !1.2 add t1n(a,i)aa <- M1(a,i)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,-One,m1,1,t13,1,rc)
end if

!2 t1n(a,i)bb <- - sum(m-b) [T1o(a,m)bb  . FII(m,i)bb]

!par
if (myRank == idaabb) then
  !2.1 mult M1(a,i) <= T1o(a,m)bb  . FII(m,i)bb
  call ccmult(wrk,wrksize,2,2,2,1,t12,1,f22,1,m1,ssc,rc)

  !2.2 add t1n(a,i)bb <- M1(a,i)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,-One,m1,1,t14,1,rc)
end if

return

end subroutine contt13

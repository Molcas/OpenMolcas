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

subroutine contf22(wrk,wrksize)
! this routine does:
! f2(m,i) <- 0.5 sum(m) [fok(m,e) . T1o(e,i)]
!
! N.B. use and destroy : M1,M2

use ccsd_global, only: f21, f22, fk3, fk4, idaabb, idbaab, m1, m2, t11, t12
use Para_Info, only: MyRank
use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize
real(kind=wp), intent(inout) :: wrk(wrksize)
integer(kind=iwp) :: post, rc, ssc

!1 f2(m,i)aa <- 0.5 sum(e-a) [ fok(e,m)aa . t1o(e,i)aa]

!par
if (myRank == idbaab) then

  !1.1 map M1(m,e) <= fok(e,m)aa
  call map(wrk,wrksize,2,2,1,0,0,fk3,1,m1,post,rc)

  !1.2 mult M2(m,i) <= M1(m,e) . T1o(e,i)aa
  call ccmult(wrk,wrksize,2,2,2,1,m1,1,t11,1,m2,ssc,rc)

  !1.3 add f2(m,i)aa <- 0.5 M2(m,i)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,Half,m2,1,f21,1,rc)

end if

!2 f2(m,i)aa <- 0.5 sum(e-b) [ fok(e,m)bb . t1o(e,i)bb]

!par
if (myRank == idaabb) then

  !2.1 map M1(m,e) <= fok(e,m)bb
  call map(wrk,wrksize,2,2,1,0,0,fk4,1,m1,post,rc)

  !2.2 mult M2(m,i) <= M1(m,e) . T1o(e,i)bb
  call ccmult(wrk,wrksize,2,2,2,1,m1,1,t12,1,m2,ssc,rc)

  !2.3 add f2(m,i)bb <- 0.5 M2(m,i)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,Half,m2,1,f22,1,rc)

end if

return

end subroutine contf22

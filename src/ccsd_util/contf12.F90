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

subroutine contf12(wrk,wrksize)
! this routine does:
! FI2   f1(a,e) <- -0.5 sum(m) [t1o(a,m) . fok(e,m)]
!
! N.B. use and destroy: M1,M2

use ccsd_global, only: f11, f12, fk3, fk4, idaabb, idbaab, m1, m2, t11, t12
use Para_Info, only: MyRank
use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: wrksize
real(kind=wp) :: wrk(wrksize)
integer(kind=iwp) :: posst, rc, ssc

!1 f1(a,e)aa <- sum(m-a) [T1o(a,m)aa . fok(e,m)aa]

!par
if (myRank == idbaab) then

  !1.1 map M1(m,e) <- fok(e,m)aa
  call map(wrk,wrksize,2,2,1,0,0,fk3%d,fk3%i,1,m1%d,m1%i,m1%pos0,posst,rc)
  !1.2 mult M2(a,e) = t1o(a,m)aa . M1(m,e)
  call mult(wrk,wrksize,2,2,2,1,t11%d,t11%i,1,m1%d,m1%i,1,m2%d,m2%i,ssc,m2%pos0,rc)
  !1.3 add f1(a,e)aa <- -0.5 M2(a,e)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,-Half,m2%d,1,f11%d,f11%i,1,rc)

end if

!2 f1(a,e)bb <- sum(m-b) [T1o(a,m)bb . fok(e,m)bb]

!par
if (myRank == idaabb) then

  !2.1 map M1(m,e) <- fok(e,m)bb
  call map(wrk,wrksize,2,2,1,0,0,fk4%d,fk4%i,1,m1%d,m1%i,m1%pos0,posst,rc)
  !2.2 mult M2(a,e) = t1o(a,m)bb . M1(m,e)
  call mult(wrk,wrksize,2,2,2,1,t12%d,t12%i,1,m1%d,m1%i,1,m2%d,m2%i,ssc,m2%pos0,rc)
  !2.3 add f1(a,e)bb <- -0.5 M2(a,e)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,-Half,m2%d,1,f12%d,f12%i,1,rc)

end if

return

end subroutine contf12

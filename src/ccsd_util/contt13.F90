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

use Para_Info, only: MyRank
#include "ccsd2.fh"
#include "parallel.fh"
#include "wrk.fh"
! help variables
integer rc, ssc

!1 t1n(a,i)aa <- - sum(m-a) [T1o(a,m)aa  . FII(m,i)aa]

!par
if (myRank == idbaab) then
  !1.1 mult M1(a,i) <= T1o(a,m)aa  . FII(m,i)aa
  call mult(wrk,wrksize,2,2,2,1,mapdt11,mapit11,1,mapdf21,mapif21,1,mapdm1,mapim1,ssc,possm10,rc)

  !1.2 add t1n(a,i)aa <- M1(a,i)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,-1.0d0,mapdm1,1,mapdt13,mapit13,1,rc)
end if

!2 t1n(a,i)bb <- - sum(m-b) [T1o(a,m)bb  . FII(m,i)bb]

!par
if (myRank == idaabb) then
  !2.1 mult M1(a,i) <= T1o(a,m)bb  . FII(m,i)bb
  call mult(wrk,wrksize,2,2,2,1,mapdt12,mapit12,1,mapdf22,mapif22,1,mapdm1,mapim1,ssc,possm10,rc)

  !2.2 add t1n(a,i)bb <- M1(a,i)
  call add(wrk,wrksize,2,2,0,0,0,0,1,1,-1.0d0,mapdm1,1,mapdt14,mapit14,1,rc)
end if

return

end subroutine contt13

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

subroutine divfok(wrk,wrksize,mapdfa,mapifa,possfa0,mapdfb,mapifb,possfb0,mapdfk1,mapifk1,possfk10,mapdfk2,mapifk2,possfk20, &
                  mapdfk3,mapifk3,possfk30,mapdfk4,mapifk4,possfk40,mapdfk5,mapifk5,possfk50,mapdfk6,mapifk6,possfk60,mapddp1, &
                  mapidp1,possdp10,mapddp2,mapidp2,possdp20,rc)
! this routine divides fok(p,q) -> fk(a,b) + fk(a,i) + f(i,j) + dp(p)
! to diagonal part and rest
!
! mapd and mapi for:
! fa,fb - fok(p,q)aa,bb
! fk1-6 - f(ab)aa,f(ab)bb,f(ai)aa,f(ai)bb,f(ij)aa,f(ij)bb
! dp1,2 - diagonal part dp(p)a,b
! rc    - return (error) code

use ccsd_global, only: noa, nob, norb, nsym, nva, nvb
use Definitions, only: wp, iwp

implicit none
!1 maps for FOKA,FOKB
!2 maps for FK
!  FK1 - f(a,b)aa
!  FK2 - f(a,b)bb
!  FK3 - f(a,i)aa
!  FK4 - f(a,i)bb
!  FK5 - f(i,j)aa
!  FK6 - f(i,j)bb
!3 maps for DP - diagonal part
!  DP1 - dp(p)a
!  DP2 - dp(p)b
integer(kind=iwp) :: wrksize, mapdfa(0:512,6), mapifa(8,8,8), possfa0, mapdfb(0:512,6), mapifb(8,8,8), possfb0, mapdfk1(0:512,6), &
                     mapifk1(8,8,8), possfk10, mapdfk2(0:512,6), mapifk2(8,8,8), possfk20, mapdfk3(0:512,6), mapifk3(8,8,8), &
                     possfk30, mapdfk4(0:512,6), mapifk4(8,8,8), possfk40, mapdfk5(0:512,6), mapifk5(8,8,8), possfk50, &
                     mapdfk6(0:512,6), mapifk6(8,8,8), possfk60, mapddp1(0:512,6), mapidp1(8,8,8), possdp10, mapddp2(0:512,6), &
                     mapidp2(8,8,8), possdp20, rc
real(kind=wp) :: wrk(wrksize)
integer(kind=iwp) :: iidp, iidpa, iidpb, iifaa, iifai, iifii, iifok, iifoka, iifokb, possdp, possdpa, possdpb, possfaa, possfai, &
                     possfii, possfok, possfoka, possfokb, rc1, symp

rc = 0

!1 define dp

do symp=1,nsym

  iidpa = mapidp1(symp,1,1)
  possdpa = mapddp1(iidpa,1)
  iidpb = mapidp2(symp,1,1)
  possdpb = mapddp2(iidpb,1)
  iifoka = mapifa(symp,1,1)
  possfoka = mapdfa(iifoka,1)
  iifokb = mapifb(symp,1,1)
  possfokb = mapdfb(iifokb,1)

  if (norb(symp) > 0) call fokunpck5(symp,wrk(possfoka),wrk(possfokb),wrk(possdpa),wrk(possdpb),norb(symp),rc1)

end do

!2 define faa,fai,fii

do symp=1,nsym
  if (norb(symp) == 0) cycle

  !2.1 alpha case

  iifok = mapifa(symp,1,1)
  iifaa = mapifk1(symp,1,1)
  iifai = mapifk3(symp,1,1)
  iifii = mapifk5(symp,1,1)
  iidp = mapidp1(symp,1,1)

  possfok = mapdfa(iifok,1)
  possfaa = mapdfk1(iifaa,1)
  possfai = mapdfk3(iifai,1)
  possfii = mapdfk5(iifii,1)
  possdp = mapddp1(iidp,1)

  call fokunpck1(wrk(possfok),wrk(possdp),norb(symp))
  if (nva(symp) > 0) call fokunpck2(wrk(possfok),wrk(possfaa),norb(symp),nva(symp),noa(symp))
  if ((noa(symp)*nva(symp)) > 0) call fokunpck3(wrk(possfok),wrk(possfai),norb(symp),nva(symp),noa(symp))
  if (noa(symp) > 0) call fokunpck4(wrk(possfok),wrk(possfii),norb(symp),noa(symp))

  !2.2 alpha case

  iifok = mapifb(symp,1,1)
  iifaa = mapifk2(symp,1,1)
  iifai = mapifk4(symp,1,1)
  iifii = mapifk6(symp,1,1)
  iidp = mapidp2(symp,1,1)

  possfok = mapdfb(iifok,1)
  possfaa = mapdfk2(iifaa,1)
  possfai = mapdfk4(iifai,1)
  possfii = mapdfk6(iifii,1)
  possdp = mapddp2(iidp,1)

  call fokunpck1(wrk(possfok),wrk(possdp),norb(symp))
  if (nvb(symp) > 0) call fokunpck2(wrk(possfok),wrk(possfaa),norb(symp),nvb(symp),nob(symp))
  if ((nob(symp)*nvb(symp)) > 0) call fokunpck3(wrk(possfok),wrk(possfai),norb(symp),nvb(symp),nob(symp))
  if (nob(symp) > 0) call fokunpck4(wrk(possfok),wrk(possfii),norb(symp),nob(symp))

end do

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(possfa0)
  call Unused_integer(possfb0)
  call Unused_integer(possfk10)
  call Unused_integer(possfk20)
  call Unused_integer(possfk30)
  call Unused_integer(possfk40)
  call Unused_integer(possfk50)
  call Unused_integer(possfk60)
  call Unused_integer(possdp10)
  call Unused_integer(possdp20)
end if

end subroutine divfok

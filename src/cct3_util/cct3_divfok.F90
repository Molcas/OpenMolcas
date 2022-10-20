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

subroutine cct3_divfok(wrk,wrksize,mapdfa,mapifa,mapdfb,mapifb,mapdfk1,mapifk1,mapdfk2,mapifk2,mapdfk3,mapifk3,mapdfk4,mapifk4, &
                       mapdfk5,mapifk5,mapdfk6,mapifk6,mapddp1,mapidp1,mapddp2,mapidp2,rc)
! this routine divides fok(p,q) -> fk(a,b) + fk(a,i) + f(i,j) + dp(p)
! to diagonal part and rest
!
! mapd and mapi for:
! fa,fb - fok(p,q)aa,bb
! fk1-6 - f(ab)aa,f(ab)bb,f(ai)aa,f(ai)bb,f(ij)aa,f(ij)bb
! dp1,2 - diagonal part dp(p)a,b
! rc    - return (error) code
!
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

use CCT3_global, only: noa, nob, norb, nsym, nva, nvb
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: wrksize, mapdfa(0:512,6), mapifa(8,8,8), mapdfb(0:512,6), mapifb(8,8,8), mapdfk1(0:512,6), mapifk1(8,8,8), &
                     mapdfk2(0:512,6), mapifk2(8,8,8), mapdfk3(0:512,6), mapifk3(8,8,8), mapdfk4(0:512,6), mapifk4(8,8,8), &
                     mapdfk5(0:512,6), mapifk5(8,8,8), mapdfk6(0:512,6), mapifk6(8,8,8), mapddp1(0:512,6), mapidp1(8,8,8), &
                     mapddp2(0:512,6), mapidp2(8,8,8), rc
real(kind=wp) :: wrk(wrksize)
integer(kind=iwp) :: iidp, iidpa, iidpb, iifaa, iifai, iifii, iifok, iifoka, iifokb, posdp, posdpa, posdpb, posfaa, posfai, &
                     posfii, posfok, posfoka, posfokb, rc1, symp

rc = 0

!1 define dp

do symp=1,nsym

  iidpa = mapidp1(symp,1,1)
  posdpa = mapddp1(iidpa,1)
  iidpb = mapidp2(symp,1,1)
  posdpb = mapddp2(iidpb,1)
  iifoka = mapifa(symp,1,1)
  posfoka = mapdfa(iifoka,1)
  iifokb = mapifb(symp,1,1)
  posfokb = mapdfb(iifokb,1)

  if (norb(symp) > 0) then
    call cct3_fokunpck5(symp,wrk(posfoka),wrk(posfokb),wrk(posdpa),wrk(posdpb),norb(symp),rc1)
  end if

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

  posfok = mapdfa(iifok,1)
  posfaa = mapdfk1(iifaa,1)
  posfai = mapdfk3(iifai,1)
  posfii = mapdfk5(iifii,1)
  posdp = mapddp1(iidp,1)

  call cct3_fokunpck1(wrk(posfok),wrk(posdp),norb(symp))
  if (nva(symp) > 0) then
    call cct3_fokunpck2(wrk(posfok),wrk(posfaa),norb(symp),nva(symp),noa(symp))
  end if
  if ((noa(symp)*nva(symp)) > 0) then
    call cct3_fokunpck3(wrk(posfok),wrk(posfai),norb(symp),nva(symp),noa(symp))
  end if
  if (noa(symp) > 0) then
    call cct3_fokunpck4(wrk(posfok),wrk(posfii),norb(symp),noa(symp))
  end if

  !2.2 alpha case

  iifok = mapifb(symp,1,1)
  iifaa = mapifk2(symp,1,1)
  iifai = mapifk4(symp,1,1)
  iifii = mapifk6(symp,1,1)
  iidp = mapidp2(symp,1,1)

  posfok = mapdfb(iifok,1)
  posfaa = mapdfk2(iifaa,1)
  posfai = mapdfk4(iifai,1)
  posfii = mapdfk6(iifii,1)
  posdp = mapddp2(iidp,1)

  call cct3_fokunpck1(wrk(posfok),wrk(posdp),norb(symp))
  if (nvb(symp) > 0) then
    call cct3_fokunpck2(wrk(posfok),wrk(posfaa),norb(symp),nvb(symp),nob(symp))
  end if
  if ((nob(symp)*nvb(symp)) > 0) then
    call cct3_fokunpck3(wrk(posfok),wrk(posfai),norb(symp),nvb(symp),nob(symp))
  end if
  if (nob(symp) > 0) then
    call cct3_fokunpck4(wrk(posfok),wrk(posfii),norb(symp),nob(symp))
  end if

end do

return

end subroutine cct3_divfok

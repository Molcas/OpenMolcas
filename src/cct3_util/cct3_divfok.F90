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

subroutine cct3_divfok(wrk,wrksize,mapdfa,mapifa,possfa0,mapdfb,mapifb,possfb0,mapdfk1,mapifk1,possfk10,mapdfk2,mapifk2,possfk20, &
                       mapdfk3,mapifk3,possfk30,mapdfk4,mapifk4,possfk40,mapdfk5,mapifk5,possfk50,mapdfk6,mapifk6,possfk60, &
                       mapddp1,mapidp1,possdp10,mapddp2,mapidp2,possdp20,rc)
! this routine divides fok(p,q) -> fk(a,b) + fk(a,i) + f(i,j) + dp(p)
! to diagonal part and rest
!
! mapd and mapi for:
! fa,fb - fok(p,q)aa,bb
! fk1-6 - f(ab)aa,f(ab)bb,f(ai)aa,f(ai)bb,f(ij)aa,f(ij)bb
! dp1,2 - diagonal part dp(p)a,b
! rc    - return (error) code

#include "t31.fh"
#include "wrk.fh"
integer rc
!1 maps for FOKA,FOKB
integer mapdfa(0:512,1:6)
integer mapifa(1:8,1:8,1:8)
integer possfa0
integer mapdfb(0:512,1:6)
integer mapifb(1:8,1:8,1:8)
integer possfb0
!2 maps for FK
!  FK1 - f(a,b)aa
!  FK2 - f(a,b)bb
!  FK3 - f(a,i)aa
!  FK4 - f(a,i)bb
!  FK5 - f(i,j)aa
!  FK6 - f(i,j)bb
integer mapdfk1(0:512,1:6)
integer mapifk1(1:8,1:8,1:8)
integer possfk10
integer mapdfk2(0:512,1:6)
integer mapifk2(1:8,1:8,1:8)
integer possfk20
integer mapdfk3(0:512,1:6)
integer mapifk3(1:8,1:8,1:8)
integer possfk30
integer mapdfk4(0:512,1:6)
integer mapifk4(1:8,1:8,1:8)
integer possfk40
integer mapdfk5(0:512,1:6)
integer mapifk5(1:8,1:8,1:8)
integer possfk50
integer mapdfk6(0:512,1:6)
integer mapifk6(1:8,1:8,1:8)
integer possfk60
!3 maps for DP - diagonal part
!  DP1 - dp(p)a
!  DP2 - dp(p)b
integer mapddp1(0:512,1:6)
integer mapidp1(1:8,1:8,1:8)
integer possdp10
integer mapddp2(0:512,1:6)
integer mapidp2(1:8,1:8,1:8)
integer possdp20
! help variables
integer symp, rc1
integer iifoka, iifokb, iifok, iifaa, iifai, iifii, iidpa, iidpb, iidp
integer possfoka, possfokb, possfok, possfaa, possfai, possfii
integer possdpa, possdpb, possdp

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

  if (norb(symp) > 0) then
    call cct3_fokunpck5(symp,wrk(possfoka),wrk(possfokb),wrk(possdpa),wrk(possdpb),norb(symp),rc1)
  end if

end do

!2 define faa,fai,fii

do symp=1,nsym
  if (norb(symp) == 0) goto 1000

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

  call cct3_fokunpck1(wrk(possfok),wrk(possdp),norb(symp))
  if (nva(symp) > 0) then
    call cct3_fokunpck2(wrk(possfok),wrk(possfaa),norb(symp),nva(symp),noa(symp))
  end if
  if ((noa(symp)*nva(symp)) > 0) then
    call cct3_fokunpck3(wrk(possfok),wrk(possfai),norb(symp),nva(symp),noa(symp))
  end if
  if (noa(symp) > 0) then
    call cct3_fokunpck4(wrk(possfok),wrk(possfii),norb(symp),noa(symp))
  end if

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

  call cct3_fokunpck1(wrk(possfok),wrk(possdp),norb(symp))
  if (nvb(symp) > 0) then
    call cct3_fokunpck2(wrk(possfok),wrk(possfaa),norb(symp),nvb(symp),nob(symp))
  end if
  if ((nob(symp)*nvb(symp)) > 0) then
    call cct3_fokunpck3(wrk(possfok),wrk(possfai),norb(symp),nvb(symp),nob(symp))
  end if
  if (nob(symp) > 0) then
    call cct3_fokunpck4(wrk(possfok),wrk(possfii),norb(symp),nob(symp))
  end if

1000 continue
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

end subroutine cct3_divfok

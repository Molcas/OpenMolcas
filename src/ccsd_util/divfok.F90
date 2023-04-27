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

subroutine divfok(wrk,wrksize,fa,fb,fk1,fk2,fk3,fk4,fk5,fk6,dp1,dp2,rc)
! this routine divides fok(p,q) -> fk(a,b) + fk(a,i) + f(i,j) + dp(p)
! to diagonal part and rest
!
! map type for:
! fa,fb - fok(p,q)aa,bb
! fk1-6 - f(ab)aa,f(ab)bb,f(ai)aa,f(ai)bb,f(ij)aa,f(ij)bb
! dp1,2 - diagonal part dp(p)a,b
! rc    - return (error) code

use ccsd_global, only: Map_Type, noa, nob, norb, nsym, nva, nvb
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
integer(kind=iwp), intent(in) :: wrksize
real(kind=wp), intent(inout) :: wrk(wrksize)
type(Map_Type), intent(in) :: fa, fb, fk1, fk2, fk3, fk4, fk5, fk6, dp1, dp2
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp) :: iidp, iidpa, iidpb, iifaa, iifai, iifii, iifok, iifoka, iifokb, posdp, posdpa, posdpb, posfaa, posfai, &
                     posfii, posfok, posfoka, posfokb, rc1, symp

rc = 0

!1 define dp

do symp=1,nsym

  iidpa = dp1%i(symp,1,1)
  posdpa = dp1%d(iidpa,1)
  iidpb = dp2%i(symp,1,1)
  posdpb = dp2%d(iidpb,1)
  iifoka = fa%i(symp,1,1)
  posfoka = fa%d(iifoka,1)
  iifokb = fb%i(symp,1,1)
  posfokb = fb%d(iifokb,1)

  if (norb(symp) > 0) call fokunpck5(symp,wrk(posfoka),wrk(posfokb),wrk(posdpa),wrk(posdpb),norb(symp),rc1)

end do

!2 define faa,fai,fii

do symp=1,nsym
  if (norb(symp) == 0) cycle

  !2.1 alpha case

  iifok = fa%i(symp,1,1)
  iifaa = fk1%i(symp,1,1)
  iifai = fk3%i(symp,1,1)
  iifii = fk5%i(symp,1,1)
  iidp = dp1%i(symp,1,1)

  posfok = fa%d(iifok,1)
  posfaa = fk1%d(iifaa,1)
  posfai = fk3%d(iifai,1)
  posfii = fk5%d(iifii,1)
  posdp = dp1%d(iidp,1)

  call fokunpck1(wrk(posfok),wrk(posdp),norb(symp))
  if (nva(symp) > 0) call fokunpck2(wrk(posfok),wrk(posfaa),norb(symp),nva(symp),noa(symp))
  if ((noa(symp)*nva(symp)) > 0) call fokunpck3(wrk(posfok),wrk(posfai),norb(symp),nva(symp),noa(symp))
  if (noa(symp) > 0) call fokunpck4(wrk(posfok),wrk(posfii),norb(symp),noa(symp))

  !2.2 alpha case

  iifok = fb%i(symp,1,1)
  iifaa = fk2%i(symp,1,1)
  iifai = fk4%i(symp,1,1)
  iifii = fk6%i(symp,1,1)
  iidp = dp2%i(symp,1,1)

  posfok = fb%d(iifok,1)
  posfaa = fk2%d(iifaa,1)
  posfai = fk4%d(iifai,1)
  posfii = fk6%d(iifii,1)
  posdp = dp2%d(iidp,1)

  call fokunpck1(wrk(posfok),wrk(posdp),norb(symp))
  if (nvb(symp) > 0) call fokunpck2(wrk(posfok),wrk(posfaa),norb(symp),nvb(symp),nob(symp))
  if ((nob(symp)*nvb(symp)) > 0) call fokunpck3(wrk(posfok),wrk(posfai),norb(symp),nvb(symp),nob(symp))
  if (nob(symp) > 0) call fokunpck4(wrk(posfok),wrk(posfii),norb(symp),nob(symp))

end do

return

end subroutine divfok

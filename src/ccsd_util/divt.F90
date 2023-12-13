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

subroutine divt(wrk,wrksize,nind,t,dp1,dp2,rc)
! this routine divides T amplitudes by denominators, i.e. differences in dp
!
! nint - number of indices in T (2 or 4) (I)
! t    - map type of T (I)
! dp1  - map type of dpa (I)
! dp2  - map type of dpb (I)
! rc   - return (error) code

use ccsd_global, only: Map_Type, noa, nob, nva, nvb
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize, nind
real(kind=wp), intent(inout) :: wrk(wrksize)
type(Map_Type), intent(in) :: t, dp1, dp2
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp) :: dima, dimab, dimb, dimi, dimij, dimj, iidp, iidpa, iidpb, iidpi, iidpj, iit, posdp, posdpa, posdpb, posdpi, &
                     posdpj, post, syma, symb, symi, symj

rc = 0

if (nind == 2) then
  !I T1 amplitudes - 2 indices

  if (t%d(0,1) == 3) then
    !I.1 T1aa case

    do iit=1,t%d(0,5)

      post = t%d(iit,1)
      syma = t%d(iit,3)
      dima = nva(syma)
      dimi = noa(syma)
      iidp = dp1%i(syma,1,1)
      posdp = dp1%d(iidp,1)

      if ((dima*dimi) > 0) call divthelp1(wrk(post),dima,dimi,wrk(posdp))

    end do

  else if (t%d(0,1) == 4) then
    !I.2 T1bb case

    do iit=1,t%d(0,5)

      post = t%d(iit,1)
      syma = t%d(iit,3)
      dima = nvb(syma)
      dimi = nob(syma)
      iidp = dp2%i(syma,1,1)
      posdp = dp2%d(iidp,1)

      if ((dima*dimi) > 0) call divthelp1(wrk(post),dima,dimi,wrk(posdp))

    end do

  else
    !I.3 invalid T1 type
    ! RC=1 : incorrect map for T1
    rc = 1
    return
  end if

else if (nind == 4) then
  !II T2 amplitudes - 4 indices

  if (t%d(0,6) == 0) then
    !II.1 T2abab case

    do iit=1,t%d(0,5)

      post = t%d(iit,1)
      syma = t%d(iit,3)
      symb = t%d(iit,4)
      symi = t%d(iit,5)
      symj = t%d(iit,6)
      dima = nva(syma)
      dimb = nvb(symb)
      dimi = noa(symi)
      dimj = nob(symj)
      iidpa = dp1%i(syma,1,1)
      iidpb = dp2%i(symb,1,1)
      iidpi = dp1%i(symi,1,1)
      iidpj = dp2%i(symj,1,1)
      posdpa = dp1%d(iidpa,1)
      posdpb = dp2%d(iidpb,1)
      posdpi = dp1%d(iidpi,1)
      posdpj = dp2%d(iidpj,1)

      if (t%d(iit,2) > 0) &
        call divthelp2(wrk(post),dima,dimb,dimi,dimj,wrk(posdpa),wrk(posdpb),wrk(posdpi),wrk(posdpj),noa(syma),nob(symb))

    end do

  else if ((t%d(0,6) == 4) .and. (t%d(0,1) == 3)) then
    !II.2 T2aaaa case

    do iit=1,t%d(0,5)

      post = t%d(iit,1)
      syma = t%d(iit,3)
      symb = t%d(iit,4)
      symi = t%d(iit,5)
      symj = t%d(iit,6)
      dima = nva(syma)
      dimb = nva(symb)
      dimi = noa(symi)
      dimj = noa(symj)
      iidpa = dp1%i(syma,1,1)
      iidpb = dp1%i(symb,1,1)
      iidpi = dp1%i(symi,1,1)
      iidpj = dp1%i(symj,1,1)
      posdpa = dp1%d(iidpa,1)
      posdpb = dp1%d(iidpb,1)
      posdpi = dp1%d(iidpi,1)
      posdpj = dp1%d(iidpj,1)

      if (t%d(iit,2) == 0) cycle

      if (syma /= symb) then
        ! different symmetries a,b; i,j
        call divthelp2(wrk(post),dima,dimb,dimi,dimj,wrk(posdpa),wrk(posdpb),wrk(posdpi),wrk(posdpj),noa(syma),noa(symb))

      else
        ! same symmetries a,b; i,j
        dimab = (dima*(dima-1))/2
        dimij = (dimi*(dimi-1))/2
        call divthelp3(wrk(post),dimab,dimij,wrk(posdpa),wrk(posdpi),dima,dimi,noa(syma))
      end if

    end do

  else if ((t%d(0,6) == 4) .and. (t%d(0,1) == 4)) then
    !II.3 T2bbbb case

    do iit=1,t%d(0,5)

      post = t%d(iit,1)
      syma = t%d(iit,3)
      symb = t%d(iit,4)
      symi = t%d(iit,5)
      symj = t%d(iit,6)
      dima = nvb(syma)
      dimb = nvb(symb)
      dimi = nob(symi)
      dimj = nob(symj)
      iidpa = dp2%i(syma,1,1)
      iidpb = dp2%i(symb,1,1)
      iidpi = dp2%i(symi,1,1)
      iidpj = dp2%i(symj,1,1)
      posdpa = dp2%d(iidpa,1)
      posdpb = dp2%d(iidpb,1)
      posdpi = dp2%d(iidpi,1)
      posdpj = dp2%d(iidpj,1)

      if (t%d(iit,2) == 0) cycle

      if (syma /= symb) then
        ! different symmetries a,b; i,j
        call divthelp2(wrk(post),dima,dimb,dimi,dimj,wrk(posdpa),wrk(posdpb),wrk(posdpi),wrk(posdpj),nob(syma),nob(symb))

      else
        ! same symmetries a,b; i,j
        dimab = (dima*(dima-1))/2
        dimij = (dimi*(dimi-1))/2
        call divthelp3(wrk(post),dimab,dimij,wrk(posdpa),wrk(posdpi),dima,dimi,nob(syma))
      end if

    end do

  else
    !II.4 RC=2 : incorrect t%d for T2
    rc = 2
    return
  end if

else
  !III invalid nind
  ! RC=3 : nind is not 2 or 4 (Stup)
  rc = 3
end if

return

end subroutine divt

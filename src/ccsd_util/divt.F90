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

subroutine divt(wrk,wrksize,nind,mapdt,mapit,mapddp1,mapidp1,mapddp2,mapidp2,rc)
! this routine divides T amplitudes by denominators, i.e. differences in dp
!
! nint    - number of indices in T (2 or 4) (I)
! mapdt   - direct map of T (I)
! mapit   - inverse map of T (I)
! mapddp1 - direct map of dpa (I)
! mapidp1 - inverse map of dpa (I)
! mapddp2 - direct map of dpb (I)
! mapidp2 - inverse map of dpb (I)
! rc      - return (error) code

#include "ccsd1.fh"
#include "wrk.fh"
integer nind, rc
integer mapdt(0:512,1:6)
integer mapit(1:8,1:8,1:8)
integer mapddp1(0:512,1:6)
integer mapidp1(1:8,1:8,1:8)
integer mapddp2(0:512,1:6)
integer mapidp2(1:8,1:8,1:8)
! help variables
integer posst, possdp, possdpa, possdpb, possdpi, possdpj
integer dimi, dimj, dima, dimb, dimab, dimij, syma, symb, symi, symj
integer iit, iidp, iidpa, iidpb, iidpi, iidpj

rc = 0

if (nind == 2) then
  !I T1 amplitudes - 2 indices

  if (mapdt(0,1) == 3) then
    !I.1 T1aa case

    do iit=1,mapdt(0,5)

      posst = mapdt(iit,1)
      syma = mapdt(iit,3)
      dima = nva(syma)
      dimi = noa(syma)
      iidp = mapidp1(syma,1,1)
      possdp = mapddp1(iidp,1)

      if ((dima*dimi) > 0) then
        call divthelp1(wrk(posst),dima,dimi,wrk(possdp))
      end if

    end do

  else if (mapdt(0,1) == 4) then
    !I.2 T1bb case

    do iit=1,mapdt(0,5)

      posst = mapdt(iit,1)
      syma = mapdt(iit,3)
      dima = nvb(syma)
      dimi = nob(syma)
      iidp = mapidp2(syma,1,1)
      possdp = mapddp2(iidp,1)

      if ((dima*dimi) > 0) then
        call divthelp1(wrk(posst),dima,dimi,wrk(possdp))
      end if

    end do

  else
    !I.3 invalid T1 type
    ! RC=1 : incorrect map for T1
    rc = 1
    return
  end if

else if (nind == 4) then
  !II T2 amplitudes - 4 indices

  if (mapdt(0,6) == 0) then
    !II.1 T2abab case

    do iit=1,mapdt(0,5)

      posst = mapdt(iit,1)
      syma = mapdt(iit,3)
      symb = mapdt(iit,4)
      symi = mapdt(iit,5)
      symj = mapdt(iit,6)
      dima = nva(syma)
      dimb = nvb(symb)
      dimi = noa(symi)
      dimj = nob(symj)
      iidpa = mapidp1(syma,1,1)
      iidpb = mapidp2(symb,1,1)
      iidpi = mapidp1(symi,1,1)
      iidpj = mapidp2(symj,1,1)
      possdpa = mapddp1(iidpa,1)
      possdpb = mapddp2(iidpb,1)
      possdpi = mapddp1(iidpi,1)
      possdpj = mapddp2(iidpj,1)

      if (mapdt(iit,2) > 0) then
        call divthelp2(wrk(posst),dima,dimb,dimi,dimj,wrk(possdpa),wrk(possdpb),wrk(possdpi),wrk(possdpj),noa(syma),nob(symb))
      end if

    end do

  else if ((mapdt(0,6) == 4) .and. (mapdt(0,1) == 3)) then
    !II.2 T2aaaa case

    do iit=1,mapdt(0,5)

      posst = mapdt(iit,1)
      syma = mapdt(iit,3)
      symb = mapdt(iit,4)
      symi = mapdt(iit,5)
      symj = mapdt(iit,6)
      dima = nva(syma)
      dimb = nva(symb)
      dimi = noa(symi)
      dimj = noa(symj)
      iidpa = mapidp1(syma,1,1)
      iidpb = mapidp1(symb,1,1)
      iidpi = mapidp1(symi,1,1)
      iidpj = mapidp1(symj,1,1)
      possdpa = mapddp1(iidpa,1)
      possdpb = mapddp1(iidpb,1)
      possdpi = mapddp1(iidpi,1)
      possdpj = mapddp1(iidpj,1)

      if (mapdt(iit,2) == 0) goto 400

      if (syma /= symb) then
        ! different symmetries a,b; i,j
        call divthelp2(wrk(posst),dima,dimb,dimi,dimj,wrk(possdpa),wrk(possdpb),wrk(possdpi),wrk(possdpj),noa(syma),noa(symb))

      else
        ! same symmetries a,b; i,j
        dimab = (dima*(dima-1))/2
        dimij = (dimi*(dimi-1))/2
        call divthelp3(wrk(posst),dimab,dimij,wrk(possdpa),wrk(possdpi),dima,dimi,noa(syma))
      end if

      400 continue
    end do

  else if ((mapdt(0,6) == 4) .and. (mapdt(0,1) == 4)) then
    !II.3 T2bbbb case

    do iit=1,mapdt(0,5)

      posst = mapdt(iit,1)
      syma = mapdt(iit,3)
      symb = mapdt(iit,4)
      symi = mapdt(iit,5)
      symj = mapdt(iit,6)
      dima = nvb(syma)
      dimb = nvb(symb)
      dimi = nob(symi)
      dimj = nob(symj)
      iidpa = mapidp2(syma,1,1)
      iidpb = mapidp2(symb,1,1)
      iidpi = mapidp2(symi,1,1)
      iidpj = mapidp2(symj,1,1)
      possdpa = mapddp2(iidpa,1)
      possdpb = mapddp2(iidpb,1)
      possdpi = mapddp2(iidpi,1)
      possdpj = mapddp2(iidpj,1)

      if (mapdt(iit,2) == 0) goto 500

      if (syma /= symb) then
        ! different symmetries a,b; i,j
        call divthelp2(wrk(posst),dima,dimb,dimi,dimj,wrk(possdpa),wrk(possdpb),wrk(possdpi),wrk(possdpj),nob(syma),nob(symb))

      else
        ! same symmetries a,b; i,j
        dimab = (dima*(dima-1))/2
        dimij = (dimi*(dimi-1))/2
        call divthelp3(wrk(posst),dimab,dimij,wrk(possdpa),wrk(possdpi),dima,dimi,nob(syma))
      end if

      500 continue
    end do

  else
    !II.4 RC=2 : incorrect mapdt for T2
    rc = 2
    return
  end if

else
  !III invalid nind
  ! RC=3 : nind is not 2 or 4 (Stup)
  rc = 3
end if

return
! Avoid unused argument warnings
if (.false.) call Unused_integer_array(mapit)

end subroutine divt

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

subroutine mktau(wrk,wrksize,mapdt2,mapit2,mapdt1a,mapit1a,mapdt1b,mapit1b,fact,rc)
! this routine does:
! t2(abij) = t2(abij) + fact. (t1(ai).t1(bj)-t1(bi).t1(aj))
! N.B. T24a,4b must be of type 4, T2abab of type 0
!
! mapdt2  - direct map of T2 (I)
! mapit2  - inverse map of T2 (I)
! mapdt1a - direct map of T1aa (I)
! mapit1a - inverse map of T1aa (I)
! mapdt1b - direct map of T1bb (I)
! mapit1b - inverse map of T1bb (I)
! fact    - numerical factor (I)
! rc      - return (error) code

use ccsd_global, only: noa, nob, nva, nvb
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: wrksize, mapdt2(0:512,6), mapit2(8,8,8), mapdt1a(0:512,6), mapit1a(8,8,8), mapdt1b(0:512,6), mapit1b(8,8,8), rc
real(kind=wp) :: wrk(wrksize), fact
integer(kind=iwp) :: dima, dimab, dimb, dimi, dimij, dimj, iit11, iit12, iit1a, iit1b, iit2, posst11, posst12, posst1a, posst1b, &
                     posst2, syma, symb, symi, symj

rc = 0

if (mapdt2(0,6) == 0) then
  !I.1 T2abab case

  do iit2=1,mapdt2(0,5)

    posst2 = mapdt2(iit2,1)
    syma = mapdt2(iit2,3)
    symb = mapdt2(iit2,4)
    symi = mapdt2(iit2,5)
    symj = mapdt2(iit2,6)
    dima = nva(syma)
    dimb = nvb(symb)
    dimi = noa(symi)
    dimj = nob(symj)
    iit1a = mapit1a(syma,1,1)
    iit1b = mapit1b(symb,1,1)
    posst1a = mapdt1a(iit1a,1)
    posst1b = mapdt1b(iit1b,1)

    if ((syma == symi) .and. (symb == symj) .and. (mapdt2(iit2,2) > 0)) &
      call mktauhelp1(wrk(posst2),wrk(posst1a),wrk(posst1b),dima,dimb,dimi,dimj,fact)

  end do

else if ((mapdt2(0,6) == 4) .and. (mapdt2(0,1) == 3)) then
  !I.2 T2aaaa case

  do iit2=1,mapdt2(0,5)

    posst2 = mapdt2(iit2,1)
    syma = mapdt2(iit2,3)
    symb = mapdt2(iit2,4)
    symi = mapdt2(iit2,5)
    symj = mapdt2(iit2,6)
    dima = nva(syma)
    dimb = nva(symb)
    dimi = noa(symi)
    dimj = noa(symj)
    iit11 = mapit1a(syma,1,1)
    iit12 = mapit1a(symb,1,1)
    posst11 = mapdt1a(iit11,1)
    posst12 = mapdt1a(iit12,1)

    if ((syma == symi) .and. (symb == symj) .and. (syma /= symj) .and. (mapdt2(iit2,2) > 0)) then
      !I.2.* case T2(sym1,sym2,sym1,sym2)

      call mktauhelp1(wrk(posst2),wrk(posst11),wrk(posst12),dima,dimb,dimi,dimj,fact)

    else if ((syma == symi) .and. (symb == symj) .and. (syma == symj) .and. (mapdt2(iit2,2) > 0)) then
      !I.2.* case T2(sym1,sym1,sym1,sym1)

      dimab = (dima*(dima-1))/2
      dimij = (dimi*(dimi-1))/2
      call mktauhelp2(wrk(posst2),wrk(posst11),dimab,dimij,dima,dimi,fact)

    end if

  end do

else if ((mapdt2(0,6) == 4) .and. (mapdt2(0,1) == 4)) then
  !I.3 T2bbbb case

  do iit2=1,mapdt2(0,5)

    posst2 = mapdt2(iit2,1)
    syma = mapdt2(iit2,3)
    symb = mapdt2(iit2,4)
    symi = mapdt2(iit2,5)
    symj = mapdt2(iit2,6)
    dima = nvb(syma)
    dimb = nvb(symb)
    dimi = nob(symi)
    dimj = nob(symj)
    iit11 = mapit1b(syma,1,1)
    iit12 = mapit1b(symb,1,1)
    posst11 = mapdt1b(iit11,1)
    posst12 = mapdt1b(iit12,1)

    if ((syma == symi) .and. (symb == symj) .and. (syma /= symj) .and. (mapdt2(iit2,2) > 0)) then
      !I.3.* case T2(sym1,sym2,sym1,sym2)

      call mktauhelp1(wrk(posst2),wrk(posst11),wrk(posst12),dima,dimb,dimi,dimj,fact)

    else if ((syma == symi) .and. (symb == symj) .and. (syma == symj) .and. (mapdt2(iit2,2) > 0)) then
      !I.3.* case T2(sym1,sym1,sym1,sym1)

      dimab = (dima*(dima-1))/2
      dimij = (dimi*(dimi-1))/2
      call mktauhelp2(wrk(posst2),wrk(posst11),dimab,dimij,dima,dimi,fact)

    end if

  end do

else
  !I.4 RC=1 : incorrect mapdt for T2
  rc = 1
  return
end if

return
! Avoid unused argument warnings
if (.false.) call Unused_integer_array(mapit2)

end subroutine mktau

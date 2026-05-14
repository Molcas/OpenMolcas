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

subroutine mktau(wrk,wrksize,t2,t1a,t1b,fact,rc)
! this routine does:
! t2(abij) = t2(abij) + fact. (t1(ai).t1(bj)-t1(bi).t1(aj))
! N.B. T24a,4b must be of type 4, T2abab of type 0
!
! t2   - map type of T2 (I)
! t1a  - map type of T1aa (I)
! t1b  - map type of T1bb (I)
! fact - numerical factor (I)
! rc   - return (error) code

use ccsd_global, only: Map_Type, noa, nob, nva, nvb
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize
real(kind=wp), intent(inout) :: wrk(wrksize)
type(Map_Type), intent(in) :: t2, t1a, t1b
real(kind=wp), intent(in) :: fact
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp) :: dima, dimab, dimb, dimi, dimij, dimj, iit11, iit12, iit1a, iit1b, iit2, post11, post12, post1a, post1b, &
                     post2, syma, symb, symi, symj

rc = 0

if (t2%d(0,6) == 0) then
  !I.1 T2abab case

  do iit2=1,t2%d(0,5)

    post2 = t2%d(iit2,1)
    syma = t2%d(iit2,3)
    symb = t2%d(iit2,4)
    symi = t2%d(iit2,5)
    symj = t2%d(iit2,6)
    dima = nva(syma)
    dimb = nvb(symb)
    dimi = noa(symi)
    dimj = nob(symj)
    iit1a = t1a%i(syma,1,1)
    iit1b = t1b%i(symb,1,1)
    post1a = t1a%d(iit1a,1)
    post1b = t1b%d(iit1b,1)

    if ((syma == symi) .and. (symb == symj) .and. (t2%d(iit2,2) > 0)) &
      call mktauhelp1(wrk(post2),wrk(post1a),wrk(post1b),dima,dimb,dimi,dimj,fact)

  end do

else if ((t2%d(0,6) == 4) .and. (t2%d(0,1) == 3)) then
  !I.2 T2aaaa case

  do iit2=1,t2%d(0,5)

    post2 = t2%d(iit2,1)
    syma = t2%d(iit2,3)
    symb = t2%d(iit2,4)
    symi = t2%d(iit2,5)
    symj = t2%d(iit2,6)
    dima = nva(syma)
    dimb = nva(symb)
    dimi = noa(symi)
    dimj = noa(symj)
    iit11 = t1a%i(syma,1,1)
    iit12 = t1a%i(symb,1,1)
    post11 = t1a%d(iit11,1)
    post12 = t1a%d(iit12,1)

    if ((syma == symi) .and. (symb == symj) .and. (syma /= symj) .and. (t2%d(iit2,2) > 0)) then
      !I.2.* case T2(sym1,sym2,sym1,sym2)

      call mktauhelp1(wrk(post2),wrk(post11),wrk(post12),dima,dimb,dimi,dimj,fact)

    else if ((syma == symi) .and. (symb == symj) .and. (syma == symj) .and. (t2%d(iit2,2) > 0)) then
      !I.2.* case T2(sym1,sym1,sym1,sym1)

      dimab = (dima*(dima-1))/2
      dimij = (dimi*(dimi-1))/2
      call mktauhelp2(wrk(post2),wrk(post11),dimab,dimij,dima,dimi,fact)

    end if

  end do

else if ((t2%d(0,6) == 4) .and. (t2%d(0,1) == 4)) then
  !I.3 T2bbbb case

  do iit2=1,t2%d(0,5)

    post2 = t2%d(iit2,1)
    syma = t2%d(iit2,3)
    symb = t2%d(iit2,4)
    symi = t2%d(iit2,5)
    symj = t2%d(iit2,6)
    dima = nvb(syma)
    dimb = nvb(symb)
    dimi = nob(symi)
    dimj = nob(symj)
    iit11 = t1b%i(syma,1,1)
    iit12 = t1b%i(symb,1,1)
    post11 = t1b%d(iit11,1)
    post12 = t1b%d(iit12,1)

    if ((syma == symi) .and. (symb == symj) .and. (syma /= symj) .and. (t2%d(iit2,2) > 0)) then
      !I.3.* case T2(sym1,sym2,sym1,sym2)

      call mktauhelp1(wrk(post2),wrk(post11),wrk(post12),dima,dimb,dimi,dimj,fact)

    else if ((syma == symi) .and. (symb == symj) .and. (syma == symj) .and. (t2%d(iit2,2) > 0)) then
      !I.3.* case T2(sym1,sym1,sym1,sym1)

      dimab = (dima*(dima-1))/2
      dimij = (dimi*(dimi-1))/2
      call mktauhelp2(wrk(post2),wrk(post11),dimab,dimij,dima,dimi,fact)

    end if

  end do

else
  !I.4 RC=1 : incorrect %d for T2
  rc = 1
  return
end if

return

end subroutine mktau

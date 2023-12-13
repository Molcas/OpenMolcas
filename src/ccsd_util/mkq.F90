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

subroutine mkq(wrk,wrksize,t2,t11,t12,fact,rc)
! this routine does:
! t2(a,b,i,j) = fact . t2(a,b,i,j) +t11(ai).t12(bj)
! for T2aaaa, T2bbbb and T2abab but they must be expanded (typ=0)
!
! t2   - map type of T2 (I)
! t11  - map type of T11 (I)
! t12  - map type of T12 (I)
! fact - numerical factor (I)
! rc   - return (error) code

use ccsd_global, only: dimm, Map_Type
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize
real(kind=wp), intent(inout) :: wrk(wrksize)
type(Map_Type), intent(in) :: t2, t11, t12
real(kind=wp), intent(in) :: fact
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp) :: dima, dimb, dimi, dimj, iit11, iit12, iit2, post11, post12, post2, syma, symb, symi, symj

rc = 0

if (t2%d(0,6) == 0) then
  !I.1 typ of t2 is 0 (T2 is expanded)

  do iit2=1,t2%d(0,5)

    post2 = t2%d(iit2,1)
    syma = t2%d(iit2,3)
    symb = t2%d(iit2,4)
    symi = t2%d(iit2,5)
    symj = t2%d(iit2,6)
    dima = dimm(t2%d(0,1),syma)
    dimb = dimm(t2%d(0,2),symb)
    dimi = dimm(t2%d(0,3),symi)
    dimj = dimm(t2%d(0,4),symj)
    iit11 = t11%i(syma,1,1)
    iit12 = t12%i(symb,1,1)
    post11 = t11%d(iit11,1)
    post12 = t12%d(iit12,1)

    if ((syma == symi) .and. (symb == symj) .and. (t2%d(iit2,2) > 0)) then
      call mkqhelp1(wrk(post2),wrk(post11),wrk(post12),dima,dimb,dimi,dimj,fact)
    else if (t2%d(iit2,2) > 0) then
      call mkqhelp2(wrk(post2),t2%d(iit2,2),t2%d(iit2,2),fact)
    end if

  end do

else
  !I.4 RC=1 : typ of T2 is not 0
  rc = 1
  return
end if

return

end subroutine mkq

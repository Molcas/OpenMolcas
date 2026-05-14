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

subroutine saamp(wrk,wrksize,key)
! this routine rearranges amplitudes to be spin adapted
! key - 0 - no adaptation
!       1 - T2 DDVV adaptation
!       2 - T2 DDVV + T1 DV adaptation
!       3 - full T1 and T2 adaptation (only for doublets)
!       4 - full T2 without SDVS (only for doublets)
!
! amplitudes T1 are in t13 - aa, t14 - bb
! T2 are in t21 - aaaa, t22 - bbbb, t23 - abab

use ccsd_global, only: dimm, mmul, nsym, t13, t14, t21, t22, t23
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: wrksize, key
real(kind=wp), intent(inout) :: wrk(wrksize)
integer(kind=iwp) :: ii, pos1, pos2, pos3, pos4, pos5, pos6, syma, symb, symi, symij, symj, syms

!0 skip this routine if SA in not turn on
if (key == 0) return

!I T1 adaptation
if ((key == 2) .or. (key == 3)) then

  !I.1 def symmetry, where S orbital is situated (only for doublet states)
  syms = 0
  do symi=1,nsym
    if (dimm(1,symi) /= dimm(2,symi)) syms = symi
  end do
  if ((key == 2) .and. (syms == 0)) syms = 1
  if (syms == 0) then
    write(u6,*) ' Full SA is turn on and there is no S orbitals'
    call Abend()
  end if

  !I.2 loop over symi
  do symi=1,nsym
    syma = symi

    ii = t13%i(syma,1,1)
    pos1 = t13%d(ii,1)
    ii = t14%i(syma,1,1)
    pos2 = t14%d(ii,1)
    ii = t23%i(syma,syms,syms)
    pos3 = t23%d(ii,1)
    call saamphlp3(wrk(pos1),wrk(pos2),wrk(pos3),dimm(1,symi),dimm(2,symi),dimm(3,symi),dimm(4,symi),dimm(1,syms),dimm(4,syms),key)

  end do
end if

!II T2 adaptation

do symi=1,nsym
  do symj=1,symi
    symij = mmul(symi,symj)

    do syma=1,nsym
      symb = mmul(symij,syma)

      ! Meggie out
      if (symb > syma) cycle

      if (symi == symj) then
        ! case si=sj, sa=sb

        ii = t21%i(syma,symb,symi)
        pos1 = t21%d(ii,1)
        ii = t22%i(syma,symb,symi)
        pos2 = t22%d(ii,1)
        ii = t23%i(syma,symb,symi)
        pos3 = t23%d(ii,1)
        call saamphlp1(wrk(pos1),wrk(pos2),wrk(pos3),dimm(1,symi),dimm(2,symi),dimm(3,syma),dimm(4,syma),key)

      else
        ! case si>sj, sa>sb

        ii = t21%i(syma,symb,symi)
        pos1 = t21%d(ii,1)
        ii = t22%i(syma,symb,symi)
        pos2 = t22%d(ii,1)
        ii = t23%i(syma,symb,symi)
        pos3 = t23%d(ii,1)
        ii = t23%i(symb,syma,symj)
        pos4 = t23%d(ii,1)

        ii = t23%i(symb,syma,symi)
        pos5 = t23%d(ii,1)
        ii = t23%i(syma,symb,symj)
        pos6 = t23%d(ii,1)
        call saamphlp2(wrk(pos1),wrk(pos2),wrk(pos3),wrk(pos4),wrk(pos5),wrk(pos6),dimm(1,symi),dimm(1,symj),dimm(2,symi), &
                       dimm(2,symj),dimm(3,syma),dimm(3,symb),dimm(4,syma),dimm(4,symb),key)

      end if

    end do
  end do
end do

return

end subroutine saamp

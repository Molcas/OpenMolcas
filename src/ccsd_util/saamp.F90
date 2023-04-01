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
! 1 - T2 DDVV adaptation
! 2 - T2 DDVV + T1 DV adaptation
! 3 - full T1 and T2 adaptation (only for doublets)
! 4 - full T2 without SDVS (only for doublets)
!
! amplitudes T1 are in t13 - aa, t14 - bb
! T2 are in t21 - aaaa, t22 - bbbb, t23 - abab

#include "ccsd1.fh"
#include "ccsd2.fh"
#include "wrk.fh"
integer key
! help variables
integer symi, symj, syma, symb, syms, symij
integer poss1, poss2, poss3, poss4, poss5, poss6, ii

!0 skip this routine if SA in not turn on
if (key == 0) then
  return
end if

!I T1 adaptation
if ((key == 2) .or. (key == 3)) then

  !I.1 def symmetry, where S orbital is situated (only for doublet states)
  syms = 0
  do symi=1,nsym
    if (dimm(1,symi) /= dimm(2,symi)) then
      syms = symi
    end if
  end do
  if ((key == 2) .and. (syms == 0)) then
    syms = 1
  end if
  if (syms == 0) then
    write(6,*) ' Full SA is turn on and there is no S orbitals'
    call Abend()
  end if

  !I.2 loop over symi
  do symi=1,nsym
    syma = symi

    ii = mapit13(syma,1,1)
    poss1 = mapdt13(ii,1)
    ii = mapit14(syma,1,1)
    poss2 = mapdt14(ii,1)
    ii = mapit23(syma,syms,syms)
    poss3 = mapdt23(ii,1)
    call saamphlp3(wrk(poss1),wrk(poss2),wrk(poss3),dimm(1,symi),dimm(2,symi),dimm(3,symi),dimm(4,symi),dimm(1,syms),dimm(4,syms), &
                   key)

  end do
end if

!II T2 adaptation

do symi=1,nsym
  do symj=1,symi
    symij = mmul(symi,symj)

    do syma=1,nsym
      symb = mmul(symij,syma)

      if (symb > syma) then
        ! Meggie out
        goto 202
      end if

      if (symi == symj) then
        ! case si=sj, sa=sb

        ii = mapit21(syma,symb,symi)
        poss1 = mapdt21(ii,1)
        ii = mapit22(syma,symb,symi)
        poss2 = mapdt22(ii,1)
        ii = mapit23(syma,symb,symi)
        poss3 = mapdt23(ii,1)
        call saamphlp1(wrk(poss1),wrk(poss2),wrk(poss3),dimm(1,symi),dimm(2,symi),dimm(3,syma),dimm(4,syma),key)

      else
        ! case si>sj, sa>sb

        ii = mapit21(syma,symb,symi)
        poss1 = mapdt21(ii,1)
        ii = mapit22(syma,symb,symi)
        poss2 = mapdt22(ii,1)
        ii = mapit23(syma,symb,symi)
        poss3 = mapdt23(ii,1)
        ii = mapit23(symb,syma,symj)
        poss4 = mapdt23(ii,1)

        ii = mapit23(symb,syma,symi)
        poss5 = mapdt23(ii,1)
        ii = mapit23(syma,symb,symj)
        poss6 = mapdt23(ii,1)
        call saamphlp2(wrk(poss1),wrk(poss2),wrk(poss3),wrk(poss4),wrk(poss5),wrk(poss6),dimm(1,symi),dimm(1,symj),dimm(2,symi), &
                       dimm(2,symj),dimm(3,syma),dimm(3,symb),dimm(4,syma),dimm(4,symb),key)

      end if

      202 continue
    end do
  end do
end do

return

end subroutine saamp

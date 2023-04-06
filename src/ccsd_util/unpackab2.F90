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

subroutine unpackab2(wrk,wrksize,n,r1,r2,r3,r4,r5,r6,ssn,key,aeqb)
! n    - map type corresponding to N (Input)
! r1   - map type corresponding to R1 (Input)
! r2   - map type corresponding to R2 (Input)
! r3   - map type corresponding to R3 (Input)
! r4   - map type corresponding to R4 (Input)
! r5   - map type corresponding to R4 (Input)
! r6   - map type corresponding to R4 (Input)
! ssn  - overall symmetry state of N (Input)
! key  - define state of _a,_b  key   a, b (Input)
!        1    S  S
!        2    V  S
!        3    S  V
!        4    V  V
! aqeb - aqeb=0 if a /= b; =1 if a=b (Input)
!
! this routine unpacks: N _a_b (p,q) = <ab|pq>             key   aeqb
! N _a_b (p,q) a>=b-aa  -> R1 _a,_b(j,e)aaaa = <ab||je>     4     0
!                       -> R2 _a,_b(j,e)bbbb = <ab||je>    1-4    0
!                       -> R3 _a,_b(j,e)abab = <ab||je>    3,4   0,1
!                       -> R4 _b,_a(j,e)abab = <ba||je>    2,4   0,1
!                       -> R5 _a,_b(j,e)abba = <ab||je>    3,4   0,1
!                       -> R6 _b,_a(j,e)abba = <ba||je>    2,4   0,1

use ccsd_global, only: dimm, Map_Type, mmul, noa, nob, nsym, nva, nvb
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize, ssn, key, aeqb
real(kind=wp), intent(inout) :: wrk(wrksize)
type(Map_Type), intent(in) :: n, r1, r2, r3, r4, r5, r6
integer(kind=iwp) :: dime, dimj, dimp, dimq, in_, inm, inp, ir1, ir2, ir3, ir4, ir5, ir6, lengthn, posn, posnm, posnp, posr1, &
                     posr2, posr3, posr4, posr5, posr6, symp, symq

do symp=1,nsym
  symq = mmul(ssn,symp)

  in_ = n%i(symp,1,1)
  lengthn = n%d(in_,2)

  if (lengthn == 0) cycle

  if (symp == symq) then
    ! symp == symq

    posn = n%d(in_,1)
    dimp = dimm(5,symp)

    !I.1 def R1 - _a_b(j,e)aaaa

    if ((key == 4) .and. (aeqb == 0)) then
      ir1 = r1%i(symp,1,1)
      posr1 = r1%d(ir1,1)
      dimj = noa(symp)
      dime = nva(symq)
      if (r1%d(ir1,2) > 0) call unpckhelp5(wrk(posn),wrk(posr1),dimp,dimj,dime,0,noa(symp),noa(symp),nva(symp))
    end if

    !I.2 def R2 - _a_b(j,e)bbbb

    if (aeqb == 0) then
      ir2 = r2%i(symp,1,1)
      posr2 = r2%d(ir2,1)
      dimj = nob(symp)
      dime = nvb(symq)
      if (r2%d(ir2,2) > 0) call unpckhelp5(wrk(posn),wrk(posr2),dimp,dimj,dime,0,nob(symp),nob(symp),nvb(symp))
    end if

    !I.3 def R3 - _a_b(j,e)abba

    if ((key == 2) .or. (key == 4)) then
      ir3 = r3%i(symp,1,1)
      posr3 = r3%d(ir3,1)
      dimj = nob(symp)
      dime = nva(symq)
      if (r3%d(ir3,2) > 0) call unpckhelp7(wrk(posn),wrk(posr3),dimp,dimp,dimj,dime,0,nob(symp),noa(symq),nva(symq))
    end if

    !I.4 def R4 - _b_a(j,e)abba

    if (((key == 3) .or. (key == 4)) .and. (aeqb == 0)) then
      ir4 = r4%i(symp,1,1)
      posr4 = r4%d(ir4,1)
      dimj = nob(symp)
      dime = nva(symq)
      if (r4%d(ir4,2) > 0) call unpckhelp6(wrk(posn),wrk(posr4),dimp,dimp,dimj,dime,0,nob(symp),noa(symq),nva(symq))
    end if

    !I.5 def R5 - _a_b (j,e)abab

    if ((key == 2) .or. (key == 4)) then
      ir5 = r5%i(symp,1,1)
      posr5 = r5%d(ir5,1)
      dimj = noa(symp)
      dime = nvb(symq)
      if (r5%d(ir5,2) > 0) call unpckhelp3(wrk(posn),wrk(posr5),dimp,dimp,dimj,dime,0,noa(symp),nob(symq),nvb(symq))
    end if

    !I.6 def R6 - _b_a(j,e)abab

    if (((key == 3) .or. (key == 4)) .and. (aeqb == 0)) then
      ir6 = r6%i(symp,1,1)
      posr6 = r6%d(ir6,1)
      dimj = noa(symp)
      dime = nvb(symq)
      if (r6%d(ir6,2) > 0) call unpckhelp4(wrk(posn),wrk(posr6),dimp,dimp,dimj,dime,0,noa(symp),nob(symq),nvb(symq))
    end if

  else
    ! symp /= symq

    inp = n%i(symp,1,1)
    inm = n%i(symq,1,1)
    posnp = n%d(inp,1)
    posnm = n%d(inm,1)
    dimp = dimm(5,symp)
    dimq = dimm(5,symq)

    !II.1 def R1 - _a_b(j,e)aaaa

    if (key == 4) then
      ir1 = r1%i(symp,1,1)
      posr1 = r1%d(ir1,1)
      dimj = noa(symp)
      dime = nva(symq)
      if (r1%d(ir1,2) > 0) call unpckhelp2(wrk(posnp),wrk(posnm),wrk(posr1),dimp,dimq,dimj,dime,0,noa(symp),noa(symq),nva(symq))
    end if

    !II.2 def R2 - _a_b(j,e)bbbb

    ir2 = r2%i(symp,1,1)
    posr2 = r2%d(ir2,1)
    dimj = nob(symp)
    dime = nvb(symq)
    if (r2%d(ir2,2) > 0) call unpckhelp2(wrk(posnp),wrk(posnm),wrk(posr2),dimp,dimq,dimj,dime,0,nob(symp),nob(symq),nvb(symq))

    !II.3 def R3 - _a_b(j,e)abba

    if ((key == 2) .or. (key == 4)) then
      ir3 = r3%i(symq,1,1)
      posr3 = r3%d(ir3,1)
      dimj = nob(symq)
      dime = nva(symp)
      if (r3%d(ir3,2) > 0) call unpckhelp7(wrk(posnp),wrk(posr3),dimp,dimq,dimj,dime,0,nob(symq),noa(symp),nva(symp))
    end if

    !II.4 def R4 - _b_a(j,e)abba

    if ((key == 3) .or. (key == 4)) then
      ir4 = r4%i(symp,1,1)
      posr4 = r4%d(ir4,1)
      dimj = nob(symp)
      dime = nva(symq)
      if (r4%d(ir4,2) > 0) call unpckhelp6(wrk(posnp),wrk(posr4),dimp,dimq,dimj,dime,0,nob(symp),noa(symq),nva(symq))
    end if

    !II.5 def R5 - _a_b(j,e)abab

    if ((key == 2) .or. (key == 4)) then
      ir5 = r5%i(symp,1,1)
      posr5 = r5%d(ir5,1)
      dimj = noa(symp)
      dime = nvb(symq)
      if (r5%d(ir5,2) > 0) call unpckhelp3(wrk(posnp),wrk(posr5),dimp,dimq,dimj,dime,0,noa(symp),nob(symq),nvb(symq))
    end if

    !II.6 def R6 - _b_a(j,e)abab

    if ((key == 3) .or. (key == 4)) then
      ir6 = r6%i(symq,1,1)
      posr6 = r6%d(ir6,1)
      dimj = noa(symq)
      dime = nvb(symp)
      if (r6%d(ir6,2) > 0) call unpckhelp4(wrk(posnp),wrk(posr6),dimp,dimq,dimj,dime,0,noa(symq),nob(symp),nvb(symp))
    end if

  end if

end do

return

end subroutine unpackab2

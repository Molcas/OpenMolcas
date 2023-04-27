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

subroutine unpackab1(wrk,wrksize,n,r1,r2,r3,r4,ssn,key,aeqb)
! n    - map type corresponding to N (Input)
! r1   - map type corresponding to R1 (Input)
! r2   - map type corresponding to R2 (Input)
! r3   - map type corresponding to R3 (Input)
! r4   - map type corresponding to R4 (Input)
! ssn  - overall symmetry state of N (Input)
! key  - define state of _a,_b  key   a, b (Input)
!        1    S  S
!        2    V  S
!        3    S  V
!        4    V  V
! aqeb - aqeb=0 if a /= b; =1 if a == b (Input)
!
! this routine unpacks: N _a_b (p,q) = <ab|pq>             key  aeqb
! N _a_b (p,q) a>=b-bb  -> R1 _a,_b(ef)aaaa = <ab||ef>      4    0
!                       -> R2 _a,_b(ef)bbbb = <ab||ef>     1-4   0
!                       -> R3 _a,_b(e,f)abab = <ab||ef>    3,4   0,1
!                       -> R4 _b,_a(e,f)abab = <ba||fe>    2,4   0,1
!
! !N.B. mylim, ze pri II.3;II.4;III.3;III.4 maju byt posn +--+ a nie ++++
! ako su terazky

use ccsd_global, only: dimm, Map_Type, mmul, noa, nob, nsym, nva, nvb
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize, ssn, key, aeqb
real(kind=wp), intent(inout) :: wrk(wrksize)
type(Map_Type), intent(in) :: n, r1, r2, r3, r4
integer(kind=iwp) :: dime, dimef, dimf, dimp, dimq, in_, inm, inp, ir1, ir2, ir3, ir4, lengthn, posn, posnm, posnp, posr1, posr2, &
                     posr3, posr4, symp, symq

do symp=1,nsym
  symq = mmul(ssn,symp)

  in_ = n%i(symp,1,1)
  lengthn = n%d(in_,2)

  if (lengthn == 0) cycle

  if (symp == symq) then
    ! symp == symq

    posn = n%d(in_,1)
    dimp = dimm(5,symp)

    !I.1 def R1 - _a_b(ef)aaaa

    if ((key == 4) .and. (aeqb == 0)) then
      ir1 = r1%i(symp,1,1)
      posr1 = r1%d(ir1,1)
      dimef = r1%d(ir1,2)
      if (r1%d(ir1,2) > 0) call unpckhelp1(wrk(posn),wrk(posr1),dimp,dimef,noa(symp),nva(symp))
    end if

    !I.2 def R2 - _a_b(ef)bbbb

    if (aeqb == 0) then
      ir2 = r2%i(symp,1,1)
      posr2 = r2%d(ir2,1)
      dimef = r2%d(ir2,2)
      if (r2%d(ir2,2) > 0) call unpckhelp1(wrk(posn),wrk(posr2),dimp,dimef,nob(symp),nvb(symp))
    end if

    !I.3  def R3 - _a_b(e,f)abab

    if ((key == 2) .or. (key == 4)) then
      ir3 = r3%i(symp,1,1)
      posr3 = r3%d(ir3,1)
      dime = dimm(3,symp)
      dimf = dimm(4,symq)
      if (r3%d(ir3,2) > 0) call unpckhelp3(wrk(posn),wrk(posr3),dimp,dimp,dime,dimf,noa(symp),nva(symp),nob(symq),nvb(symq))
    end if

    !I.4 def R4 - _b_a(e,f)abab

    if (((key == 3) .or. (key == 4)) .and. (aeqb == 0)) then
      ir4 = r4%i(symp,1,1)
      posr4 = r4%d(ir4,1)
      dime = dimm(3,symp)
      dimf = dimm(4,symq)
      if (r4%d(ir4,2) > 0) call unpckhelp4(wrk(posn),wrk(posr4),dimp,dimp,dime,dimf,noa(symp),nva(symp),nob(symq),nvb(symq))
    end if

  else if (symp > symq) then
    ! symp > symq

    inp = n%i(symp,1,1)
    inm = n%i(symq,1,1)
    posnp = n%d(inp,1)
    posnm = n%d(inm,1)
    dimp = dimm(5,symp)
    dimq = dimm(5,symq)

    !II.1 def R1 - _a_b(ef)aaaa

    if (key == 4) then
      ir1 = r1%i(symp,1,1)
      posr1 = r1%d(ir1,1)
      dime = dimm(3,symp)
      dimf = dimm(3,symq)
      if (r1%d(ir1,2) > 0) &
        call unpckhelp2(wrk(posnp),wrk(posnm),wrk(posr1),dimp,dimq,dime,dimf,noa(symp),nva(symp),noa(symq),nva(symq))
    end if

    !II.2 def R2 - _a_b(ef)bbbb

    ir2 = r2%i(symp,1,1)
    posr2 = r2%d(ir2,1)
    dime = dimm(4,symp)
    dimf = dimm(4,symq)
    if (r2%d(ir2,2) > 0) &
      call unpckhelp2(wrk(posnp),wrk(posnm),wrk(posr2),dimp,dimq,dime,dimf,nob(symp),nvb(symp),nob(symq),nvb(symq))

    !II.3 def R3 - _a_b(e,f)abab

    if ((key == 2) .or. (key == 4)) then
      ir3 = r3%i(symp,1,1)
      posr3 = r3%d(ir3,1)
      dime = dimm(3,symp)
      dimf = dimm(4,symq)
      if (r3%d(ir3,2) > 0) call unpckhelp3(wrk(posnp),wrk(posr3),dimp,dimq,dime,dimf,noa(symp),nva(symp),nob(symq),nvb(symq))
    end if

    !II.4 def R4 - _b_a(e,f)abab

    if ((key == 3) .or. (key == 4)) then
      ir4 = r4%i(symq,1,1)
      posr4 = r4%d(ir4,1)
      dime = dimm(3,symq)
      dimf = dimm(4,symp)
      if (r4%d(ir4,2) > 0) call unpckhelp4(wrk(posnp),wrk(posr4),dimp,dimq,dime,dimf,noa(symq),nva(symq),nob(symp),nvb(symp))
    end if

  else
    ! symp < symq

    inp = n%i(symp,1,1)
    inm = n%i(symq,1,1)
    posnp = n%d(inp,1)
    posnm = n%d(inm,1)
    dimp = dimm(5,symp)
    dimq = dimm(5,symq)

    !III.1 def R1 - _a_b(ef)aaaa -no operation
    !III.2 def R2 - _a_b(ef)bbbb -no operation

    !III.3 def R3 - _a_b(e,f)abab

    if ((key == 2) .or. (key == 4)) then
      ir3 = r3%i(symp,1,1)
      posr3 = r3%d(ir3,1)
      dime = dimm(3,symp)
      dimf = dimm(4,symq)
      if (r3%d(ir3,2) > 0) call unpckhelp3(wrk(posnp),wrk(posr3),dimp,dimq,dime,dimf,noa(symp),nva(symp),nob(symq),nvb(symq))
    end if

    !III.4 def R4 - _b_a(e,f)abab

    if ((key == 3) .or. (key == 4)) then
      ir4 = r4%i(symq,1,1)
      posr4 = r4%d(ir4,1)
      dime = dimm(3,symq)
      dimf = dimm(4,symp)
      if (r4%d(ir4,2) > 0) call unpckhelp4(wrk(posnp),wrk(posr4),dimp,dimq,dime,dimf,noa(symq),nva(symq),nob(symp),nvb(symp))
    end if

  end if

end do

return

end subroutine unpackab1

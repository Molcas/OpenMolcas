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

subroutine unpackab3(wrk,wrksize,mapdn,mapin,mapdr1,mapir1,ssn,nabnow,possab0,lentotab,key)
! mapdn    - direct map matrix corresponding to N (Input)
! mapin    - direct map matrix corresponding to N (Input)
! mapdr1   - inverse map matrix corresponding to R1 (Input)
! mapir1   - inverse map matrix corresponding to R1 (Input)
! ssn      - overall symmetry state of N (Input)
! nabnow   - #of b indices stored in stack (dimension of stack)
! possab0  - real position of AB_stack (I)
! lentotab - length of one record in stack (per one b ind) (I)
! key      - key, specifying which integrals need to be defined (I)
!
! this routine unpacks: N _a_b (p,q) = <ab|pq>             key
! N _a_b (p,q) a>=b-bb  -> R1 _a,_b(ef)aaaa = <ab||ef>      1
!                       -> R1 _a,_b(ef)bbbb = <ab||ef>      2
!                       -> R1 _a,_b(e,f)abab = <ab||ef>     3
!                       -> R1 _b,_a(e,f)abab = <ba||fe>     4
!
! !N.B. mylim, ze pri II.3;II.4;III.3;III.4 maju byt possn +--+ a nie ++++
! ako su terazky

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: wrksize, mapdn(0:512,6), mapin(8,8,8), mapdr1(0:512,6), mapir1(8,8,8), ssn, nabnow, possab0, lentotab, key
real(kind=wp) :: wrk(wrksize)
#include "ccsd1.fh"
integer(kind=iwp) :: bb, dime, dimef, dimf, dimp, dimq, in_, inm, inp, ir1, lengthn, possn, possnm, possnp, possr1, symp, symq

do bb=1,nabnow

  do symp=1,nsym
    symq = mmul(ssn,symp)

    in_ = mapin(symp,1,1)
    lengthn = mapdn(in_,2)

    if (lengthn == 0) cycle

    if (symp == symq) then
      ! symp == symq

      !old possn = mapdn(in_,1)
      possn = possab0+(mapdn(in_,1)-mapdn(1,1))+(bb-1)*lentotab
      dimp = dimm(5,symp)

      !I.1 def R1 - _a_b(ef)aaaa

      if (key == 1) then
        ir1 = mapir1(symp,1,1)
        possr1 = mapdr1(ir1,1)
        dimef = (dimm(3,symp)*(dimm(3,symp)-1))/2
        if (mapdr1(ir1,2) > 0) call unpckhelp8(wrk(possn),wrk(possr1),dimp,dimef,noa(symp),nva(symp),bb,nabnow)
      end if

      !I.2 def R1 - _a_b(ef)bbbb

      if (key == 2) then
        ir1 = mapir1(symp,1,1)
        possr1 = mapdr1(ir1,1)
        dimef = (dimm(4,symp)*(dimm(4,symp)-1))/2
        if (mapdr1(ir1,2) > 0) call unpckhelp8(wrk(possn),wrk(possr1),dimp,dimef,nob(symp),nvb(symp),bb,nabnow)
      end if

      !I.3 def R1 - _a_b(e,f)abab

      if (key == 3) then
        ir1 = mapir1(symp,1,1)
        possr1 = mapdr1(ir1,1)
        dime = dimm(3,symp)
        dimf = dimm(4,symq)
        if (mapdr1(ir1,2) > 0) call unpckhelp10(wrk(possn),wrk(possr1),dimp,dimp,dime,dimf,noa(symp),nva(symp),nob(symq), &
                                                nvb(symq),bb,nabnow)
      end if

      !I.4 def R1 - _b_a(e,f)abab

      if (key == 4) then
        ir1 = mapir1(symp,1,1)
        possr1 = mapdr1(ir1,1)
        dime = dimm(3,symp)
        dimf = dimm(4,symq)
        if (mapdr1(ir1,2) > 0) call unpckhelp11(wrk(possn),wrk(possr1),dimp,dimp,dime,dimf,noa(symp),nva(symp),nob(symq), &
                                                nvb(symq),bb,nabnow)
      end if

    else if (symp > symq) then
      ! symp > symq

      inp = mapin(symp,1,1)
      inm = mapin(symq,1,1)
      !old possnp = mapdn(inp,1)
      possnp = possab0+(mapdn(inp,1)-mapdn(1,1))+(bb-1)*lentotab
      !old possnm = mapdn(inm,1)
      possnm = possab0+(mapdn(inm,1)-mapdn(1,1))+(bb-1)*lentotab
      dimp = dimm(5,symp)
      dimq = dimm(5,symq)

      !II.1 def R1 - _a_b(ef)aaaa

      if (key == 1) then
        ir1 = mapir1(symp,1,1)
        possr1 = mapdr1(ir1,1)
        dime = dimm(3,symp)
        dimf = dimm(3,symq)
        if (mapdr1(ir1,2) > 0) call unpckhelp9(wrk(possnp),wrk(possnm),wrk(possr1),dimp,dimq,dime,dimf,noa(symp),nva(symp), &
                                               noa(symq),nva(symq),bb,nabnow)
      end if

      !II.2 def R1 - _a_b(ef)bbbb

      if (key == 2) then
        ir1 = mapir1(symp,1,1)
        possr1 = mapdr1(ir1,1)
        dime = dimm(4,symp)
        dimf = dimm(4,symq)
        if (mapdr1(ir1,2) > 0) call unpckhelp9(wrk(possnp),wrk(possnm),wrk(possr1),dimp,dimq,dime,dimf,nob(symp),nvb(symp), &
                                               nob(symq),nvb(symq),bb,nabnow)
      end if

      !II.3 def R1 - _a_b(e,f)abab

      if (key == 3) then
        ir1 = mapir1(symp,1,1)
        possr1 = mapdr1(ir1,1)
        dime = dimm(3,symp)
        dimf = dimm(4,symq)
        if (mapdr1(ir1,2) > 0) call unpckhelp10(wrk(possnp),wrk(possr1),dimp,dimq,dime,dimf,noa(symp),nva(symp),nob(symq), &
                                                nvb(symq),bb,nabnow)
      end if

      !II.4 def R4 - _b_a(e,f)abab

      if (key == 4) then
        ir1 = mapir1(symq,1,1)
        possr1 = mapdr1(ir1,1)
        dime = dimm(3,symq)
        dimf = dimm(4,symp)
        if (mapdr1(ir1,2) > 0) call unpckhelp11(wrk(possnp),wrk(possr1),dimp,dimq,dime,dimf,noa(symq),nva(symq),nob(symp), &
                                                nvb(symp),bb,nabnow)
      end if

    else
      ! symp < symq

      inp = mapin(symp,1,1)
      inm = mapin(symq,1,1)
      !old possnp = mapdn(inp,1)
      !old possnm = mapdn(inm,1)
      possnp = possab0+(mapdn(inp,1)-mapdn(1,1))+(bb-1)*lentotab
      possnm = possab0+(mapdn(inm,1)-mapdn(1,1))+(bb-1)*lentotab
      dimp = dimm(5,symp)
      dimq = dimm(5,symq)

      !III.1 def R1 - _a_b(ef)aaaa -no operation
      !III.2 def R1 - _a_b(ef)bbbb -no operation

      !III.3 def R1 - _a_b(e,f)abab

      if (key == 3) then
        ir1 = mapir1(symp,1,1)
        possr1 = mapdr1(ir1,1)
        dime = dimm(3,symp)
        dimf = dimm(4,symq)
        if (mapdr1(ir1,2) > 0) call unpckhelp10(wrk(possnp),wrk(possr1),dimp,dimq,dime,dimf,noa(symp),nva(symp),nob(symq), &
                                                nvb(symq),bb,nabnow)
      end if

      !III.4 def R3 - _b_a(e,f)abab

      if (key == 4) then
        ir1 = mapir1(symq,1,1)
        possr1 = mapdr1(ir1,1)
        dime = dimm(3,symq)
        dimf = dimm(4,symp)
        if (mapdr1(ir1,2) > 0) call unpckhelp11(wrk(possnp),wrk(possr1),dimp,dimq,dime,dimf,noa(symq),nva(symq),nob(symp), &
                                                nvb(symp),bb,nabnow)
      end if

    end if

  end do
end do

return

end subroutine unpackab3

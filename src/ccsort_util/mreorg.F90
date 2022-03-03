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

subroutine mreorg(wrk,wrksize,symp,symq,symr,typp,typq,typr,pospv2,posqv2,posrv2,typpv2,typqv2,typrv2,typv3,posv20,posv30,fact)
! this routine is up level routine for mreorg1 (also more detailed
! description can be found there).
! #2 must be of type 0, #3 can be 0, and 2
! this routine only prepares some constants, required by ireorg1,
! that can be deduced form input data - dimp,dimqr,dimt-dimv
!
! symp-r     - symmetries of p-r (I)
! typp-r     - types of indices p-r in V2 (I)
! posp-rv2   - positions of p-r ind. in V2 (I)
! typp-rv2   - types of indices, corresponding to p-r in V2 (I)
! typv3      - type of V3 (0,2) (I)
! posv20,30  - initial positions of V2 and V3 in wrk (I)
! fact       - multiplication factors (usually +-1.0) (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize, symp, symq, symr, typp, typq, typr, pospv2, posqv2, posrv2, typpv2, typqv2, typrv2, &
                                 typv3, posv20, posv30
real(kind=wp), intent(inout) :: wrk(wrksize)
real(kind=wp), intent(in) :: fact
integer(kind=iwp) :: dimp, dimqr, ind(4), mhelp, nhelp, rc

! define dimensions of V2

call ireorg2(symp,typpv2,nhelp,rc)
ind(pospv2) = nhelp
call ireorg2(symq,typqv2,nhelp,rc)
ind(posqv2) = nhelp
call ireorg2(symr,typrv2,nhelp,rc)
ind(posrv2) = nhelp

! def dimp,dimqr

call ireorg2(symp,typp,dimp,rc)

call ireorg2(symq,typq,nhelp,rc)
call ireorg2(symr,typr,mhelp,rc)

if ((typv3 == 2) .and. (symq == symr)) then
  dimqr = (nhelp*(nhelp-1))/2
else
  dimqr = nhelp*mhelp
end if

! use mreorg1

call mreorg1(symp,symq,symr,typp,typq,typr,pospv2,posqv2,posrv2,typpv2,typqv2,typrv2,typv3,wrk(posv20),wrk(posv30),fact,dimp, &
             dimqr,ind(1),ind(2),ind(3))

return

end subroutine mreorg

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

subroutine ireorg(wrk,wrksize,symp,symq,symr,syms,typp,typq,typr,typs,pospv1,posqv1,posrv1,possv1,typpv1,typqv1,typrv1,typsv1, &
                  typv2,posv10,posv20,fact)
! this routine is up level routine for ireorg1 (also more detailed
! description can be found there).
! v1 must be of type 0, v2 can be 0,1,3 and 4
! this routine only prepares some constants, required by ireorg1,
! that can be deduced form input data - dimpq,dimrs,dimt-dimx
!
! symp-s     - symmetries of p-s (I)
! typp-s     - types of indices p-s in V2 (I)
! posp-sv1   - positions of p-s ind. in v1 (I)
! typp-sv1   - types of indices, corresponding to p-s in V1 (I)
! typv2      - type of V2 (0,1,2,4) (I)
! posv10,20  - initial positions of V1 and V2 in wrk (I)
! fact       - multiplication factors (usually +-1.0) (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize, symp, symq, symr, syms, typp, typq, typr, typs, pospv1, posqv1, posrv1, possv1, typpv1, &
                                 typqv1, typrv1, typsv1, typv2, posv10, posv20
real(kind=wp), intent(inout) :: wrk(wrksize)
real(kind=wp), intent(in) :: fact
integer(kind=iwp) :: dimpq, dimrs, ind(4), mhelp = -1, nhelp = -1, rc

! define dimensions of V1

call ireorg2(symp,typpv1,nhelp,rc)
ind(pospv1) = nhelp
call ireorg2(symq,typqv1,nhelp,rc)
ind(posqv1) = nhelp
call ireorg2(symr,typrv1,nhelp,rc)
ind(posrv1) = nhelp
call ireorg2(syms,typsv1,nhelp,rc)
ind(possv1) = nhelp

! def dimpq,dimrs

call ireorg2(symp,typp,nhelp,rc)
call ireorg2(symq,typq,mhelp,rc)

if (((typv2 == 1) .or. (typv2 == 4)) .and. (symp == symq)) then
  dimpq = (nhelp*(nhelp-1))/2
else
  dimpq = nhelp*mhelp
end if

call ireorg2(symr,typr,nhelp,rc)
call ireorg2(syms,typs,mhelp,rc)

if (((typv2 == 3) .or. (typv2 == 4)) .and. (symr == syms)) then
  dimrs = (nhelp*(nhelp-1))/2
else
  dimrs = nhelp*mhelp
end if

! use ireorg1

call ireorg1(symp,symq,symr,syms,typp,typq,typr,typs,pospv1,posqv1,posrv1,possv1,typpv1,typqv1,typrv1,typsv1,typv2,wrk(posv10), &
             wrk(posv20),fact,dimpq,dimrs,ind(1),ind(2),ind(3),ind(4))

return

end subroutine ireorg

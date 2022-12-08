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

subroutine ireorg1(symp,symq,symr,syms,typp,typq,typr,typs,pospv1,posqv1,posrv1,possv1,typpv1,typqv1,typrv1,typsv1,typv2,v1, &
                   v2,fact,dimpq,dimrs,dimt,dimu,dimv,dimx)
! this routine maps v2(pq,rs) <+ fact . v1 (t,u,v,z)
! v2 may be of type 0,1,3,4 (typv2) while type v1 is always 0
! symp-syms and typp-typs are symmetries and types of p-s indices
! pospv1-possv1 are corresponding positions of p-s indices in v1
!
! symp-s     - symmetries of p-s (I)
! typp-s     - types of indices p-s in V2 (I)
! posp-sv1   - positions of p-s ind. in v1 (I)
! typp-sv1   - types of indices, corresponding to p-s in V1 (I)
! typv2      - type of V2 (0,1,2,4) (I)
! v1,v2      - arrays V1 and V2  (I,O)
! fact       - multiplication factors (usually +-1.0) (I)
! dimpq,rs   - dimensions of V2 (I)
! dimt-x     - dimensions of V1 (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: symp, symq, symr, syms, typp, typq, typr, typs, pospv1, posqv1, posrv1, possv1, typpv1, typqv1, &
                                 typrv1, typsv1, typv2, dimpq, dimrs, dimt, dimu, dimv, dimx
real(kind=wp), intent(in) :: v1(dimt,dimu,dimv,dimx), fact
real(kind=wp), intent(inout) :: v2(dimpq,dimrs)
integer(kind=iwp) :: ind(4), p, paddv1 = -1, pq, pqyes, pup = 0, q, qaddv1 = -1, qup = 0, r, raddv1 = -1, rc, rs, rsyes, rup = 0, &
                     s, saddv1 = -1, sup = 0

! def additive constants

call ireorg3(symp,typp,typpv1,paddv1,rc)
call ireorg3(symq,typq,typqv1,qaddv1,rc)
call ireorg3(symr,typr,typrv1,raddv1,rc)
call ireorg3(syms,typs,typsv1,saddv1,rc)

! def summation limits

call ireorg2(symp,typp,pup,rc)
call ireorg2(symq,typq,qup,rc)
call ireorg2(symr,typr,rup,rc)
call ireorg2(syms,typs,sup,rc)

! def pqyes, rsyes (i.e. if there is a reduced summation)

if ((typv2 == 1) .or. (typv2 == 4)) then
  if (symp == symq) then
    pqyes = 1
  else
    pqyes = 0
  end if
else
  pqyes = 0
end if

if ((typv2 == 3) .or. (typv2 == 4)) then
  if (symr == syms) then
    rsyes = 1
  else
    rsyes = 0
  end if
else
  rsyes = 0
end if

if ((pqyes == 1) .and. (rsyes == 1)) then

  ! case p>q, r>s

  rs = 0
  do r=2,rup
    ind(posrv1) = raddv1+r
    do s=1,r-1
      ind(possv1) = saddv1+s
      rs = rs+1

      pq = 0
      do p=2,pup
        ind(pospv1) = paddv1+p
        do q=1,p-1
          ind(posqv1) = qaddv1+q
          pq = pq+1

          v2(pq,rs) = v2(pq,rs)+fact*v1(ind(1),ind(2),ind(3),ind(4))

        end do
      end do
    end do
  end do

else if (pqyes == 1) then

  ! case p>q, r,s

  rs = 0
  do s=1,sup
    ind(possv1) = saddv1+s
    do r=1,rup
      ind(posrv1) = raddv1+r
      rs = rs+1

      pq = 0
      do p=2,pup
        ind(pospv1) = paddv1+p
        do q=1,p-1
          ind(posqv1) = qaddv1+q
          pq = pq+1

          v2(pq,rs) = v2(pq,rs)+fact*v1(ind(1),ind(2),ind(3),ind(4))

        end do
      end do
    end do
  end do

else if (rsyes == 1) then

  ! case p,q, r>s

  rs = 0
  do r=2,rup
    ind(posrv1) = raddv1+r
    do s=1,r-1
      ind(possv1) = saddv1+s
      rs = rs+1

      pq = 0
      do q=1,qup
        ind(posqv1) = qaddv1+q
        do p=1,pup
          ind(pospv1) = paddv1+p
          pq = pq+1

          v2(pq,rs) = v2(pq,rs)+fact*v1(ind(1),ind(2),ind(3),ind(4))

        end do
      end do
    end do
  end do

else

  ! case p,q, r,s

  rs = 0
  do s=1,sup
    ind(possv1) = saddv1+s
    do r=1,rup
      ind(posrv1) = raddv1+r
      rs = rs+1

      pq = 0
      do q=1,qup
        ind(posqv1) = qaddv1+q
        do p=1,pup
          ind(pospv1) = paddv1+p
          pq = pq+1

          v2(pq,rs) = v2(pq,rs)+fact*v1(ind(1),ind(2),ind(3),ind(4))

        end do
      end do
    end do
  end do

end if

return

end subroutine ireorg1

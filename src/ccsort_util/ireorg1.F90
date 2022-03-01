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

subroutine ireorg1(symp,symq,symr,syms,typp,typq,typr,typs,posspv1,possqv1,possrv1,posssv1,typpv1,typqv1,typrv1,typsv1,typv2,v1, &
                   v2,fact,dimpq,dimrs,dimt,dimu,dimv,dimx)
! this routine maps v2(pq,rs) <+ fact . v1 (t,u,v,z)
! v2 may be of type 0,1,3,4 (typv2) while type v1 is always 0
! symp-syms and typp-typs are symmetries and types of p-s indices
! posspv1-posssv1 are corresponding positions of p-s indices in v1
!
! symp-s     - symmetries of p-s (I)
! typp-s     - types of indices p-s in V2 (I)
! possp-sv1  - positions of p-s ind. in v1 (I)
! typp-sv1   - types of indices, corresponding to p-s in V1 (I)
! typv2      - type of V2 (0,1,2,4) (I)
! v1,v2      - arrays V1 and V2  (I,O)
! fact       - multiplication factors (usually +-1.0d0) (I)
! dimpq,rs   - dimensions of V2 (I)
! dimt-x     - dimensions of V1 (I)
!
! reorg.fh may not be included

#include "ccsort.fh"
integer symp, symq, symr, syms, typp, typq, typr, typs
integer posspv1, possqv1, possrv1, posssv1
integer typpv1, typqv1, typrv1, typsv1, typv2
integer dimpq, dimrs, dimt, dimu, dimv, dimx
real*8 v2(1:dimpq,1:dimrs)
real*8 v1(1:dimt,1:dimu,1:dimv,1:dimx)
real*8 fact
! help variables
integer p, q, r, s, pq, rs, rc, pqyes, rsyes
integer :: pup = 0, qup = 0, rup = 0, sup = 0
integer :: paddv1 = -1, qaddv1 = -1, raddv1 = -1, saddv1 = -1
integer ind(1:4)

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
    ind(possrv1) = raddv1+r
    do s=1,r-1
      ind(posssv1) = saddv1+s
      rs = rs+1

      pq = 0
      do p=2,pup
        ind(posspv1) = paddv1+p
        do q=1,p-1
          ind(possqv1) = qaddv1+q
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
    ind(posssv1) = saddv1+s
    do r=1,rup
      ind(possrv1) = raddv1+r
      rs = rs+1

      pq = 0
      do p=2,pup
        ind(posspv1) = paddv1+p
        do q=1,p-1
          ind(possqv1) = qaddv1+q
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
    ind(possrv1) = raddv1+r
    do s=1,r-1
      ind(posssv1) = saddv1+s
      rs = rs+1

      pq = 0
      do q=1,qup
        ind(possqv1) = qaddv1+q
        do p=1,pup
          ind(posspv1) = paddv1+p
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
    ind(posssv1) = saddv1+s
    do r=1,rup
      ind(possrv1) = raddv1+r
      rs = rs+1

      pq = 0
      do q=1,qup
        ind(possqv1) = qaddv1+q
        do p=1,pup
          ind(posspv1) = paddv1+p
          pq = pq+1

          v2(pq,rs) = v2(pq,rs)+fact*v1(ind(1),ind(2),ind(3),ind(4))

        end do
      end do
    end do
  end do

end if

return

end subroutine ireorg1

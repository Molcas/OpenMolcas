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

subroutine mreorg1(symp,symq,symr,typp,typq,typr,pospv1,posqv1,posrv1,typpv1,typqv1,typrv1,typv2,v1,v2,fact,dimp,dimqr,dimt, &
                   dimu,dimv)
! this routine maps v2(p,qr) <+ fact . v1 (t,u,v)
! v2 may be of type 0,2 (typv2) while type v1 is always 0
! symp-symr and typp-typr are symmetries and types of p-r indices
! pospv1-posrv1 are corresponding positions of p-r indices in v1
! N.B. v1 and v2 have no direct relation to #1 or #2, since they
! are not imported, v1,v2 corresponds to arbitrary matrices
!
! symp-r     - symmetries of p-r (I)
! typp-r     - types of indices p-r in V2 (I)
! posp-rv1   - positions of p-r ind. in v1 (I)
! typp-rv1   - types of indices, corresponding to p-r in V1 (I)
! typv2      - type of V2 (0,1,2,4) (I)
! v1,v2      - arrays V1 and V2  (I,O)
! fact       - multiplication factors (usually +-1.0) (I)
! dimp,qr    - dimensions of V2 (I)
! dimt-s     - dimensions of V1 (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: symp, symq, symr, typp, typq, typr, pospv1, posqv1, posrv1, typpv1, typqv1, typrv1, typv2, dimp, &
                                 dimqr, dimt, dimu, dimv
real(kind=wp), intent(in) :: v1(dimt,dimu,dimv), fact
real(kind=wp), intent(inout) :: v2(dimp,dimqr)
integer(kind=iwp) :: ind(4), p, paddv1, pup, q, qaddv1, qr, qryes, qup, r, raddv1, rc, rup

! def additional constants

call ireorg3(symp,typp,typpv1,paddv1,rc)
call ireorg3(symq,typq,typqv1,qaddv1,rc)
call ireorg3(symr,typr,typrv1,raddv1,rc)

! def sumation limits

call ireorg2(symp,typp,pup,rc)
call ireorg2(symq,typq,qup,rc)
call ireorg2(symr,typr,rup,rc)

! def qryes, rsyes (i.e. if there is a reduced sumations)

if (typv2 == 2) then
  if (symq == symr) then
    qryes = 1
  else
    qryes = 0
  end if
else
  qryes = 0
end if

if (qryes == 1) then

  ! case p, q>s

  qr = 0
  do q=2,qup
    ind(posqv1) = qaddv1+q
    do r=1,q-1
      ind(posrv1) = raddv1+r
      qr = qr+1

      do p=1,pup
        ind(pospv1) = paddv1+p

        v2(p,qr) = v2(p,qr)+fact*v1(ind(1),ind(2),ind(3))

      end do
    end do
  end do

else

  ! case p q,r

  qr = 0
  do r=1,rup
    ind(posrv1) = raddv1+r
    do q=1,qup
      ind(posqv1) = qaddv1+q
      qr = qr+1

      do p=1,pup
        ind(pospv1) = paddv1+p

        v2(p,qr) = v2(p,qr)+fact*v1(ind(1),ind(2),ind(3))

      end do
    end do
  end do

end if

return

end subroutine mreorg1

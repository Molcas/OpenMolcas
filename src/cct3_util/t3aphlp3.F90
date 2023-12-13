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

subroutine t3aphlp3(a1,a2,a3,b,dimp,dimq,dimr,dimqr,ns,szkey)
! this routine realizes following procedure
! for symp/=symq=symr
!
! B(p,qr)=B(p,qr)+ns.(A1(qr,p)-A2(p,r,q)+A3(p,q,r))
!
! a1    - A1 matrix (I)
! a2    - A2 matrix (I)
! a3    - A3 matrix (I)
! b     - B  matrix (I/O)
! dimp  - dimension of p index (I)
! dimq  - dimension of q index (I)
! dimr  - dimension of r index (I)
! dimqr - dimension of qr (I)
! ns    - signum of the permutation (+-1) (I)
! szkey - set zero key (I)
!         = 0 no vanishing
!         = 1 set B=0 at the beginning

use CCT3_global, only: nshf
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimp, dimq, dimr, dimqr, ns, szkey
real(kind=wp), intent(in) :: a1(dimqr,dimp), a2(dimp,dimr,dimq), a3(dimp,dimq,dimr)
real(kind=wp), intent(inout) :: b(dimp,dimqr)
integer(kind=iwp) :: q, qr, qr0

if (szkey == 1) b(:,:) = Zero

if (ns == 1) then
  ! phase +1

  do q=2,dimq
    qr0 = nshf(q)
    b(:,qr0+1:qr0+q-1) = b(:,qr0+1:qr0+q-1)+a3(:,q,1:q-1)
  end do

  do q=2,dimq
    qr0 = nshf(q)
    b(:,qr0+1:qr0+q-1) = b(:,qr0+1:qr0+q-1)-a2(:,1:q-1,q)
  end do

  do qr=1,dimqr
    b(:,qr) = b(:,qr)+a1(qr,:)
  end do

else
  ! phase -1

  do q=2,dimq
    qr0 = nshf(q)
    b(:,qr0+1:qr0+q-1) = b(:,qr0+1:qr0+q-1)-a3(:,q,1:q-1)
  end do

  do q=2,dimq
    qr0 = nshf(q)
    b(:,qr0+1:qr0+q-1) = b(:,qr0+1:qr0+q-1)+a2(:,1:q-1,q)
  end do

  do qr=1,dimqr
    b(:,qr) = b(:,qr)-a1(qr,:)
  end do

end if

return

end subroutine t3aphlp3

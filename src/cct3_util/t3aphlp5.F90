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

subroutine t3aphlp5(a1,a2,b,dimp,dimq,dimr,ns,szkey)
! this routine realizes following procedure
! for symp/=symq/=symr
!
! B(p,q,r)=B(p,q,r)+ns*(A1(q,r,p)-A2(p,r,q))
!
! a1    - A1 matrix (I)
! a2    - A2 matrix (I)
! b     - B  matrix (I/O)
! dimp  - dimension of p index (I)
! dimq  - dimension of q index (I)
! dimr  - dimension of r index (I)
! ns    - signum of the permutation (+-1) (I)
! szkey - set zero key (I)
!         = 0 no vanishing
!         = 1 set B=0 at the beginning

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimp, dimq, dimr, ns, szkey
real(kind=wp), intent(in) :: a1(dimq,dimr,dimp), a2(dimp,dimr,dimq)
real(kind=wp), intent(inout) :: b(dimp,dimq,dimr)
integer(kind=iwp) :: p, r

if (szkey == 1) b(:,:,:) = Zero

if (ns == 1) then
  ! phase +1

  do r=1,dimr
    b(:,:,r) = b(:,:,r)-a2(:,r,:)
  end do

  do p=1,dimp
    b(p,:,:) = b(p,:,:)+a1(:,:,p)
  end do

else
  ! phase -1

  do r=1,dimr
    b(:,:,r) = b(:,:,r)+a2(:,r,:)
  end do

  do p=1,dimp
    b(p,:,:) = b(p,:,:)-a1(:,:,p)
  end do

end if

return

end subroutine t3aphlp5

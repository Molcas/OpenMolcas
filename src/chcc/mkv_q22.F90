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

subroutine MkV_Q22(W2,W1,dima)
! this routine does
! W1(j,u,i,a') = 2W2(i,u,j,a')-W2(j,u,i,a')
!
! N.B. Kvajt odflaknute

use chcc_global, only: no
use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima
real(kind=wp), intent(in) :: W2(no,no,no,dima)
real(kind=wp), intent(out) :: W1(no,no,no,dima)
integer(kind=iwp) :: i, u

do i=1,no
  do u=1,no
    W1(:,u,i,:) = W2(:,u,i,:)-Two*W2(i,u,:,:)
  end do
end do

return

end subroutine MkV_Q22

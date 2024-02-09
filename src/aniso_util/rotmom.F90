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

subroutine rotmom(MOM,N,R,MOMR)
! inverse rotation

use Constants, only: cZero, cOne
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N
complex(kind=wp), intent(in) :: MOM(3,N,N) ! initial momentum matrix
real(kind=wp), intent(in) :: R(3,3) ! rotation matrix
complex(kind=wp), intent(out) :: MOMR(3,N,N) ! rotated momentum matrix
integer(kind=iwp) :: i, j, k
complex(kind=wp) :: RC(3,3)

! rotate the matrix

MOMR(:,:,:) = cZero

RC(:,:) = R(:,:)*cOne

do i=1,N
  do j=1,N
    do k=1,3
      MOMR(:,i,j) = MOMR(:,i,j)+RC(:,k)*MOM(k,i,j)
    end do
  end do
end do

return

end subroutine rotmom

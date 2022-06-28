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

subroutine daxpint(from,to,fact,ndim1,ndim2,ndim3,ndim4)
!bs subroutine similar to daxpy with interchange of two indices
!bs change from physicists notation to chemists notaion
!bs to(i,j,k,l)=to(i,j,k,l)+fact*from(i,k,j,l)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ndim1, ndim2, ndim3, ndim4
real(kind=wp), intent(in) :: from(ndim1,ndim2,ndim3,ndim4), fact
real(kind=wp), intent(inout) :: to(ndim1,ndim3,ndim2,ndim4)
integer(kind=iwp) :: irun2, irun3

if (fact /= Zero) then
  do irun3=1,ndim3
    do irun2=1,ndim2
      to(:,irun3,irun2,:) = to(:,irun3,irun2,:)+fact*from(:,irun2,irun3,:)
    end do
  end do
end if

return

end subroutine daxpint

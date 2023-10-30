!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine schmidtd2_cvb(c1,sxc1,nvec1,c2,nvec2,n)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nvec1, nvec2, n
real(kind=wp), intent(in) :: c1(n,nvec1), sxc1(n,nvec1)
real(kind=wp), intent(inout) :: c2(n,nvec2)
integer(kind=iwp) :: i, j
real(kind=wp), external :: ddot_

do i=1,nvec2
  do j=1,nvec1
    c2(:,i) = c2(:,i)-ddot_(n,c2(:,i),1,sxc1(:,j),1)*c1(:,j)
  end do
end do

return

end subroutine schmidtd2_cvb

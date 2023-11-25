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

subroutine saoon_cvb(c1,c2,n2,s,n,metr)
! Put matrix product S*C1 in C2:

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n2, n, metr
real(kind=wp), intent(in) :: c1(n,n2), s(*)
real(kind=wp), intent(out) :: c2(n,n2)
integer(kind=iwp) :: ik, j, k

if (metr == 0) then
  c2(:,:) = c1(:,:)
else if (metr == 1) then
  call mxatb_cvb(s,c1,n,n,n2,c2)
else if (metr == 2) then
  c2(:,:) = Zero
  do j=1,n2
    ik = 0
    do k=1,n
      c2(1:k-1,j) = c2(1:k-1,j)+s(ik+1:ik+k-1)*c1(k,j)
      c2(k,j) = c2(k,j)+sum(s(ik+1:ik+k)*c1(1:k,j))
      ik = ik+k
    end do
  end do
end if

return

end subroutine saoon_cvb

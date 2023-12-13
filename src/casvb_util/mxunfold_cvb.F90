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

subroutine mxunfold_cvb(avec,a,n)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: avec(n*(n-1))
real(kind=wp), intent(out) :: a(n,n)
integer(kind=iwp) :: i, iprm, j

a(:,:) = Zero
iprm = 0
do i=1,n
  do j=1,n
    if (j /= i) then
      iprm = iprm+1
      a(j,i) = avec(iprm)
    end if
  end do
end do

return

end subroutine mxunfold_cvb

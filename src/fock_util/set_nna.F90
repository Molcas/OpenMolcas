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

subroutine set_nnA(nSym,nAorb,nnA)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nSym, nAorb(8)
integer(kind=iwp), intent(inout) :: nnA(8,8)
integer(kind=iwp) :: i, j

do j=1,nSym
  do i=1,j-1
    nnA(i,j) = nAorb(i)*nAorb(j)
    nnA(j,i) = nnA(i,j)
  end do
  nnA(j,j) = nAorb(j)*(nAorb(j)+1)/2
end do

return

end subroutine set_nnA

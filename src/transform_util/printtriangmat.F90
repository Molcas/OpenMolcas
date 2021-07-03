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

subroutine PrintTriangMat(nRow,A)
! Prints a triangular matrix A(nRow,nRow)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nRow
real(kind=wp), intent(in) :: A(nRow*(nRow+1)/2)
integer(kind=iwp) :: i, iCount

if (nRow > 8) return
iCount = 1
do i=1,nRow
  write(u6,'(8F10.6)') A(iCount:iCount+i-1)
  iCount = iCount+i
end do

return

end subroutine PrintTriangMat

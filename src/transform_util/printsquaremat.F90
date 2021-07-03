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

subroutine PrintSquareMat(nRow,A)
! Prints a square matrix A(nRow,nRow)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nRow
real(kind=wp), intent(in) :: A(nRow*(nRow+1)/2)
integer(kind=iwp) :: i, iCount
character(len=8) :: Frmt

if (nRow > 8) return
write(Frmt,'("(",i1,"F10.6)")') nRow
iCount = 1
do i=1,nRow
  write(u6,Frmt) A(iCount:iCount+nRow-1)
  iCount = iCount+nRow
end do

return

end subroutine PrintSquareMat

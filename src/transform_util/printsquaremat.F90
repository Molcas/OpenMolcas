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

implicit real*8(a-h,o-z)
implicit integer(i-n)
dimension A(nRow**2)

if (nRow > 8) return
iCount = 0
do i=1,nRow
  if (nRow == 1) write(6,'(1F10.6)')(A(iCount+k),k=1,nRow)
  if (nRow == 2) write(6,'(2F10.6)')(A(iCount+k),k=1,nRow)
  if (nRow == 3) write(6,'(3F10.6)')(A(iCount+k),k=1,nRow)
  if (nRow == 4) write(6,'(4F10.6)')(A(iCount+k),k=1,nRow)
  if (nRow == 5) write(6,'(5F10.6)')(A(iCount+k),k=1,nRow)
  if (nRow == 6) write(6,'(6F10.6)')(A(iCount+k),k=1,nRow)
  if (nRow == 7) write(6,'(7F10.6)')(A(iCount+k),k=1,nRow)
  if (nRow == 8) write(6,'(8F10.6)')(A(iCount+k),k=1,nRow)
  iCount = iCount+nRow
end do

return

end subroutine PrintSquareMat

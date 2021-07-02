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

subroutine PrintDiagMat(nRow,A)
! Prints a diagonal matrix A(nRow,nRow)

implicit real*8(a-h,o-z)
implicit integer(i-n)
dimension A(nRow*(nRow+1))

if (nRow > 8) return
iCount = 0
do i=1,nRow
  write(6,'(8F10.6)')(A(iCount+k),k=1,i)
  iCount = iCount+i
end do

return

end subroutine PrintDiagMat

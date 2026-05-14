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

function FindDetR(matrix,n)
! Function to find the determinant of a square matrix
!
! Description: The Subroutine is based on two key points:
!
! 1]  A determinant is unaltered when row operations are performed:
!     Hence, using this principle, row operations (column operations
!     would work as well) are used to convert the matrix the matrix
!     into upper triangular form
! 2]  The determinant of a triangular matrix is obtained by finding
!     the product of the diagonal elements

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: FindDetR
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(inout) :: matrix(N,N)
integer(kind=iwp) :: i, j, k
real(kind=wp) :: m, temp
logical(kind=iwp) :: DetExists
real(kind=wp), parameter :: MINIMAL_REAL = tiny(MINIMAL_REAL)

DetExists = .true.
temp = 0
! Convert to upper triangular form
do k=1,N-1
  if (abs(matrix(k,k)) < MINIMAL_REAL) then
    DetExists = .false.
    do i=k+1,N
      if (abs(matrix(i,k)) > MINIMAL_REAL) then
        do j=1,N
          temp = matrix(i,j)
          matrix(i,j) = matrix(k,j)
          matrix(k,j) = temp
        end do
        DetExists = .true.
      end if
    end do
    if (.not. DetExists) then
      FindDetR = 0
      return
    end if
  end if
  do j=k+1,N
    m = matrix(j,k)/matrix(k,k)
    matrix(j,k+1:N) = matrix(j,k+1:N)-m*matrix(k,k+1:N)
  end do
end do ! k
! Evaluate determinant by finding product of diagonal elements
FindDetR = One
do i=1,N
  FindDetR = FindDetR*matrix(i,i)
end do

return

end function FindDetR

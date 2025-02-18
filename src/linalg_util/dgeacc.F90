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

!--------------------------------------------------------
! PURPOSE: Adds matrix B to matrix A element by element
!
! INPUT:
! A(*)    - Target matrix stored in column-major format
! LDA     - Leading dimension of matrix A
! B(*)    - Source matrix stored in column-major format
! LDB     - Leading dimension of matrix B
! M       - Number of rows in both matrices
! N       - Number of columns in both matrices
!
! OUTPUT:
! A(*)    - Updated matrix containing A + B
!--------------------------------------------------------
subroutine dgeAcc(A, LDA, B, LDB, M, N)
  use Definitions, only: wp, iwp

  implicit none
  real(kind=wp), intent(inout) :: A(*)
  real(kind=wp), intent(in) :: B(*)
  integer(kind=iwp), intent(in) :: LDA, LDB, M, N
  integer(kind=iwp) :: iCol, iRow

  do iRow = 0, m-1
    do iCol = 0, n-1
      a(iRow+iCol*lda+1) = a(iRow+iCol*lda+1)+b(iRow+iCol*ldb+1)
    end do
  end do

  return
end subroutine dgeAcc


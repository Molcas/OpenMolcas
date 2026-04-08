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
subroutine dgeAcc(A,LDA,B,LDB,M,N)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LDA, LDB, M, N
real(kind=wp), intent(inout) :: A(LDA,N)
real(kind=wp), intent(in) :: B(LDB,N)

A(1:M,:) = A(1:M,:)+B(1:M,:)

end subroutine dgeAcc

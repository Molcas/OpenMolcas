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
! Copyright (C) Anders Bernhardsson                                    *
!***********************************************************************

!--------------------------------------------------------
! PURPOSE: Adds matrix B to matrix A element by element
!
! INPUT:
! r       - Factor to multiply matrix A
! A(*)    - Source matrix stored in column-major format
! LDA     - Leading dimension of matrix A
! FORMA   - 'T': transpose A; 'N': do not transpose A
! B(*)    - Target matrix stored in column-major format
! LDB     - Leading dimension of matrix B
! M       - Number of rows in both matrices
! N       - Number of columns in both matrices
!
! OUTPUT:
! B(*)    - Updated matrix containing B + r*op(A)
!--------------------------------------------------------
subroutine dgeAcc(r,A,LDA,FORMA,B,LDB,M,N)
! MATRIX Addition (accumulation) FOR GENERAL MATRICES

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LDA, LDB, M, N
real(kind=wp), intent(in) :: r, A(LDA,*)
character, intent(in) :: FORMA
real(kind=wp), intent(inout) :: B(LDB,N)
integer(kind=iwp) :: iCol

if (FORMA == 'N') then
  b(1:m,:) = b(1:m,:)+r*a(1:m,1:n)
else if (FORMA == 'T') then
  do iCol=1,n
    b(1:m,iCol) = b(1:m,iCol)+r*a(iCol,1:m)
  end do
else
  write(u6,*) FORMA
  call Abend()
end if

end subroutine dgeAcc

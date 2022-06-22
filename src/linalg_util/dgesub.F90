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

subroutine dgeSub(A,LDA,FORMA,B,LDB,FORMB,C,LDC,M,N)
! MATRIX subtraction FOR GENERAL MATRICES

#include "intent.fh"

use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: A(*), B(*)
character, intent(in) :: FORMA, FORMB
integer(kind=iwp), intent(in) :: LDA, LDB, LDC, M, N
real(kind=wp), intent(_OUT_) :: C(*)
integer(kind=iwp) :: iCol, iRow

if ((FORMA == 'N') .and. (FORMB == 'N')) then
  do iRow=0,m-1
    do iCol=0,n-1
      c(iRow+iCol*ldc+1) = a(iRow+iCol*lda+1)-b(iRow+iCol*ldb+1)
    end do
  end do
else if ((FORMA == 'T') .and. (FORMB == 'N')) then
  do iRow=0,m-1
    do iCol=0,n-1
      c(iRow+iCol*ldc+1) = a(iCol+iRow*lda+1)-b(iRow+iCol*ldb+1)
    end do
  end do
else if ((FORMA == 'N') .and. (FORMB == 'T')) then
  do iRow=0,m-1
    do iCol=0,n-1
      c(iRow+iCol*ldc+1) = a(iRow+iCol*lda+1)-b(iCol+iRow*ldb+1)
    end do
  end do
else if ((FORMA == 'T') .and. (FORMB == 'T')) then
  do iRow=0,m-1
    do iCol=0,n-1
      c(iRow+iCol*ldc+1) = a(iCol+iRow*lda+1)-b(iCol+iRow*ldb+1)
    end do
  end do
else
  write(u6,*) 'Error when calling DGESUB, forma=',FORMA,'   formb=',FORMB
  call Abend()
end if

return

end subroutine dgeSub

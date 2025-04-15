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

subroutine dgeAdd2(r,A,LDA,FORMA,B,LDB,FORMB,C,LDC,M,N)
! MATRIX Addition FOR GENERAL MATRICES

use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: r, A(*), B(*), C(*)
integer(kind=iwp) :: LDA, LDB, LDC, M, N
character :: FORMA, FORMB
integer(kind=iwp) :: iCol

if ((FORMA == 'N') .and. (FORMB == 'N')) then
  do iCol=0,n-1
    c(iCol*ldc+1:iCol*ldc+m) = r*a(iCol*lda+1:iCol*lda+m)+b(iCol*ldb+1:iCol*ldb+m)
  end do
else if ((FORMA == 'T') .and. (FORMB == 'N')) then
  do iCol=0,n-1
    c(iCol*ldc+1:iCol*ldc+m) = r*a(iCol+1:iCol+m*lda:lda)+b(iCol*ldb+1:iCol*ldb+m)
  end do
else if ((FORMA == 'N') .and. (FORMB == 'T')) then
  do iCol=0,n-1
    c(iCol*ldc+1:iCol*ldc+m) = r*a(iCol*lda+1:iCol*lda+m)+b(iCol+1:iCol+m*ldb:ldb)
  end do
else if ((FORMA == 'T') .and. (FORMB == 'T')) then
  do iCol=0,n-1
    c(iCol*ldc+1:iCol*ldc+m) = r*a(iCol+1:iCol+m*lda:lda)+b(iCol+1:iCol+m*ldb:ldb)
  end do
else
  write(u6,*) FORMA,FORMB
  call Abend()
end if

return

end subroutine dgeAdd2

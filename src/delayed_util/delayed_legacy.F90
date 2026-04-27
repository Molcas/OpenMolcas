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

subroutine dgetf2(m,n,a,lda,ipiv,info)
use link_blas, only: lb_dgetf2
use Definitions, only: BLASInt, BLASR8
implicit none
integer(kind=BLASInt) :: info, lda, m, n
integer(kind=BLASInt) :: ipiv(*)
real(kind=BLASR8) :: a(lda,*)
call lb_dgetf2(m,n,a,lda,ipiv,info)
end subroutine dgetf2

subroutine dpotf2(uplo,n,a,lda,info)
use link_blas, only: lb_dpotf2
use Definitions, only: BLASInt, BLASR8
implicit none
character :: uplo
integer(kind=BLASInt) :: info, lda, n
real(kind=BLASR8) :: a(lda,*)
call lb_dpotf2(uplo,n,a,lda,info)
end subroutine dpotf2


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

subroutine dgetf2( m, n, a, lda, ipiv, info )
  use link_blas
  implicit none
  integer :: info, lda, m, n
  integer :: ipiv( * )
  real*8 :: a( lda, * )
  call lb_dgetf2( m, n, a, lda, ipiv, info )
end subroutine dgetf2

subroutine dpotf2( uplo, n, a, lda, info )
  use link_blas
  implicit none
  character :: uplo
  integer :: info, lda, n
  real*8 :: a( lda, * )
  call lb_dpotf2( uplo, n, a, lda, info )
end subroutine dpotf2


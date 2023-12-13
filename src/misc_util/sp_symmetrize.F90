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
! Copyright (C) 2014, Ignacio Fdez. Galvan                             *
!***********************************************************************
!  Sp_Symmetrize
!
!> @ingroup Sparse
!> @brief Converts a sparse matrix to symmetric format
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Converts a matrix \p A to a matrix \p B, stored in symmetric mode.
!> Only the lower triangle of \p A is stored.
!>
!> @param[in]  n   Size of the matrix
!> @param[in]  A   Input matrix, in sparse format
!> @param[in]  ija Index vector of matrix \p A
!> @param[out] B   Output matrix, in sparse format
!> @param[out] ijb Index vector of matrix \p B
!***********************************************************************

subroutine Sp_Symmetrize(n,A,ija,B,ijb)

use Constants, only: One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: n, ija(*)
real(kind=wp), intent(in) :: A(*)
real(kind=wp), intent(_OUT_) :: B(*)
integer(kind=iwp), intent(_OUT_) :: ijb(*)
integer(kind=iwp) :: i, j, k, nijb

ijb(1) = n+2
nijb = n+1
do i=1,n
  B(i) = A(i)
  do k=ija(i),ija(i+1)-1
    j = ija(k)
    if (j < i) then
      nijb = nijb+1
      B(nijb) = A(k)
      ijb(nijb) = ija(k)
    end if
  end do
  ijb(i+1) = nijb+1
end do
B(n+1) = One

end subroutine Sp_Symmetrize

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
!  Sp_Unpack
!
!> @ingroup Sparse
!> @brief
!>   Converts a sparse matrix to a dense (conventional) one
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Unpacks the sparse matrix \p Sp into a conventional dense matrix \p A.
!>
!> @param[in]  n     Size of the matrix
!> @param[in]  Sp    Input matrix, in sparse format
!> @param[in]  ij_Sp Index vector of matrix \p Sp
!> @param[out] A     Output matrix, in dense (conventional) format
!***********************************************************************

subroutine Sp_Unpack(n,Sp,ij_Sp,A)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n, ij_Sp(*)
real(kind=wp), intent(in) :: Sp(*)
real(kind=wp), intent(out) :: A(n,n)
integer(kind=iwp) :: i, j, k

A(:,:) = Zero
if (Sp(n+1) > Zero) then
  do i=1,n
    A(i,i) = Sp(i)
    do k=ij_Sp(i),ij_Sp(i+1)-1
      j = ij_Sp(k)
      A(i,j) = Sp(k)
      A(j,i) = Sp(k)
    end do
  end do
else
  do i=1,n
    A(i,i) = Sp(i)
    do k=ij_Sp(i),ij_Sp(i+1)-1
      j = ij_Sp(k)
      A(i,j) = Sp(k)
    end do
  end do
end if

end subroutine Sp_Unpack

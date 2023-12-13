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
!  Sp_Transpose
!
!> @ingroup Sparse
!> @brief
!>   Transpose a matrix in sparse format
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Saves in \p B the transpose of the input matrix \p A, both in sparse format.
!>
!> @param[in]  n   Size of the matrices
!> @param[in]  A   Input matrix, in sparse format
!> @param[in]  ija Index vector of matrix \p A
!> @param[out] B   Output matrix, in sparse format
!> @param[out] ijb Index vector of matrix \p B
!> @param[in]  nij Length of the vectors
!***********************************************************************

subroutine Sp_Transpose(n,A,ija,B,ijb,nij)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: n, ija(*), nij
real(kind=wp), intent(in) :: A(*)
real(kind=wp), intent(_OUT_) :: B(*)
integer(kind=iwp), intent(_OUT_) :: ijb(*)
integer(kind=iwp) :: i, j, k, kk
integer(kind=iwp), allocatable :: ia(:)

if (A(n+1) > Zero) then
  B(1:nij) = A(1:nij)
  ijb(1:nij) = ija(1:nij)
else
  call mma_allocate(ia,nij)

  ! Create an index of the rows in A
  do i=1,n
    B(i) = A(i)
    do k=ija(i),ija(i+1)-1
      ia(k) = i
    end do
  end do

  ! Lookup each column in A, save in B as a row
  ijb(1) = n+2
  kk = ijb(1)
  do j=1,n
    do k=ija(1),nij
      if (ija(k) == j) then
        ijb(kk) = ia(k)
        B(kk) = A(k)
        kk = kk+1
      end if
    end do
    ijb(j+1) = kk
  end do
  B(n+1) = Zero
  call mma_deallocate(ia)
end if

end subroutine Sp_Transpose

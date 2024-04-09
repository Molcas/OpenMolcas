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
! Copyright (C) 2019, Stefano Battaglia                                *
!***********************************************************************

subroutine transmat(A,U,N)
!**
! This subroutine carries out the following transformation:
!       U^T * A * U
! where A and U are NxN matrices
!**

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(inout) :: A(N,N)
real(kind=wp), intent(in) :: U(N,N)
real(kind=wp), allocatable :: B(:,:)

! Allocate temporary array B
call mma_allocate(B,N,N,Label='B')

! B = U^T * A
call dgemm_('T','N',N,N,N,One,U,N,A,N,Zero,B,N)
! A = B * U
call dgemm_('N','N',N,N,N,One,B,N,U,N,Zero,A,N)

call mma_deallocate(B)

end subroutine transmat

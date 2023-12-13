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

subroutine MatInvert(A,n)
! Inverts a square matrix

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(inout) :: A(n,n)
integer(kind=iwp) :: err, nw
real(kind=wp) :: dum(1)
integer(kind=iwp), allocatable :: ipiv(:)
real(kind=wp), allocatable :: wrk(:)

call mma_allocate(ipiv,n)
call dGeTRF_(n,n,A,n,ipiv,err)
call dGeTRI_(n,A,n,ipiv,dum,-1,err)
nw = int(dum(1))
call mma_allocate(wrk,nw)
call dGeTRI_(n,A,n,ipiv,wrk,nw,err)
call mma_deallocate(ipiv)
call mma_deallocate(wrk)

return

end subroutine MatInvert

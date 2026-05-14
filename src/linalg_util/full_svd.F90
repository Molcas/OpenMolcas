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

subroutine full_svd(m,n,amat,umat,vmat,svals)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: m, n
real(kind=wp), intent(inout) :: amat(m,*)
real(kind=wp), intent(_OUT_) :: umat(m,*), vmat(n,*), svals(*)
integer(kind=iwp) :: info, lwork
real(kind=wp) :: wrk1_lapack(1)
real(kind=wp), allocatable :: lapckwrk(:)

! Note that dgesvd returns V**T, not V.
!write(u6,*) ' In full_svd. Calling dgesvd:'
call dgesvd_('A','A',m,n,amat,m,svals,umat,m,vmat,n,wrk1_lapack,-1,info)
!write(u6,*) ' full_svd back from dgesvd'
lwork = int(wrk1_lapack(1))
!write(u6,*) ' lwork:',lwork
call mma_allocate(lapckwrk,lwork,label='lapckwrk')
!write(u6,*) ' Calling dgesvd again:'
call dgesvd_('A','A',m,n,amat,m,svals,umat,m,vmat,n,lapckwrk,lwork,info)
!write(u6,*) ' full_svd back from dgesvd'
call mma_deallocate(lapckwrk)

return

end subroutine full_svd

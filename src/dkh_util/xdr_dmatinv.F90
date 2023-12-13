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

subroutine XDR_dmatinv(a,n)
! Invert a real square matrix

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(inout) :: a(n*n)
integer(kind=iwp) :: info1, info2
#ifdef _MOLCAS_MPP_
logical(kind=iwp) :: isCloneQ
#endif
integer(kind=iwp), allocatable :: piv(:)
real(kind=wp), allocatable :: Tmp(:)

#ifdef _MOLCAS_MPP_
call check_parallel_data(a,n*n,isCloneQ,'C')
#endif
call mma_allocate(piv,n,label='piv')
call mma_allocate(Tmp,n,label='tmp')
call dgetrf_(n,n,a(1),n,piv,info1)
call dgetri_(n,a(1),n,piv,Tmp,n,info2)
call mma_deallocate(piv)
call mma_deallocate(Tmp)
#ifdef _MOLCAS_MPP_
if (isCloneQ) call check_parallel_data(a,n*n,isCloneQ,'S')
#endif

return

end subroutine XDR_dmatinv

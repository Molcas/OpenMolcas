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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine mxdiag_cvb(a,eigval,n)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(inout) :: a(n,n)
real(kind=wp), intent(out) :: eigval(n)
integer(kind=iwp) :: ierr
real(kind=wp), allocatable :: tmp(:)

call mma_allocate(tmp,n*3,label='tmp')
call dsyev_('V','L',n,a,n,eigval,tmp,n*3,ierr)
call mma_deallocate(tmp)
if (ierr /= 0) then
  write(u6,*) ' Fatal error in mxdiag, ierr :',ierr
  call abend_cvb()
end if

return

end subroutine mxdiag_cvb

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

subroutine sorti_cvb(n,arrin)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: n
integer(kind=iwp), intent(inout) :: arrin(n)
integer(kind=iwp) :: i
integer(kind=iwp), allocatable :: indx(:), tmp(:)

call mma_allocate(indx,n,label='indx')
call sortindxi_cvb(n,arrin,indx)
call mma_allocate(tmp,n,label='tmp')
do i=1,n
  tmp(i) = arrin(indx(i))
end do
arrin(:) = tmp(:)
call mma_deallocate(indx)
call mma_deallocate(tmp)

return

end subroutine sorti_cvb

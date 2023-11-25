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

subroutine prgrad_cvb(grad,n)

use casvb_global, only: ipr, norb, nprorb
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: grad(n)
real(kind=wp), allocatable :: tmp(:,:)

if (ipr(3) < 2) return
call mma_allocate(tmp,norb,norb,label='tmp')
call mxunfold_cvb(grad,tmp,norb)
write(u6,'(/,a)') ' Orbital gradient :'
call mxprint_cvb(tmp,norb,norb,0)
if (n-nprorb > 0) then
  write(u6,'(a)') ' Structure coefficient gradient :'
  call mxprint_cvb(grad(nprorb+1),1,n-nprorb,0)
end if
call mma_deallocate(tmp)

return

end subroutine prgrad_cvb

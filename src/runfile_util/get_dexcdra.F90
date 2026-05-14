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

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

subroutine Get_dExcdRa(dExcdRa,ndExcdRa)

use stdalloc, only: mma_allocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp), allocatable, intent(out) :: dExcdRa(:)
integer(kind=iwp), intent(out) :: ndExcdRa
logical(kind=iwp) :: Found
character(len=*), parameter :: Label = 'dExcdRa'

call qpg_dArray(Label,Found,ndExcdRa)
if ((.not. Found) .or. (ndExcdRa == 0)) call SysAbendmsg('Get_dExcdRa','Did not find:',Label)
call mma_allocate(dExcdRa,ndExcdRa,label='dExcdRa')
call Get_dArray(Label,dExcdRa,ndExcdRa)

return

end subroutine Get_dExcdRa

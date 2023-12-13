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

subroutine citrans_free

use citrans, only: ndet_group, ncsf_group, ndoc_group, nsoc_group, spintabs, spintabs_free
use stdalloc, only: mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: i

call mma_deallocate(ndoc_group)
call mma_deallocate(nsoc_group)
call mma_deallocate(ndet_group)
call mma_deallocate(ncsf_group)
do i=lbound(spintabs,1),ubound(spintabs,1)
  call mma_deallocate(spintabs(i)%coef)
end do
call spintabs_free()

end subroutine citrans_free

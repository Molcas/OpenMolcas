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

subroutine ZTRNSF(N,UR,UI,AR,AI)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(in) :: UR(N,N), UI(N,N)
real(kind=wp), intent(inout) :: AR(N,N), AI(N,N)
real(kind=wp), allocatable :: CI(:,:), CR(:,:)

call mma_allocate(CR,N,N,Label='CR')
call mma_allocate(CI,N,N,Label='CI')
call DGEMM_('N','N',N,N,N,One,AR,N,UR,N,Zero,CR,N)
call DGEMM_('N','N',N,N,N,-One,AI,N,UI,N,One,CR,N)
call DGEMM_('N','N',N,N,N,One,AR,N,UI,N,Zero,CI,N)
call DGEMM_('N','N',N,N,N,One,AI,N,UR,N,One,CI,N)
call DGEMM_('T','N',N,N,N,One,UR,N,CR,N,Zero,AR,N)
call DGEMM_('T','N',N,N,N,One,UI,N,CI,N,One,AR,N)
call DGEMM_('T','N',N,N,N,One,UR,N,CI,N,Zero,AI,N)
call DGEMM_('T','N',N,N,N,-One,UI,N,CR,N,One,AI,N)
call mma_deallocate(CI)
call mma_deallocate(CR)

end subroutine ZTRNSF

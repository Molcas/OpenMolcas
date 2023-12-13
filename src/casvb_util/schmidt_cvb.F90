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

subroutine schmidt_cvb(c,nvec,sao,n,metr)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nvec, n, metr
real(kind=wp), intent(inout) :: c(n,nvec)
real(kind=wp), intent(in) :: sao(*)
real(kind=wp), allocatable :: c2(:,:)

if (metr == 0) then
  call schmidt2_cvb(c,c,nvec,sao,n,metr)
else
  call mma_allocate(c2,n,nvec,label='c2')
  call schmidt2_cvb(c,c2,nvec,sao,n,metr)
  call mma_deallocate(c2)
end if

return

end subroutine schmidt_cvb

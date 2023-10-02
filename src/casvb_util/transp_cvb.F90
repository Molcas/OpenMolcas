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

subroutine transp_cvb(a,b,n1,n2)
! Transposes matrix A; A and B may share memory.
! (FIXME: That goes against the Fortran standard)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: n1, n2
real(kind=wp) :: a(n1,n2), b(n2,n1)
integer(kind=iwp) :: i, j
real(kind=wp), allocatable :: tmp(:,:)

call mma_allocate(tmp,n2,n1,label='tmp')
do i=1,n1
  do j=1,n2
    tmp(j,i) = a(i,j)
  end do
end do
call fmove_cvb(tmp,b,n2*n1)
call mma_deallocate(tmp)

return

end subroutine transp_cvb

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

function detm_cvb(a,n)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Ten
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: detm_cvb
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: a(n,n)
integer(kind=iwp) :: ierr
real(kind=wp) :: det(2)
integer(kind=iwp), allocatable :: itmp(:)
real(kind=wp), allocatable :: tmp1(:,:), tmp2(:,:)

! start linpack_determinant
if (n == 0) then
  detm_cvb = One
  return
end if
call mma_allocate(tmp1,n,n,label='tmp1')
call mma_allocate(itmp,n,label='itmp')
ierr = 0
tmp1(:,:) = a(:,:)
call dgetrf_(n,n,tmp1,n,itmp,ierr)
! start linpack_determinant
!call dgefa(tmp1,n,n,itmp,ierr)
if (ierr /= 0) then
  detm_cvb = Zero
else
  call mma_allocate(tmp2,n,n,label='tmp2')
  call dgedi(tmp1,n,n,itmp,det,tmp2,10)
  detm_cvb = det(1)*Ten**det(2)
end if

call mma_deallocate(tmp1)
call mma_deallocate(itmp)
call mma_deallocate(tmp2)

return

end function detm_cvb

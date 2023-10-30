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

subroutine stringen_cvb(norb,nel,locc,lunocc)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: norb, nel
integer(kind=iwp), intent(inout) :: locc(*), lunocc(*)
integer(kind=iwp) :: i_locc, i_lunocc, indx, iorb, rc
integer(kind=iwp), allocatable :: nk(:), nkmax(:), nkmin(:)

call mma_allocate(nk,[0,norb],label='nk')
call mma_allocate(nkmax,[0,norb],label='nkmax')
call mma_allocate(nkmin,[0,norb],label='nkmin')
! Spin string loop initialization (use xdet as graph storage):
do iorb=0,norb
  nkmin(iorb) = max(iorb-norb+nel,0)
  nkmax(iorb) = min(iorb,nel)
end do
nk(:) = nkmax(:)
! Spin string loop starts here:
indx = 0
do
  indx = indx+1
  i_locc = (indx-1)*nel+1
  i_lunocc = (indx-1)*(norb-nel)+1
  call occupy_cvb(nk,norb,locc(i_locc),lunocc(i_lunocc))
  call loop_cvb(norb,nk,nkmin,nkmax,rc)
  if (rc == 0) exit
end do
call mma_deallocate(nk)
call mma_deallocate(nkmax)
call mma_deallocate(nkmin)

return

end subroutine stringen_cvb

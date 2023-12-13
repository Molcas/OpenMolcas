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

subroutine report_cvb(orbs,norb)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: norb
real(kind=wp), intent(in) :: orbs(norb,norb)
real(kind=wp), allocatable :: tmp(:,:)

write(u6,'(/,a)') ' Orbital coefficients :'
write(u6,'(a)') ' ----------------------'
call mxprint_cvb(orbs,norb,norb,0)
write(u6,'(/,a)') ' Overlap between orbitals :'
write(u6,'(a)') ' --------------------------'

call mma_allocate(tmp,norb,norb,label='tmp')
call mxattb_cvb(orbs,orbs,norb,norb,norb,tmp)
call mxprint_cvb(tmp,norb,norb,0)
call mma_deallocate(tmp)

return

end subroutine report_cvb

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

subroutine ao2mo_cvb(orbsao,orbs,norb1)

use casvb_global, only: nbas_mo, norb
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: norb1
real(kind=wp), intent(in) :: orbsao(nbas_mo,norb1)
real(kind=wp), intent(out) :: orbs(norb,norb1)
real(kind=wp), allocatable :: tmp(:,:)

if (norb1 == 0) return
call mma_allocate(tmp,nbas_mo,norb,label='tmp')
call getmo_cvb(tmp,3)
call mxattb_cvb(tmp,orbsao,norb,nbas_mo,norb1,orbs)
call mma_deallocate(tmp)

return

end subroutine ao2mo_cvb

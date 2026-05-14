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
! Copyright (C) 2021, Vladislav Kochetov                               *
!***********************************************************************

subroutine cut_matrices()
!***********************************************************************
! cut dynamics matrices from size (lrootstot,lrootstot) to size (d,d)
!***********************************************************************

use rhodyn_data, only: CSF2SO, d, density0, dipole_basis, dysamp_bas, flag_dyson, hamiltonian, istates, lrootstot, SO_CI, U_CI
use rhodyn_utils, only: removeColumn, removeLineAndColumn
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp

implicit none
complex(kind=wp), allocatable :: d1(:,:), d2(:,:), d3(:,:)

! cut dynamics matrices
call removeLineAndColumn(hamiltonian,istates)
call removeLineAndColumn(density0,istates)
if (flag_dyson) call removeLineAndColumn(dysamp_bas,istates)
! cut dipole moment
call mma_allocate(d1,lrootstot,lrootstot,label='d1')
call mma_allocate(d2,lrootstot,lrootstot,label='d2')
call mma_allocate(d3,lrootstot,lrootstot,label='d3')
d1(:,:) = dipole_basis(:,:,1)
d2(:,:) = dipole_basis(:,:,2)
d3(:,:) = dipole_basis(:,:,3)
call removeLineAndColumn(d1,istates)
call removeLineAndColumn(d2,istates)
call removeLineAndColumn(d3,istates)
call mma_deallocate(dipole_basis)
call mma_allocate(dipole_basis,d,d,3)
dipole_basis(:,:,1) = d1
dipole_basis(:,:,2) = d2
dipole_basis(:,:,3) = d3
call mma_deallocate(d1)
call mma_deallocate(d2)
call mma_deallocate(d3)

! cut transform marices
call removeColumn(U_CI,istates)
call removeColumn(CSF2SO,istates)
call removeLineAndColumn(SO_CI,istates)

end subroutine cut_matrices

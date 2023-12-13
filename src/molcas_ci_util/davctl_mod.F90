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

module davctl_mod
! Parameters and variables for the Davidson diagonalization scheme

use Definitions, only: wp, iwp

implicit none
private

! flags indicating the storage mode used by the Davidson diagonalization scheme
integer(kind=iwp), parameter :: in_core = 0 , on_disk = 2, mixed_mode_1 = 3, mixed_mode_2 = 4, &
                                llab = 16 ! label length
integer(kind=iwp) :: istart, mxDiskStk, mxMemStk, n_Roots, nDiskStk, nkeep, nMemStk, nvec, save_mode ! flag for the storage handling
logical(kind=iwp) :: save_in_memory
integer(kind=iwp), allocatable :: disk_address(:) ! disk address table
real(kind=wp), allocatable :: memory_vectors(:,:)
character(len=llab), allocatable :: LblStk(:) ! stack of labels

public :: disk_address, in_core, istart, llab, LblStk, memory_vectors, mixed_mode_1, mixed_mode_2, mxDiskStk, mxMemStk, n_Roots, &
          nDiskStk, nkeep, nMemStk, nvec, on_disk, save_in_memory, save_mode

end module davctl_mod

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

subroutine ALLOC_LUCIA()
! Dimensions and
! Allocation of static memory
!
! =====
! Input
! =====
!
! KFREE : Pointer to first element of free space
! Information in /LUCINP/,/ORBINP/,/CSYM/
!
! ======
! Output
! ======
! KFREE : First array of free space after allocation of static memory
! /GLBBAS/,/CDIM/
!
! =======
! Version
! =======
!
! Modified Jan 1997
!           Fall 97 (PGINT1 added)
!           Spring 99

use lucia_data, only: INT1, INT1O, KINH1, KINH1_NOCCSYM, LSM1, LSM2, NBINT1, NBINT2, NSMOB, NTOOB, PGINT1, PGINT1A, PINT1, PINT2, &
                      RHO1, SRHO1
use stdalloc, only: mma_allocate
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ISM

! 1 : One electron integrals (Complete matrix allocated)
call mma_allocate(INT1,NTOOB**2,Label='INT1')
! A copy of the original UNMODIFIED 1-elecs ints
call mma_allocate(INT1O,NTOOB**2,Label='Int1O')
! Zero to avoid problems with elements that will not
! be initialized
INT1(:) = Zero
INT1O(:) = Zero
! 2 : Two electron integrals
! Pointers to symmetry block of integrals
call mma_allocate(PINT1,NBINT1,Label='PINT1')
call mma_allocate(PINT2,NBINT2,Label='PINT2')
! Pointers to nonsymmetric one-electron integrals
do ISM=1,NSMOB
  ! triangular packed
  call mma_allocate(PGINT1(ISM)%A,NSMOB,Label='PGINT1(ISM)%I')
  ! no packing
  call mma_allocate(PGINT1A(ISM)%A,NSMOB,Label='PGINT1A(ISM)%I')
end do
! Symmetry of last index as a function of initial index
call mma_allocate(LSM1,NBINT1,Label='LSM1')
call mma_allocate(LSM2,NBINT2,Label='LSM2')
! 3 One-body density
call mma_allocate(RHO1,NTOOB**2,Label='RHO1')
! 3.1 : One-body spin density
call mma_allocate(SRHO1,NTOOB**2,Label='SRHO1')
! indices for pair of orbitals symmetry ordered
! Lower half packed
call mma_allocate(KINH1,NTOOB*NTOOB,Label='KINTH1')
! Complete form
call mma_allocate(KINH1_NOCCSYM,NTOOB*NTOOB,Label='KINTH1_NOCCSYM')

end subroutine ALLOC_LUCIA

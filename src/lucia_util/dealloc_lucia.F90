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

subroutine DEALLOC_LUCIA()
! Deallocate memory allocated during alloc_lucia

use stdalloc, only: mma_deallocate
use GLBBAS, only: INT1, INT1O, PINT1, PINT2, PGINT1, PGINT1A, LSM1, LSM2, RHO1, SRHO1, KINH1_NOCCSYM, KINH1, CI_VEC, SIGMA_VEC
use lucia_data, only: MXSOOB, XISPSM
use lucia_data, only: LCSBLK
use lucia_data, only: IREFSM, PSSIGN
use lucia_data, only: NSMOB
use Constants, only: Zero, Two

implicit none
! Input
integer ISM, LBLOCK

! 1 : One electron integrals( Complete matrix allocated )
call mma_deallocate(INT1)
! A copy of the original UNMODIFIED 1-elecs ints
call mma_deallocate(INT1O)
! 2 : Two electron integrals
!  Pointers to symmetry block of integrals
call mma_deallocate(PINT1)
call mma_deallocate(PINT2)
!  Pointers to nonsymmetric one-electron integrals
do ISM=1,NSMOB
  ! triangular packed
  call mma_deallocate(PGINT1(ISM)%I)
  ! no packing
  call mma_deallocate(PGINT1A(ISM)%I)
end do
! Symmetry of last index as a function of initial index
call mma_deallocate(LSM1)
call mma_deallocate(LSM2)
! 3 One-body density
call mma_deallocate(RHO1)
! 3.1 : One-body spin density
call mma_deallocate(SRHO1)
! indices for pair of orbitals symmetry ordered
! Lower half packed
call mma_deallocate(KINH1)
! Complete form
call mma_deallocate(KINH1_NOCCSYM)

! arrays allocated at end of lucia
LBLOCK = MXSOOB
LBLOCK = max(LBLOCK,LCSBLK)
LBLOCK = max(int(XISPSM(IREFSM,1)),MXSOOB)
if (PSSIGN /= Zero) LBLOCK = int(Two*XISPSM(IREFSM,1))
call mma_deallocate(CI_VEC)
call mma_deallocate(SIGMA_VEC)

end subroutine DEALLOC_LUCIA

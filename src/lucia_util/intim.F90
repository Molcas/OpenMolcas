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
! Copyright (C) 1991,1997, Jeppe Olsen                                 *
!***********************************************************************

subroutine INTIM()
! Interface to external integrals
!
! If NOINT /= 0, only pointers are constructed
! Jeppe Olsen, Winter of 1991
!
! Version : Fall 97

use lucia_data, only: ECORE, ECORE_HEX, ECORE_ORIG, INT1, INT1O, KINH1, KINH1_NOCCSYM, LSM1, LSM2, NOINT, NSMOB, NTOOB, NTOOBS, &
                      PINT1, PINT2
use Constants, only: Zero

implicit none

! Pointers for symmetry blocks of integrals
call INTPNT(PINT1,LSM1,PINT2,LSM2)

! Pointer for orbital indices for symmetry blocked matrices
call ORBINH1(KINH1,KINH1_NOCCSYM,NTOOBS,NTOOB,NSMOB)

! Change one-electron integrals to inactive fock matrix
if (NOINT == 0) then
  INT1O(:) = INT1(:)
  ECORE_HEX = Zero
end if
ECORE_ORIG = ECORE
ECORE = ECORE+ECORE_HEX

end subroutine INTIM

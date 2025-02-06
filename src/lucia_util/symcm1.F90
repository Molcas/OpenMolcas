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

subroutine SYMCM1(ITASK,I1,I2,I12)
! Symmetries I1,I2,I12 are related as
! I1*I2 = I12
! IF (ITASK == 1) I2 and I12 are known, find I1
! IF (ITASK == 2) I1 and I12 are known, find I1
! IF (ITASK == 3) I1 and I2 are known, find I12
!
! D2h version, written for compatibility with general symmetry

use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ITASK, I1, I2, I12

if (ITASK == 1) then
  I1 = Mul(I2,I12)

  !if ((I12 > 8) .or. (I2 > 8) .or. (I12 <= 0) .or. (I2 <= 0)) write(u6,*) ' I12 and I2 = ',I12,I2

else if (ITASK == 2) then

  !if ((I12 > 8) .or. (I1 > 8) .or. (I12 <= 0) .or. (I1 <= 0)) write(u6,*) ' I12 and I1 = ',I12,I1

  I2 = Mul(I1,I12)
else if (ITASK == 3) then

  !if ((I2 > 8) .or. (I1 > 8) .or. (I2 <= 0) .or. (I1 <= 0)) write(u6,*) ' I2 and I1 = ',I2,I1

  I12 = Mul(I1,I2)

end if

end subroutine SYMCM1

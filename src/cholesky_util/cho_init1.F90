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

subroutine CHO_INIT1()
!
! Purpose: initialize counter arrays.

use Cholesky, only: Cho_Real_Par, InfRed, InfVec, myNumCho, nSym, NumCho, NumChT, RstCho, XnPass

implicit none

if (RSTCHO) then

  ! Read restart information.
  ! -------------------------

  call CHO_GETRSTC()
  NUMCHT = sum(NUMCHO(1:NSYM))

else

  ! Initialize vector info and counters.
  ! ------------------------------------

  INFVEC(:,:,:) = 0
  NUMCHO(1:NSYM) = 0
  NUMCHT = 0

  ! Initialize reduced set info.
  ! ----------------------------

  INFRED(:) = 0

  ! Initialize global integral pass counter.
  ! ----------------------------------------

  XNPASS = 0

end if

! Parallel init.
! --------------

if (Cho_Real_Par) MYNUMCHO(1:NSYM) = 0

end subroutine CHO_INIT1

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

use ChoSwp, only: InfRed, InfVec
use Definitions, only: iwp

implicit none
#include "cholesky.fh"
#include "choglob.fh"
#include "cho_para_info.fh"
integer(kind=iwp), external :: CHO_ISUMELM

if (RSTCHO) then

  ! Read restart information.
  ! -------------------------

  call CHO_GETRSTC()
  NUMCHT = CHO_ISUMELM(NUMCHO,NSYM)

else

  ! Initialize vector info and counters.
  ! ------------------------------------

  call IZERO(INFVEC,size(INFVEC))
  call IZERO(NUMCHO,NSYM)
  NUMCHT = 0

  ! Initialize reduced set info.
  ! ----------------------------

  call IZERO(INFRED,size(INFRED))

  ! Initialize global integral pass counter.
  ! ----------------------------------------

  XNPASS = 0

end if

! Parallel init.
! --------------

if (Cho_Real_Par) call IZERO(MYNUMCHO,NSYM)

end subroutine CHO_INIT1

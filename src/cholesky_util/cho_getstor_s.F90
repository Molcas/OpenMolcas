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

subroutine CHO_GETSTOR_S(VCSTOR,ISYM)
!
! Purpose: get total vector storage (in words), symmetry ISYM.

use ChoArr, only: nDimRS
use ChoSwp, only: InfVec

implicit real*8(a-h,o-z)
#include "cholesky.fh"

if (NUMCHO(ISYM) < 1) then
  VCSTOR = 0.0d0
else if (.not. allocated(nDimRS)) then
  IRED = INFVEC(NUMCHO(ISYM),2,ISYM)
  JRED = 3
  call CHO_GETRED(IRED,JRED,.false.)
  call CHO_SETREDIND(JRED)
  VCSTOR = dble(INFVEC(NUMCHO(ISYM),4,ISYM))+dble(NNBSTR(ISYM,JRED))
else
  IRED = INFVEC(NUMCHO(ISYM),2,ISYM)
  VCSTOR = dble(INFVEC(NUMCHO(ISYM),4,ISYM))+dble(NDIMRS(ISYM,IRED))
end if

end subroutine CHO_GETSTOR_S

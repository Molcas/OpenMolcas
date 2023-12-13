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

subroutine CHO_DBGINT()
!
! Purpose: regenerate and check integrals as specified in input.

use Cholesky, only: CHO_MINCHK, ICHKQ, NCHKQ

implicit none

if (CHO_MINCHK) then ! minimal check
  call CHO_MCA_DBGINT_S(ICHKQ,NCHKQ,.true.)
else ! check all (or the number of col. from input)
  call CHO_MCA_DBGINT_A()
end if

end subroutine CHO_DBGINT

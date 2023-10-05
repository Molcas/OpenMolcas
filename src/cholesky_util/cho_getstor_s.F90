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

use Cholesky, only: InfVec, nDimRS, nnBstR, NumCho
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: VCSTOR
integer(kind=iwp), intent(in) :: ISYM
integer(kind=iwp) :: IRED, JRED

if (NUMCHO(ISYM) < 1) then
  VCSTOR = Zero
else if (.not. allocated(nDimRS)) then
  IRED = INFVEC(NUMCHO(ISYM),2,ISYM)
  JRED = 3
  call CHO_GETRED(IRED,JRED,.false.)
  call CHO_SETREDIND(JRED)
  VCSTOR = real(INFVEC(NUMCHO(ISYM),4,ISYM),kind=wp)+real(NNBSTR(ISYM,JRED),kind=wp)
else
  IRED = INFVEC(NUMCHO(ISYM),2,ISYM)
  VCSTOR = real(INFVEC(NUMCHO(ISYM),4,ISYM),kind=wp)+real(NDIMRS(ISYM,IRED),kind=wp)
end if

end subroutine CHO_GETSTOR_S

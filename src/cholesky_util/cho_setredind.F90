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

subroutine CHO_SETREDIND(IRED)
!
! Purpose: set index arrays for reduced set IRED. The counter
!          array NNBSTRSH must be set on entry.

use Cholesky, only: iiBstR, iiBstRSh, nnBstR, nnBstRSh, nnBstRT, nnShl, nSym
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: IRED
integer(kind=iwp) :: ISHLAB, ISYM, J
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: mmShl, mSym
#endif
character(len=*), parameter :: SECNAM = 'CHO_SETREDIND'

J = IRED

#ifdef _DEBUGPRINT_
MSYM = size(iiBstRSh,1)
MMSHL = size(iiBstRSh,2)
if ((NNSHL /= MMSHL) .or. (NSYM /= MSYM)) call CHO_QUIT('[1] Dimension error in '//SECNAM,104)
if ((J < 1) .or. (J > 3)) call CHO_QUIT('[2] Dimension error in '//SECNAM,104)
#endif

if (NNSHL < 1) then ! may occur in parallel runs
  NNBSTRT(J) = 0
  IIBSTR(1:NSYM,J) = 0
  NNBSTR(1:NSYM,J) = 0
  return
end if

NNBSTRT(J) = 0
do ISYM=1,NSYM
  IIBSTRSH(ISYM,1,J) = 0
  NNBSTR(ISYM,J) = NNBSTRSH(ISYM,1,J)
  do ISHLAB=2,NNSHL
    IIBSTRSH(ISYM,ISHLAB,J) = NNBSTR(ISYM,J)
    NNBSTR(ISYM,J) = NNBSTR(ISYM,J)+NNBSTRSH(ISYM,ISHLAB,J)
  end do
  IIBSTR(ISYM,J) = NNBSTRT(J)
  NNBSTRT(J) = NNBSTRT(J)+NNBSTR(ISYM,J)
end do

end subroutine CHO_SETREDIND

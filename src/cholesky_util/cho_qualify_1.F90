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

subroutine CHO_QUALIFY_1(DIAG,ISYM,ISHLAB,MEM,MEM0,LEFT)
!
! Purpose: qualify diagonals ("qualify until full").

use Cholesky, only: DiaMin, iiBstR, iiBstRSh, IndRed, iOffq, iQuAB, MaxQual, nnBstR, nnBstRSh, nQual
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: Diag(*)
integer(kind=iwp), intent(in) :: ISYM, ISHLAB, MEM
integer(kind=iwp), intent(inout) :: MEM0, LEFT
integer(kind=iwp) :: I, I2, J, MAXQ, NUMQ

if (NNBSTRSH(ISYM,ISHLAB,2) > 0) then
  I = IIBSTR(ISYM,2)+IIBSTRSH(ISYM,ISHLAB,2)
  I2 = I+NNBSTRSH(ISYM,ISHLAB,2)
  MAXQ = min(MAXQUAL-NQUAL(ISYM),LEFT/NNBSTR(ISYM,2))
  NUMQ = 0
  do while ((I < I2) .and. (NUMQ < MAXQ))
    I = I+1
    J = INDRED(I,2)
    if (DIAG(J) >= DIAMIN(ISYM)) then
      NUMQ = NUMQ+1
      iQuAB(IOFFQ(ISYM)+NUMQ,ISYM) = I
    end if
  end do
  NQUAL(ISYM) = NQUAL(ISYM)+NUMQ
  MEM0 = MEM0+NUMQ*NNBSTR(ISYM,2)
  LEFT = MEM-MEM0
end if

end subroutine CHO_QUALIFY_1

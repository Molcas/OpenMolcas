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

subroutine CHO_SETMAXSHL(DIAG,DIASH,ISYSH,IRED)
!
! Purpose: set max. shell pair data for selection procedure.

use ChoArr, only: iSP2F, iAtomShl
use ChoSwp, only: nnBstRSh, iiBstRSh, IndRed
use Constants

implicit real*8(a-h,o-z)
real*8 DIAG(*), DIASH(*)
integer ISYSH(*)
#include "cholesky.fh"
character*13 SECNAM
parameter(SECNAM='CHO_SETMAXSHL')

! Initialize the largest diagonal in each shell pair.
! ---------------------------------------------------

DIASH(1:NNSHL) = Zero
ISYSH(1:NNSHL) = 0

! Find largest diagonal in each shell pair. Loop only
! over those that are included in the reduced set at hand.
! --------------------------------------------------------

if (IRED == 1) then
  do ISYMAB=1,NSYM
    do ISHLAB=1,NNSHL
      IAB1 = IIBSTR(ISYMAB,IRED)+IIBSTRSH(ISYMAB,ISHLAB,IRED)+1
      IAB2 = IAB1+NNBSTRSH(ISYMAB,ISHLAB,IRED)-1
      do IAB=IAB1,IAB2
        DIASH(ISHLAB) = max(DIASH(ISHLAB),DIAG(IAB))
        if (DIASH(ISHLAB) == DIAG(IAB)) ISYSH(ISHLAB) = ISYMAB
      end do
    end do
  end do
else if ((IRED == 2) .or. (IRED == 3)) then
  do ISYMAB=1,NSYM
    do ISHLAB=1,NNSHL
      JAB1 = IIBSTR(ISYMAB,IRED)+IIBSTRSH(ISYMAB,ISHLAB,IRED)+1
      JAB2 = JAB1+NNBSTRSH(ISYMAB,ISHLAB,IRED)-1
      do JAB=JAB1,JAB2
        IAB = INDRED(JAB,IRED)
        DIASH(ISHLAB) = max(DIASH(ISHLAB),DIAG(IAB))
        if (DIASH(ISHLAB) == DIAG(IAB)) ISYSH(ISHLAB) = ISYMAB
      end do
    end do
  end do
else
  write(LUPRI,*) SECNAM,': unknown reduced set, IRED = ',IRED
  call CHO_QUIT('Unknown reduced set in '//SECNAM,104)
end if

! Exclude 2-center diagonals (if requested).
! The effect of this is that 2-center diagonals can never be
! qualified; they may still be included in the vectors, though.
! If CHO_NO2CENTER=T, the 2-center diagonals are removed from the
! initial diagonal and we need not worry here.
! ---------------------------------------------------------------

if (CHO_1CENTER .and. (.not. CHO_NO2CENTER)) then
  do ISAB=1,NNSHL
    ISHLAB = ISP2F(ISAB)
    call CHO_INVPCK(ISHLAB,ISHLA,ISHLB,.true.)
    if (IATOMSHL(ISHLA) /= IATOMSHL(ISHLB)) DIASH(ISAB) = 0.0d0
  end do
end if

end subroutine CHO_SETMAXSHL

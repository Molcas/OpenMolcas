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

subroutine STEPVEC(ICL,IOP,NCL,NOP,ISPIN,NORB,IWALK)
! PURPOSE: A SPIN-COUPLED CSF IS SPECIFIED BY NCL CLOSED SHELL
!          AND NOP OPEN SHELL AND OCCUPATION VECTORS ICL AND IOP,
!          RESPECTIVELY. THE SPIN COUPLING IS STORED IN THE
!          VECTOR ISPIN. TRANSLATE THESE DATA INTO THE
!          CORRESPONDING STEP VECTOR.

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ICL(*), IOP(*), NCL, NOP, ISPIN(*), NORB, IWALK(*)
integer(kind=iwp) :: IDELSP, IORB, NXTCL, NXTOP
logical(kind=iwp) :: Test

NXTCL = 1
NXTOP = 1
do IORB=1,NORB

  Test = NXTOP <= NOP
  if (Test) Test = IORB == IOP(NXTOP)

  if ((NXTCL <= NCL) .and. (IORB == ICL(NXTCL))) then

    IWALK(IORB) = 3
    NXTCL = NXTCL+1

  else if (Test) then

    if (ISPIN(NXTOP) == 1) then
      IDELSP = 1
    else
      IDELSP = -1
    end if
    if (IDELSP == 1) then
      IWALK(IORB) = 1
    else if (IDELSP == -1) then
      IWALK(IORB) = 2
    end if
    NXTOP = NXTOP+1

  else

    IWALK(IORB) = 0

  end if

end do

! EXIT

return

end subroutine STEPVEC

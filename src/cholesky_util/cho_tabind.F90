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

function CHO_TABIND(TABLE,LKEY,NTABLE,EOINP,LEOINP,NEOINP,WORD)
!
! Purpose: table lookup.
!          First, try to find WORD in TABLE. If success, return ID,
!          else, check if WORD is a special string
!          in EOINP (if any supplied). If success, return NTABLE+1,
!          else, return -1.

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: CHO_TABIND
integer(kind=iwp), intent(in) :: LKEY, NTABLE, LEOINP, NEOINP
character(len=LKEY), intent(in) :: TABLE(NTABLE), WORD
character(len=LEOINP), intent(in) :: EOINP(NEOINP)
integer(kind=iwp) :: IJUMP, LCMP
logical(kind=iwp) :: Test

! Find entry.
! -----------

if ((LKEY > 0) .and. (NTABLE > 0)) then
  IJUMP = 1
  Test = IJUMP <= NTABLE
  if (Test) Test = TABLE(IJUMP) /= WORD
  do while (Test)
    IJUMP = IJUMP+1
    Test = IJUMP <= NTABLE
    if (Test) Test = TABLE(IJUMP) /= WORD
  end do
  if (IJUMP > NTABLE) then
    if ((LEOINP > 0) .and. (NEOINP > 0)) then
      LCMP = min(LEOINP,LKEY)
      IJUMP = 1
      Test = (IJUMP <= NEOINP)
      if (Test) Test = EOINP(IJUMP)(1:LCMP) /= WORD(1:LCMP)
      do while (Test)
        IJUMP = IJUMP+1
        Test = IJUMP <= NTABLE
        if (Test) Test = EOINP(IJUMP)(1:LCMP) /= WORD(1:LCMP)
      end do
      if (IJUMP > NEOINP) then
        CHO_TABIND = -1
      else
        CHO_TABIND = NTABLE+1
      end if
    else
      CHO_TABIND = -1
    end if
  else
    CHO_TABIND = IJUMP
  end if
else
  CHO_TABIND = -1
end if

end function CHO_TABIND

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

function PIKNAM(LINED,KEYW)
! picks a character item up from a line like
! ITEM_NAME  value  COMMENTS
!
! returns ' ' when it fails to pick the item up.
!
!   LINED        input line to be searched
!   KEYW         keyword

use Definitions, only: iwp

implicit none
character(len=40) :: PIKNAM
character(len=80), intent(in) :: LINED
character(len=40), intent(in) :: KEYW
character(len=80) LINE
integer(kind=iwp) :: I, IK, IL, ISTRT, LK, LL
logical(kind=iwp) :: FOUND

IK = 0
LL = 0
IL = 0

LINE = LINED
LK = 1
do I=40,1,-1
  if (KEYW(I:I) /= ' ') then
    LK = I
    exit
  end if
end do
do I=1,40
  if (KEYW(I:I) /= ' ') then
    IK = I
    exit
  end if
end do

ISTRT = index(LINE,KEYW(IK:LK))
if (ISTRT == 0) then
  ! the keyword has not been found in the line.
  PIKNAM = ' '
  return
end if
ISTRT = ISTRT+LK-IK+1
if (LINE(ISTRT:ISTRT) /= ' ') then
  ! the keyword is not followed by a blank character.
  PIKNAM = ' '
  return
else if (ISTRT >= 80) then
  ! the line finishes after the keyword.
  PIKNAM = ' '
  return
else
  ISTRT = ISTRT+1
end if

FOUND = .false.
do I=ISTRT,80
  if (FOUND) then
    if (LINE(I:I) == ' ') then
      LL = I-1
      exit
    end if
  else
    if (LINE(I:I) /= ' ') then
      FOUND = .true.
      IL = I
    end if
  end if
end do

if (FOUND) then
  PIKNAM = LINE(IL:LL)
  write(LINE,601) LK,LL-IL+1
  !write(u6,LINE) KEYW(:LK),PIKNAM
else
  PIKNAM = ' '
end if

return

601 format("(' PIKNAM:         ',A",I2.2,',1X,A',I2.2,')')

end function PIKNAM

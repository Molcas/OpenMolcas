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

character*40 function PIKNAM(LINED,KEYW)
! picks a character item up from a line like
! ITEM_NAME  value  COMMENTS
!
! returns ' ' when it fails to pick the item up.
!
!   LINED        input line to be searched
!   KEYW         keyword

implicit real*8(A-H,O-Z)
character*40 KEYW
character*80 LINE, LINED
logical FOUND

IK = 0
LL = 0
IL = 0

LINE = LINED
LK = 1
do I=40,1,-1
  if (KEYW(I:I) /= ' ') then
    LK = I
    GO TO 11
  end if
end do
11 continue
do I=1,40
  if (KEYW(I:I) /= ' ') then
    IK = I
    GO TO 13
  end if
end do
13 continue

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
      GO TO 21
    end if
  else
    if (LINE(I:I) /= ' ') then
      FOUND = .true.
      IL = I
    end if
  end if
end do
21 continue

if (FOUND) then
  PIKNAM = LINE(IL:LL)
  write(LINE,601) LK,LL-IL+1
  !write(6,LINE) KEYW(:LK),PIKNAM
else
  PIKNAM = ' '
end if

return
601 format('('' PIKNAM:         '',A',I2.2,',1X,A',I2.2,')')

end function PIKNAM

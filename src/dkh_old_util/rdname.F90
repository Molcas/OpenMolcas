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

character*40 function RDNAME(LUT,KEYW)
! reads a character item up wherever it is in a file in the form
! ITEM_NAME  value
!
! Returns ' ' when it fails to read the item up.
!
!   LUT       logical unit to be searched
!             LUT>0 for formatted files, LUT<0 for unformatted files
!   KEYW      keyword

integer LU, LUT
character*40 KEYW, value, PIKNAM
character*80 LINE, BLANK
logical UNFORM
external PIKNAM

BLANK = '                                        '
BLANK = BLANK(:40)//'                                        '

if (LUT > 0) then
  UNFORM = .false.
  LU = LUT
else if (LUT < 0) then
  UNFORM = .true.
  LU = -LUT
  ! if not IBM: no error messages for reading short records
else
  LU = 0
  UNFORM = .true.
  write(6,*) 'RdName: LUT=0!'
  call Abend()
end if

rewind LU
1 continue
if (UNFORM) then
  ! unformatted files
  LINE = BLANK
  read(LU,end=999) LINE
else
  ! formatted files
  read(LU,'(A80)',end=999) LINE
end if
value = PIKNAM(LINE,KEYW)
if (value == ' ') GO TO 1

! the item has been found
RDNAME = value
return

! the item has not been found in the file
999 continue
RDNAME = value

return

end function RDNAME

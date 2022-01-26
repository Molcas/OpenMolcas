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

function RDNAME(LUT,KEYW)
! reads a character item up wherever it is in a file in the form
! ITEM_NAME  value
!
! Returns ' ' when it fails to read the item up.
!
!   LUT       logical unit to be searched
!             LUT>0 for formatted files, LUT<0 for unformatted files
!   KEYW      keyword

use Definitions, only: iwp, u6

implicit none
character(len=40) :: RDNAME
integer(kind=iwp), intent(in) :: LUT
character(len=40), intent(in) :: KEYW
integer(kind=iwp) :: LU
logical(kind=iwp) :: UNFORM
character(len=80) :: LINE
character(len=40) :: VAL
character(len=40), external :: PIKNAM

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
  write(u6,*) 'RdName: LUT=0!'
  call Abend()
end if

rewind LU
1 continue
if (UNFORM) then
  ! unformatted files
  LINE = ''
  read(LU,end=999) LINE
else
  ! formatted files
  read(LU,'(A80)',end=999) LINE
end if
VAL = PIKNAM(LINE,KEYW)
if (VAL == ' ') GO TO 1

! the item has been found
RDNAME = VAL
return

! the item has not been found in the file
999 continue
RDNAME = VAL

return

end function RDNAME

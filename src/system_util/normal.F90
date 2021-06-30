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

subroutine normal(line)

use Definitions, only: iwp

implicit none
character(len=*), intent(inout) :: line
character :: ch
integer(kind=iwp), save :: ifirst = 1, inew(0:255)
integer(kind=iwp) :: i, iflag, iold, ipos
character(len=*), parameter :: upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', &
                               lower = 'abcdefghijklmnopqrstuvwxyz'

if (ifirst == 1) then
  ifirst = 0
  do i=0,255
    inew(i) = i
  end do
  do i=1,26
    iold = ichar(lower(i:i))
    inew(iold) = ichar(upper(i:i))
  end do
end if

ipos = 0
iflag = 1
do i=1,len(line)
  ch = line(i:i)
  if (ch /= ' ') then
    ipos = ipos+1
    iflag = 0
    line(ipos:ipos) = char(inew(ichar(ch)))
  else
    if (iflag == 0) then
      ipos = ipos+1
      line(ipos:ipos) = ' '
      iflag = 1
    end if
  end if
end do
do i=ipos+1,len(line)
  line(i:i) = ' '
end do

return

end subroutine normal

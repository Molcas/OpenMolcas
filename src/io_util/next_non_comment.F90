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

function next_non_comment(lu,line)

use Definitions, only: iwp

implicit none
logical(kind=iwp) :: next_non_comment
integer(kind=iwp), intent(in) :: lu
character(len=*), intent(out) :: line
integer(kind=iwp) :: stat

next_non_comment = .false.
do
  read(lu,'(A)',iostat=stat) line
  if (stat < 0) return
  if (stat > 0) call abend()
  line = adjustl(line)
  if ((line(1:1) /= '*') .and. (line /= ' ')) exit
end do
next_non_comment = .true.

return

end function next_non_comment

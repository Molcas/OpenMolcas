!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine strip_blanks_cvb(line,blanks,nblank,blankdelim)
! If BLANKDELIM, a given number of blanks (not leading or trailing)
! will be subsituted by a single blank. Otherwise all blanks are stripped.

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
character(len=*), intent(inout) :: line
integer(kind=iwp), intent(in) :: nblank
character, intent(in) :: blanks(nblank)
logical(kind=iwp), intent(in) :: blankdelim
integer(kind=iwp) :: iblank, ich, ich2, lenline
integer(kind=iwp), allocatable :: lv(:)

lenline = len_trim(line)
do iblank=1,nblank
  if (blanks(iblank) /= ' ') then
    do ich=1,lenline
      if (line(ich:ich) == blanks(iblank)) line(ich:ich) = ' '
    end do
  end if
end do
call mma_allocate(lv,lenline,label='lv')
ich2 = 0
do ich=1,lenline
  if (line(ich:ich) /= ' ') then
    ich2 = ich2+1
    lv(ich2) = ich
  else if (blankdelim .and. (ich >= 2) .and. (line(ich-1:ich-1) /= ' ') .and. (ich2 /= 0)) then
    ! (Final condition eliminates leading blanks:)
    ich2 = ich2+1
    lv(ich2) = ich
  end if
end do
do ich=1,ich2
  line(ich:ich) = line(lv(ich):lv(ich))
end do
lenline = ich2
call mma_deallocate(lv)

return

end subroutine strip_blanks_cvb

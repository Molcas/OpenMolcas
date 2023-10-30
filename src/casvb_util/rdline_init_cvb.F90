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

subroutine rdline_init_cvb(variat)

use casvb_global, only: inp, lenline, line
use Definitions, only: iwp, u6

implicit none
logical(kind=iwp), intent(in) :: variat
integer(kind=iwp) :: istatus
logical(kind=iwp), parameter :: blankdelim = .true. ! BLANKDELIM signifies whether blanks are used to delimit fields
integer(kind=iwp), parameter :: nblank = 2
character(len=*), parameter :: blanks(nblank) = [' ',',']

if (variat) return
rewind(inp)
do
  read(inp,'(a)',iostat=istatus) line
  if (istatus < 0) then
    write(u6,*) ' WARNING: Initiation string not found in input file.'
    return
  end if
  call strip_blanks_cvb(line,blanks,nblank,blankdelim)
  call upcase(line)
  lenline = len_trim(line)
  if (line(1:6) == '&CASVB') exit
  !if (((.not. blankdelim) .and. (line(1:10) == '&CASVB&END')) .or. (blankdelim .and. (line(1:11) == '&CASVB &END'))) exit
end do

return

end subroutine rdline_init_cvb

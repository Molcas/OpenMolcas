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

subroutine PrintLine(unt,line,length,isBinLuscus)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************

use grid_it_globals, only: iBinary, isLuscus
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: unt, length
logical(kind=iwp), intent(in) :: isBinLuscus
character(len=128) :: line
#include "WrkSpc.fh"
integer(kind=iwp) :: ll, li

if (isLuscus) then
  ll = length
  li = merge(1,0,isBinLuscus)
  !write(u6,*) 'before pl ',line,' ',ll,li
  call prt_lusc(unt,line,ll,li)
else
  if (iBinary == 1) then
    write(unt) line(1:length)
  else
    write(unt,'(A)') line(1:length)
  end if
end if

return

end subroutine PrintLine

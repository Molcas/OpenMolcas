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
! Copyright (C) Valera Veryazov                                        *
!***********************************************************************

subroutine NoBlanks(strout,strin)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************

use Definitions, only: iwp

implicit none
character(len=*), intent(out) :: strout
character(len=*), intent(in) :: strin
integer(kind=iwp) :: flag, i, j

flag = -1
j = 0
do i=1,len(strin)
  if ((flag == -1) .and. (strin(i:i) == ' ')) cycle
  flag = 0
  if (i <= len(strin)-1) then
    if (strin(i:i+1) == '  ') cycle
  end if
  j = j+1
  strout(j:j) = strin(i:i)
end do
strout(j+1:) = ' '

return

end subroutine NoBlanks

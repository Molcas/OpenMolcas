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
! Copyright (C) 1990, IBM                                              *
!***********************************************************************

function iCLast(KWord,iChrct)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iCLast
character(len=*), intent(in) :: KWord
integer(kind=iwp), intent(in) :: iChrct
integer(kind=iwp) :: i

iCLast = 0
do i=iChrct,1,-1
  if (KWord(i:i) /= ' ') then
    iCLast = i
    return
  end if
end do

end function iCLast

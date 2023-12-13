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
! Copyright (C) 2004, Giovanni Ghigo                                   *
!***********************************************************************
!***********************************************************************
! This routine checks if -string- is a Real*8 number. If it is a number*
! NaN=.False. and -Number- contains the number, otherwise, NaN=.True.  *
!----------------------------------------------------------------------*
! Author:   Giovanni Ghigo - November 2004 - Lund(SE)                  *
! Author:   Giovanni Ghigo                                             *
!***********************************************************************

subroutine Get_dNumber(string,dNumber,iErr)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
character(len=*), intent(in) :: string
real(kind=wp), intent(out) :: dNumber
integer(kind=iwp), intent(out) :: iErr
integer(kind=iwp) :: i, j
logical(kind=iwp) :: NaN
character(len=*), parameter :: Chars = ' +-1234567890.'

iErr = 0
dNumber = Zero
NaN = .true.
do i=1,len(string)
  NaN = .true.
  do j=1,len(Chars)
    if (string(i:i) == Chars(j:j)) NaN = .false.
  end do
  if (NaN) exit
end do
if (NaN) then
  iErr = 1
else
  read(string,*) dNumber
end if

return

end subroutine Get_dNumber

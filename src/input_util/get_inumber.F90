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
! This routine checks if -string- is an integer. If it is a number     *
! NaN=.False. and -iNumber- contains the number, otherwise, NaN=.True. *
!----------------------------------------------------------------------*
! Author:   Giovanni Ghigo - November 2004 - Lund(SE)                  *
! Author:   Giovanni Ghigo                                             *
!***********************************************************************

subroutine Get_iNumber(string,iNumber,iErr)

implicit integer(i-n)
character*(*) string
character*11 Chars
data Chars/' 1234567890'/
logical NaN

iErr = 0
iNumber = 0
NaN = .true.
do i=1,len(string)
  NaN = .true.
  do j=1,11
    if (string(i:i) == Chars(j:j)) NaN = .false.
  end do
  if (NaN) goto 10
end do
10 if (.not. NaN) read(string,*) iNumber
if (NaN) iErr = 1

return

end subroutine Get_iNumber

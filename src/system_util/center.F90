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
! Copyright (C) 1993, Markus P. Fuelscher                              *
!***********************************************************************
!  Center
!
!> @brief
!>   Center a string for printing purpose
!> @author M. P. F&uuml;lscher
!> @author P. O. Widmark
!>
!> @details
!> Add spaces to the beginning and the end of a string
!>
!> @param[in,out] String a string
!***********************************************************************

subroutine Center(String)

use Definitions, only: iwp

implicit none
character(len=*), intent(inout) :: String
integer(kind=iwp) :: i, lLeft, lRight, lShift, lString

!----------------------------------------------------------------------*
! get the length of the line                                           *
!----------------------------------------------------------------------*
lString = len(String)
!----------------------------------------------------------------------*
! get the number of leading blanks                                     *
!----------------------------------------------------------------------*
lLeft = 0
do i=1,lString
  if (String(i:i) /= ' ') then
    lLeft = i-1
    exit
  end if
end do
!----------------------------------------------------------------------*
! get the number of trailing blanks                                    *
!----------------------------------------------------------------------*
lRight = 0
do i=lString,1,-1
  if (String(i:i) /= ' ') then
    lRight = lString-i
    exit
  endif
end do
!----------------------------------------------------------------------*
! shift the line                                                       *
!----------------------------------------------------------------------*
if (lLeft+lRight /= 0) then
  lShift = (lRight-lLeft)/2
  if (lShift > 0) then
    do i=lString,lShift+1,-1
      String(i:i) = String(i-lShift:i-lShift)
    end do
    do i=1,lLeft+lShift
      String(i:i) = ' '
    end do
  else if (lShift < 0) then
    do i=1,lString-lShift
      String(i:i) = String(i-lShift:i-lShift)
    end do
    do i=lString,lString-lRight-lShift,-1
      String(i:i) = ' '
    end do
  end if
end if
!----------------------------------------------------------------------*
! normal termination                                                   *
!----------------------------------------------------------------------*
return

end subroutine Center

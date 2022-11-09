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
!               1993, Per-Olof Widmark                                 *
!***********************************************************************

function StrnLn(String)
!***********************************************************************
!                                                                      *
!     Determine the position of the last nonblank character in         *
!     the input string.                                                *
!                                                                      *
!     calling arguments:                                               *
!     String : Character string                                        *
!                                                                      *
!     return value:                                                    *
!     StrnLn  : Integer, position of last nonblank character           *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher and P.O. Widmark                                  *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: StrnLn
character(len=*), intent(in) :: String
integer(kind=iwp) :: i, lString
character, parameter :: cBlank = char(32), cNull = char(0)

!----------------------------------------------------------------------*
! get the length of the line                                           *
!----------------------------------------------------------------------*
lString = len(String)
!----------------------------------------------------------------------*
! loop over the string                                                 *
!----------------------------------------------------------------------*
StrnLn = 0
do i=1,lString
  if ((String(i:i) /= cBlank) .and. (String(i:i) /= cNull)) StrnLn = i
end do

!----------------------------------------------------------------------*
! Normal exit                                                          *
!----------------------------------------------------------------------*
return

end function StrnLn

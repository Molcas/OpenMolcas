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

module getline_mod

private

integer igetline, myunit
integer, parameter :: Len_Line = 180
character(len=Len_Line) Line
integer, parameter :: mxn = Len_Line/2+1
integer iend(mxn), istrt(mxn), ncol
logical Quit_On_Error

public :: igetline, myunit, Line, mxn, iend, istrt, ncol, Quit_On_Error

end module getline_mod

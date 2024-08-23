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
Module getline_mod
Private
Integer igetline,myunit

integer, parameter:: Len_Line=180
character(Len=Len_Line) Line
integer, parameter:: mxn=Len_Line/2+1
integer iend(mxn), istrt(mxn), ncol

Logical Quit_On_Error

Public :: igetline, myunit, Line, mxn, iend, istrt, ncol, Quit_On_Error
End Module getline_mod

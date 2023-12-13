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

subroutine RdNLst(iUnit,NameIn)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iUnit
character(len=*), intent(in) :: NameIn
logical(kind=iwp) :: No_Input_OK

No_Input_OK = .false.
call RdNLst_(iUnit,NameIn,No_Input_OK)

return

end subroutine RdNLst

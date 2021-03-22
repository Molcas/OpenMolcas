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
! Copyright (C) 2013, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine RPA_ReadIntegerInput(Key,nInp,Lu,iVal,n)

implicit none
character*4 Key
integer nInp
integer Lu
integer n
integer iVal(n)
character*180 Line

character*180 Get_Ln
external Get_Ln

if (n >= nInp) then
  Line = Get_Ln(Lu)
  call Get_I(1,iVal,nInp)
else
  ! insufficent memory for reading (fix in calling routine)
  call RPA_Warn(3,'Integer read problem for keyword '//Key)
end if

return
#ifdef _WARNING_WORKAROUND_
if (.false.) call Unused_character(Line)
#endif

end subroutine RPA_ReadIntegerInput

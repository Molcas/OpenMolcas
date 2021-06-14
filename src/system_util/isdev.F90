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
! Copyright (C) 2000-2016, Valera Veryazov                             *
!***********************************************************************
!*****************************************************************************
!                                                                            *
! Author:   Valera Veryazov 2000-2016                                        *
!           Theoretical Chemistry                                            *
!           Lund University                                                  *
!           Sweden                                                           *
!                                                                            *
!*****************************************************************************
!
!   Use this routine wisely!
!
!   in case if you would like to block some dev. code
!   from usage, or apply some conditions for this:
!
!   add a call
!     Call OnlyIMayUseIt("name")

subroutine OnlyIMayUseIt(Name)

character*(*) Name
character*256 value
character*1024 message

value = ' '
call getenvf('MOLCAS_ISDEV',value)
if (value == 'PRODUCTION') then
  return
end if
if (value == ' ' .or. value /= Name) then
  message = 'This code is for testing purpose only;if you want to use it by other means, ;e.g. to publish the results -- you must &
            &contact;'//Name//' to find the conditions applied;;Alternatively you may either use the production version;'
  call WarningMessage(2,message)
  return
end if

end subroutine OnlyIMayUseIt

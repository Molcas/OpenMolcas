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
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_XCV_TmpFiles(irc,iOpt)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: iOpt

irc = 0
if (iOpt == 1) then
  call Cho_XCV_OpenTmpFiles()
else if (iOpt == 2) then
  call Cho_XCV_CloseAndKeepTmpFiles()
else if (iOpt == 3) then
  call Cho_XCV_CloseAndEraseTmpFiles()
else
  irc = 1
end if

return

end subroutine Cho_XCV_TmpFiles

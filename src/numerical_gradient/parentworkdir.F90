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

subroutine ParentWorkDir()

use subdirs, only: f_setsubdir, Sub, OldWorkDir, NewWorkDir
use filesystem, only: remove_, chdir_
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: i

Sub = ''
call f_setsubdir(Sub)
call chdir_(OldWorkDir)

! Remove the subdirectory,
! try to make sure it is indeed a subdirectory
if (index(NewWorkDir,trim(OldWorkDir)) == 1) then
  i = len_trim(OldWorkDir)
  if ((len_trim(NewWorkDir) >= i+2) .and. (NewWorkDir(i+1:i+1) == '/') .and. (NewWorkDir(i+2:i+2) /= '/')) then
    call remove_(NewWorkDir)
  end if
end if

return

end subroutine ParentWorkDir

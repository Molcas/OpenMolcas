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

subroutine fileorb(filein,fileout)

use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: filein
character(len=*), intent(out) :: fileout
character(len=256) :: tmp
logical(kind=iwp) :: Exists

if (index(filein,'/') /= 0) then
  fileout = filein
  goto 100
end if
tmp = ' '
call getenvf('MOLCAS_SUBMIT_DIR',tmp)
if (tmp /= ' ') then
  fileout = trim(tmp)//'/'//filein
  !write(u6,*) 'vv',fileout
  call f_inquire(fileout,Exists)
  if (Exists) goto 100
end if
fileout = filein
call f_inquire(fileout,Exists)
if (.not. Exists) then
  tmp = 'file '//trim(fileout)//' not found'
  call WarningMessage(2,tmp)
  call Quit_OnUserError()
end if
100 continue
!write(u6,*) 'INPORB file=',fileout
!call fcopy(trim(fileout),'INPORB')

return

end subroutine fileorb

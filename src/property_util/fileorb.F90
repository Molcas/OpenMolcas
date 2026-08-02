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
else
  Exists = .false.
  tmp = ' '
  call getenvf('MOLCAS_SUBMIT_DIR',tmp)
  if (tmp /= ' ') then
    fileout = trim(tmp)//'/'//filein
    !write(u6,*) 'vv',fileout
    call f_inquire(fileout,Exists)
  end if
  if (.not. Exists) then
    fileout = filein
    call f_inquire(fileout,Exists)
    if (.not. Exists) then
      tmp = 'file '//trim(fileout)//' not found'
      call WarningMessage(2,tmp)
      call Quit_OnUserError()
    end if
  end if
end if
!write(u6,*) 'INPORB file=',fileout
!call fcopy(trim(fileout),'INPORB')

return

end subroutine fileorb

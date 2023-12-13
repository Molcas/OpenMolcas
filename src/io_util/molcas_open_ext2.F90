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
! Copyright (C) 2001-2005, Valera Veryazov                             *
!***********************************************************************

subroutine molcas_open_ext2(Lu,f_Name,f_access,f_form,f_iostat,is_recl,f_recl,f_status,is_error)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: Lu, f_recl
integer(kind=iwp), intent(out) :: f_iostat
character(len=*), intent(in) :: f_Name, f_access, f_form, f_status
logical(kind=iwp), intent(in) :: is_recl
logical(kind=iwp), intent(out) :: is_error
integer(kind=iwp) :: lRealName
character(len=4096) :: RealName

is_error = .false.
call prgmtranslate(f_Name,RealName,lRealName)
if (index(RealName,'UNK_VAR') /= 0) then
  write(u6,*) '*** attempt to open ',RealName(1:lRealName)
  RealName = f_Name
  lRealName = index(RealName,' ')
end if

if (is_recl) then
  open(unit=Lu,file=Realname(1:lRealName),status=f_status,access=f_access,form=f_form,iostat=f_iostat,recl=f_recl)
else
  open(unit=Lu,file=RealName(1:lRealName),status=f_status,access=f_access,form=f_form,iostat=f_iostat)
end if
!write(u6,*) 'DEBUG open ',RealName(1:lRealName)
!write(u6,*) 'Unit ', Lu

if (f_iostat /= 0) is_error = .true.

return

end subroutine molcas_open_ext2

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

subroutine molcas_open(Lu,FileName)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: Lu
character(len=*), intent(in) :: FileName
integer(kind=iwp) :: f_recl, f_iostat
character(len=10) :: f_access, f_form, f_status
logical(kind=iwp) :: is_recl, is_error

f_recl = 1
f_iostat = 100
f_access = 'SEQUENTIAL'
f_form = 'FORMATTED'
f_status = 'UNKNOWN'
is_recl = .false.
is_error = .false.

call molcas_open_ext2(Lu,trim(FileName),f_access,f_form,f_iostat,is_recl,f_recl,f_status,is_error)
if (f_iostat /= 0) then
  write(u6,*)
  write(u6,'(3a)') 'molcas_open: Error opening file "',trim(FileName),'"'
  write(u6,'(a,i9)') '   iostat is',f_iostat
  write(u6,'(a)') '   Aborting'
  write(u6,*)
  call Abend()
end if

return

end subroutine molcas_open

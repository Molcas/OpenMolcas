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

subroutine vf(fname,lun)
! this routine opens file vanisf file with a given name
! fname- name of the vanished file (I)
! lun  - lun number with which file will be opened (I)

use Definitions, only: iwp

implicit none
character(len=8), intent(in) :: fname
integer(kind=iwp), intent(in) :: lun

call molcas_open(lun,fname)
!open(unit=lun,file=fname)
write(lun,*) ' File scratched'
close(lun)

return

end subroutine vf

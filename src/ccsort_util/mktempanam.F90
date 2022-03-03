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

subroutine mktempanam()
! this routine prepares names for TEMP and files as
! TEMP001 - TEMPmbas and stores them into
! tmpnam and tmanam arrays (mbas-maximum number of basis functions)
!
! variables used:
! tmpnam - array of TEMP file names (imported from ccsort_global)
! this routine (I)

use ccsort_global, only: lunpublic, mbas, tmpnam
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: itemp, k1, lun

lun = lunpublic
call molcas_open(lun,'TEMP000')
!open(unit=lun,file='TEMP000')

do k1=1,mbas
  if (k1 < 10) then
    write(lun,99) k1
  else if (k1 < 100) then
    write(lun,199) k1
  else
    write(lun,299) k1
  end if
end do

rewind(lun)

do itemp=1,mbas
  read(lun,599) tmpnam(itemp)
end do

rewind(lun)
write(lun,*) ' File scratched'
close(lun)

return

99 format('TEMP00',i1)
199 format('TEMP0',i2)
299 format('TEMP',i3)
599 format(a7)

end subroutine mktempanam

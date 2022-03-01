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
! tmp-anam - array of TEMP file names (Transported through common /tmnames/)
! this routine (I)

implicit real*8(a-h,o-z)
#include "reorg.fh"
integer lun, itemp, k1

lun = lunpublic
call molcas_open(lun,'TEMP000')
!open(unit=lun,file='TEMP000')

itemp = 0
do k1=1,9
  itemp = itemp+1
  if (itemp > mbas) goto 500
  write(lun,99) k1
99 format(6hTEMP00,i1)
end do

do k1=10,99
  itemp = itemp+1
  if (itemp > mbas) goto 500
  write(lun,199) k1
199 format(5hTEMP0,i2)
end do

do k1=100,mbas
  itemp = itemp+1
  if (itemp > mbas) goto 500
  write(lun,299) k1
299 format(4hTEMP,i3)
end do

500 rewind(lun)

do itemp=1,mbas
  read(lun,599) tmpnam(itemp)
599 format(a7)
end do

rewind(lun)
write(lun,*) ' File scratched'
close(lun)

return

end subroutine mktempanam

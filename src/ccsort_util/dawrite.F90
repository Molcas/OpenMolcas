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

subroutine dawrite(lun,irec0,vector,length,recl)
! this routine writes vector with required length to
! open direct access file lun starting from record number
! irec0
!
! lun   - logical unit of direct access file (I)
! irec0 - initial record number (I)
! vector- vector (I)
! length- number of R8 data to be read (I)
! recl  - length of one record in lun in R8 (I)

real*8 vector(1:length)
integer lun, irec0, length, recl
! help variables
integer ilow, iup, need, irec, i

if (length == 0) return

! def need,ilow,iup,irec

need = length
ilow = 1
irec = irec0
iup = 0

do
  if (recl >= need) then
    iup = iup+need
  else
    iup = iup+recl
  end if

  write(lun,rec=irec) (vector(i),i=ilow,iup)

  need = need-(iup-ilow+1)
  irec = irec+1
  ilow = ilow+recl

  if (need <= 0) exit
end do

return

end subroutine dawrite

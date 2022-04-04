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

subroutine daread(lun,irec0,vector,length,reclen)
! this routine reads vector with required length from
! open direct access file lun starting from record number
! irec0
! lun   - logical unit of direct access file (I)
! irec0 - initial record number (I)
! vector- vector (O)
! length- number of R8 data to be read (I)
! reclen- length of one record in lun in R8 (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lun, irec0, length, reclen
real(kind=wp), intent(out) :: vector(length)
integer(kind=iwp) :: ilow, irec, iup, need

if (length == 0) return

! def need,ilow,iup,irec

need = length
ilow = 1
iup = 0
irec = irec0

do
  if (reclen >= need) then
    iup = iup+need
  else
    iup = iup+reclen
  end if

  read(lun,rec=irec) vector(ilow:iup)

  need = need-(iup-ilow+1)
  irec = irec+1
  ilow = ilow+reclen

  if (need <= 0) exit
end do

return

end subroutine daread

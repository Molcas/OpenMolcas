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

subroutine cct3_getmap(lun,med,length,rc)
! this routine reads med%d and med%i of the given mediate
! from lun and reconstructs med%d to actual positions med%pos0
!
! lun    - Logical unit number of file, where mediate is stored (Input)
! med    - mediate (Input/Output)
! length - overall length of mediate (Output)
! rc     - return (error) code (Output)
!
! N.B.
! all mediates are stored as follows
! 1 - %d, %i
! 2 - one record with complete mediate

use CCT3_global, only: daddr, iokey, Map_Type
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lun
type(Map_Type), intent(inout) :: med
integer(kind=iwp), intent(out) :: length, rc
integer(kind=iwp) :: im, pos

rc = 0

!1 read med%d

if (iokey == 1) then
  ! Fortran IO
  read(lun) med%d,med%i

else
  ! MOLCAS IO
  call idafile(lun,2,med%d,513*6,daddr(lun))
  call idafile(lun,2,med%i,8*8*8,daddr(lun))
end if

!2 change positions in med%d to proper one and calculate overall length

pos = med%pos0
length = 0

do im=1,med%d(0,5)

  med%d(im,1) = pos
  pos = pos+med%d(im,2)
  length = length+med%d(im,2)

end do

return

end subroutine cct3_getmap

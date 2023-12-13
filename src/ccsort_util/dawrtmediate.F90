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

subroutine dawrtmediate(wrk,wrksize,lun,map,rc)
! this routine writes required mediate to open unformatted file
! with number lun
! it also stores %d and %i of the given mediate
!
! lun - Logical unit number of file, where mediate will be stored (Input)
! map - map type corresponding to given mediate (Input)
! rc  - return (error) code (Output)
!
! N.B.
! all mediates are stored as follows
! 1 - Map_Type
! 2 - one record with complete mediate

use ccsort_global, only: Map_Type
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize, lun
real(kind=wp), intent(_IN_) :: wrk(wrksize)
type(Map_Type), intent(_IN_) :: map
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp) :: im, length, pos0

rc = 0

!1 write map%d

call dawrtmap(lun,map,rc)

!2 calculate overall length

length = 0

do im=1,map%d(0,5)
  length = length+map%d(im,2)
end do

! write mediate in one block

if (length == 0) then
  ! RC=1 : there is nothing to write, length of mediate is 0
  rc = 1
else
  pos0 = map%d(1,1)
  call dawri(lun,length,wrk(pos0))
end if

return

end subroutine dawrtmediate

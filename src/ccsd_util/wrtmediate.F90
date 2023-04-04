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

subroutine wrtmediate(wrk,wrksize,lun,map,rc)
! this routine writes required mediate to opened unformatted file
! with number lun
! it also stores %d and %i of the given mediate
!
! lun - Logical unit number of file, where mediate will be stored (Input)
! map - given mediate (Input)
! rc  - return (error) code (Output)
!
! N.B.
! all mediates are stored as follows
! 1 - Map_Type
! 2 - one record with complete mediate

use ccsd_global, only: Map_Type
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: wrksize, lun, rc
real(kind=wp) :: wrk(wrksize)
type(Map_Type) :: map
integer(kind=iwp) :: im, length, poss0, rc1

rc = 0

!1 write map

call wrtmap(lun,map,rc1)

!2 calculate overall length

length = 0

do im=1,map%d(0,5)
  length = length+map%d(im,2)
end do

! write mediate in one block

if (length == 0) then
  ! RC=1 : there is nothing to write, length of mediate is 0
  rc = 1
  return
end if

poss0 = map%d(1,1)
call wri(lun,length,wrk(poss0))

return

end subroutine wrtmediate

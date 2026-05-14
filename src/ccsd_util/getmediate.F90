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

subroutine getmediate(wrk,wrksize,lun,map,rc)
! this routine reads required mediate from opened unformatted file
! with number lun, and places it starting with the %pos0
! it also reads %d and %i of the given mediate, and reconstructs
! %d to actual positions
!
! lun - Logical unit number of file, where mediate is stored (Input)
! map - map type corresponding to given mediate (Output)
! rc  - return (error) code (Output)
!
! N.B.
! all mediates are stored as follows
! 1 - Map_Type
! 2 - one record with complete mediate

use ccsd_global, only: Map_Type
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize, lun
real(kind=wp), intent(inout) :: wrk(wrksize)
type(Map_Type), intent(inout) :: map
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp) :: length, rc1

rc = 0

!1 read %d

call getmap(lun,length,map,rc1)

!2 read mediate in one block

if (length == 0) then
  ! RC=1 : there is nothing to read, length of mediate is 0
  rc = 1
  return
end if

call rea(lun,length,wrk(map%pos0))

return

end subroutine getmediate

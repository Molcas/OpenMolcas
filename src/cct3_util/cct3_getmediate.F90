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

subroutine cct3_getmediate(wrk,wrksize,lun,med,rc)
! this routine reads required mediate from opened unformatted file
! with number lun, and places it starting with the %pos0
! it also reads %d and %i of the given mediate, and reconstructs
! %d to actual positions
!
! lun  - Logical unit number of file, where mediate is stored (Input)
! med  - mediate (Input/Output)
! rc   - return (error) code (Output)
!
! N.B.
! all mediates are stored as follows
! 1 - %d, %i
! 2 - one record with complete mediate

use CCT3_global, only: Map_Type
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize, lun
real(kind=wp), intent(_OUT_) :: wrk(wrksize)
type(Map_Type), intent(inout) :: med
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp) :: length, rc1

rc = 0

!1 read med%d

call cct3_getmap(lun,med,length,rc1)

!2 read mediate in one block

if (length == 0) then
  ! RC=1 : there is nothing to read, length of mediate is 0
  rc = 1
  return
end if

call cct3_rea(lun,length,wrk(med%pos0))

return

end subroutine cct3_getmediate

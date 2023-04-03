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

subroutine getmediate(wrk,wrksize,lun,poss0,mapd,mapi,rc)
! this routine reads required mediate from opened unformatted file
! with number lun, and places it starting with the poss0
! it also reads mapd and mapi of the given mediade, and reconstructs
! mapd to actual positions
!
! lun   - Logical unit number of file, where mediate is stored (Input)
! poss0 - initial position in WRK, where mediate will be stored (Input)
! mapd  - direct map matrix corresponding to given mediate (Output)
! mapi  - inverse map matrix corresponding to given mediate (Output)
! rc    - return (error) code (Output)
!
! N.B.
! all mediates are stored as follows
! 1 - mapd, mapi
! 2 - one record with complete mediate

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: wrksize, lun, poss0, mapd(0:512,6), mapi(8,8,8), rc
real(kind=wp) :: wrk(wrksize)
integer(kind=iwp) :: length, rc1

rc = 0

!1 read mapd

call getmap(lun,poss0,length,mapd,mapi,rc1)

!2 read mediate in one block

if (length == 0) then
  ! RC=1 : there is nothing to read, length of mediate is 0
  rc = 1
  return
end if

call rea(lun,length,wrk(poss0))

return

end subroutine getmediate

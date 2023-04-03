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

subroutine wrtmap(lun,mapd,mapi,rc)
! this routine writes required mapd and mapi to opened unformatted file
! with number lun
!
! lun   - Logical unit number of file, where mediate will be stored (Input)
! mapd  - direct map matrix corresponding to given mediate (Input)
! mapi  - inverse map matrix corresponding to given mediate (Input)
! rc    - return (error) code (Output)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: lun, mapd(0:512,6), mapi(8,8,8), rc
#include "filemgr.fh"
#include "ccsd1.fh"

rc = 0

!1 write mapd

if (iokey == 1) then
  ! Fortran IO
  write(lun) mapd,mapi

else
  ! MOLCAS IO
  call idafile(lun,1,mapd,513*6,daddr(lun))
  call idafile(lun,1,mapi,8*8*8,daddr(lun))
end if

return

end subroutine wrtmap

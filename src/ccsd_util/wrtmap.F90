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

subroutine wrtmap(lun,map,rc)
! this routine writes required map to opened unformatted file
! with number lun
!
! lun - Logical unit number of file, where mediate will be stored (Input)
! map - given mediate (Input)
! rc  - return (error) code (Output)

use ccsd_global, only: daddr, iokey, Map_Type
use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: lun
type(Map_Type), intent(_IN_) :: map
integer(kind=iwp), intent(out) :: rc

rc = 0

!1 write map

if (iokey == 1) then
  ! Fortran IO
  write(lun) map%d,map%i

else
  ! MOLCAS IO
  call idafile(lun,1,map%d,size(map%d),daddr(lun))
  call idafile(lun,1,map%i,size(map%i),daddr(lun))
end if

return

end subroutine wrtmap

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

subroutine dawri(lun,length,vector)
! this routine writes length-R8 numbers to open unformatted file
! with number lun at the given position as one record
!
! lun    - Logical unit number of file, where mediate will be stored (Input)
! length - # of R8 numbers to be written (Input)
! vector - space, where numbers are stored (Input)

use ccsort_global, only: daddr, iokey
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: lun, length
real(kind=wp), intent(_IN_) :: vector(length)

if (iokey == 1) then
  ! Fortran IO
  write(lun) vector

else
  ! MOLCAS IO
  call ddafile(lun,1,vector,length,daddr(lun))
end if

return

end subroutine dawri

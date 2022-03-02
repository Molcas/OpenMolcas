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

subroutine ccsort_wri(lun,length,vector)
! this routine writes length-R8 numbers to open unformatted file
! with number lun at the given position as one record
!
! lun    - Logical unit number of file, where mediate will be stored (Input)
! length - # of R8 numbers to be written (Input)
! vector - space, where numbers are stored Input)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: lun, length
real(kind=wp) :: vector(length)

write(lun) vector

return

end subroutine ccsort_wri

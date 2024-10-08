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

subroutine Get_I1(icol,ival)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iCol
integer(kind=iwp), intent(out) :: ival
integer(kind=iwp) :: idum(1)

call Get_I(icol,idum,1)
ival = idum(1)

end subroutine Get_I1

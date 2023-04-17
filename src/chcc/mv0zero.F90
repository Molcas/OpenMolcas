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

subroutine mv0zero(dd,length,mat)
! mat = 0

use chcc_global, only: mhkey
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: dd, length
real(kind=wp) :: mat(dd)

if (mhkey == 1) then
  ! ESSL

  call dcopy_(length,[Zero],0,mat,1)

else
  ! Fortran matrix handling

  mat(1:length) = Zero

end if

return

end subroutine mv0zero

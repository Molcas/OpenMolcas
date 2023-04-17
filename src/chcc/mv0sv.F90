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

subroutine mv0sv(dd,length,mat,f)
! A = A .f

use chcc_global, only: mhkey
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: dd, length
real(kind=wp) :: mat(dd), f
integer(kind=iwp) :: init

if (mhkey == 1) then
  ! ESSL

  do init=1,length
    mat(init) = mat(init)*f
  end do

else
  ! Fortran matrix handling

  do init=1,length
    mat(init) = mat(init)*f
  end do

end if

return

end subroutine mv0sv

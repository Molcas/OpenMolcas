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

subroutine sleepf(sec)

use, intrinsic :: iso_c_binding, only: c_int
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: sec
integer(kind=c_int) :: r
interface
  function sleepc(seconds) bind(C,name='sleep')
    import :: c_int
    integer(c_int) :: sleepc
    integer(c_int), intent(in), value :: seconds
  end function sleepc
end interface

r = sleepc(int(sec,kind=c_int))

#include "macros.fh"
unused_var(r)

end subroutine sleepf

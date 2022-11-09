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

#include "macros.fh"

subroutine Unused_real(r)

use Definitions, only: wp

implicit none
real(kind=wp), intent(in) :: r

unused_var(r)

end subroutine Unused_real

!***********************************************************************

subroutine Unused_real_array(r)

use Definitions, only: wp

implicit none
real(kind=wp), intent(in) :: r(*)

unused_var(r(1))

end subroutine Unused_real_array

!***********************************************************************

subroutine Unused_integer(i)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: i

unused_var(i)

end subroutine Unused_integer

!***********************************************************************

subroutine Unused_integer_array(i)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: i(*)

unused_var(i(1))

end subroutine Unused_integer_array

!***********************************************************************

subroutine Unused_logical(l)

use Definitions, only: iwp

implicit none
logical(kind=iwp), intent(in) :: l

unused_var(l)

end subroutine Unused_logical

!***********************************************************************

subroutine Unused_character(c)

implicit none
character, intent(in) :: c(*)

unused_var(c(1))

end subroutine Unused_character

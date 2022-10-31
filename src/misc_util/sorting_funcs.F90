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

module sorting_funcs

use Definitions, only: wp, iwp

implicit none
private

public :: geq_c, geq_i, geq_r, leq_c, leq_i, leq_r

contains

pure function leq_i(x,y)
  logical(kind=iwp) :: leq_i
  integer(kind=iwp), intent(in) :: x, y
  leq_i = x <= y
end function leq_i

pure function leq_r(x,y)
  logical(kind=iwp) :: leq_r
  real(kind=wp), intent(in) :: x, y
  leq_r = x <= y
end function leq_r

pure function leq_c(x,y)
  logical(kind=iwp) :: leq_c
  complex(kind=wp), intent(in) :: x, y
  leq_c = real(x)**2+aimag(x)**2 <= real(y)**2+aimag(y)**2
end function leq_c

pure function geq_i(x,y)
  logical(kind=iwp) :: geq_i
  integer(kind=iwp), intent(in) :: x, y
  geq_i = x >= y
end function geq_i

pure function geq_r(x,y)
  logical(kind=iwp) :: geq_r
  real(kind=wp), intent(in) :: x, y
  geq_r = x >= y
end function geq_r

pure function geq_c(x,y)
  logical(kind=iwp) :: geq_c
  complex(kind=wp), intent(in) :: x, y
  geq_c = real(x)**2+aimag(x)**2 >= real(y)**2+aimag(y)**2
end function geq_c

end module sorting_funcs

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
implicit none
public

interface le
  module procedure le_i, le_r, le_c
end interface

interface ge
  module procedure ge_i, ge_r, ge_c
end interface

contains

logical pure function le_i(x, y)
  integer, intent(in) :: x, y
  le_i = x <= y
end function

logical pure function le_r(x, y)
  real*8, intent(in) :: x, y
  le_r = x <= y
end function

logical pure function le_c(x, y)
  complex*16, intent(in) :: x, y
  le_c = real(x)**2+aimag(x)**2 <= real(y)**2+aimag(y)**2
end function

logical pure function ge_i(x, y)
  integer, intent(in) :: x, y
  ge_i = x >= y
end function

logical pure function ge_r(x, y)
  real*8, intent(in) :: x, y
  ge_r = x >= y
end function

logical pure function ge_c(x, y)
  complex*16, intent(in) :: x, y
  ge_c = real(x)**2+aimag(x)**2 >= real(y)**2+aimag(y)**2
end function

end module sorting_funcs

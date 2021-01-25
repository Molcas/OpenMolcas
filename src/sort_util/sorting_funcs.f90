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

contains

logical pure function leq_i(x, y)
  integer, intent(in) :: x, y
  leq_i = x <= y
end function

logical pure function leq_r(x, y)
  real*8, intent(in) :: x, y
  leq_r = x <= y
end function

logical pure function leq_c(x, y)
  complex*16, intent(in) :: x, y
  leq_c = real(x)**2+aimag(x)**2 <= real(y)**2+aimag(y)**2
end function

logical pure function geq_i(x, y)
  integer, intent(in) :: x, y
  geq_i = x >= y
end function

logical pure function geq_r(x, y)
  real*8, intent(in) :: x, y
  geq_r = x >= y
end function

logical pure function geq_c(x, y)
  complex*16, intent(in) :: x, y
  geq_c = real(x)**2+aimag(x)**2 >= real(y)**2+aimag(y)**2
end function

end module sorting_funcs

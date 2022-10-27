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

subroutine Unused_real(r)

implicit none
real*8 r, r2

#ifdef _WARNING_WORKAROUND_
if (.false.) then
  r2 = r
  r = r2
end if
#endif

end subroutine Unused_real

!***********************************************************************

subroutine Unused_real_array(r)

implicit none
real*8 r(*), r2

#ifdef _WARNING_WORKAROUND_
if (.false.) then
  r2 = r(1)
  r(1) = r2
end if
#endif

end subroutine Unused_real_array

!***********************************************************************

subroutine Unused_integer(i)

implicit none
integer i, i2

#ifdef _WARNING_WORKAROUND_
if (.false.) then
  i2 = i
  i = i2
end if
#endif

end subroutine Unused_integer

!***********************************************************************

subroutine Unused_integer_array(i)

implicit none
integer i(*), i2

#ifdef _WARNING_WORKAROUND_
if (.false.) then
  i2 = i(1)
  i(1) = i2
end if
#endif

end subroutine Unused_integer_array

!***********************************************************************

subroutine Unused_logical(l)

implicit none
logical l, l2

#ifdef _WARNING_WORKAROUND_
if (.false.) then
  l2 = l
  l = l2
end if
#endif

end subroutine Unused_logical

!***********************************************************************

subroutine Unused_character(c)

implicit none
character c(*), c2

#ifdef _WARNING_WORKAROUND_
if (.false.) then
  c2 = c(1)
  c(1) = c2
end if
#endif

end subroutine Unused_character

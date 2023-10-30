!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

function isitanint_cvb(a)

use Definitions, only: iwp

implicit none
logical(kind=iwp) :: isitanint_cvb
character(len=*), intent(in) :: a
integer(kind=iwp) :: ich, j
logical(kind=iwp) :: done
integer(kind=iwp), parameter :: nallowed = 12
character, parameter :: allowedchars(nallowed) = ['+','-','0','1','2','3','4','5','6','7','8','9']

done = .false.
do ich=1,len_trim(a)
  do j=1,nallowed
    if (a(ich:ich) == allowedchars(j)) then
      done = .true.
      exit
    end if
  end do
  if (.not. done) then
    isitanint_cvb = .false.
    return
  end if
end do
isitanint_cvb = .true.

return

end function isitanint_cvb

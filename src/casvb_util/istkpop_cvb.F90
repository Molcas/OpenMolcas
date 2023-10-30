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

subroutine istkpop_cvb(iarr,ival)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: iarr(*)
integer(kind=iwp), intent(out) :: ival

if (iarr(2) == 2) then
  write(u6,*) ' Trying to pop off empty stack!'
  call abend_cvb()
end if
ival = iarr(iarr(2))
iarr(2) = iarr(2)-1

return

end subroutine istkpop_cvb

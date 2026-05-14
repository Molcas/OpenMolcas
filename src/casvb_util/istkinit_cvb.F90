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
!******************************
!** Simple stack of integers **
!******************************

subroutine istkinit_cvb(iarr,n)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n
integer(kind=iwp), intent(out) :: iarr(2)

if (n < 2) then
  write(u6,*) ' Too small dimension in ISTKINIT_CVB :',n
  call abend_cvb()
end if
iarr(1) = n
iarr(2) = 2

return

end subroutine istkinit_cvb

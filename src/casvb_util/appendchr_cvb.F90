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

subroutine appendchr_cvb(c,string,iskip)

use Definitions, only: iwp

implicit none
character(len=*), intent(inout) :: c
character(len=*), intent(in) :: string
integer(kind=iwp), intent(in) :: iskip
integer(kind=iwp) :: ibegin, iend

ibegin = len_trim(c)+1+iskip
iend = min(len(c),ibegin+len_trim(string)-1)
c(ibegin:iend) = trim(string)

return

end subroutine appendchr_cvb

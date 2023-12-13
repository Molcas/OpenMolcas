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
!**************************************
!** Low-level input parsing routines **
!**************************************

subroutine rdstring_cvb(string,ierr)
! Check if field is applicable:

use casvb_global, only: ifield, nfield
use Definitions, only: wp, iwp

implicit none
character(len=8), intent(inout) :: string
integer(kind=iwp), intent(out) :: ierr
integer(kind=iwp) :: idi
real(kind=wp) :: rdr

ierr = 0
if (nfield == -1) ierr = 1
if (ifield > nfield) ierr = 2
if (ierr /= 0) then
  string = '        '
else
  !call gtstring_cvb(string,ifield)
  call gtany_cvb(string,idi,rdr,1,ifield,ierr)
end if

return

end subroutine rdstring_cvb

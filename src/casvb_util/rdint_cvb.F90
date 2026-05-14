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

subroutine rdint_cvb(intval,ierr)
! Check if field is applicable:

use casvb_global, only: ifield, nfield
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(inout) :: intval
integer(kind=iwp), intent(out) :: ierr
integer(kind=iwp) :: jerr
real(kind=wp) :: rdr
character(len=8) :: string

ierr = 0
if (nfield == -1) ierr = 1
if (ifield > nfield) ierr = 2
if (ierr /= 0) return
call gtany_cvb(string,intval,rdr,2,ifield,jerr)
if (jerr == 1) then
  if (ifield == 1) ierr = 3
  if (ifield /= 1) ierr = 4
end if

return

end subroutine rdint_cvb

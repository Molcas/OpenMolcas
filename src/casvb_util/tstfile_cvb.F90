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

function tstfile_cvb(fileid)

use casvb_global, only: ifilio, filename
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp) :: tstfile_cvb
real(kind=wp), intent(in) :: fileid
integer(kind=iwp) :: ibf
logical(kind=iwp) :: ex

if (fileid < 1.0e-2_wp) then
  tstfile_cvb = .false.
  return
end if
call mkfn_cvb(fileid,ibf)
if (ifilio(ibf) == 0) then
  call f_inquire(filename(ibf),ex)
  tstfile_cvb = ex
else
  tstfile_cvb = .true.
end if

return

end function tstfile_cvb

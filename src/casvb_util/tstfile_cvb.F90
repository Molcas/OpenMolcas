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

logical function tstfile_cvb(fileid)

implicit real*8(a-h,o-z)
logical ex
#include "io_cvb.fh"

if (fileid < 1.d-2) then
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

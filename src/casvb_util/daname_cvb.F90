!s***********************************************************************
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

subroutine daname_cvb(lu,fname)

implicit real*8(a-h,o-z)
character*(*) fname
#include "dancom_cvb.fh"
integer find_lu, isfreeunit, i_open
logical find_unused, is_opened
external find_lu, isfreeunit, is_opened

ilu = find_lu(fname)
if (ilu > 0) then
  lu = ilu
else
  find_unused = (lu < 1)
  if (.not. find_unused) find_unused = is_opened(lu)
  if (find_unused) ilu = isfreeunit(10)
end if

if (is_opened(lu)) then
  i_open = 1
else
  i_open = 0
  lu = 10  ! initialize
end if
call istkpush_cvb(idan,i_open)
if (i_open == 0) call daname(lu,fname)

return

end subroutine daname_cvb
